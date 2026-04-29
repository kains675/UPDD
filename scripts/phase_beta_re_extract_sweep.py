#!/usr/bin/env python
"""scripts/phase_beta_re_extract_sweep.py — Phase β L387 v54 patch re-extraction.

Re-extracts the n=25 snapshot ensemble from each target system's restrained
DCD using the L387 v54 patched ``utils.extract_snapshots.save_snapshots``
(Layer 0 ``_ensure_intrachain_peptide_bonds`` + Layer 1 ``image_molecules``
``make_whole=True`` + Layer 2 ``_unwrap_binder_manual``).

This is a CPU-only **re-extraction** — no MD re-run, no MM-GBSA. The
sweep writes to a NEW directory ``snapshots_n25_postl387_patch/`` per
system, leaving the original ``snapshots_n25/`` untouched (read-only).

Frame selection policy:
    - If ``snapshots_n25/`` already contains ``..._snap*_f<idx>.pdb``
      filenames, the frame indices are reused (preserves the original
      K-Means selection so the post-patch ensemble is statistically
      comparable to the pre-patch ensemble).
    - Otherwise, ``cluster_and_select(traj, 25)`` is invoked fresh.

Per-system errors are caught + logged; the sweep continues.
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import os
import re
import sys
import time
import warnings
from pathlib import Path
from typing import List, Optional, Tuple

# Silence mdtraj warnings before import.
warnings.filterwarnings("ignore")
logging.getLogger("mdtraj").setLevel(logging.ERROR)
logging.getLogger("openmm").setLevel(logging.ERROR)

# Make ``utils`` importable.
_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO))
sys.path.insert(0, str(_REPO / "utils"))

import mdtraj as md  # noqa: E402

try:
    md.set_logger_level(logging.ERROR)  # type: ignore[attr-defined]
except Exception:
    pass

from utils.extract_snapshots import (  # noqa: E402
    cluster_and_select,
    save_snapshots,
)


OUTPUTS_ROOT = _REPO / "outputs"
NEW_SUBDIR = "snapshots_n25_postl387_patch"
ORIG_SUBDIR = "snapshots_n25"
N_SNAPSHOTS = 25


# ---------------------------------------------------------------------------
# Target list — explicit (system_name, seed) tuples.
# ---------------------------------------------------------------------------

TARGETS: List[Tuple[str, str]] = [
    # Phase α primary (regression preservation) — 1EBP_MTR13 5 seeds
    ("1EBP_MTR13", "s7"),
    ("1EBP_MTR13", "s19"),
    ("1EBP_MTR13", "s23"),
    ("1EBP_MTR13", "s42"),
    ("1EBP_MTR13", "s101"),

    # Phase α primary — 2QKI_Cp4 8 completed seeds (skip s127, s19, s199)
    ("2QKI_Cp4", "s7"),
    ("2QKI_Cp4", "s19_reseed55"),
    ("2QKI_Cp4", "s23"),
    ("2QKI_Cp4", "s42"),
    ("2QKI_Cp4", "s83"),
    ("2QKI_Cp4", "s101"),
    ("2QKI_Cp4", "s163"),
    ("2QKI_Cp4", "s251"),

    # Phase α primary — 2QKH_MTR25 reseed (3 seeds)
    ("2QKH_MTR25", "s7"),
    ("2QKH_MTR25", "s19"),
    ("2QKH_MTR25", "s42"),

    # Phase α primary — 7TL8_MTR6 reseed (2 seeds; s7 already done in #81)
    ("7TL8_MTR6", "s19"),
    ("7TL8_MTR6", "s42"),

    # Phase β BROKEN — 3IOL_NML20 (4 seeds; mechanism (c) box-padding)
    ("3IOL_NML20", "s19"),
    ("3IOL_NML20", "s23"),
    ("3IOL_NML20", "s42"),
    ("3IOL_NML20", "s101"),

    # AMBIGUOUS — 1YCR_NML22 (5 seeds)
    ("1YCR_NML22", "s7"),
    ("1YCR_NML22", "s19"),
    ("1YCR_NML22", "s23"),
    ("1YCR_NML22", "s42"),
    ("1YCR_NML22", "s101"),

    # Already done in #81 — sweep verifies but does not re-extract.
    ("7TL8_MTR6", "s7"),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_FRAME_RE = re.compile(r"_snap\d+_f(\d+)\.pdb$")


def _system_dir(system: str, seed: str) -> Path:
    return OUTPUTS_ROOT / f"{system}_calib_{seed}"


def _find_dcd(sysdir: Path) -> Optional[Path]:
    for c in sysdir.glob("mdresult/*_restrained.dcd"):
        return c
    for c in sysdir.glob("mdresult/*.dcd"):
        return c
    return None


def _find_topology(sysdir: Path) -> Optional[Path]:
    """Pick the topology PDB the original extraction would have used.

    Preference order: mdresult/*_final.pdb → _md_input/*_renum.pdb →
    _md_input/*.pdb → mdresult/*.pdb.
    """
    candidates: List[Path] = []
    for pat in ("mdresult/*_final.pdb",):
        candidates.extend(sysdir.glob(pat))
    for pat in ("_md_input/*_renum.pdb", "_md_input/*.pdb"):
        candidates.extend(sysdir.glob(pat))
    for pat in ("mdresult/*.pdb",):
        candidates.extend(sysdir.glob(pat))
    seen = set()
    for c in candidates:
        if c in seen:
            continue
        seen.add(c)
        try:
            md.load(str(c))
            return c
        except Exception:
            continue
    return None


def _existing_basename_and_frames(orig_dir: Path) -> Tuple[Optional[str], List[int]]:
    """Inspect the original ``snapshots_n25/`` and recover (basename, frames).

    Basename is the prefix shared by ``<basename>_snap01_fXX.pdb`` filenames
    so the re-extraction emits an identical filename pattern.

    Returns (None, []) if the directory is missing or has no snapshots.
    """
    if not orig_dir.is_dir():
        return None, []
    pdbs = sorted(orig_dir.glob("*_snap*_f*.pdb"))
    if not pdbs:
        return None, []
    frames: List[int] = []
    basenames: List[str] = []
    for p in pdbs:
        m = _FRAME_RE.search(p.name)
        if not m:
            continue
        frames.append(int(m.group(1)))
        # Strip "_snapNN_fMMM.pdb" tail to get basename
        basenames.append(p.name[: m.start()])
    if not frames:
        return None, []
    # Basename should be unique per directory; take majority value.
    basename = max(set(basenames), key=basenames.count)
    return basename, sorted(set(frames))


def _verify_existing_patched(patched_dir: Path) -> bool:
    """Check the existing patched dir is fully populated (25 PDBs)."""
    if not patched_dir.is_dir():
        return False
    n = len(list(patched_dir.glob("*_snap*_f*.pdb")))
    return n >= N_SNAPSHOTS


# ---------------------------------------------------------------------------
# Main per-system worker
# ---------------------------------------------------------------------------

def reextract_one(system: str, seed: str) -> dict:
    """Re-extract n=25 snapshots for one (system, seed). Returns status dict."""
    record = {
        "system": system,
        "seed": seed,
        "tag": f"{system}_calib_{seed}",
        "status": "PENDING",
        "frames": [],
        "n_saved": 0,
        "skip_reason": None,
        "elapsed_s": None,
    }
    t0 = time.time()
    sysdir = _system_dir(system, seed)
    if not sysdir.is_dir():
        record["status"] = "MISSING_DIR"
        record["skip_reason"] = f"no dir: {sysdir}"
        return record

    patched_dir = sysdir / NEW_SUBDIR

    # Already-done short-circuit (only for 7TL8_MTR6_s7 in this sweep,
    # but also generic for any system with a complete patched dir).
    if _verify_existing_patched(patched_dir):
        n = len(list(patched_dir.glob("*_snap*_f*.pdb")))
        record["status"] = "ALREADY_DONE"
        record["n_saved"] = n
        record["elapsed_s"] = round(time.time() - t0, 2)
        return record

    dcd = _find_dcd(sysdir)
    if dcd is None:
        record["status"] = "SKIPPED"
        record["skip_reason"] = "no DCD"
        record["elapsed_s"] = round(time.time() - t0, 2)
        return record

    top = _find_topology(sysdir)
    if top is None:
        record["status"] = "SKIPPED"
        record["skip_reason"] = "no topology"
        record["elapsed_s"] = round(time.time() - t0, 2)
        return record

    orig_dir = sysdir / ORIG_SUBDIR
    basename, frames = _existing_basename_and_frames(orig_dir)

    # Load trajectory
    try:
        traj = md.load(str(dcd), top=str(top))
    except Exception as exc:
        record["status"] = "LOAD_FAILED"
        record["skip_reason"] = f"{type(exc).__name__}: {exc}"
        record["elapsed_s"] = round(time.time() - t0, 2)
        return record

    # If existing snapshots reference frames > traj.n_frames, refuse them.
    if frames and max(frames) >= traj.n_frames:
        # The original ensemble was extracted from a *trimmed* DCD (e.g.
        # 2QKI_Cp4_s7 uses the trimmed sub-trajectory). Fall back to fresh
        # K-Means on the current DCD to avoid IndexError.
        frames = []
        basename = None

    if not frames:
        # Fresh K-Means selection
        try:
            frames = cluster_and_select(traj, N_SNAPSHOTS, binder_chain="B")
        except Exception as exc:
            record["status"] = "SELECT_FAILED"
            record["skip_reason"] = f"{type(exc).__name__}: {exc}"
            record["elapsed_s"] = round(time.time() - t0, 2)
            return record

    if basename is None:
        # Standard convention used by extract_snapshots.main(): topology stem.
        basename = top.stem.replace("_final", "")

    # If the original snapshot count is < N_SNAPSHOTS, trim/pad.
    frames = sorted(set(int(f) for f in frames))[:N_SNAPSHOTS]

    record["frames"] = frames
    record["basename"] = basename

    # Run the patched save_snapshots
    patched_dir.mkdir(parents=True, exist_ok=True)
    try:
        saved = save_snapshots(traj, frames, str(patched_dir), basename)
    except Exception as exc:
        record["status"] = "SAVE_FAILED"
        record["skip_reason"] = f"{type(exc).__name__}: {exc}"
        record["elapsed_s"] = round(time.time() - t0, 2)
        return record

    record["status"] = "OK"
    record["n_saved"] = len(saved)
    record["elapsed_s"] = round(time.time() - t0, 2)
    return record


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Phase β re-extraction sweep using L387 v54 patched save_snapshots.",
    )
    parser.add_argument(
        "--out-json", default="/tmp/phase_beta_re_extract_status_20260429.json",
        help="Path to write per-system status JSON.",
    )
    parser.add_argument(
        "--filter", default=None,
        help="Optional substring filter on '<system>_calib_<seed>' tag.",
    )
    args = parser.parse_args(argv)

    targets = TARGETS
    if args.filter:
        targets = [t for t in targets if args.filter in f"{t[0]}_calib_{t[1]}"]

    print(f"Phase β re-extraction sweep — {len(targets)} systems")
    print(f"L387 v54 patched extract_snapshots.py (Layer 0+1+2)")
    print("=" * 72)

    records = []
    t_sweep = time.time()
    for i, (system, seed) in enumerate(targets, 1):
        tag = f"{system}_calib_{seed}"
        print(f"\n[{i}/{len(targets)}] {tag}")
        try:
            rec = reextract_one(system, seed)
        except Exception as exc:  # last-resort isolation
            rec = {
                "system": system, "seed": seed, "tag": tag,
                "status": "UNHANDLED",
                "skip_reason": f"{type(exc).__name__}: {exc}",
                "n_saved": 0, "elapsed_s": None,
            }
        status_short = rec["status"]
        if rec.get("skip_reason"):
            print(f"   → {status_short}: {rec['skip_reason']}")
        else:
            print(f"   → {status_short} (n_saved={rec.get('n_saved', 0)}, "
                  f"elapsed={rec.get('elapsed_s')}s)")
        records.append(rec)

    elapsed = round(time.time() - t_sweep, 1)
    n_ok = sum(1 for r in records if r["status"] == "OK")
    n_already = sum(1 for r in records if r["status"] == "ALREADY_DONE")
    n_skip = len(records) - n_ok - n_already
    print("\n" + "=" * 72)
    print(f"Sweep complete: OK={n_ok}, ALREADY_DONE={n_already}, "
          f"FAILED/SKIPPED={n_skip}, total_elapsed={elapsed}s")

    # Persist status
    Path(args.out_json).write_text(json.dumps({
        "schema": "phase_beta_re_extract/0.1",
        "n_targets": len(records),
        "n_ok": n_ok,
        "n_already_done": n_already,
        "n_skipped": n_skip,
        "elapsed_s": elapsed,
        "records": records,
    }, indent=2))
    print(f"Status JSON: {args.out_json}")

    return 0 if n_skip == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
