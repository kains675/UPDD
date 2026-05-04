#!/usr/bin/env python
"""scripts/phase_beta_repbsa_v2_reextract.py — Phase β re-PBSA v2 re-extraction.

Re-extracts the n=25 snapshot ensemble for the 10 systems that previously
failed in #84 due to CONECT wraparound, using the CONECT v55 patched
``utils.extract_snapshots.save_snapshots`` (atom-index Part A + hex-extended
Part B re-serialization at >99,999 atoms).

Output dir: ``snapshots_n25_postl387_patch_v2/`` per system. The previous
``snapshots_n25_postl387_patch/`` (pre-v55 broken-CONECT) is preserved
read-only for forensic.

Frame selection policy (preserves K-Means selection from #83):
    1. Read frame indices from ``snapshots_n25_postl387_patch/`` filenames
       (``..._snap*_f<idx>.pdb``).
    2. Fallback: read from ``snapshots_n25/``.
    3. Fallback: fresh ``cluster_and_select(traj, 25)``.
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import sys
import time
import warnings
from pathlib import Path
from typing import List, Optional, Tuple

warnings.filterwarnings("ignore")
logging.getLogger("mdtraj").setLevel(logging.ERROR)
logging.getLogger("openmm").setLevel(logging.ERROR)

_REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO))
sys.path.insert(0, str(_REPO / "utils"))

import mdtraj as md  # noqa: E402

try:
    md.set_logger_level(logging.ERROR)
except Exception:
    pass

from utils.extract_snapshots import cluster_and_select, save_snapshots  # noqa: E402

OUTPUTS_ROOT = _REPO / "outputs"
NEW_SUBDIR = "snapshots_n25_postl387_patch_v2"
PRIOR_SUBDIRS = ("snapshots_n25_postl387_patch", "snapshots_n25")
N_SNAPSHOTS = 25

# Phase β re-PBSA v2 — 10 target systems (CONECT wraparound failures from #84).
TARGETS: List[Tuple[str, str]] = [
    ("2QKI_Cp4", "s7"),
    ("2QKI_Cp4", "s19_reseed55"),
    ("2QKI_Cp4", "s23"),
    ("2QKI_Cp4", "s42"),
    ("2QKI_Cp4", "s83"),
    ("2QKI_Cp4", "s101"),
    ("2QKI_Cp4", "s163"),
    ("2QKI_Cp4", "s251"),
    ("7TL8_MTR6", "s19"),
    ("7TL8_MTR6", "s42"),
]

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


def _existing_basename_and_frames(sysdir: Path) -> Tuple[Optional[str], List[int]]:
    """Recover (basename, frames) from prior snapshot dirs (preserve K-Means)."""
    for sub in PRIOR_SUBDIRS:
        d = sysdir / sub
        if not d.is_dir():
            continue
        pdbs = sorted(d.glob("*_snap*_f*.pdb"))
        if not pdbs:
            continue
        frames: List[int] = []
        basenames: List[str] = []
        for p in pdbs:
            m = _FRAME_RE.search(p.name)
            if not m:
                continue
            frames.append(int(m.group(1)))
            basenames.append(p.name[: m.start()])
        if frames:
            basename = max(set(basenames), key=basenames.count)
            return basename, sorted(set(frames))
    return None, []


def _verify_existing(patched_dir: Path) -> bool:
    if not patched_dir.is_dir():
        return False
    n = len(list(patched_dir.glob("*_snap*_f*.pdb")))
    return n >= N_SNAPSHOTS


def reextract_one(system: str, seed: str) -> dict:
    record = {
        "system": system, "seed": seed, "tag": f"{system}_calib_{seed}",
        "status": "PENDING", "frames": [], "n_saved": 0,
        "skip_reason": None, "elapsed_s": None,
        "n_conect_per_pdb": [],
    }
    t0 = time.time()
    sysdir = _system_dir(system, seed)
    if not sysdir.is_dir():
        record["status"] = "MISSING_DIR"
        record["skip_reason"] = f"no dir: {sysdir}"
        return record

    patched_dir = sysdir / NEW_SUBDIR

    if _verify_existing(patched_dir):
        n = len(list(patched_dir.glob("*_snap*_f*.pdb")))
        record["status"] = "ALREADY_DONE"
        record["n_saved"] = n
        record["elapsed_s"] = round(time.time() - t0, 2)
        # Audit CONECT counts even on already-done
        for p in sorted(patched_dir.glob("*_snap*_f*.pdb")):
            n_c = sum(1 for line in open(p) if line.startswith("CONECT"))
            record["n_conect_per_pdb"].append(n_c)
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

    basename, frames = _existing_basename_and_frames(sysdir)

    try:
        traj = md.load(str(dcd), top=str(top))
    except Exception as exc:
        record["status"] = "LOAD_FAILED"
        record["skip_reason"] = f"{type(exc).__name__}: {exc}"
        record["elapsed_s"] = round(time.time() - t0, 2)
        return record

    if frames and max(frames) >= traj.n_frames:
        frames = []
        basename = None

    if not frames:
        try:
            frames = cluster_and_select(traj, N_SNAPSHOTS, binder_chain="B")
        except Exception as exc:
            record["status"] = "SELECT_FAILED"
            record["skip_reason"] = f"{type(exc).__name__}: {exc}"
            record["elapsed_s"] = round(time.time() - t0, 2)
            return record

    if basename is None:
        basename = top.stem.replace("_final", "")

    frames = sorted(set(int(f) for f in frames))[:N_SNAPSHOTS]
    record["frames"] = frames
    record["basename"] = basename

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

    # Audit CONECT counts
    for p in sorted(patched_dir.glob("*_snap*_f*.pdb")):
        n_c = sum(1 for line in open(p) if line.startswith("CONECT"))
        record["n_conect_per_pdb"].append(n_c)

    return record


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Phase β re-PBSA v2 re-extraction (CONECT v55 patched).",
    )
    parser.add_argument(
        "--out-json",
        default="/tmp/phase_beta_repbsa_v2_reextract_status.json",
    )
    parser.add_argument("--filter", default=None)
    args = parser.parse_args(argv)

    targets = TARGETS
    if args.filter:
        targets = [t for t in targets if args.filter in f"{t[0]}_calib_{t[1]}"]

    print(f"Phase β re-PBSA v2 re-extraction — {len(targets)} systems")
    print(f"Source: extract_snapshots.py CONECT v55 patched")
    print(f"Output: outputs/<sys>/{NEW_SUBDIR}/")
    print("=" * 72)

    records = []
    t_sweep = time.time()
    for i, (system, seed) in enumerate(targets, 1):
        tag = f"{system}_calib_{seed}"
        print(f"\n[{i}/{len(targets)}] {tag}")
        try:
            rec = reextract_one(system, seed)
        except Exception as exc:
            rec = {
                "system": system, "seed": seed, "tag": tag,
                "status": "UNHANDLED",
                "skip_reason": f"{type(exc).__name__}: {exc}",
                "n_saved": 0, "elapsed_s": None,
                "n_conect_per_pdb": [],
            }
        n_conect = rec.get("n_conect_per_pdb", [])
        if n_conect:
            mn, mx = min(n_conect), max(n_conect)
            n_zero = sum(1 for x in n_conect if x == 0)
            cstat = f"CONECT min={mn} max={mx} n_zero={n_zero}"
        else:
            cstat = "CONECT —"
        if rec.get("skip_reason"):
            print(f"   → {rec['status']}: {rec['skip_reason']}")
        else:
            print(f"   → {rec['status']} (n_saved={rec.get('n_saved')}, "
                  f"elapsed={rec.get('elapsed_s')}s, {cstat})")
        records.append(rec)

    elapsed = round(time.time() - t_sweep, 1)
    n_ok = sum(1 for r in records if r["status"] == "OK")
    n_already = sum(1 for r in records if r["status"] == "ALREADY_DONE")
    n_fail = len(records) - n_ok - n_already
    print("\n" + "=" * 72)
    print(f"Sweep: OK={n_ok}, ALREADY_DONE={n_already}, "
          f"FAILED={n_fail}, total={elapsed}s")

    Path(args.out_json).write_text(json.dumps({
        "schema": "phase_beta_repbsa_v2_reextract/0.1",
        "n_targets": len(records),
        "n_ok": n_ok,
        "n_already_done": n_already,
        "n_failed": n_fail,
        "elapsed_s": elapsed,
        "records": records,
    }, indent=2))
    print(f"Status JSON: {args.out_json}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
