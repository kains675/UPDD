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
TARGETS_PHASE_BETA: List[Tuple[str, str]] = [
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

# T1 N→30 expansion — Phase 1 (2026-05-18) — Cp4 MD-ready seeds that need
# v2 snapshot extraction. 7 seeds already have mdresult/ but no snapshots
# under either v1 or v2 paths. Reextract under CONECT v55 (Cp4 atom count
# 123,083 ≥ 99,999 → Case B wraparound trigger, v2 mandatory).
TARGETS_T1_PHASE_1_CP4: List[Tuple[str, str]] = [
    ("2QKI_Cp4", "s127"),
    ("2QKI_Cp4", "s19"),
    ("2QKI_Cp4", "s199"),
    ("2QKI_Cp4", "s433"),
    ("2QKI_Cp4", "s523"),
    ("2QKI_Cp4", "s773"),
    ("2QKI_Cp4", "s991"),
]

# T1 N→30 expansion — Phase 1.5 (2026-05-18) — Cp4 fresh GAFF2 MD seeds.
# Option β: 4 new seeds picked from free integers (no current or archive
# collision). MD launched with UNPATCHED params (from
# outputs/_archive/pre_amb14_patch_20260427/) for cohort parity with
# baseline 11.
TARGETS_T1_PHASE_15_CP4: List[Tuple[str, str]] = [
    ("2QKI_Cp4", "s53"),
    ("2QKI_Cp4", "s89"),
    ("2QKI_Cp4", "s113"),
    ("2QKI_Cp4", "s173"),
]

# T1 N→30 expansion — Phase 1.5 (2026-05-18) — 2QKI_WT v2 reextract for
# ALL 10 WT seeds (all have completed DCDs per preflight inventory).
# Includes pilot 2 (s7, s23) for idempotent skip + 8 new (s19, s42, s83,
# s101, s127, s163, s199, s251) for v2 reextract + MMPBSA. Matched-integer
# pairing with Cp4 baseline cohort for branched ΔΔG aggregation.
TARGETS_T1_PHASE_15_WT_FULL: List[Tuple[str, str]] = [
    ("2QKI_WT", "s19"),
    ("2QKI_WT", "s42"),
    ("2QKI_WT", "s83"),
    ("2QKI_WT", "s101"),
    ("2QKI_WT", "s127"),
    ("2QKI_WT", "s163"),
    ("2QKI_WT", "s199"),
    ("2QKI_WT", "s251"),
]

# T1 N→30 expansion — Phase 1 (2026-05-18) — 2QKH_MTR25 MD-ready seeds.
# Case A confirmed (45,828 atoms < 99,999) → v1 and v2 are numerically
# equivalent for this system. For symmetry with the existing 15 reps
# under snapshots_n25_postl387_patch/, write to the SAME v1 path.
# The CONECT v55 patch is harmless at this atom count.
# NOTE: this script writes to snapshots_n25_postl387_patch_v2/ by default
# (NEW_SUBDIR); for these 3 seeds run with --out-subdir snapshots_n25_postl387_patch
# to write into the v1 path.
TARGETS_T1_PHASE_1_MTR25: List[Tuple[str, str]] = [
    ("2QKH_MTR25", "s479"),
    ("2QKH_MTR25", "s997"),
    ("2QKH_MTR25", "s1549"),
]

# Phase 1C — 2QKI_WT pilot under postl387_v2 (Cp4 reference branch).
# Two pilot seeds matched to Cp4 baseline (s7, s23) to validate the v2
# snapshot extractor on the WT MD trajectories before committing to full
# N=30 in Phase 3.
TARGETS_T1_PHASE_1_WT_PILOT: List[Tuple[str, str]] = [
    ("2QKI_WT", "s7"),
    ("2QKI_WT", "s23"),
]

TARGETS: List[Tuple[str, str]] = (
    TARGETS_PHASE_BETA
    + TARGETS_T1_PHASE_1_CP4
    + TARGETS_T1_PHASE_1_MTR25
    + TARGETS_T1_PHASE_1_WT_PILOT
    + TARGETS_T1_PHASE_15_CP4
    + TARGETS_T1_PHASE_15_WT_FULL
)

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


def reextract_one(system: str, seed: str, out_subdir: str = NEW_SUBDIR) -> dict:
    record = {
        "system": system, "seed": seed, "tag": f"{system}_calib_{seed}",
        "status": "PENDING", "frames": [], "n_saved": 0,
        "skip_reason": None, "elapsed_s": None,
        "n_conect_per_pdb": [],
        "out_subdir": out_subdir,
    }
    t0 = time.time()
    sysdir = _system_dir(system, seed)
    if not sysdir.is_dir():
        record["status"] = "MISSING_DIR"
        record["skip_reason"] = f"no dir: {sysdir}"
        return record

    patched_dir = sysdir / out_subdir

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
    parser.add_argument(
        "--out-subdir", default=NEW_SUBDIR,
        help=(
            "Subdir name under each outputs/<system>_calib_<seed>/ to write "
            f"snapshots into (default: {NEW_SUBDIR}). Use "
            "'snapshots_n25_postl387_patch' to write into the v1 path for "
            "Case-A systems (atom count < 99,999) that stay numerically "
            "equivalent under the v55 patch."
        ),
    )
    args = parser.parse_args(argv)

    targets = TARGETS
    if args.filter:
        targets = [t for t in targets if args.filter in f"{t[0]}_calib_{t[1]}"]

    print(f"Phase β re-PBSA v2 re-extraction — {len(targets)} systems")
    print(f"Source: extract_snapshots.py CONECT v55 patched")
    print(f"Output: outputs/<sys>/{args.out_subdir}/")
    print("=" * 72)

    records = []
    t_sweep = time.time()
    for i, (system, seed) in enumerate(targets, 1):
        tag = f"{system}_calib_{seed}"
        print(f"\n[{i}/{len(targets)}] {tag}")
        try:
            rec = reextract_one(system, seed, out_subdir=args.out_subdir)
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
