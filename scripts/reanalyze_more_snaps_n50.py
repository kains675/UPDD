#!/usr/bin/env python3
"""Re-extract snapshots (n=50) and re-run MM-GBSA on existing UPDD calibration trajectories.

Rationale
---------
Derived from `reanalyze_more_snaps.py` (n=25). Gate 2 (ΔΔG z>2) verdict on
2026-04-22 showed 0/5 pairs passing at n=25 → σ_wi dominates. SciVal §A.3
sampling expansion axis (orthogonal to bias-axis §A.2 Fix 4) calls for an
n=50 pass to test whether doubling snapshots reduces σ_wi toward Gate 2
threshold. Output dirs: `snapshots_n50/` and `mmgbsa_results_n50/` — parallel
to n=25, forensic preservation of all prior tiers (R-7).

This is the v0.6.7 3-traj protocol baseline (no Fix 4). Fix 4 is Coder's
parallel in-place task and MUST NOT be invoked in this run (see task brief).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import subprocess
import sys
import time
from typing import Dict, List, Optional, Tuple

PY = "/home/san/miniconda3/envs/md_simulation/bin/python"
PROJ = "/home/san/UPDD_proj"
DCD_PATTERN = os.path.join(PROJ, "outputs", "*_calib_s*", "mdresult", "*_restrained.dcd")
STATUS_OUT = os.path.join(PROJ, "outputs", "analysis", "reanalyze_more_snaps_n50_status_20260422.json")
N_SNAPS = 50


def derive_target_id(run_tag: str) -> str:
    """Run dir name '2QKI_WT_calib_s42' -> target_id '2QKI'."""
    return run_tag.split("_", 1)[0]


def discover_runs() -> List[Tuple[str, str]]:
    """Return list of (run_dir, tag) for each matching DCD trajectory."""
    dcds = sorted(glob.glob(DCD_PATTERN))
    runs: List[Tuple[str, str]] = []
    for dcd in dcds:
        run_dir = os.path.dirname(os.path.dirname(dcd))
        tag = os.path.basename(run_dir)
        runs.append((run_dir, tag))
    return runs


def build_commands(run_dir: str, tag: str, platform: str) -> Tuple[List[str], List[str], Dict[str, str]]:
    target_id = derive_target_id(tag)
    snap_dir = os.path.join(run_dir, "snapshots_n50")
    mmgbsa_dir = os.path.join(run_dir, "mmgbsa_results_n50")
    md_dir = os.path.join(run_dir, "mdresult")

    extract_cmd = [
        PY, os.path.join(PROJ, "utils", "extract_snapshots.py"),
        "--md_dir", md_dir,
        "--outputdir", snap_dir,
        "--n_snapshots", str(N_SNAPS),
        "--binder_chain", "B",
        "--target_id", target_id,
    ]
    mmgbsa_cmd = [
        PY, os.path.join(PROJ, "utils", "run_mmgbsa.py"),
        "--md_dir", snap_dir,
        "--outputdir", mmgbsa_dir,
        "--ncaa_elem", "none",
        "--receptor_chain", "A",
        "--binder_chain", "B",
        "--target_id", target_id,
    ]
    env_overlay = {"UPDD_MMGBSA_PLATFORM": platform}
    return extract_cmd, mmgbsa_cmd, env_overlay


def run_cmd(cmd: List[str], env_overlay: Dict[str, str], log_fp) -> int:
    env = os.environ.copy()
    env.update(env_overlay)
    log_fp.write("$ " + " ".join(cmd) + "\n")
    log_fp.flush()
    proc = subprocess.run(cmd, stdout=log_fp, stderr=subprocess.STDOUT, env=env)
    return proc.returncode


def load_summary(mmgbsa_dir: str) -> Tuple[Optional[float], Optional[float]]:
    summary_path = os.path.join(mmgbsa_dir, "mmgbsa_summary.json")
    if not os.path.isfile(summary_path):
        return None, None
    try:
        with open(summary_path) as fh:
            d = json.load(fh)
        mean_dg = d.get("mean_dg")
        ie = list(d.get("interaction_entropy_per_design", {}).values())
        sigma = ie[0].get("std_delta_g_kcal") if ie else None
        return mean_dg, sigma
    except Exception:
        return None, None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--platform", choices=["CPU", "CUDA"], default="CPU",
                        help="MM-GBSA platform (default CPU for GPU courtesy)")
    parser.add_argument("--max_runs", type=int, default=0,
                        help="Cap number of runs processed (0 = all)")
    parser.add_argument("--only", default="",
                        help="Substring filter on run tag")
    parser.add_argument("--dry_run", action="store_true",
                        help="Print commands without executing")
    args = parser.parse_args()

    runs = discover_runs()
    if args.only:
        runs = [r for r in runs if args.only in r[1]]
    if args.max_runs > 0:
        runs = runs[: args.max_runs]

    print(f"[orchestrator] detected {len(runs)} run(s) after filters (platform={args.platform}, dry_run={args.dry_run})")
    status: Dict[str, Dict] = {}
    elapsed_samples: List[float] = []

    for idx, (run_dir, tag) in enumerate(runs, 1):
        mmgbsa_dir = os.path.join(run_dir, "mmgbsa_results_n50")
        summary_path = os.path.join(mmgbsa_dir, "mmgbsa_summary.json")
        if os.path.isfile(summary_path):
            print(f"[{idx}/{len(runs)}] SKIP {tag} (mmgbsa_results_n50/mmgbsa_summary.json exists)")
            mean_dg, sigma = load_summary(mmgbsa_dir)
            status[tag] = {"state": "skipped", "elapsed_s": 0.0,
                           "mean_dg": mean_dg, "sigma_wi": sigma}
            continue

        extract_cmd, mmgbsa_cmd, env_overlay = build_commands(run_dir, tag, args.platform)
        if args.dry_run:
            print(f"[{idx}/{len(runs)}] DRY {tag}")
            print("  extract: " + " ".join(extract_cmd))
            print(f"  mmgbsa : UPDD_MMGBSA_PLATFORM={args.platform} " + " ".join(mmgbsa_cmd))
            status[tag] = {"state": "dry", "elapsed_s": 0.0, "mean_dg": None, "sigma_wi": None}
            continue

        log_path = f"/tmp/reanalyze_more_snaps_n50_{tag}.log"
        t0 = time.time()
        with open(log_path, "w") as log_fp:
            log_fp.write(f"# reanalyze_more_snaps_n50 tag={tag} platform={args.platform}\n")
            rc_e = run_cmd(extract_cmd, {}, log_fp)
            rc_m = run_cmd(mmgbsa_cmd, env_overlay, log_fp) if rc_e == 0 else -1
        elapsed = time.time() - t0
        state = "done" if (rc_e == 0 and rc_m == 0 and os.path.isfile(summary_path)) else "failed"
        mean_dg, sigma = load_summary(mmgbsa_dir) if state == "done" else (None, None)
        status[tag] = {"state": state, "elapsed_s": round(elapsed, 1),
                       "mean_dg": mean_dg, "sigma_wi": sigma, "log": log_path,
                       "rc_extract": rc_e, "rc_mmgbsa": rc_m}
        elapsed_samples.append(elapsed)
        print(f"[{idx}/{len(runs)}] {state.upper()} {tag} elapsed={elapsed:.1f}s "
              f"⟨ΔG⟩={mean_dg} σ_wi={sigma}")
        if idx == 2 and elapsed_samples:
            avg = sum(elapsed_samples) / len(elapsed_samples)
            remaining = len(runs) - idx
            eta_min = (avg * remaining) / 60.0
            print(f"[orchestrator] ETA estimate: avg {avg:.1f}s/run × {remaining} remaining ≈ {eta_min:.1f} min")

    os.makedirs(os.path.dirname(STATUS_OUT), exist_ok=True)
    with open(STATUS_OUT, "w") as fh:
        json.dump({"n_runs": len(runs), "platform": args.platform,
                   "n_snapshots": N_SNAPS, "results": status}, fh, indent=2)
    print(f"[orchestrator] status written: {STATUS_OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
