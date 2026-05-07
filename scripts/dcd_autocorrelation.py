#!/usr/bin/env python3
"""Compute autocorrelation time τ of DCD trajectories for Gate 1 / R9 mitigation.

For each outputs/*_calib_s*/mdresult/*_restrained.dcd:
  1. Load trajectory (mdtraj)
  2. Compute RMSD time series (backbone atoms vs first frame)
  3. Autocorrelation function (normalized) of that time series
  4. τ = first lag where ACF drops below 1/e (~0.368)
  5. Compare τ to snap spacing for n=10/25/50
     - ratio = spacing / τ
     - ratio > 1: snaps independent (good)
     - ratio < 1: snaps correlated (bad — effective N < nominal N)

Output: outputs/analysis/dcd_autocorrelation_20260422.json
CPU-only (numpy + mdtraj), parallel-safe with GPU MM-GBSA run.
"""
from __future__ import annotations

import glob
import json
import os
import sys
import time
from typing import Dict, List, Optional

import numpy as np

PROJ = "/home/san/UPDD_proj"
DCD_PATTERN = os.path.join(PROJ, "outputs", "*_calib_s*", "mdresult", "*_restrained.dcd")
OUT_PATH = os.path.join(PROJ, "outputs", "analysis", "dcd_autocorrelation_20260422.json")


def find_top_for_dcd(dcd_path: str) -> Optional[str]:
    """Locate a corresponding topology file (prmtop or pdb) for the DCD."""
    mdresult = os.path.dirname(dcd_path)
    # Try common candidates in preference order
    candidates: List[str] = []
    candidates.extend(sorted(glob.glob(os.path.join(mdresult, "*.prmtop"))))
    candidates.extend(sorted(glob.glob(os.path.join(mdresult, "*.parm7"))))
    # Parent run dir
    run_dir = os.path.dirname(mdresult)
    candidates.extend(sorted(glob.glob(os.path.join(run_dir, "*.prmtop"))))
    candidates.extend(sorted(glob.glob(os.path.join(run_dir, "*.parm7"))))
    # Fallback: first-frame PDB if available
    candidates.extend(sorted(glob.glob(os.path.join(mdresult, "*_init.pdb"))))
    candidates.extend(sorted(glob.glob(os.path.join(mdresult, "*.pdb"))))
    candidates.extend(sorted(glob.glob(os.path.join(run_dir, "*.pdb"))))
    for c in candidates:
        if os.path.isfile(c):
            return c
    return None


def compute_rmsd_series(dcd_path: str, top_path: str) -> np.ndarray:
    import mdtraj as md
    traj = md.load(dcd_path, top=top_path)
    # Select backbone atoms (C, CA, N) — robust selection that works for peptide+protein
    top = traj.topology
    bb_idx = top.select("backbone")
    if bb_idx.size < 4:
        bb_idx = top.select("name CA")
    if bb_idx.size < 4:
        bb_idx = np.arange(min(traj.n_atoms, 100))
    ref = traj.slice(0)
    rmsd = md.rmsd(traj, ref, atom_indices=bb_idx)
    return np.asarray(rmsd, dtype=np.float64)


def autocorr_tau(series: np.ndarray) -> Dict:
    """Compute normalized autocorrelation function and return τ (in frames).

    τ = first lag where ACF drops below 1/e (~0.3679).
    If ACF never crosses threshold within half the trajectory, return -1 and
    flag as inconclusive (trajectory likely shorter than τ).
    """
    n = series.size
    if n < 8:
        return {"tau_frames": -1, "n_frames": n, "flag": "too_short"}

    x = series - series.mean()
    var = np.dot(x, x) / n
    if var <= 0:
        return {"tau_frames": 0, "n_frames": n, "flag": "zero_variance"}

    # Normalized ACF via FFT (fast, standard)
    max_lag = n // 2
    # Use np.correlate full then slice, or FFT. FFT is faster for long series.
    fft_size = 1
    while fft_size < 2 * n:
        fft_size *= 2
    f = np.fft.rfft(x, n=fft_size)
    acf_full = np.fft.irfft(f * np.conj(f), n=fft_size)[:max_lag]
    # Normalize: acf(0) should be var * n
    acf = acf_full / (var * n)

    threshold = 1.0 / np.e  # ~0.3679
    # Find first index where acf drops below threshold
    below = np.where(acf < threshold)[0]
    if below.size == 0:
        return {"tau_frames": -1, "n_frames": n, "max_lag": int(max_lag),
                "acf_at_max_lag": float(acf[-1]), "flag": "no_decorrelation_within_half"}
    tau = int(below[0])
    return {"tau_frames": tau, "n_frames": n, "max_lag": int(max_lag),
            "acf_at_max_lag": float(acf[-1]), "flag": "ok"}


def analyze_one(dcd_path: str) -> Dict:
    run_dir = os.path.dirname(os.path.dirname(dcd_path))
    tag = os.path.basename(run_dir)
    top_path = find_top_for_dcd(dcd_path)
    if top_path is None:
        return {"tag": tag, "error": "no_topology_found"}
    t0 = time.time()
    try:
        series = compute_rmsd_series(dcd_path, top_path)
    except Exception as e:
        return {"tag": tag, "error": f"rmsd_failed: {type(e).__name__}: {e}",
                "top": top_path}
    stats = autocorr_tau(series)
    n_frames = int(series.size)
    elapsed = time.time() - t0

    # Snap spacing for n=10/25/50
    adequacy = {}
    for n_snap in (10, 25, 50):
        spacing = n_frames / n_snap if n_snap > 0 else 0
        tau = stats.get("tau_frames", -1)
        if tau is None or tau < 0:
            ratio = None
            flag = "tau_indeterminate"
        elif tau == 0:
            ratio = float("inf")
            flag = "ok"
        else:
            ratio = float(spacing / tau)
            flag = "ok" if ratio >= 1.0 else "correlated"
        adequacy[f"n{n_snap}"] = {
            "spacing_frames": round(spacing, 2),
            "ratio_spacing_over_tau": None if ratio is None else round(ratio, 3),
            "flag": flag,
        }

    return {
        "tag": tag,
        "dcd": dcd_path,
        "top": top_path,
        "n_frames": n_frames,
        "rmsd_mean_nm": float(series.mean()),
        "rmsd_std_nm": float(series.std()),
        "autocorr": stats,
        "adequacy": adequacy,
        "elapsed_s": round(elapsed, 2),
    }


def main() -> int:
    dcds = sorted(glob.glob(DCD_PATTERN))
    print(f"[autocorr] found {len(dcds)} DCD trajectories")
    results: List[Dict] = []
    t_overall = time.time()
    for idx, dcd in enumerate(dcds, 1):
        r = analyze_one(dcd)
        tag = r.get("tag", "?")
        if "error" in r:
            print(f"[{idx}/{len(dcds)}] ERROR {tag}: {r['error']}")
        else:
            ac = r["autocorr"]
            print(f"[{idx}/{len(dcds)}] {tag} n={r['n_frames']} "
                  f"tau={ac.get('tau_frames')} "
                  f"n25_ratio={r['adequacy']['n25']['ratio_spacing_over_tau']} "
                  f"n50_ratio={r['adequacy']['n50']['ratio_spacing_over_tau']} "
                  f"elapsed={r['elapsed_s']}s")
        results.append(r)

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w") as fh:
        json.dump({
            "n_trajectories": len(dcds),
            "threshold": "1/e",
            "units": {"tau": "frames", "spacing": "frames"},
            "results": results,
            "total_elapsed_s": round(time.time() - t_overall, 1),
        }, fh, indent=2)
    print(f"[autocorr] output: {OUT_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
