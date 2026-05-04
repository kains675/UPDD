#!/usr/bin/env python3
"""DEPRECATED (2026-04-22): superseded by utils/mmgbsa_outlier_screen.py
(production helper, MAD-any-of v1 with sanitization_policy schema field).

Use scripts/sanitize_summaries.py for batch sanitization or import
mmgbsa_outlier_screen.sanitize_mmgbsa_summary() directly.

Original functionality retained here for reference + Quick Win 3 provenance
(R-7: no deletion of historical artefacts).

Screen MM-GBSA per-snap outliers via MAD-based robust z-scores.

Rationale (Hampel FR, 1974, J. Am. Stat. Assoc. 69(346), 383-393):
    Classical z-score (x - mean)/stdev is non-robust: a single pathological
    outlier inflates stdev and masks its own deviation. The median absolute
    deviation (MAD = median(|x_i - median(x)|)) provides a breakdown-point
    50% robust scale estimator. Scaling MAD / 0.6745 makes it a consistent
    estimator of Gaussian stdev (since 0.6745 = Phi^{-1}(0.75)). The 3-sigma
    threshold corresponds to the standard 0.3% two-tailed Gaussian tail.

Motivating case: 2QKI_Cp4_s7 snap08 has e_ligand ~= +26415 kcal/mol whereas
its neighbours span -550 to -604 kcal/mol. This is a bound-geometry strain
artefact of 3-trajectory MM-GBSA split_complex: the isolated-ligand energy
is evaluated at a strained bound conformation. We flag such snaps and
recompute a sanitised mean_dg over the retained set.
"""

import glob
import json
import math
import os
import statistics
from typing import Any, Dict, List, Optional, Tuple

GLOB_PATTERN = "/home/san/UPDD_proj/outputs/*_calib_s*/mmgbsa_results/mmgbsa_summary.json"
OUT_JSON = "/home/san/UPDD_proj/outputs/analysis/outlier_screen_20260422.json"
COMPONENTS = ("e_complex_kcal", "e_receptor_kcal", "e_ligand_kcal", "delta_g_kcal")
Z_THRESHOLD = 3.0
MAD_CONSISTENCY = 0.6745  # median(|x - med|) / 0.6745 ~ stdev for Gaussian


def robust_scale(xs: List[float]) -> Tuple[float, float, bool]:
    """Return (median, stdev-equivalent scale, mad_degenerate)."""
    med = statistics.median(xs)
    mad = statistics.median([abs(x - med) for x in xs])
    if mad == 0.0:
        scale = statistics.pstdev(xs) if len(xs) > 1 else 0.0
        return med, scale, True
    return med, mad / MAD_CONSISTENCY, False


def z_scores(xs: List[float]) -> Tuple[List[float], bool]:
    med, scale, degen = robust_scale(xs)
    if scale == 0.0:
        return [0.0] * len(xs), degen
    return [(x - med) / scale for x in xs], degen


def process_run(json_path: str) -> Optional[Dict[str, Any]]:
    with open(json_path) as f:
        data = json.load(f)
    results = data.get("results", [])
    if not results:
        return None
    tag = os.path.basename(os.path.dirname(os.path.dirname(json_path)))
    comp_z: Dict[str, List[float]] = {}
    comp_degen: Dict[str, bool] = {}
    for key in COMPONENTS:
        xs = [float(r[key]) for r in results]
        zs, degen = z_scores(xs)
        comp_z[key] = zs
        comp_degen[key] = degen
    flagged: List[Dict[str, Any]] = []
    flag_mask = [False] * len(results)
    for i, r in enumerate(results):
        reasons = [k for k in ("e_complex_kcal", "e_ligand_kcal", "delta_g_kcal")
                   if abs(comp_z[k][i]) > Z_THRESHOLD]
        if reasons:
            flag_mask[i] = True
            flagged.append({
                "run_tag": tag,
                "snap_name": r["snapshot"],
                "z_e_complex": round(comp_z["e_complex_kcal"][i], 3),
                "z_e_receptor": round(comp_z["e_receptor_kcal"][i], 3),
                "z_e_ligand": round(comp_z["e_ligand_kcal"][i], 3),
                "z_delta_g": round(comp_z["delta_g_kcal"][i], 3),
                "e_ligand_kcal": r["e_ligand_kcal"],
                "delta_g_kcal": r["delta_g_kcal"],
                "flag_reason": ",".join(reasons),
                "mad_degenerate": {k: comp_degen[k] for k in reasons},
            })
    dgs = [float(r["delta_g_kcal"]) for r in results]
    dgs_clean = [dg for dg, bad in zip(dgs, flag_mask) if not bad]
    mean_dg_orig = statistics.fmean(dgs)
    mean_dg_clean = statistics.fmean(dgs_clean) if dgs_clean else math.nan
    sigma_dg_clean = statistics.pstdev(dgs_clean) if len(dgs_clean) > 1 else 0.0
    return {
        "tag": tag,
        "n_snap_original": len(results),
        "n_snap_retained": len(dgs_clean),
        "mean_dg_original": round(mean_dg_orig, 4),
        "mean_dg_clean": round(mean_dg_clean, 4) if not math.isnan(mean_dg_clean) else None,
        "sigma_dg_clean": round(sigma_dg_clean, 4),
        "delta_mean": round(mean_dg_clean - mean_dg_orig, 4) if not math.isnan(mean_dg_clean) else None,
        "flagged": flagged,
    }


def main() -> None:
    paths = sorted(glob.glob(GLOB_PATTERN))
    per_run = [r for r in (process_run(p) for p in paths) if r is not None]
    all_flagged = [f for r in per_run for f in r["flagged"]]
    all_flagged_sorted = sorted(all_flagged, key=lambda f: abs(f["z_e_ligand"]), reverse=True)
    total_snaps = sum(r["n_snap_original"] for r in per_run)
    shifts = [abs(r["delta_mean"]) for r in per_run if r["delta_mean"] is not None]
    avg_shift = statistics.fmean(shifts) if shifts else 0.0
    print(f"Runs scanned: {len(per_run)}  total snaps: {total_snaps}  flagged: {len(all_flagged)}")
    print(f"Avg |delta_mean_dg| per run: {avg_shift:.4f} kcal/mol")
    print("\nFlagged snaps (sorted by |z_e_ligand|):")
    print(f"{'run_tag':<28} {'snap':<48} {'z_eL':>8} {'z_eC':>8} {'z_dG':>8} {'e_lig_kcal':>14}  reason")
    for f in all_flagged_sorted:
        print(f"{f['run_tag']:<28} {f['snap_name']:<48} {f['z_e_ligand']:>8.2f} "
              f"{f['z_e_complex']:>8.2f} {f['z_delta_g']:>8.2f} {f['e_ligand_kcal']:>14.3f}  {f['flag_reason']}")
    out = {
        "threshold_z": Z_THRESHOLD,
        "mad_consistency_constant": MAD_CONSISTENCY,
        "n_runs": len(per_run),
        "n_snaps_total": total_snaps,
        "n_flagged": len(all_flagged),
        "avg_abs_delta_mean_dg": round(avg_shift, 4),
        "per_run": per_run,
        "flagged_sorted_by_z_e_ligand": all_flagged_sorted,
    }
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\nJSON report -> {OUT_JSON}")


if __name__ == "__main__":
    main()
