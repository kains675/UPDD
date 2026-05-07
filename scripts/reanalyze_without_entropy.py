#!/usr/bin/env python3
"""Entropy-free reanalysis of MM-GBSA summaries (R-11 ranking-only regime).

Purpose
-------
Produce an entropy-free ranking table by removing the Duan 2016 interaction-entropy
correction for runs whose per-design σ_wi exceeds the Jiao-Lin 2021 reliability
threshold (σ_wi < 3.6 kcal/mol). When IE averaging is unreliable
(σ_wi ≫ 3.6 kcal/mol), the exp(+βΔG) weighting is dominated by a single outlier
frame, producing a spurious -TΔS term that swamps the ΔΔG ranking signal.

Regime
------
R-11 ranking-only. Absolute ΔG values are not physically calibrated; only
inter-design ΔΔG ordering is reported.

Jiao-Lin 2021 rule
------------------
IE is usable iff (σ_wi ≤ 3.6 kcal/mol) AND (reliable flag from upstream IE
computation is True). When not usable, fall back to uncorrected mean ΔG
(ΔG_preferred = mean_dg), matching the standard MM-GBSA endpoint method.

Inputs
------
Globs outputs/*_calib_s*/mmgbsa_results/mmgbsa_summary.json

Outputs
-------
- stdout: pretty per-run table + per-target ΔΔG table
- outputs/analysis/ranking_entropy_free_20260422.json: structured dump
"""

import glob
import json
import os
import statistics
from typing import Any, Dict, List, Optional, Tuple

PROJECT = "/home/san/UPDD_proj"
GLOB_PATTERN = os.path.join(
    PROJECT, "outputs", "*_calib_s*", "mmgbsa_results", "mmgbsa_summary.json"
)
OUT_DIR = os.path.join(PROJECT, "outputs", "analysis")
OUT_JSON = os.path.join(OUT_DIR, "ranking_entropy_free_20260422.json")
SIGMA_THRESHOLD = 3.6  # Jiao-Lin 2021 kcal/mol


def parse_tag(run_dir: str) -> Tuple[str, Optional[int]]:
    """Split run dir like '2QKI_Cp4_calib_s7' into ('2QKI_Cp4_calib', 7)."""
    tag = os.path.basename(run_dir)
    if "_calib_s" not in tag:
        return tag, None
    left, seed_str = tag.rsplit("_calib_s", 1)
    try:
        seed = int(seed_str)
    except ValueError:
        return tag, None
    return left + "_calib", seed


def extract_run(summary_path: str) -> Dict[str, Any]:
    """Pull per-run fields; tolerate absent IE dicts."""
    with open(summary_path, "r") as fh:
        data = json.load(fh)
    run_dir = os.path.dirname(os.path.dirname(summary_path))
    tag = os.path.basename(run_dir)
    group, seed = parse_tag(run_dir)
    mean_dg = data.get("mean_dg")
    n_snap = data.get("n_calc")
    ie_map = data.get("interaction_entropy_per_design") or {}
    # IE dict has one entry per design; use the first (runs are single-design).
    sigma_wi = None
    minus_tds = None
    dg_corrected = None
    reliable_flag = False
    if ie_map:
        first_key = next(iter(ie_map))
        entry = ie_map[first_key] or {}
        sigma_wi = entry.get("std_delta_g_kcal")
        minus_tds = entry.get("minus_TdS_kcal")
        dg_corrected = entry.get("delta_g_corrected_kcal")
        reliable_flag = bool(entry.get("reliable", False))
    ie_usable = (
        sigma_wi is not None
        and sigma_wi <= SIGMA_THRESHOLD
        and reliable_flag
    )
    dg_preferred = dg_corrected if ie_usable else mean_dg
    return {
        "tag": tag,
        "group": group,
        "seed": seed,
        "n_snap": n_snap,
        "mean_dg": mean_dg,
        "sigma_wi": sigma_wi,
        "minus_TdS": minus_tds,
        "delta_g_corrected": dg_corrected,
        "reliable_flag": reliable_flag,
        "ie_usable": ie_usable,
        "delta_g_preferred": dg_preferred,
    }


def aggregate_group(runs: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Compute inter-replicate aggregates for one target_variant group."""
    dgs = [r["delta_g_preferred"] for r in runs if r["delta_g_preferred"] is not None]
    sigmas = [r["sigma_wi"] for r in runs if r["sigma_wi"] is not None]
    mean_dg = statistics.fmean(dgs) if dgs else None
    sigma_between = statistics.stdev(dgs) if len(dgs) >= 2 else 0.0
    sigma_within_avg = statistics.fmean(sigmas) if sigmas else None
    any_ie_usable = any(r["ie_usable"] for r in runs)
    all_ie_usable = all(r["ie_usable"] for r in runs) if runs else False
    return {
        "n_rep": len(runs),
        "mean_dg": mean_dg,
        "sigma_between": sigma_between,
        "sigma_within_avg": sigma_within_avg,
        "any_ie_usable": any_ie_usable,
        "all_ie_usable": all_ie_usable,
        "seeds": sorted(r["seed"] for r in runs if r["seed"] is not None),
    }


def split_target_variant(group: str) -> Tuple[str, str]:
    """'2QKI_Cp4_calib' -> ('2QKI', 'Cp4'); 'X_Y_Z_calib' -> ('X', 'Y_Z')."""
    core = group[: -len("_calib")] if group.endswith("_calib") else group
    parts = core.split("_", 1)
    if len(parts) != 2:
        return core, ""
    return parts[0], parts[1]


def build_ddg_rows(
    group_agg: Dict[str, Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """For each target with both WT and ncAA variant, compute ΔΔG."""
    by_target: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for group, agg in group_agg.items():
        target, variant = split_target_variant(group)
        by_target.setdefault(target, {})[variant] = agg
    rows = []
    for target, variants in sorted(by_target.items()):
        wt = variants.get("WT")
        ncaa_variants = {k: v for k, v in variants.items() if k != "WT"}
        if wt is None or not ncaa_variants:
            continue
        for ncaa_label, nc in sorted(ncaa_variants.items()):
            if wt["mean_dg"] is None or nc["mean_dg"] is None:
                continue
            ddg = nc["mean_dg"] - wt["mean_dg"]
            # σ_within_avg may be None if no IE data at all; treat as 0.
            sw_wt = wt["sigma_within_avg"] or 0.0
            sw_nc = nc["sigma_within_avg"] or 0.0
            sigma_ddg = (sw_wt ** 2 + sw_nc ** 2) ** 0.5
            sigma_ddg_between = (
                wt["sigma_between"] ** 2 + nc["sigma_between"] ** 2
            ) ** 0.5
            rows.append({
                "target": target,
                "ncaa_variant": ncaa_label,
                "dg_wt": wt["mean_dg"],
                "sigma_wt_between": wt["sigma_between"],
                "dg_ncaa": nc["mean_dg"],
                "sigma_ncaa_between": nc["sigma_between"],
                "ddg": ddg,
                "sigma_ddg_within": sigma_ddg,
                "sigma_ddg_between": sigma_ddg_between,
                "ie_usable_wt": wt["all_ie_usable"],
                "ie_usable_ncaa": nc["all_ie_usable"],
                "n_rep_wt": wt["n_rep"],
                "n_rep_ncaa": nc["n_rep"],
            })
    return rows


def fmt(x: Optional[float], width: int = 10, prec: int = 2) -> str:
    if x is None:
        return "NA".rjust(width)
    return f"{x:{width}.{prec}f}"


def print_per_run_table(runs: List[Dict[str, Any]]) -> None:
    header = (
        f"{'tag':<30} {'n':>3} {'mean_dg':>10} {'sigma_wi':>10} "
        f"{'-TdS':>10} {'dG_corr':>10} {'usable':>6} {'dG_pref':>10}"
    )
    print("=" * len(header))
    print("PER-RUN TABLE (R-11 ranking-only)")
    print("=" * len(header))
    print(header)
    print("-" * len(header))
    for r in sorted(runs, key=lambda x: x["tag"]):
        print(
            f"{r['tag']:<30} {r['n_snap']:>3} "
            f"{fmt(r['mean_dg'])} {fmt(r['sigma_wi'])} "
            f"{fmt(r['minus_TdS'])} {fmt(r['delta_g_corrected'])} "
            f"{str(r['ie_usable']):>6} {fmt(r['delta_g_preferred'])}"
        )


def print_ddg_table(rows: List[Dict[str, Any]]) -> None:
    header = (
        f"{'target':<6} {'ncaa':<8} {'<dG_WT>':>10} {'sig_bWT':>8} "
        f"{'<dG_nc>':>10} {'sig_bnc':>8} {'DDG':>10} "
        f"{'sig_DDGw':>10} {'sig_DDGb':>10} {'IEwt':>5} {'IEnc':>5}"
    )
    print()
    print("=" * len(header))
    print("PER-TARGET DDG TABLE (ncAA - WT)")
    print("=" * len(header))
    print(header)
    print("-" * len(header))
    for r in sorted(rows, key=lambda x: (x["target"], x["ncaa_variant"])):
        print(
            f"{r['target']:<6} {r['ncaa_variant']:<8} "
            f"{fmt(r['dg_wt'])} {fmt(r['sigma_wt_between'], 8)} "
            f"{fmt(r['dg_ncaa'])} {fmt(r['sigma_ncaa_between'], 8)} "
            f"{fmt(r['ddg'])} {fmt(r['sigma_ddg_within'])} "
            f"{fmt(r['sigma_ddg_between'])} "
            f"{str(r['ie_usable_wt']):>5} {str(r['ie_usable_ncaa']):>5}"
        )


def main() -> None:
    paths = sorted(glob.glob(GLOB_PATTERN))
    runs = [extract_run(p) for p in paths]

    # Group by target_variant (e.g., '2QKI_Cp4_calib')
    groups: Dict[str, List[Dict[str, Any]]] = {}
    for r in runs:
        groups.setdefault(r["group"], []).append(r)
    group_agg = {g: aggregate_group(rs) for g, rs in groups.items()}
    ddg_rows = build_ddg_rows(group_agg)

    print_per_run_table(runs)
    print_ddg_table(ddg_rows)

    # Summary line
    n_usable = sum(1 for r in runs if r["ie_usable"])
    print()
    print(f"Summary: {len(runs)} runs processed, "
          f"ie_usable={n_usable} vs not_usable={len(runs) - n_usable} "
          f"(threshold sigma_wi <= {SIGMA_THRESHOLD} kcal/mol)")

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        json.dump({
            "regime": "R-11 ranking-only",
            "sigma_threshold_kcal": SIGMA_THRESHOLD,
            "source": "Jiao-Lin 2021",
            "n_runs": len(runs),
            "n_ie_usable": n_usable,
            "per_run": runs,
            "per_group": group_agg,
            "per_target_ddg": ddg_rows,
        }, fh, indent=2)
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
