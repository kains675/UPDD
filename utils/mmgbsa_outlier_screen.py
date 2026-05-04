"""
MAD-based outlier screening for MM-GBSA per-snap results.

Production helper for Path-recommended 'MAD-any-of' sanitization policy.
Supersedes Quick Win 3 ad-hoc script (scripts/screen_outliers.py).

Reference:
    Hampel FR (1974) "The Influence Curve and its Role in Robust Estimation."
    J. Am. Stat. Assoc. 69(346), 383-393. doi:10.1080/01621459.1974.10482962

Background
----------
Classical z-score (x - mean)/stdev is non-robust: a single pathological
outlier inflates stdev and masks its own deviation. The median absolute
deviation (MAD = median(|x_i - median(x)|)) is a 50%-breakdown-point
robust scale estimator. Scaling MAD / 0.6745 produces a consistent
Gaussian-stdev estimator (0.6745 = Phi^{-1}(0.75)), so the familiar
3-sigma cutoff maps to ~0.3% two-tailed Gaussian tail probability.

Motivating cases (from Path verdicts, 2026-04-22):
  - 2QKI_Cp4_s7  snap08_f176  : e_ligand ~ +26415 kcal/mol (neighbours -550..-604)
  - 2QKI_Cp4_s23 snap17_f303  : delta_g -3985, e_ligand +3638 (|z_eL|=179.7)
  - 3IOL_NML20_s101 snap25_f490 : delta_g +1160 (|z_dg|=264.9) flips ensemble sign

Policy: 'MAD-any-of' -- flag snap if any of (e_complex, e_ligand, delta_g)
has |robust-z| > threshold (default 3.0).

Schema contract:
  - Sanitized summary JSON gets top-level `sanitization_policy: "mad_any_of_v1"`.
  - Sanitized summary JSON gets top-level `sanitization` block with audit fields.
  - Original `results[]` is preserved unchanged (R-7).
  - Original summary JSONs are never modified in place.
"""

from __future__ import annotations

import json
import math
import os
import statistics
from typing import Any, Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
MAD_CONSISTENCY = 0.6745  # Phi^{-1}(0.75); MAD / 0.6745 ~= Gaussian stdev.
DEFAULT_FIELDS: Tuple[str, ...] = (
    "e_complex_kcal",
    "e_ligand_kcal",
    "delta_g_kcal",
)
DEFAULT_Z_THRESHOLD: float = 3.0
SANITIZATION_POLICY_ID: str = "mad_any_of_v1"


# ---------------------------------------------------------------------------
# Core robust statistics
# ---------------------------------------------------------------------------
def robust_z_scores(
    values: List[float],
) -> Tuple[float, float, List[float], bool]:
    """Compute robust z-scores using MAD-derived Gaussian-equivalent scale.

    Parameters
    ----------
    values : list of float
        Univariate sample.

    Returns
    -------
    median : float
        Sample median.
    mad_scaled : float
        MAD / 0.6745 (approx. Gaussian stdev). Falls back to pstdev if
        MAD == 0 (degenerate case).
    z_scores : list of float
        (x_i - median) / mad_scaled; zero-vector if scale is zero.
    mad_degenerate : bool
        True if the raw MAD was zero (pstdev fallback used).
    """
    if not values:
        return 0.0, 0.0, [], False
    med = statistics.median(values)
    mad_raw = statistics.median([abs(v - med) for v in values])
    if mad_raw == 0.0:
        # Degenerate — majority of values collapse onto median. Fall back
        # to pstdev so a genuine single outlier is still detectable.
        scale = statistics.pstdev(values) if len(values) > 1 else 0.0
        if scale == 0.0:
            return med, 0.0, [0.0] * len(values), True
        return med, scale, [(v - med) / scale for v in values], True
    scale = mad_raw / MAD_CONSISTENCY
    return med, scale, [(v - med) / scale for v in values], False


# ---------------------------------------------------------------------------
# Screening
# ---------------------------------------------------------------------------
def mad_any_of_screen(
    snap_results: List[Dict[str, Any]],
    fields: Tuple[str, ...] = DEFAULT_FIELDS,
    z_threshold: float = DEFAULT_Z_THRESHOLD,
) -> Dict[str, Any]:
    """Apply 'MAD-any-of' outlier screen across chosen energy fields.

    A snap is flagged if |robust-z| > z_threshold for AT LEAST ONE of the
    requested fields. This catches pathologies that manifest in delta_g
    but are hidden in e_ligand alone (and vice versa).

    Parameters
    ----------
    snap_results : list of dict
        Each dict is one entry from a mmgbsa_summary.json `results[]`.
    fields : tuple of str
        Energy fields to screen. Default screens complex, ligand, delta_g.
    z_threshold : float
        |z| cutoff. Default 3.0 (Hampel-standard 3-sigma equivalent).

    Returns
    -------
    dict
        Report with per-field z arrays, flagged snap detail, retained
        index list, and policy identifier.
    """
    n = len(snap_results)
    per_field_stats: Dict[str, Dict[str, Any]] = {}
    per_field_z: Dict[str, List[float]] = {}
    for field in fields:
        values = [float(r[field]) for r in snap_results]
        med, scale, zs, degen = robust_z_scores(values)
        per_field_z[field] = zs
        per_field_stats[field] = {
            "median": med,
            "scale_gaussian_stdev": scale,
            "mad_degenerate": degen,
        }

    flagged: List[Dict[str, Any]] = []
    retained_indices: List[int] = []
    flagged_snap_names: List[str] = []
    for i, r in enumerate(snap_results):
        z_per_field = {f: per_field_z[f][i] for f in fields}
        offending = [f for f, z in z_per_field.items() if abs(z) > z_threshold]
        if offending:
            max_field = max(offending, key=lambda f: abs(z_per_field[f]))
            flagged.append(
                {
                    "snap": r.get("snapshot", f"index_{i}"),
                    "index": i,
                    "field_max_z": max_field,
                    "max_abs_z": abs(z_per_field[max_field]),
                    "z_per_field": {f: round(z, 4) for f, z in z_per_field.items()},
                    "offending_fields": offending,
                    "values": {f: float(r[f]) for f in fields},
                }
            )
            flagged_snap_names.append(r.get("snapshot", f"index_{i}"))
        else:
            retained_indices.append(i)

    return {
        "fields_screened": list(fields),
        "z_threshold": z_threshold,
        "n_snap_original": n,
        "n_snap_flagged": len(flagged),
        "n_snap_retained": len(retained_indices),
        "flagged_snaps": flagged,
        "flagged_snap_names": flagged_snap_names,
        "retained_indices": retained_indices,
        "per_field_stats": per_field_stats,
        "sanitization_policy": SANITIZATION_POLICY_ID,
    }


# ---------------------------------------------------------------------------
# Summary-level sanitization
# ---------------------------------------------------------------------------
def _mean_or_none(xs: List[float]) -> Optional[float]:
    if not xs:
        return None
    return float(statistics.fmean(xs))


def _pstdev_or_zero(xs: List[float]) -> float:
    if len(xs) < 2:
        return 0.0
    return float(statistics.pstdev(xs))


def _recompute_interaction_entropy(
    delta_g_values: List[float],
    temperature_K: float = 300.0,
) -> Optional[Dict[str, Any]]:
    """Duan 2016 JCTC -TΔS from ΔG fluctuations.

    -TΔS_IE = (1 / β) * ln<exp(β * (ΔG_i - <ΔG>))>, β = 1 / (kB * T).
    Uses kcal/mol units; kB*T at 300 K ~= 0.5961 kcal/mol.

    Returns None if fewer than 2 samples (IE ill-defined).
    """
    if len(delta_g_values) < 2:
        return None
    kBT = 0.0019872041 * temperature_K  # kcal/mol
    mean_dg = statistics.fmean(delta_g_values)
    # Numerically-stable mean of exp((dg - mean) / kBT).
    exps = [math.exp((dg - mean_dg) / kBT) for dg in delta_g_values]
    mean_exp = statistics.fmean(exps)
    minus_T_dS = kBT * math.log(mean_exp)
    return {
        "mean_delta_g_kcal": mean_dg,
        "std_delta_g_kcal": _pstdev_or_zero(delta_g_values),
        "minus_TdS_kcal": minus_T_dS,
        "delta_g_corrected_kcal": mean_dg + minus_T_dS,
        "n_frames": len(delta_g_values),
        "reliable": True,
        "temperature_K": temperature_K,
    }


def sanitize_mmgbsa_summary(
    summary_path: str,
    fields: Optional[Tuple[str, ...]] = None,
    z_threshold: float = DEFAULT_Z_THRESHOLD,
    write_sanitized: bool = False,
    sanitized_path: Optional[str] = None,
) -> Dict[str, Any]:
    """Read a mmgbsa_summary.json, apply MAD-any-of screen, and return a
    sanitized summary dict.

    Parameters
    ----------
    summary_path : str
        Path to original mmgbsa_summary.json.
    fields : tuple of str, optional
        Energy fields to screen (defaults to DEFAULT_FIELDS).
    z_threshold : float
        |z| cutoff (default 3.0).
    write_sanitized : bool
        If True, write the sanitized dict to sanitized_path.
    sanitized_path : str, optional
        Output path. Defaults to summary_path with '.json' -> '_sanitized.json'.

    Returns
    -------
    dict
        Sanitized summary. Contains:
          * original top-level keys (gb_model, ncaa_elem, si_radius, n_calc,
            results, etc.) preserved verbatim
          * mean_dg and best_dg recomputed over retained snaps
          * sanitization_policy == "mad_any_of_v1" (top-level)
          * sanitization : {policy, z_threshold, fields_screened,
            n_snap_original, n_snap_retained, flagged_snap_names,
            mean_dg_original, mean_dg_sanitized,
            sigma_wi_original, sigma_wi_sanitized, delta_mean_dg,
            per_field_stats}
          * interaction_entropy_per_design recomputed over retained snaps
            (if originally present)
    """
    if fields is None:
        fields = DEFAULT_FIELDS

    with open(summary_path, "r") as f:
        summary = json.load(f)

    results = summary.get("results", [])
    if not isinstance(results, list):
        raise ValueError(
            f"{summary_path}: 'results' missing or not a list."
        )

    report = mad_any_of_screen(results, fields=fields, z_threshold=z_threshold)

    # Recompute aggregate values over retained snaps.
    dgs_original = [float(r["delta_g_kcal"]) for r in results]
    retained = [results[i] for i in report["retained_indices"]]
    dgs_sanitized = [float(r["delta_g_kcal"]) for r in retained]

    mean_dg_original = _mean_or_none(dgs_original)
    mean_dg_sanitized = _mean_or_none(dgs_sanitized)
    sigma_wi_original = _pstdev_or_zero(dgs_original)
    sigma_wi_sanitized = _pstdev_or_zero(dgs_sanitized)

    best_dg_sanitized = min(dgs_sanitized) if dgs_sanitized else None

    sanitized = dict(summary)  # shallow copy preserves original results[]
    sanitized["sanitization_policy"] = SANITIZATION_POLICY_ID
    sanitized["sanitization"] = {
        "policy": SANITIZATION_POLICY_ID,
        "z_threshold": z_threshold,
        "fields_screened": list(fields),
        "n_snap_original": report["n_snap_original"],
        "n_snap_retained": report["n_snap_retained"],
        "n_snap_flagged": report["n_snap_flagged"],
        "flagged_snap_names": report["flagged_snap_names"],
        "flagged_snaps": report["flagged_snaps"],
        "mean_dg_original": mean_dg_original,
        "mean_dg_sanitized": mean_dg_sanitized,
        "sigma_wi_original": sigma_wi_original,
        "sigma_wi_sanitized": sigma_wi_sanitized,
        "delta_mean_dg": (
            (mean_dg_sanitized - mean_dg_original)
            if (mean_dg_sanitized is not None and mean_dg_original is not None)
            else None
        ),
        "per_field_stats": report["per_field_stats"],
        "source_summary_path": os.path.abspath(summary_path),
    }
    sanitized["mean_dg"] = mean_dg_sanitized
    sanitized["best_dg"] = best_dg_sanitized

    # Recompute interaction_entropy_per_design if the original had it.
    ie_block = summary.get("interaction_entropy_per_design")
    if isinstance(ie_block, dict):
        new_ie: Dict[str, Any] = {}
        # Group retained snaps by prefix-stripped design key matching
        # original grouping (key == snapshot prefix up to "_snapXX_*").
        for design_key, orig_entry in ie_block.items():
            temperature_K = float(orig_entry.get("temperature_K", 300.0))
            # Select retained snaps whose snapshot starts with design_key.
            group_dgs = [
                float(r["delta_g_kcal"])
                for r in retained
                if str(r.get("snapshot", "")).startswith(design_key)
            ]
            if not group_dgs:
                # Preserve original-shape placeholder but mark unreliable.
                new_ie[design_key] = {
                    "mean_delta_g_kcal": None,
                    "std_delta_g_kcal": None,
                    "minus_TdS_kcal": None,
                    "delta_g_corrected_kcal": None,
                    "n_frames": 0,
                    "reliable": False,
                    "temperature_K": temperature_K,
                    "note": "all snaps flagged by mad_any_of_v1",
                }
                continue
            recomputed = _recompute_interaction_entropy(
                group_dgs, temperature_K=temperature_K
            )
            if recomputed is None:
                new_ie[design_key] = {
                    "mean_delta_g_kcal": statistics.fmean(group_dgs),
                    "std_delta_g_kcal": _pstdev_or_zero(group_dgs),
                    "minus_TdS_kcal": None,
                    "delta_g_corrected_kcal": None,
                    "n_frames": len(group_dgs),
                    "reliable": False,
                    "temperature_K": temperature_K,
                    "note": "fewer than 2 retained snaps; IE undefined",
                }
            else:
                new_ie[design_key] = recomputed
        sanitized["interaction_entropy_per_design"] = new_ie

    if write_sanitized:
        if sanitized_path is None:
            base, ext = os.path.splitext(summary_path)
            sanitized_path = f"{base}_sanitized{ext or '.json'}"
        os.makedirs(os.path.dirname(sanitized_path) or ".", exist_ok=True)
        with open(sanitized_path, "w") as f:
            json.dump(sanitized, f, indent=2)
        sanitized["sanitization"]["sanitized_path"] = os.path.abspath(
            sanitized_path
        )

    return sanitized


__all__ = [
    "MAD_CONSISTENCY",
    "DEFAULT_FIELDS",
    "DEFAULT_Z_THRESHOLD",
    "SANITIZATION_POLICY_ID",
    "robust_z_scores",
    "mad_any_of_screen",
    "sanitize_mmgbsa_summary",
]
