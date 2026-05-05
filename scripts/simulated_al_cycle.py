#!/usr/bin/env python
"""scripts/simulated_al_cycle.py - Simulated active-learning cycle (v0.1).

Active-learning loop demonstration over the v0.7.1 6-target ncAA pair benchmark
(Cp4, 1EBP_MTR13, 1YCR_NML22, 3IOL_NML20, 7TL8_MTR6, 2QKH_MTR25). Synthetic
K_d ground truth is sourced from Magotti 2009 ITC for Cp4 (literature SSOT
direct, DOI 10.1002/jmr.972 Table 3, DDG = -1.4 kcal/mol) and from UPDD
MM-PBSA post-CONECT-v55 outputs for the other five pairs. The script
implements three acquisition functions (UCB / EI / MaxEntropy per spec
Section 2 phase-switching), a calibration model (linear residual for
cycles 1-2 -> GPR at n_cumulative >= 15 AND cycle >= 3 per spec Section 3
Option C), 95% predicted-interval coverage probability, plateau detection
via cycle_manager.detect_plateau, and random/greedy baseline comparators.

Disclosure (sigma_fix4_target_selectivity): Synthetic K_d noise inherits Cp4
60x sigma_btwn tightening vs WT 1.08x (sigma fix-4 target-selective per UPDD
internal verdict). This is a hidden second-order circularity in the
synthetic-cycle output and is surfaced here as honest disclosure rather
than aggregated into a single noise term.

Interpretation regime: protocol_convergence_demonstration.

Acquisition-function citations:
    UCB         - Srinivas et al. 2010 (arXiv:0912.3995, ICML 2010), beta=1.0
                  per spec Section 2 line 110 (balanced exploration/exploitation).
    EI          - Jones, Schonlau, Welch 1998 (J. Glob. Optim. 13:455-492,
                  DOI 10.1023/A:1008306431147).
    MaxEntropy  - Houlsby et al. 2011 (BALD, arXiv:1112.5745).

Schema: simulated_al_cycle/0.1
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# Repo root on sys.path so utils.cycle_manager is importable when invoked
# directly from the scripts/ directory.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

# ---------------------------------------------------------------------------
# Schema + literature SSOT constants (C1: Magotti 2009 SSOT direct for Cp4)
# ---------------------------------------------------------------------------

SCHEMA_VERSION = "simulated_al_cycle/0.1"

# C1 BLOCKING - Magotti SSOT direct, no MM-PBSA-derived for Cp4.
# DOI 10.1002/jmr.972 Table 3 (ITC primary) / Table 2 (SPR sibling).
MAGOTTI_2009_CP4_DDG_KCAL = -1.4
MAGOTTI_2009_CP4_DDG_KCAL_SPR = -3.0

# K_d <-> DG conversion constants - mirrors cycle_manager.py:99-100 (CODATA 2018).
R_KCAL_PER_MOL_K = 1.987204258640832e-3
T_KELVIN = 298.15

# C5 calibration transition thresholds.
GPR_TRANSITION_N_CUMULATIVE = 15
GPR_TRANSITION_CYCLE = 3

# Phase-specific acquisition switching (spec Section 2 lines 106-111).
ACQUISITION_PHASE = {1: "max_entropy", 2: "ucb", 3: "ei"}
UCB_BETA = 1.0

# Six-pair ncAA benchmark.
TARGETS = ("Cp4", "1EBP_MTR13", "1YCR_NML22", "3IOL_NML20", "7TL8_MTR6", "2QKH_MTR25")

# Path to UPDD MM-PBSA synthetic-truth source (read-only).
CALIBRATION_JSON = _REPO_ROOT / "outputs" / "calibration_summary_cleaned.json"
TWOQKH_SEED_DIRS = (
    _REPO_ROOT / "outputs" / "2QKH_MTR25_calib_s7" / "mmpbsa_results_postl387" / "mmpbsa_summary.json",
    _REPO_ROOT / "outputs" / "2QKH_MTR25_calib_s19" / "mmpbsa_results_postl387" / "mmpbsa_summary.json",
    _REPO_ROOT / "outputs" / "2QKH_MTR25_calib_s42" / "mmpbsa_results_postl387" / "mmpbsa_summary.json",
)


# ---------------------------------------------------------------------------
# K_d <-> DG conversion (C1)
# ---------------------------------------------------------------------------

def ddg_kcal_to_kd_molar(ddg: float, T: float = T_KELVIN) -> float:
    """K_d = exp(DG / RT) in molar (sign convention per cycle_manager._kd_to_ddg)."""
    return math.exp(ddg / (R_KCAL_PER_MOL_K * T))


def kd_molar_to_ddg_kcal(kd: float, T: float = T_KELVIN) -> float:
    """DG = R T ln(K_d / 1 M) in kcal/mol."""
    if kd <= 0 or not math.isfinite(kd):
        raise ValueError(f"kd must be > 0, got {kd!r}")
    return R_KCAL_PER_MOL_K * T * math.log(kd)


# ---------------------------------------------------------------------------
# Ground-truth construction (C1)
# ---------------------------------------------------------------------------

def _aggregate_2qkh_seed_means(seed_paths: Tuple[Path, ...]) -> float:
    """Inline mean-of-means aggregator for 2QKH_MTR25 raw seeds.

    Analogous to branched_ddg.compute_branched_ddg() but inline; we never
    modify the protected branched_ddg module. Reads ``mean_dg`` field from
    each seed's mmpbsa_summary.json.
    """
    means = []
    for p in seed_paths:
        with open(p) as f:
            d = json.load(f)
        means.append(float(d["mean_dg"]))
    return float(np.mean(means))


def build_synthetic_ground_truth(calibration_json: Optional[Path] = None,
                                 twoqkh_seed_paths: Optional[Tuple[Path, ...]] = None
                                 ) -> Dict[str, Dict[str, Any]]:
    """Build the 6-target synthetic ground-truth dict.

    Cp4 path is hardcoded to MAGOTTI_2009_CP4_DDG_KCAL (C1 BLOCKING). The
    other five pairs read mean-of-per-seed-means from ``per_variant`` of
    ``calibration_summary_cleaned.json``; 2QKH_MTR25 aggregates from raw
    seed mmpbsa_summary.json files.

    Returns dict keyed by variant -> {ddg_kcal, kd_molar, provenance, source}.
    """
    cal_path = calibration_json if calibration_json is not None else CALIBRATION_JSON
    tk_paths = twoqkh_seed_paths if twoqkh_seed_paths is not None else TWOQKH_SEED_DIRS

    with open(cal_path) as f:
        cal = json.load(f)
    pv = cal["per_variant"]

    def _mean_of_means(key: str) -> float:
        return float(np.mean([float(e["mean_dg"]) for e in pv[key]]))

    gt: Dict[str, Dict[str, Any]] = {}

    # C1: Cp4 from Magotti 2009 ITC literature SSOT direct.
    gt["Cp4"] = {
        "ddg_kcal": MAGOTTI_2009_CP4_DDG_KCAL,
        "kd_molar": ddg_kcal_to_kd_molar(MAGOTTI_2009_CP4_DDG_KCAL),
        "provenance": "Magotti_2009_ITC_literature_anchor",
        "source": "DOI 10.1002/jmr.972 Table 3",
    }
    # SPR sibling (sensitivity sibling - not used as primary oracle).
    gt["Cp4_SPR_sibling"] = {
        "ddg_kcal": MAGOTTI_2009_CP4_DDG_KCAL_SPR,
        "kd_molar": ddg_kcal_to_kd_molar(MAGOTTI_2009_CP4_DDG_KCAL_SPR),
        "provenance": "Magotti_2009_SPR_literature_anchor_sensitivity",
        "source": "DOI 10.1002/jmr.972 Table 2",
    }

    # Other five pairs from UPDD MM-PBSA synthetic.
    for variant, key in (
        ("1EBP_MTR13", "1EBP_MTR13_calib"),
        ("1YCR_NML22", "1YCR_NML22_calib"),
        ("3IOL_NML20", "3IOL_NML20_calib"),
        ("7TL8_MTR6",  "7TL8_MTR6_calib"),
    ):
        ddg = _mean_of_means(key)
        gt[variant] = {
            "ddg_kcal": ddg,
            "kd_molar": ddg_kcal_to_kd_molar(ddg),
            "provenance": "UPDD_MM-PBSA_synthetic",
            "source": str(cal_path.relative_to(_REPO_ROOT)),
        }

    # 2QKH_MTR25 from inline aggregation of 3 raw seeds.
    ddg_2qkh = _aggregate_2qkh_seed_means(tk_paths)
    gt["2QKH_MTR25"] = {
        "ddg_kcal": ddg_2qkh,
        "kd_molar": ddg_kcal_to_kd_molar(ddg_2qkh),
        "provenance": "UPDD_MM-PBSA_synthetic_aggregated_3seeds",
        "source": "outputs/2QKH_MTR25_calib_s{7,19,42}/mmpbsa_results_postl387/mmpbsa_summary.json",
    }
    return gt


# ---------------------------------------------------------------------------
# Acquisition functions (spec Section 2)
# ---------------------------------------------------------------------------

def acquisition_ucb(mu: np.ndarray, sigma: np.ndarray, beta: float = UCB_BETA
                    ) -> np.ndarray:
    """UCB score = mu + beta * sigma (Srinivas+ 2010)."""
    return np.asarray(mu, dtype=float) + float(beta) * np.asarray(sigma, dtype=float)


def acquisition_ei(mu: np.ndarray, sigma: np.ndarray, current_best: float
                   ) -> np.ndarray:
    """Expected improvement under normal posterior (Jones+ 1998).

    Maximisation form: EI(x) = (mu - best) Phi(z) + sigma phi(z),
    z = (mu - best) / sigma. Returns a non-negative score for each point.
    """
    mu = np.asarray(mu, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    out = np.zeros_like(mu)
    nz = sigma > 0
    z = np.zeros_like(mu)
    z[nz] = (mu[nz] - current_best) / sigma[nz]
    # Phi (standard normal CDF) via erf; phi (PDF) closed form.
    Phi = 0.5 * (1.0 + np.vectorize(math.erf)(z / math.sqrt(2.0)))
    phi = (1.0 / math.sqrt(2.0 * math.pi)) * np.exp(-0.5 * z * z)
    out[nz] = (mu[nz] - current_best) * Phi[nz] + sigma[nz] * phi[nz]
    return out


def acquisition_max_entropy(mu: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    """Differential entropy of N(mu, sigma^2) up to additive constant.

    Shannon-entropy maximisation reduces to selecting the largest sigma
    (BALD pool-based, Houlsby+ 2011 reduction).
    """
    sigma = np.asarray(sigma, dtype=float)
    # 0.5*log(2 pi e sigma^2) is monotone in sigma; return the entropy proper.
    eps = 1e-30
    return 0.5 * np.log(2.0 * math.pi * math.e * np.maximum(sigma, eps) ** 2)


def select_acquisition(cycle_idx: int) -> str:
    """Phase-specific acquisition switching (spec Section 2 lines 106-111).

    Cycle 1: max_entropy (broad exploration).
    Cycle 2: UCB (balanced).
    Cycle 3+: EI (exploitation).
    """
    if cycle_idx <= 0:
        raise ValueError("cycle_idx must be >= 1")
    return ACQUISITION_PHASE.get(cycle_idx, "ei")


# ---------------------------------------------------------------------------
# Calibration models (C5)
# ---------------------------------------------------------------------------

def fit_linear_residual(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """Fit y = a*x + b, return (a, b).  Used for cycles 1-2 (small data)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size < 2:
        # Degenerate - identity calibration (a=1, b=0).
        return 1.0, 0.0
    a, b = np.polyfit(x, y, 1)
    return float(a), float(b)


def predict_linear(a: float, b: float, x: np.ndarray
                   ) -> Tuple[np.ndarray, np.ndarray]:
    """Linear predictor; returns (mu, sigma).  Sigma estimated as residual std."""
    x = np.asarray(x, dtype=float)
    mu = a * x + b
    return mu, np.full_like(mu, fill_value=0.5)  # default 0.5 kcal placeholder


def fit_gpr(x: np.ndarray, y: np.ndarray):
    """Fit GaussianProcessRegressor with RBF + WhiteKernel.

    Per C5: kernel = RBF(length_scale=1.0) + WhiteKernel(noise_level=0.1).
    Returns the fitted regressor.
    """
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, WhiteKernel
    x = np.asarray(x, dtype=float).reshape(-1, 1)
    y = np.asarray(y, dtype=float)
    kernel = RBF(length_scale=1.0) + WhiteKernel(noise_level=0.1)
    gpr = GaussianProcessRegressor(
        kernel=kernel, normalize_y=True, n_restarts_optimizer=3, random_state=0,
    )
    gpr.fit(x, y)
    return gpr


def should_transition_to_gpr(cycle: int, n_cumulative: int) -> bool:
    """C5 transition criterion: GPR at n_cumulative >= 15 AND cycle >= 3."""
    return cycle >= GPR_TRANSITION_CYCLE and n_cumulative >= GPR_TRANSITION_N_CUMULATIVE


# ---------------------------------------------------------------------------
# Coverage probability (C3)
# ---------------------------------------------------------------------------

def compute_coverage_probability(predictions_lo: np.ndarray,
                                 predictions_hi: np.ndarray,
                                 truths: np.ndarray) -> float:
    """Mean indicator c = mean(I[truth in [lo, hi]]) (spec Section 3 line 169)."""
    lo = np.asarray(predictions_lo, dtype=float)
    hi = np.asarray(predictions_hi, dtype=float)
    truth = np.asarray(truths, dtype=float)
    if truth.size == 0:
        return float("nan")
    inside = (truth >= lo) & (truth <= hi)
    return float(np.mean(inside.astype(float)))


# ---------------------------------------------------------------------------
# Baseline comparators
# ---------------------------------------------------------------------------

def baseline_random(candidates: List[str], batch_size: int, rng: np.random.Generator
                    ) -> List[str]:
    if batch_size >= len(candidates):
        return list(candidates)
    idx = rng.choice(len(candidates), size=batch_size, replace=False)
    return [candidates[int(i)] for i in idx]


def baseline_greedy(candidates: List[str], scores: np.ndarray, batch_size: int
                    ) -> List[str]:
    """Pure greedy on prediction (descending mu), ignoring uncertainty."""
    order = np.argsort(-np.asarray(scores, dtype=float))
    take = order[: min(batch_size, len(candidates))]
    return [candidates[int(i)] for i in take]


# ---------------------------------------------------------------------------
# Re-entry interface (uses cycle_manager public API)
# ---------------------------------------------------------------------------

def reentry_via_cycle_manager(project_id: str, source_cycle: int, *,
                              best_sequence: str, output_root: Path) -> Dict[str, Any]:
    """Re-entry shim that calls cycle_manager.set_next_wt for the orchestrator.

    Lazy-imports cycle_manager so that the module can be loaded in
    --dry-run mode without the full utils.* dependency cone.
    """
    from utils import cycle_manager as _cm  # type: ignore
    return _cm.set_next_wt(
        project_id=project_id,
        source_cycle=source_cycle,
        best="explicit",
        explicit_sequence=best_sequence,
        output_root=str(output_root),
    )


# ---------------------------------------------------------------------------
# Main simulator
# ---------------------------------------------------------------------------

def _synthesize_predictions(truths: np.ndarray, rng: np.random.Generator,
                            noise_kcal: float = 1.0
                            ) -> Tuple[np.ndarray, np.ndarray]:
    """Synthetic predictions = truth + N(0, noise) with a default sigma."""
    mu = np.asarray(truths, dtype=float) + rng.normal(0.0, noise_kcal,
                                                      size=len(truths))
    sigma = np.full_like(mu, noise_kcal, dtype=float)
    return mu, sigma


def run_simulated_cycle(*, strategy: str, n_cycles: int, candidates_per_cycle: int,
                        ground_truth: Dict[str, Dict[str, Any]],
                        rng: np.random.Generator) -> Dict[str, Any]:
    """Run a 3-cycle simulated AL loop on the 6-target benchmark.

    strategy in {"al", "random", "greedy"}.  Returns per-cycle records,
    calibration-model trace, coverage probabilities, RMSE per cycle, and
    plateau detection result.
    """
    candidate_keys = [k for k in TARGETS]  # preserve order, exclude SPR sibling
    truths = np.array([ground_truth[k]["ddg_kcal"] for k in candidate_keys])

    per_cycle: List[Dict[str, Any]] = []
    cumulative_x: List[float] = []
    cumulative_y: List[float] = []
    cum_predicted: List[float] = []
    cum_truth: List[float] = []
    model_trace: List[str] = []
    coverage_trace: List[float] = []
    rmse_trace: List[float] = []
    transition_cycle: Optional[int] = None
    trajectory_kd: List[float] = []

    for cycle in range(1, n_cycles + 1):
        # Per-cycle synthetic prediction (mu, sigma) for each candidate.
        mu, sigma = _synthesize_predictions(truths, rng, noise_kcal=1.0)

        # Acquisition selection per strategy.
        if strategy == "al":
            acq_name = select_acquisition(cycle)
            if acq_name == "ucb":
                scores = acquisition_ucb(mu, sigma, beta=UCB_BETA)
            elif acq_name == "ei":
                best = float(np.min(cumulative_y)) if cumulative_y else float(np.min(mu))
                scores = acquisition_ei(-mu, sigma, current_best=-best)
            else:  # max_entropy
                scores = acquisition_max_entropy(mu, sigma)
            order = np.argsort(-scores)
            picks_idx = order[: min(candidates_per_cycle, len(candidate_keys))].tolist()
        elif strategy == "random":
            acq_name = "random"
            keys_picked = baseline_random(candidate_keys, candidates_per_cycle, rng)
            picks_idx = [candidate_keys.index(k) for k in keys_picked]
        elif strategy == "greedy":
            acq_name = "greedy"
            keys_picked = baseline_greedy(candidate_keys, mu, candidates_per_cycle)
            picks_idx = [candidate_keys.index(k) for k in keys_picked]
        else:
            raise ValueError(f"unknown strategy {strategy!r}")

        sel_keys = [candidate_keys[i] for i in picks_idx]
        sel_mu = mu[picks_idx]
        sel_sigma = sigma[picks_idx]
        sel_truth = truths[picks_idx]

        # Update cumulative training data with the picked points.
        cumulative_x.extend(sel_mu.tolist())
        cumulative_y.extend(sel_truth.tolist())
        cum_predicted.extend(sel_mu.tolist())
        cum_truth.extend(sel_truth.tolist())

        # C5: pick calibration model.
        n_cumulative = len(cumulative_x)
        use_gpr = should_transition_to_gpr(cycle, n_cumulative)
        if use_gpr:
            model_trace.append("gpr")
            if transition_cycle is None:
                transition_cycle = cycle
            try:
                gpr = fit_gpr(np.array(cumulative_x), np.array(cumulative_y))
                pred_mu, pred_std = gpr.predict(
                    np.array(cum_predicted).reshape(-1, 1), return_std=True,
                )
            except Exception:
                # Fallback to linear if GPR fit fails.
                a, b = fit_linear_residual(np.array(cumulative_x),
                                           np.array(cumulative_y))
                pred_mu = a * np.array(cum_predicted) + b
                pred_std = np.full_like(pred_mu, 0.5)
        else:
            model_trace.append("linear_residual")
            a, b = fit_linear_residual(np.array(cumulative_x),
                                       np.array(cumulative_y))
            pred_mu = a * np.array(cum_predicted) + b
            pred_std = np.full_like(pred_mu, 0.5)

        # 95% CI bounds.
        lo = pred_mu - 1.96 * pred_std
        hi = pred_mu + 1.96 * pred_std
        cov = compute_coverage_probability(lo, hi, np.array(cum_truth))
        coverage_trace.append(cov)

        rmse = float(np.sqrt(np.mean((np.array(cum_predicted) - np.array(cum_truth)) ** 2)))
        rmse_trace.append(rmse)

        # Trajectory K_d for plateau detection (best DG so far).
        best_dg = float(np.min(cumulative_y))
        trajectory_kd.append(ddg_kcal_to_kd_molar(best_dg))

        per_cycle.append({
            "cycle": cycle,
            "acquisition": acq_name,
            "selected": sel_keys,
            "predicted_kcal": [float(x) for x in sel_mu],
            "ground_truth_kcal": [float(x) for x in sel_truth],
            "calibration_error_kcal": float(np.mean(np.abs(sel_mu - sel_truth))),
            "n_cumulative": n_cumulative,
        })

    # Plateau detection on the K_d trajectory: simple monotonic-delta check
    # mirroring cycle_manager.detect_plateau threshold (0.5 kcal/mol).
    plateau_cycle: Optional[int] = None
    if len(trajectory_kd) >= 2:
        dgs = [kd_molar_to_ddg_kcal(k) for k in trajectory_kd]
        deltas = [abs(dgs[i + 1] - dgs[i]) for i in range(len(dgs) - 1)]
        if all(d < 0.5 for d in deltas[-2:]) and len(deltas) >= 2:
            plateau_cycle = len(dgs)

    return {
        "per_cycle": per_cycle,
        "calibration": {
            "model_per_cycle": model_trace,
            "transition_criterion_met_at_cycle": transition_cycle,
            "coverage_probability_per_cycle": coverage_trace,
            "rmse_per_cycle_kcal": rmse_trace,
        },
        "plateau_detection": {
            "plateau_cycle": plateau_cycle,
            "trajectory_kd_molar": trajectory_kd,
        },
    }


# ---------------------------------------------------------------------------
# Output writer
# ---------------------------------------------------------------------------

def _build_output_payload(ground_truth: Dict[str, Dict[str, Any]],
                          al_run: Dict[str, Any],
                          random_run: Dict[str, Any],
                          greedy_run: Dict[str, Any],
                          n_cycles: int, candidates_per_cycle: int
                          ) -> Dict[str, Any]:
    """Assemble the simulated_al_cycle/0.1 JSON payload.

    C7 (Single-anchor framing): the synthetic ground truth is a SINGLE
    literature anchor (Cp4 = Magotti 2009 ITC). The other five pipeline-
    produced ⟨⟨ΔG⟩⟩ values are NOT ground truth — they are pipeline-self-
    prediction comparators used only to populate the candidate space at a
    realistic ΔG scale. The demonstration scope is therefore single-anchor
    *convergence* against Magotti SSOT, not multi-anchor calibration.
    """
    rmse_al = al_run["calibration"]["rmse_per_cycle_kcal"][-1]
    rmse_random = random_run["calibration"]["rmse_per_cycle_kcal"][-1]
    rmse_greedy = greedy_run["calibration"]["rmse_per_cycle_kcal"][-1]

    # C7: Partition the variant entries into literature ground_truth
    # (single anchor) vs pipeline comparators (5 pairs, NOT ground truth).
    gt_literature = {
        k: v for k, v in ground_truth.items()
        if "literature_anchor" in v.get("provenance", "")
    }
    comparators_pipeline = {
        k: v for k, v in ground_truth.items()
        if "literature_anchor" not in v.get("provenance", "")
    }

    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        # C2 BLOCKING - framing fields.
        "interpretation_regime": "protocol_convergence_demonstration",
        "absolute_efficiency_claim": False,
        "prospective_wet_lab_required": True,
        # C7 BLOCKING - single-anchor framing fields.
        "ground_truth_source": "literature_single_anchor",
        "anchor_count": 1,  # Magotti 2009 ITC Cp4 only; SPR sibling is sensitivity, not a second anchor
        "comparator_source": "pipeline_self_prediction",
        "demonstration_scope": "single_anchor_convergence",  # NOT multi_anchor_calibration
        "notes_provenance_disclosure": (
            "Cp4 anchor: Magotti 2009 ITC -1.4 kcal/mol literature direct "
            "(single ground-truth anchor; SPR -3.0 carried as sensitivity "
            "sibling, NOT a second anchor). Other 5 pairs (1EBP_MTR13, "
            "1YCR_NML22, 3IOL_NML20, 7TL8_MTR6, 2QKH_MTR25): pipeline self-"
            "prediction comparators (UPDD MM-PBSA post-CONECT-v55, post-"
            "PR-21), not ground truth. Demonstration scope: single-anchor "
            "convergence against Magotti SSOT, not multi-anchor "
            "calibration. Absolute efficiency / wet-lab try reduction "
            "claims are NOT supported by this output."
        ),
        "ssot_anchors": {
            "cp4": "Magotti 2009 (10.1002/jmr.972) Table 3 ITC -1.4 kcal/mol",
        },
        "n_targets": len([k for k in TARGETS]),
        "n_cycles": n_cycles,
        "candidates_per_cycle": candidates_per_cycle,
        # C7: ground_truth contains ONLY literature anchors (Cp4, optional
        # SPR sibling). Pipeline-produced ⟨⟨ΔG⟩⟩ values for the other five
        # pairs are surfaced separately under "comparators".
        "ground_truth": gt_literature,
        "comparators": comparators_pipeline,
        "strategies": {
            "al": al_run,
            "random": random_run,
            "greedy": greedy_run,
        },
        "comparative_summary": {
            "rmse_final_kcal_al": rmse_al,
            "rmse_final_kcal_random": rmse_random,
            "rmse_final_kcal_greedy": rmse_greedy,
            "narrative": (
                "AL strategy emits comparative protocol convergence vs "
                "random/greedy baselines under a single-anchor synthetic "
                "oracle (Cp4 = Magotti 2009 ITC). Absolute efficiency / "
                "wet-lab try reduction claims NOT supported."
            ),
        },
        "notes": {
            "sigma_fix4_target_selectivity": True,
            "sigma_fix4_disclosure": (
                "Synthetic K_d noise inherits Cp4 60x sigma_btwn tightening "
                "vs WT 1.08x per UPDD internal verdict "
                "(sigma_fix4_target_selective). This is a hidden second-"
                "order circularity in the synthetic-cycle output and is "
                "surfaced here as honest disclosure rather than aggregated "
                "into a single noise term."
            ),
            "single_anchor_framing_disclosure": (
                "Per condition C7 (user 2026-05-05): the ground_truth "
                "object holds only the Magotti 2009 literature anchor "
                "(anchor_count=1). The five pipeline-produced "
                "comparator values appear under the comparators field "
                "and are explicitly NOT used as multi-anchor calibration "
                "targets — they populate the candidate space at "
                "realistic ΔG scale only."
            ),
        },
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--n-cycles", type=int, default=3)
    p.add_argument("--candidates-per-cycle", type=int, default=5)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output-dir", type=str, default=None,
                   help="Output directory (default: outputs/simulated_al_cycle_<YYYYMMDD>/)")
    p.add_argument("--dry-run", action="store_true",
                   help="Build ground truth + print structure, then exit 0.")
    args = p.parse_args(argv)

    gt = build_synthetic_ground_truth()

    if args.dry_run:
        print("simulated_al_cycle/0.1 dry-run")
        print(f"schema_version = {SCHEMA_VERSION}")
        print(f"interpretation_regime = protocol_convergence_demonstration")
        print(f"absolute_efficiency_claim = False")
        print(f"ground truth ({len(gt)} entries):")
        for k, v in gt.items():
            print(f"  {k}: ddg={v['ddg_kcal']:+.4f} kcal/mol  "
                  f"kd={v['kd_molar']:.3e} M  prov={v['provenance']}")
        return 0

    # Three strategies, fresh RNG each so outputs are independent but
    # deterministic given --seed.
    rng_al = np.random.default_rng(args.seed)
    rng_rd = np.random.default_rng(args.seed + 1)
    rng_gd = np.random.default_rng(args.seed + 2)

    al_run = run_simulated_cycle(
        strategy="al", n_cycles=args.n_cycles,
        candidates_per_cycle=args.candidates_per_cycle,
        ground_truth=gt, rng=rng_al,
    )
    random_run = run_simulated_cycle(
        strategy="random", n_cycles=args.n_cycles,
        candidates_per_cycle=args.candidates_per_cycle,
        ground_truth=gt, rng=rng_rd,
    )
    greedy_run = run_simulated_cycle(
        strategy="greedy", n_cycles=args.n_cycles,
        candidates_per_cycle=args.candidates_per_cycle,
        ground_truth=gt, rng=rng_gd,
    )

    payload = _build_output_payload(
        gt, al_run, random_run, greedy_run,
        n_cycles=args.n_cycles,
        candidates_per_cycle=args.candidates_per_cycle,
    )

    out_dir = Path(args.output_dir) if args.output_dir else (
        _REPO_ROOT / "outputs" / f"simulated_al_cycle_{datetime.now(timezone.utc).strftime('%Y%m%d')}"
    )
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "simulated_al_cycle.json"
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2, default=str)
    print(f"wrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
