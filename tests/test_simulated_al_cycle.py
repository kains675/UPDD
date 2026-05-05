"""tests/test_simulated_al_cycle.py - SciVal C6 (>=8 functions) coverage for
``scripts/simulated_al_cycle.py`` (schema ``simulated_al_cycle/0.1``).

Test layout (per SciVal verdict 2026-05-05 conditions):
    1. test_synthetic_kd_generation_determinism   - same seed -> same output
    2. test_cp4_anchor_is_magotti_ssot            - C1 BLOCKING enforcement
    3. test_acquisition_ucb_correctness           - UCB analytical case
    4. test_acquisition_ei_correctness            - EI analytical case
    5. test_acquisition_max_entropy_correctness   - MaxEntropy correctness
    6. test_calibration_linear_residual           - sklearn parity
    7. test_calibration_gpr_transition            - n_cumulative + cycle gate
    8. test_baseline_random_greedy_comparators    - schema validity
    9. test_plateau_triggering_monotone_sequence  - monotone-delta plateau
    10.test_output_schema_regression              - schema_version stable
"""

from __future__ import annotations

import importlib.util as _ilu
import math
import sys
from pathlib import Path

import numpy as np
import pytest

# Mirror conftest pattern - put repo root + utils on sys.path so the script's
# lazy imports work, then load scripts/simulated_al_cycle.py by filespec.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

_SCRIPT_PATH = _REPO_ROOT / "scripts" / "simulated_al_cycle.py"
_spec = _ilu.spec_from_file_location("scripts_simulated_al_cycle", _SCRIPT_PATH)
_mod = _ilu.module_from_spec(_spec)  # type: ignore[arg-type]
sys.modules["scripts_simulated_al_cycle"] = _mod
assert _spec is not None and _spec.loader is not None
_spec.loader.exec_module(_mod)  # type: ignore[union-attr]
sal = _mod


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def ground_truth():
    return sal.build_synthetic_ground_truth()


# ---------------------------------------------------------------------------
# Test 1: determinism under seed
# ---------------------------------------------------------------------------

def test_synthetic_kd_generation_determinism(ground_truth):
    """Same seed -> identical per-cycle picks + RMSE trace."""
    rng_a = np.random.default_rng(42)
    rng_b = np.random.default_rng(42)
    out_a = sal.run_simulated_cycle(
        strategy="al", n_cycles=3, candidates_per_cycle=5,
        ground_truth=ground_truth, rng=rng_a,
    )
    out_b = sal.run_simulated_cycle(
        strategy="al", n_cycles=3, candidates_per_cycle=5,
        ground_truth=ground_truth, rng=rng_b,
    )
    a_picks = [tuple(c["selected"]) for c in out_a["per_cycle"]]
    b_picks = [tuple(c["selected"]) for c in out_b["per_cycle"]]
    assert a_picks == b_picks
    assert out_a["calibration"]["rmse_per_cycle_kcal"] == \
        out_b["calibration"]["rmse_per_cycle_kcal"]


# ---------------------------------------------------------------------------
# Test 2: C1 BLOCKING - Cp4 from Magotti SSOT, no 2QKI_Cp4_calib lookup
# ---------------------------------------------------------------------------

def test_cp4_anchor_is_magotti_ssot(ground_truth):
    """Hardcoded constant + no 2QKI_Cp4_calib in Cp4 path."""
    assert sal.MAGOTTI_2009_CP4_DDG_KCAL == -1.4
    assert sal.MAGOTTI_2009_CP4_DDG_KCAL_SPR == -3.0
    cp4 = ground_truth["Cp4"]
    assert cp4["ddg_kcal"] == -1.4
    assert cp4["provenance"] == "Magotti_2009_ITC_literature_anchor"
    assert "10.1002/jmr.972" in cp4["source"]
    # Anti-pattern guard: no 2QKI_Cp4_calib literal anywhere in script.
    src = (_REPO_ROOT / "scripts" / "simulated_al_cycle.py").read_text()
    assert "2QKI_Cp4_calib" not in src
    assert "Cp4_calib" not in src


# ---------------------------------------------------------------------------
# Test 3: UCB on a 3-candidate analytical case
# ---------------------------------------------------------------------------

def test_acquisition_ucb_correctness():
    mu = np.array([1.0, 2.0, 3.0])
    sigma = np.array([0.5, 0.5, 0.5])
    scores = sal.acquisition_ucb(mu, sigma, beta=1.0)
    # mu + beta*sigma = [1.5, 2.5, 3.5]
    assert np.allclose(scores, np.array([1.5, 2.5, 3.5]))
    assert int(np.argmax(scores)) == 2


# ---------------------------------------------------------------------------
# Test 4: EI on a known closed-form normal-posterior case
# ---------------------------------------------------------------------------

def test_acquisition_ei_correctness():
    """For mu = best, sigma > 0, EI reduces to sigma * phi(0) = sigma / sqrt(2 pi).

    Also verifies EI is monotone increasing in mu when mu > best with sigma fixed.
    """
    mu = np.array([0.0, 1.0, 2.0])
    sigma = np.array([1.0, 1.0, 1.0])
    best = 0.0
    ei = sal.acquisition_ei(mu, sigma, current_best=best)
    # At mu=best=0, sigma=1: EI = 1 * (1/sqrt(2 pi)) since z=0, Phi(0)=0.5,
    # phi(0)=1/sqrt(2 pi); contribution (mu-best)*Phi = 0 -> EI = 1/sqrt(2 pi).
    assert math.isclose(ei[0], 1.0 / math.sqrt(2.0 * math.pi), rel_tol=1e-9)
    # Monotone increasing in mu when sigma fixed and mu > best.
    assert ei[1] < ei[2]
    assert ei[0] < ei[1]


# ---------------------------------------------------------------------------
# Test 5: MaxEntropy reduces to log(sigma) ranking
# ---------------------------------------------------------------------------

def test_acquisition_max_entropy_correctness():
    mu = np.array([0.0, 0.0, 0.0])
    sigma = np.array([0.1, 1.0, 5.0])
    ent = sal.acquisition_max_entropy(mu, sigma)
    # Entropy of N(0, sigma^2) is monotone in sigma.
    assert ent[0] < ent[1] < ent[2]
    # Closed form: 0.5 log(2 pi e sigma^2).
    expected = 0.5 * np.log(2.0 * math.pi * math.e * sigma ** 2)
    assert np.allclose(ent, expected)


# ---------------------------------------------------------------------------
# Test 6: linear-residual fit parity with sklearn LinearRegression
# ---------------------------------------------------------------------------

def test_calibration_linear_residual():
    from sklearn.linear_model import LinearRegression
    rng = np.random.default_rng(0)
    x = rng.normal(size=20)
    y = 1.7 * x + 0.3 + rng.normal(scale=0.05, size=20)
    a, b = sal.fit_linear_residual(x, y)
    lr = LinearRegression().fit(x.reshape(-1, 1), y)
    assert math.isclose(a, float(lr.coef_[0]), abs_tol=1e-6)
    assert math.isclose(b, float(lr.intercept_), abs_tol=1e-6)


# ---------------------------------------------------------------------------
# Test 7: GPR transition criterion (cycle >= 3 AND n_cumulative >= 15)
# ---------------------------------------------------------------------------

def test_calibration_gpr_transition():
    # Cycle < 3: never transition, regardless of n_cumulative.
    assert sal.should_transition_to_gpr(cycle=1, n_cumulative=100) is False
    assert sal.should_transition_to_gpr(cycle=2, n_cumulative=100) is False
    # n_cumulative < 15: never transition, regardless of cycle.
    assert sal.should_transition_to_gpr(cycle=10, n_cumulative=14) is False
    # Both conditions met -> transition.
    assert sal.should_transition_to_gpr(cycle=3, n_cumulative=15) is True
    assert sal.should_transition_to_gpr(cycle=4, n_cumulative=20) is True


# ---------------------------------------------------------------------------
# Test 8: baseline random + greedy emit valid output schemas
# ---------------------------------------------------------------------------

def test_baseline_random_greedy_comparators(ground_truth):
    rng = np.random.default_rng(7)
    rnd_out = sal.run_simulated_cycle(
        strategy="random", n_cycles=3, candidates_per_cycle=3,
        ground_truth=ground_truth, rng=rng,
    )
    grd_rng = np.random.default_rng(8)
    grd_out = sal.run_simulated_cycle(
        strategy="greedy", n_cycles=3, candidates_per_cycle=3,
        ground_truth=ground_truth, rng=grd_rng,
    )
    for out in (rnd_out, grd_out):
        assert "per_cycle" in out
        assert "calibration" in out
        assert "plateau_detection" in out
        assert len(out["per_cycle"]) == 3
        for entry in out["per_cycle"]:
            for k in ("cycle", "acquisition", "selected", "predicted_kcal",
                      "ground_truth_kcal", "calibration_error_kcal", "n_cumulative"):
                assert k in entry
        assert "rmse_per_cycle_kcal" in out["calibration"]
        assert "coverage_probability_per_cycle" in out["calibration"]
        assert len(out["calibration"]["rmse_per_cycle_kcal"]) == 3
    # Direct comparator helpers.
    cands = ["a", "b", "c", "d", "e"]
    rng2 = np.random.default_rng(0)
    pick_rnd = sal.baseline_random(cands, batch_size=2, rng=rng2)
    assert len(pick_rnd) == 2
    pick_gr = sal.baseline_greedy(cands, scores=np.array([1, 5, 3, 4, 2]),
                                  batch_size=2)
    # Greedy descending: indices 1 (5) and 3 (4) -> ['b', 'd'].
    assert pick_gr == ["b", "d"]


# ---------------------------------------------------------------------------
# Test 9: plateau detection on monotone-improving DG sequence (-2, -2.4, -2.8)
# ---------------------------------------------------------------------------

def test_plateau_triggering_monotone_sequence():
    """detect_plateau must fire when last 2 abs(delta) < 0.5 kcal."""
    from utils import cycle_manager as cm
    # Each step DDG ~ 0.4 kcal, all < 0.5 threshold (lookback default 3
    # requires lookback+1 = 4 wet-lab cycles; we instead use lookback=2 and
    # 3 wet-lab cycles to confirm the underlying mechanism).
    history = cm.CycleHistory(schema=cm.SCHEMA_VERSION, project_id="t",
                              target_pdb="t.pdb", created_at="2026-05-05")
    for i, dg in enumerate([-2.0, -2.4, -2.8]):
        kd = math.exp(dg / (cm.R_KCAL_PER_MOL_K * cm.DEFAULT_TEMPERATURE_K))
        ent = cm.CycleEntry(cycle_n=i, wt_sequence="X", variants=["X"])
        ent.wet_lab_results.append(
            cm.WetLabResult(sequence="X", kd_molar=kd, assay="ITC")
        )
        ent.best_validated_sequence = "X"
        history.cycles.append(ent)
    assert cm.detect_plateau(history, lookback=2, threshold_kcal=0.5) is True
    # And does NOT fire with the default 0.3 kcal threshold (deltas ~ 0.4):
    assert cm.detect_plateau(history, lookback=2, threshold_kcal=0.3) is False


# ---------------------------------------------------------------------------
# Test 10: output schema regression (top-level fields stable)
# ---------------------------------------------------------------------------

def test_output_schema_regression(ground_truth):
    rng_al = np.random.default_rng(11)
    rng_rd = np.random.default_rng(12)
    rng_gd = np.random.default_rng(13)
    al = sal.run_simulated_cycle(strategy="al", n_cycles=3,
                                 candidates_per_cycle=3,
                                 ground_truth=ground_truth, rng=rng_al)
    rd = sal.run_simulated_cycle(strategy="random", n_cycles=3,
                                 candidates_per_cycle=3,
                                 ground_truth=ground_truth, rng=rng_rd)
    gd = sal.run_simulated_cycle(strategy="greedy", n_cycles=3,
                                 candidates_per_cycle=3,
                                 ground_truth=ground_truth, rng=rng_gd)
    payload = sal._build_output_payload(ground_truth, al, rd, gd,
                                        n_cycles=3, candidates_per_cycle=3)
    assert payload["schema_version"] == "simulated_al_cycle/0.1"
    # C2 framing - mandatory fields.
    assert payload["interpretation_regime"] == "protocol_convergence_demonstration"
    assert payload["absolute_efficiency_claim"] is False
    assert payload["prospective_wet_lab_required"] is True
    assert "Cp4 anchor: Magotti 2009 ITC -1.4" in payload["notes_provenance_disclosure"]
    # C4 sigma fix-4 disclosure flag.
    assert payload["notes"]["sigma_fix4_target_selectivity"] is True
    # C3 coverage probability surfaced for each strategy.
    for strat in ("al", "random", "greedy"):
        cov = payload["strategies"][strat]["calibration"]["coverage_probability_per_cycle"]
        assert isinstance(cov, list) and len(cov) == 3
    # C1 - ground truth is propagated.
    assert payload["ground_truth"]["Cp4"]["ddg_kcal"] == -1.4
    # C7 - single-anchor framing fields.
    assert payload["ground_truth_source"] == "literature_single_anchor"
    assert payload["anchor_count"] == 1
    assert payload["comparator_source"] == "pipeline_self_prediction"
    assert payload["demonstration_scope"] == "single_anchor_convergence"
    # C7 - ground_truth holds ONLY literature anchor entries.
    for k, v in payload["ground_truth"].items():
        assert "literature_anchor" in v["provenance"], (
            f"ground_truth must hold only literature anchors; "
            f"{k} provenance={v['provenance']}"
        )
    # C7 - comparators key holds the 5 pipeline-self-prediction pairs.
    assert "comparators" in payload, "C7: comparators field mandatory"
    pipeline_pairs = {"1EBP_MTR13", "1YCR_NML22", "3IOL_NML20",
                     "7TL8_MTR6", "2QKH_MTR25"}
    assert pipeline_pairs.issubset(set(payload["comparators"].keys())), (
        f"C7: comparators must contain all 5 pipeline pairs; "
        f"got {set(payload['comparators'].keys())}"
    )
    for k in pipeline_pairs:
        assert "literature_anchor" not in payload["comparators"][k]["provenance"], (
            f"C7: comparator {k} must NOT carry literature_anchor provenance"
        )
    # C7 disclosure text in notes.
    assert "single_anchor_framing_disclosure" in payload["notes"]
    # No prohibited efficiency claims in the narrative.
    nar = payload["comparative_summary"]["narrative"].lower()
    assert "30%" not in nar
    # If "wet-lab try reduction" appears, it must be in the negation context.
    if "wet-lab try reduction" in nar:
        assert "not supported" in nar
