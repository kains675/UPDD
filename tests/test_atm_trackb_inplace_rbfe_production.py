# -*- coding: utf-8 -*-
"""Tests for the in-place residue-4 RBFE per-direction PRODUCTION launcher.

Covers the launcher's pure / logic surface WITHOUT a real OpenMM build or GPU:
  - the C8 SIGN-critical bound-leg decouple-direction gate;
  - the matched-seed availability gate;
  - the C7 Welch-Satterthwaite quadrature SE + z_SE sign rule;
  - the COMBINED + single-direction cntl writers (UWHAM-consumable / SSOT);
  - the C11 pre-registration of the 4 outcomes;
  - the soft-core usc helper (the .out pertE must be the CAPPED perturbation,
    not the raw u1-u0 that would overflow the WHAM solve).

Uses the importlib spec-loader pattern (no sys.modules pollution; the launcher +
bridge are openmm-only and the qmmm env has no atom_openmm).
"""

import importlib.util
import math
import os
import sys

import pytest


_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ = os.path.dirname(_HERE)
_UTILS = os.path.join(_PROJ, "utils")
_SCRIPTS = os.path.join(_PROJ, "scripts")


def _load(name, relpath):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_PROJ, relpath))
    if spec is None or spec.loader is None:
        pytest.skip("could not locate %s" % relpath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def prod():
    return _load("trackb_inplace_rbfe_production",
                 "scripts/trackb_inplace_rbfe_production.py")


@pytest.fixture(scope="module")
def rbfe():
    if _UTILS not in sys.path:
        sys.path.insert(0, _UTILS)
    return _load("atm_trackB_inplace_rbfe", "utils/atm_trackB_inplace_rbfe.py")


# --------------------------- C8 decouple-direction gate --------------------
def test_c8_free_leg_none_is_ok(prod):
    r = prod.gate_decouple_direction("free", None, raise_on_fail=True)
    assert r["passed"] is True
    assert "free_leg" in r["reason"]


def test_c8_bound_leg_none_fails_loud(prod):
    with pytest.raises(RuntimeError):
        prod.gate_decouple_direction("bound", None, raise_on_fail=True)
    # non-raising form reports the violation.
    r = prod.gate_decouple_direction("bound", None, raise_on_fail=False)
    assert r["passed"] is False and "C8 VIOLATION" in r["reason"]


def test_c8_bound_leg_unit_vector_passes(prod):
    # An outward unit vector (the kind compute_decouple_direction returns).
    v = (0.0, 0.0, 1.0)
    r = prod.gate_decouple_direction("bound", v, raise_on_fail=True)
    assert r["passed"] is True
    assert abs(r["magnitude"] - 1.0) < 1e-6


def test_c8_bound_leg_non_unit_fails(prod):
    with pytest.raises(RuntimeError):
        prod.gate_decouple_direction("bound", (0.0, 0.0, 2.0), raise_on_fail=True)


def test_c8_bound_leg_non_finite_fails(prod):
    with pytest.raises(RuntimeError):
        prod.gate_decouple_direction("bound", (float("nan"), 0.0, 1.0),
                                     raise_on_fail=True)


def test_c8_bound_leg_wrong_length_fails(prod):
    with pytest.raises(RuntimeError):
        prod.gate_decouple_direction("bound", (0.0, 1.0), raise_on_fail=True)


# --------------------------- seed-availability gate ------------------------
def test_seed_gate_default_seeds_present(prod):
    # The default matched seeds (s7, s101, s127) must exist for both endpoints
    # in this repo — that is the basis of the paired ddG. If this fails the
    # cohort drifted and the default must be updated (a real, not silent, gap).
    r = prod.gate_seed_availability(list(prod.ENDPOINTS), list(prod.DEFAULT_SEEDS))
    assert r["passed"] is True, r["missing"]


def test_seed_gate_missing_seed_fails(prod):
    r = prod.gate_seed_availability(["cp4", "wt"], ["s_does_not_exist_999"])
    assert r["passed"] is False
    assert len(r["missing"]) == 2  # both endpoints missing this seed


# --------------------------- C7 quadrature / z_SE --------------------------
def test_quadrature_se_and_z(prod):
    q = prod._welch_satterthwaite_quadrature(10.0, 1.0, 7.0, 1.0)
    assert abs(q["ddg"] - 3.0) < 1e-12
    assert abs(q["se_quadrature"] - math.sqrt(2.0)) < 1e-9
    assert abs(q["z_se"] - 3.0 / math.sqrt(2.0)) < 1e-9


def test_quadrature_zero_se_returns_none_z(prod):
    q = prod._welch_satterthwaite_quadrature(1.0, 0.0, 1.0, 0.0)
    assert q["z_se"] is None


def test_z_se_threshold_is_three(prod):
    # The pre-registration sign rule: z_SE >= 3 AND n >= 3 (R-18 honest).
    assert prod.Z_SE_SIGN_THRESHOLD == 3.0
    assert prod.MIN_REPLICATES_FOR_SIGN == 3


# --------------------------- cntl writers ----------------------------------
def test_combined_cntl_has_both_directions(prod, rbfe, tmp_path):
    sched = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2,
                                   single_direction=None)
    cntl = str(tmp_path / "trackb_asyncre.cntl")
    prod._write_combined_cntl(cntl, sched, "free", 250, 1.0)
    text = open(cntl).read()
    assert "MODE = 'INPLACE_RBFE'" in text
    # The combined cntl is the SSOT: DIRECTION must carry BOTH a +1 and a -1
    # block (UWHAM requires both legs); LIGAND_ATOMS / DISPLACEMENT must NOT
    # appear (that is the ABFE whole-binder path, not the in-place swap).
    dir_line = [ln for ln in text.splitlines() if ln.startswith("DIRECTION")][0]
    assert "1" in dir_line and "-1" in dir_line
    assert "LIGAND_ATOMS" not in text and "DISPLACEMENT" not in text


def test_direction_cntl_single_direction_only(prod, rbfe, tmp_path):
    sched = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2,
                                   single_direction="forward")
    cntl = str(tmp_path / "trackb_dplus_asyncre.cntl")
    prod._write_direction_cntl(cntl, sched, "dplus", 250, 1.0)
    text = open(cntl).read()
    assert "BASENAME = 'trackb_dplus'" in text
    dir_line = [ln for ln in text.splitlines() if ln.startswith("DIRECTION")][0]
    # forward standalone -> only +1, no -1.
    assert "-1" not in dir_line


# --------------------------- C11 pre-registration --------------------------
def test_pre_registration_has_four_outcomes(prod, tmp_path):
    path = prod.write_pre_registration(str(tmp_path), {"leg": "free"})
    import json
    payload = json.load(open(path))
    ids = [o["id"] for o in payload["outcomes"]]
    assert len(ids) == 4
    assert "O1_favorable_sign_resolved" in ids
    assert "O3_sign_undetermined" in ids
    assert "O4_non_convergent" in ids
    # The sign-claim policy must state z_SE >= 3 AND n >= 3 (no silent assume).
    assert "z_SE >= 3" in payload["sign_claim_policy"]


# --------------------------- soft-core usc helper --------------------------
def test_softcore_usc_caps_large_perturbation(rbfe):
    # A raw u1-u0 of ~50,000 kcal/mol (the single-shared-core clash) must be
    # capped to < umax. The .out pertE column carries usc, NOT the raw value —
    # otherwise UWHAM's _bias_fcn (no cap of its own) overflows the WHAM solve.
    kcal = 4.184
    umax_kj = 200.0 * kcal
    ubcore_kj = 100.0 * kcal
    acore = 0.0625
    big_pert_kj = 50000.0 * kcal
    base, usc = rbfe._atm_softcore_components_kj(
        big_pert_kj, -17000.0 * kcal, 1.0, umax_kj, ubcore_kj, acore)
    usc_kcal = usc / kcal
    assert usc_kcal < 200.0, usc_kcal
    assert usc_kcal > 100.0  # above ubcore (it was a large perturbation)


def test_softcore_usc_small_perturbation_passthrough(rbfe):
    kcal = 4.184
    umax_kj = 200.0 * kcal
    ubcore_kj = 100.0 * kcal
    acore = 0.0625
    # A small (< ubcore) perturbation passes through uncapped.
    small_kj = 5.0 * kcal
    base, usc = rbfe._atm_softcore_components_kj(
        small_kj, 0.0, 1.0, umax_kj, ubcore_kj, acore)
    assert abs(usc / kcal - 5.0) < 1e-9


def test_softcore_usc_backward_sign(rbfe):
    # Backward (Direction<0) flips the sign of the perturbation + bases on u1.
    kcal = 4.184
    umax_kj = 200.0 * kcal
    ubcore_kj = 100.0 * kcal
    base_f, usc_f = rbfe._atm_softcore_components_kj(
        10.0 * kcal, 0.0, 1.0, umax_kj, ubcore_kj, 0.0625)
    base_b, usc_b = rbfe._atm_softcore_components_kj(
        10.0 * kcal, 0.0, -1.0, umax_kj, ubcore_kj, 0.0625)
    # forward base = u0 = 0; backward base = u1 = u0 + pert.
    assert abs(base_f - 0.0) < 1e-9
    assert abs(base_b - 10.0 * kcal) < 1e-6
    # usc sign flips.
    assert usc_f > 0 and usc_b < 0


def test_state_energy_uses_softcore_components(rbfe):
    # _atm_state_energy_kj must equal base + bias_fcn(usc) so the .out potE +
    # capped pertE reconstruct in UWHAM (e0 = potE - bias_fcn(usc) = base).
    kcal = 4.184
    umax_kj = 200.0 * kcal
    ubcore_kj = 100.0 * kcal
    acore = 0.0625
    pert_kj = 8.0 * kcal
    u0_kj = -100.0 * kcal
    lambda1, lambda2 = 0.3, 0.4
    alpha_per_kj = 0.10 / kcal
    uh_kj = 30.0 * kcal
    w0_kj = 0.0
    base, usc = rbfe._atm_softcore_components_kj(
        pert_kj, u0_kj, 1.0, umax_kj, ubcore_kj, acore)
    U = rbfe._atm_state_energy_kj(
        pert_kj, u0_kj, lambda1, lambda2, alpha_per_kj, uh_kj, w0_kj,
        1.0, umax_kj, ubcore_kj, acore)
    # Reconstruct bias_fcn(usc) in kJ and check U == base + bias.
    bias = (((lambda2 - lambda1) / alpha_per_kj)
            * math.log(1.0 + math.exp(-alpha_per_kj * (usc - uh_kj)))
            + lambda2 * usc + w0_kj)
    assert abs(U - (base + bias)) < 1e-6


# --------------------------- parse helpers ---------------------------------
def test_parse_directions_default_both(prod):
    assert prod._parse_directions(None) == ["dplus", "dminus"]
    assert prod._parse_directions("dplus") == ["dplus"]
    with pytest.raises(ValueError):
        prod._parse_directions("dbogus")


def test_parse_seeds_default_and_custom(prod):
    assert prod._parse_seeds(None) == list(prod.DEFAULT_SEEDS)
    assert prod._parse_seeds("sA,sB") == ["sA", "sB"]
    with pytest.raises(ValueError):
        prod._parse_seeds(",")


def test_max_concurrent_units_packs_to_vram(prod):
    # 14 GiB free, 2 GiB/unit, 1.5 GiB headroom -> floor(12.5/2) = 6 units.
    assert prod._max_concurrent_units(14.0, 2.0) == 6
    # bound: 30 GiB free, 6.5 GiB/unit -> floor(28.5/6.5) = 4.
    assert prod._max_concurrent_units(30.0, 6.5) == 4
    # never below 1.
    assert prod._max_concurrent_units(1.0, 6.5) == 1
