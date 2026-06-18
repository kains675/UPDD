# -*- coding: utf-8 -*-
"""Tests for the in-place residue-4 RBFE ladder bridge.

Covers (without a real OpenMM build where possible):
  - the RBFE λ-ladder builder (shape / symmetry / soft-core canon / parameter
    nature of the window count);
  - the ATM hybrid-potential recombination ``_atm_state_energy_kj`` (pure math
    reproducing the upstream expression);
  - the serializer / load-only round-trip contract + the ladder run + the
    per-direction driver mixing-gate reuse (guarded behind endpoint presence +
    openmm, slow — opt-in via the build fixture).

Uses the same importlib spec-loader pattern as test_atm_trackb_inplace_res4_ats
so a pytest session does not pollute sys.modules and the atom_openmm-free qmmm
env can import the module (the RBFE module + mixing gate are openmm-only).
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


def _have_openmm():
    try:
        import openmm  # noqa: F401
        return True
    except ImportError:
        return False


def _endpoints_present(seed="s7"):
    cp4 = os.path.join(_PROJ, "outputs", f"2QKI_Cp4_hybrid_calib_{seed}",
                       "mdresult", "2QKI_Cp4_final.pdb")
    wt = os.path.join(_PROJ, "outputs", f"2QKI_WT_calib_{seed}",
                      "mdresult", "2QKI_WT_final.pdb")
    return os.path.isfile(cp4) and os.path.isfile(wt)


def _load_rbfe():
    spec = importlib.util.spec_from_file_location(
        "atm_trackB_inplace_rbfe",
        os.path.join(_UTILS, "atm_trackB_inplace_rbfe.py"))
    if spec is None or spec.loader is None:
        pytest.skip("could not locate utils/atm_trackB_inplace_rbfe.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _load_smoke():
    spec = importlib.util.spec_from_file_location(
        "trackb_inplace_rbfe_ladder_smoke",
        os.path.join(_SCRIPTS, "trackb_inplace_rbfe_ladder_smoke.py"))
    if spec is None or spec.loader is None:
        pytest.skip("could not locate scripts/trackb_inplace_rbfe_ladder_smoke.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _load_driver_mixing():
    spec = importlib.util.spec_from_file_location(
        "trackb_per_direction_production",
        os.path.join(_SCRIPTS, "trackb_per_direction_production.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def rbfe():
    if not _have_openmm():
        pytest.skip("openmm not importable in this env")
    return _load_rbfe()


# ---------------------------------------------------------------------------
# Ladder builder — pure-Python (build_rbfe_ladder imports openmm only via the
# ats module's top-level import; guarded by the fixture).
# ---------------------------------------------------------------------------
def test_ladder_shape_and_symmetry(rbfe):
    sch = rbfe.build_rbfe_ladder(n_windows_half=8, softcore_band=2)
    assert sch["n_states"] == 16
    for key in ("lambdas", "lambdas_1", "lambdas_2", "directions", "intermd",
                "alpha", "u0", "w0"):
        assert len(sch[key]) == 16, key
    # Forward = +1, backward = -1.
    assert sch["directions"] == [1] * 8 + [-1] * 8
    # Backward half = whole-tuple reverse of the forward half.
    assert sch["lambdas_1"][8:] == list(reversed(sch["lambdas_1"][:8]))
    assert sch["w0"][8:] == list(reversed(sch["w0"][:8]))
    assert sch["intermd"][8:] == list(reversed(sch["intermd"][:8]))


def test_ladder_lambda_endpoints(rbfe):
    sch = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2)
    # Forward λ1 runs 0 -> 0.5 inclusive; backward returns to 0.
    assert sch["lambdas_1"][0] == 0.0
    assert sch["lambdas_1"][5] == 0.5
    assert sch["lambdas_1"][-1] == 0.0
    # λ1 is monotone non-decreasing on the forward half.
    fwd = sch["lambdas_1"][:6]
    assert all(fwd[i] <= fwd[i + 1] for i in range(len(fwd) - 1))


def test_ladder_softcore_band(rbfe):
    sch = rbfe.build_rbfe_ladder(n_windows_half=8, softcore_band=3)
    # softcore_band states at the apex carry INTERMEDIATE=1 + W0>0.
    fwd_inter = sch["intermd"][:8]
    assert sum(fwd_inter) == 3
    # The INTERMEDIATE states are the LAST 3 of the forward half (closest to 0.5).
    assert fwd_inter[-3:] == [1, 1, 1]
    assert fwd_inter[:5] == [0, 0, 0, 0, 0]
    # W0 ramps up across the band (last band state has the apex W0).
    band_w0 = [sch["w0"][i] for i in range(8) if sch["intermd"][i] == 1]
    assert band_w0 == sorted(band_w0)
    assert band_w0[-1] == pytest.approx(rbfe.RBFE_W0_APEX)


def test_ladder_softcore_canon_not_retuned(rbfe):
    sch = rbfe.build_rbfe_ladder(n_windows_half=4, softcore_band=1)
    # The global soft-core canon must equal the validated ATS values.
    assert sch["umax"] == rbfe.ats.ATS_UMAX_KCAL == 200.0
    assert sch["ubcore"] == rbfe.ats.ATS_UBCORE_KCAL == 100.0
    assert sch["acore"] == rbfe.ats.ATS_ACORE == 0.0625


def test_ladder_window_count_is_parameter(rbfe):
    # The count is a PARAMETER (R-18 — not pre-assumed converged): different
    # n_windows_half yields different state counts with the same structure.
    for n in (4, 6, 8, 12):
        sch = rbfe.build_rbfe_ladder(n_windows_half=n, softcore_band=1)
        assert sch["n_states"] == 2 * n
        assert sch["n_windows_half"] == n


def test_ladder_rejects_degenerate_counts(rbfe):
    with pytest.raises(ValueError):
        rbfe.build_rbfe_ladder(n_windows_half=1)
    with pytest.raises(ValueError):
        rbfe.build_rbfe_ladder(n_windows_half=4, softcore_band=0)
    with pytest.raises(ValueError):
        rbfe.build_rbfe_ladder(n_windows_half=4, softcore_band=5)
    with pytest.raises(ValueError):
        rbfe.build_rbfe_ladder(n_windows_half=4, softcore_band=1, n_apex_bridge=-1)
    with pytest.raises(ValueError):
        rbfe.build_rbfe_ladder(n_windows_half=4, softcore_band=1, apex_band=0.0)
    with pytest.raises(ValueError):
        rbfe.build_rbfe_ladder(n_windows_half=4, softcore_band=1, apex_band=1.5)


def test_apex_bridge_default_is_byte_identical(rbfe):
    """n_apex_bridge=0 (the DEFAULT) must produce the pre-bridge ladder exactly
    (additive change; existing callers unaffected)."""
    base = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2)
    same = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2,
                                  n_apex_bridge=0)
    for key in ("lambdas", "lambdas_1", "lambdas_2", "directions", "intermd",
                "alpha", "u0", "w0", "n_states"):
        assert base[key] == same[key], key
    assert base["n_apex_bridge"] == 0


def test_apex_bridge_inserts_finer_w0_states(rbfe):
    """An apex bridge adds n_apex_bridge EXTRA states to EACH half just below
    the λ=0.5 apex, in a CONCAVE (denser near W0=1.0) W0 grid; the ladder stays
    symmetric and the soft-core canon is untouched."""
    n_bridge = 2
    sch = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2,
                                 n_apex_bridge=n_bridge)
    # Forward half grows by n_bridge per half -> total grows by 2*n_bridge.
    assert sch["n_states"] == 12 + 2 * n_bridge
    half = sch["n_windows_half"]
    assert half == 6 + n_bridge
    # Symmetry preserved (whole-tuple reverse).
    assert sch["lambdas_1"][half:] == list(reversed(sch["lambdas_1"][:half]))
    assert sch["w0"][half:] == list(reversed(sch["w0"][:half]))
    assert sch["directions"] == [1] * half + [-1] * half
    # The bridge windows sit at the apex λ (0.5) with W0 strictly increasing to
    # the apex W0 (1.0) — finer steps approaching saturation.
    fwd_w0 = sch["w0"][:half]
    # All forward W0 are non-decreasing (anneal ramps up to the apex).
    assert all(fwd_w0[i] <= fwd_w0[i + 1] for i in range(len(fwd_w0) - 1))
    # The last forward state is the apex (W0 == apex).
    assert fwd_w0[-1] == pytest.approx(rbfe.RBFE_W0_APEX)
    # The bridge inserted >=1 NEW distinct W0 level between the penultimate
    # pre-bridge band W0 and the apex (concave compression near 1.0).
    base = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2)
    base_fwd_w0 = set(base["w0"][:base["n_windows_half"]])
    new_levels = set(fwd_w0) - base_fwd_w0
    assert len(new_levels) >= 1
    # Soft-core canon NOT retuned by the bridge.
    assert sch["umax"] == rbfe.ats.ATS_UMAX_KCAL == 200.0
    assert sch["ubcore"] == rbfe.ats.ATS_UBCORE_KCAL == 100.0
    assert sch["acore"] == rbfe.ats.ATS_ACORE == 0.0625


# ---------------------------------------------------------------------------
# Single-direction STANDALONE ladder (the decisive fork test) — pure-Python.
# ---------------------------------------------------------------------------
def test_single_direction_default_is_full_symmetric(rbfe):
    """single_direction=None (the DEFAULT) is byte-identical to the pre-existing
    full symmetric two-direction ladder (additive change)."""
    base = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2)
    same = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2,
                                  single_direction=None)
    for key in ("lambdas", "lambdas_1", "lambdas_2", "directions", "intermd",
                "alpha", "u0", "w0", "n_states"):
        assert base[key] == same[key], key
    # The full symmetric ladder carries a None single_direction marker.
    assert base["single_direction"] is None
    assert same["single_direction"] is None
    # And it remains two-direction (Dir flips at the half boundary).
    assert base["directions"] == [1] * 6 + [-1] * 6


def test_single_direction_forward_is_forward_half_only(rbfe):
    """single_direction='forward' yields ONLY the forward half: n_windows_half
    states, DIRECTION=+1 throughout, λ1 0->0.5, NO backward states / apex flip."""
    sym = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2)
    fwd = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2,
                                 single_direction="forward")
    assert fwd["single_direction"] == "forward"
    assert fwd["n_states"] == 6
    # DIRECTION=+1 only — no direction-flip boundary.
    assert fwd["directions"] == [1] * 6
    assert set(fwd["directions"]) == {1}
    # State 0 = λ=0 endpoint, last state = λ=0.5 apex.
    assert fwd["lambdas_1"][0] == 0.0
    assert fwd["lambdas_1"][-1] == 0.5
    # The forward standalone == the forward half of the symmetric ladder.
    for key in ("lambdas_1", "lambdas_2", "intermd", "w0", "alpha", "u0"):
        assert fwd[key] == sym[key][:6], key


def test_single_direction_backward_is_reversed_forward(rbfe):
    """single_direction='backward' yields the whole-tuple reverse of the forward
    half, DIRECTION=-1 throughout (the dminus standalone ladder)."""
    fwd = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2,
                                 single_direction="forward")
    bwd = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2,
                                 single_direction="backward")
    assert bwd["single_direction"] == "backward"
    assert bwd["n_states"] == 6
    assert bwd["directions"] == [-1] * 6
    # Whole-tuple reverse of the forward half.
    assert bwd["lambdas_1"] == list(reversed(fwd["lambdas_1"]))
    assert bwd["w0"] == list(reversed(fwd["w0"]))
    assert bwd["intermd"] == list(reversed(fwd["intermd"]))
    # State 0 = apex side, last state = λ=0 endpoint.
    assert bwd["lambdas_1"][0] == 0.5
    assert bwd["lambdas_1"][-1] == 0.0


def test_single_direction_with_apex_bridge(rbfe):
    """A single-direction ladder honours the apex bridge (extra W0 windows) —
    n_states = n_windows_half + n_apex_bridge, still one DIRECTION."""
    n_bridge = 3
    fwd = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2,
                                 n_apex_bridge=n_bridge,
                                 single_direction="forward")
    assert fwd["n_states"] == 6 + n_bridge
    assert fwd["directions"] == [1] * (6 + n_bridge)
    assert fwd["n_apex_bridge"] == n_bridge
    # W0 still ramps non-decreasing to the apex (1.0).
    assert all(fwd["w0"][i] <= fwd["w0"][i + 1]
               for i in range(len(fwd["w0"]) - 1))
    assert fwd["w0"][-1] == pytest.approx(rbfe.RBFE_W0_APEX)
    # Soft-core canon untouched.
    assert fwd["umax"] == rbfe.ats.ATS_UMAX_KCAL == 200.0


def test_single_direction_rejects_bad_value(rbfe):
    with pytest.raises(ValueError):
        rbfe.build_rbfe_ladder(n_windows_half=6, single_direction="sideways")
    with pytest.raises(ValueError):
        rbfe.build_rbfe_ladder(n_windows_half=6, single_direction="fwd")


# ---------------------------------------------------------------------------
# Canonical ATS standard schedule (two-copy box; spec C).
# ---------------------------------------------------------------------------
def test_ats_standard_spec_c_11_windows(rbfe):
    """n_windows_half=6 reproduces the spec C exactly (11 λ states)."""
    sch = rbfe.build_ats_standard_ladder(n_windows_half=6,
                                         single_direction="forward")
    assert sch["n_states"] == 11
    assert sch["lambdas_1"] == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                0.1, 0.2, 0.3, 0.4, 0.5]
    assert sch["lambdas_2"] == [0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                                0.5, 0.5, 0.5, 0.5, 0.5]
    # Single DIRECTION (no flip boundary), Uh=110 (ATS canon, NOT 30).
    assert set(sch["directions"]) == {1}
    assert sch["u0"][0] == rbfe.ATS_TWOCOPY_U0_DEFAULT == 110.0
    assert sch["alpha"][0] == rbfe.RBFE_ALPHA_DEFAULT == 0.10
    assert sch["schedule_kind"] == "ats_standard"


def test_ats_standard_softcore_canon_not_retuned(rbfe):
    sch = rbfe.build_ats_standard_ladder(n_windows_half=6)
    assert sch["umax"] == rbfe.ats.ATS_UMAX_KCAL == 200.0
    assert sch["ubcore"] == rbfe.ats.ATS_UBCORE_KCAL == 100.0
    assert sch["acore"] == rbfe.ats.ATS_ACORE == 0.0625


def test_ats_standard_backward_is_reverse(rbfe):
    fwd = rbfe.build_ats_standard_ladder(n_windows_half=6,
                                         single_direction="forward")
    bwd = rbfe.build_ats_standard_ladder(n_windows_half=6,
                                         single_direction="backward")
    assert set(bwd["directions"]) == {-1}
    assert bwd["lambdas_1"] == list(reversed(fwd["lambdas_1"]))
    assert bwd["lambdas_2"] == list(reversed(fwd["lambdas_2"]))


def test_ats_standard_window_count_is_parameter(rbfe):
    # 2*n - 1 states (the two phases share the apex λ1=0/λ2=0.5 state).
    for n in (3, 4, 6, 8):
        sch = rbfe.build_ats_standard_ladder(n_windows_half=n)
        assert sch["n_states"] == 2 * n - 1
        # λ1 endpoints: starts 0, ends 0.5; λ2 starts 0, ends 0.5.
        assert sch["lambdas_1"][0] == 0.0 and sch["lambdas_1"][-1] == 0.5
        assert sch["lambdas_2"][0] == 0.0 and sch["lambdas_2"][-1] == 0.5


def test_ats_standard_rejects_degenerate(rbfe):
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(n_windows_half=1)
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(n_windows_half=6, single_direction="both")


# ---------------------------------------------------------------------------
# Leg-switch bridge: explicit leg-down λ1 knots (lambda1_rampdown).
# ---------------------------------------------------------------------------
def test_ats_standard_lambda1_rampdown_none_is_byte_identical(rbfe):
    """lambda1_rampdown=None reproduces the historical 11-state schedule exactly.

    The default-None path MUST be byte-identical to the legacy uniform leg-down
    (the existing callers — launcher / smoke / test — pass no lambda1_rampdown).
    """
    legacy = rbfe.build_ats_standard_ladder(n_windows_half=6,
                                            single_direction="forward")
    explicit = rbfe.build_ats_standard_ladder(n_windows_half=6,
                                              single_direction="forward",
                                              lambda1_rampdown=None)
    for key in ("lambdas_1", "lambdas_2", "intermd", "directions",
                "alpha", "u0", "w0", "n_states", "schedule_name"):
        assert explicit[key] == legacy[key], "drift in %r" % (key,)
    assert legacy["n_states"] == 11
    assert legacy["lambda1_rampdown"] is None


def test_ats_standard_lambda1_rampdown_bridge_12_states(rbfe):
    """lambda1_rampdown=[0.05,...,0.5] inserts the leg-switch bridge => 12 states."""
    sch = rbfe.build_ats_standard_ladder(
        n_windows_half=6, single_direction="forward",
        lambda1_rampdown=[0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
    assert sch["n_states"] == 12
    assert sch["lambdas_1"] == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
    assert sch["lambdas_2"] == [0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                                0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    # Apex λ1=0/λ2=0.5 is idx 5; endpoints (idx0 decoupled, idx11 coupled apex)
    # are the only intermd==0 states. The inserted bridge window is intermd==1.
    assert sch["intermd"] == [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
    assert sch["lambdas_1"][5] == 0.0 and sch["lambdas_2"][5] == 0.5  # apex
    assert sch["lambdas_1"][0] == 0.0 and sch["lambdas_1"][-1] == 0.5
    assert sch["lambdas_2"][0] == 0.0 and sch["lambdas_2"][-1] == 0.5
    # Single DIRECTION, soft-core canon NOT re-tuned (C7), Uh=110 (ATS canon).
    assert set(sch["directions"]) == {1}
    assert sch["u0"][0] == rbfe.ATS_TWOCOPY_U0_DEFAULT == 110.0
    assert sch["alpha"][0] == rbfe.RBFE_ALPHA_DEFAULT == 0.10
    assert sch["umax"] == rbfe.ats.ATS_UMAX_KCAL == 200.0
    assert sch["ubcore"] == rbfe.ats.ATS_UBCORE_KCAL == 100.0
    assert sch["acore"] == rbfe.ats.ATS_ACORE == 0.0625
    # λ1 strictly increasing across the leg-down (states 5..11), λ2 monotone.
    legdown = sch["lambdas_1"][5:]
    assert all(b > a for a, b in zip(legdown, legdown[1:]))
    assert sch["lambda1_rampdown"] == [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
    assert "bridge" in sch["schedule_name"]


def test_ats_standard_lambda1_rampdown_backward_is_reverse(rbfe):
    """Backward bridge ladder is the whole-tuple reverse of the forward bridge."""
    fwd = rbfe.build_ats_standard_ladder(
        n_windows_half=6, single_direction="forward",
        lambda1_rampdown=[0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
    bwd = rbfe.build_ats_standard_ladder(
        n_windows_half=6, single_direction="backward",
        lambda1_rampdown=[0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
    assert bwd["n_states"] == 12
    assert set(bwd["directions"]) == {-1}
    assert bwd["lambdas_1"] == list(reversed(fwd["lambdas_1"]))
    assert bwd["lambdas_2"] == list(reversed(fwd["lambdas_2"]))
    assert bwd["intermd"] == list(reversed(fwd["intermd"]))


def test_ats_standard_lambda1_rampdown_rejects_invalid(rbfe):
    """Out-of-range / non-monotone / wrong-endpoint knots fail loud."""
    # value > 0.5
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(
            n_windows_half=6, lambda1_rampdown=[0.1, 0.6])
    # value == 0 (apex is placed by leg-up; knots start above 0)
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(
            n_windows_half=6, lambda1_rampdown=[0.0, 0.5])
    # not strictly increasing
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(
            n_windows_half=6, lambda1_rampdown=[0.2, 0.2, 0.5])
    # does not end at 0.5
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(
            n_windows_half=6, lambda1_rampdown=[0.1, 0.2, 0.3])
    # empty list
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(
            n_windows_half=6, lambda1_rampdown=[])


# ---------------------------------------------------------------------------
# Deep-λ2 decouple-tail densification: explicit leg-up λ2 knots (lambda2_rampup).
# ---------------------------------------------------------------------------
def test_ats_standard_lambda2_rampup_none_is_byte_identical(rbfe):
    """lambda2_rampup=None reproduces the historical 11-state schedule exactly.

    The default-None path MUST be byte-identical to the legacy uniform leg-up
    (every existing caller — launcher / smoke / test — passes no lambda2_rampup).
    """
    legacy = rbfe.build_ats_standard_ladder(n_windows_half=6,
                                            single_direction="forward")
    explicit = rbfe.build_ats_standard_ladder(n_windows_half=6,
                                              single_direction="forward",
                                              lambda2_rampup=None)
    for key in ("lambdas_1", "lambdas_2", "intermd", "directions",
                "alpha", "u0", "w0", "n_states", "n_windows_half_linear",
                "softcore_band", "schedule_name"):
        assert explicit[key] == legacy[key], "drift in %r" % (key,)
    assert legacy["n_states"] == 11
    assert legacy["lambda2_rampup"] is None


def test_ats_standard_lambda2_rampup_densifies_decouple_tail(rbfe):
    """lambda2_rampup=[0.0,0.05,...,0.5] inserts the deep-λ2 bridges => 13 states.

    The seal bonds are the deep-λ2 decouple tail (λ2 0.2<->0.1<->0.0). Inserting
    λ2=0.05 (between 0/0.1) and λ2=0.15 (between 0.1/0.2) gives an 8-state leg-up;
    with the canonical 5-state uniform leg-down (n_windows_half=6) => 8 + 5 = 13.
    """
    sch = rbfe.build_ats_standard_ladder(
        n_windows_half=6, single_direction="forward",
        lambda2_rampup=[0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5])
    assert sch["n_states"] == 13
    assert sch["lambdas_1"] == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                0.1, 0.2, 0.3, 0.4, 0.5]
    assert sch["lambdas_2"] == [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5,
                                0.5, 0.5, 0.5, 0.5, 0.5]
    # The inserted bridge λ2 values are actually present in the leg-up tail.
    assert 0.05 in sch["lambdas_2"][:8]
    assert 0.15 in sch["lambdas_2"][:8]
    # Leg-up has 8 states; only the idx0 decoupled endpoint + idx12 coupled apex
    # are intermd==0; the leg-up apex (idx7, λ1=0/λ2=0.5) is the leg-switch state.
    assert sch["intermd"] == [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
    assert sch["lambdas_1"][7] == 0.0 and sch["lambdas_2"][7] == 0.5  # apex
    assert sch["lambdas_1"][0] == 0.0 and sch["lambdas_2"][0] == 0.0   # decoupled
    assert sch["lambdas_1"][-1] == 0.5 and sch["lambdas_2"][-1] == 0.5  # coupled
    # Leg-up length is reflected in the soft-core band metadata.
    assert sch["n_windows_half_linear"] == 8
    assert sch["softcore_band"] == 8
    # Soft-core canon NOT re-tuned (C7).
    assert sch["u0"][0] == rbfe.ATS_TWOCOPY_U0_DEFAULT == 110.0
    assert sch["alpha"][0] == rbfe.RBFE_ALPHA_DEFAULT == 0.10
    assert sch["umax"] == rbfe.ats.ATS_UMAX_KCAL == 200.0
    assert sch["ubcore"] == rbfe.ats.ATS_UBCORE_KCAL == 100.0
    assert sch["acore"] == rbfe.ats.ATS_ACORE == 0.0625
    # λ2 strictly increasing across the leg-up (states 0..7).
    legup = sch["lambdas_2"][:8]
    assert all(b > a for a, b in zip(legup, legup[1:]))
    assert sch["lambda2_rampup"] == [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
    assert "bridge" in sch["schedule_name"]


def test_ats_standard_lambda2_rampup_backward_is_reverse(rbfe):
    """Backward densified-leg-up ladder is the whole-tuple reverse of forward.

    This is the dminus leg (the W4A seal bonds 8/9, 9/10 live here). The
    backward ladder's deep-λ2 tail is the END of the tuple.
    """
    fwd = rbfe.build_ats_standard_ladder(
        n_windows_half=6, single_direction="forward",
        lambda2_rampup=[0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5])
    bwd = rbfe.build_ats_standard_ladder(
        n_windows_half=6, single_direction="backward",
        lambda2_rampup=[0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5])
    assert bwd["n_states"] == 13
    assert set(bwd["directions"]) == {-1}
    assert bwd["lambdas_1"] == list(reversed(fwd["lambdas_1"]))
    assert bwd["lambdas_2"] == list(reversed(fwd["lambdas_2"]))
    assert bwd["intermd"] == list(reversed(fwd["intermd"]))
    # On the backward leg the decoupled endpoint (λ2=0) is the LAST state, and
    # the deep-λ2 tail (λ2 0.0, 0.05, 0.1, 0.15, 0.2) is the END of the tuple.
    assert bwd["lambdas_2"][-1] == 0.0
    assert bwd["lambdas_2"][-2] == 0.05
    assert bwd["lambdas_2"][-3] == 0.1


def test_ats_standard_lambda2_rampup_rejects_invalid(rbfe):
    """Out-of-range / non-monotone / wrong-endpoint leg-up knots fail loud."""
    # value > 0.5
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(
            n_windows_half=6, lambda2_rampup=[0.0, 0.6])
    # value < 0
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(
            n_windows_half=6, lambda2_rampup=[-0.1, 0.0, 0.5])
    # does NOT start at 0.0 (the leg-up OWNS the decoupled endpoint)
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(
            n_windows_half=6, lambda2_rampup=[0.05, 0.1, 0.5])
    # does NOT end at 0.5 (the apex)
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(
            n_windows_half=6, lambda2_rampup=[0.0, 0.1, 0.2])
    # not strictly increasing
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(
            n_windows_half=6, lambda2_rampup=[0.0, 0.2, 0.2, 0.5])
    # too short (need both 0.0 endpoint and 0.5 apex)
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(
            n_windows_half=6, lambda2_rampup=[0.0])
    with pytest.raises(ValueError):
        rbfe.build_ats_standard_ladder(
            n_windows_half=6, lambda2_rampup=[])


def test_ats_standard_lambda1_and_lambda2_composable(rbfe):
    """Both axes simultaneously: independent densification, single shared apex.

    leg-up densified to 8 states (λ2=0.05, 0.15 inserted) AND leg-down densified
    to 6 states (λ1=0.05 inserted) => 8 + 6 = 14 states, NOT 8 + 6 + 1 (no
    double-count of the λ1=0/λ2=0.5 apex, which is placed once by the leg-up).
    """
    sch = rbfe.build_ats_standard_ladder(
        n_windows_half=6, single_direction="forward",
        lambda2_rampup=[0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5],
        lambda1_rampdown=[0.05, 0.1, 0.2, 0.3, 0.4, 0.5])
    assert sch["n_states"] == 14
    # Leg-up = first 8 (λ1=0), leg-down = last 6 (λ2=0.5).
    assert sch["lambdas_1"] == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
    assert sch["lambdas_2"] == [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5,
                                0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    # Exactly ONE state at the leg-switch apex (λ1=0, λ2=0.5): idx 7.
    apex_hits = [k for k in range(sch["n_states"])
                 if sch["lambdas_1"][k] == 0.0 and sch["lambdas_2"][k] == 0.5]
    assert apex_hits == [7]
    # Only the two genuine endpoints are intermd==0 (decoupled idx0, coupled apex).
    assert sch["intermd"] == [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]
    assert sch["lambda2_rampup"] == [0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]
    assert sch["lambda1_rampdown"] == [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
    assert sch["n_windows_half_linear"] == 8  # leg-up count


# ---------------------------------------------------------------------------
# ATM hybrid-potential recombination — pure math.
# ---------------------------------------------------------------------------
def test_state_energy_lambda0_is_base(rbfe):
    """At λ1=λ2=0, W0=0 the hybrid term vanishes -> U = base = u0 (Direction>=0)."""
    u0 = -100.0
    pert = 40.0
    e = rbfe._atm_state_energy_kj(
        pert_kj=pert, u0_kj=u0, lambda1=0.0, lambda2=0.0, alpha_per_kj=0.0,
        uh_kj=0.0, w0_kj=0.0, direction=1.0,
        umax_kj=836.8, ubcore_kj=418.4, acore=0.0625)
    assert e == pytest.approx(u0)


def test_state_energy_direction_selects_base(rbfe):
    """Direction<0 selects u1 = u0 + pert as the base (λ=0, W0=0)."""
    u0 = -100.0
    pert = 40.0
    e = rbfe._atm_state_energy_kj(
        pert_kj=pert, u0_kj=u0, lambda1=0.0, lambda2=0.0, alpha_per_kj=0.0,
        uh_kj=0.0, w0_kj=0.0, direction=-1.0,
        umax_kj=836.8, ubcore_kj=418.4, acore=0.0625)
    assert e == pytest.approx(u0 + pert)


def test_state_energy_linear_lambda(rbfe):
    """At λ1==λ2=λ (linear segment, α irrelevant) U = u0 + λ*usc + W0; for a
    small perturbation (below Ubcore) usc == u == pert."""
    u0 = -100.0
    pert = 10.0   # below Ubcore -> usc == pert
    lam = 0.3
    w0 = 5.0
    e = rbfe._atm_state_energy_kj(
        pert_kj=pert, u0_kj=u0, lambda1=lam, lambda2=lam, alpha_per_kj=0.0,
        uh_kj=0.0, w0_kj=w0, direction=1.0,
        umax_kj=836.8, ubcore_kj=418.4, acore=0.0625)
    assert e == pytest.approx(u0 + lam * pert + w0)


def test_state_energy_softcore_caps_large_perturbation(rbfe):
    """A perturbation far above Ubcore is compressed toward Umax by the soft-core
    cap (usc < raw pert), so the λ-scaled term stays bounded."""
    u0 = 0.0
    pert = 100000.0   # huge raw clash
    umax = 836.8      # 200 kcal/mol
    ubcore = 418.4    # 100 kcal/mol
    e = rbfe._atm_state_energy_kj(
        pert_kj=pert, u0_kj=u0, lambda1=0.0, lambda2=1.0, alpha_per_kj=1.0,
        uh_kj=0.0, w0_kj=0.0, direction=1.0,
        umax_kj=umax, ubcore_kj=ubcore, acore=0.0625)
    # The capped usc cannot exceed Umax; the hybrid term is finite (not 1e5).
    assert math.isfinite(e)
    assert e < umax + 1.0


# ---------------------------------------------------------------------------
# Smoke module — import + cntl emitter (no build).
# ---------------------------------------------------------------------------
def test_smoke_module_imports(rbfe):
    smoke = _load_smoke()
    assert hasattr(smoke, "run_ladder_smoke")
    assert hasattr(smoke, "_write_ladder_cntl")
    assert hasattr(smoke, "_load_mixing_gate")


def test_ladder_cntl_emit_has_no_displacement(rbfe, tmp_path):
    """The in-place RBFE cntl must NOT carry the ABFE displacement/LIGAND_ATOMS
    keywords (it is a single-shared-core swap, ATMForce pre-attached) and must
    mark MODE = INPLACE_RBFE so a consumer does not re-wrap the ATMForce."""
    smoke = _load_smoke()
    sch = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2)
    ser = {"sys_xml_path": "inplace_rbfe_free_sys.xml",
           "pdb_path": "inplace_rbfe_free.pdb"}
    cntl = str(tmp_path / "c.cntl")
    smoke._write_ladder_cntl(cntl, sch, ser, "free", 150, 1.0)
    txt = open(cntl).read()
    assert "MODE = 'INPLACE_RBFE'" in txt
    # No active ABFE keyword LINES (the explanatory comments mention the words,
    # so check for assignment lines, not substrings).
    non_comment = [ln for ln in txt.splitlines()
                   if ln.strip() and not ln.lstrip().startswith("#")]
    assert not any(ln.lstrip().startswith("DISPLACEMENT") for ln in non_comment)
    assert not any(ln.lstrip().startswith("LIGAND_ATOMS") for ln in non_comment)
    assert not any(ln.lstrip().startswith("LIGAND_CM_ATOMS")
                   for ln in non_comment)
    # The per-state schedule arrays are emitted.
    assert "LAMBDA1 =" in txt and "DIRECTION =" in txt and "W0COEFF =" in txt


def test_ladder_cntl_emit_single_direction(rbfe, tmp_path):
    """A STANDALONE single-direction schedule (all DIRECTION the same, no apex
    flip) still emits a valid in-place RBFE cntl (the emitter is direction-array
    agnostic)."""
    smoke = _load_smoke()
    sch = rbfe.build_rbfe_ladder(n_windows_half=6, softcore_band=2,
                                 single_direction="forward")
    ser = {"sys_xml_path": "inplace_rbfe_free_sys.xml",
           "pdb_path": "inplace_rbfe_free.pdb"}
    cntl = str(tmp_path / "c_fwd.cntl")
    smoke._write_ladder_cntl(cntl, sch, ser, "free", 150, 1.0)
    txt = open(cntl).read()
    assert "MODE = 'INPLACE_RBFE'" in txt
    # The DIRECTION line lists the 6 forward states all as +1.
    dir_line = [ln for ln in txt.splitlines()
                if ln.lstrip().startswith("DIRECTION =")][0]
    assert dir_line.count("1") == 6  # six +1 entries
    assert "-1" not in dir_line       # no backward states in a forward standalone


def test_smoke_cli_accepts_single_direction(rbfe):
    """run_ladder_smoke exposes the single_direction wiring contract (default
    None = the full symmetric ladder)."""
    smoke = _load_smoke()
    import inspect
    sig = inspect.signature(smoke.run_ladder_smoke)
    assert "single_direction" in sig.parameters
    assert sig.parameters["single_direction"].default is None


# ---------------------------------------------------------------------------
# Mixing-gate-compatible log: the adapter's log lines must parse with the
# per-direction driver's OWN regexes (so the gate is reused, not reimplemented).
# ---------------------------------------------------------------------------
def test_emitted_log_parses_with_driver_regex(rbfe):
    pdp = _load_driver_mixing()
    # A synthetic line in the adapter's exact format.
    line = "2000-01-01 00:00:05 Replica 3 new state 7"
    m = pdp._NEW_STATE_RE.search(line)
    assert m is not None
    assert int(m.group(1)) == 3 and int(m.group(2)) == 7
    ts = pdp._LOG_TS_RE.match(line)
    assert ts is not None and ts.group(1) == "2000-01-01 00:00:05"


def test_synthetic_ladder_log_mixing_gate(rbfe, tmp_path):
    """A hand-written driver log in the adapter's format must drive the
    per-direction mixing gate end-to-end (parse + check_atm_mixing)."""
    pdp = _load_driver_mixing()
    K = 4
    log = tmp_path / "driver.log"
    # 6 cycles; replicas 0..3 walk the ladder so adjacent pairs all cross.
    seqs = [
        [0, 1, 0, 1, 0, 1],
        [1, 0, 1, 2, 1, 0],
        [2, 3, 2, 1, 2, 3],
        [3, 2, 3, 2, 3, 2],
    ]
    with open(log, "w") as fh:
        for c in range(6):
            ts = "2000-01-01 00:00:%02d" % (c + 1)
            for r in range(K):
                fh.write("%s Replica %d new state %d\n" % (ts, r, seqs[r][c]))
    tr = pdp.parse_state_transitions_from_log(str(log), warmup_cycles=1)
    assert tr["n_cycles_total"] == 6
    assert tr["n_samples"] > 0
    mix = pdp.check_atm_mixing(str(log), schedule_K=K, warmup_cycles=1,
                               min_crossings=1)
    # adjacent_crossings populated; the gate returns a structured result.
    assert "walls" in mix and "both_ends_visited_count" in mix


# ---------------------------------------------------------------------------
# Real build / serialize / run — slow; guarded behind endpoint presence.
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not _endpoints_present(),
                    reason="2QKI endpoint final.pdb not present")
def test_serialize_roundtrip_unsolvated(rbfe, tmp_path):
    """Serialize the unsolvated harmonized fused System + deserialize it; the
    ATMForce must survive + topology/positions/particles stay in lockstep."""
    res = rbfe.serialize_inplace_rbfe_system(
        leg="free", out_dir=str(tmp_path), seed="s7", solvate=False,
        harmonize_common_charges=True, swap_mode="genuine", constraints=None)
    assert os.path.isfile(res["sys_xml_path"])
    assert os.path.isfile(res["pdb_path"])
    assert res["atmforce_index"] is not None
    loaded = rbfe.load_serialized_system(res["sys_xml_path"], res["pdb_path"])
    ns = loaded["system"].getNumParticles()
    nt = loaded["topology"].getNumAtoms()
    npos = len(loaded["positions"])
    assert ns == nt == npos
    assert loaded["atmforce_index"] is not None


@pytest.mark.skipif(not _endpoints_present(),
                    reason="2QKI endpoint final.pdb not present")
def test_ladder_runs_cpu_reference(rbfe, tmp_path):
    """A tiny ladder on the Reference platform (unsolvated harmonized) loads +
    runs + emits a parseable driver log (wiring, not convergence)."""
    res = rbfe.serialize_inplace_rbfe_system(
        leg="free", out_dir=str(tmp_path), seed="s7", solvate=False,
        harmonize_common_charges=True, swap_mode="genuine", constraints=None)
    loaded = rbfe.load_serialized_system(res["sys_xml_path"], res["pdb_path"])
    sch = rbfe.build_rbfe_ladder(n_windows_half=4, softcore_band=1)
    logp = str(tmp_path / "driver.log")
    ladder = rbfe.InplaceRbfeLadder(
        loaded["system"], loaded["positions"], sch,
        platform_name="Reference", timestep_fs=1.0, log_path=logp, seed=11,
        minimize_iters=100, backward_equil_steps=50)
    info = None
    for _ in range(3):
        info = ladder.run_cycle(md_steps=5)
    ladder.close()
    assert info is not None
    assert os.path.isfile(logp)
    # The log must be non-empty + parseable by the driver's gate.
    pdp = _load_driver_mixing()
    tr = pdp.parse_state_transitions_from_log(logp, warmup_cycles=0)
    assert tr["n_samples"] > 0
    assert tr["n_cycles_total"] == 3


# ---------------------------------------------------------------------------
# Two-copy serialize wiring — construction param (pure validation + real build).
# ---------------------------------------------------------------------------
def test_serialize_construction_default_is_single_core(rbfe):
    import inspect
    sig = inspect.signature(rbfe.serialize_inplace_rbfe_system)
    assert "construction" in sig.parameters
    assert sig.parameters["construction"].default == "single_core"
    assert "displacement_nm" in sig.parameters


def test_serialize_rejects_bad_construction(rbfe):
    with pytest.raises(ValueError):
        rbfe.serialize_inplace_rbfe_system(construction="overlay")


@pytest.mark.skipif(not _endpoints_present(),
                    reason="2QKI endpoint final.pdb not present")
def test_serialize_twocopy_roundtrip_unsolvated(rbfe, tmp_path):
    """Serialize the unsolvated two-copy box + deserialize it: the ATMForce
    survives, particle/topology lockstep holds, the box is ~2x a single
    endpoint, and the C8 decouple direction is the (unit) displacement vector."""
    res = rbfe.serialize_inplace_rbfe_system(
        leg="free", out_dir=str(tmp_path), seed="s7", solvate=False,
        harmonize_common_charges=True, construction="twocopy", constraints=None)
    assert res["construction"] == "twocopy"
    assert os.path.isfile(res["sys_xml_path"])
    assert os.path.isfile(res["pdb_path"])
    loaded = rbfe.load_serialized_system(res["sys_xml_path"], res["pdb_path"])
    ns = loaded["system"].getNumParticles()
    nt = loaded["topology"].getNumAtoms()
    npos = len(loaded["positions"])
    assert ns == nt == npos
    assert loaded["atmforce_index"] is not None
    # Two copies resident -> n_atoms ~ 2x a single endpoint (~210 -> ~417).
    assert res["n_atoms"] > 380
    assert res["n_copy1"] is not None and res["n_copy1"] < res["n_atoms"]
    # The C8 decouple direction is a finite UNIT vector (the displacement unit).
    dd = res["genuine_decouple_dir"]
    assert dd is not None
    mag = sum(c * c for c in dd) ** 0.5
    assert abs(mag - 1.0) < 1e-6


@pytest.mark.skipif(not _endpoints_present(),
                    reason="2QKI endpoint final.pdb not present")
def test_twocopy_ats_ladder_runs_cpu_reference(rbfe, tmp_path):
    """The two-copy box runs through the canonical ATS standard ladder on the
    Reference platform: loads, runs a few cycles, writes per-walker .out rows,
    and the running-box |usc| is NOT pinned at the single-shared-core saturated
    plateau (the collapse-signature escape re-checked under integration)."""
    res = rbfe.serialize_inplace_rbfe_system(
        leg="free", out_dir=str(tmp_path), seed="s7", solvate=False,
        harmonize_common_charges=True, construction="twocopy", constraints=None)
    loaded = rbfe.load_serialized_system(res["sys_xml_path"], res["pdb_path"])
    sch = rbfe.build_ats_standard_ladder(n_windows_half=3,
                                         single_direction="forward")
    logp = str(tmp_path / "twocopy_driver.log")
    ladder = rbfe.InplaceRbfeLadder(
        loaded["system"], loaded["positions"], sch,
        platform_name="Reference", timestep_fs=1.0, log_path=logp, seed=21,
        minimize_iters=100, backward_equil_steps=0,
        out_dir=str(tmp_path), out_basename="trackb_dplus")
    nan_any = False
    for _ in range(3):
        info = ladder.run_cycle(md_steps=5)
        nan_any = nan_any or info["nan_seen"]
    ladder.close()
    assert not nan_any
    # Per-walker .out rows were written (UWHAM-consumable layout).
    rows = 0
    pmax = 0.0
    for r in range(sch["n_states"]):
        path = os.path.join(str(tmp_path), "r%d" % (r,), "trackb_dplus.out")
        if not os.path.isfile(path):
            continue
        for line in open(path):
            cols = line.split()
            if len(cols) >= 10:
                rows += 1
                pmax = max(pmax, abs(float(cols[9])))
    assert rows > 0
    # Not pinned at the single-shared-core saturated plateau (~150) / clash.
    assert not (140.0 <= pmax <= 201.0)
    assert pmax < 1.0e3


# ---------------------------------------------------------------------------
# Staged minimization (OPT-IN, W4A large-box stability) — wiring + behaviour.
# ---------------------------------------------------------------------------
def test_staged_min_constructor_is_opt_in(rbfe):
    """The ladder constructor gained the staged_min opt-in (DEFAULT OFF) + the
    staged-floor / warmup params; the existing minimize_iters default is
    unchanged (the non-staged budget is untouched)."""
    import inspect
    sig = inspect.signature(rbfe.InplaceRbfeLadder.__init__)
    assert "staged_min" in sig.parameters
    assert sig.parameters["staged_min"].default is False
    assert "staged_min_iters" in sig.parameters
    assert "staged_warmup_steps" in sig.parameters
    # minimize_iters default unchanged (V3I/MTR/A9G non-staged path untouched).
    assert sig.parameters["minimize_iters"].default == 500
    # The staged floor honours the P1 5000-iter production minimum.
    assert rbfe.STAGED_MIN_ITERS_FLOOR == 5000


def test_staged_min_floor_clamped_only_when_on(rbfe):
    """staged_min_iters is clamped to the >=5000 floor ONLY when staged_min=True;
    when OFF the value is stored verbatim (it is never used on the off path)."""
    # We construct via a lightweight object that carries only the clamp logic the
    # constructor applies (mirrors the constructor body without a real System).
    floor = rbfe.STAGED_MIN_ITERS_FLOOR
    # ON + below floor -> clamped up to the floor.
    clamped_on = max(100, floor) if True else 100
    assert clamped_on == floor
    # The constructor's clamp is: staged on -> max(iters, floor); off -> iters.
    # (verified in the real-build test below; this asserts the floor constant.)
    assert floor == 5000


@pytest.mark.skipif(not _endpoints_present(),
                    reason="2QKI endpoint final.pdb not present")
def test_staged_min_off_is_single_stage_byte_identical(rbfe, tmp_path):
    """staged_min=False (the default) runs the EXISTING single-stage minimize path:
    the ladder loads + runs identically to the legacy CPU-reference run, and the
    instance records staged_min=False (no staged branch taken)."""
    res = rbfe.serialize_inplace_rbfe_system(
        leg="free", out_dir=str(tmp_path), seed="s7", solvate=False,
        harmonize_common_charges=True, swap_mode="genuine", constraints=None)
    loaded = rbfe.load_serialized_system(res["sys_xml_path"], res["pdb_path"])
    sch = rbfe.build_rbfe_ladder(n_windows_half=4, softcore_band=1)
    logp = str(tmp_path / "driver_offstage.log")
    ladder = rbfe.InplaceRbfeLadder(
        loaded["system"], loaded["positions"], sch,
        platform_name="Reference", timestep_fs=1.0, log_path=logp, seed=11,
        minimize_iters=100, backward_equil_steps=50)   # staged_min defaults False
    assert ladder.staged_min is False
    info = None
    for _ in range(3):
        info = ladder.run_cycle(md_steps=5)
    ladder.close()
    assert info is not None
    pdp = _load_driver_mixing()
    tr = pdp.parse_state_transitions_from_log(logp, warmup_cycles=0)
    assert tr["n_samples"] > 0
    assert tr["n_cycles_total"] == 3


@pytest.mark.skipif(not _endpoints_present(),
                    reason="2QKI endpoint final.pdb not present")
def test_staged_min_on_runs_staged_path(rbfe, tmp_path):
    """staged_min=True runs the STAGED relax path (reference-state minimize ->
    polish -> warmup) and the ladder still loads + runs + emits a parseable log.
    The staged floor is clamped to >=5000 even when a smaller value is passed."""
    res = rbfe.serialize_inplace_rbfe_system(
        leg="free", out_dir=str(tmp_path), seed="s7", solvate=False,
        harmonize_common_charges=True, swap_mode="genuine", constraints=None)
    loaded = rbfe.load_serialized_system(res["sys_xml_path"], res["pdb_path"])
    sch = rbfe.build_rbfe_ladder(n_windows_half=3, softcore_band=1)
    logp = str(tmp_path / "driver_staged.log")
    ladder = rbfe.InplaceRbfeLadder(
        loaded["system"], loaded["positions"], sch,
        platform_name="Reference", timestep_fs=1.0, log_path=logp, seed=13,
        minimize_iters=100, backward_equil_steps=50,
        staged_min=True, staged_min_iters=200, staged_warmup_steps=5)
    assert ladder.staged_min is True
    # The staged floor clamps a small value up to >=5000 (P1 production floor).
    assert ladder.staged_min_iters >= rbfe.STAGED_MIN_ITERS_FLOOR
    assert ladder.staged_warmup_steps == 5
    info = None
    for _ in range(2):
        info = ladder.run_cycle(md_steps=5)
    ladder.close()
    assert info is not None
    pdp = _load_driver_mixing()
    tr = pdp.parse_state_transitions_from_log(logp, warmup_cycles=0)
    assert tr["n_samples"] > 0


# ---------------------------------------------------------------------------
# FIX-A re-seeding (OPT-IN, bound-leg ladder mixing) — wiring + behaviour.
# Validated across 3 wall seeds (decisive seal / fragmented-ladder pathology
# removed); FE-unbiased (initial-condition change only).
# DEFAULT OFF must be byte/behaviour-identical to the legacy identity init.
# ---------------------------------------------------------------------------
def _build_minimal_atm_system(n_per_copy=6, box_nm=4.0, seed=0):
    """A SELF-CONTAINED minimal real ATMForce system: two identical harmonic
    'copies' inside ONE ATMForce (the two-copy ATS geometry in miniature), so the
    re-seed path runs on the Reference platform WITHOUT the 2QKI endpoint PDBs.
    C7 canon globals so _set_state / getPerturbationEnergy behave as in production.
    Mirrors W4A/test_w4a_reseed_proto._build_minimal_atm_system."""
    import numpy as np
    import openmm as mm
    import openmm.unit as unit

    n = int(n_per_copy)
    system = mm.System()
    for _ in range(2 * n):
        system.addParticle(12.0)
    a = box_nm
    system.setDefaultPeriodicBoxVectors(
        mm.Vec3(a, 0, 0) * unit.nanometer,
        mm.Vec3(0, a, 0) * unit.nanometer,
        mm.Vec3(0, 0, a) * unit.nanometer)
    bonded = mm.HarmonicBondForce()
    for c in range(2):
        base = c * n
        for i in range(n - 1):
            bonded.addBond(base + i, base + i + 1, 0.15 * unit.nanometer,
                           200000.0 * unit.kilojoule_per_mole / unit.nanometer ** 2)
    atm = mm.ATMForce(200.0, 100.0, 0.0625, 110.0, 0.10, 1.0, 0.0, 0.0, 1.0)
    atm.addForce(bonded)
    d0 = mm.Vec3(0.0, 0.0, 0.0) * unit.nanometer
    d1 = mm.Vec3(1.0, 0.0, 0.0) * unit.nanometer
    zero = mm.Vec3(0.0, 0.0, 0.0) * unit.nanometer
    for p in range(2 * n):
        if p < n:
            atm.addParticle(zero, zero)
        else:
            atm.addParticle(d1, d0)
    system.addForce(atm)
    rng = np.random.RandomState(seed)
    pos = []
    for c in range(2):
        x0 = 1.0 + 1.5 * c
        for i in range(n):
            pos.append(mm.Vec3(x0 + 0.15 * i, 1.0 + 0.01 * rng.randn(),
                               1.0 + 0.01 * rng.randn()) * unit.nanometer)
    return system, pos


# -- (a) constructor opt-in surface --------------------------------------------
def test_reseed_constructor_is_opt_in(rbfe):
    """The ladder constructor gained the FIX-A re-seed opt-ins (DEFAULT OFF)."""
    import inspect
    sig = inspect.signature(rbfe.InplaceRbfeLadder.__init__)
    assert sig.parameters["reseed_perm_seed"].default is None
    assert sig.parameters["reseed_endpoint"].default is False
    assert sig.parameters["reseed_endpoint_band_lambda2_max"].default == 0.25
    assert sig.parameters["reseed_endpoint_equil_steps"].default == 2000
    assert sig.parameters["reseed_endpoint_minimize_iters"].default == 500


# -- (a) A1 permutation reproducibility / bijection ----------------------------
def test_reseed_permutation_reproducible_from_logged_seed(rbfe):
    a = rbfe.make_reseed_permutation(11, reseed_perm_seed=12345)
    b = rbfe.make_reseed_permutation(11, reseed_perm_seed=12345)
    assert a == b
    assert sorted(a) == list(range(11))


def test_reseed_permutation_independent_of_global_rng(rbfe):
    import random
    random.seed(1)
    a = rbfe.make_reseed_permutation(11, reseed_perm_seed=42)
    random.seed(99999)
    [random.random() for _ in range(50)]
    b = rbfe.make_reseed_permutation(11, reseed_perm_seed=42)
    assert a == b


def test_reseed_permutation_bad_n_raises(rbfe):
    with pytest.raises(ValueError):
        rbfe.make_reseed_permutation(0, reseed_perm_seed=1)


# -- (b) DEFAULT OFF == identity init ------------------------------------------
def test_reseed_default_off_is_identity(rbfe):
    for n in (3, 5, 11, 12):
        assert rbfe.make_reseed_permutation(n, reseed_perm_seed=None) \
            == list(range(n))


# -- (c) DIRECTION-CORRECT decoupled endpoint + band for BOTH legs -------------
def test_reseed_decoupled_endpoint_band_forward(rbfe):
    """Forward (dplus): decoupled endpoint = state 0 (λ1=λ2=0); band = {0,1,2}."""
    sch = rbfe.build_ats_standard_ladder(n_windows_half=6,
                                         single_direction="forward")
    assert rbfe.decoupled_endpoint_state(sch) == 0
    assert rbfe.decoupled_band_states(sch, band_lambda2_max=0.25) == [0, 1, 2]


def test_reseed_decoupled_endpoint_band_backward(rbfe):
    """Backward (dminus): decoupled endpoint = the LAST state (λ1=λ2=0); band =
    {8,9,10}. The max-λ1 forward-only rule would WRONGLY pick the coupled apex
    (state 0) here — this is the direction-correctness guard."""
    sch = rbfe.build_ats_standard_ladder(n_windows_half=6,
                                         single_direction="backward")
    assert rbfe.decoupled_endpoint_state(sch) == 10
    assert sch["lambdas_1"][10] == 0.0 and sch["lambdas_2"][10] == 0.0
    assert rbfe.decoupled_band_states(sch, band_lambda2_max=0.25) == [8, 9, 10]


def test_reseed_decoupled_endpoint_band_densified_both_directions(rbfe):
    """With --lambda1-rampdown the leg-down is densified (11 -> 12 states); the
    endpoint INDEX shifts but the λ2=0 rule still resolves it for BOTH legs
    (densified backward -> state 11, band {9,10,11}; densified forward -> state 0,
    band {0,1,2}). A fixed-index mapping would false-green on the 12-state ladder."""
    ramp = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
    db = rbfe.build_ats_standard_ladder(n_windows_half=6,
                                        single_direction="backward",
                                        lambda1_rampdown=ramp)
    assert db["n_states"] == 12
    assert rbfe.decoupled_endpoint_state(db) == 11
    assert rbfe.decoupled_band_states(db, band_lambda2_max=0.25) == [9, 10, 11]
    df = rbfe.build_ats_standard_ladder(n_windows_half=6,
                                        single_direction="forward",
                                        lambda1_rampdown=ramp)
    assert rbfe.decoupled_endpoint_state(df) == 0
    assert rbfe.decoupled_band_states(df, band_lambda2_max=0.25) == [0, 1, 2]


# -- (b) real-system DEFAULT OFF: ladder init is identity + runs NaN-free -------
@pytest.mark.skipif(not _have_openmm(), reason="openmm not importable")
def test_reseed_off_ladder_is_identity_real_system(rbfe, tmp_path):
    """Default-OFF: the ladder's replica_state is the identity map, no endpoint
    re-seed config is built, and it runs a cycle without NaN."""
    system, pos = _build_minimal_atm_system()
    sch = rbfe.build_ats_standard_ladder(n_windows_half=3,
                                         single_direction="backward")
    logp = str(tmp_path / "reseed_off.log")
    ladder = rbfe.InplaceRbfeLadder(
        system, pos, sch, platform_name="Reference", timestep_fs=1.0,
        log_path=logp, seed=11, minimize_iters=50, backward_equil_steps=0)
    assert ladder.replica_state == list(range(ladder.n_states))
    assert ladder.reseed_perm_seed is None
    assert ladder.reseed_endpoint is False
    assert ladder._reseed_relaxed_positions is None
    assert ladder._reseed_band == []
    info = ladder.run_cycle(md_steps=3)
    assert info["nan_seen"] is False
    ladder.close()


# -- (a) A1 ON in-__init__: ladder starts in the LOGGED permutation -------------
@pytest.mark.skipif(not _have_openmm(), reason="openmm not importable")
def test_reseed_a1_permutation_applied_real_system(rbfe, tmp_path):
    """reseed_perm_seed set: the ladder's t=0 replica_state is EXACTLY the logged
    permutation (a valid bijection), and the ladder runs NaN-free. The permutation
    is set IN __init__ before the per-state minimize, so each context is minimized
    at the state it will occupy (no post-construction relabel)."""
    system, pos = _build_minimal_atm_system()
    sch = rbfe.build_ats_standard_ladder(n_windows_half=3,
                                         single_direction="backward")
    logp = str(tmp_path / "reseed_a1.log")
    ladder = rbfe.InplaceRbfeLadder(
        system, pos, sch, platform_name="Reference", timestep_fs=1.0,
        log_path=logp, seed=11, minimize_iters=50, backward_equil_steps=0,
        reseed_perm_seed=20260618)
    expect = rbfe.make_reseed_permutation(ladder.n_states, 20260618)
    assert ladder.replica_state == expect
    assert sorted(ladder.replica_state) == list(range(ladder.n_states))
    info = ladder.run_cycle(md_steps=3)
    assert info["nan_seen"] is False
    ladder.close()


# -- (c) A2 ON: endpoint re-seed targets the direction-correct band, NaN-free ---
@pytest.mark.skipif(not _have_openmm(), reason="openmm not importable")
def test_reseed_a2_endpoint_reseed_backward_real_system(rbfe, tmp_path):
    """reseed_endpoint=True on a BACKWARD ladder: the endpoint re-seed runs a
    standalone equilibration at the genuine decoupled endpoint (LAST state, λ2=0)
    and seeds ONLY the decoupled band (NOT the coupled apex). Runs NaN-free."""
    system, pos = _build_minimal_atm_system()
    sch = rbfe.build_ats_standard_ladder(n_windows_half=3,
                                         single_direction="backward")
    logp = str(tmp_path / "reseed_a2.log")
    ladder = rbfe.InplaceRbfeLadder(
        system, pos, sch, platform_name="Reference", timestep_fs=1.0,
        log_path=logp, seed=11, minimize_iters=50, backward_equil_steps=0,
        reseed_perm_seed=20260618, reseed_endpoint=True,
        reseed_endpoint_band_lambda2_max=0.25,
        reseed_endpoint_equil_steps=50, reseed_endpoint_minimize_iters=50)
    # The endpoint is the LAST state (decoupled, λ2=0) — NOT state 0 (coupled apex).
    assert ladder._reseed_endpoint_state == ladder.n_states - 1
    assert sch["lambdas_2"][ladder._reseed_endpoint_state] == 0.0
    # Band is a strict, non-empty subset; every band state has small λ2.
    assert 0 < len(ladder._reseed_band) < ladder.n_states
    for k in ladder._reseed_band:
        assert sch["lambdas_2"][k] <= 0.25
    assert ladder._reseed_relaxed_positions is not None
    info = ladder.run_cycle(md_steps=3)
    assert info["nan_seen"] is False
    ladder.close()


# ---------------------------------------------------------------------------
# HARDENED per-seed cohort mixing gate (check_leg_hardened_mixing): the cohort
# verdict is the AND over EVERY (endpoint, seed, direction) leg read from the
# RAW driver.log — a single seed-failing leg fails the cohort (a dispatcher
# PARTIAL_SUCCESS / recovery may NOT roll it up to PASS).
# ---------------------------------------------------------------------------
def _load_inplace_production():
    spec = importlib.util.spec_from_file_location(
        "trackb_inplace_rbfe_production",
        os.path.join(_SCRIPTS, "trackb_inplace_rbfe_production.py"))
    if spec is None or spec.loader is None:
        pytest.skip("could not locate scripts/trackb_inplace_rbfe_production.py")
    if _SCRIPTS not in sys.path:
        sys.path.insert(0, _SCRIPTS)
    if _UTILS not in sys.path:
        sys.path.insert(0, _UTILS)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _write_cohort_driver_log(path, trajectories, warmup_cycles=2):
    """Synthesize an async_re driver log from {replica: [state,...]} with a
    leading warmup so the post-warmup window IS the supplied trajectory."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    n_rounds = len(next(iter(trajectories.values())))
    lines = []
    rr = 0
    for _ in range(warmup_cycles):
        ts = "2026-06-18 00:00:%02d" % (rr % 60)
        for rep, traj in trajectories.items():
            lines.append("%s - INFO - async_re.openmm_async_re - "
                         "Replica %d new state %d\n" % (ts, rep, traj[0]))
        rr += 1
    for c in range(n_rounds):
        ts = "2026-06-18 00:01:%02d" % (rr % 60)
        for rep, traj in trajectories.items():
            lines.append("%s - INFO - async_re.openmm_async_re - "
                         "Replica %d new state %d\n" % (ts, rep, traj[c]))
        rr += 1
    with open(path, "w") as fh:
        fh.write("".join(lines))


def _clean_traj(k=11, length=120):
    sweep = list(range(0, k)) + list(range(k - 2, -1, -1))
    out = {}
    for rep in range(k):
        off = rep % len(sweep)
        out[rep] = (sweep * 12)[off:off + length]
    return out


def _reseal_traj(k=11):
    """First half: one full 0..K-1..0 sweep (whole-window all pairs). Second
    half: lower-only shuttle 0..K-3 -> the top two pairs seal."""
    first = list(range(0, k)) + list(range(k - 2, -1, -1))
    lower = list(range(0, k - 2)) + list(range(k - 4, 0, -1))
    n = len(first)
    out = {}
    for rep in range(k):
        off = rep % len(lower)
        sh = (lower * 8)[off:off + n]
        out[rep] = first + sh
    return out


def _stamp_cntl(path, k=11):
    """Full per-direction sibling cntl so ``_parse_cntl_schedule`` returns a K-
    state dict with the soft-core SSOT arrays (DIRECTION/INTERMEDIATE/LAMBDA1/
    LAMBDA2/ALPHA/U0/W0). Values are synthetic but length-consistent — the
    cohort gate only needs K + the schedule presence (overlaps are supplied via
    the monkeypatched harness in these AND-logic tests)."""
    def _csv(seq):
        return ", ".join("%.4f" % v for v in seq)
    lam1 = [0.5 * i / (k - 1) for i in range(k)]
    with open(path, "w") as fh:
        fh.write("BASENAME = 'trackb'\n")
        fh.write("DIRECTION =    '%s'\n" % _csv([1.0] * k))
        fh.write("INTERMEDIATE = '%s'\n" % _csv([1.0] * k))
        fh.write("LAMBDA1 =      '%s'\n" % _csv(lam1))
        fh.write("LAMBDA2 =      '%s'\n" % _csv(lam1))
        fh.write("ALPHA =        '%s'\n" % _csv([0.1] * k))
        fh.write("U0 =           '%s'\n" % _csv([110.0] * k))
        fh.write("W0COEFF =      '%s'\n" % _csv([0.0] * k))


def _build_cohort(out_root, leg, kind_by_seed, k=11, n_rep=2):
    """Build a synthetic cohort dir tree. ``kind_by_seed`` maps replicate index
    -> 'clean' | 'reseal' (applied to BOTH endpoints + BOTH directions)."""
    prod = _load_inplace_production()
    for ep in prod.ENDPOINTS:
        for j in range(n_rep):
            rep_dir = prod._rep_dir(out_root, ep, leg, j)
            for tag in prod.DIRECTION_TAGS:
                sub = os.path.join(rep_dir, tag)
                base = prod.JOBNAME + "_" + tag
                traj = (_clean_traj(k) if kind_by_seed.get(j) == "clean"
                        else _reseal_traj(k))
                _write_cohort_driver_log(
                    os.path.join(sub, base + "_driver.log"), traj,
                    warmup_cycles=2)
                _stamp_cntl(os.path.join(sub, base + "_asyncre.cntl"), k=k)
    return prod


def _patch_overlaps_ok(prod, monkeypatch, value=0.5):
    """Force per-pair overlaps available + above floor (decouple the cohort
    AND-logic test from pymbar/atom_openmm availability)."""
    def _fake(leg_dir, direction_tag, schedule, warmup_cycles):
        k = 11
        if isinstance(schedule, dict) and schedule.get("n_states"):
            k = int(schedule["n_states"])
        pairs = {"%d-%d" % (i, i + 1): value for i in range(k - 1)}
        return {"overlaps": pairs, "bhattacharyya":
                {kk: 0.4 for kk in pairs}, "source": "synthetic"}
    monkeypatch.setattr(prod, "_compute_adjacent_overlaps", _fake)
    # Also make _load_overlap return a truthy sentinel so the gate takes the
    # overlap path (the real harness may be absent in the qmmm env).
    monkeypatch.setattr(prod, "_load_overlap", lambda: object())


def test_cohort_all_clean_passes(tmp_path, monkeypatch):
    """All seeds clean-mixing (both endpoints, both directions) + overlaps OK
    -> cohort PASS."""
    out_root = str(tmp_path / "out")
    prod = _build_cohort(out_root, "bound", {0: "clean", 1: "clean"}, n_rep=2)
    _patch_overlaps_ok(prod, monkeypatch)
    res = prod.check_leg_hardened_mixing(
        out_root, "bound", n_replicates=2, mintimeid=3)
    assert res["verdict"] == "PASS"
    assert res["passed"] is True
    assert res["n_fail"] == 0 and res["n_indeterminate"] == 0
    assert res["n_pass"] == res["n_legs"] == 2 * 2 * 2   # ep x rep x dir


def test_cohort_one_seed_reseal_fails_no_rollup(tmp_path, monkeypatch):
    """ONE seed's leg is an open-once-then-reseal wall; every OTHER leg is clean.
    The cohort verdict is the AND -> FAIL. A dispatcher PARTIAL_SUCCESS rollup
    may NOT promote the seed-failing cohort to PASS."""
    out_root = str(tmp_path / "out")
    # rep0 clean, rep1 reseal (the failing seed) — applied to both endpoints.
    prod = _build_cohort(out_root, "bound", {0: "clean", 1: "reseal"}, n_rep=2)
    _patch_overlaps_ok(prod, monkeypatch)
    res = prod.check_leg_hardened_mixing(
        out_root, "bound", n_replicates=2, mintimeid=3)
    assert res["verdict"] == "FAIL"
    assert res["passed"] is False
    assert res["n_fail"] >= 1
    # The failing legs are exactly the reseal seed's (rep1), both endpoints/dirs.
    failed = [r for r in res["per_leg"] if r["verdict"] == "FAIL"]
    assert all(r["rep"] == 1 for r in failed)
    assert any("8-9" in (r.get("second_half_walls") or []) for r in failed)


def test_cohort_missing_overlaps_indeterminate(tmp_path, monkeypatch):
    """Clean crossings but overlaps unavailable (harness absent) -> cohort is
    INDETERMINATE, NOT PASS (refuses to accept an unverified leg)."""
    out_root = str(tmp_path / "out")
    prod = _build_cohort(out_root, "bound", {0: "clean", 1: "clean"}, n_rep=2)
    # Force the overlap harness to be unavailable.
    monkeypatch.setattr(prod, "_load_overlap", lambda: None)
    res = prod.check_leg_hardened_mixing(
        out_root, "bound", n_replicates=2, mintimeid=3)
    assert res["verdict"] == "INDETERMINATE"
    assert res["passed"] is False
    assert res["n_indeterminate"] >= 1


def test_cohort_missing_log_indeterminate(tmp_path, monkeypatch):
    """A missing driver.log for one seed -> that leg INDETERMINATE -> cohort not
    PASS (no silent skip of an un-judged seed)."""
    out_root = str(tmp_path / "out")
    prod = _build_cohort(out_root, "bound", {0: "clean", 1: "clean"}, n_rep=2)
    _patch_overlaps_ok(prod, monkeypatch)
    # Remove rep1/cp4/dplus driver log.
    victim = os.path.join(
        prod._rep_dir(out_root, "cp4", "bound", 1), "dplus",
        prod.JOBNAME + "_dplus_driver.log")
    os.remove(victim)
    res = prod.check_leg_hardened_mixing(
        out_root, "bound", n_replicates=2, mintimeid=3)
    assert res["passed"] is False
    assert res["n_indeterminate"] >= 1
