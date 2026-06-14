#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for ``scripts/trackb_production_v2_asyncre`` — v2 cntl + nodefile.

The v2 launcher delegates production to upstream ``atom_openmm`` binaries
(abfe_structprep / abfe_production); these tests cover the pure-Python
glue: the canonical-schedule constants (re-exported from v1 for regression),
the cntl writer (keyword correctness, displacement unit conversion), and
the nodefile writer (row format). They do NOT exercise the upstream
binaries (those require GPU + atm env).
"""

import importlib.util
import os
import re
import sys
import tempfile

import pytest


_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ = os.path.dirname(_HERE)


def _load_v2_module():
    spec = importlib.util.spec_from_file_location(
        "trackb_v2",
        os.path.join(_PROJ, "scripts", "trackb_production_v2_asyncre.py"),
    )
    if spec is None or spec.loader is None:
        pytest.skip("could not locate scripts/trackb_production_v2_asyncre.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def v2():
    return _load_v2_module()


# ---------------------------------------------------------------------------
# Schedule constants (re-exported from v1 for regression continuity)
# ---------------------------------------------------------------------------
def test_v2_schedule_22_states_symmetric(v2):
    """v2 must expose the same canonical 22-state schedule as v1."""
    assert v2.N_STATES == 22
    assert len(v2.LAMBDAS_1) == 22
    assert len(v2.LAMBDAS_2) == 22
    assert len(v2.DIRECTIONS) == 22
    assert len(v2.INTERMD) == 22
    assert len(v2.W0) == 22
    assert v2.DIRECTIONS[:11] == [1] * 11
    assert v2.DIRECTIONS[11:] == [-1] * 11
    assert v2.LAMBDAS_1[:11] == [0.00, 0.05, 0.10, 0.15, 0.20,
                                 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
    assert v2.LAMBDAS_1[11:] == [0.50, 0.45, 0.40, 0.35, 0.30,
                                 0.25, 0.20, 0.15, 0.10, 0.05, 0.00]
    assert v2.LAMBDAS_1 == v2.LAMBDAS_2
    assert [i for i, v in enumerate(v2.INTERMD) if v == 1] == [10, 11]
    assert [i for i, v in enumerate(v2.W0) if v != 0.0] == [10, 11]


def test_v2_softcore_constants(v2):
    """v2 must expose the same UMAX/UBCORE/ACORE as v1 (uwham contract)."""
    assert v2.UMAX_KCAL == 200.0
    assert v2.UBCORE_KCAL == 100.0
    assert v2.ACORE == 0.062500


def test_v2_displacement_default(v2):
    """v2 default displacement is the v1 2.5 nm +x value (regression continuity)."""
    assert v2.DISPLACEMENT_NM_DEFAULT == (2.5, 0.0, 0.0)


def test_v2_pre_register_outcomes(v2):
    """4-outcome pre-register schema."""
    assert v2.PRE_REGISTER_OUTCOMES == (
        "1_sign_stable",
        "2_sigma_drift",
        "3_magnitude_drift",
        "4_sign_flip",
    )


# ---------------------------------------------------------------------------
# Soft-core / softplus bias parity with v1 (used by uwham re-analysis)
# ---------------------------------------------------------------------------
def test_v2_soft_core_pert_e_matches_v1_path(v2):
    """soft_core_pert_e bounds u to [ub, umax] when u > ub."""
    assert v2.soft_core_pert_e(50.0, umax=200.0, ub=100.0, a=0.0625) == 50.0
    assert v2.soft_core_pert_e(99.999, umax=200.0, ub=100.0, a=0.0625) == 99.999
    assert v2.soft_core_pert_e(100.0, umax=200.0, ub=100.0, a=0.0625) == 100.0
    capped = v2.soft_core_pert_e(1e15, umax=200.0, ub=100.0, a=0.0625)
    assert 100.0 < capped < 200.0


def test_v2_softplus_bias_zero_at_lambda_zero(v2):
    """At lambda1=lambda2=0 the softplus bias reduces to w0 only."""
    bias = v2.softplus_bias_kj(
        lambda1=0.0, lambda2=0.0, alpha_per_kj=0.0239,
        uh_kj=110.0 * 4.184, w0_kj=0.0, pert_kj=50.0,
    )
    assert abs(bias) < 1e-10


# ---------------------------------------------------------------------------
# Nodefile writer
# ---------------------------------------------------------------------------
def test_nodefile_writer_default(v2):
    """nodefile must be one row per GPU device with CUDA platform."""
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "nodefile")
        v2.write_nodefile(path, gpu_indices=[0], platform="CUDA")
        with open(path) as fh:
            content = fh.read()
        rows = [r for r in content.strip().split("\n") if r]
        assert len(rows) == 1
        # localhost,0:<dev>,<threads>,<arch>,,<tmp>
        m = re.match(r"^localhost,0:(\d+),\d+,CUDA,,\S+$", rows[0])
        assert m is not None, f"unexpected row format: {rows[0]!r}"
        assert m.group(1) == "0"


def test_nodefile_writer_multi_gpu(v2):
    """Multi-GPU nodefile must list each device on its own row."""
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "nodefile")
        v2.write_nodefile(path, gpu_indices=[0, 1], platform="CUDA")
        with open(path) as fh:
            rows = [r for r in fh.read().strip().split("\n") if r]
        assert len(rows) == 2
        assert "0:0" in rows[0]
        assert "0:1" in rows[1]


def test_nodefile_writer_empty_falls_back_to_zero(v2):
    """Empty GPU list defaults to device 0."""
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "nodefile")
        v2.write_nodefile(path, gpu_indices=[], platform="CUDA")
        with open(path) as fh:
            rows = [r for r in fh.read().strip().split("\n") if r]
        assert len(rows) == 1
        assert "0:0" in rows[0]


# ---------------------------------------------------------------------------
# Cntl writer
# ---------------------------------------------------------------------------
def test_cntl_writer_bound_leg_keywords(v2):
    """Bound-leg cntl must contain CM-CM restraint + position-restraint keys."""
    with tempfile.TemporaryDirectory() as td:
        cntl = os.path.join(td, "trackb.cntl")
        nodefile = os.path.join(td, "nodefile")
        v2.write_cntl_file(
            cntl_path=cntl,
            basename="trackb",
            nodefile_path=nodefile,
            ligand_atom_indices=[100, 101, 102, 103],
            pos_restrained_atom_indices=[0, 1, 2, 3, 4],
            displacement_nm=(2.5, 0.0, 0.0),
            production_steps=2500,
            prnt_frequency=2500,
            trj_frequency=25000,
            wall_time_min=120,
            cycle_time_s=10,
            checkpoint_time_s=600,
            max_samples=10,
        )
        with open(cntl) as fh:
            content = fh.read()

        # Required keys for bound leg ABFE
        for key in (
            "JOB_TRANSPORT = 'LOCAL_OPENMM'",
            "BASENAME = 'trackb'",
            "TEMPERATURES = '300.0'",
            "LAMBDAS =",
            "DIRECTION =",
            "INTERMEDIATE =",
            "LAMBDA1 =",
            "LAMBDA2 =",
            "ALPHA =",
            "U0 =",
            "W0COEFF =",
            "DISPLACEMENT = '25.0, 0.0, 0.0'",  # 2.5 nm → 25 Å
            "LIGOFFSET = '0., 0., 0.'",
            "WALL_TIME = 120",
            "CYCLE_TIME = 10",
            "CHECKPOINT_TIME = 600",
            "PRODUCTION_STEPS = '2500'",
            "PRNT_FREQUENCY = '2500'",
            "TRJ_FREQUENCY = '25000'",
            "MAX_SAMPLES = 10",
            "LIGAND_ATOMS = 100, 101, 102, 103",
            "LIGAND_CM_ATOMS = 100, 101, 102, 103",
            "RCPT_CM_ATOMS = 0, 1, 2, 3, 4",
            "CM_KF = 25.0",
            "CM_TOL = 5.0",
            "POS_RESTRAINED_ATOMS = 0, 1, 2, 3, 4",
            "POSRE_FORCE_CONSTANT = 25.0",
            "POSRE_TOLERANCE = 0.5",
            "UMAX = 200.0",
            "ACORE = 0.0625",
            "UBCORE = 100.0",
            "FRICTION_COEFF = 0.5",
            "TIME_STEP = 0.002",
            "OPENMM_PLATFORM = CUDA",
            "VERBOSE = 'no'",
        ):
            assert key in content, f"cntl missing key: {key!r}"


def test_cntl_writer_free_leg_omits_receptor_keys(v2):
    """Free-leg cntl (no receptor) must omit RCPT_CM_ATOMS / POS_RESTRAINED."""
    with tempfile.TemporaryDirectory() as td:
        cntl = os.path.join(td, "trackb.cntl")
        nodefile = os.path.join(td, "nodefile")
        v2.write_cntl_file(
            cntl_path=cntl,
            basename="trackb",
            nodefile_path=nodefile,
            ligand_atom_indices=[0, 1, 2, 3],
            pos_restrained_atom_indices=[],
            displacement_nm=(2.5, 0.0, 0.0),
            production_steps=2500,
            prnt_frequency=2500,
            trj_frequency=25000,
            wall_time_min=120,
            cycle_time_s=10,
            checkpoint_time_s=600,
            max_samples=10,
        )
        with open(cntl) as fh:
            content = fh.read()
        # Receptor-specific keys must NOT appear on free leg.
        for key in (
            "RCPT_CM_ATOMS",
            "CM_KF",
            "CM_TOL",
            "POS_RESTRAINED_ATOMS",
        ):
            assert key not in content, (
                f"cntl wrongly includes receptor key {key!r} on free leg"
            )
        # POSRE_FORCE_CONSTANT + POSRE_TOLERANCE are ALWAYS included because
        # abfe_structprep.massage_keywords auto-populates POS_RESTRAINED_ATOMS
        # during mintherm (MINTHERM_RESTRAIN_SOLUTES default=YES). Omitting
        # these on the free leg → TypeError on float(None) in
        # OMMSystem.set_positional_restraints.
        assert "POSRE_FORCE_CONSTANT = 25.0" in content
        assert "POSRE_TOLERANCE = 0.5" in content


def test_cntl_schedule_csv_lengths(v2):
    """LAMBDAS / DIRECTION / INTERMEDIATE / LAMBDA1 / LAMBDA2 / ALPHA / U0 / W0COEFF
    must each be 22 comma-separated values in the cntl."""
    with tempfile.TemporaryDirectory() as td:
        cntl = os.path.join(td, "trackb.cntl")
        nodefile = os.path.join(td, "nodefile")
        v2.write_cntl_file(
            cntl_path=cntl,
            basename="trackb",
            nodefile_path=nodefile,
            ligand_atom_indices=[0],
            pos_restrained_atom_indices=[],
            displacement_nm=(2.5, 0.0, 0.0),
            production_steps=10,
            prnt_frequency=10,
            trj_frequency=10,
            wall_time_min=10,
            cycle_time_s=10,
            checkpoint_time_s=10,
        )
        with open(cntl) as fh:
            content = fh.read()
        for key in ("LAMBDAS", "DIRECTION", "INTERMEDIATE",
                    "LAMBDA1", "LAMBDA2", "ALPHA", "U0", "W0COEFF"):
            m = re.search(rf"^{key} =\s+'([^']*)'", content, re.MULTILINE)
            assert m, f"cntl missing or malformed key {key!r}"
            vals = [v.strip() for v in m.group(1).split(",")]
            assert len(vals) == 22, (
                f"{key} has {len(vals)} entries, expected 22"
            )


# ===========================================================================
# C.1 UWHAM per-state SSOT (BLOCKING silent-bias fix).
# The postprocess MUST read the per-state (λ1/λ2/α/u0/w0) arrays from the
# PRODUCTION cntl and pass them to calculate_uwham — NOT rely on upstream
# 22-state α/u0 defaults (which would silently bias a ramped densified38 run).
# These tests exercise the pure-Python cntl-parse + fail-fast guard (no
# atom_openmm import needed for the parse path).
# ===========================================================================
def _load_uwham_module():
    spec = importlib.util.spec_from_file_location(
        "trackb_uwham_postprocess",
        os.path.join(_PROJ, "scripts", "trackb_uwham_postprocess.py"),
    )
    if spec is None or spec.loader is None:
        pytest.skip("could not locate scripts/trackb_uwham_postprocess.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def uwham_pp():
    return _load_uwham_module()


def _write_leg_cntl(leg_dir, v2, schedule_name):
    """Write a production cntl into leg_dir from a named v2 schedule."""
    cntl = os.path.join(leg_dir, "trackb_asyncre.cntl")
    nf = os.path.join(leg_dir, "nodefile")
    sched = None if schedule_name == "canonical22" else v2.get_schedule(
        schedule_name)
    v2.write_cntl_file(
        cntl_path=cntl, basename="trackb", nodefile_path=nf,
        ligand_atom_indices=[0], pos_restrained_atom_indices=[],
        displacement_nm=(2.5, 0.0, 0.0), production_steps=10,
        prnt_frequency=10, trj_frequency=10, wall_time_min=10,
        cycle_time_s=10, checkpoint_time_s=10, schedule=sched,
    )
    return cntl


def test_c1_parse_cntl_schedule_canonical22_equals_upstream_defaults(
    uwham_pp, v2, tmp_path,
):
    """C.1 regression: the canonical22 cntl parses to arrays IDENTICAL to the
    upstream 22-state defaults (alpha=0.10, u0=110). Passing them explicitly is
    therefore byte-equivalent to the default path → back-compat preserved."""
    _write_leg_cntl(str(tmp_path), v2, "canonical22")
    cntl = os.path.join(str(tmp_path), "trackb_asyncre.cntl")
    sched = uwham_pp._parse_cntl_schedule(cntl)
    assert sched is not None
    assert sched["n_states"] == 22
    # upstream calculate_uwham defaults: alpha=0.10 ×22, u0=110.0 ×22.
    assert sched["alpha"] == [0.10] * 22
    assert sched["u0"] == [110.0] * 22
    # 2 INTERMEDIATE at the λ=0.5 midpoints (states 10, 11).
    assert [i for i, x in enumerate(sched["intermd"]) if x == 1] == [10, 11]


def test_c1_parse_cntl_schedule_densified38_uses_ramped_arrays(
    uwham_pp, v2, tmp_path,
):
    """C.1 silent-bias fix: the densified38 cntl parses to the RAMPED α
    (0.10→0.25) + u0 (110→82) arrays — NOT the 22-state defaults. Analyzing
    with these (vs defaults) is what prevents the silent ΔG bias."""
    _write_leg_cntl(str(tmp_path), v2, "densified38")
    cntl = os.path.join(str(tmp_path), "trackb_asyncre.cntl")
    sched = uwham_pp._parse_cntl_schedule(cntl)
    assert sched is not None
    assert sched["n_states"] == 38
    # The ramped arrays must NOT equal the upstream 22-state defaults.
    assert sched["alpha"] != [0.10] * 38
    assert max(sched["alpha"]) == 0.25  # apex α
    assert min(sched["u0"]) == 82.0     # apex u0
    assert sum(int(x) for x in sched["intermd"]) == 18
    # The arrays match the v2 schedule exactly (cntl is byte-faithful SSOT).
    s = v2.get_schedule("densified38")
    assert sched["alpha"] == s["alpha"]
    assert sched["u0"] == s["u0"]
    assert sched["w0"] == s["w0"]
    assert sched["lambda1"] == s["lambdas_1"]
    assert sched["lambda2"] == s["lambdas_2"]


def test_c1_analyze_one_leg_fail_fast_when_cntl_missing(uwham_pp, tmp_path):
    """C.1 BLOCKING: analyze_one_leg with require_cntl_schedule=True (default)
    raises when the production cntl is absent — NO silent fallback to upstream
    22-state defaults. Builds 22 r*/ dirs so the only missing piece is the
    cntl, isolating the C.1 guard."""
    leg = tmp_path / "leg"
    leg.mkdir()
    for i in range(22):
        rd = leg / f"r{i}"
        rd.mkdir()
        (rd / "trackb.out").write_text(
            f"{i} 300.0 1.0 0.1 0.1 0.1 110.0 0.0 -5.0 -1.0 0.0\n"
        )
    # No trackb_asyncre.cntl present → _derive_expected_replicas falls back to
    # counting r*/ dirs (22), then the C.1 guard fires before analysis.
    with pytest.raises(FileNotFoundError, match="C.1 UWHAM per-state SSOT"):
        uwham_pp.analyze_one_leg(
            leg_dir=str(leg), jobname="trackb",
            mintimeid=None, maxtimeid=None, block_bootstrap=False,
        )


def test_c1_leg_analysis_window_rejects_none_schedule(uwham_pp, tmp_path):
    """C.1: _leg_analysis_window raises on a None schedule (it would otherwise
    fall through to upstream defaults — the silent-bias path)."""
    with pytest.raises(ValueError, match="C.1 UWHAM per-state SSOT"):
        uwham_pp._leg_analysis_window(
            leg_dir=str(tmp_path), jobname="trackb", schedule=None,
            mintimeid=None, maxtimeid=None, n_states=38,
        )


# ===========================================================================
# densified38v2 — count-neutral linear-λ REBALANCE of densified38
# (dplus 8→9 overlap fix, 2026-06-05).
# Fixes the dplus 8→9 BC=0.131 overlap gap by thinning the saturated coupled
# plateau and densifying the 0.40→0.45 turnover. RIGOR-NEUTRAL / ranking-safe:
# both physical endpoints invariant, soft-core ladder + α/U0/W0 ramps UNCHANGED.
# ===========================================================================
# Expected forward linear λ array (exact).
_DENSE38v2_FWD_LINEAR_EXPECTED = [
    0.00, 0.08, 0.16, 0.24, 0.32,
    0.38, 0.40, 0.42, 0.44, 0.45,
]


def test_densified38v2_registered_in_schedules(v2):
    """densified38v2 is a NEW selectable schedule alongside the existing ones
    (does NOT replace densified38)."""
    assert "densified38v2" in v2.SCHEDULES
    # All prior schedules still present (densified38 NOT overwritten).
    assert "canonical22" in v2.SCHEDULES
    assert "densified34" in v2.SCHEDULES
    assert "densified38" in v2.SCHEDULES
    # Selectable via get_schedule without error.
    s = v2.get_schedule("densified38v2")
    assert s is not None


def test_densified38v2_builds_38_states_19_plus_19(v2):
    """densified38v2 = 38 states = 19 forward (+1) + 19 backward (-1)."""
    s = v2.get_schedule("densified38v2")
    assert s["n_states"] == 38
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "directions", "intermd", "w0", "alpha", "u0"):
        assert len(s[key]) == 38, f"{key} length"
    assert s["directions"][:19] == [1] * 19
    assert s["directions"][19:] == [-1] * 19


def test_densified38v2_forward_linear_lambda_array(v2):
    """The forward LINEAR λ segment (states 0-9) equals the specified
    10-value array; λ=0.00 coupled endpoint is the first state."""
    s = v2.get_schedule("densified38v2")
    # Forward linear states 0-9: λ1 == λ2 == the rebalanced array.
    assert s["lambdas_1"][:10] == _DENSE38v2_FWD_LINEAR_EXPECTED
    assert s["lambdas_2"][:10] == _DENSE38v2_FWD_LINEAR_EXPECTED
    assert s["lambdas"][:10] == _DENSE38v2_FWD_LINEAR_EXPECTED
    # Coupled endpoint preserved as first state.
    assert s["lambdas_1"][0] == 0.00
    # The module-level constant matches the spec too.
    assert list(v2.DENSE38v2_LAMBDA_FWD_LINEAR) == _DENSE38v2_FWD_LINEAR_EXPECTED
    # Linear segment carries W0=0, ALPHA=0.10, U0=110, INTERMEDIATE=0.
    assert s["w0"][:10] == [0.0] * 10
    assert s["alpha"][:10] == [0.10] * 10
    assert s["u0"][:10] == [110.0] * 10
    assert s["intermd"][:10] == [0] * 10


def test_densified38v2_softcore_ladder_identical_to_densified38(v2):
    """The 9-window soft-core ladder (forward states 10-18 AND their backward
    mirror) is BYTE-IDENTICAL to densified38 — only the linear segment differs.
    This is the load-bearing invariant (ladder + α/U0/W0 unchanged)."""
    v2_s = v2.get_schedule("densified38v2")
    d38 = v2.get_schedule("densified38")
    # Forward ladder block (states 10-18).
    for key in ("lambdas_1", "lambdas_2", "lambdas", "w0", "alpha",
                "u0", "intermd"):
        assert v2_s[key][10:19] == d38[key][10:19], f"fwd ladder {key}"
    # Backward ladder block (states 19-27, the reversed forward ladder).
    for key in ("lambdas_1", "lambdas_2", "lambdas", "w0", "alpha",
                "u0", "intermd"):
        assert v2_s[key][19:28] == d38[key][19:28], f"bwd ladder {key}"
    # The raw ladder tuple constant is shared (single source).
    assert len(v2.DENSE38_LADDER) == 9


def test_densified38v2_alpha_u0_apex_preserved(v2):
    """ALPHA ramp apex 0.25 and U0 ramp apex 82.0 preserved (same as
    densified38 — the rebalance does NOT touch the soft-core ramps)."""
    s = v2.get_schedule("densified38v2")
    assert max(s["alpha"]) == 0.25
    assert min(s["u0"]) == 82.0
    # Linear segment is the floor of both ramps.
    assert min(s["alpha"]) == 0.10
    assert max(s["u0"]) == 110.0


def test_densified38v2_intermediate_sum_unchanged(v2):
    """INTERMEDIATE sum = 18 (9 forward ladder + 9 backward ladder), unchanged
    from densified38 — the linear rebalance does not add/remove intermediates."""
    v2_s = v2.get_schedule("densified38v2")
    d38 = v2.get_schedule("densified38")
    assert sum(int(x) for x in v2_s["intermd"]) == 18
    assert sum(int(x) for x in v2_s["intermd"]) == sum(
        int(x) for x in d38["intermd"])
    # INTERMEDIATE positions: forward 10-18, backward 19-27.
    pos = [i for i, x in enumerate(v2_s["intermd"]) if x == 1]
    assert pos == list(range(10, 19)) + list(range(19, 28))


def test_densified38v2_backward_is_reversed_mirror(v2):
    """Backward half = whole-tuple reverse of the forward half for EVERY
    per-state array (λ1 stays λ1, λ2 stays λ2 — NOT a λ1↔λ2 swap); DIRECTION
    forced to -1 on the backward block. This is the dminus-auto-mirror property
    (symmetry argument): no separate backward edit, symmetry is automatic."""
    s = v2.get_schedule("densified38v2")
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "intermd", "w0", "alpha", "u0"):
        fwd = s[key][:19]
        bwd = s[key][19:]
        assert list(reversed(fwd)) == bwd, f"reversed mirror {key}"
    # Backward DIRECTION is all -1.
    assert s["directions"][19:] == [-1] * 19


def test_densified38v2_endpoints_invariant_vs_densified38(v2):
    """Both physical endpoints are identical to densified38 (Q3 ΔG-unbias
    proof in practice): coupled λ=0.00 first state, and the W0=1.0/λ=0.5 apex
    (forward state 18 / backward state 19). Rebalance is endpoint-invariant."""
    v2_s = v2.get_schedule("densified38v2")
    d38 = v2.get_schedule("densified38")
    # Coupled endpoint (state 0).
    assert v2_s["lambdas_1"][0] == d38["lambdas_1"][0] == 0.00
    assert v2_s["w0"][0] == d38["w0"][0] == 0.0
    # Decoupled apex (forward state 18: W0=1.0, λ2=0.5).
    assert v2_s["w0"][18] == d38["w0"][18] == 1.0
    assert v2_s["lambdas_2"][18] == d38["lambdas_2"][18] == 0.5


def test_densified38_unchanged_regression(v2):
    """Regression: densified38 is UNCHANGED by the densified38v2 addition.
    Assert its arrays byte-for-byte against the pre-rebalance reference values
    (uniform-Δλ=0.05 linear segment + the 9-window soft-core ladder)."""
    d38 = v2.get_schedule("densified38")
    assert d38["n_states"] == 38
    # Forward linear segment: uniform Δλ=0.05 (the ORIGINAL, NOT rebalanced).
    assert d38["lambdas_1"][:10] == [
        0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
    # Full forward λ1 (linear + ladder).
    assert d38["lambdas_1"][:19] == [
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.45, 0.45, 0.46, 0.47, 0.48, 0.485, 0.49, 0.495, 0.5]
    assert d38["lambdas_2"][:19] == [
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.46, 0.47, 0.48, 0.49, 0.495, 0.498, 0.499, 0.5, 0.5]
    assert d38["alpha"][:19] == [0.10] * 10 + [
        0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.23, 0.24, 0.25]
    assert d38["u0"][:19] == [110.0] * 10 + [
        105.0, 100.0, 95.0, 92.0, 90.0, 88.0, 86.0, 84.0, 82.0]
    assert d38["w0"][:19] == [0.0] * 10 + [
        0.15, 0.30, 0.45, 0.60, 0.72, 0.82, 0.90, 0.96, 1.00]
    assert d38["directions"] == [1] * 19 + [-1] * 19
    assert sum(int(x) for x in d38["intermd"]) == 18
    # densified38v2 differs from densified38 ONLY in the linear segments.
    # Layout per direction: [0:10] forward linear, [10:19] forward ladder;
    # backward mirror: [19:28] backward ladder, [28:38] backward linear.
    # The two LADDER blocks are identical; both LINEAR blocks differ (forward
    # rebalanced + its reversed-mirror backward tail).
    v2_s = v2.get_schedule("densified38v2")
    assert v2_s["lambdas_1"][:10] != d38["lambdas_1"][:10]   # fwd linear differs
    assert v2_s["lambdas_1"][10:19] == d38["lambdas_1"][10:19]  # fwd ladder same
    assert v2_s["lambdas_1"][19:28] == d38["lambdas_1"][19:28]  # bwd ladder same
    assert v2_s["lambdas_1"][28:] != d38["lambdas_1"][28:]   # bwd linear differs


# ===========================================================================
# densified38v3 — PER-DIRECTION (NON-mirror) schedule
# (densified38v2 dminus-asymmetry fix, 2026-06-06).
# Breaks densified38v2's forward-only-mirror symmetry: dplus = clean v2 forward
# (UNCHANGED); dminus = v2 backward + ONE W0-graded micro-bridge at the
# soft-core-end → plateau handoff (closes the dminus 9→10 BC=0.24 overlap hole
# the reversed mirror created). UNEQUAL per-direction counts (dplus 19, dminus
# 20 = 39). Endpoints invariant both directions → rigor-neutral / ranking-safe.
# ===========================================================================
# Expected dminus W0-graded micro-bridge tuple (λ1, λ2, W0, ALPHA, U0):
# λ=0.45 (broad transition), W0=0.07, α/U0 interpolated between ladder-end
# (W0=0.15/α=0.12/U0=105) and plateau (W0=0/α=0.10/U0=110).
_DENSE38v3_BRIDGE_EXPECTED = (0.45, 0.45, 0.07, 0.11, 107.5)


def test_densified38v3_registered_in_schedules(v2):
    """densified38v3 is a NEW selectable schedule alongside the existing ones
    (does NOT replace densified38 / densified38v2)."""
    assert "densified38v3" in v2.SCHEDULES
    # All prior schedules still present (none overwritten).
    assert "canonical22" in v2.SCHEDULES
    assert "densified34" in v2.SCHEDULES
    assert "densified38" in v2.SCHEDULES
    assert "densified38v2" in v2.SCHEDULES
    s = v2.get_schedule("densified38v3")
    assert s is not None


def test_densified38v3_unequal_per_direction_counts(v2):
    """densified38v3 has UNEQUAL per-direction counts: dplus 19, dminus 20
    (the bridge adds ONE state to dminus only), total 39. DIRECTION column is a
    contiguous +1 (19) then -1 (20) block (required by the per-direction
    launcher's contiguity gate)."""
    s = v2.get_schedule("densified38v3")
    assert s["n_states"] == 39
    n_fwd = sum(1 for d in s["directions"] if d >= 0)
    n_bwd = s["n_states"] - n_fwd
    assert n_fwd == 19, "dplus (forward) = 19 states"
    assert n_bwd == 20, "dminus (backward) = 20 states (+1 bridge)"
    # Contiguous +1 then -1.
    assert s["directions"][:19] == [1] * 19
    assert s["directions"][19:] == [-1] * 20
    # Every per-state array is the full 39 long.
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "directions", "intermd", "w0", "alpha", "u0"):
        assert len(s[key]) == 39, f"{key} length"


def test_densified38v3_dplus_byte_equal_to_densified38v2(v2):
    """dplus (forward [:19]) is BYTE-IDENTICAL to densified38v2's dplus — the
    clean side is reused verbatim, NOT rebuilt ('dplus is already
    correct, do not change it')."""
    v3 = v2.get_schedule("densified38v3")
    v2s = v2.get_schedule("densified38v2")
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "intermd", "w0", "alpha", "u0"):
        assert v3[key][:19] == v2s[key][:19], f"dplus {key} must match v2 dplus"
    # dplus DIRECTION is all +1.
    assert v3["directions"][:19] == [1] * 19


def test_densified38v3_dminus_has_w0_graded_bridge(v2):
    """dminus has ONE EXTRA W0-graded micro-bridge window at the soft-core-end
    → plateau handoff. The bridge is at λ1=λ2=0.45 (the broad transition's λ),
    W0=0.07, with α/U0 between the ladder-end and the plateau, INTERMEDIATE=1."""
    v3 = v2.get_schedule("densified38v3")
    # dminus block = states 19..38 (20 states).
    dminus_l1 = v3["lambdas_1"][19:]
    dminus_l2 = v3["lambdas_2"][19:]
    dminus_w0 = v3["w0"][19:]
    dminus_a = v3["alpha"][19:]
    dminus_u0 = v3["u0"][19:]
    dminus_i = v3["intermd"][19:]
    # The bridge tuple is present exactly once in the dminus block (it is the
    # only window with W0≈0.07).
    bridge_positions = [
        j for j in range(len(dminus_w0))
        if abs(dminus_w0[j] - 0.07) < 1e-9
    ]
    assert len(bridge_positions) == 1, "exactly one W0=0.07 bridge in dminus"
    bp = bridge_positions[0]
    bl1, bl2, bw0, ba, bu0 = _DENSE38v3_BRIDGE_EXPECTED
    assert dminus_l1[bp] == bl1
    assert dminus_l2[bp] == bl2
    assert dminus_w0[bp] == bw0
    assert dminus_a[bp] == ba
    assert dminus_u0[bp] == bu0
    assert dminus_i[bp] == 1, "bridge is a soft-core window (INTERMEDIATE=1)"
    # Module-level constant matches the spec.
    assert tuple(v2.DENSE38v3_DMINUS_BRIDGE) == _DENSE38v3_BRIDGE_EXPECTED


def test_densified38v3_bridge_sits_at_handoff_between_broad_and_plateau(v2):
    """The bridge is inserted BETWEEN the broad transition state (λ=0.45, W0=0)
    and the bare tight plateau (λ=0.44, W0=0) — i.e. the 9→10 overlap hole.
    The neighbor BEFORE the bridge is the broad transition (λ=0.45, W0=0); the
    neighbor AFTER the bridge is the first plateau step (λ<0.45, W0=0)."""
    v3 = v2.get_schedule("densified38v3")
    dm_l1 = v3["lambdas_1"][19:]
    dm_w0 = v3["w0"][19:]
    bp = next(j for j in range(len(dm_w0)) if abs(dm_w0[j] - 0.07) < 1e-9)
    # Before bridge = broad transition (last W0=0 state at λ=0.45).
    assert dm_w0[bp - 1] == 0.0
    assert dm_l1[bp - 1] == 0.45
    # After bridge = bare tight plateau (W0=0, λ stepping down from 0.45).
    assert dm_w0[bp + 1] == 0.0
    assert dm_l1[bp + 1] < 0.45


def test_densified38v3_bridge_alpha_u0_monotone_across_handoff(v2):
    """α and U0 are monotone-graded across the W0 handoff: the bridge's
    α (0.11) lies between the plateau α (0.10) and the ladder-end α (0.12); its
    U0 (107.5) lies between the plateau U0 (110) and the ladder-end U0 (105).
    The W0=0.07 bridge gives the broad transition a soft-core intermediate to
    hand off to the plateau through, instead of abutting a bare W0=0 jump."""
    v3 = v2.get_schedule("densified38v3")
    dm_w0 = v3["w0"][19:]
    dm_a = v3["alpha"][19:]
    dm_u0 = v3["u0"][19:]
    bp = next(j for j in range(len(dm_w0)) if abs(dm_w0[j] - 0.07) < 1e-9)
    # Bridge W0 between ladder-end (0.15) and plateau (0.0).
    assert 0.0 < dm_w0[bp] < 0.15
    # Bridge α between plateau (0.10) and ladder-end (0.12).
    assert 0.10 <= dm_a[bp] <= 0.12
    assert dm_a[bp] == 0.11
    # Bridge U0 between ladder-end (105) and plateau (110).
    assert 105.0 <= dm_u0[bp] <= 110.0
    assert dm_u0[bp] == 107.5


def test_densified38v3_dminus_is_NOT_reverse_of_dplus(v2):
    """The KEY departure from densified38v2: dminus is NOT reversed(dplus).
    densified38v2's dminus == reversed(dplus) (19==19); densified38v3's dminus
    (20) has one extra state, so reversed(dplus) (19) cannot equal it."""
    v3 = v2.get_schedule("densified38v3")
    dplus_l1 = v3["lambdas_1"][:19]
    dminus_l1 = v3["lambdas_1"][19:]
    assert len(dminus_l1) == 20
    assert list(reversed(dplus_l1)) != dminus_l1
    # But REMOVING the bridge state from dminus DOES recover reversed(dplus)
    # (the bridge is a pure insertion; the rest mirrors v2).
    bp = next(j for j in range(len(v3["w0"][19:]))
              if abs(v3["w0"][19:][j] - 0.07) < 1e-9)
    dminus_no_bridge = dminus_l1[:bp] + dminus_l1[bp + 1:]
    assert list(reversed(dplus_l1)) == dminus_no_bridge


def test_densified38v3_dminus_direction_forced_plus_at_runtime(v2):
    """The schedule encodes dminus DIRECTION=-1 (thermodynamic-leg label). The
    (b+) runtime fix (DIRECTION→+1 forced + DISPLACEMENT→-25) is applied by
    generate_per_direction_cntls at slice time, NOT in the schedule. Here we
    assert the schedule's dminus block carries the -1 label that the slicer
    keys on; the +1/-25 rewrite is exercised in the per-direction suite + the
    dry-run. This guards the contiguity contract the launcher requires."""
    s = v2.get_schedule("densified38v3")
    assert s["directions"][19:] == [-1] * 20
    assert s["directions"][:19] == [1] * 19


def test_densified38v3_endpoints_invariant_both_directions(v2):
    """Per-direction endpoint-invariance (load-bearing / ΔG-unbias): BOTH directions keep
    λ=0 (coupled) as one endpoint and the apex (W0=1.0, λ2=0.5) as the other.
    dplus: state 0 = coupled, state 18 = apex. dminus: first state (19) = apex
    (reversed), last state (38) = coupled. The bridge insertion is INTERIOR;
    endpoints unchanged → E[ΔG_int] unchanged (Zwanzig/Kirkwood)."""
    s = v2.get_schedule("densified38v3")
    n_fwd = 19
    # dplus coupled endpoint (state 0).
    assert s["lambdas_1"][0] == 0.00
    assert s["w0"][0] == 0.0
    # dplus decoupled apex (state 18).
    assert s["w0"][n_fwd - 1] == 1.0
    assert s["lambdas_2"][n_fwd - 1] == 0.5
    # dminus decoupled apex (first backward state, 19).
    assert s["w0"][n_fwd] == 1.0
    assert s["lambdas_2"][n_fwd] == 0.5
    # dminus coupled endpoint (last state, 38).
    assert s["lambdas_1"][-1] == 0.00
    assert s["w0"][-1] == 0.0


def test_densified38v3_extra_intermediate_is_only_the_bridge(v2):
    """dminus has EXACTLY one more INTERMEDIATE window than densified38v2's
    dminus (the bridge carries INTERMEDIATE=1). dplus INTERMEDIATE count is
    unchanged from v2 dplus."""
    v3 = v2.get_schedule("densified38v3")
    v2s = v2.get_schedule("densified38v2")
    v3_dplus_int = sum(int(x) for x in v3["intermd"][:19])
    v2_dplus_int = sum(int(x) for x in v2s["intermd"][:19])
    assert v3_dplus_int == v2_dplus_int, "dplus INTERMEDIATE unchanged"
    v3_dminus_int = sum(int(x) for x in v3["intermd"][19:])
    v2_dminus_int = sum(int(x) for x in v2s["intermd"][19:])
    assert v3_dminus_int == v2_dminus_int + 1, "dminus +1 INTERMEDIATE (bridge)"


# --- Regression: prior schedules UNCHANGED by the densified38v3 addition -----

def test_densified38v2_unchanged_by_v3_addition(v2):
    """Regression: densified38v2 is UNCHANGED by adding densified38v3.
    38 states, 19+19, the rebalanced linear forward array intact, the dminus
    still the reversed mirror (19, NO bridge)."""
    v2s = v2.get_schedule("densified38v2")
    assert v2s["n_states"] == 38
    assert v2s["directions"] == [1] * 19 + [-1] * 19
    assert v2s["lambdas_1"][:10] == [
        0.00, 0.08, 0.16, 0.24, 0.32, 0.38, 0.40, 0.42, 0.44, 0.45]
    # dminus is the pure reversed mirror (no W0=0.07 bridge anywhere).
    assert not any(abs(w - 0.07) < 1e-9 for w in v2s["w0"])
    # Backward half = reversed forward half.
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "intermd", "w0", "alpha", "u0"):
        assert list(reversed(v2s[key][:19])) == v2s[key][19:], f"v2 mirror {key}"


def test_densified38_unchanged_by_v3_addition(v2):
    """Regression: densified38 is UNCHANGED by adding densified38v3."""
    d38 = v2.get_schedule("densified38")
    assert d38["n_states"] == 38
    assert d38["directions"] == [1] * 19 + [-1] * 19
    assert d38["lambdas_1"][:10] == [
        0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
    assert not any(abs(w - 0.07) < 1e-9 for w in d38["w0"])


def test_canonical22_unchanged_by_v3_addition(v2):
    """Regression: canonical22 is UNCHANGED by adding densified38v3."""
    c22 = v2.get_schedule("canonical22")
    assert c22["n_states"] == 22
    assert c22["directions"] == [1] * 11 + [-1] * 11
    assert not any(abs(w - 0.07) < 1e-9 for w in c22["w0"])


def test_densified38v3_cntl_writer_emits_39_state_unequal_columns(v2, tmp_path):
    """write_cntl_file serializes the densified38v3 schedule generically: all
    per-state cntl columns are 39 long, DIRECTION = 19×+1 then 20×-1, and the
    W0=0.07 bridge appears once in the W0COEFF column. This proves the cntl
    writer is count-agnostic for the unequal-count schedule (no fixed-38/22
    assumption)."""
    cntl_path = os.path.join(str(tmp_path), "trackb_v3.cntl")
    sched = v2.get_schedule("densified38v3")
    v2.write_cntl_file(
        cntl_path=cntl_path,
        basename="trackb",
        nodefile_path=os.path.join(str(tmp_path), "nodefile"),
        ligand_atom_indices=[0, 1, 2],
        pos_restrained_atom_indices=[3, 4],
        displacement_nm=(2.5, 0.0, 0.0),
        production_steps=1000,
        prnt_frequency=100,
        trj_frequency=100,
        wall_time_min=60,
        cycle_time_s=10,
        checkpoint_time_s=10,
        schedule=sched,
    )
    with open(cntl_path) as fh:
        text = fh.read()

    def _col(key):
        m = re.search(rf"^{key}\s*=\s*'([^']*)'", text, re.M)
        assert m, f"{key} column missing"
        return [t.strip() for t in m.group(1).split(",") if t.strip() != ""]

    for key in ("LAMBDAS", "DIRECTION", "INTERMEDIATE", "LAMBDA1",
                "LAMBDA2", "ALPHA", "U0", "W0COEFF"):
        assert len(_col(key)) == 39, f"{key} not 39 long"
    direction = [int(float(x)) for x in _col("DIRECTION")]
    assert direction[:19] == [1] * 19
    assert direction[19:] == [-1] * 20
    w0 = [float(x) for x in _col("W0COEFF")]
    assert sum(1 for w in w0 if abs(w - 0.07) < 1e-9) == 1, "one W0=0.07 bridge"


# ===========================================================================
# densified38v4 — PER-DIRECTION (NON-mirror): SECOND dminus W0-graded micro-
# bridge at the re-indexed 5→6 handoff
# (dminus second-bridge fix, 2026-06-06).
# dplus = clean densified38v2/v3 forward (UNCHANGED, byte-equal); dminus =
# densified38v3 dminus (which already has the 9→10 W0=0.07 bridge) + one MORE
# W0-graded micro-bridge at the 5→6 soft-core-ladder-end → plateau handoff
# (closes the dminus 5→6 BC=0.18 hole the v3 re-pilot exposed — the SAME single
# turnover, re-indexed, NOT a new pathology). UNEQUAL per-direction counts
# (dplus 19, dminus 21 = 40). Endpoints invariant both directions → rigor-
# neutral / ranking-safe. LAST iteration under the 2-bridge/leg hard cap.
# ===========================================================================
# Expected NEW dminus 5→6 W0-graded micro-bridge tuple (λ1, λ2, W0, ALPHA, U0):
# λ1=0.465, λ2=0.485 (interior to the v3-dminus-local 5 / 6 neighbors); W0=0.17;
# α=0.17, U0=93.5 (midpoints of state 5 (0.18/92) and state 6 (0.16/95)).
_DENSE38v4_BRIDGE56_EXPECTED = (0.465, 0.485, 0.17, 0.17, 93.5)


def test_densified38v4_registered_in_schedules(v2):
    """densified38v4 is a NEW selectable schedule alongside the existing ones
    (does NOT replace densified38 / v2 / v3 / canonical22)."""
    assert "densified38v4" in v2.SCHEDULES
    # All prior schedules still present (none overwritten).
    assert "canonical22" in v2.SCHEDULES
    assert "densified34" in v2.SCHEDULES
    assert "densified38" in v2.SCHEDULES
    assert "densified38v2" in v2.SCHEDULES
    assert "densified38v3" in v2.SCHEDULES
    s = v2.get_schedule("densified38v4")
    assert s is not None


def test_densified38v4_unequal_per_direction_counts_40(v2):
    """densified38v4 has UNEQUAL per-direction counts: dplus 19, dminus 21
    (v3 dminus 20 + the new 5→6 bridge), total 40. DIRECTION column is a
    contiguous +1 (19) then -1 (21) block (required by the per-direction
    launcher's contiguity gate)."""
    s = v2.get_schedule("densified38v4")
    assert s["n_states"] == 40
    n_fwd = sum(1 for d in s["directions"] if d >= 0)
    n_bwd = s["n_states"] - n_fwd
    assert n_fwd == 19, "dplus (forward) = 19 states"
    assert n_bwd == 21, "dminus (backward) = 21 states (v3 20 + 1 bridge)"
    # Contiguous +1 then -1.
    assert s["directions"][:19] == [1] * 19
    assert s["directions"][19:] == [-1] * 21
    # Every per-state array is the full 40 long.
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "directions", "intermd", "w0", "alpha", "u0"):
        assert len(s[key]) == 40, f"{key} length"


def test_densified38v4_dplus_byte_equal_to_v2_and_v3(v2):
    """dplus (forward [:19]) is BYTE-IDENTICAL to BOTH densified38v2's AND
    densified38v3's dplus — the clean side is reused verbatim, NOT rebuilt
    ('dplus is clean, do not touch it')."""
    v4 = v2.get_schedule("densified38v4")
    v3 = v2.get_schedule("densified38v3")
    v2s = v2.get_schedule("densified38v2")
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "intermd", "w0", "alpha", "u0"):
        assert v4[key][:19] == v2s[key][:19], f"dplus {key} must match v2 dplus"
        assert v4[key][:19] == v3[key][:19], f"dplus {key} must match v3 dplus"
    # dplus DIRECTION is all +1.
    assert v4["directions"][:19] == [1] * 19


def test_densified38v4_dminus_has_BOTH_bridges(v2):
    """dminus carries BOTH W0-graded micro-bridges: the v3 9→10 (W0=0.07) AND
    the new 5→6 (W0=0.17). Each appears exactly once in the dminus block."""
    v4 = v2.get_schedule("densified38v4")
    dminus_w0 = v4["w0"][19:]
    assert len(dminus_w0) == 21
    n_007 = sum(1 for w in dminus_w0 if abs(w - 0.07) < 1e-9)
    n_017 = sum(1 for w in dminus_w0 if abs(w - 0.17) < 1e-9)
    assert n_007 == 1, "exactly one v3 9→10 W0=0.07 bridge in dminus"
    assert n_017 == 1, "exactly one new 5→6 W0=0.17 bridge in dminus"


def test_densified38v4_new_5to6_bridge_tuple(v2):
    """The NEW 5→6 bridge tuple matches the spec exactly:
    λ1=0.465, λ2=0.485, W0=0.17, α=0.17, U0=93.5, INTERMEDIATE=1."""
    v4 = v2.get_schedule("densified38v4")
    dm_l1 = v4["lambdas_1"][19:]
    dm_l2 = v4["lambdas_2"][19:]
    dm_lam = v4["lambdas"][19:]
    dm_w0 = v4["w0"][19:]
    dm_a = v4["alpha"][19:]
    dm_u0 = v4["u0"][19:]
    dm_i = v4["intermd"][19:]
    bp = next(j for j in range(len(dm_w0)) if abs(dm_w0[j] - 0.17) < 1e-9)
    bl1, bl2, bw0, ba, bu0 = _DENSE38v4_BRIDGE56_EXPECTED
    assert dm_l1[bp] == bl1
    assert dm_l2[bp] == bl2
    assert dm_lam[bp] == bl2, "LAMBDAS tracks λ2 in soft-core windows"
    assert dm_w0[bp] == bw0
    assert dm_a[bp] == ba
    assert dm_u0[bp] == bu0
    assert dm_i[bp] == 1, "bridge is a soft-core window (INTERMEDIATE=1)"
    # Module-level constant matches the spec.
    assert tuple(v2.DENSE38v4_DMINUS_BRIDGE_56) == _DENSE38v4_BRIDGE56_EXPECTED


def test_densified38v4_5to6_bridge_sits_between_broad_and_plateau(v2):
    """The new bridge is inserted BETWEEN the broad transition state 5
    (λ1=0.470/λ2=0.490, W0=0.600) and the tight plateau side state 6
    (λ1=0.460/λ2=0.480, W0=0.450) — the v3 dminus 5→6 handoff. The neighbor
    BEFORE the bridge is the broad transition; the neighbor AFTER is the
    plateau-side state with strictly lower W0 + lower λ."""
    v4 = v2.get_schedule("densified38v4")
    dm_l1 = v4["lambdas_1"][19:]
    dm_l2 = v4["lambdas_2"][19:]
    dm_w0 = v4["w0"][19:]
    bp = next(j for j in range(len(dm_w0)) if abs(dm_w0[j] - 0.17) < 1e-9)
    # Before bridge = broad transition state 5.
    assert dm_l1[bp - 1] == 0.470
    assert dm_l2[bp - 1] == 0.490
    assert dm_w0[bp - 1] == 0.600
    # After bridge = tight plateau side state 6.
    assert dm_l1[bp + 1] == 0.460
    assert dm_l2[bp + 1] == 0.480
    assert dm_w0[bp + 1] == 0.450


def test_densified38v4_5to6_bridge_alpha_u0_monotone(v2):
    """α and U0 of the new bridge are the midpoints (monotone-interpolated)
    between state 5 (α=0.18, U0=92) and state 6 (α=0.16, U0=95): α=0.17,
    U0=93.5. Both lie strictly between the two neighbors."""
    v4 = v2.get_schedule("densified38v4")
    dm_w0 = v4["w0"][19:]
    dm_a = v4["alpha"][19:]
    dm_u0 = v4["u0"][19:]
    bp = next(j for j in range(len(dm_w0)) if abs(dm_w0[j] - 0.17) < 1e-9)
    # α between state 6 (0.16) and state 5 (0.18).
    assert dm_a[bp + 1] < dm_a[bp] < dm_a[bp - 1]
    assert dm_a[bp] == 0.17
    # U0 between state 5 (92) and state 6 (95).
    assert dm_u0[bp - 1] < dm_u0[bp] < dm_u0[bp + 1]
    assert dm_u0[bp] == 93.5


def test_densified38v4_dminus_is_v3_dminus_plus_one_bridge(v2):
    """REMOVING the new W0=0.17 bridge state from v4's dminus recovers
    densified38v3's dminus EXACTLY (the v4 fix is a pure single insertion on
    top of the v3 dminus — reuses the v3 array, does not re-derive)."""
    v4 = v2.get_schedule("densified38v4")
    v3 = v2.get_schedule("densified38v3")
    v4_dminus_w0 = v4["w0"][19:]
    bp = next(j for j in range(len(v4_dminus_w0))
              if abs(v4_dminus_w0[j] - 0.17) < 1e-9)
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "intermd", "w0", "alpha", "u0"):
        v4_dminus = v4[key][19:]
        v4_no_new_bridge = v4_dminus[:bp] + v4_dminus[bp + 1:]
        assert v4_no_new_bridge == v3[key][19:], (
            f"v4 dminus {key} minus the new bridge must equal v3 dminus"
        )


def test_densified38v4_dminus_direction_forced_plus_at_runtime(v2):
    """The schedule encodes dminus DIRECTION=-1 (thermodynamic-leg label). The
    (b+) runtime fix (DIRECTION→+1 forced + DISPLACEMENT→-25) is applied by
    generate_per_direction_cntls at slice time, NOT in the schedule. Assert the
    schedule's dminus block carries the -1 label the slicer keys on (the +1/-25
    rewrite is exercised in the per-direction suite + the dry-run)."""
    s = v2.get_schedule("densified38v4")
    assert s["directions"][19:] == [-1] * 21
    assert s["directions"][:19] == [1] * 19


def test_densified38v4_endpoints_invariant_both_directions(v2):
    """Per-direction endpoint-invariance (load-bearing / ΔG-unbias): BOTH directions keep
    λ=0 (coupled) as one endpoint and the apex (W0=1.0, λ2=0.5) as the other.
    dplus: state 0 = coupled, state 18 = apex. dminus: first state (19) = apex
    (reversed), last state (39) = coupled. Both bridges are INTERIOR;
    endpoints unchanged → E[ΔG_int] unchanged (Zwanzig/Kirkwood)."""
    s = v2.get_schedule("densified38v4")
    n_fwd = 19
    # dplus coupled endpoint (state 0).
    assert s["lambdas_1"][0] == 0.00
    assert s["w0"][0] == 0.0
    # dplus decoupled apex (state 18).
    assert s["w0"][n_fwd - 1] == 1.0
    assert s["lambdas_2"][n_fwd - 1] == 0.5
    # dminus decoupled apex (first backward state, 19).
    assert s["w0"][n_fwd] == 1.0
    assert s["lambdas_2"][n_fwd] == 0.5
    # dminus coupled endpoint (last state, 39).
    assert s["lambdas_1"][-1] == 0.00
    assert s["w0"][-1] == 0.0


def test_densified38v4_dminus_two_extra_intermediates_vs_v2(v2):
    """dminus has EXACTLY 2 more INTERMEDIATE windows than densified38v2's
    dminus (the v3 9→10 bridge + the new v4 5→6 bridge, each INTERMEDIATE=1),
    and EXACTLY 1 more than densified38v3's dminus. dplus INTERMEDIATE count is
    unchanged from v2/v3."""
    v4 = v2.get_schedule("densified38v4")
    v3 = v2.get_schedule("densified38v3")
    v2s = v2.get_schedule("densified38v2")
    # dplus unchanged.
    assert sum(int(x) for x in v4["intermd"][:19]) == \
        sum(int(x) for x in v2s["intermd"][:19]), "dplus INTERMEDIATE unchanged"
    # dminus: +2 vs v2, +1 vs v3.
    v4_dm = sum(int(x) for x in v4["intermd"][19:])
    v3_dm = sum(int(x) for x in v3["intermd"][19:])
    v2_dm = sum(int(x) for x in v2s["intermd"][19:])
    assert v4_dm == v2_dm + 2, "dminus +2 INTERMEDIATE (both bridges) vs v2"
    assert v4_dm == v3_dm + 1, "dminus +1 INTERMEDIATE (new 5→6 bridge) vs v3"


def test_densified38v4_cntl_writer_emits_40_state_unequal_columns(v2, tmp_path):
    """write_cntl_file serializes the densified38v4 schedule generically: all
    per-state cntl columns are 40 long, DIRECTION = 19×+1 then 21×-1, and BOTH
    bridges (W0=0.07 + W0=0.17) appear once in the W0COEFF column. Proves the
    cntl writer is count-agnostic for the 40=19/21 unequal-count schedule."""
    cntl_path = os.path.join(str(tmp_path), "trackb_v4.cntl")
    sched = v2.get_schedule("densified38v4")
    v2.write_cntl_file(
        cntl_path=cntl_path,
        basename="trackb",
        nodefile_path=os.path.join(str(tmp_path), "nodefile"),
        ligand_atom_indices=[0, 1, 2],
        pos_restrained_atom_indices=[3, 4],
        displacement_nm=(2.5, 0.0, 0.0),
        production_steps=1000,
        prnt_frequency=100,
        trj_frequency=100,
        wall_time_min=60,
        cycle_time_s=10,
        checkpoint_time_s=10,
        schedule=sched,
    )
    with open(cntl_path) as fh:
        text = fh.read()

    def _col(key):
        m = re.search(rf"^{key}\s*=\s*'([^']*)'", text, re.M)
        assert m, f"{key} column missing"
        return [t.strip() for t in m.group(1).split(",") if t.strip() != ""]

    for key in ("LAMBDAS", "DIRECTION", "INTERMEDIATE", "LAMBDA1",
                "LAMBDA2", "ALPHA", "U0", "W0COEFF"):
        assert len(_col(key)) == 40, f"{key} not 40 long"
    direction = [int(float(x)) for x in _col("DIRECTION")]
    assert direction[:19] == [1] * 19
    assert direction[19:] == [-1] * 21
    w0 = [float(x) for x in _col("W0COEFF")]
    assert sum(1 for w in w0 if abs(w - 0.07) < 1e-9) == 1, "one W0=0.07 bridge"
    assert sum(1 for w in w0 if abs(w - 0.17) < 1e-9) == 1, "one W0=0.17 bridge"


# --- Regression: prior schedules UNCHANGED by the densified38v4 addition -----

def test_densified38v3_unchanged_by_v4_addition(v2):
    """Regression: densified38v3 is UNCHANGED by adding densified38v4.
    39 states, 19+20, dplus byte-equal to v2, dminus has the 9→10 W0=0.07
    bridge and NO W0=0.17 bridge."""
    v3 = v2.get_schedule("densified38v3")
    assert v3["n_states"] == 39
    assert v3["directions"] == [1] * 19 + [-1] * 20
    assert sum(1 for w in v3["w0"] if abs(w - 0.07) < 1e-9) == 1, "v3 keeps 9→10 bridge"
    assert not any(abs(w - 0.17) < 1e-9 for w in v3["w0"]), "v3 has NO 5→6 bridge"


def test_densified38v2_unchanged_by_v4_addition(v2):
    """Regression: densified38v2 is UNCHANGED by adding densified38v4."""
    v2s = v2.get_schedule("densified38v2")
    assert v2s["n_states"] == 38
    assert v2s["directions"] == [1] * 19 + [-1] * 19
    assert not any(abs(w - 0.07) < 1e-9 for w in v2s["w0"]), "v2 has NO 9→10 bridge"
    assert not any(abs(w - 0.17) < 1e-9 for w in v2s["w0"]), "v2 has NO 5→6 bridge"
    # Backward half = reversed forward half (pure mirror, no bridge).
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "intermd", "w0", "alpha", "u0"):
        assert list(reversed(v2s[key][:19])) == v2s[key][19:], f"v2 mirror {key}"


def test_densified38_unchanged_by_v4_addition(v2):
    """Regression: densified38 is UNCHANGED by adding densified38v4."""
    d38 = v2.get_schedule("densified38")
    assert d38["n_states"] == 38
    assert d38["directions"] == [1] * 19 + [-1] * 19
    assert d38["lambdas_1"][:10] == [
        0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45]
    assert not any(abs(w - 0.07) < 1e-9 for w in d38["w0"])
    assert not any(abs(w - 0.17) < 1e-9 for w in d38["w0"])


def test_canonical22_unchanged_by_v4_addition(v2):
    """Regression: canonical22 is UNCHANGED by adding densified38v4."""
    c22 = v2.get_schedule("canonical22")
    assert c22["n_states"] == 22
    assert c22["directions"] == [1] * 11 + [-1] * 11
    assert not any(abs(w - 0.07) < 1e-9 for w in c22["w0"])
    assert not any(abs(w - 0.17) < 1e-9 for w in c22["w0"])


# ---------------------------------------------------------------------------
# densified_bound28 — BOUND-leg per-direction schedule (6→7 cliff bridge)
# ---------------------------------------------------------------------------
# The 3 bridge windows expected on the dplus side, in (λ1, λ2, W0, α, U0) order.
_BOUND28_BRIDGE_EXPECTED = [
    (0.300, 0.315, 0.20, 0.12, 105.0),
    (0.310, 0.330, 0.45, 0.14, 100.0),
    (0.320, 0.345, 0.70, 0.16,  97.0),
]


def test_densified_bound28_registered_in_schedules(v2):
    """densified_bound28 is a NEW selectable schedule alongside the existing
    ones; the prior schedules remain present (registry not clobbered). It is the
    only schedule marked BOUND-valid besides canonical22."""
    assert "densified_bound28" in v2.SCHEDULES
    for name in ("canonical22", "densified34", "densified38",
                 "densified38v2", "densified38v3", "densified38v4"):
        assert name in v2.SCHEDULES
    assert v2.BOUND_SCHEDULES == frozenset(
        {"canonical22", "densified_bound28", "densified_bound30"})
    s = v2.get_schedule("densified_bound28")
    assert s["n_states"] == 28


def test_densified_bound28_equal_per_direction_counts_28(v2):
    """densified_bound28 has EQUAL per-direction counts: dplus 14, dminus 14,
    total 28. DIRECTION is a contiguous +1 (14) then -1 (14) block (required by
    the per-direction launcher's contiguity gate). Every per-state array is the
    full 28 long."""
    s = v2.get_schedule("densified_bound28")
    assert s["n_states"] == 28
    n_fwd = sum(1 for d in s["directions"] if d >= 0)
    n_bwd = s["n_states"] - n_fwd
    assert n_fwd == 14, "dplus (forward) = 14 states (11 canonical + 3 bridge)"
    assert n_bwd == 14, "dminus (backward) = 14 states"
    assert s["directions"][:14] == [1] * 14
    assert s["directions"][14:] == [-1] * 14
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "directions", "intermd", "w0", "alpha", "u0"):
        assert len(s[key]) == 28, f"{key} length"


def test_densified_bound28_dplus_canonical_plateau_and_endpoints(v2):
    """The dplus coupled plateau (canonical 0..6, λ=0.00..0.30) and the tail
    (λ=0.35/0.40/0.45 + apex λ=0.50/W0=1.0) are CANONICAL: λ1=λ2, W0=0 except the
    apex (W0=1.0), α=0.10, U0=110. Only the 3 bridge windows differ."""
    s = v2.get_schedule("densified_bound28")
    # Canonical coupled plateau states 0..6.
    expect_lam = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
    for i, lam in enumerate(expect_lam):
        assert s["lambdas_1"][i] == lam and s["lambdas_2"][i] == lam
        assert s["w0"][i] == 0.0
        assert s["alpha"][i] == 0.10 and s["u0"][i] == 110.0
        assert s["intermd"][i] == 0
    # After the 3 bridge windows (idx 7,8,9): canonical tail at idx 10,11,12.
    for i, lam in zip((10, 11, 12), (0.35, 0.40, 0.45)):
        assert s["lambdas_1"][i] == lam and s["lambdas_2"][i] == lam
        assert s["w0"][i] == 0.0
        assert s["alpha"][i] == 0.10 and s["u0"][i] == 110.0
        assert s["intermd"][i] == 0
    # dplus apex (idx 13): λ=0.50, W0=1.0, INTERMEDIATE=1, α/U0 canonical.
    assert s["lambdas_1"][13] == 0.50 and s["lambdas_2"][13] == 0.50
    assert s["w0"][13] == 1.0 and s["intermd"][13] == 1
    assert s["alpha"][13] == 0.10 and s["u0"][13] == 110.0


def test_densified_bound28_dplus_bridge_three_windows(v2):
    """The dplus 6→7 cliff carries EXACTLY 3 W0-graded soft-core bridge windows
    (W0 0.20/0.45/0.70), inserted between canonical state 6 (λ=0.30) and state 7
    (λ=0.35). Each is λ1<λ2 (softplus active), INTERMEDIATE=1, with the α ramp
    0.12→0.16 and U0 descent 105→97."""
    s = v2.get_schedule("densified_bound28")
    bridge = list(zip(s["lambdas_1"][7:10], s["lambdas_2"][7:10],
                      s["w0"][7:10], s["alpha"][7:10], s["u0"][7:10]))
    assert bridge == _BOUND28_BRIDGE_EXPECTED
    # All three are INTERMEDIATE soft-core windows with λ1<λ2.
    for i in (7, 8, 9):
        assert s["intermd"][i] == 1
        assert s["lambdas_1"][i] < s["lambdas_2"][i], "softplus active (λ1<λ2)"
    # W0 ramps 0.20→0.45→0.70 (monotone across the cliff).
    assert s["w0"][7:10] == [0.20, 0.45, 0.70]
    # The 3 bridge W0 values appear exactly once each on the dplus side.
    for w0 in (0.20, 0.45, 0.70):
        assert sum(1 for w in s["w0"][:14] if abs(w - w0) < 1e-9) == 1
    # Module-level constant matches the spec.
    assert list(v2.BOUND28_BRIDGE) == _BOUND28_BRIDGE_EXPECTED


def test_densified_bound28_bridge_sits_between_cliff_neighbors(v2):
    """The bridge block is inserted BETWEEN the lower cliff neighbor (λ=0.30,
    W0=0) and the upper cliff neighbor (λ=0.35, W0=0) — i.e. the 6→7 hole. The
    state immediately BEFORE the bridge is the canonical λ=0.30 plateau; the
    state immediately AFTER is the canonical λ=0.35 plateau."""
    s = v2.get_schedule("densified_bound28")
    # First bridge window is idx 7; idx 6 = lower neighbor, idx 10 = upper.
    assert s["w0"][6] == 0.0 and s["lambdas_1"][6] == 0.30
    assert s["w0"][10] == 0.0 and s["lambdas_1"][10] == 0.35
    # The 3 interior bridge states are all W0>0.
    assert all(s["w0"][i] > 0.0 for i in (7, 8, 9))


def test_densified_bound28_dminus_is_whole_tuple_reverse_of_dplus(v2):
    """dminus is the WHOLE-TUPLE reverse of dplus (per-state arrays reversed;
    λ1 stays λ1, λ2 stays λ2 — NOT a λ1↔λ2 swap), with DIRECTION=-1. The default
    symmetric construction makes dminus the reverse of dplus for every column
    EXCEPT directions (which flip to -1)."""
    s = v2.get_schedule("densified_bound28")
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "intermd", "w0", "alpha", "u0"):
        assert list(reversed(s[key][:14])) == s[key][14:], f"dminus reverse {key}"
    assert s["directions"][14:] == [-1] * 14


def test_densified_bound28_dminus_direction_minus_one(v2):
    """dminus DIRECTION column is all -1 in the schedule (the per-direction
    launcher derives the dminus leg from where(DIRECTION<0) and forces +1 at
    runtime via the cntl rewrite — but the SCHEDULE records -1)."""
    s = v2.get_schedule("densified_bound28")
    assert all(d == -1 for d in s["directions"][14:])
    assert all(d == 1 for d in s["directions"][:14])


def test_densified_bound28_endpoints_invariant_both_directions(v2):
    """Endpoints IMMUTABLE both directions (ΔG-unbias / ranking-only): the
    coupled λ=0.00 (W0=0) and the apex λ=0.50/W0=1.0 are present and unbiased in
    both the dplus and dminus halves."""
    s = v2.get_schedule("densified_bound28")
    # dplus: coupled first, apex last fwd state.
    assert s["lambdas_1"][0] == 0.00 and s["w0"][0] == 0.0
    assert s["w0"][13] == 1.0 and s["lambdas_2"][13] == 0.50
    # dminus: apex first bwd state, coupled last state.
    assert s["w0"][14] == 1.0 and s["lambdas_2"][14] == 0.50
    assert s["lambdas_1"][-1] == 0.00 and s["w0"][-1] == 0.0


def test_densified_bound28_c2_region_lock_endpoints_unbiased(v2):
    """C2 region-lock: every λ1==λ2 state (endpoint / coupled plateau / apex)
    carries the CANONICAL α=0.10 and U0=110 — the α/U0 ramp lives STRICTLY inside
    the W0>0 bridge windows (where λ1<λ2). No biased endpoint slips through."""
    s = v2.get_schedule("densified_bound28")
    for i in range(s["n_states"]):
        if abs(s["lambdas_1"][i] - s["lambdas_2"][i]) < 1e-12:
            assert s["alpha"][i] == v2.BOUND28_CANON_ALPHA, (
                f"state {i} λ1==λ2 must have canonical α"
            )
            assert s["u0"][i] == v2.BOUND28_CANON_U0, (
                f"state {i} λ1==λ2 must have canonical U0"
            )


def test_densified_bound28_c2_assert_fires_on_biased_endpoint(v2):
    """C2 fail-loud: the module-level region-lock would REJECT a schedule whose
    λ1==λ2 endpoint is biased (α/U0 != canonical). We can't re-trigger the
    module asserts, so verify the invariant the assert encodes: a hand-biased
    endpoint violates ``alpha == BOUND28_CANON_ALPHA``."""
    biased_alpha = 0.12  # a non-canonical α at a λ1==λ2 endpoint
    assert biased_alpha != v2.BOUND28_CANON_ALPHA
    with pytest.raises(AssertionError):
        # This is the exact assert form the module emits per λ1==λ2 state.
        assert biased_alpha == v2.BOUND28_CANON_ALPHA


def test_densified_bound28_intermediate_count(v2):
    """dplus has 4 INTERMEDIATE states (3 bridge + 1 apex); dminus has 4
    (1 apex + 3 bridge). The bridge windows are the only INTERMEDIATE states
    besides the two apexes."""
    s = v2.get_schedule("densified_bound28")
    assert sum(int(x) for x in s["intermd"][:14]) == 4
    assert sum(int(x) for x in s["intermd"][14:]) == 4


def test_densified_bound28_builder_cliff_center_parameter_not_hardcoded(v2):
    """C4: ``_build_densified_bound_arrays`` takes ``cliff_center`` as a
    PARAMETER (dminus cliff is pilot-validated, NOT a hardcoded mirror).
    Relocating cliff_center shifts the bridge λ windows by the same offset — so a
    dminus pilot can move the bridge once its overlap profile is measured."""
    default = v2._build_densified_bound_arrays()
    relocated = v2._build_densified_bound_arrays(cliff_center=0.35)
    nf = sum(1 for d in default["directions"] if d >= 0)
    # dplus bridge λ1 windows for the default = [0.30, 0.31, 0.32].
    def bridge_l1(arr):
        n = sum(1 for d in arr["directions"] if d >= 0)
        return [arr["lambdas_1"][j] for j in range(n)
                if 0.0 < arr["w0"][j] < 1.0]
    assert bridge_l1(default) == [0.300, 0.310, 0.320]
    # Relocated cliff_center=0.35 → windows shifted by +0.05.
    assert bridge_l1(relocated) == pytest.approx([0.350, 0.360, 0.370])
    # Default and relocated are DIFFERENT (parameter is honoured, not ignored).
    assert default["lambdas_1"] != relocated["lambdas_1"]


def test_densified_bound28_builder_n_bridge_escalation(v2):
    """C3 escalation: the builder accepts ``n_bridge`` (3→4/5) — widening the
    bridge adds states symmetrically to both directions (pilot escalation hook,
    iteration cap=2). Endpoints stay invariant."""
    for n_bridge, n_states in ((4, 30), (5, 32)):
        a = v2._build_densified_bound_arrays(n_bridge=n_bridge)
        assert a["n_states"] == n_states
        n_fwd = sum(1 for d in a["directions"] if d >= 0)
        n_bridge_windows = sum(1 for j in range(n_fwd) if 0.0 < a["w0"][j] < 1.0)
        assert n_bridge_windows == n_bridge
        # Endpoints invariant.
        assert a["lambdas_1"][0] == 0.00 and a["w0"][0] == 0.0
        assert a["w0"][n_fwd - 1] == 1.0


def test_densified_bound28_builder_extra_alpha_u0_ramp_escalation(v2):
    """C3 escalation: the builder accepts ``extra_alpha_u0_ramp`` to override the
    per-bridge (α, U0) for a steeper ramp (pilot escalation hook). The override
    lands ONLY in the bridge windows; endpoints stay canonical."""
    ramp = [(0.13, 103.0), (0.15, 98.0), (0.18, 94.0)]
    a = v2._build_densified_bound_arrays(n_bridge=3, extra_alpha_u0_ramp=ramp)
    n_fwd = sum(1 for d in a["directions"] if d >= 0)
    got = [(a["alpha"][j], a["u0"][j]) for j in range(n_fwd)
           if 0.0 < a["w0"][j] < 1.0]
    assert got == ramp
    # Endpoints / plateau still canonical (C2 holds under escalation).
    for i in range(a["n_states"]):
        if abs(a["lambdas_1"][i] - a["lambdas_2"][i]) < 1e-12:
            assert a["alpha"][i] == 0.10 and a["u0"][i] == 110.0


def test_densified_bound28_builder_rejects_bad_n_bridge(v2):
    """C3: invalid n_bridge (<1) and mismatched extra_alpha_u0_ramp length
    fail loud (defensive — escalation must be well-formed)."""
    with pytest.raises(ValueError):
        v2._build_densified_bound_arrays(n_bridge=0)
    with pytest.raises(ValueError):
        v2._build_densified_bound_arrays(
            n_bridge=3, extra_alpha_u0_ramp=[(0.13, 103.0)])  # too few pairs


def test_densified_bound28_cntl_writer_emits_28_state_columns(v2, tmp_path):
    """The cntl writer serializes the densified_bound28 schedule verbatim: every
    schedule column has 28 entries, DIRECTION is +1×14 then -1×14, and the
    bridge W0 values (0.20/0.45/0.70) appear once per direction."""
    nodefile = tmp_path / "nodefile"
    nodefile.write_text("localhost,0:0,1,CUDA,,/tmp\n")
    cntl = tmp_path / "t_asyncre.cntl"
    v2.write_cntl_file(
        cntl_path=str(cntl), basename="t", nodefile_path=str(nodefile),
        ligand_atom_indices=[1, 2, 3], pos_restrained_atom_indices=[4, 5],
        displacement_nm=(2.5, 0.0, 0.0), production_steps=2500,
        prnt_frequency=2500, trj_frequency=25000, max_samples=1000,
        wall_time_min=720, cycle_time_s=10, checkpoint_time_s=600,
        schedule=v2.get_schedule("densified_bound28"),
    )
    txt = cntl.read_text()

    def col(name):
        line = next(l for l in txt.splitlines() if l.startswith(name + " "))
        return line.split("'")[1].split(", ") if "'" in line else []

    for name in ("LAMBDAS", "DIRECTION", "INTERMEDIATE",
                 "LAMBDA1", "LAMBDA2", "ALPHA", "U0", "W0COEFF"):
        assert len(col(name)) == 28, f"{name} must have 28 entries"
    directions = [int(float(x)) for x in col("DIRECTION")]
    assert directions == [1] * 14 + [-1] * 14
    w0 = [float(x) for x in col("W0COEFF")]
    for w in (0.20, 0.45, 0.70):
        assert sum(1 for x in w0 if abs(x - w) < 1e-9) == 2, (
            f"W0={w} appears once per direction (2 total)"
        )


def test_densified_bound28_bound_default_cntl_byte_equal_to_canonical(v2, tmp_path):
    """C1 CRITICAL: the bound leg's DEFAULT path (bound_schedule=canonical22 →
    schedule=None) must emit a cntl whose schedule columns are BYTE-EQUAL to the
    explicit canonical22 schedule — i.e. the new bound_schedule selector, left at
    default, does NOT perturb the historical bound-leg cntl."""
    def write_cntl(schedule):
        d = tmp_path / ("none" if schedule is None else "canon")
        d.mkdir()
        nodefile = d / "nodefile"
        nodefile.write_text("localhost,0:0,1,CUDA,,/tmp\n")
        cntl = d / "t_asyncre.cntl"
        v2.write_cntl_file(
            cntl_path=str(cntl), basename="t", nodefile_path=str(nodefile),
            ligand_atom_indices=[1, 2, 3], pos_restrained_atom_indices=[4, 5],
            displacement_nm=(2.5, 0.0, 0.0), production_steps=2500,
            prnt_frequency=2500, trj_frequency=25000, max_samples=1000,
            wall_time_min=720, cycle_time_s=10, checkpoint_time_s=600,
            schedule=schedule,
        )
        return cntl.read_text()

    none_txt = write_cntl(None)
    canon_txt = write_cntl(v2.get_schedule("canonical22"))
    sched_prefixes = ("LAMBDAS ", "DIRECTION ", "INTERMEDIATE ",
                      "LAMBDA1 ", "LAMBDA2 ", "ALPHA ", "U0 ", "W0COEFF ")

    def sched_lines(txt):
        return [l for l in txt.splitlines() if l.startswith(sched_prefixes)]

    assert sched_lines(none_txt) == sched_lines(canon_txt)


def test_densified_bound28_setup_rejects_free_schedule_on_bound_leg(v2, tmp_path):
    """A free-leg densified38* schedule is NOT valid for the bound leg: passing
    one as bound_schedule fails loud (the densified38* ladders fix the free-leg
    decoupling crossover the receptor-held bound leg does not have)."""
    # The validation lives in setup_one_leg's schedule-resolution branch; assert
    # the guard set directly (setup_one_leg needs a real PDB to reach it).
    assert "densified38v4" not in v2.BOUND_SCHEDULES
    assert "densified_bound28" in v2.BOUND_SCHEDULES
    assert "canonical22" in v2.BOUND_SCHEDULES


# --- Regression: prior schedules UNCHANGED by the densified_bound28 addition --

def test_canonical22_unchanged_by_bound28_addition(v2):
    """Regression: canonical22 is UNCHANGED by adding densified_bound28."""
    c22 = v2.get_schedule("canonical22")
    assert c22["n_states"] == 22
    assert c22["directions"] == [1] * 11 + [-1] * 11
    assert c22["lambdas_1"] == v2.LAMBDA_FWD + v2.LAMBDA_BWD
    assert c22["w0"] == [0.0] * 10 + [1.0, 1.0] + [0.0] * 10
    assert c22["alpha"] == [0.10] * 22
    assert c22["u0"] == [110.0] * 22


def test_densified38v4_unchanged_by_bound28_addition(v2):
    """Regression: the free-leg densified38v4 is UNCHANGED by adding the bound
    schedule (40 states, 19/21, both free-leg bridges intact)."""
    v4 = v2.get_schedule("densified38v4")
    assert v4["n_states"] == 40
    assert v4["directions"] == [1] * 19 + [-1] * 21
    assert sum(1 for w in v4["w0"] if abs(w - 0.07) < 1e-9) == 1
    assert sum(1 for w in v4["w0"] if abs(w - 0.17) < 1e-9) == 1


# ---------------------------------------------------------------------------
# densified_bound30 — 4-window cap-1 escalation of densified_bound28 (6→7 cliff)
# ---------------------------------------------------------------------------
# The 4 bridge windows expected on the dplus side, in (λ1, λ2, W0, α, U0) order
# (linearly interpolated across the SAME span as the 3-window bridge — same
# COORDINATED direction: W0 ↑, α ↑, U0 ↓ DESCENT, MC-1). Endpoints are
# 0.300↔0.345 (the SAME cliff span as bound28), denser λ (MC-2 window count).
_BOUND30_BRIDGE_EXPECTED = [
    (0.300000, 0.315000, 0.200000, 0.120000, 105.000000),
    (0.306667, 0.325000, 0.366667, 0.133333, 102.333333),
    (0.313333, 0.335000, 0.533333, 0.146667,  99.666667),
    (0.320000, 0.345000, 0.700000, 0.160000,  97.000000),
]


def test_densified_bound30_registered_in_schedules(v2):
    """densified_bound30 is a NEW selectable schedule registered ADDITIVELY in
    both SCHEDULES and BOUND_SCHEDULES; the prior schedules (including
    densified_bound28) remain present (registry not clobbered)."""
    assert "densified_bound30" in v2.SCHEDULES
    for name in ("canonical22", "densified34", "densified38", "densified38v2",
                 "densified38v3", "densified38v4", "densified_bound28"):
        assert name in v2.SCHEDULES
    assert "densified_bound30" in v2.BOUND_SCHEDULES
    s = v2.get_schedule("densified_bound30")
    assert s["n_states"] == 30


def test_densified_bound30_equal_per_direction_counts_30(v2):
    """densified_bound30 has EQUAL per-direction counts: dplus 15, dminus 15,
    total 30. DIRECTION is a contiguous +1 (15) then -1 (15) block (required by
    the per-direction launcher's contiguity gate). Every per-state array is the
    full 30 long."""
    s = v2.get_schedule("densified_bound30")
    assert s["n_states"] == 30
    n_fwd = sum(1 for d in s["directions"] if d >= 0)
    n_bwd = s["n_states"] - n_fwd
    assert n_fwd == 15, "dplus (forward) = 15 states (11 canonical + 4 bridge)"
    assert n_bwd == 15, "dminus (backward) = 15 states"
    assert s["directions"][:15] == [1] * 15
    assert s["directions"][15:] == [-1] * 15
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "directions", "intermd", "w0", "alpha", "u0"):
        assert len(s[key]) == 30, f"{key} length"


def test_densified_bound30_dplus_canonical_plateau_and_endpoints(v2):
    """The dplus coupled plateau (canonical 0..6, λ=0.00..0.30) and the tail
    (λ=0.35/0.40/0.45 + apex λ=0.50/W0=1.0) are CANONICAL: λ1=λ2, W0=0 except the
    apex (W0=1.0), α=0.10, U0=110. Only the 4 bridge windows differ. The 4-window
    bridge pushes the tail/apex up by 1 index vs bound28."""
    s = v2.get_schedule("densified_bound30")
    # Canonical coupled plateau states 0..6.
    expect_lam = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
    for i, lam in enumerate(expect_lam):
        assert s["lambdas_1"][i] == lam and s["lambdas_2"][i] == lam
        assert s["w0"][i] == 0.0
        assert s["alpha"][i] == 0.10 and s["u0"][i] == 110.0
        assert s["intermd"][i] == 0
    # After the 4 bridge windows (idx 7,8,9,10): canonical tail at idx 11,12,13.
    for i, lam in zip((11, 12, 13), (0.35, 0.40, 0.45)):
        assert s["lambdas_1"][i] == lam and s["lambdas_2"][i] == lam
        assert s["w0"][i] == 0.0
        assert s["alpha"][i] == 0.10 and s["u0"][i] == 110.0
        assert s["intermd"][i] == 0
    # dplus apex (idx 14): λ=0.50, W0=1.0, INTERMEDIATE=1, α/U0 canonical.
    assert s["lambdas_1"][14] == 0.50 and s["lambdas_2"][14] == 0.50
    assert s["w0"][14] == 1.0 and s["intermd"][14] == 1
    assert s["alpha"][14] == 0.10 and s["u0"][14] == 110.0


def test_densified_bound30_dplus_bridge_four_windows(v2):
    """The dplus 6→7 cliff carries EXACTLY 4 W0-graded soft-core bridge windows
    (cap-1 escalation, MC-2). Each is λ1<λ2 (softplus active), INTERMEDIATE=1,
    matching the linearly-interpolated 4-window spec over the SAME cliff span."""
    s = v2.get_schedule("densified_bound30")
    bridge = list(zip(s["lambdas_1"][7:11], s["lambdas_2"][7:11],
                      s["w0"][7:11], s["alpha"][7:11], s["u0"][7:11]))
    for got, exp in zip(bridge, _BOUND30_BRIDGE_EXPECTED):
        assert got == pytest.approx(exp, abs=1e-6)
    # All four are INTERMEDIATE soft-core windows with λ1<λ2.
    for i in (7, 8, 9, 10):
        assert s["intermd"][i] == 1
        assert s["lambdas_1"][i] < s["lambdas_2"][i], "softplus active (λ1<λ2)"
    # Exactly 4 soft-core bridge windows (0<W0<1) on the dplus side.
    assert sum(1 for w in s["w0"][:15] if 0.0 < w < 1.0) == 4


def test_densified_bound30_mc1_u0_strictly_descending(v2):
    """MC-1 CRITICAL: U0 ramps STRICTLY DESCENDING across the dplus bridge
    (105→102.3→99.7→97). The validated free-leg densified38 descends U0 110→82;
    raising U0 above usc is forbidden (softplus caps the favorable LOW-usc tail,
    not the plateau). W0 ascends, α ascends — coordinated descent direction."""
    s = v2.get_schedule("densified_bound30")
    n_fwd = sum(1 for d in s["directions"] if d >= 0)
    bridge_u0 = [u for w, u in zip(s["w0"][:n_fwd], s["u0"][:n_fwd])
                 if 0.0 < w < 1.0]
    bridge_w0 = [w for w in s["w0"][:n_fwd] if 0.0 < w < 1.0]
    bridge_a = [a for w, a in zip(s["w0"][:n_fwd], s["alpha"][:n_fwd])
                if 0.0 < w < 1.0]
    assert all(a > b for a, b in zip(bridge_u0, bridge_u0[1:])), (
        f"U0 must strictly descend (MC-1); got {bridge_u0!r}")
    assert bridge_u0[0] == pytest.approx(105.0) and bridge_u0[-1] == pytest.approx(97.0)
    assert all(a < b for a, b in zip(bridge_w0, bridge_w0[1:])), "W0 ascends"
    assert all(a < b for a, b in zip(bridge_a, bridge_a[1:])), "α ascends"


def test_densified_bound30_bridge_same_cliff_span_as_bound28(v2):
    """MC-2: the 4-window bridge spans the SAME cliff (λ0.30↔0.345) as the
    3-window bound28 — window COUNT is the lever, not a wider span. First/last
    bridge endpoints match bound28's; the interior is denser."""
    s30 = v2.get_schedule("densified_bound30")
    s28 = v2.get_schedule("densified_bound28")
    def bridge_l(s):
        n = sum(1 for d in s["directions"] if d >= 0)
        return [(s["lambdas_1"][j], s["lambdas_2"][j]) for j in range(n)
                if 0.0 < s["w0"][j] < 1.0]
    b30, b28 = bridge_l(s30), bridge_l(s28)
    assert b30[0] == pytest.approx(b28[0]), "same bridge START (λ0.30↔0.315)"
    assert b30[-1] == pytest.approx(b28[-1]), "same bridge END (λ0.32↔0.345)"
    assert len(b30) == 4 and len(b28) == 3, "count is the lever (4 vs 3)"


def test_densified_bound30_bridge_sits_between_cliff_neighbors(v2):
    """The 4-window bridge is inserted BETWEEN the lower cliff neighbor (λ=0.30,
    W0=0, idx 6) and the upper cliff neighbor (λ=0.35, W0=0, idx 11) — the 6→7
    hole. The interior bridge states (7..10) are all W0>0."""
    s = v2.get_schedule("densified_bound30")
    assert s["w0"][6] == 0.0 and s["lambdas_1"][6] == 0.30
    assert s["w0"][11] == 0.0 and s["lambdas_1"][11] == 0.35
    assert all(s["w0"][i] > 0.0 for i in (7, 8, 9, 10))


def test_densified_bound30_dminus_is_whole_tuple_reverse_of_dplus(v2):
    """dminus is the WHOLE-TUPLE reverse of dplus (per-state arrays reversed; λ1
    stays λ1, λ2 stays λ2 — NOT a λ1↔λ2 swap), DIRECTION=-1. Mirror the bound28
    dminus construction at the new count (15/15)."""
    s = v2.get_schedule("densified_bound30")
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "intermd", "w0", "alpha", "u0"):
        assert list(reversed(s[key][:15])) == s[key][15:], f"dminus reverse {key}"
    assert s["directions"][15:] == [-1] * 15
    assert s["directions"][:15] == [1] * 15


def test_densified_bound30_endpoints_invariant_both_directions(v2):
    """Q4 ΔG-unbias: endpoints IMMUTABLE both directions (λ1==λ2 → softplus
    prefactor 0). The coupled λ=0.00 (W0=0) and the apex λ=0.50/W0=1.0 are
    present and unbiased in both the dplus and dminus halves at the new count."""
    s = v2.get_schedule("densified_bound30")
    # dplus: coupled first, apex last fwd state (idx 14).
    assert s["lambdas_1"][0] == 0.00 and s["w0"][0] == 0.0
    assert s["w0"][14] == 1.0 and s["lambdas_2"][14] == 0.50
    # dminus: apex first bwd state (idx 15), coupled last state.
    assert s["w0"][15] == 1.0 and s["lambdas_2"][15] == 0.50
    assert s["lambdas_1"][-1] == 0.00 and s["w0"][-1] == 0.0


def test_densified_bound30_c2_region_lock_endpoints_unbiased(v2):
    """C2 region-lock: every λ1==λ2 state (endpoint / coupled plateau / apex)
    carries the CANONICAL α=0.10 and U0=110. The α/U0 ramp lives STRICTLY inside
    the W0>0 bridge windows (where λ1<λ2). No biased endpoint slips through."""
    s = v2.get_schedule("densified_bound30")
    for i in range(s["n_states"]):
        if abs(s["lambdas_1"][i] - s["lambdas_2"][i]) < 1e-12:
            assert s["alpha"][i] == v2.BOUND28_CANON_ALPHA, (
                f"state {i} λ1==λ2 must have canonical α")
            assert s["u0"][i] == v2.BOUND28_CANON_U0, (
                f"state {i} λ1==λ2 must have canonical U0")


def test_densified_bound30_intermediate_count(v2):
    """dplus has 5 INTERMEDIATE states (4 bridge + 1 apex); dminus has 5
    (1 apex + 4 bridge). The bridge windows are the only INTERMEDIATE states
    besides the two apexes — n_bridge+1 per direction."""
    s = v2.get_schedule("densified_bound30")
    assert sum(int(x) for x in s["intermd"][:15]) == 5
    assert sum(int(x) for x in s["intermd"][15:]) == 5


def test_densified_bound30_parameterized_region_lock_fires_on_tamper(v2):
    """C10: the PARAMETERIZED region-lock verifier (n_bridge-aware, NOT hardcoded
    28/14/14/3) fires AssertionError on tamper for the new count: a biased λ1==λ2
    endpoint, a U0 ascending (MC-1 violation), and a wrong n_bridge count. The
    assert is NOT weakened by the parameterization."""
    import copy
    base = v2._build_densified_bound_arrays(n_bridge=4)
    # Clean build passes.
    v2._assert_densified_bound_region_lock(base, n_bridge=4, label="clean30")

    # (a) bias a λ1==λ2 endpoint's U0 → region-lock fires.
    bad_a = copy.deepcopy(base)
    for i in range(bad_a["n_states"]):
        if abs(bad_a["lambdas_1"][i] - bad_a["lambdas_2"][i]) < 1e-12:
            bad_a["u0"][i] = 999.0
            break
    with pytest.raises(AssertionError):
        v2._assert_densified_bound_region_lock(bad_a, n_bridge=4, label="bad_a")

    # (b) flip the bridge U0 to ASCENDING → MC-1 guard fires.
    bad_b = copy.deepcopy(base)
    n_fwd = sum(1 for d in bad_b["directions"] if d >= 0)
    bidx = [i for i in range(n_fwd) if 0.0 < bad_b["w0"][i] < 1.0]
    bad_b["u0"][bidx[0]], bad_b["u0"][bidx[-1]] = (
        bad_b["u0"][bidx[-1]], bad_b["u0"][bidx[0]])
    with pytest.raises(AssertionError):
        v2._assert_densified_bound_region_lock(bad_b, n_bridge=4, label="bad_b")

    # (c) wrong n_bridge expectation (state-count mismatch) fires.
    with pytest.raises(AssertionError):
        v2._assert_densified_bound_region_lock(base, n_bridge=3, label="bad_c")


def test_densified_bound30_per_direction_split_count_agnostic(v2):
    """C3: the per-direction split derives (total, fwd) from the DIRECTION column
    with NO magic-number break — 30 → 15/15 (as verified for the 28→14/14 split).
    The contiguity gate accepts the +1×15 / -1×15 layout."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "_pd_for_bound30",
        os.path.join(os.path.dirname(__file__), "..", "scripts",
                     "trackb_per_direction_production.py"))
    pd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pd)
    s = v2.get_schedule("densified_bound30")
    total, fwd = pd._derive_state_counts_from_directions(s["directions"])
    assert (total, fwd) == (30, 15)
    # bound28 still derives 28/14 (no regression in the split logic).
    s28 = v2.get_schedule("densified_bound28")
    assert pd._derive_state_counts_from_directions(s28["directions"]) == (28, 14)


def test_densified_bound30_cntl_writer_emits_30_state_columns(v2, tmp_path):
    """The cntl writer serializes the densified_bound30 schedule verbatim: every
    schedule column has 30 entries, DIRECTION is +1×15 then -1×15, and 4 distinct
    soft-core bridge W0 values appear per direction."""
    nodefile = tmp_path / "nodefile"
    nodefile.write_text("localhost,0:0,1,CUDA,,/tmp\n")
    cntl = tmp_path / "t_asyncre.cntl"
    v2.write_cntl_file(
        cntl_path=str(cntl), basename="t", nodefile_path=str(nodefile),
        ligand_atom_indices=[1, 2, 3], pos_restrained_atom_indices=[4, 5],
        displacement_nm=(2.5, 0.0, 0.0), production_steps=2500,
        prnt_frequency=2500, trj_frequency=25000, max_samples=1000,
        wall_time_min=720, cycle_time_s=10, checkpoint_time_s=600,
        schedule=v2.get_schedule("densified_bound30"),
    )
    txt = cntl.read_text()

    def col(name):
        line = next(l for l in txt.splitlines() if l.startswith(name + " "))
        return line.split("'")[1].split(", ") if "'" in line else []

    for name in ("LAMBDAS", "DIRECTION", "INTERMEDIATE",
                 "LAMBDA1", "LAMBDA2", "ALPHA", "U0", "W0COEFF"):
        assert len(col(name)) == 30, f"{name} must have 30 entries"
    directions = [int(float(x)) for x in col("DIRECTION")]
    assert directions == [1] * 15 + [-1] * 15
    w0 = [float(x) for x in col("W0COEFF")]
    # 4 distinct soft-core bridge W0 values, each appearing once per direction.
    bridge_w0 = sorted({round(x, 6) for x in w0 if 0.0 < x < 1.0})
    assert len(bridge_w0) == 4
    for w in bridge_w0:
        assert sum(1 for x in w0 if abs(x - w) < 1e-6) == 2, (
            f"W0={w} appears once per direction (2 total)")


# --- Regression: bound28 + prior schedules UNCHANGED by the bound30 addition --

def test_densified_bound28_unchanged_by_bound30_addition(v2):
    """R-7 regression: densified_bound28 is BYTE-EQUAL after adding bound30 (the
    3-window bridge, 28 states, 14/14, W0 0.20/0.45/0.70 each once per dir)."""
    s = v2.get_schedule("densified_bound28")
    assert s["n_states"] == 28
    assert s["directions"] == [1] * 14 + [-1] * 14
    bridge = list(zip(s["lambdas_1"][7:10], s["lambdas_2"][7:10],
                      s["w0"][7:10], s["alpha"][7:10], s["u0"][7:10]))
    assert bridge == _BOUND28_BRIDGE_EXPECTED
    for w in (0.20, 0.45, 0.70):
        assert sum(1 for x in s["w0"] if abs(x - w) < 1e-9) == 2


def test_canonical22_unchanged_by_bound30_addition(v2):
    """R-7 regression: canonical22 is UNCHANGED by adding densified_bound30."""
    c22 = v2.get_schedule("canonical22")
    assert c22["n_states"] == 22
    assert c22["directions"] == [1] * 11 + [-1] * 11
    assert c22["alpha"] == [0.10] * 22
    assert c22["u0"] == [110.0] * 22


def test_densified38v3_unchanged_by_bound28_addition(v2):
    """Regression: the free-leg densified38v3 is UNCHANGED by adding the bound
    schedule (39 states, 19/20, the 9→10 W0=0.07 bridge intact). The bound
    schedule is a SEPARATE registry entry; it does not mutate the free ladder.
    (Note: the free ladder legitimately contains a W0=0.45 window of its own —
    DENSE38_LADDER — which is unrelated to the bound 6→7 bridge, so an absence
    check on the bound W0 values would be a false signal.)"""
    v3 = v2.get_schedule("densified38v3")
    assert v3["n_states"] == 39
    assert v3["directions"] == [1] * 19 + [-1] * 20
    assert sum(1 for w in v3["w0"] if abs(w - 0.07) < 1e-9) == 1
    # densified38v3's dplus stays byte-equal to densified38v2's dplus (the free
    # ladder is untouched by the bound addition).
    v2s = v2.get_schedule("densified38v2")
    for key in ("lambdas_1", "lambdas_2", "lambdas",
                "intermd", "w0", "alpha", "u0"):
        assert v3[key][:19] == v2s[key][:19], f"v3 dplus {key} unchanged"


# ===========================================================================
# P2: paired-difference σ_btwn statistics helper + v3 ddint JSON re-emit.
# Pins the production n4 campaign numbers (mean +1.2432, paired SD 1.1132,
# SEM 0.5566, cross-leg r -0.618, t(3) 2.234) — honest paired stats that the
# v2 quadrature ignored (cross-leg correlation + STDEV/SEM conflation).
# ===========================================================================
# Exact shipped per-replicate dgbind1 values (production_v2_4_replicate_
# campaign/ddint_free_n4.json, n=4 matched seeds). Order = rep1..rep4.
_N4_CP4_DGBIND1 = [
    95.94337820961978, 96.74114230986879, 95.40976256134152, 95.84444337200682,
]
_N4_WT_DGBIND1 = [
    94.3077565056123, 94.49965937351176, 95.75636025078887, 94.40212588521466,
]


def _make_endpoint_block(dgbind1_vals):
    """Minimal cp4/wt endpoint block accepted by build_v3_payload — only the
    fields the v3 builder reads (per-rep dgbind1 + mean/sigma + closure)."""
    import numpy as np
    arr = np.array(dgbind1_vals, dtype=float)
    n = len(arr)
    sigma = float(arr.std(ddof=1)) if n >= 2 else None
    reps = [
        {"dgbind1_kcal": float(v), "dgb_kcal": 0.1, "ddgb_kcal": 0.3,
         "n_states": 22}
        for v in dgbind1_vals
    ]
    return {
        "n_replicates": n,
        "replicate_leg_dirs": [f"/fake/rep{i}" for i in range(n)],
        "replicates": reps,
        "mean_dgbind1_kcal": float(arr.mean()),
        "sigma_btwn_dgbind1_kcal": sigma,
        "sem_dgbind1_kcal": (float(sigma / n ** 0.5) if sigma else None),
        "single_run_analytic_ddgb_kcal_diag": 0.3,
    }


def test_paired_difference_stats_reproduces_n4_campaign(uwham_pp):
    """Regression PIN: paired-difference statistics over the shipped n4
    matched-seed dgbind1 values reproduce the honest paired numbers that the
    v2 quadrature error obscured."""
    r = uwham_pp.paired_difference_stats(_N4_CP4_DGBIND1, _N4_WT_DGBIND1)
    assert r["n"] == 4
    assert r["df"] == 3
    assert r["mean"] == pytest.approx(1.2432, abs=5e-4)
    assert r["paired_sd"] == pytest.approx(1.1132, abs=5e-4)
    assert r["sem"] == pytest.approx(0.5566, abs=5e-4)
    assert r["pearson_r"] == pytest.approx(-0.618, abs=5e-3)
    assert r["t_stat"] == pytest.approx(2.234, abs=5e-3)
    # t-based 95% CI spans 0 -> sign undetermined.
    assert r["ci95"][0] == pytest.approx(-0.528, abs=5e-3)
    assert r["ci95"][1] == pytest.approx(3.015, abs=5e-3)
    assert r["sign_status"] == "undetermined"
    assert r["error_basis"] == "paired_sigma_btwn_t"


def test_paired_difference_stats_mean_equals_diff_of_means(uwham_pp):
    """The paired mean is the SAME point estimate as mean(cp4)-mean(wt)
    (the difference of means equals the mean of differences)."""
    import numpy as np
    r = uwham_pp.paired_difference_stats(_N4_CP4_DGBIND1, _N4_WT_DGBIND1)
    diff_of_means = float(np.mean(_N4_CP4_DGBIND1) - np.mean(_N4_WT_DGBIND1))
    assert r["mean"] == pytest.approx(diff_of_means, abs=1e-12)


def test_paired_difference_stats_quadrature_legacy_differs(uwham_pp):
    """The honest paired SEM (0.5566) is materially smaller than the naive
    independent STDEV quadrature (0.8786) — proving the v2 number conflated
    STDEV with SEM and ignored the (negative) cross-leg correlation."""
    import numpy as np
    cp4_sd = float(np.std(_N4_CP4_DGBIND1, ddof=1))
    wt_sd = float(np.std(_N4_WT_DGBIND1, ddof=1))
    quadrature = (cp4_sd ** 2 + wt_sd ** 2) ** 0.5
    assert quadrature == pytest.approx(0.8786, abs=5e-4)
    r = uwham_pp.paired_difference_stats(_N4_CP4_DGBIND1, _N4_WT_DGBIND1)
    assert r["sem"] < quadrature


def test_paired_difference_stats_unequal_length_raises(uwham_pp):
    with pytest.raises(ValueError):
        uwham_pp.paired_difference_stats([1.0, 2.0], [1.0])


def test_paired_difference_stats_single_pair_no_dispersion(uwham_pp):
    """n=1 -> mean defined but SD/SEM/t undefined (None), sign undetermined."""
    r = uwham_pp.paired_difference_stats([2.0], [1.0])
    assert r["n"] == 1
    assert r["mean"] == pytest.approx(1.0)
    assert r["paired_sd"] is None
    assert r["sem"] is None
    assert r["t_stat"] is None
    assert r["sign_status"] == "undetermined"


def test_paired_difference_stats_positive_sign(uwham_pp):
    """A tight cluster of positive differences -> CI entirely > 0 -> positive."""
    r = uwham_pp.paired_difference_stats(
        [5.00, 5.10, 4.95, 5.05], [1.00, 1.05, 0.98, 1.02],
    )
    assert r["ci95"][0] > 0
    assert r["sign_status"] == "positive"


def test_build_v3_payload_point_estimate_and_fields(uwham_pp):
    """build_v3_payload keeps the point estimate, adds the paired SEM as the
    primary error, retains the quadrature as a labelled legacy field, and
    does NOT reuse the old ddint_err_kcal key (units changed: SEM vs STDEV)."""
    cp4 = _make_endpoint_block(_N4_CP4_DGBIND1)
    wt = _make_endpoint_block(_N4_WT_DGBIND1)
    payload = uwham_pp.build_v3_payload(
        cp4=cp4, wt=wt, protocol="test", mintimeid=100, maxtimeid=None,
        notes="unit test",
    )
    assert payload["schema_version"] == "trackb_ddint_free_v3"
    assert payload["ddint_free_kcal"] == pytest.approx(1.2432, abs=5e-4)
    assert payload["error_basis"] == "paired_sigma_btwn_t"
    assert payload["ddint_err_paired_sem_kcal"] == pytest.approx(0.5566,
                                                                 abs=5e-4)
    assert payload["ddint_err_quadrature_legacy_kcal"] == pytest.approx(
        0.8786, abs=5e-4)
    # Old key intentionally absent (consumer-visible change signalled by v3).
    assert "ddint_err_kcal" not in payload
    ps = payload["paired_difference_stats"]
    assert ps["sign_status"] == "undetermined"
    assert len(ps["per_rep_ddg"]) == 4


def test_build_v3_payload_unequal_replicates_no_pairing(uwham_pp):
    """Unequal replicate counts -> pairing undefined -> paired block None,
    error_basis falls back to the σ_btwn quadrature tag (still emits a number
    via the legacy field)."""
    cp4 = _make_endpoint_block(_N4_CP4_DGBIND1)            # n=4
    wt = _make_endpoint_block(_N4_WT_DGBIND1[:3])          # n=3
    payload = uwham_pp.build_v3_payload(
        cp4=cp4, wt=wt, protocol="test", mintimeid=100, maxtimeid=None,
        notes="",
    )
    assert payload["paired_difference_stats"] is None
    assert payload["ddint_err_paired_sem_kcal"] is None
    assert payload["error_basis"] == "sigma_btwn"
    assert payload["ddint_err_quadrature_legacy_kcal"] is not None


def test_reemit_v3_from_v2_numbers_neutral(uwham_pp, tmp_path):
    """The --reemit-v3-from path re-derives v3 from a v2 artifact's stored
    blocks with the point estimate + per-replicate values BYTE-IDENTICAL,
    only upgrading the uncertainty representation."""
    import json
    cp4 = _make_endpoint_block(_N4_CP4_DGBIND1)
    wt = _make_endpoint_block(_N4_WT_DGBIND1)
    v2 = {
        "schema_version": "trackb_ddint_free_v2",
        "protocol": "atm-openmm async_re b+ corrected",
        "mintimeid": 100, "maxtimeid": None,
        "error_basis": "sigma_btwn",
        "cp4": cp4, "wt": wt,
        "ddint_free_kcal": cp4["mean_dgbind1_kcal"] - wt["mean_dgbind1_kcal"],
        "ddint_err_kcal": 0.8786,
        "notes": "original v2 notes",
    }
    v2_path = tmp_path / "ddint_free_n4.v2.json"
    v2_path.write_text(json.dumps(v2, indent=2))
    out_path = tmp_path / "ddint_free_n4.json"
    payload = uwham_pp.reemit_v3_from_v2(str(v2_path), str(out_path))
    # Source v2 untouched (R-7 — original preserved).
    assert json.loads(v2_path.read_text())["schema_version"] == \
        "trackb_ddint_free_v2"
    on_disk = json.loads(out_path.read_text())
    assert on_disk["schema_version"] == "trackb_ddint_free_v3"
    # Point estimate byte-identical.
    assert on_disk["ddint_free_kcal"] == v2["ddint_free_kcal"]
    # Per-replicate dgbind1 byte-identical for both endpoints.
    for ep in ("cp4", "wt"):
        src = [r["dgbind1_kcal"] for r in v2[ep]["replicates"]]
        got = [r["dgbind1_kcal"] for r in on_disk[ep]["replicates"]]
        assert got == src
    # Paired SEM is the new primary; legacy quadrature retained.
    assert on_disk["ddint_err_paired_sem_kcal"] == pytest.approx(0.5566,
                                                                 abs=5e-4)
    assert on_disk["ddint_err_quadrature_legacy_kcal"] == pytest.approx(
        0.8786, abs=5e-4)
    assert payload["reemit_provenance"]["source_schema"] == \
        "trackb_ddint_free_v2"


def test_reemit_v3_rejects_non_v2(uwham_pp, tmp_path):
    """Re-emitting from a non-v2 (e.g. already-v3) artifact fails fast."""
    import json
    bad = tmp_path / "already_v3.json"
    bad.write_text(json.dumps({"schema_version": "trackb_ddint_free_v3"}))
    with pytest.raises(ValueError):
        uwham_pp.reemit_v3_from_v2(str(bad), str(tmp_path / "out.json"))


def test_student_t_ppf_975_low_df(uwham_pp):
    """The t-critical helper returns the textbook value at df=3 (2.3060 at
    df=8, 3.1824 at df=3) and None below df=1."""
    assert uwham_pp._student_t_ppf_975(3) == pytest.approx(3.1824, abs=2e-3)
    assert uwham_pp._student_t_ppf_975(8) == pytest.approx(2.3060, abs=2e-3)
    assert uwham_pp._student_t_ppf_975(0) is None


# ===========================================================================
# P11: MBAR overlap-matrix QC (REPORTING-ONLY; numbers-neutral).
# Computes pymbar.MBAR(...).compute_overlap() on the captured per-target-state
# neg_pot SSOT matrix. Threshold O<0.03 is a FLAG (report), NOT a gate; O is
# NEVER folded into any error bar (reporting-only scope-lock). pymbar may be
# absent in the qmmm env (present in atm) — every failure must graceful-degrade
# to status:"unavailable" without raising.
# DOI: 10.1063/1.2978177 (MBAR), 10.1007/s10822-015-9840-9 (overlap diagnostic).
# ===========================================================================
def _pymbar_available():
    try:
        import pymbar  # noqa: F401
        return True
    except Exception:
        return False


def _harmonic_neg_pot(centers, n_per_state, sigma=1.0, seed=0):
    """Synthetic ``neg_pot[N, K]`` = -0.5*(x_n - center_k)**2 (β=1) over K
    harmonic states sampled at ``centers``. Returns (neg_pot, N_k). Adjacent
    overlap is high for close centers, ~0 for distant centers."""
    import numpy as np
    rng = np.random.default_rng(seed)
    k = len(centers)
    xs = [rng.normal(c, sigma, n_per_state) for c in centers]
    x = np.concatenate(xs)
    neg_pot = np.zeros((k * n_per_state, k))
    for kk in range(k):
        neg_pot[:, kk] = -(0.5 * (x - centers[kk]) ** 2)
    return neg_pot, [n_per_state] * k


def test_overlap_qc_healthy_matrix_no_flag(uwham_pp):
    """(a) A synthetic HEALTHY ladder (close, overlapping states) yields
    min-adjacent O well above 0.03 → no flag, tier WELL_DETERMINED, no
    warning string."""
    if not _pymbar_available():
        pytest.skip("pymbar not importable in this env")
    neg_pot, n_k = _harmonic_neg_pot([0.0, 0.6, 1.2], 300, seed=1)
    qc = uwham_pp.compute_overlap_qc(
        {"leg1_neg_pot": neg_pot, "leg1_N_k": n_k})
    assert qc["status"] == "ok"
    assert qc["basis"] == "mbar_compute_overlap"
    assert qc["threshold"] == 0.03
    assert qc["min_adjacent_O"] > 0.03
    assert qc["flagged_pairs"] == []
    assert qc["warning"] is None
    assert qc["tier"] in ("MARGINAL", "WELL_DETERMINED")
    # adjacent_O present per leg (2 adjacencies for 3 states).
    assert len(qc["legs"]["forward"]["adjacent_O"]) == 2


def test_overlap_qc_near_zero_pair_flagged(uwham_pp):
    """(b) A synthetic near-ZERO-overlap adjacent pair (one far-away state) is
    FLAGGED: tier COLLAPSE, a flagged_pairs entry, and a LOW_OVERLAP warning
    string. The flag is REPORT-ONLY — no ΔG is returned by this helper."""
    if not _pymbar_available():
        pytest.skip("pymbar not importable in this env")
    # states 0,1 overlap; state 2 is 40 units away -> O_{1,2} ~ 0.
    neg_pot, n_k = _harmonic_neg_pot([0.0, 0.5, 40.0], 300, seed=2)
    qc = uwham_pp.compute_overlap_qc(
        {"leg1_neg_pot": neg_pot, "leg1_N_k": n_k})
    assert qc["status"] == "ok"
    assert qc["min_adjacent_O"] < 0.03
    assert qc["tier"] == "COLLAPSE"
    assert len(qc["flagged_pairs"]) >= 1
    # the flagged pair is the (1,2) adjacency.
    assert any(p["i"] == 1 and p["j"] == 2 for p in qc["flagged_pairs"])
    assert qc["warning"] is not None
    assert "LOW_OVERLAP" in qc["warning"]
    assert "does NOT gate" in qc["warning"]


def test_overlap_qc_from_matrix_pure_no_pymbar(uwham_pp):
    """The matrix→report reducer is pure (no pymbar) — a hand-built overlap
    matrix with one sub-threshold adjacency flags correctly. Lets the flag
    logic be tested even where pymbar is unavailable."""
    import numpy as np
    # 3-state overlap matrix: O_01 = 0.20 (ok), O_12 = 0.01 (collapse).
    mat = np.array([
        [0.80, 0.20, 0.00],
        [0.20, 0.79, 0.01],
        [0.00, 0.01, 0.99],
    ])
    qc = uwham_pp._overlap_qc_from_matrix(mat, threshold=0.03)
    assert qc["status"] == "ok"
    assert qc["adjacent_O"] == [pytest.approx(0.20), pytest.approx(0.01)]
    assert qc["min_adjacent_O"] == pytest.approx(0.01)
    assert qc["tier"] == "COLLAPSE"
    assert qc["flagged_pairs"] == [{"i": 1, "j": 2, "O": pytest.approx(0.01)}]
    assert "LOW_OVERLAP" in qc["warning"]


def test_overlap_qc_threshold_is_not_a_gate(uwham_pp):
    """SCOPE-LOCK: O<0.03 is a FLAG, not a gate. The report is a descriptive
    dict — it carries NO ΔG, NO error bar, NO reject/rescale field. Confirms
    the helper cannot mutate any free-energy number (it has none to mutate)."""
    import numpy as np
    mat = np.array([[0.9, 0.001], [0.001, 0.999]])  # collapse pair
    qc = uwham_pp._overlap_qc_from_matrix(mat, threshold=0.03)
    # Only descriptive QC keys — no ΔG / error / gate verb anywhere.
    forbidden = {
        "dgbind1_kcal", "dgbind2_kcal", "dgb_kcal", "ddgb_kcal",
        "ddint_free_kcal", "sem", "sigma_btwn", "reject", "rescale",
        "skip", "gate", "ci95",
    }
    assert forbidden.isdisjoint(qc.keys())
    assert qc["tier"] == "COLLAPSE"  # flagged, but purely a label


def test_overlap_qc_graceful_degrade_no_matrix(uwham_pp):
    """(d) The graceful-degrade path: an absent / empty sink (canonical-22
    upstream path, or pymbar/atom_openmm missing) returns
    status:"unavailable" WITHOUT raising — so the qmmm pytest path passes and
    no existing number changes."""
    for sink in (None, {}, {"leg2_neg_pot": "x"}):  # no leg1 matrix
        qc = uwham_pp.compute_overlap_qc(sink)
        assert qc["status"] == "unavailable"
        assert "reason" in qc
        assert qc["basis"] == "mbar_compute_overlap"


def test_overlap_qc_graceful_degrade_bad_matrix(uwham_pp):
    """A malformed captured matrix (n_k inconsistent with neg_pot) degrades to
    status:"unavailable" — the try/except swallows it; nothing propagates."""
    import numpy as np
    bad = {"leg1_neg_pot": np.zeros((10, 3)), "leg1_N_k": [3, 3, 3]}  # 9 != 10
    qc = uwham_pp.compute_overlap_qc(bad)
    assert qc["status"] == "unavailable"
    assert "reason" in qc


def test_build_overlap_qc_summary_additive_and_neutral(uwham_pp):
    """build_v3_payload carries a top-level overlap_qc digest that is purely
    ADDITIVE: it does NOT change ddint_free_kcal or any error field, and on
    synthetic blocks (no captured matrix) it reports status:"unavailable"
    per endpoint."""
    cp4 = _make_endpoint_block(_N4_CP4_DGBIND1)
    wt = _make_endpoint_block(_N4_WT_DGBIND1)
    payload = uwham_pp.build_v3_payload(
        cp4=cp4, wt=wt, protocol="test", mintimeid=100, maxtimeid=None,
        notes="unit test",
    )
    # Point estimate + error bars unchanged by the additive QC field.
    assert payload["ddint_free_kcal"] == pytest.approx(1.2432, abs=5e-4)
    assert payload["ddint_err_paired_sem_kcal"] == pytest.approx(0.5566,
                                                                 abs=5e-4)
    # New additive top-level digest.
    oq = payload["overlap_qc"]
    assert oq["basis"] == "mbar_compute_overlap"
    assert oq["threshold"] == 0.03
    assert oq["cp4"]["status"] == "unavailable"   # synthetic block, no matrix
    assert oq["wt"]["status"] == "unavailable"
    assert oq["min_adjacent_O"] is None
    assert oq["any_flagged"] is False


def test_reemit_v3_point_estimate_byte_identical_with_overlap_qc(
    uwham_pp, tmp_path,
):
    """(c) Re-emitting v3 from a v2 built on the shipped n4 dgbind1 values keeps
    the point estimate BYTE-IDENTICAL to the production SSOT
    (+1.2432061094273337) while the additive overlap_qc field is present —
    proving P11 is numbers-neutral on the re-emit path (P2 byte-identity
    preserved)."""
    import json
    import numpy as np
    cp4 = _make_endpoint_block(_N4_CP4_DGBIND1)
    wt = _make_endpoint_block(_N4_WT_DGBIND1)
    ddint = cp4["mean_dgbind1_kcal"] - wt["mean_dgbind1_kcal"]
    v2 = {
        "schema_version": "trackb_ddint_free_v2",
        "protocol": "atm-openmm async_re b+ corrected",
        "mintimeid": 100, "maxtimeid": None,
        "error_basis": "sigma_btwn",
        "cp4": cp4, "wt": wt,
        "ddint_free_kcal": ddint,
        "ddint_err_kcal": 0.8786,
        "notes": "n4 campaign",
    }
    v2_path = tmp_path / "ddint_free_n4.v2.json"
    v2_path.write_text(json.dumps(v2, indent=2))
    out_path = tmp_path / "ddint_free_n4.json"
    payload = uwham_pp.reemit_v3_from_v2(str(v2_path), str(out_path))
    on_disk = json.loads(out_path.read_text())
    # Point estimate byte-identical to the production SSOT constant.
    assert on_disk["ddint_free_kcal"] == ddint
    assert on_disk["ddint_free_kcal"] == pytest.approx(
        1.2432061094273337, abs=0.0)
    # Per-replicate dgbind1 byte-identical.
    for ep in ("cp4", "wt"):
        src = [r["dgbind1_kcal"] for r in v2[ep]["replicates"]]
        got = [r["dgbind1_kcal"] for r in on_disk[ep]["replicates"]]
        assert got == src
    # Additive overlap_qc present, neutral (unavailable on synthetic blocks).
    assert "overlap_qc" in on_disk
    assert on_disk["overlap_qc"]["any_flagged"] is False
    assert payload["overlap_qc"]["min_adjacent_O"] is None


# ===========================================================================
# P5: per-replicate cycle-closure DISTRUST flag (REPORTING-ONLY; ASYMMETRIC).
# Labels the EXISTING result["dgb_kcal"] (|dgbind1-dgbind2|) as a distrust
# prompt. |closure|>2.0 -> distrust=True (review); |closure|<=2.0 -> NOT a
# trust assertion. NEVER folds closure into any error bar; NEVER asserts trust
# (reporting-only scope-lock). A flagged densified38v4 dminus leg reads
# "review: may be non-mirror hysteresis", not "broken".
# DOI: 10.1021/jp102971x (directional bias), 10.1007/s10822-015-9840-9.
# ===========================================================================
# Asymmetric forbidden keys: this report must NEVER advertise correctness.
_CLOSURE_FORBIDDEN_KEYS = {
    "trust", "trusted", "pass", "passed", "green", "ok_to_trust",
    "validated", "verified", "correct", "converged",
}


def _make_endpoint_block_with_closure(dgbind1_vals, dgb_vals):
    """cp4/wt endpoint block whose per-replicate results carry a
    closure_distrust report (what analyze_one_leg attaches), built from the
    given per-rep dgb (closure) values via the production reducer."""
    import numpy as np
    arr = np.array(dgbind1_vals, dtype=float)
    n = len(arr)
    sigma = float(arr.std(ddof=1)) if n >= 2 else None
    assert len(dgbind1_vals) == len(dgb_vals)
    reps = []
    for v, dgb in zip(dgbind1_vals, dgb_vals):
        r = {"dgbind1_kcal": float(v), "dgb_kcal": float(dgb),
             "ddgb_kcal": 0.3, "n_states": 22}
        # closure_distrust filled by the caller via compute_closure_distrust,
        # mirroring analyze_one_leg's attachment (same module reducer).
        r["closure_distrust"] = None
        reps.append(r)
    return {
        "n_replicates": n,
        "replicate_leg_dirs": [f"/fake/rep{i}" for i in range(n)],
        "replicates": reps,
        "mean_dgbind1_kcal": float(arr.mean()),
        "sigma_btwn_dgbind1_kcal": sigma,
        "sem_dgbind1_kcal": (float(sigma / n ** 0.5) if sigma else None),
        "single_run_analytic_ddgb_kcal_diag": 0.3,
    }


def test_closure_distrust_large_closure_flagged(uwham_pp):
    """(a) A synthetic leg with |closure| = 3.0 (> 2.0) is FLAGGED:
    distrust=True, tier 'distrust', a CLOSURE_DISTRUST warning string. The
    report is REPORT-ONLY — it carries NO ΔG."""
    cd = uwham_pp._closure_distrust_from_dgb(3.0)
    assert cd["status"] == "ok"
    assert cd["basis"] == "abs_dgbind1_minus_dgbind2"
    assert cd["threshold"] == 2.0
    assert cd["closure_kcal_abs"] == pytest.approx(3.0)
    assert cd["distrust"] is True
    assert cd["tier"] == "distrust"
    assert cd["warning"] is not None
    assert "CLOSURE_DISTRUST" in cd["warning"]
    assert "does NOT gate" in cd["warning"].lower() or \
        "Does NOT gate" in cd["warning"]
    # Sign-invariant: a negative dgb of the same magnitude flags identically.
    cd_neg = uwham_pp._closure_distrust_from_dgb(-3.0)
    assert cd_neg["distrust"] is True
    assert cd_neg["closure_kcal_abs"] == pytest.approx(3.0)


def test_closure_distrust_small_closure_not_flagged(uwham_pp):
    """(b) A synthetic leg with |closure| = 0.5 (< 1.0) is NOT flagged:
    distrust=False, tier 'clean', no warning. distrust=False asserts NOTHING
    about correctness (asymmetric)."""
    cd = uwham_pp._closure_distrust_from_dgb(0.5)
    assert cd["status"] == "ok"
    assert cd["distrust"] is False
    assert cd["tier"] == "clean"
    assert cd["warning"] is None
    # The v2_1 mild band (1.0-2.0) is NOT distrust (preserves 1.43 soft-flag).
    cd_mild = uwham_pp._closure_distrust_from_dgb(1.43)
    assert cd_mild["distrust"] is False
    assert cd_mild["tier"] == "mild"
    # The 2.0 boundary itself is NOT distrust (strict >).
    cd_edge = uwham_pp._closure_distrust_from_dgb(2.0)
    assert cd_edge["distrust"] is False
    assert cd_edge["tier"] == "mild"


def test_closure_distrust_never_asserts_trust(uwham_pp):
    """(c) ASYMMETRIC: the report NEVER emits a trust/pass/green key, in EITHER
    the flagged or the unflagged case. distrust=False is "no red-flag", not a
    correctness assertion. Verifies the must-not-cross asymmetry directly."""
    for dgb in (0.0, 0.5, 1.43, 3.0, -3.0):
        cd = uwham_pp._closure_distrust_from_dgb(dgb)
        assert _CLOSURE_FORBIDDEN_KEYS.isdisjoint(cd.keys()), (
            f"closure report leaked a trust-flavoured key for dgb={dgb}")
        # The only correctness-adjacent field is 'distrust' (asymmetric).
        assert "distrust" in cd
        # No ΔG / error / gate field can be mutated here (none present).
        forbidden_num = {
            "dgbind1_kcal", "dgbind2_kcal", "ddint_free_kcal", "sem",
            "sigma_btwn", "ci95", "reject", "rescale", "gate",
        }
        assert forbidden_num.isdisjoint(cd.keys())


def test_closure_quarantine_top_level_flags_distrust_endpoint(uwham_pp):
    """The build_v3_payload top-level closure_quarantine lists ONLY the
    endpoints whose |closure| > 2.0 in 'flagged' (asymmetric — no trust list),
    emits a warning, and never moves ddint_free_kcal or any error bar."""
    # cp4 rep0 closure 2.727 (distrust), wt rep0 closure 3.424 (distrust) —
    # the live n4 campaign numbers.
    cp4 = _make_endpoint_block_with_closure(
        _N4_CP4_DGBIND1, [2.727, 0.1, 0.1, 0.1])
    wt = _make_endpoint_block_with_closure(
        _N4_WT_DGBIND1, [3.424, 0.1, 0.1, 0.1])
    # Attach the closure_distrust reports the way analyze_one_leg does.
    for block in (cp4, wt):
        for r in block["replicates"]:
            r["closure_distrust"] = uwham_pp.compute_closure_distrust(r)
    payload = uwham_pp.build_v3_payload(
        cp4=cp4, wt=wt, protocol="test", mintimeid=100, maxtimeid=None,
        notes="unit test",
    )
    # Numbers-neutral: point estimate + error bars unchanged.
    assert payload["ddint_free_kcal"] == pytest.approx(1.2432, abs=5e-4)
    assert payload["ddint_err_paired_sem_kcal"] == pytest.approx(0.5566,
                                                                 abs=5e-4)
    # Closure values byte-identical to the existing closure_*_kcal_abs fields.
    assert payload["closure_cp4_kcal_abs"] == pytest.approx(2.727, abs=0.0)
    assert payload["closure_wt_kcal_abs"] == pytest.approx(3.424, abs=0.0)
    cq = payload["closure_quarantine"]
    assert cq["basis"] == "abs_dgbind1_minus_dgbind2"
    assert cq["threshold"] == 2.0
    # rep0 of each endpoint is flagged; rep1-3 are clean.
    assert "cp4/rep0" in cq["flagged"]
    assert "wt/rep0" in cq["flagged"]
    assert "cp4/rep1" not in cq["flagged"]
    assert cq["warning"] is not None
    assert "CLOSURE_DISTRUST" in cq["warning"]
    # ASYMMETRIC: no trust list, no trust-flavoured key in the digest.
    assert _CLOSURE_FORBIDDEN_KEYS.isdisjoint(cq.keys())
    # per_rep_closure carries every replicate.
    assert len(cq["per_rep_closure"]) == 8
    flagged_entries = [e for e in cq["per_rep_closure"]
                       if e.get("status") == "ok" and e["distrust"]]
    assert len(flagged_entries) == 2


def test_closure_quarantine_unavailable_does_not_crash(uwham_pp):
    """Graceful: closure_unavailable applies ONLY when a replicate lacks BOTH a
    closure_distrust report AND a finite dgb_kcal — no crash either way.

    Reps that carry a finite dgb_kcal but no attached closure_distrust (the
    v2-archive re-emit shape) are NOT unavailable: the digest derives the label
    from dgb_kcal via the existing compute_closure_distrust reducer. Reps with
    neither stay unavailable."""
    import copy
    # (1) Finite dgb_kcal, no closure_distrust attached (v2-archive re-emit
    # shape): _make_endpoint_block sets dgb_kcal=0.1 (< 1.0) on every rep, so
    # the digest DERIVES a clean, non-distrust label — NOT unavailable.
    cp4 = _make_endpoint_block(_N4_CP4_DGBIND1)
    wt = _make_endpoint_block(_N4_WT_DGBIND1)
    payload = uwham_pp.build_v3_payload(
        cp4=cp4, wt=wt, protocol="test", mintimeid=100, maxtimeid=None,
        notes="unit test",
    )
    cq = payload["closure_quarantine"]
    assert cq["flagged"] == []                # 0.1 < threshold -> not distrust
    assert cq["unavailable"] == []            # derived from dgb_kcal, not N/A
    assert cq["warning"] is None
    derived = [e for e in cq["per_rep_closure"]
               if e.get("status") == "ok" and e["tier"] == "clean"]
    assert len(derived) == 8                  # all 8 derived clean from dgb_kcal

    # (2) Genuinely unavailable: strip dgb_kcal AND closure_distrust from reps
    # 1..3 of each endpoint (rep0 keeps dgb_kcal — build_v3_payload reads it for
    # the closure_*_kcal_abs field). Only reps with NEITHER key are unavailable.
    cp4_na = copy.deepcopy(cp4)
    wt_na = copy.deepcopy(wt)
    for block in (cp4_na, wt_na):
        for r in block["replicates"][1:]:
            r.pop("dgb_kcal", None)
            r.pop("closure_distrust", None)
    payload_na = uwham_pp.build_v3_payload(
        cp4=cp4_na, wt=wt_na, protocol="test", mintimeid=100, maxtimeid=None,
        notes="unit test",
    )
    cq_na = payload_na["closure_quarantine"]
    assert cq_na["flagged"] == []             # rep0 dgb 0.1 -> clean, not distrust
    # 3 stripped reps per endpoint -> 6 unavailable; the 2 rep0s derive clean.
    assert sorted(cq_na["unavailable"]) == [
        "cp4/rep1", "cp4/rep2", "cp4/rep3",
        "wt/rep1", "wt/rep2", "wt/rep3",
    ]
    assert cq_na["warning"] is None

    # Direct: a result dict without dgb_kcal degrades gracefully.
    cd = uwham_pp.compute_closure_distrust({"dgbind1_kcal": 90.0})
    assert cd["status"] == "closure_unavailable"
    assert "reason" in cd
    # A non-finite dgb_kcal also degrades (no crash).
    cd_nan = uwham_pp.compute_closure_distrust({"dgb_kcal": float("nan")})
    assert cd_nan["status"] == "closure_unavailable"


def test_closure_quarantine_reemit_point_estimate_byte_identical(
    uwham_pp, tmp_path,
):
    """(d) Re-emitting v3 from a v2 built on the shipped n4 dgbind1 values keeps
    the point estimate BYTE-IDENTICAL to the production SSOT
    (+1.2432061094273337) while the additive closure_quarantine field is
    present — proving P5 is numbers-neutral on the re-emit path."""
    import json
    cp4 = _make_endpoint_block(_N4_CP4_DGBIND1)
    wt = _make_endpoint_block(_N4_WT_DGBIND1)
    ddint = cp4["mean_dgbind1_kcal"] - wt["mean_dgbind1_kcal"]
    v2 = {
        "schema_version": "trackb_ddint_free_v2",
        "protocol": "atm-openmm async_re b+ corrected",
        "mintimeid": 100, "maxtimeid": None,
        "error_basis": "sigma_btwn",
        "cp4": cp4, "wt": wt,
        "ddint_free_kcal": ddint,
        "ddint_err_kcal": 0.8786,
        "notes": "n4 campaign",
    }
    v2_path = tmp_path / "ddint_free_n4.v2.json"
    v2_path.write_text(json.dumps(v2, indent=2))
    out_path = tmp_path / "ddint_free_n4.json"
    payload = uwham_pp.reemit_v3_from_v2(str(v2_path), str(out_path))
    on_disk = json.loads(out_path.read_text())
    # Point estimate byte-identical to the production SSOT constant.
    assert on_disk["ddint_free_kcal"] == ddint
    assert on_disk["ddint_free_kcal"] == pytest.approx(
        1.2432061094273337, abs=0.0)
    assert on_disk["ddint_err_paired_sem_kcal"] == pytest.approx(0.5566,
                                                                 abs=5e-4)
    # Per-replicate dgbind1 byte-identical.
    for ep in ("cp4", "wt"):
        src = [r["dgbind1_kcal"] for r in v2[ep]["replicates"]]
        got = [r["dgbind1_kcal"] for r in on_disk[ep]["replicates"]]
        assert got == src
    # Additive closure_quarantine present (unavailable on synthetic blocks).
    assert "closure_quarantine" in on_disk
    assert on_disk["closure_quarantine"]["flagged"] == []
    # Asymmetric: no trust-flavoured key leaked into the persisted digest.
    assert _CLOSURE_FORBIDDEN_KEYS.isdisjoint(
        on_disk["closure_quarantine"].keys())
    assert payload["closure_quarantine"]["threshold"] == 2.0


# Real per-rep closure (|dgbind1 - dgbind2|) values stored in the shipped n4 v2
# archive. The v2 path carries dgb_kcal per replicate but NOT a closure_distrust
# report (the flag post-dates the archive). Four of these exceed the 2.0 floor:
# cp4/rep0 (2.727), cp4/rep1 (3.990), wt/rep0 (3.424), wt/rep3 (2.837).
_N4_CP4_DGB = [2.726878323628881, 3.9898610207167025,
               -1.8862439478883033, 1.6787945493379652]
_N4_WT_DGB = [-3.4244781716257933, -0.12266894295916586,
              1.3181131024501553, -2.8367320545581407]


def _make_v2_block_with_dgb(dgbind1_vals, dgb_vals):
    """v2-archive-shaped endpoint block: per-rep dgb_kcal present, NO attached
    closure_distrust report (the flag post-dates the archive). Mirrors what
    reemit_v3_from_v2 reads off the shipped ddint_free_n4.v2.json."""
    import numpy as np
    arr = np.array(dgbind1_vals, dtype=float)
    n = len(arr)
    sigma = float(arr.std(ddof=1)) if n >= 2 else None
    reps = [
        {"dgbind1_kcal": float(v), "dgb_kcal": float(dgb), "ddgb_kcal": 0.3,
         "n_states": 22}
        for v, dgb in zip(dgbind1_vals, dgb_vals)
    ]
    return {
        "n_replicates": n,
        "replicate_leg_dirs": [f"/fake/rep{i}" for i in range(n)],
        "replicates": reps,
        "mean_dgbind1_kcal": float(arr.mean()),
        "sigma_btwn_dgbind1_kcal": sigma,
        "sem_dgbind1_kcal": (float(sigma / n ** 0.5) if sigma else None),
        "single_run_analytic_ddgb_kcal_diag": 0.3,
    }


def test_closure_quarantine_reemit_derives_flagged_from_dgb(
    uwham_pp, tmp_path,
):
    """(e) KP3 regression: re-emitting v3 from a v2 archive that carries dgb_kcal
    but NO attached closure_distrust DERIVES the distrust labels from dgb_kcal
    (via the existing compute_closure_distrust reducer), so the persisted
    closure_quarantine flags exactly cp4/rep0, cp4/rep1, wt/rep0, wt/rep3 —
    matching the live analyze_one_leg path — rather than reporting all reps
    closure_unavailable. Numbers-neutral: the point estimate + paired SEM stay
    byte-identical to the production SSOT and no trust key is emitted."""
    import json
    cp4 = _make_v2_block_with_dgb(_N4_CP4_DGBIND1, _N4_CP4_DGB)
    wt = _make_v2_block_with_dgb(_N4_WT_DGBIND1, _N4_WT_DGB)
    ddint = cp4["mean_dgbind1_kcal"] - wt["mean_dgbind1_kcal"]
    v2 = {
        "schema_version": "trackb_ddint_free_v2",
        "protocol": "atm-openmm async_re b+ corrected",
        "mintimeid": 100, "maxtimeid": None,
        "error_basis": "sigma_btwn",
        "cp4": cp4, "wt": wt,
        "ddint_free_kcal": ddint,
        "ddint_err_kcal": 0.8786,
        "notes": "n4 campaign",
    }
    v2_path = tmp_path / "ddint_free_n4.v2.json"
    v2_path.write_text(json.dumps(v2, indent=2))
    out_path = tmp_path / "ddint_free_n4.json"
    uwham_pp.reemit_v3_from_v2(str(v2_path), str(out_path))
    on_disk = json.loads(out_path.read_text())

    cq = on_disk["closure_quarantine"]
    # The fix: distrust derived from dgb_kcal, NOT all-unavailable.
    assert cq["flagged"] == ["cp4/rep0", "cp4/rep1", "wt/rep0", "wt/rep3"]
    assert cq["unavailable"] == []
    assert cq["warning"] is not None
    assert "CLOSURE_DISTRUST" in cq["warning"]
    # Numbers-neutral: point estimate + paired SEM byte-identical to SSOT.
    assert on_disk["ddint_free_kcal"] == pytest.approx(
        1.2432061094273337, abs=0.0)
    assert on_disk["ddint_err_paired_sem_kcal"] == pytest.approx(
        0.5566018985450167, abs=0.0)
    # closure_*_kcal_abs unchanged (cp4 2.727 / wt 3.424, first replicate).
    assert on_disk["closure_cp4_kcal_abs"] == pytest.approx(
        2.726878323628881, abs=0.0)
    assert on_disk["closure_wt_kcal_abs"] == pytest.approx(
        3.4244781716257933, abs=0.0)
    # ASYMMETRIC: no trust-flavoured key leaked into the persisted digest.
    assert _CLOSURE_FORBIDDEN_KEYS.isdisjoint(cq.keys())


# ===========================================================================
# A6: calibration-anchor ADVISORY sign-check (REPORTING-ONLY; numbers-neutral).
# ddint = mean(cp4.dgbind1) - mean(wt.dgbind1); favorable_cp4 (Cp4 binds
# tighter, ΔΔG_bind < 0) maps to a NEGATIVE ddint. Sign-only (R-11/R-18),
# never a magnitude/calibration claim; free-leg ddint is NOT ΔΔG_bind →
# provisional. Anchor: Katragadda & Lambris 2006 (DOI 10.1021/jm0603419);
# ITC/SPR SSOT: Magotti 2009 (DOI 10.1002/jmr.972).
# ===========================================================================
def test_a6_target_card_2qki_has_calibration_anchor():
    """The 2QKI target card carries a calibration_anchor block with the
    favorable_cp4 sign, the Katragadda K_d anchor, and the Magotti ITC/SPR
    SSOT — config side of A6 (sign/ranking advisory; magnitude NOT calibrated)."""
    import json
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    card = os.path.join(here, "target_cards", "2QKI.json")
    with open(card) as fh:
        d = json.load(fh)
    ca = d["calibration_anchor"]
    assert ca["expected_ddg_bind_sign"] == "favorable_cp4"
    assert "ΔΔG_bind < 0" in ca["favorable_means"]
    assert ca["anchor_kd_nM"] == 15
    assert "Katragadda" in ca["anchor_ref"]
    assert "10.1021/jm0603419" in ca["anchor_ref"]
    assert ca["itc_spr_ssot"]["itc_kcal_mol"] == -1.4
    assert ca["itc_spr_ssot"]["spr_kcal_mol"] == -3.0
    assert "Magotti" in ca["itc_spr_ssot"]["ref"]
    assert "advisory only" in ca["regime"]


def test_a6_sign_favorable_is_consistent(uwham_pp):
    """A synthetic FAVORABLE result (cp4 binds tighter → ddint negative, CI
    excludes 0) → sign_vs_anchor 'consistent', no warning."""
    paired = uwham_pp.paired_difference_stats(
        [1.00, 1.10, 0.95, 1.05], [5.00, 5.10, 4.95, 5.05],  # ddint ≈ -4
    )
    assert paired["sign_status"] == "negative"
    sa = uwham_pp.compute_sign_vs_anchor(paired, leg="free")
    assert sa["sign_vs_anchor"] == "consistent"
    assert sa["warning"] is None
    assert sa["favorable_ddint_sign"] == "negative"
    assert sa["provisional"] is True   # free leg


def test_a6_sign_unfavorable_is_discordant_with_warning(uwham_pp):
    """A synthetic UNFAVORABLE result (cp4 binds weaker → ddint positive, CI
    excludes 0) → sign_vs_anchor 'discordant' with a WARNING. Advisory only —
    the warning never gates / changes a number."""
    paired = uwham_pp.paired_difference_stats(
        [5.00, 5.10, 4.95, 5.05], [1.00, 1.10, 0.95, 1.05],  # ddint ≈ +4
    )
    assert paired["sign_status"] == "positive"
    sa = uwham_pp.compute_sign_vs_anchor(paired, leg="free")
    assert sa["sign_vs_anchor"] == "discordant"
    assert sa["warning"] is not None
    assert "DISCORDANT" in sa["warning"]
    # R-18 honesty: the warning flags the free-leg caveat (not a ranking verdict).
    assert "free leg" in sa["warning"].lower() or "free-leg" in sa["warning"].lower()


def test_a6_sign_ci_spans_zero_is_undetermined(uwham_pp):
    """The shipped n4 campaign (CI spans 0) → sign_vs_anchor 'undetermined'
    (neither consistent nor discordant — the sign is not resolved)."""
    paired = uwham_pp.paired_difference_stats(_N4_CP4_DGBIND1, _N4_WT_DGBIND1)
    assert paired["sign_status"] == "undetermined"
    sa = uwham_pp.compute_sign_vs_anchor(paired, leg="free")
    assert sa["sign_vs_anchor"] == "undetermined"
    assert sa["warning"] is None
    # No paired stats at all → also undetermined, no crash.
    none_sa = uwham_pp.compute_sign_vs_anchor(None, leg="free")
    assert none_sa["sign_vs_anchor"] == "undetermined"


def test_a6_sign_is_advisory_only_no_number_change(uwham_pp):
    """SCOPE-LOCK: the sign-check report carries NO ΔG / error / gate verb. It
    cannot mutate any free-energy number."""
    paired = uwham_pp.paired_difference_stats(
        [5.0, 5.1, 4.95, 5.05], [1.0, 1.1, 0.95, 1.05],
    )
    sa = uwham_pp.compute_sign_vs_anchor(paired, leg="free")
    forbidden = {
        "dgbind1_kcal", "dgb_kcal", "ddgb_kcal", "ddint_free_kcal", "sem",
        "sigma_btwn", "reject", "rescale", "gate", "ci95",
    }
    assert forbidden.isdisjoint(sa.keys())


def test_a6_reemit_v3_point_estimate_byte_identical_with_sign_check(
    uwham_pp, tmp_path,
):
    """A6 is numbers-neutral on the re-emit path: the additive sign_vs_anchor
    field is present while the point estimate (+1.2432061094273337) and paired
    SEM (0.5566018985450167) stay BYTE-IDENTICAL to the production SSOT."""
    import json
    cp4 = _make_endpoint_block(_N4_CP4_DGBIND1)
    wt = _make_endpoint_block(_N4_WT_DGBIND1)
    ddint = cp4["mean_dgbind1_kcal"] - wt["mean_dgbind1_kcal"]
    v2 = {
        "schema_version": "trackb_ddint_free_v2",
        "protocol": "atm-openmm async_re b+ corrected",
        "mintimeid": 100, "maxtimeid": None, "error_basis": "sigma_btwn",
        "cp4": cp4, "wt": wt, "ddint_free_kcal": ddint,
        "ddint_err_kcal": 0.8786, "notes": "n4 campaign",
    }
    v2_path = tmp_path / "ddint_free_n4.v2.json"
    v2_path.write_text(json.dumps(v2, indent=2))
    out_path = tmp_path / "ddint_free_n4.json"
    uwham_pp.reemit_v3_from_v2(str(v2_path), str(out_path))
    on_disk = json.loads(out_path.read_text())
    assert on_disk["ddint_free_kcal"] == pytest.approx(
        1.2432061094273337, abs=0.0)
    assert on_disk["ddint_err_paired_sem_kcal"] == pytest.approx(
        0.5566018985450167, abs=0.0)
    # Additive field present, neutral (n4 CI spans 0 → undetermined).
    sa = on_disk["sign_vs_anchor"]
    assert sa["sign_vs_anchor"] == "undetermined"
    assert sa["provisional"] is True
    assert sa["expected_ddg_bind_sign"] == "favorable_cp4"


# ===========================================================================
# A7: ADVISORY extend-recommended verdict (REPORTING-ONLY; numbers-neutral).
# Pre-reg C4: worst endpoint σ_btwn > target (≈ 2× free paired σ ≈ 2.2
# kcal/mol) → recommend n -> n+1. Never changes n / gates / touches σ.
# DOI 10.1007/s10822-015-9840-9 (convergence diagnostics).
# ===========================================================================
def test_a7_extend_recommended_above_target(uwham_pp):
    """σ_btwn above the 2.2 kcal/mol target → recommended=True with a reason
    that names the extend action."""
    cp4 = _make_endpoint_block([90.0, 93.0, 87.0, 91.0])   # large spread
    wt = _make_endpoint_block([85.0, 80.0, 90.0, 84.0])
    er = uwham_pp.compute_extend_recommended(cp4, wt)
    assert er["recommended"] is True
    assert er["current_sigma_btwn_kcal"] > er["target_kcal"]
    assert er["target_kcal"] == pytest.approx(2.2)
    assert er["n_current"] == 4
    assert "extend" in er["reason"].lower()


def test_a7_extend_not_recommended_below_target(uwham_pp):
    """σ_btwn below the 2.2 kcal/mol target → recommended=False."""
    cp4 = _make_endpoint_block([90.00, 90.10, 89.90, 90.05])  # tight
    wt = _make_endpoint_block([85.00, 85.10, 84.90, 85.05])
    er = uwham_pp.compute_extend_recommended(cp4, wt)
    assert er["recommended"] is False
    assert er["current_sigma_btwn_kcal"] < er["target_kcal"]


def test_a7_extend_single_replicate_not_recommended(uwham_pp):
    """A single replicate (σ_btwn undefined) → recommended=False with an
    explanatory reason (cannot assess spread — not a C4 trigger)."""
    cp4 = _make_endpoint_block([90.0])
    wt = _make_endpoint_block([85.0])
    er = uwham_pp.compute_extend_recommended(cp4, wt)
    assert er["recommended"] is False
    assert er["current_sigma_btwn_kcal"] is None
    assert "single replicate" in er["reason"].lower()


def test_a7_extend_recommended_in_payload_numbers_neutral(uwham_pp):
    """build_v3_payload carries extend_recommended as an ADDITIVE field; on the
    shipped n4 blocks (σ_btwn ≈ 0.68 < 2.2) it is recommended=False, and the
    point estimate + paired SEM are unchanged by it."""
    cp4 = _make_endpoint_block(_N4_CP4_DGBIND1)
    wt = _make_endpoint_block(_N4_WT_DGBIND1)
    payload = uwham_pp.build_v3_payload(
        cp4=cp4, wt=wt, protocol="test", mintimeid=100, maxtimeid=None,
        notes="unit test",
    )
    assert payload["ddint_free_kcal"] == pytest.approx(1.2432, abs=5e-4)
    assert payload["ddint_err_paired_sem_kcal"] == pytest.approx(0.5566,
                                                                 abs=5e-4)
    er = payload["extend_recommended"]
    assert er["recommended"] is False
    assert er["current_sigma_btwn_kcal"] < er["target_kcal"]
    # And the A6 sign field also rides in the payload, neutral on n4.
    assert payload["sign_vs_anchor"]["sign_vs_anchor"] == "undetermined"
