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
