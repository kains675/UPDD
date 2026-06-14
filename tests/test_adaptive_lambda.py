#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the Track B adaptive λ-scheduling engine MVP.

Headline acceptance proofs (must pass):
  1. overlap reproduces the manual numbers on the EXISTING densified38
     pilot .out: dplus state 8→9 BC ≈ 0.131 and dminus state 9→10 BC ≈ 0.364
     (tolerance ±0.03).
  2. prescribe reproduces DENSE38v2_LAMBDA_FWD_LINEAR (the manual densified38→v2
     step) from the densified38 pilot overlaps — byte-equal to the registry
     densified38v2 schedule.
  3. constraints guards: endpoint move → raises; soft-core-ladder mutation →
     raises; symmetry desync → raises; count change under "fixed" → raises.
  4. schedule_io round-trip + --free-schedule-file loads a JSON dict and the
     launcher dry-run (cpu) generates cntls with those λ.
  5. estimate_overlap is GPU-free (reads .out only) and the package exposes NO
     auto-launch path (no engine.py iterate-loop in the MVP).

GPU-free throughout (reads only the existing pilot .out text). pymbar optional.

Run with the qmmm env:
    /home/san/miniconda3/envs/qmmm/bin/python -m pytest tests/test_adaptive_lambda.py
"""

import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pytest


_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_PILOT_RUN_DIR = os.path.join(
    _REPO_ROOT, "outputs", "_trackb",
    "production_v2_3_densify38_pilot", "cp4", "free",
)
_PILOT_COMBINED_CNTL = os.path.join(_PILOT_RUN_DIR, "trackb_asyncre.cntl")

# Reference numbers (manual densified38 overlap re-derivation, acceptance #1).
_BC_DPLUS_8_9 = 0.131
_BC_DMINUS_9_10 = 0.364
_BC_TOL = 0.03

# The manual densified38→v2 linear segment (asyncre DENSE38v2_LAMBDA_FWD_LINEAR).
_DENSE38v2_FWD_LINEAR = [0.00, 0.08, 0.16, 0.24, 0.32, 0.38, 0.40, 0.42, 0.44, 0.45]


def _load_module(name: str, path: str):
    """Project test-isolation loader (tests/test_atm_trackb_per_direction.py
    reference) — spec_from_file_location + utils/scripts on sys.path."""
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, os.path.join(_REPO_ROOT, "scripts"))
    sys.path.insert(0, os.path.join(_REPO_ROOT, "utils"))
    sys.path.insert(0, _REPO_ROOT)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def adaptive():
    """Import the adaptive_lambda package submodules under test."""
    sys.path.insert(0, _REPO_ROOT)
    sys.path.insert(0, os.path.join(_REPO_ROOT, "utils"))
    from utils.adaptive_lambda import overlap, judge, prescribe, constraints
    from utils.adaptive_lambda import schedule_io
    return {
        "overlap": overlap,
        "judge": judge,
        "prescribe": prescribe,
        "constraints": constraints,
        "schedule_io": schedule_io,
    }


@pytest.fixture(scope="module")
def asyncre():
    """The asyncre launcher module (densified38 / densified38v2 schedules)."""
    return _load_module(
        "_asyncre_for_test",
        os.path.join(_REPO_ROOT, "scripts", "trackb_production_v2_asyncre.py"),
    )


@pytest.fixture(scope="module")
def seed_schedule(asyncre):
    return asyncre.get_schedule("densified38")


@pytest.fixture(scope="module")
def overlap_report(adaptive, seed_schedule):
    return adaptive["overlap"].estimate_overlap(
        _PILOT_RUN_DIR, jobname="trackb", schedule=seed_schedule,
        use_mbar="never", warmup_cycles=5,
    )


_PILOT_PRESENT = os.path.isdir(_PILOT_RUN_DIR) and any(
    os.path.isfile(os.path.join(_PILOT_RUN_DIR, "dplus", f"r{i}", "trackb_dplus.out"))
    for i in range(19)
)
_needs_pilot = pytest.mark.skipif(
    not _PILOT_PRESENT,
    reason="densified38 pilot .out files not present (archived?)",
)


# ---------------------------------------------------------------------------
# P6 (true soft-core reduced potential) — bit-equivalence fixtures.
# The bit-equivalence reference requires atom_openmm.uwham (atm conda env only).
# In qmmm pytest atom_openmm is absent, so the bit-equiv test SKIPs gracefully
# and the SYNTHETIC test (which monkeypatches the bias helper) carries the
# deterministic wiring proof.
# ---------------------------------------------------------------------------
_REP1_FREE = {
    "cp4": os.path.join(
        _REPO_ROOT, "outputs", "_trackb",
        "production_v2_4_replicate_campaign", "rep1", "cp4", "cp4", "free",
    ),
    "wt": os.path.join(
        _REPO_ROOT, "outputs", "_trackb",
        "production_v2_4_replicate_campaign", "rep1", "wt", "wt", "free",
    ),
}


def _atom_openmm_importable():
    try:
        import atom_openmm.uwham  # noqa: F401
        return True
    except Exception:
        return False


def _rep1_free_present(tag):
    leg = _REP1_FREE[tag]
    cntl = os.path.join(leg, "trackb_asyncre.cntl")
    return os.path.isfile(cntl) and os.path.isfile(
        os.path.join(leg, "r0", "trackb.out")
    )


_needs_bit_equiv = pytest.mark.skipif(
    not (_atom_openmm_importable()
         and _rep1_free_present("cp4") and _rep1_free_present("wt")),
    reason="bit-equivalence needs atom_openmm (atm env) + rep1 free-leg data",
)


# ===========================================================================
# ACCEPTANCE #1 — BC reproduces the verdict numbers
# ===========================================================================
@_needs_pilot
def test_acceptance_1_bc_dplus_8_9(overlap_report):
    pairs = {(p.i, p.j): p for p in overlap_report.per_direction["dplus"]}
    assert (8, 9) in pairs, "dplus 8→9 pair missing from overlap report"
    bc = pairs[(8, 9)].bc
    assert abs(bc - _BC_DPLUS_8_9) <= _BC_TOL, (
        f"dplus 8→9 BC={bc:.4f} not within {_BC_TOL} of verdict "
        f"{_BC_DPLUS_8_9}"
    )


@_needs_pilot
def test_acceptance_1_bc_dminus_9_10(overlap_report):
    pairs = {(p.i, p.j): p for p in overlap_report.per_direction["dminus"]}
    assert (9, 10) in pairs, "dminus 9→10 pair missing from overlap report"
    bc = pairs[(9, 10)].bc
    assert abs(bc - _BC_DMINUS_9_10) <= _BC_TOL, (
        f"dminus 9→10 BC={bc:.4f} not within {_BC_TOL} of verdict "
        f"{_BC_DMINUS_9_10}"
    )


@_needs_pilot
def test_extract_pertE_by_state_groups_by_col0(adaptive):
    by_state = adaptive["overlap"].extract_pertE_by_state(
        _PILOT_RUN_DIR, "trackb", "dplus", warmup_cycles=5,
    )
    # 19-state forward leg; states 0..18 should be present.
    assert set(by_state.keys()) <= set(range(19))
    assert 8 in by_state and 9 in by_state
    # pertE arrays are 1-D float.
    assert by_state[8].ndim == 1 and by_state[8].dtype == float


def test_bhattacharyya_closed_form():
    """The Gaussian BC closed form: identical samples → ~1; well-separated → ~0."""
    bc_mod = pytest.importorskip("utils.adaptive_lambda.overlap")
    rng = np.random.default_rng(0)
    a = rng.normal(0.0, 1.0, 500)
    b = rng.normal(0.0, 1.0, 500)
    near = bc_mod.bhattacharyya_coefficient(a, b)
    assert near > 0.9, f"same-distribution BC should be high, got {near}"
    c = rng.normal(50.0, 1.0, 500)
    far = bc_mod.bhattacharyya_coefficient(a, c)
    assert far < 0.05, f"far-separated BC should be ~0, got {far}"


def test_bhattacharyya_requires_two_samples():
    from utils.adaptive_lambda.overlap import bhattacharyya_coefficient
    with pytest.raises(ValueError):
        bhattacharyya_coefficient(np.array([1.0]), np.array([1.0, 2.0]))


# ===========================================================================
# ACCEPTANCE #2 — prescribe reproduces the densified38→v2 step
# ===========================================================================
@_needs_pilot
def test_acceptance_2_reproduces_dense38v2_linear(
    adaptive, asyncre, seed_schedule, overlap_report,
):
    cfg = adaptive["constraints"].AdaptiveConfig()
    jr = adaptive["judge"].classify(overlap_report, seed_schedule, cfg)
    proposed = adaptive["prescribe"].propose_rebalance(
        overlap_report, jr, seed_schedule, cfg,
    )
    new_linear = proposed["_provenance"]["new_linear_segment"]
    assert new_linear == _DENSE38v2_FWD_LINEAR, (
        f"prescribe did not reproduce DENSE38v2 linear: {new_linear}"
    )


@_needs_pilot
def test_acceptance_2_proposed_equals_registry_v2(
    adaptive, asyncre, seed_schedule, overlap_report,
):
    cfg = adaptive["constraints"].AdaptiveConfig()
    jr = adaptive["judge"].classify(overlap_report, seed_schedule, cfg)
    proposed = adaptive["prescribe"].propose_rebalance(
        overlap_report, jr, seed_schedule, cfg,
    )
    v2 = asyncre.get_schedule("densified38v2")
    for key in ("lambdas_1", "lambdas_2", "w0", "alpha", "u0",
                "intermd", "directions", "lambdas"):
        assert np.allclose(
            [float(x) for x in proposed[key]],
            [float(x) for x in v2[key]],
        ), f"proposed[{key}] != registry densified38v2[{key}]"
    assert proposed["n_states"] == 38


@_needs_pilot
def test_acceptance_2_count_neutral_and_endpoints_pinned(
    adaptive, seed_schedule, overlap_report,
):
    cfg = adaptive["constraints"].AdaptiveConfig()
    jr = adaptive["judge"].classify(overlap_report, seed_schedule, cfg)
    proposed = adaptive["prescribe"].propose_rebalance(
        overlap_report, jr, seed_schedule, cfg,
    )
    # count-neutral
    assert proposed["n_states"] == seed_schedule["n_states"]
    # coupled endpoint pinned
    assert float(proposed["lambdas_1"][0]) == 0.0


@_needs_pilot
def test_judge_finds_bottleneck_and_redundant_plateau(
    adaptive, seed_schedule, overlap_report,
):
    cfg = adaptive["constraints"].AdaptiveConfig()
    jr = adaptive["judge"].classify(overlap_report, seed_schedule, cfg)
    # The single genuine bottleneck is dplus 8→9.
    bn = {(p.i, p.j) for p in jr.bottlenecks.get("dplus", [])}
    assert (8, 9) in bn
    # The plateau (0,1)..(7,8) is redundant.
    rd = {(p.i, p.j) for p in jr.redundant_runs.get("dplus", [])}
    assert (0, 1) in rd and (7, 8) in rd
    # The schedule should NOT pass (it has a bottleneck).
    assert jr.all_pass is False


@_needs_pilot
def test_judge_region_derived_from_w0_intermd_not_hardcoded(
    adaptive, seed_schedule, overlap_report,
):
    """region_of must come from W0/INTERMEDIATE, not a hardcoded 'state 10'."""
    cfg = adaptive["constraints"].AdaptiveConfig()
    jr = adaptive["judge"].classify(overlap_report, seed_schedule, cfg)
    pairs = {(p.i, p.j): p for p in overlap_report.per_direction["dplus"]}
    # (8,9): both W0==0, INTERMEDIATE==0 → linear.
    assert jr.region_of("dplus", pairs[(8, 9)]) == "linear"
    # (10,11): W0!=0 / INTERMEDIATE==1 → softcore.
    if (10, 11) in pairs:
        assert jr.region_of("dplus", pairs[(10, 11)]) == "softcore"


# ===========================================================================
# ACCEPTANCE #3 — constraint guards fail-loud
# ===========================================================================
def test_acceptance_3_endpoint_move_raises(adaptive, asyncre):
    c = adaptive["constraints"]
    seed = asyncre.get_schedule("densified38")
    bad = asyncre.get_schedule("densified38")
    bad["lambdas_1"][0] = 0.01  # move the coupled endpoint λ=0
    with pytest.raises(c.EndpointInvariantError):
        c.endpoint_invariant(seed, bad)


def test_acceptance_3_apex_mutation_raises(adaptive, asyncre):
    c = adaptive["constraints"]
    seed = asyncre.get_schedule("densified38")
    bad = asyncre.get_schedule("densified38")
    # Mutate the apex (W0==1.0) state's ALPHA — apex must be byte-identical.
    apex = [i for i, w in enumerate(bad["w0"]) if float(w) == 1.0][0]
    bad["alpha"][apex] = bad["alpha"][apex] + 0.5
    with pytest.raises(c.EndpointInvariantError):
        c.endpoint_invariant(seed, bad)


def test_acceptance_3_softcore_ladder_mutation_raises(adaptive, asyncre):
    c = adaptive["constraints"]
    seed = asyncre.get_schedule("densified38")
    bad = asyncre.get_schedule("densified38")
    # Mutate a soft-core ladder state (find a non-apex W0!=0 state).
    ladder = [
        i for i in range(len(bad["w0"]))
        if float(bad["w0"][i]) not in (0.0, 1.0)
        and int(round(float(bad["intermd"][i]))) == 1
    ][0]
    bad["u0"][ladder] = float(bad["u0"][ladder]) + 5.0
    with pytest.raises(c.RegionLockedError):
        c.region_locked(seed, bad)


def test_acceptance_3_symmetry_desync_raises(adaptive, asyncre):
    c = adaptive["constraints"]
    bad = asyncre.get_schedule("densified38")
    # Break the whole-tuple-reverse mirror: perturb ONE backward-half lambda.
    n = bad["n_states"]
    fwd = sum(1 for d in bad["directions"] if float(d) >= 0)
    bad["lambdas_1"][fwd] = float(bad["lambdas_1"][fwd]) + 0.07
    with pytest.raises(c.SymmetryMirrorError):
        c.symmetry_mirror(bad)


def test_acceptance_3_count_change_fixed_raises(adaptive, asyncre):
    c = adaptive["constraints"]
    seed = asyncre.get_schedule("densified38")
    bad = asyncre.get_schedule("densified38")
    # Drop one state from every per-state array → count changes.
    for k in ("lambdas", "lambdas_1", "lambdas_2", "directions",
              "intermd", "w0", "alpha", "u0"):
        bad[k] = list(bad[k])[:-1]
    bad["n_states"] = len(bad["directions"])
    with pytest.raises(c.CountPolicyError):
        c.count_policy(seed, bad, policy="fixed")


def test_constraints_pass_on_valid_v2(adaptive, asyncre):
    """The valid densified38→v2 rebalance must clear ALL four guards."""
    c = adaptive["constraints"]
    seed = asyncre.get_schedule("densified38")
    v2 = asyncre.get_schedule("densified38v2")
    assert c.assert_all_constraints(seed, v2) is True


def test_adaptive_config_rejects_bad_thresholds(adaptive):
    c = adaptive["constraints"]
    with pytest.raises(ValueError):
        c.AdaptiveConfig(bottleneck_floor=0.9, redundant_ceiling=0.5)
    with pytest.raises(ValueError):
        c.AdaptiveConfig(count_policy="grow-a-lot")


def test_adaptive_config_owned_values(adaptive):
    """The default config values must be exactly as specified."""
    cfg = adaptive["constraints"].AdaptiveConfig()
    assert cfg.bottleneck_floor == 0.30
    assert cfg.redundant_ceiling == 0.85
    assert cfg.target_margin == 0.40
    assert cfg.max_iterations == 4
    assert cfg.count_policy == "fixed"
    assert cfg.damping_max_moves == 2
    assert cfg.pilot_max_samples == 25
    assert cfg.region == "linear"


# ===========================================================================
# ACCEPTANCE #4 — schedule_io round-trip + --free-schedule-file launcher wiring
# ===========================================================================
def test_acceptance_4_registry_dict_roundtrip(adaptive, asyncre):
    io = adaptive["schedule_io"]
    v2 = asyncre.get_schedule("densified38v2")
    rt = io.to_registry_dict(v2)
    back = io.from_registry_dict(rt)
    for k in ("lambdas_1", "lambdas_2", "w0", "alpha", "u0"):
        assert np.allclose([float(x) for x in rt[k]], [float(x) for x in back[k]])
    assert rt["directions"] == back["directions"]
    assert rt["intermd"] == back["intermd"]
    assert rt["n_states"] == 38


def test_acceptance_4_emit_validated_schedule(adaptive, asyncre, tmp_path):
    io = adaptive["schedule_io"]
    seed = asyncre.get_schedule("densified38")
    v2 = asyncre.get_schedule("densified38v2")
    path = str(tmp_path / "adaptive_schedule.json")
    env = io.emit_validated_schedule(
        v2, path, seed=seed,
        overlap_table=[{"i": 8, "j": 9, "bc": 0.137}],
    )
    assert os.path.isfile(path)
    assert env["provenance"]["regime"] == "ranking_only"
    # endpoints invariant between seed and v2 (the empirical ranking-unbias check).
    assert env["provenance"]["endpoints_invariant"] is True
    assert "proposed_endpoint_hashes" in env["provenance"]
    assert "overlap_table" in env["provenance"]
    # The JSON file is loadable back into a schedule dict.
    loaded = io.load_schedule_dict(path)
    assert loaded["n_states"] == 38


def test_acceptance_4_load_schedule_dict_bare_and_envelope(adaptive, asyncre, tmp_path):
    io = adaptive["schedule_io"]
    v2 = asyncre.get_schedule("densified38v2")
    # bare dict
    bare = str(tmp_path / "bare.json")
    with open(bare, "w") as fh:
        json.dump(io.to_registry_dict(v2), fh)
    assert io.load_schedule_dict(bare)["n_states"] == 38
    # envelope
    env = str(tmp_path / "env.json")
    io.emit_validated_schedule(v2, env)
    assert io.load_schedule_dict(env)["n_states"] == 38


@_needs_pilot
def test_acceptance_4_free_schedule_file_launcher_dryrun(
    adaptive, asyncre, tmp_path,
):
    """--free-schedule-file loads a JSON dict and the launcher dry-run (cpu)
    generates per-direction cntls with those λ (SSOT single write point)."""
    io = adaptive["schedule_io"]
    v2 = asyncre.get_schedule("densified38v2")
    # Build the adaptive_schedule.json.
    sched_json = str(tmp_path / "adaptive_schedule.json")
    io.emit_validated_schedule(v2, sched_json)

    # Stage a COPY of the pilot combined cntl into a temp v21-out-root (do NOT
    # touch the pilot dir). The launcher reads <v21-out-root>/cp4/free/
    # trackb_asyncre.cntl.
    v21_root = tmp_path / "v21root"
    leg_dir = v21_root / "cp4" / "free"
    leg_dir.mkdir(parents=True)
    shutil.copy2(_PILOT_COMBINED_CNTL, str(leg_dir / "trackb_asyncre.cntl"))

    # Apply the schedule file directly via the launcher helper (the same code
    # path main() runs), then verify the per-direction slices inherit the λ.
    pd = _load_module(
        "_pd_for_free_sched_test",
        os.path.join(_REPO_ROOT, "scripts",
                     "trackb_per_direction_production.py"),
    )
    audit = pd.apply_free_schedule_file_to_combined_cntl(
        combined_cntl_path=str(leg_dir / "trackb_asyncre.cntl"),
        schedule_file=sched_json,
    )
    assert audit["n_states"] == 38
    assert set(audit["rewritten_keys"]) >= {
        "LAMBDAS", "DIRECTION", "INTERMEDIATE", "LAMBDA1", "LAMBDA2",
        "ALPHA", "U0", "W0COEFF",
    }

    # The combined cntl now carries the v2 linear segment.
    lambda1_line = None
    with open(leg_dir / "trackb_asyncre.cntl") as fh:
        for line in fh:
            if line.startswith("LAMBDA1"):
                lambda1_line = line
                break
    assert lambda1_line is not None
    for v in (0.08, 0.16, 0.24, 0.32, 0.38, 0.42, 0.44):
        assert str(v) in lambda1_line, (
            f"λ={v} (v2 linear) missing from rewritten combined cntl LAMBDA1"
        )

    # Slicing the rewritten combined cntl produces per-direction cntls with the
    # loaded λ — the downstream consumer inherits the override (no post-hoc
    # per-direction patch).
    gen = pd.generate_per_direction_cntls(
        leg_dir=str(leg_dir), jobname="trackb",
    )
    dplus_cntl = os.path.join(str(leg_dir), "dplus",
                              "trackb_dplus_asyncre.cntl")
    assert os.path.isfile(dplus_cntl)
    dplus_l1 = None
    with open(dplus_cntl) as fh:
        for line in fh:
            if line.startswith("LAMBDA1"):
                dplus_l1 = line
                break
    assert dplus_l1 is not None
    # The dplus slice = forward 19 states; its linear head carries the v2 λ.
    for v in (0.08, 0.16, 0.24, 0.32, 0.38, 0.42, 0.44):
        assert str(v) in dplus_l1


def test_acceptance_4_launcher_dryrun_subprocess(tmp_path):
    """End-to-end: invoke the launcher with --free-schedule-file + --dry-run on
    cpu and confirm exit 0 (no GPU, no launch)."""
    if not _PILOT_PRESENT:
        pytest.skip("pilot .out absent")
    py = sys.executable
    # Build the schedule json via the package.
    sys.path.insert(0, _REPO_ROOT)
    from utils.adaptive_lambda import schedule_io
    asyncre = _load_module(
        "_asyncre_sub",
        os.path.join(_REPO_ROOT, "scripts", "trackb_production_v2_asyncre.py"),
    )
    v2 = asyncre.get_schedule("densified38v2")
    sched_json = str(tmp_path / "adaptive_schedule.json")
    schedule_io.emit_validated_schedule(v2, sched_json)

    v21_root = tmp_path / "v21root"
    leg_dir = v21_root / "cp4" / "free"
    leg_dir.mkdir(parents=True)
    shutil.copy2(_PILOT_COMBINED_CNTL, str(leg_dir / "trackb_asyncre.cntl"))
    out_root = tmp_path / "out"
    out_root.mkdir()

    cmd = [
        py,
        os.path.join(_REPO_ROOT, "scripts",
                     "trackb_per_direction_production.py"),
        "--free-pilot",
        "--legs", "free",
        "--endpoints", "cp4",
        "--gpu-host", "cpu",
        "--dry-run",
        "--free-schedule-file", sched_json,
        "--v21-out-root",
        os.path.relpath(str(v21_root), _REPO_ROOT),
        "--out-root", os.path.relpath(str(out_root), _REPO_ROOT),
    ]
    proc = subprocess.run(
        cmd, cwd=_REPO_ROOT, capture_output=True, text=True, timeout=180,
    )
    assert proc.returncode == 0, (
        f"launcher dry-run exit {proc.returncode}\nSTDOUT:\n{proc.stdout}\n"
        f"STDERR:\n{proc.stderr}"
    )
    assert "combined cntl re-emitted" in proc.stdout, (
        f"free-schedule-file application not reported:\n{proc.stdout}"
    )


# ===========================================================================
# ACCEPTANCE #5 — GPU-free + no auto-launch path in the MVP
# ===========================================================================
def test_acceptance_5_estimate_overlap_is_gpu_free(adaptive):
    """estimate_overlap must NOT import openmm / trigger any GPU path — it reads
    .out text only. We assert openmm is not pulled in by the call."""
    # Snapshot whether openmm was already imported (other tests may have).
    pre_openmm = "openmm" in sys.modules
    if not _PILOT_PRESENT:
        pytest.skip("pilot .out absent")
    sys.path.insert(0, _REPO_ROOT)
    asyncre = _load_module(
        "_asyncre_gpu_free",
        os.path.join(_REPO_ROOT, "scripts", "trackb_production_v2_asyncre.py"),
    )
    seed = asyncre.get_schedule("densified38")
    rep = adaptive["overlap"].estimate_overlap(
        _PILOT_RUN_DIR, jobname="trackb", schedule=seed,
        use_mbar="never", warmup_cycles=5,
    )
    assert rep.per_direction  # produced something
    # estimate_overlap itself must not have introduced openmm.
    if not pre_openmm:
        assert "openmm" not in sys.modules, (
            "estimate_overlap imported openmm — it must be GPU/OpenMM-free"
        )


def test_acceptance_5_no_engine_module_in_mvp():
    """The MVP must NOT ship engine.py (the iterate/auto-launch loop is Stage 2)."""
    engine_path = os.path.join(
        _REPO_ROOT, "utils", "adaptive_lambda", "engine.py"
    )
    assert not os.path.isfile(engine_path), (
        "engine.py present — the auto-iterate/auto-launch loop is OUT OF SCOPE "
        "for the MVP (Stage 2)"
    )


def test_acceptance_5_package_exposes_no_launch_symbol(adaptive):
    """The public package surface must expose diagnose/judge/prescribe/io only —
    no launch / subprocess / iterate symbol."""
    import utils.adaptive_lambda as pkg
    banned = ("launch", "iterate", "run_campaign", "auto", "subprocess",
              "spawn", "engine")
    for name in pkg.__all__:
        low = name.lower()
        assert not any(b in low for b in banned), (
            f"package __all__ exposes a launch-like symbol {name!r}"
        )


def test_estimate_overlap_use_mbar_modes(adaptive):
    """use_mbar='never' yields None matrices; 'require' raises if pymbar absent."""
    if not _PILOT_PRESENT:
        pytest.skip("pilot .out absent")
    asyncre = _load_module(
        "_asyncre_mbar",
        os.path.join(_REPO_ROOT, "scripts", "trackb_production_v2_asyncre.py"),
    )
    seed = asyncre.get_schedule("densified38")
    rep = adaptive["overlap"].estimate_overlap(
        _PILOT_RUN_DIR, jobname="trackb", schedule=seed,
        use_mbar="never", warmup_cycles=5,
    )
    assert all(m is None for m in rep.mbar_matrices.values())


def test_mbar_overlap_matrix_optional(adaptive):
    """mbar_overlap_matrix degrades to None gracefully and NEVER hard-depends on
    pymbar / atom_openmm. P6 changed its contract: it consumes the FULL-column
    sample dict ({state: array[n,7]}) + a per-state soft-core schedule and builds
    the TRUE soft-core operator (no linear-λ surrogate). All graceful paths below
    return None (atom_openmm absent in qmmm), and a malformed OLD-API 1-D dict
    must degrade — not raise."""
    overlap = adaptive["overlap"]
    rng = np.random.default_rng(1)

    # (a) OLD-API shape (bare 1-D pertE dict, no schedule) → graceful None.
    old_shape = {
        0: rng.normal(190.0, 4.0, 25),
        1: rng.normal(188.0, 4.0, 25),
        2: rng.normal(50.0, 60.0, 25),
    }
    assert overlap.mbar_overlap_matrix(old_shape) is None

    # (b) NEW-API full-column samples + schedule. When atom_openmm is importable
    # AND pymbar present it returns a square matrix; otherwise None (graceful).
    def _row(l1, l2, a, u0, w0, potE, pertE):
        return [l1, l2, a, u0, w0, potE, pertE]
    samples = {
        0: np.array([_row(0.0, 0.0, 0.1, 110.0, 0.0, -14288.0 + rng.normal(0, 2),
                          186.0 + rng.normal(0, 4)) for _ in range(25)]),
        1: np.array([_row(0.45, 0.46, 0.12, 105.0, 0.15,
                          -14270.0 + rng.normal(0, 2),
                          95.0 + rng.normal(0, 4)) for _ in range(25)]),
        2: np.array([_row(0.5, 0.5, 0.25, 82.0, 1.0,
                          -14250.0 + rng.normal(0, 2),
                          50.0 + rng.normal(0, 6)) for _ in range(25)]),
    }
    schedule = {
        "lambda1": [0.0, 0.45, 0.5], "lambda2": [0.0, 0.46, 0.5],
        "alpha": [0.1, 0.12, 0.25], "u0": [110.0, 105.0, 82.0],
        "w0": [0.0, 0.15, 1.0],
    }
    mat = overlap.mbar_overlap_matrix(samples, schedule=schedule)
    if mat is not None:
        assert mat.shape == (3, 3)


def test_histogram_intersection_is_advisory_only(adaptive):
    """Histogram intersection is computed but is NEVER consulted by judge.classify
    (banned as gate). Sanity: it returns a value in [0,1]."""
    rng = np.random.default_rng(2)
    a = rng.normal(0.0, 1.0, 100)
    b = rng.normal(0.0, 1.0, 100)
    val = adaptive["overlap"].histogram_intersection_advisory(a, b)
    assert 0.0 <= val <= 1.0


@_needs_pilot
def test_propose_rebalance_refuses_when_no_bottleneck(adaptive, asyncre):
    """If the seed already clears the overlap floor (no bottleneck), propose_
    rebalance REFUSES (raises ValueError) rather than silently emitting the seed
    (verdict refuse-and-escalate)."""
    # densified38v2 already fixed the bottleneck; on its OWN overlaps it should
    # have no linear bottleneck. We approximate by feeding the v2 schedule with a
    # synthetic all-pass report.
    overlap = adaptive["overlap"]
    judge = adaptive["judge"]
    prescribe = adaptive["prescribe"]
    cfg = adaptive["constraints"].AdaptiveConfig()
    v2 = asyncre.get_schedule("densified38v2")
    # Synthetic report: all pairs high BC (no bottleneck).
    rep = overlap.OverlapReport(run_dir="x", jobname="trackb", warmup_cycles=5)
    pairs = []
    n = v2["n_states"]
    for i in range(n - 1):
        pairs.append(overlap.PairOverlap(
            i=i, j=i + 1, bc=0.9, n_i=25, n_j=25,
            mean_i=100.0, mean_j=101.0, sd_i=3.0, sd_j=3.0,
        ))
    rep.per_direction["dplus"] = pairs
    jr = judge.classify(rep, v2, cfg)
    assert jr.all_pass is True
    with pytest.raises(ValueError):
        prescribe.propose_rebalance(rep, jr, v2, cfg)


# ===========================================================================
# P6 — TRUE ATM soft-core reduced potential in overlap.py (P6, UNWIRED)
# ===========================================================================
#
# overlap.mbar_overlap_matrix previously built u_kn from a LINEAR-λ surrogate
# (u = β·λ·pertE), which omits the ATM soft-core bias precisely in the bridge
# region where the bias DEFINES the states — the physics-wrong operator, not an
# approximation (P6 / Klimovich-Shirts-Mobley 2015). P6 replaces it with
# the TRUE soft-core reduced potential reconstructed via build_softcore_neg_pot,
# which REUSES atom_openmm.uwham._bias_fcn / _npot_fcn (the SAME helpers the
# postprocess UWHAM estimator trusts — NOT a re-implementation).
#
# It stays NUMBERS-NEUTRAL because estimate_overlap is UNWIRED from every
# production gate (no production caller imports it — only __init__ + tests).


def _ref_soft_core_helpers():
    """Reference (test-double) ATM soft-core bias / negative-reduced-potential.

    VERBATIM transcription of atom_openmm.uwham._bias_fcn / _npot_fcn (the
    upstream form) used ONLY to drive the SYNTHETIC deterministic test in the
    qmmm env, where atom_openmm is NOT importable. Production code never uses
    this — build_softcore_neg_pot reuses the REAL atom_openmm helpers. The two
    must compute the identical operator (the bit-equivalence test, run in the
    atm env, proves the real path matches the postprocess SSOT)."""
    def _bias_fcn(epert, lam1, lam2, alpha, u0, w0):
        ebias1 = np.zeros_like(epert)
        if alpha > 0:
            ee = 1 + np.exp(-alpha * (epert - u0))
            ebias1 = (lam2 - lam1) * np.log(ee) / alpha
        return ebias1 + lam2 * epert + w0

    def _npot_fcn(e0, epert, bet, lam1, lam2, alpha, u0, w0):
        return -bet * (e0 + _bias_fcn(epert, lam1, lam2, alpha, u0, w0))

    return _bias_fcn, _npot_fcn


@_needs_bit_equiv
def test_p6_bit_equivalence_vs_postprocess_neg_pot(adaptive):
    """BC1 (NON-NEGOTIABLE): the overlap.py soft-core reconstruction is
    bit-equivalent (atol=1e-9) to the postprocess UWHAM SSOT neg_pot on the
    free-leg n=4 cp4 AND wt rep1 reduced potentials.

    Compared after a canonical row-lexsort: postprocess builds rows in
    walker-file (r0..rN) order; build_softcore_neg_pot groups rows by state. The
    row SET (and every per-sample reduced-potential value) must be identical —
    row order is a construction artifact, irrelevant to the overlap matrix
    (row-permutation invariant). A wrong column / β / target-vs-walker indexing
    shifts values by ≫ 1e-9, so atol=1e-9 on the sorted matrices is the
    bug-detector (BC1)."""
    import importlib.util as _ilu
    import pandas as pd
    from atom_openmm import uwham

    overlap = adaptive["overlap"]
    spec = _ilu.spec_from_file_location(
        "_pp_for_p6",
        os.path.join(_REPO_ROOT, "scripts", "trackb_uwham_postprocess.py"),
    )
    pp = _ilu.module_from_spec(spec)
    sys.path.insert(0, os.path.join(_REPO_ROOT, "scripts"))
    spec.loader.exec_module(pp)

    columns = ["stateid", "temperature", "direction", "lambda1", "lambda2",
               "alpha", "u0", "w0", "potE", "pertE", "trash"]
    _F = overlap._SAMPLE_FIELDS  # (l1,l2,alpha,u0,w0,potE,pertE)

    def _canon(M):
        # Stable lexsort over all columns → order-independent comparison.
        return M[np.lexsort(M.T[::-1])]

    for tag in ("cp4", "wt"):
        leg = _REP1_FREE[tag]
        sched = pp._parse_cntl_schedule(
            os.path.join(leg, "trackb_asyncre.cntl"))
        assert sched is not None, f"{tag}: cntl parse failed"
        directions = np.array(sched["directions"])
        nstates = len(directions)
        lambda1 = np.array(sched["lambda1"]); lambda2 = np.array(sched["lambda2"])
        alpha = np.array(sched["alpha"]); u0 = np.array(sched["u0"])
        w0 = np.array(sched["w0"])
        bet = 1.0 / (0.001986209 * 300.0)
        leg1istate = int(np.where(directions >= 0)[0][-1])
        leg2istate = int(np.where(directions < 0)[0][0])

        dfs = []
        for i in range(nstates):
            df = pd.read_csv(os.path.join(leg, f"r{i}", "trackb.out"),
                             sep=r"\s+", header=None, names=columns,
                             index_col=False)
            df["timeid"] = np.arange(1, len(df) + 1)
            dfs.append(df)
        data = pd.concat(dfs, ignore_index=True)
        mask = data["timeid"] >= 100  # postprocess default mintimeid

        # ---- forward leg (states 0..leg1istate) ----
        d1 = data[mask & (data["stateid"] <= leg1istate)]
        ids1 = np.arange(leg1istate + 1)
        e0 = d1["potE"].values.copy()
        for i in range(len(d1)):
            e0[i] -= uwham._bias_fcn(
                d1["pertE"].iloc[i], d1["lambda1"].iloc[i],
                d1["lambda2"].iloc[i], d1["alpha"].iloc[i],
                d1["u0"].iloc[i], d1["w0"].iloc[i])
        ref1 = np.zeros((len(d1), len(ids1)))
        for c, be in enumerate(ids1):
            ref1[:, c] = uwham._npot_fcn(
                e0, d1["pertE"].values, bet,
                lambda1[be], lambda2[be], alpha[be], u0[be], w0[be])
        sbs1 = {}
        for st in ids1:
            sub = d1[d1["stateid"] == st]
            sbs1[int(st)] = np.column_stack(
                [sub[_F[k]].values for k in range(len(_F))])
        mine1 = overlap.build_softcore_neg_pot(
            sbs1, target_states=[int(s) for s in ids1],
            schedule=sched, temperature_K=300.0)

        assert ref1.shape == mine1.shape, f"{tag}/fwd shape mismatch"
        assert np.allclose(_canon(ref1), _canon(mine1), rtol=0, atol=1e-9), (
            f"{tag}/fwd NOT bit-equivalent to postprocess neg_pot "
            f"(maxdiff={np.max(np.abs(_canon(ref1)-_canon(mine1))):.3e})"
        )

        # ---- backward leg (states leg2istate..nstates-1, reversed cols) ----
        d2 = data[mask & (data["stateid"] >= leg2istate)]
        ids2 = np.arange(leg2istate, nstates)[::-1]
        e0b = d2["potE"].values.copy()
        for i in range(len(d2)):
            e0b[i] -= uwham._bias_fcn(
                d2["pertE"].iloc[i], d2["lambda1"].iloc[i],
                d2["lambda2"].iloc[i], d2["alpha"].iloc[i],
                d2["u0"].iloc[i], d2["w0"].iloc[i])
        ref2 = np.zeros((len(d2), len(ids2)))
        for c, be in enumerate(ids2):
            ref2[:, c] = uwham._npot_fcn(
                e0b, d2["pertE"].values, bet,
                lambda1[be], lambda2[be], alpha[be], u0[be], w0[be])
        sbs2 = {}
        for st in np.arange(leg2istate, nstates):
            sub = d2[d2["stateid"] == st]
            sbs2[int(st)] = np.column_stack(
                [sub[_F[k]].values for k in range(len(_F))])
        mine2 = overlap.build_softcore_neg_pot(
            sbs2, target_states=[int(s) for s in ids2],
            schedule=sched, temperature_K=300.0)

        assert ref2.shape == mine2.shape, f"{tag}/bwd shape mismatch"
        assert np.allclose(_canon(ref2), _canon(mine2), rtol=0, atol=1e-9), (
            f"{tag}/bwd NOT bit-equivalent to postprocess neg_pot "
            f"(maxdiff={np.max(np.abs(_canon(ref2)-_canon(mine2))):.3e})"
        )


def test_p6_synthetic_soft_core_reconstruction(adaptive, monkeypatch):
    """SYNTHETIC deterministic proof (qmmm-runnable, no atom_openmm) that
    build_softcore_neg_pot wires the soft-core operator correctly:
      * subtracts the per-sample OWN bias (e0 = potE − bias(own params)),
      * builds each column from the per-TARGET-state SCHEDULE params (not walker
        rows) — BC2,
      * uses β = 1/(_UWHAM_KB_KCAL·T) — BC4.

    Monkeypatches overlap._import_uwham_bias to return a reference (verbatim)
    soft-core implementation, then asserts the reconstruction equals a hand-rolled
    expected matrix to machine precision (atol=1e-12)."""
    overlap = adaptive["overlap"]
    _bias_fcn, _npot_fcn = _ref_soft_core_helpers()
    monkeypatch.setattr(overlap, "_import_uwham_bias",
                        lambda: (_bias_fcn, _npot_fcn))

    # Two states, deterministic samples. Each sample carries its OWN sampled
    # (l1,l2,alpha,u0,w0,potE,pertE) — DELIBERATELY different from the schedule
    # SSOT so a target-vs-walker indexing bug would surface.
    F = overlap._SAMPLE_FIELDS
    # state 0 samples (own params: alpha 0.1, u0 110, w0 0)
    s0 = np.array([
        [0.10, 0.10, 0.10, 110.0, 0.0, -14288.0, 186.0],
        [0.10, 0.10, 0.10, 110.0, 0.0, -14290.0, 150.0],
        [0.10, 0.10, 0.10, 110.0, 0.0, -14285.0, 200.0],
    ])
    # state 1 samples (own params: alpha 0.12, u0 105, w0 0.15 — soft-core)
    s1 = np.array([
        [0.45, 0.46, 0.12, 105.0, 0.15, -14270.0, 95.0],
        [0.45, 0.46, 0.12, 105.0, 0.15, -14275.0, 80.0],
        [0.45, 0.46, 0.12, 105.0, 0.15, -14272.0, 110.0],
    ])
    samples_by_state = {0: s0, 1: s1}

    # Schedule SSOT (per GLOBAL state). Column k uses schedule[k], NOT walker row.
    schedule = {
        "lambda1": [0.0, 0.45],
        "lambda2": [0.0, 0.46],
        "alpha":   [0.10, 0.12],
        "u0":      [110.0, 105.0],
        "w0":      [0.0, 0.15],
    }
    target_states = [0, 1]
    T = 300.0
    beta = 1.0 / (overlap._UWHAM_KB_KCAL * T)

    got = overlap.build_softcore_neg_pot(
        samples_by_state, target_states=target_states,
        schedule=schedule, temperature_K=T)

    # Hand-rolled expected: rows = s0 then s1 (sorted state order), e0 from OWN
    # bias, columns from SCHEDULE params.
    data = np.concatenate([s0, s1], axis=0)
    fi = {n: i for i, n in enumerate(F)}
    pertE = data[:, fi["pertE"]]; potE = data[:, fi["potE"]]
    e0 = potE.copy()
    for i in range(len(data)):
        e0[i] -= _bias_fcn(pertE[i], data[i, fi["lambda1"]], data[i, fi["lambda2"]],
                           data[i, fi["alpha"]], data[i, fi["u0"]], data[i, fi["w0"]])
    expected = np.zeros((len(data), len(target_states)))
    for c, st in enumerate(target_states):
        expected[:, c] = _npot_fcn(
            e0, pertE, beta,
            schedule["lambda1"][st], schedule["lambda2"][st],
            schedule["alpha"][st], schedule["u0"][st], schedule["w0"][st])

    assert got.shape == (6, 2)
    assert np.allclose(got, expected, rtol=0, atol=1e-12), (
        f"synthetic reconstruction mismatch maxdiff="
        f"{np.max(np.abs(got-expected)):.3e}"
    )

    # Negative control: a LINEAR-λ surrogate (the deleted operator) would give a
    # materially different matrix — confirm the soft-core operator is NOT the
    # surrogate (otherwise the fix is a no-op).
    surrogate = np.zeros_like(expected)
    lambdas = [0.0, 1.0]
    all_pert = pertE
    for c in range(2):
        surrogate[:, c] = beta * lambdas[c] * all_pert
    assert not np.allclose(got, surrogate, rtol=0, atol=1e-3), (
        "soft-core reconstruction collapsed to the linear-λ surrogate"
    )


def test_p6_no_surrogate_fallback_when_atom_openmm_absent(adaptive, monkeypatch):
    """The deleted linear-λ surrogate must NOT come back as a fallback. When
    atom_openmm is unavailable, build_softcore_neg_pot RAISES (no surrogate) and
    mbar_overlap_matrix degrades to None (BC stays the gate proxy) — P6:
    a surrogate O is the physics-wrong operator, so emitting NO number is correct
    where the true operator is unavailable."""
    overlap = adaptive["overlap"]
    monkeypatch.setattr(overlap, "_import_uwham_bias", lambda: None)

    sbs = {
        0: np.tile([0.0, 0.0, 0.1, 110.0, 0.0, -14288.0, 186.0], (5, 1)),
        1: np.tile([0.45, 0.46, 0.12, 105.0, 0.15, -14270.0, 95.0], (5, 1)),
    }
    schedule = {
        "lambda1": [0.0, 0.45], "lambda2": [0.0, 0.46],
        "alpha": [0.1, 0.12], "u0": [110.0, 105.0], "w0": [0.0, 0.15],
    }
    with pytest.raises(RuntimeError):
        overlap.build_softcore_neg_pot(sbs, [0, 1], schedule)

    # mbar_overlap_matrix swallows the RuntimeError and returns None (graceful,
    # never a surrogate number).
    mat = overlap.mbar_overlap_matrix(sbs, schedule=schedule)
    assert mat is None


def test_p6_mbar_overlap_matrix_returns_none_without_schedule(adaptive):
    """Without the per-state soft-core SSOT arrays, mbar_overlap_matrix returns
    None (it refuses to emit a surrogate O) — BC remains the gate."""
    overlap = adaptive["overlap"]
    sbs = {
        0: np.tile([0.0, 0.0, 0.1, 110.0, 0.0, -14288.0, 186.0], (5, 1)),
        1: np.tile([0.45, 0.46, 0.12, 105.0, 0.15, -14270.0, 95.0], (5, 1)),
    }
    assert overlap.mbar_overlap_matrix(sbs, schedule=None) is None


def test_p6_estimate_overlap_unwired_from_production_gates():
    """THE dormancy invariant (P6 BC3): estimate_overlap /
    mbar_overlap_matrix / judge.classify are NOT wired into any production gate.
    The free-leg numbers come from trackb_uwham_postprocess.py → atom_openmm.uwham,
    a path overlap.py never touches; the occupancy/crossing production gate
    (trackb_per_direction_production.py) does NOT import these symbols. P6 keeps
    the surrogate-replacement UNWIRED — it gates nothing today.

    This test greps the production launchers for any import/call of the overlap
    diagnostic symbols and asserts there is none (a future Stage-2 wiring would
    require its OWN verdict — BC5)."""
    prod_files = [
        os.path.join(_REPO_ROOT, "scripts",
                     "trackb_per_direction_production.py"),
        os.path.join(_REPO_ROOT, "scripts",
                     "trackb_production_v2_asyncre.py"),
        os.path.join(_REPO_ROOT, "scripts", "trackb_uwham_postprocess.py"),
    ]
    # ``estimate_overlap`` is the unique gate-entry symbol of the adaptive
    # overlap diagnostic — it must appear in NO production launcher. (NOTE: a
    # bare-substring scan of ``mbar_overlap_matrix`` would FALSE-POSITIVE on the
    # postprocess's OWN distinct P11 function ``_compute_mbar_overlap_matrix``,
    # which is unrelated; so we key on ``estimate_overlap`` + the import of the
    # adaptive overlap/judge MODULES, not on shared substrings.)
    for pf in prod_files:
        if not os.path.isfile(pf):
            continue
        with open(pf) as fh:
            text = fh.read()
        assert "estimate_overlap" not in text, (
            f"{os.path.basename(pf)} references estimate_overlap — it must stay "
            f"UNWIRED from production (P6 BC3). Wiring it is a separate "
            f"Stage-2 decision (BC5)."
        )
        # The ONLY adaptive_lambda symbol production may import is schedule_io
        # (the validated-schedule SSOT loader). It must NOT import the overlap or
        # judge diagnostic modules (those carry the MBAR-O gate path).
        for banned_import in (
            "adaptive_lambda.overlap", "adaptive_lambda.judge",
            "adaptive_lambda import overlap", "adaptive_lambda import judge",
            "from .overlap import", "from .judge import",
            "judge.classify", "build_softcore_neg_pot",
        ):
            assert banned_import not in text, (
                f"{os.path.basename(pf)} imports/calls {banned_import!r} — the "
                f"adaptive overlap/judge gate must stay UNWIRED from production "
                f"(P6 BC3). Only schedule_io is permitted."
            )


# ---------------------------------------------------------------------------
# Robust MBAR solver + convergence guard (Track B overlap-harness fix).
#
# The pymbar 4.2 DEFAULT solver STALLS on the huge-dynamic-range two-copy ATS
# reduced potential and returns a NON-PHYSICAL overlap matrix (diagonal > 1,
# row-sum != 1) WITHOUT raising. The harness used to swallow that and emit a
# clean adjacent overlap = 0.0 → phantom "lambda-cliff" misdiagnosis. The fix:
#   (1) construct MBAR with solver_protocol="robust" + maximum_iterations>=20000
#       + relative_tolerance=1e-8 (TypeError fallback for pymbar 3.x);
#   (2) a convergence guard (diag <= 1+tol AND |row-sum - 1| <= tol) rejects the
#       non-physical garbage — overlap.py -> None, postprocess -> raise.
# These tests stub pymbar (no real solve / no atom_openmm needed) so they run in
# the qmmm env.
# ---------------------------------------------------------------------------
class _FakeMBAR:
    """Records the ctor kwargs and returns a caller-controlled overlap matrix.

    Class attributes ``last_kwargs`` / ``overlap_matrix`` are set by each test
    before the call site constructs the instance.
    """

    last_args = None
    last_kwargs = None
    overlap_matrix = None

    def __init__(self, u_kn, N_k, **kwargs):
        type(self).last_args = (u_kn, N_k)
        type(self).last_kwargs = dict(kwargs)

    def compute_overlap(self):
        return {"matrix": type(self).overlap_matrix}


def _install_fake_pymbar(monkeypatch, overlap_matrix):
    """Inject a fake ``pymbar`` module exposing ``_FakeMBAR`` as ``MBAR``."""
    import types as _types

    _FakeMBAR.last_args = None
    _FakeMBAR.last_kwargs = None
    _FakeMBAR.overlap_matrix = np.asarray(overlap_matrix, dtype=float)
    fake = _types.ModuleType("pymbar")
    fake.MBAR = _FakeMBAR
    monkeypatch.setitem(sys.modules, "pymbar", fake)
    return _FakeMBAR


def _overlap_test_inputs():
    """Minimal valid (samples, schedule) for mbar_overlap_matrix's pre-MBAR
    validation (2 states, full-column rows). The neg_pot build is monkeypatched,
    so the actual values only need the right SHAPE."""
    def _row(l1, l2, a, u0, w0, potE, pertE):
        return [l1, l2, a, u0, w0, potE, pertE]
    samples = {
        0: np.array([_row(0.0, 0.0, 0.1, 110.0, 0.0, -14288.0, 186.0)
                     for _ in range(8)]),
        1: np.array([_row(0.5, 0.5, 0.25, 82.0, 1.0, -14250.0, 50.0)
                     for _ in range(8)]),
    }
    schedule = {
        "lambda1": [0.0, 0.5], "lambda2": [0.0, 0.5],
        "alpha": [0.1, 0.25], "u0": [110.0, 82.0], "w0": [0.0, 1.0],
    }
    return samples, schedule


def test_overlap_mbar_passes_robust_solver_kwargs(adaptive, monkeypatch):
    """overlap.mbar_overlap_matrix constructs MBAR with the robust solver kwargs
    (solver_protocol='robust', maximum_iterations>=20000, relative_tolerance)."""
    overlap = adaptive["overlap"]
    samples, schedule = _overlap_test_inputs()
    # neg_pot shape [N_samples, K_states] = [16, 2]; values irrelevant (MBAR fake)
    monkeypatch.setattr(
        overlap, "build_softcore_neg_pot",
        lambda sbs, target_states, schedule, temperature_K: np.zeros((16, 2)),
    )
    # A physical 2x2 overlap matrix so the guard PASSES.
    fake = _install_fake_pymbar(monkeypatch, [[0.7, 0.3], [0.3, 0.7]])

    mat = overlap.mbar_overlap_matrix(samples, schedule=schedule)
    assert mat is not None and mat.shape == (2, 2)
    kw = fake.last_kwargs
    assert kw is not None, "MBAR was never constructed"
    assert kw.get("solver_protocol") == "robust"
    assert kw.get("maximum_iterations", 0) >= 20000
    assert "relative_tolerance" in kw


def test_overlap_guard_rejects_diag_gt_one(adaptive, monkeypatch):
    """A non-physical matrix with a diagonal > 1 (stalled default solver) is
    rejected as None — never a clean 0.0 (phantom lambda-cliff guard)."""
    overlap = adaptive["overlap"]
    samples, schedule = _overlap_test_inputs()
    monkeypatch.setattr(
        overlap, "build_softcore_neg_pot",
        lambda sbs, target_states, schedule, temperature_K: np.zeros((16, 2)),
    )
    # diag 8.81 > 1 (mirrors the observed config-A garbage); off-diag 0.0 → the
    # adjacent overlap would be a phantom 0.0 if the guard did not fire.
    _install_fake_pymbar(monkeypatch, [[8.81, 0.0], [0.0, 8.81]])
    assert overlap.mbar_overlap_matrix(samples, schedule=schedule) is None


def test_overlap_guard_rejects_bad_row_sum(adaptive, monkeypatch):
    """A matrix whose rows do NOT sum to ~1 (non-convergence) is rejected."""
    overlap = adaptive["overlap"]
    samples, schedule = _overlap_test_inputs()
    monkeypatch.setattr(
        overlap, "build_softcore_neg_pot",
        lambda sbs, target_states, schedule, temperature_K: np.zeros((16, 2)),
    )
    # diag <= 1 but row-sum = 9.35 (config-A signature) → guard must fire.
    _install_fake_pymbar(monkeypatch, [[0.5, 8.85], [8.85, 0.5]])
    assert overlap.mbar_overlap_matrix(samples, schedule=schedule) is None


def test_overlap_guard_accepts_physical_matrix(adaptive, monkeypatch):
    """A converged (doubly-stochastic) matrix passes the guard unchanged."""
    overlap = adaptive["overlap"]
    samples, schedule = _overlap_test_inputs()
    monkeypatch.setattr(
        overlap, "build_softcore_neg_pot",
        lambda sbs, target_states, schedule, temperature_K: np.zeros((16, 2)),
    )
    physical = [[0.79, 0.21], [0.065, 0.935]]
    _install_fake_pymbar(monkeypatch, physical)
    mat = overlap.mbar_overlap_matrix(samples, schedule=schedule)
    assert mat is not None
    assert np.allclose(mat, np.asarray(physical), atol=0, rtol=0)


def _load_postprocess_module():
    """Isolation-load trackb_uwham_postprocess (spec_from_file_location)."""
    import importlib.util as _ilu
    spec = _ilu.spec_from_file_location(
        "_pp_overlap_guard",
        os.path.join(_REPO_ROOT, "scripts", "trackb_uwham_postprocess.py"),
    )
    pp = _ilu.module_from_spec(spec)
    sys.path.insert(0, os.path.join(_REPO_ROOT, "scripts"))
    spec.loader.exec_module(pp)
    return pp


def test_postprocess_mbar_passes_robust_solver_kwargs(monkeypatch):
    """_compute_mbar_overlap_matrix (P11 REPORTING) uses the same robust kwargs."""
    pp = _load_postprocess_module()
    neg_pot = np.zeros((4, 2))          # [N_samples=4, K_states=2]
    n_k = [2, 2]
    fake = _install_fake_pymbar(monkeypatch, [[0.7, 0.3], [0.3, 0.7]])
    mat = pp._compute_mbar_overlap_matrix(neg_pot, n_k)
    assert mat.shape == (2, 2)
    kw = fake.last_kwargs
    assert kw is not None
    assert kw.get("solver_protocol") == "robust"
    assert kw.get("maximum_iterations", 0) >= 20000
    assert "relative_tolerance" in kw


def test_postprocess_guard_raises_on_nonphysical_matrix(monkeypatch):
    """The postprocess guard RAISES on a non-physical matrix (caller wraps it in
    graceful-degrade → reported 'unavailable', never a phantom 0.0). diag>1 and
    bad row-sum are both rejected."""
    pp = _load_postprocess_module()
    neg_pot = np.zeros((4, 2))
    n_k = [2, 2]

    _install_fake_pymbar(monkeypatch, [[8.81, 0.0], [0.0, 8.81]])
    with pytest.raises(ValueError):
        pp._compute_mbar_overlap_matrix(neg_pot, n_k)

    _install_fake_pymbar(monkeypatch, [[0.5, 8.85], [8.85, 0.5]])
    with pytest.raises(ValueError):
        pp._compute_mbar_overlap_matrix(neg_pot, n_k)


def test_postprocess_guard_accepts_physical_matrix(monkeypatch):
    """A converged matrix passes the postprocess guard and is returned."""
    pp = _load_postprocess_module()
    neg_pot = np.zeros((4, 2))
    n_k = [2, 2]
    physical = [[0.79, 0.21], [0.065, 0.935]]
    _install_fake_pymbar(monkeypatch, physical)
    mat = pp._compute_mbar_overlap_matrix(neg_pot, n_k)
    assert np.allclose(mat, np.asarray(physical), atol=0, rtol=0)
