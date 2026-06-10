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


def test_adaptive_config_scival_owned_values(adaptive):
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
    """mbar_overlap_matrix returns None gracefully when pymbar is unavailable;
    when available it returns a square matrix. Either outcome is acceptable
    (MVP must not hard-depend on pymbar)."""
    rng = np.random.default_rng(1)
    by_state = {
        0: rng.normal(190.0, 4.0, 25),
        1: rng.normal(188.0, 4.0, 25),
        2: rng.normal(50.0, 60.0, 25),
    }
    mat = adaptive["overlap"].mbar_overlap_matrix(by_state)
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
