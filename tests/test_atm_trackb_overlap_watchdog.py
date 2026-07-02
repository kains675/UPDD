# -*- coding: utf-8 -*-
"""Tests for the Track B per-direction MBAR adjacent-overlap WATCHDOG.

Covers the watchdog's PURE logic surface WITHOUT a real MBAR solve, pymbar, or
atom_openmm (those live only in the atm env; the live end-to-end run is
validated separately on the completed campaign):
  - the three-tier per-pair classifier (HARD_ALERT / THIN_FLAG / OK) incl. the
    CI-contains-0 hard rule and the cry-wolf guard (THIN never drives exit 1);
  - across-rep pooling (skips unavailable reps, max-K-1 alignment);
  - the per-pair bootstrap CI (deterministic with a seeded Generator);
  - the top-level run/overall verdict via a monkeypatched harness (no GPU);
  - stdout RESULT line + exit-code mapping;
  - CLI parsing + the HARKing guard (no densify/rerun symbols in the module).

Uses the importlib spec-loader pattern (no sys.modules pollution).
"""

import importlib.util
import json
import os
import sys

import numpy as np
import pytest


_HERE = os.path.dirname(os.path.abspath(__file__))
_PROJ = os.path.dirname(_HERE)


def _load(name, relpath):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(_PROJ, relpath))
    if spec is None or spec.loader is None:
        pytest.skip("could not locate %s" % relpath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def wd():
    return _load("trackb_overlap_watchdog",
                 "scripts/trackb_overlap_watchdog.py")


# --------------------------- per-pair classifier ---------------------------
def test_classify_ok_above_thin(wd):
    assert wd.classify_pair(0.20, 0.18, 0.03, 0.10) == wd.PAIR_OK


def test_classify_thin_between_floor_and_flag(wd):
    # 0.03 <= pooled < 0.10 with a CI strictly above 0 -> advisory THIN_FLAG.
    assert wd.classify_pair(0.096, 0.094, 0.03, 0.10) == wd.PAIR_THIN_FLAG


def test_classify_hard_below_floor(wd):
    assert wd.classify_pair(0.01, 0.005, 0.03, 0.10) == wd.PAIR_HARD_ALERT


def test_classify_hard_when_ci_contains_zero(wd):
    # Pooled above the thin flag but CI lower bound <= 0 -> HARD (mixing not
    # statistically distinguishable from a broken pair).
    assert wd.classify_pair(0.12, -0.01, 0.03, 0.10) == wd.PAIR_HARD_ALERT
    assert wd.classify_pair(0.12, 0.0, 0.03, 0.10) == wd.PAIR_HARD_ALERT


def test_classify_nonfinite_is_unavailable(wd):
    assert wd.classify_pair(float("nan"), float("nan"), 0.03, 0.10) == \
        wd.PAIR_UNAVAILABLE


# ------------------------------- pooling -----------------------------------
def test_pool_skips_unavailable_reps(wd):
    per_rep = [[0.2, 0.3], None, [0.25, 0.35]]
    n_pairs, by_pair = wd.pool_pairs_across_reps(per_rep)
    assert n_pairs == 2
    assert by_pair[0] == [0.2, 0.25]
    assert by_pair[1] == [0.3, 0.35]


def test_pool_all_unavailable_returns_empty(wd):
    n_pairs, by_pair = wd.pool_pairs_across_reps([None, None])
    assert n_pairs == 0
    assert by_pair == []


def test_pool_skips_none_and_nonfinite_entries(wd):
    per_rep = [[0.2, None], [float("nan"), 0.3]]
    n_pairs, by_pair = wd.pool_pairs_across_reps(per_rep)
    assert n_pairs == 2
    assert by_pair[0] == [0.2]
    assert by_pair[1] == [0.3]


# ----------------------------- bootstrap CI --------------------------------
def test_bootstrap_single_rep_collapses_to_point(wd):
    rng = np.random.default_rng(0)
    pooled, lo, hi = wd.bootstrap_pair_ci([0.2], 2000, rng)
    assert pooled == 0.2 and lo == 0.2 and hi == 0.2


def test_bootstrap_ci_brackets_mean_and_is_deterministic(wd):
    vals = [0.18, 0.20, 0.22]
    p1, lo1, hi1 = wd.bootstrap_pair_ci(vals, 2000, np.random.default_rng(0))
    p2, lo2, hi2 = wd.bootstrap_pair_ci(vals, 2000, np.random.default_rng(0))
    assert p1 == p2 == pytest.approx(0.20)
    # Same seed -> identical CI (reproducible report).
    assert (lo1, hi1) == (lo2, hi2)
    assert lo1 <= 0.20 <= hi1


def test_bootstrap_empty_is_nan(wd):
    rng = np.random.default_rng(0)
    p, lo, hi = wd.bootstrap_pair_ci([], 100, rng)
    assert np.isnan(p) and np.isnan(lo) and np.isnan(hi)


# ------------------- top-level run via fake harness ------------------------
class _FakeOverlapMod:
    """Stand-in for utils.adaptive_lambda.overlap that returns canned matrices
    keyed by (rep_dir, direction) so the watchdog's run logic can be exercised
    with no pymbar / atom_openmm / GPU."""

    def __init__(self, matrices_by_key):
        self._m = matrices_by_key

    def extract_samples_by_state(self, rep_dir, jobname, direction,
                                 warmup_cycles=0):
        # Non-empty sentinel; the fake mbar_overlap_matrix ignores its content.
        key = (rep_dir, direction)
        return {0: np.zeros((1, 7))} if key in self._m else {}

    def mbar_overlap_matrix(self, samples_by_state, schedule=None,
                            temperature_K=300.0):
        # The caller passes the rep/direction via the harness call sequence; we
        # cannot see it here, so this fake is wired by monkeypatching
        # adjacent_overlaps_for_rep instead (see test below). Kept for API shape.
        return None


def test_run_watchdog_healthy_and_alert(wd, monkeypatch, tmp_path):
    # Build a fake layout: 2 endpoints x 1 direction x 2 reps. We monkeypatch
    # the import + the per-rep extractor so no real harness / GPU is needed.
    out_root = str(tmp_path)
    for ep in ("cp4", "wt"):
        for rep in (0, 1):
            os.makedirs(os.path.join(out_root, ep, "bound",
                                     "rep%d" % rep, "dplus"))

    # cp4 healthy (all ~0.2); wt has a broken pair 2 (~0.01) -> HARD_ALERT.
    canned = {
        "cp4": [0.20, 0.21, 0.22],
        "wt": [0.20, 0.21, 0.01],
    }

    def fake_adj(overlap_mod, parse_cntl_schedule, rep_dir, direction,
                 mintimeid, maxtimeid):
        ep = "cp4" if (os.sep + "cp4" + os.sep) in (rep_dir + os.sep) else "wt"
        return list(canned[ep])

    monkeypatch.setattr(wd, "_import_overlap_harness",
                        lambda: (object(), lambda p: {}))
    monkeypatch.setattr(wd, "adjacent_overlaps_for_rep", fake_adj)

    report = wd.run_watchdog(
        out_root=out_root, leg="bound",
        endpoints=["cp4", "wt"], directions=["dplus"],
        mintimeid=100, maxtimeid=None,
        hard_floor=0.03, thin_flag=0.10, n_boot=500, seed=0)

    assert report["overall"] == wd.VERDICT_ALERT  # wt's broken pair drives it
    assert report["n_hard_alert"] == 1
    cp4 = next(g for g in report["groups"] if g["endpoint"] == "cp4")
    wt = next(g for g in report["groups"] if g["endpoint"] == "wt")
    assert cp4["n_hard_alert"] == 0
    assert wt["n_hard_alert"] == 1
    # the broken pair is 2-3 in wt.
    bad = [p for p in wt["per_pair"] if p["verdict"] == wd.PAIR_HARD_ALERT]
    assert len(bad) == 1 and bad[0]["i"] == 2 and bad[0]["j"] == 3


def test_run_watchdog_thin_does_not_alert(wd, monkeypatch, tmp_path):
    out_root = str(tmp_path)
    for rep in (0, 1, 2):
        os.makedirs(os.path.join(out_root, "cp4", "bound",
                                 "rep%d" % rep, "dplus"))

    def fake_adj(*a, **k):
        # one thin pair (0.096) but everything >= hard_floor -> HEALTHY.
        return [0.20, 0.096, 0.22]

    monkeypatch.setattr(wd, "_import_overlap_harness",
                        lambda: (object(), lambda p: {}))
    monkeypatch.setattr(wd, "adjacent_overlaps_for_rep", fake_adj)
    report = wd.run_watchdog(
        out_root=out_root, leg="bound", endpoints=["cp4"],
        directions=["dplus"], mintimeid=100, maxtimeid=None,
        hard_floor=0.03, thin_flag=0.10, n_boot=300, seed=0)
    assert report["overall"] == wd.VERDICT_HEALTHY
    assert report["n_hard_alert"] == 0
    assert report["n_thin"] == 1


def test_run_watchdog_unavailable_reps(wd, monkeypatch, tmp_path):
    out_root = str(tmp_path)
    os.makedirs(os.path.join(out_root, "cp4", "bound", "rep0", "dplus"))

    monkeypatch.setattr(wd, "_import_overlap_harness",
                        lambda: (object(), lambda p: {}))
    monkeypatch.setattr(wd, "adjacent_overlaps_for_rep",
                        lambda *a, **k: None)  # harness produced nothing
    report = wd.run_watchdog(
        out_root=out_root, leg="bound", endpoints=["cp4"],
        directions=["dplus"], mintimeid=100, maxtimeid=None,
        hard_floor=0.03, thin_flag=0.10, n_boot=100, seed=0)
    # No available reps -> no pairs, no hard alerts -> HEALTHY (alert-only must
    # not false-alarm on missing data; that is the orchestrator's completeness
    # check, not the overlap watchdog's).
    assert report["overall"] == wd.VERDICT_HEALTHY
    grp = report["groups"][0]
    assert grp["n_reps_available"] == 0
    assert grp["n_pairs"] == 0


# ----------------------------- stdout / render -----------------------------
def test_render_result_line_healthy(wd):
    report = {
        "leg": "bound", "out_root": "X", "mintimeid": 100, "maxtimeid": None,
        "hard_floor": 0.03, "thin_flag": 0.10, "n_boot": 2000,
        "groups": [], "n_hard_alert": 0, "n_thin": 0, "overall": "HEALTHY",
    }
    txt = wd.render_stdout(report)
    assert "RESULT: HEALTHY" in txt
    assert "RESULT: ALERT" not in txt


def test_render_result_line_alert_has_recommendation(wd):
    report = {
        "leg": "bound", "out_root": "X", "mintimeid": 100, "maxtimeid": None,
        "hard_floor": 0.03, "thin_flag": 0.10, "n_boot": 2000,
        "groups": [{
            "endpoint": "wt", "direction": "dplus",
            "n_reps_available": 3, "n_reps_total": 3,
            "n_hard_alert": 1, "n_thin": 0,
            "per_pair": [{
                "i": 2, "j": 3, "pooled_overlap": 0.01,
                "ci_lo": 0.0, "ci_hi": 0.02, "n_rep": 3,
                "verdict": wd.PAIR_HARD_ALERT,
            }],
        }],
        "n_hard_alert": 1, "n_thin": 0, "overall": "ALERT",
    }
    txt = wd.render_stdout(report)
    assert "RESULT: ALERT" in txt
    assert "lambda-bridge" in txt  # the lambda-spacing-only recommendation


# ------------------------------- CLI ---------------------------------------
def test_parse_csv(wd):
    assert wd._parse_csv("cp4, wt") == ["cp4", "wt"]
    assert wd._parse_csv("dplus") == ["dplus"]
    assert wd._parse_csv(None) == []
    assert wd._parse_csv(" , ") == []


def test_cli_defaults(wd):
    a = wd.build_arg_parser().parse_args(["--out-root", "X"])
    assert a.leg == "bound"
    assert a.endpoints == "cp4,wt"
    assert a.directions == "dplus,dminus"
    assert a.mintimeid == wd.DEFAULT_MINTIMEID
    assert a.hard_floor == wd.DEFAULT_HARD_FLOOR
    assert a.thin_flag == wd.DEFAULT_THIN_FLAG
    assert a.n_boot == wd.DEFAULT_N_BOOT


def test_cli_rejects_unknown_direction(wd):
    rc = wd.main(["--out-root", "X", "--directions", "dbogus"])
    assert rc == 2


def test_main_writes_json_and_maps_exit(wd, monkeypatch, tmp_path):
    out_root = str(tmp_path / "camp")
    os.makedirs(os.path.join(out_root, "cp4", "bound", "rep0", "dplus"))

    monkeypatch.setattr(wd, "_import_overlap_harness",
                        lambda: (object(), lambda p: {}))
    # One broken pair -> ALERT -> exit 1 + JSON written.
    monkeypatch.setattr(wd, "adjacent_overlaps_for_rep",
                        lambda *a, **k: [0.2, 0.005])
    json_path = str(tmp_path / "sub" / "wd.json")
    rc = wd.main(["--out-root", out_root, "--leg", "bound",
                  "--endpoints", "cp4", "--directions", "dplus",
                  "--json", json_path, "--n-boot", "100"])
    assert rc == 1
    assert os.path.isfile(json_path)
    rep = json.load(open(json_path))
    assert rep["overall"] == "ALERT"


# --------------------------- HARKing guard ---------------------------------
def test_module_has_no_rerun_or_densify_calls():
    # The watchdog is ALERT-ONLY: it must NOT import or call any production /
    # densify / rerun entry point, nor shell out, nor push notifications
    # (HARKing guard). Source-level assertion on actual CALL/IMPORT patterns —
    # the words may legitimately appear in prose (the docstring describes what
    # the module deliberately does NOT do), so we forbid invocation syntax only.
    src_path = os.path.join(_PROJ, "scripts", "trackb_overlap_watchdog.py")
    with open(src_path) as fh:
        src = fh.read()
    forbidden_patterns = (
        "run_pool_local(", "run_leg(", "run_one_replicate(",
        "subprocess.run(", "subprocess.Popen(", "subprocess.call(",
        "subprocess.check_output(", "os.system(", "os.popen(",
        "PushNotification(",
        "import subprocess",
        "import trackb_inplace_rbfe_production",  # never import the launcher
        "from trackb_inplace_rbfe_production",
    )
    for pat in forbidden_patterns:
        assert pat not in src, (
            "watchdog must not contain %r (alert-only, no auto-action / no "
            "shell-out / no launcher import)" % pat)
    # No densify entry point invoked (a densify *call*, not the word in prose).
    for densify_call in ("densify(", "_densify", "densify_pilot("):
        assert densify_call not in src, (
            "watchdog must not invoke densify (alert + recommend only)")
