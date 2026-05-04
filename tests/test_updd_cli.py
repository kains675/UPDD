"""tests/test_updd_cli.py — PR-NEW-D Day 1 User CLI tests (T1-T12).

Covers ``scripts/updd_cli.py`` (schema ``updd_cli/0.1``).

Test layout:
    T1:  ``updd check`` end-to-end (mock pipeline)
    T2:  ``updd compare`` 3-sequence ranking
    T3:  ``updd refine`` Cycle 1+ workflow (WT + 2 variants)
    T4:  ``updd discover`` file-based input (Day 1 scope)
    T5:  ``updd report`` markdown rendering from a fixture branched_ddg JSON
    T6:  Subcommand dispatch + argparse validation
    T7:  Common-options propagation (CLIContext)
    T8:  Sequence parser integration via ``--sequence "SLL(MTR)GE"``
    T9:  Branched ΔΔG integration shape on ``--wt`` / ``--variants``
    T10: Error handling — invalid ncAA + missing target file
    T11: Dry-run mode — no execution + status="dry_run"
    T12: Output format consistency between json + markdown
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

# Make ``scripts.updd_cli`` importable just like utils.* — mirror conftest's
# repo-root-on-sys.path pattern (conftest already adds utils/).
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

# scripts/ is not a package by default; load by filespec under the
# stable name ``updd_cli`` so attribute access (CLIContext, main, ...) is
# clean.
import importlib.util as _ilu

_CLI_PATH = _REPO_ROOT / "scripts" / "updd_cli.py"
_spec = _ilu.spec_from_file_location("scripts_updd_cli", _CLI_PATH)
_mod = _ilu.module_from_spec(_spec)  # type: ignore[arg-type]
sys.modules["scripts_updd_cli"] = _mod
assert _spec is not None and _spec.loader is not None
_spec.loader.exec_module(_mod)  # type: ignore[union-attr]
updd_cli = _mod

main = updd_cli.main
SCHEMA_VERSION = updd_cli.SCHEMA_VERSION
CLIContext = updd_cli.CLIContext


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def dummy_target(tmp_path: Path) -> Path:
    """Empty placeholder PDB — exists, content does not matter for Day-1 mock."""
    p = tmp_path / "target.pdb"
    p.write_text("HEADER  dummy fixture\n")
    return p


@pytest.fixture()
def out_dir(tmp_path: Path) -> Path:
    return tmp_path / "cli_out"


def _read_json(p: Path):
    with open(p, "r", encoding="utf-8") as f:
        return json.load(f)


# ---------------------------------------------------------------------------
# T1: `updd check` end-to-end (mock pipeline)
# ---------------------------------------------------------------------------

def test_t1_check_end_to_end(dummy_target: Path, out_dir: Path) -> None:
    rc = main([
        "check",
        "--target", str(dummy_target),
        "--sequence", "ICVVQDWGHHRCT",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "check.json")
    assert payload["schema"] == SCHEMA_VERSION
    assert payload["subcommand"] == "check"
    assert payload["sequence"] == "ICVVQDWGHHRCT"
    assert payload["parsed"]["length"] == 13
    # Mock ΔG must be a finite float, SE positive.
    assert isinstance(payload["dG_estimate"], (int, float))
    assert payload["dG_se"] > 0
    assert payload["binder_verdict"] in {"yes", "no", "undetermined"}
    md_text = (out_dir / "check.md").read_text()
    assert "UPDD `check`" in md_text
    assert "Binder verdict" in md_text


# ---------------------------------------------------------------------------
# T2: `updd compare` 3-sequence ranking
# ---------------------------------------------------------------------------

def test_t2_compare_three_sequences(dummy_target: Path, out_dir: Path) -> None:
    rc = main([
        "compare",
        "--target", str(dummy_target),
        "--sequences", "ICVVQDWGHHRCT,ICVVQDW,ICVVQ(MTR)WGHHRCT",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "compare.json")
    assert payload["subcommand"] == "compare"
    assert len(payload["sequences"]) == 3
    ranked = payload["ranked_entries"]
    assert len(ranked) == 3
    # Ranks must be 1..3, no duplicates.
    assert sorted(e["rank"] for e in ranked) == [1, 2, 3]
    # Mock ΔG ascending: rank 1 should have smallest dG_estimate.
    dGs = [e["pipeline"]["dG_estimate"] for e in ranked]
    assert dGs == sorted(dGs)
    # Pairwise vs top should produce 2 entries (n-1).
    assert len(payload["pairwise_vs_top"]) == 2
    for pw in payload["pairwise_vs_top"]:
        for k in ("ddG", "SE_combined", "z_SE", "tier"):
            assert k in pw


# ---------------------------------------------------------------------------
# T3: `updd refine` Cycle 1+ workflow
# ---------------------------------------------------------------------------

def test_t3_refine_wt_two_variants(dummy_target: Path, out_dir: Path) -> None:
    rc = main([
        "refine",
        "--target", str(dummy_target),
        "--wt", "ICVVQDWGHHRCT",
        "--variants", "ICVVQDW(MTR)HHRCT,ICVVQD(NMA)WGHHRCT",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "refine.json")
    assert payload["subcommand"] == "refine"
    assert payload["wt"]["sequence"] == "ICVVQDWGHHRCT"
    variants = payload["variants"]
    assert len(variants) == 2
    # Each variant must carry ΔΔG, SE, z_SE, tier.
    for v in variants:
        for key in ("ddG", "SE_combined", "z_SE", "tier", "rank"):
            assert key in v
        assert v["tier"] in {"X.A", "X.B", "X.C", "X.D", None}
    # Sign convention propagated for `compute_branched_ddg` parity.
    assert payload["sign_convention"].startswith("ddG = mean(variant)")
    assert payload["compute_branched_ddg_compatible"] is True


# ---------------------------------------------------------------------------
# T4: `updd discover` file-based input (Day 1)
# ---------------------------------------------------------------------------

def test_t4_discover_file_input(dummy_target: Path, out_dir: Path,
                                tmp_path: Path) -> None:
    seq_file = tmp_path / "candidates.txt"
    seq_file.write_text(
        "# header comment\n"
        "ICVVQDWGHHRCT\n"
        "SLL(MTR)GEYLL\n"
        "\n"
        "ABCDE\n"  # invalid → parse_error
    )
    rc = main([
        "discover",
        "--target", str(dummy_target),
        "--sequence-file", str(seq_file),
        "--output", str(out_dir),
        "--ncaa-budget", "5",
    ])
    assert rc == 0
    payload = _read_json(out_dir / "discover.json")
    assert payload["subcommand"] == "discover"
    assert payload["n_input"] == 3
    assert payload["ncaa_budget"] == 5
    cands = payload["candidates"]
    assert len(cands) == 3
    # At least one mock-status candidate and one parse_error.
    statuses = [c.get("status") for c in cands]
    assert "mock" in statuses
    assert "parse_error" in statuses


# ---------------------------------------------------------------------------
# T5: `updd report` markdown from branched_ddg fixture
# ---------------------------------------------------------------------------

def _fixture_branched_ddg() -> dict:
    """Minimal branched_ddg/0.1 result fixture."""
    return {
        "schema": "branched_ddg/0.1",
        "tool": "utils.branched_ddg",
        "generated_at": "2026-04-28T12:00:00Z",
        "target_id": "1EBP",
        "solvent_model": "pbsa",
        "df_strategy": "conservative",
        "branches": {
            "wt": {
                "label": "wt", "branch_dir": "/tmp/wt_fixture",
                "n_seed": 3, "degenerate": False,
                "mean_of_seeds": -50.123, "sigma_btwn": 0.5,
                "sigma_w_median": 1.2, "SE": 0.3,
                "solvent_model": "pbsa", "charge_audit": "PASS",
            },
            "variant": {
                "label": "variant", "branch_dir": "/tmp/var_fixture",
                "n_seed": 3, "degenerate": False,
                "mean_of_seeds": -52.456, "sigma_btwn": 0.6,
                "sigma_w_median": 1.3, "SE": 0.35,
                "solvent_model": "pbsa", "charge_audit": "PASS",
            },
        },
        "pairwise": {
            "ddG": -2.333, "SE": 0.46, "df": 4, "t_crit": 2.776,
            "CI95": [-3.61, -1.06], "z_SE": -5.07, "tier": "X.A",
            "sign_convention": "ddG = mean(variant) − mean(wt)", "note": None,
        },
        "charge_audit": "PASS",
    }


def test_t5_report_markdown(tmp_path: Path) -> None:
    fixture_path = tmp_path / "fix_branched_ddg.json"
    fixture_path.write_text(json.dumps(_fixture_branched_ddg()))
    out = tmp_path / "rpt"
    rc = main([
        "report",
        "--result", str(fixture_path),
        "--output", str(out),
        "--report-format", "markdown",
    ])
    assert rc == 0
    md = (out / "report.md").read_text()
    assert "Branched ΔΔG Summary" in md
    assert "tier" in md
    # Source ΔΔG should be rendered verbatim with sign.
    assert "-2.333" in md or "−2.333" in md or "-2.333" in md  # signed minus


# ---------------------------------------------------------------------------
# T6: Subcommand dispatch + argparse validation
# ---------------------------------------------------------------------------

def test_t6_argparse_missing_subcommand() -> None:
    with pytest.raises(SystemExit) as ei:
        main([])
    assert ei.value.code != 0


def test_t6_argparse_invalid_subcommand() -> None:
    with pytest.raises(SystemExit):
        main(["foobar"])


def test_t6_check_missing_target() -> None:
    with pytest.raises(SystemExit):
        main(["check", "--sequence", "ICVVQDWGHHRCT"])


def test_t6_check_missing_sequence(tmp_path: Path) -> None:
    target = tmp_path / "t.pdb"
    target.write_text("HEADER\n")
    with pytest.raises(SystemExit):
        main(["check", "--target", str(target)])


# ---------------------------------------------------------------------------
# T7: Common-options propagation via CLIContext
# ---------------------------------------------------------------------------

def test_t7_cli_context_propagates_options(dummy_target: Path,
                                           out_dir: Path) -> None:
    import argparse
    ns = argparse.Namespace(
        target=str(dummy_target),
        output=str(out_dir),
        solvent="pbsa",
        df_strategy="min",
        strict_charge_audit=True,
        verbose=True,
        quiet=False,
        dry_run=True,
        format="json",
    )
    ctx = CLIContext.from_args(ns)
    assert ctx.solvent == "pbsa"
    assert ctx.df_strategy == "min"
    assert ctx.strict_charge_audit is True
    assert ctx.verbose is True
    assert ctx.dry_run is True
    assert ctx.fmt == "json"
    assert str(ctx.target) == str(dummy_target)
    assert ctx.output_dir == out_dir
    # Logger setup must not raise; level reflects verbose=True (DEBUG).
    import logging as _logging
    log = ctx.setup_logging()
    assert log.level == _logging.DEBUG


def test_t7_cli_context_default_output(tmp_path: Path,
                                       monkeypatch: pytest.MonkeyPatch) -> None:
    import argparse
    monkeypatch.chdir(tmp_path)
    ns = argparse.Namespace(target=None, output=None,
                            solvent="auto", df_strategy="conservative",
                            strict_charge_audit=False, verbose=False,
                            quiet=False, dry_run=False, format="both")
    ctx = CLIContext.from_args(ns)
    # Default output_dir starts with the documented prefix.
    assert str(ctx.output_dir).startswith(updd_cli.DEFAULT_OUTPUT_PREFIX)


# ---------------------------------------------------------------------------
# T8: Sequence parser integration — parens ncAA "SLL(MTR)GE"
# ---------------------------------------------------------------------------

def test_t8_parser_integration_parens_ncaa(dummy_target: Path,
                                           out_dir: Path) -> None:
    rc = main([
        "check",
        "--target", str(dummy_target),
        "--sequence", "SLL(MTR)GE",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "check.json")
    parsed = payload["parsed"]
    assert parsed["format_detected"] == "parenthesis"
    assert parsed["length"] == 6
    ncaa = parsed["ncaa_positions"]
    assert len(ncaa) == 1
    assert ncaa[0]["code"] == "MTR"
    assert ncaa[0]["position"] == 4


# ---------------------------------------------------------------------------
# T9: Branched ΔΔG integration — refine output is compute_branched_ddg shape
# ---------------------------------------------------------------------------

def test_t9_refine_branched_ddg_shape(dummy_target: Path,
                                      out_dir: Path) -> None:
    rc = main([
        "refine",
        "--target", str(dummy_target),
        "--wt", "ICVVQDWGHHRCT",
        "--variants", "ICVVQ(MTR)WGHHRCT",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "refine.json")
    # The schema must mark this payload as compute_branched_ddg-compatible
    # so downstream consumers (or PR-NEW-E HTML report) can ingest it.
    assert payload["compute_branched_ddg_compatible"] is True
    assert "sign_convention" in payload
    v = payload["variants"][0]
    # The tier vocabulary matches branched_ddg.tier_classify exactly.
    assert v["tier"] in {"X.A", "X.B", "X.C", "X.D"}
    # ΔΔG / SE / CI95 fields are mock but typed.
    assert isinstance(v["ddG"], (int, float))
    assert isinstance(v["SE_combined"], (int, float))
    assert v["SE_combined"] > 0
    assert isinstance(v["CI95"], list)
    assert len(v["CI95"]) == 2


# ---------------------------------------------------------------------------
# T10: Error handling
# ---------------------------------------------------------------------------

def test_t10_invalid_ncaa_returns_error(dummy_target: Path,
                                        out_dir: Path,
                                        capsys: pytest.CaptureFixture) -> None:
    rc = main([
        "check",
        "--target", str(dummy_target),
        "--sequence", "SLL(NOTAREALCODEXYZ)GE",
        "--output", str(out_dir),
    ])
    # NcAAValidationError → exit code 4
    assert rc == 4
    captured = capsys.readouterr()
    assert "NcAAValidationError" in captured.err


def test_t10_missing_target_file(tmp_path: Path, out_dir: Path) -> None:
    rc = main([
        "check",
        "--target", str(tmp_path / "does_not_exist.pdb"),
        "--sequence", "ICVVQDWGHHRCT",
        "--output", str(out_dir),
    ])
    # FileNotFoundError → exit code 5
    assert rc == 5


# ---------------------------------------------------------------------------
# T11: Dry-run mode
# ---------------------------------------------------------------------------

def test_t11_dry_run_marks_status(dummy_target: Path, out_dir: Path) -> None:
    rc = main([
        "check",
        "--target", str(dummy_target),
        "--sequence", "ICVVQDWGHHRCT",
        "--dry-run",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "check.json")
    assert payload["pipeline"]["status"] == "dry_run"
    assert payload["pipeline"]["would_run"] is True
    assert payload["dG_estimate"] is None
    # Verdict must be 'undetermined' when no ΔG was produced.
    assert payload["binder_verdict"] == "undetermined"


def test_t11_dry_run_allows_missing_target(tmp_path: Path,
                                           out_dir: Path) -> None:
    """--dry-run lets the target path be missing (warning only)."""
    rc = main([
        "check",
        "--target", str(tmp_path / "phantom.pdb"),
        "--sequence", "ICVVQDWGHHRCT",
        "--dry-run",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "check.json")
    assert payload["pipeline"]["status"] == "dry_run"


# ---------------------------------------------------------------------------
# T12: Output format consistency between json and markdown
# ---------------------------------------------------------------------------

def test_t12_format_json_only(dummy_target: Path, out_dir: Path) -> None:
    rc = main([
        "compare",
        "--target", str(dummy_target),
        "--sequences", "ICVVQDWGHHRCT,SLLGEYLL",
        "--output", str(out_dir),
        "--format", "json",
    ])
    assert rc == 0
    assert (out_dir / "compare.json").exists()
    assert not (out_dir / "compare.md").exists()


def test_t12_format_markdown_only(dummy_target: Path, out_dir: Path) -> None:
    rc = main([
        "compare",
        "--target", str(dummy_target),
        "--sequences", "ICVVQDWGHHRCT,SLLGEYLL",
        "--output", str(out_dir),
        "--format", "markdown",
    ])
    assert rc == 0
    assert (out_dir / "compare.md").exists()
    assert not (out_dir / "compare.json").exists()


def test_t12_format_both_semantic_equivalence(dummy_target: Path,
                                              out_dir: Path) -> None:
    """JSON + markdown carry equivalent ΔG values for the same run."""
    rc = main([
        "compare",
        "--target", str(dummy_target),
        "--sequences", "ICVVQDWGHHRCT,ICVVQ(MTR)WGHHRCT",
        "--output", str(out_dir),
        "--format", "both",
    ])
    assert rc == 0
    payload = _read_json(out_dir / "compare.json")
    md = (out_dir / "compare.md").read_text()
    # Every JSON-reported ΔG (rounded to 2 decimals) shows up verbatim in the
    # markdown table.
    for e in payload["ranked_entries"]:
        dG = e["pipeline"]["dG_estimate"]
        if isinstance(dG, (int, float)):
            assert f"{dG:+.2f}" in md


# ===========================================================================
# Day 2 — Real Stage 4 hookup tests (T13-T15)
# ===========================================================================
#
# These tests exercise the read-existing path (Mode A) of PR-NEW-D Day 2.
# They depend on Phase α post-patch fixtures already present in the repo
# under outputs/1EBP_{WT,MTR13}_calib_s*/mmpbsa_results_fix4_postpatch/. If
# those dirs are absent (fresh checkout / scrubbed working copy), the
# tests skip with an informative message rather than failing.

_REPO = Path(__file__).resolve().parents[1]
_WT_GLOB = (
    "outputs/1EBP_WT_calib_s*/mmpbsa_results_fix4_postpatch"
)
_MTR13_GLOB = (
    "outputs/1EBP_MTR13_calib_s*/mmpbsa_results_fix4_postpatch"
)
_2QKI_GLOB = (
    "outputs/2QKI_Cp4_calib_s*/mmpbsa_results_fix4_postpatch"
)


def _have_phase_alpha_fixtures() -> bool:
    """True iff at least one 1EBP WT and one 1EBP MTR13 mmpbsa_summary.json
    exists in the repo. Used to guard real-data tests on fresh checkouts."""
    import glob as _g
    wt_hits = _g.glob(str(_REPO / _WT_GLOB / "mmpbsa_summary.json"))
    var_hits = _g.glob(str(_REPO / _MTR13_GLOB / "mmpbsa_summary.json"))
    return bool(wt_hits and var_hits)


# ---------------------------------------------------------------------------
# T13 — End-to-end 1EBP_MTR13 regression (Phase α anchor reproduction)
# ---------------------------------------------------------------------------
#
# Reproduces the Phase α track_b post-patch numbers via the CLI layer:
#
#     ΔΔG ≈ +3.37 kcal/mol
#     z_SE ≈ +2.09
#     tier == "X.B"
#
# Two sub-tests: (a) direct --wt-results-dir / --variant-results-dir paths
# and (b) --wt-system-name / --system-name auto-resolution. Both must
# produce identical numeric results.

@pytest.mark.skipif(not _have_phase_alpha_fixtures(),
                    reason="Phase α 1EBP fixtures not present in this checkout")
def test_t13_1ebp_mtr13_regression_direct_paths(out_dir: Path) -> None:
    rc = main([
        "refine",
        "--target", "1EBP",
        "--wt", "ICVVQDWGHHRCT",
        "--variants", "ICVVQ(MTR)WGHHRCT",
        "--wt-results-dir", str(_REPO / _WT_GLOB),
        "--variant-results-dir", str(_REPO / _MTR13_GLOB),
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "refine.json")
    assert payload["mode"] == "real"
    assert payload["compute_branched_ddg_compatible"] is True
    v = payload["variants"][0]
    assert v["mode"] == "real"
    # Phase α anchor reproduction (PR-NEW-B compute_branched_ddg).
    assert abs(v["ddG"] - 3.373) < 0.01, \
        f"ΔΔG drift: got {v['ddG']:.4f}, expected ~+3.37"
    assert abs(v["z_SE"] - 2.088) < 0.05, \
        f"z_SE drift: got {v['z_SE']:.4f}, expected ~+2.09"
    assert v["tier"] == "X.B", \
        f"tier mismatch: got {v['tier']!r}, expected 'X.B'"
    # CI95 must include 0 (X.B by definition: sign-significant only).
    ci_lo, ci_hi = v["CI95"]
    assert ci_lo < 0.0 < ci_hi, \
        f"X.B implies CI95 includes 0; got [{ci_lo:.3f}, {ci_hi:.3f}]"
    # The full PR-NEW-B branched_ddg result must be embedded for downstream
    # consumers (PR-NEW-E HTML report etc).
    assert "branched_ddg" in v
    # Engine schema bumped per release:
    #   v0.1 → v0.2 (Light D Hybrid, 2026-04-29)
    #   v0.2 → v0.3 (PR-21 Auto-Diagnostic Level 1, 2026-05-04)
    # Anchor on the live SCHEMA_VERSION import so this stays in lockstep
    # with branched_ddg.py without a manual update on every bump.
    from utils.branched_ddg import SCHEMA_VERSION as _BR_DDG_SCHEMA  # noqa: E402
    assert v["branched_ddg"]["schema"] == _BR_DDG_SCHEMA
    assert v["branched_ddg"]["pairwise"]["tier"] == "X.B"
    # No mock contamination.
    assert v["pipeline"].get("mock") is False
    assert payload["wt"]["mode"] == "real"


@pytest.mark.skipif(not _have_phase_alpha_fixtures(),
                    reason="Phase α 1EBP fixtures not present in this checkout")
def test_t13_1ebp_mtr13_regression_system_name(out_dir: Path) -> None:
    rc = main([
        "refine",
        "--target", "1EBP",
        "--wt", "ICVVQDWGHHRCT",
        "--variants", "ICVVQ(MTR)WGHHRCT",
        "--wt-system-name", "1EBP_WT",
        "--system-name", "1EBP_MTR13",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "refine.json")
    assert payload["mode"] == "real"
    v = payload["variants"][0]
    # Same numeric reproduction whether we go via direct paths or
    # system-name resolution.
    assert abs(v["ddG"] - 3.373) < 0.01
    assert abs(v["z_SE"] - 2.088) < 0.05
    assert v["tier"] == "X.B"


# ---------------------------------------------------------------------------
# T14 — Real Stage 4 hookup smoke test (mode/data_source provenance)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not _have_phase_alpha_fixtures(),
                    reason="Phase α 1EBP fixtures not present in this checkout")
def test_t14_real_mode_provenance_check(out_dir: Path) -> None:
    """`updd check` against an existing single branch reports mode='real',
    data_source path, and no 'mock': True flag."""
    rc = main([
        "check",
        "--target", "1EBP",
        "--sequence", "ICVVQ(MTR)WGHHRCT",
        "--system-name", "1EBP_MTR13",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "check.json")
    assert payload["mode"] == "real"
    # data_source must be a path string, not "mock".
    assert isinstance(payload["data_source"], str)
    assert payload["data_source"] != "mock"
    assert "1EBP_MTR13" in payload["data_source"]
    # The pipeline payload must NOT carry a true mock flag in real mode.
    assert payload["pipeline"].get("mock") is False
    assert payload["pipeline"].get("status") == "real"
    # Real ΔG should be in a physically meaningful range for 1EBP_MTR13
    # (post-patch mean_of_seeds ≈ -20.5 kcal/mol; loose envelope here).
    assert -50.0 < payload["dG_estimate"] < 0.0


def test_t14_mock_flag_forces_mock_path(dummy_target: Path,
                                        out_dir: Path) -> None:
    """`--mock` is a hard override: even if results-dir is supplied (or
    would resolve), mock pipeline runs and 'mock': true is preserved."""
    # Use a real-looking system-name argument with --mock to force mock path.
    rc = main([
        "check",
        "--target", str(dummy_target),
        "--sequence", "ICVVQDWGHHRCT",
        "--system-name", "1EBP_MTR13",  # would resolve, but --mock overrides
        "--mock",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "check.json")
    assert payload["mode"] == "mock"
    assert payload["data_source"] == "mock"
    assert payload["pipeline"]["mock"] is True


# ---------------------------------------------------------------------------
# T15 — Error handling for Mode A
# ---------------------------------------------------------------------------

def test_t15_missing_results_dir_falls_back(dummy_target: Path,
                                            out_dir: Path) -> None:
    """Non-existent --system-name falls back to mock with a logged warning;
    exit code remains 0 (graceful degradation when --mock not explicit)."""
    rc = main([
        "check",
        "--target", str(dummy_target),
        "--sequence", "ICVVQDWGHHRCT",
        "--system-name", "DOES_NOT_EXIST_X9Z",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "check.json")
    # Resolution failed → mock fallback (system-name resolution is best-effort
    # when --execute-pipeline is not set; this matches the documented Day-2
    # behaviour for read-existing missing data).
    assert payload["mode"] == "mock"


def test_t15_corrupt_summary_json_raises(tmp_path: Path,
                                         out_dir: Path) -> None:
    """A results dir with a corrupt mmpbsa_summary.json (empty JSON object)
    raises a BranchedDDGError -> exit code 2 (uncaught generic error)."""
    bogus = tmp_path / "fake_results"
    bogus.mkdir()
    # Empty JSON dict — no 'results' key, no per-snap ΔGs.
    (bogus / "mmpbsa_summary.json").write_text("{}\n")
    rc = main([
        "check",
        "--target", "1EBP",
        "--sequence", "ICVVQDWGHHRCT",
        "--results-dir", str(bogus),
        "--output", str(out_dir),
    ])
    # Generic exception path → exit code 2 per Day-1 main() error chain.
    assert rc == 2


def test_t15_ambiguous_system_name_picks_postpatch(out_dir: Path) -> None:
    """When multiple mmpbsa_results_* variants exist for one system
    (e.g. mmpbsa_results_fix4 + mmpbsa_results_fix4_postpatch), the
    resolver prefers the most-postpatched suffix.

    Behaviour definition: Day-2 _resolve_results_dir picks the first
    matching suffix in priority order
    [fix*_postpatch, *_postpatch, fix*, mmpbsa_results].
    """
    if not _have_phase_alpha_fixtures():
        pytest.skip("Phase α 1EBP fixtures absent")
    rc = main([
        "check",
        "--target", "1EBP",
        "--sequence", "ICVVQDWGHHRCT",
        "--system-name", "1EBP_WT",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "check.json")
    # Resolved data_source must include "_postpatch" — the highest-priority
    # variant. The bare mmpbsa_results / fix4 variants must not be picked.
    ds = payload["data_source"]
    assert "_postpatch" in ds, \
        f"Expected postpatch variant in data_source; got {ds!r}"


# ===========================================================================
# Day 3 — Stage 1-2 wrappers + de-novo discovery tests (T16-T22)
# ===========================================================================
#
# Coverage:
#     T16  RFdiffusion wrapper invocation (subprocess mocked + mock fallback)
#     T17  ProteinMPNN wrapper invocation (subprocess mocked + mock fallback)
#     T18  AF2 wrapper invocation (subprocess mocked + mock fallback)
#     T19  Full Stage 1-2 chain end-to-end via `discover --de-novo`
#     T20  --check semantics for each wrapper (CLI invocation)
#     T21  Failure handling — tool not found with allow_mock=False
#     T22  Mock-mode fallback — no tools installed; entire chain mocks

from unittest.mock import patch, MagicMock

# Wrapper modules — load once, reuse across the Day 3 tests.
from utils.rfdiffusion_wrapper import (
    RFdiffusionWrapper, RFdiffusionResult,
    RFdiffusionNotFoundError, RFdiffusionInvocationError,
)
from utils.proteinmpnn_wrapper import (
    ProteinMPNNWrapper, ProteinMPNNResult,
    ProteinMPNNNotFoundError,
)
from utils.af2_wrapper import (
    AF2Wrapper, AF2Result, AF2NotFoundError,
)


# ---------------------------------------------------------------------------
# T16 — RFdiffusion wrapper invocation
# ---------------------------------------------------------------------------

def test_t16_rfdiffusion_mock_fallback(tmp_path: Path) -> None:
    """When no binary is found, run_diffusion synthesises mock backbones."""
    # Patch _resolve_binary to return None (forces mock fallback path).
    w = RFdiffusionWrapper(allow_mock=True)
    with patch.object(w, "_resolve_binary", return_value=None):
        target = tmp_path / "target.pdb"
        target.write_text("HEADER mock\n")
        res = w.run_diffusion(
            target_pdb=str(target), contigs="[A1-50/0 30-50]",
            num_designs=3, output_dir=str(tmp_path / "rf_out"),
        )
    assert res.mock is True
    assert len(res.backbones) == 3
    for bb in res.backbones:
        assert Path(bb).exists()
    assert res.binary is None
    assert res.schema == "stage12_pipeline/0.1"


def test_t16_rfdiffusion_subprocess_mocked(tmp_path: Path) -> None:
    """A real binary path triggers subprocess.run; the wrapper parses the
    designs/ output directory afterwards."""
    target = tmp_path / "target.pdb"
    target.write_text("HEADER mock\n")
    out = tmp_path / "rf_real_out"
    out.mkdir()

    # Pre-create the expected output PDBs the wrapper will scan after the
    # subprocess "succeeds".
    (out / "design_0.pdb").write_text("ATOM      1  CA  GLY A   1\n")
    (out / "design_1.pdb").write_text("ATOM      1  CA  GLY A   1\n")

    fake_cp = MagicMock()
    fake_cp.returncode = 0
    fake_cp.stdout = "ok\n"
    fake_cp.stderr = ""

    w = RFdiffusionWrapper(binary="/fake/rfdiffusion")
    with patch.object(w, "_resolve_binary",
                      return_value="/fake/rfdiffusion"), \
         patch("utils.rfdiffusion_wrapper.subprocess.run",
               return_value=fake_cp) as mock_run:
        res = w.run_diffusion(
            target_pdb=str(target), contigs="[A1-50/0 30-50]",
            num_designs=2, output_dir=str(out),
        )
    assert mock_run.called
    assert res.mock is False
    assert res.binary == "/fake/rfdiffusion"
    assert len(res.backbones) >= 2


def test_t16_rfdiffusion_invalid_args() -> None:
    """num_designs <= 0 raises ValueError."""
    w = RFdiffusionWrapper(allow_mock=True)
    with pytest.raises(ValueError):
        w.run_diffusion(target_pdb="/tmp/x.pdb",
                        contigs="[A]", num_designs=0)


# ---------------------------------------------------------------------------
# T17 — ProteinMPNN wrapper invocation
# ---------------------------------------------------------------------------

def test_t17_proteinmpnn_mock_fallback(tmp_path: Path) -> None:
    """No binary present → mock sequences synthesised matching CA count."""
    bb = tmp_path / "design.pdb"
    # 5 CA atoms → mock sequences must be length 5.
    pdb_lines = [
        "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 30.00           C",
        "ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00 30.00           C",
        "ATOM      3  CA  GLY A   3       7.600   0.000   0.000  1.00 30.00           C",
        "ATOM      4  CA  GLY A   4      11.400   0.000   0.000  1.00 30.00           C",
        "ATOM      5  CA  GLY A   5      15.200   0.000   0.000  1.00 30.00           C",
    ]
    bb.write_text("\n".join(pdb_lines) + "\n")

    w = ProteinMPNNWrapper(allow_mock=True)
    with patch.object(w, "_resolve_binary", return_value=None):
        res = w.run_mpnn(backbone_pdb=str(bb), num_sequences=4)

    assert res.mock is True
    assert len(res.sequences) == 4
    for s in res.sequences:
        assert isinstance(s["seq"], str)
        # Every sequence is alphabetic and L=5.
        assert s["seq"].isalpha()
        assert len(s["seq"]) == 5
        assert isinstance(s["score"], (int, float))
    assert res.backbone_length == 5
    assert res.schema == "stage12_pipeline/0.1"


def test_t17_proteinmpnn_subprocess_mocked(tmp_path: Path) -> None:
    """Subprocess invocation parses upstream FASTA output."""
    bb = tmp_path / "design.pdb"
    bb.write_text(
        "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 30.00           C\n"
    )
    out = tmp_path / "mpnn_out"
    (out / "seqs").mkdir(parents=True)
    fa = out / "seqs" / "design.fa"
    fa.write_text(
        ">design\nMKTGLAEDKK\n"
        ">T=0.1, sample=1, score=0.8124, seq_recovery=0.42\n"
        "MKTAAAEDKL\n"
        ">T=0.1, sample=2, score=0.7351, seq_recovery=0.38\n"
        "MKVALAEDKR\n"
    )

    fake_cp = MagicMock()
    fake_cp.returncode = 0
    fake_cp.stdout = ""
    fake_cp.stderr = ""

    w = ProteinMPNNWrapper(binary="/fake/mpnn")
    with patch.object(w, "_resolve_binary", return_value="/fake/mpnn"), \
         patch("utils.proteinmpnn_wrapper.subprocess.run",
               return_value=fake_cp):
        res = w.run_mpnn(backbone_pdb=str(bb), num_sequences=2,
                         output_dir=str(out))

    assert res.mock is False
    assert res.binary == "/fake/mpnn"
    assert len(res.sequences) == 2
    # Score parsing must pick up the value from the header.
    assert any(abs((s.get("score") or 0) - 0.8124) < 1e-3
               for s in res.sequences)


# ---------------------------------------------------------------------------
# T18 — AF2 wrapper invocation
# ---------------------------------------------------------------------------

def test_t18_af2_mock_fallback(tmp_path: Path) -> None:
    """Mock AF2 emits a tiny placeholder PDB + deterministic pLDDT/pTM."""
    w = AF2Wrapper(allow_mock=True)
    with patch.object(w, "_resolve_binary", return_value=None):
        res = w.run_af2(
            sequence="MKTLAEEDKL", recycles=3,
            output_dir=str(tmp_path / "af2_out"),
        )
    assert res.mock is True
    assert Path(res.pdb_path).exists()
    assert 50.0 <= res.plddt <= 95.0
    # pTM ∈ [0.5, 0.9] per mock formula.
    assert res.ptm is not None
    assert 0.45 <= res.ptm <= 0.95
    assert len(res.plddt_per_residue) == len("MKTLAEEDKL")
    assert res.schema == "stage12_pipeline/0.1"


def test_t18_af2_subprocess_mocked(tmp_path: Path) -> None:
    """Real AF2 path: subprocess returncode=0 + scores.json is parsed."""
    out = tmp_path / "af2_real"
    out.mkdir()
    pdb = out / "rank_1_model.pdb"
    pdb.write_text("ATOM      1  CA  GLY A   1\n")
    scores = out / "rank_1_model.json"
    scores.write_text(json.dumps({
        "plddt": [80.0, 81.0, 82.0, 83.0],
        "ptm": 0.72,
        "pae": [[0, 1], [1, 0]],
    }))

    fake_cp = MagicMock()
    fake_cp.returncode = 0
    fake_cp.stdout = ""
    fake_cp.stderr = ""

    w = AF2Wrapper(binary="/fake/af2")
    with patch.object(w, "_resolve_binary", return_value="/fake/af2"), \
         patch("utils.af2_wrapper.subprocess.run", return_value=fake_cp):
        res = w.run_af2(sequence="MKTL", recycles=2, output_dir=str(out))

    assert res.mock is False
    assert res.binary == "/fake/af2"
    assert abs(res.plddt - 81.5) < 1e-6
    assert abs((res.ptm or 0.0) - 0.72) < 1e-6
    assert len(res.pae_matrix) == 2


def test_t18_af2_invalid_sequence() -> None:
    """Non-alphabetic sequence raises ValueError."""
    w = AF2Wrapper(allow_mock=True)
    with pytest.raises(ValueError):
        w.run_af2(sequence="MKTL123", recycles=1)


# ---------------------------------------------------------------------------
# T19 — Full Stage 1-2 chain end-to-end via `discover --de-novo`
# ---------------------------------------------------------------------------

def test_t19_de_novo_full_chain(dummy_target: Path, out_dir: Path) -> None:
    """`discover --de-novo` runs RFdiffusion + ProteinMPNN + AF2 (all mocks
    when binaries absent) and produces a ranked candidate list."""
    rc = main([
        "discover",
        "--target", str(dummy_target),
        "--de-novo",
        "--num-designs", "2",
        "--num-sequences-per-backbone", "3",
        "--af2-recycles", "2",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "discover.json")
    assert payload["subcommand"] == "discover"
    assert payload["mode"] == "de_novo"
    assert payload["stage12"]["stage12_schema"] == "stage12_pipeline/0.1"
    # All 3 wrappers should have flagged mock_used (no binaries on this CI).
    for nm in ("rfdiffusion", "proteinmpnn", "af2"):
        assert payload["stage12"][nm]["mock_used"] is True
    # Two designs, three sequences each → up to 6 ranked entries (some may
    # be filtered out by pLDDT/pTM threshold; require >= 1).
    assert payload["n_designs"] == 2
    assert payload["n_ranked_candidates"] >= 1
    # Ranks must be 1..N with no duplicates.
    cands = payload["candidates"]
    ranks = [c["rank"] for c in cands]
    assert ranks == list(range(1, len(cands) + 1))
    # First ranked entry has the best (lowest) ΔG_mock.
    if len(cands) >= 2:
        dGs = [c["pipeline"]["dG_estimate"] for c in cands]
        assert dGs == sorted(dGs)
    md = (out_dir / "discover.md").read_text()
    assert "discover --de-novo" in md
    assert "Stage availability" in md


# ---------------------------------------------------------------------------
# T20 — --check-availability flags via direct module CLI
# ---------------------------------------------------------------------------

def test_t20_check_availability_all_three(tmp_path: Path) -> None:
    """Each wrapper's check_availability() returns the expected schema-keyed
    dict shape, both for "no binary" and patched-available cases."""
    # No binary case (default in CI).
    rf = RFdiffusionWrapper().check_availability()
    mp = ProteinMPNNWrapper().check_availability()
    af = AF2Wrapper().check_availability()
    for info in (rf, mp, af):
        assert set(info.keys()) >= {"available", "version", "path",
                                     "schema", "tool"}
        assert info["schema"] == "stage12_pipeline/0.1"
    # When binary exists, check_availability still returns a usable dict.
    fake_bin = tmp_path / "fake_rfdiff"
    fake_bin.write_text("#!/bin/sh\nexit 0\n")
    fake_bin.chmod(0o755)
    info2 = RFdiffusionWrapper(binary=str(fake_bin)).check_availability()
    assert info2["available"] is True
    assert info2["path"] == str(fake_bin)


# ---------------------------------------------------------------------------
# T21 — Failure handling: tool not found w/ allow_mock=False
# ---------------------------------------------------------------------------

def test_t21_rfdiffusion_strict_no_mock_raises(tmp_path: Path) -> None:
    w = RFdiffusionWrapper(allow_mock=False)
    with patch.object(w, "_resolve_binary", return_value=None):
        with pytest.raises(RFdiffusionNotFoundError):
            w.run_diffusion(target_pdb=str(tmp_path / "x.pdb"),
                            contigs="[A1-10/0 5-10]", num_designs=1)


def test_t21_proteinmpnn_strict_no_mock_raises(tmp_path: Path) -> None:
    bb = tmp_path / "bb.pdb"
    bb.write_text("ATOM      1  CA  GLY A   1\n")
    w = ProteinMPNNWrapper(allow_mock=False)
    with patch.object(w, "_resolve_binary", return_value=None):
        with pytest.raises(ProteinMPNNNotFoundError):
            w.run_mpnn(backbone_pdb=str(bb), num_sequences=1)


def test_t21_af2_strict_no_mock_raises(tmp_path: Path) -> None:
    w = AF2Wrapper(allow_mock=False)
    with patch.object(w, "_resolve_binary", return_value=None):
        with pytest.raises(AF2NotFoundError):
            w.run_af2(sequence="MKL", recycles=1)


def test_t21_subprocess_failure_raises(tmp_path: Path) -> None:
    """Real subprocess returning non-zero raises *InvocationError."""
    target = tmp_path / "x.pdb"
    target.write_text("HEADER\n")
    out = tmp_path / "out"
    out.mkdir()
    fake_cp = MagicMock()
    fake_cp.returncode = 1
    fake_cp.stdout = ""
    fake_cp.stderr = "boom\n"

    w = RFdiffusionWrapper(binary="/fake/rfdiffusion")
    with patch.object(w, "_resolve_binary",
                      return_value="/fake/rfdiffusion"), \
         patch("utils.rfdiffusion_wrapper.subprocess.run",
               return_value=fake_cp):
        with pytest.raises(RFdiffusionInvocationError):
            w.run_diffusion(target_pdb=str(target),
                            contigs="[A]", num_designs=1,
                            output_dir=str(out))


# ---------------------------------------------------------------------------
# T22 — Mock-mode fallback for demo (no tools installed)
# ---------------------------------------------------------------------------

def test_t22_demo_dry_run_no_subprocess(dummy_target: Path,
                                        out_dir: Path) -> None:
    """`--de-novo --dry-run` must NOT invoke subprocess.run on any wrapper —
    the dry-run branch returns a placeholder candidate immediately."""
    with patch("utils.rfdiffusion_wrapper.subprocess.run") as rf_run, \
         patch("utils.proteinmpnn_wrapper.subprocess.run") as mp_run, \
         patch("utils.af2_wrapper.subprocess.run") as af_run:
        rc = main([
            "discover",
            "--target", str(dummy_target),
            "--de-novo",
            "--num-designs", "3",
            "--dry-run",
            "--output", str(out_dir),
        ])
        # subprocess.run may be called by check_availability paths
        # (e.g. probing for --version) but never for the actual stage runs.
        # The tell-tale stage-1 invocation includes "inference.input_pdb=";
        # assert no such call is recorded.
        for run_mock in (rf_run, mp_run, af_run):
            for call in run_mock.call_args_list:
                argv = call.args[0] if call.args else []
                joined = " ".join(map(str, argv))
                assert "inference.input_pdb" not in joined
                assert "--num_seq_per_target" not in joined
                assert "--num-recycle" not in joined
    assert rc == 0
    payload = _read_json(out_dir / "discover.json")
    assert payload["mode"] == "de_novo"
    assert payload["n_designs"] >= 1
    # Dry-run candidates carry status="dry_run" and no ranked sequences.
    assert payload["designs"][0]["status"] == "dry_run"
    assert payload["n_ranked_candidates"] == 0


def test_t22_de_novo_chain_mock_flags_propagated(dummy_target: Path,
                                                 out_dir: Path) -> None:
    """End-to-end mock chain: each ranked candidate must report
    af2.mock=True, stage1_mock=True, stage1p5_mock=True so PR-NEW-E /
    downstream consumers can flag the demo nature unambiguously."""
    rc = main([
        "discover",
        "--target", str(dummy_target),
        "--de-novo",
        "--num-designs", "1",
        "--num-sequences-per-backbone", "2",
        "--output", str(out_dir),
    ])
    assert rc == 0
    payload = _read_json(out_dir / "discover.json")
    assert payload["stage12"]["rfdiffusion"]["mock_used"] is True
    designs = payload["designs"]
    assert designs and designs[0]["stage1_mock"] is True
    seqs = designs[0]["sequences"]
    assert seqs
    for s in seqs:
        # Every per-sequence record carries 1.5/2 mock flags.
        assert s.get("stage1p5_mock") is True
        if "af2" in s:
            assert s["af2"]["mock"] is True

