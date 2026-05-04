"""tests/test_html_report.py — PR-NEW-E HTML Report Generator tests (T1-T8).

Per ``utils/html_report.py`` (schema ``html_report/0.1``).

Test layout:
    T1: Single comparison HTML generation (synthetic BranchedDDGResult)
    T2: Multi-comparison ranking (3 variants → ranked table + matrix)
    T3: σ decomposition SVG chart validity
    T4: HTML basic validity (balanced tags, DOCTYPE, no Jinja2 fragments)
    T5: Theme parameterization (light vs dark CSS differs)
    T6: Output file write — write_report atomic + absolute path
    T7: Standalone (embed_assets=True) is the only Day-1 mode
    T8: PR-NEW-D integration (deferred — see Day-2 plan)
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

# Make ``utils.html_report`` importable.
_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from utils.html_report import (  # noqa: E402
    SCHEMA_VERSION,
    HTMLReportError,
    InputValidationError,
    RenderError,
    render_single_comparison,
    render_multi_ranking,
    render_cycle_progression,
    write_report,
    _validate_branched_ddg_input,
    _render_sigma_chart_svg,
    _inline_css,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_branched_ddg(
    *, target: str = "1EBP",
    var_label: str = "MTR13",
    wt_label: str = "WT",
    ddg: float = -2.5,
    se: float = 1.0,
    tier: str = "X.B",
    ci: tuple = (-4.5, 0.5),
    z: float = -2.5,
    sigma_btwn_var: float = 1.2,
    sigma_btwn_wt: float = 1.5,
    sigma_w_var: float = 5.0,
    sigma_w_wt: float = 5.5,
    n_seed: int = 5,
) -> dict:
    """Synthetic ``branched_ddg/0.1`` result dict for testing."""
    def _branch(label: str, mean: float, sb: float, sw: float) -> dict:
        return {
            "label": label,
            "branch_dir": f"/tmp/{label}",
            "per_seed": {},
            "mean_of_seeds": mean,
            "sigma_btwn": sb,
            "sigma_w_median": sw,
            "sigma_w_max": sw * 1.2,
            "SE": sb / (n_seed ** 0.5),
            "n_seed": n_seed,
            "n_snap_per_seed": 25,
            "solvent_model": "pbsa",
            "degenerate": False,
            "charge_audit": "PASS",
            "charge_audit_detail": {},
        }

    return {
        "schema": "branched_ddg/0.1",
        "tool": "utils.branched_ddg",
        "generated_at": "2026-04-28T12:00:00Z",
        "target_id": target,
        "solvent_model": "pbsa",
        "df_strategy": "conservative",
        "branches": {
            "wt": _branch(wt_label, -10.0, sigma_btwn_wt, sigma_w_wt),
            "variant": _branch(var_label, -10.0 + ddg,
                               sigma_btwn_var, sigma_w_var),
        },
        "pairwise": {
            "ddG": ddg,
            "SE": se,
            "df": 8,
            "t_crit": 2.306,
            "CI95": list(ci),
            "z_SE": z,
            "tier": tier,
            "sign_convention": "ddG = mean(variant) − mean(wt)",
            "note": None,
        },
        "charge_audit": "PASS",
        "charge_audit_strict": False,
        "wt_bit_identical": None,
        "wt_bit_identical_detail": None,
    }


# ---------------------------------------------------------------------------
# T1 — Single comparison HTML
# ---------------------------------------------------------------------------

def test_t1_single_comparison_html_generation():
    """T1: Single comparison HTML generation with synthetic input.

    Asserts: valid HTML doctype, contains <title>, contains <table>, contains
    ΔΔG numeric value, contains tier badge with correct CSS class, contains
    inline SVG.
    """
    r = _make_branched_ddg(ddg=-2.5, tier="X.B")
    out = render_single_comparison(r)
    assert isinstance(out, str)
    assert out.startswith("<!DOCTYPE html>"), \
        "missing HTML5 doctype"
    assert "<title>" in out and "</title>" in out
    assert "<table" in out and "</table>" in out
    # ΔΔG -2.5 should round-trip into the table.
    assert "-2.500" in out, "ΔΔG value not rendered with 3-digit format"
    # Tier badge for X.B → tier-B class.
    assert 'tier-B' in out, "tier-B CSS class not present"
    assert "X.B" in out
    # Inline SVG chart present.
    assert "<svg" in out and "</svg>" in out
    # Schema attribution.
    assert SCHEMA_VERSION in out
    # Branch labels present.
    assert "MTR13" in out and "WT" in out


def test_t1b_tier_color_class_each_tier():
    """T1b: Each tier maps to its corresponding CSS class."""
    for tier, cls in (("X.A", "tier-A"), ("X.B", "tier-B"),
                      ("X.C", "tier-C"), ("X.D", "tier-D")):
        r = _make_branched_ddg(tier=tier)
        out = render_single_comparison(r)
        assert cls in out, f"tier {tier} did not produce CSS class {cls}"


def test_t1c_dataclass_to_dict_input():
    """T1c: render_single_comparison accepts an object with .to_dict()."""
    r_dict = _make_branched_ddg(ddg=-1.0, tier="X.C")

    class _FakeResult:
        def __init__(self, d):
            self._d = d

        def to_dict(self):
            return self._d

    out = render_single_comparison(_FakeResult(r_dict))
    assert "<!DOCTYPE html>" in out
    assert "X.C" in out


# ---------------------------------------------------------------------------
# T2 — Multi-comparison ranking
# ---------------------------------------------------------------------------

def test_t2_multi_comparison_ranking():
    """T2: Multi-comparison ranking with 3 variants."""
    rs = [
        _make_branched_ddg(var_label="VarA", ddg=-3.0, tier="X.A",
                           ci=(-4.5, -1.5), z=-3.5),
        _make_branched_ddg(var_label="VarB", ddg=-1.0, tier="X.B",
                           ci=(-2.5, 0.5), z=-2.1),
        _make_branched_ddg(var_label="VarC", ddg=+0.5, tier="X.C",
                           ci=(-1.5, 2.5), z=+0.5),
    ]
    out = render_multi_ranking(rs)
    assert "<!DOCTYPE html>" in out
    # All 3 variants rendered.
    for label in ("VarA", "VarB", "VarC"):
        assert label in out, f"variant {label} not rendered"
    # Pairwise matrix table present (i \ j header sentinel).
    assert "i \\ j" in out, "pairwise matrix header missing"
    # Synthesis priority list present (ordered list).
    assert "<ol>" in out and "</ol>" in out
    # Best variant should be ranked #1 — VarA (X.A, most negative).
    rank_position = out.find("VarA")
    other_position = out.find("VarB")
    assert rank_position < other_position, \
        "VarA (X.A) should rank ahead of VarB (X.B)"
    # Tier badges all present.
    for cls in ("tier-A", "tier-B", "tier-C"):
        assert cls in out


def test_t2b_multi_empty_raises():
    """T2b: Empty results list raises InputValidationError."""
    with pytest.raises(InputValidationError):
        render_multi_ranking([])


# ---------------------------------------------------------------------------
# T3 — σ decomposition SVG chart
# ---------------------------------------------------------------------------

def test_t3_sigma_chart_svg_validity():
    """T3: SVG bar chart has 4 bars, axis labels, valid viewBox."""
    r = _make_branched_ddg(
        sigma_btwn_var=1.2, sigma_btwn_wt=1.5,
        sigma_w_var=5.0, sigma_w_wt=5.5,
    )
    svg = _render_sigma_chart_svg(r)
    assert svg.startswith('<svg')
    assert "viewBox=" in svg
    # Should contain 4 <rect> bars.
    n_rects = svg.count("<rect")
    assert n_rects == 4, f"expected 4 bars, got {n_rects}"
    # Axis labels and bar labels.
    assert "σ_btwn" in svg or "σ_w_med" in svg
    assert "kcal/mol" in svg
    # Numeric labels (each bar shows its value).
    assert "1.20" in svg or "1.50" in svg
    # SVG closed.
    assert svg.endswith("</svg>")


def test_t3b_sigma_chart_handles_missing_sigma_btwn():
    """T3b: SVG chart handles None sigma_btwn (degenerate branch) gracefully."""
    r = _make_branched_ddg()
    r["branches"]["variant"]["sigma_btwn"] = None  # degenerate
    svg = _render_sigma_chart_svg(r)
    assert "<svg" in svg
    # Missing-value bar shows "n/a" label.
    assert "n/a" in svg


# ---------------------------------------------------------------------------
# T4 — HTML structural validity
# ---------------------------------------------------------------------------

def test_t4_html_structural_validity():
    """T4: Generated HTML has balanced major tags, no Jinja2 fragments."""
    r = _make_branched_ddg()
    out = render_single_comparison(r)
    # No Jinja2 markers (we use f-strings only).
    for marker in ("{{", "}}", "{%", "%}"):
        assert marker not in out, f"unexpected Jinja2 marker {marker!r}"
    # Balance major structural tags.
    for tag in ("html", "head", "body", "header", "footer", "table",
                "section", "style"):
        opens = len(re.findall(rf"<{tag}\b", out))
        closes = len(re.findall(rf"</{tag}>", out))
        assert opens == closes, (
            f"unbalanced <{tag}>: {opens} opens vs {closes} closes"
        )


# ---------------------------------------------------------------------------
# T5 — Theme parameterization
# ---------------------------------------------------------------------------

def test_t5_theme_parameterization():
    """T5: Light and dark themes produce different CSS outputs."""
    css_light = _inline_css("light")
    css_dark = _inline_css("dark")
    assert css_light != css_dark, "themes are identical"
    # Light theme should have a near-white background, dark a near-black one.
    assert "#ffffff" in css_light
    assert "#1c1f24" in css_dark
    # End-to-end theme switch reflected in output.
    r = _make_branched_ddg()
    out_light = render_single_comparison(r, theme="light")
    out_dark = render_single_comparison(r, theme="dark")
    assert out_light != out_dark
    assert 'data-theme="light"' in out_light
    assert 'data-theme="dark"' in out_dark


def test_t5b_unknown_theme_raises():
    """T5b: Unknown theme name raises RenderError."""
    with pytest.raises(RenderError):
        _inline_css("solarized")


# ---------------------------------------------------------------------------
# T6 — Output file write
# ---------------------------------------------------------------------------

def test_t6_write_report_atomic_returns_absolute(tmp_path):
    """T6: write_report writes file at path, returns absolute, file readable."""
    html = "<!DOCTYPE html><html><body>x</body></html>"
    out_path = tmp_path / "subdir" / "report.html"
    abs_returned = write_report(html, str(out_path))
    assert Path(abs_returned).is_absolute()
    assert Path(abs_returned).exists()
    with open(abs_returned, "r", encoding="utf-8") as f:
        content = f.read()
    assert content == html
    # Subdir should have been created.
    assert (tmp_path / "subdir").is_dir()


def test_t6b_write_report_overwrites_existing(tmp_path):
    """T6b: Atomic overwrite on existing file."""
    out_path = tmp_path / "report.html"
    write_report("<html>v1</html>", str(out_path))
    assert out_path.read_text(encoding="utf-8") == "<html>v1</html>"
    write_report("<html>v2</html>", str(out_path))
    assert out_path.read_text(encoding="utf-8") == "<html>v2</html>"


# ---------------------------------------------------------------------------
# T7 — embed_assets behaviour
# ---------------------------------------------------------------------------

def test_t7_embed_assets_required_day1():
    """T7: embed_assets=True is the only Day-1 mode (False raises)."""
    r = _make_branched_ddg()
    # Default is True → succeeds.
    out = render_single_comparison(r)
    assert "<style>" in out  # CSS embedded inline
    assert "<svg" in out      # SVG embedded inline
    # Explicit False is reserved for Day 2 → RenderError.
    with pytest.raises(RenderError):
        render_single_comparison(r, embed_assets=False)
    with pytest.raises(RenderError):
        render_multi_ranking([r], embed_assets=False)


# ---------------------------------------------------------------------------
# T8 — PR-NEW-D integration (deferred Day 2)
# ---------------------------------------------------------------------------

@pytest.mark.skip(
    reason="PR-NEW-D scripts/updd_cli.py `report --report-format html` "
    "currently falls back to markdown by design (Day-1 hook). Wiring "
    "`--report-format html` into html_report.render_single_comparison is "
    "deferred to PR-NEW-E Day 2."
)
def test_t8_pr_new_d_integration_html_format(tmp_path):
    """T8: scripts/updd_cli.py report --format html (deferred)."""
    pass


# ---------------------------------------------------------------------------
# Schema validation tests (input_validation)
# ---------------------------------------------------------------------------

def test_schema_validation_rejects_wrong_schema():
    """Schema mismatch is caught up-front."""
    bad = {"schema": "something_else/0.1"}
    with pytest.raises(InputValidationError):
        _validate_branched_ddg_input(bad)


def test_schema_validation_rejects_missing_branches():
    """Missing 'branches' key raises."""
    bad = {"schema": "branched_ddg/0.1"}
    with pytest.raises(InputValidationError):
        _validate_branched_ddg_input(bad)


def test_schema_validation_accepts_full_record():
    """Full synthetic record validates cleanly."""
    r = _make_branched_ddg()
    _validate_branched_ddg_input(r)  # must not raise


# ---------------------------------------------------------------------------
# Cycle progression stub (smoke)
# ---------------------------------------------------------------------------

def test_cycle_progression_stub_renders():
    """Day-2 stub renders without raising and is valid HTML."""
    history = {
        "schema": "cycle_manager/0.1",
        "project_id": "test_proj",
        "target_pdb": "/tmp/2QKI.pdb",
        "created_at": "2026-04-28T00:00:00Z",
        "cycles": [
            {"cycle_n": 0, "variants": [], "predicted_ddg": [],
             "wet_lab_results": [], "synthesized": False},
            {"cycle_n": 1, "variants": [], "predicted_ddg": [],
             "wet_lab_results": [], "synthesized": False},
        ],
    }
    out = render_cycle_progression(history)
    assert "<!DOCTYPE html>" in out
    assert "test_proj" in out
    assert "Cycle Progression" in out
    # Closing tags balanced.
    assert out.count("<html") == out.count("</html>")


# ---------------------------------------------------------------------------
# CLI driver smoke (uses subprocess on the module)
# ---------------------------------------------------------------------------

def test_cli_single_mode(tmp_path):
    """CLI driver writes a valid HTML file."""
    r = _make_branched_ddg()
    in_json = tmp_path / "input.json"
    out_html = tmp_path / "out.html"
    with open(in_json, "w", encoding="utf-8") as f:
        json.dump(r, f)
    cmd = [
        sys.executable, "-m", "utils.html_report",
        "--result", str(in_json),
        "--output", str(out_html),
        "--mode", "single",
        "--theme", "light",
    ]
    proc = subprocess.run(
        cmd, cwd=str(_REPO_ROOT),
        capture_output=True, text=True, timeout=30,
    )
    assert proc.returncode == 0, (
        f"CLI failed: stdout={proc.stdout!r}, stderr={proc.stderr!r}"
    )
    assert out_html.is_file()
    text = out_html.read_text(encoding="utf-8")
    assert "<!DOCTYPE html>" in text
    assert "X.B" in text  # default tier from synthetic record


# ---------------------------------------------------------------------------
# T9-T14 — PR-NEW-E Day 2 additions
# ---------------------------------------------------------------------------

def _make_cycle_history(*, schema: str = "cycle_manager/0.2") -> dict:
    """Synthetic 3-cycle history with K_d 100nM → 20nM → 5nM (no plateau)."""
    return {
        "schema": schema,
        "project_id": "test_proj",
        "target_pdb": "/tmp/2QKI.pdb",
        "created_at": "2026-04-29T00:00:00Z",
        "metadata": {"operator": "san", "ncaa_budget": 5},
        "cycles": [
            {
                "cycle_n": 0,
                "wt_sequence": None,
                "variants": ["ICVVQDWGHHRCT"],
                "predicted_ddg": [
                    {"sequence": "ICVVQDWGHHRCT", "ddg": -1.5, "se": 0.4,
                     "tier": "X.B", "z_se": -3.75, "ci95": [-2.5, -0.5]},
                ],
                "wet_lab_results": [
                    {"sequence": "ICVVQDWGHHRCT", "kd_molar": 1e-7,
                     "assay": "ITC", "ddg_calc": -9.57,
                     "temperature_k": 298.15, "confidence": None,
                     "notes": None},
                ],
                "best_validated_sequence": "ICVVQDWGHHRCT",
                "synthesized": True,
                "notes": None,
                "created_at": "2026-04-29T00:00:00Z",
            },
            {
                "cycle_n": 1,
                "wt_sequence": "ICVVQDWGHHRCT",
                "variants": ["ICVVQDW(MTR)HHRCT", "ICVVQDW(NMA)HHRCT"],
                "predicted_ddg": [
                    {"sequence": "ICVVQDW(MTR)HHRCT", "ddg": -2.5, "se": 0.5,
                     "tier": "X.A", "z_se": -5.0, "ci95": [-3.6, -1.4]},
                    {"sequence": "ICVVQDW(NMA)HHRCT", "ddg": -1.0, "se": 0.6,
                     "tier": "X.B", "z_se": -1.67, "ci95": [-2.3, 0.3]},
                ],
                "wet_lab_results": [
                    {"sequence": "ICVVQDW(MTR)HHRCT", "kd_molar": 2e-8,
                     "assay": "ITC", "ddg_calc": -10.52,
                     "temperature_k": 298.15, "confidence": None,
                     "notes": None},
                ],
                "best_validated_sequence": "ICVVQDW(MTR)HHRCT",
                "synthesized": True,
                "notes": None,
                "created_at": "2026-04-29T00:00:00Z",
            },
            {
                "cycle_n": 2,
                "wt_sequence": "ICVVQDW(MTR)HHRCT",
                "variants": ["ICVVQDW(MTR)H(MEH)HRCT", "ICVVQDW(MTR)HHR(NMA)CT"],
                "predicted_ddg": [
                    {"sequence": "ICVVQDW(MTR)H(MEH)HRCT", "ddg": -1.2,
                     "se": 0.3, "tier": "X.A", "z_se": -4.0,
                     "ci95": [-1.8, -0.6]},
                    {"sequence": "ICVVQDW(MTR)HHR(NMA)CT", "ddg": -0.8,
                     "se": 0.4, "tier": "X.B", "z_se": -2.0,
                     "ci95": [-1.6, 0.0]},
                ],
                "wet_lab_results": [
                    {"sequence": "ICVVQDW(MTR)H(MEH)HRCT", "kd_molar": 5e-9,
                     "assay": "SPR", "ddg_calc": -11.34,
                     "temperature_k": 298.15, "confidence": None,
                     "notes": None},
                ],
                "best_validated_sequence": "ICVVQDW(MTR)H(MEH)HRCT",
                "synthesized": True,
                "notes": None,
                "created_at": "2026-04-29T00:00:00Z",
            },
        ],
        "trajectory": {
            "per_cycle": [
                {"cycle_n": 0, "kd_molar": 1e-7, "ddg": -9.57,
                 "sequence": "ICVVQDWGHHRCT"},
                {"cycle_n": 1, "kd_molar": 2e-8, "ddg": -10.52,
                 "sequence": "ICVVQDW(MTR)HHRCT"},
                {"cycle_n": 2, "kd_molar": 5e-9, "ddg": -11.34,
                 "sequence": "ICVVQDW(MTR)H(MEH)HRCT"},
            ],
            "delta_kd_kcal": [-0.95, -0.82],
            "ratio_per_cycle": [0.2, 0.25],
            "plateau_detected": False,
            "lookback": 3,
            "threshold_kcal": 0.5,
            "recommendation": "continue",
            "final_ddg_kcal": -11.34,
            "computed_at": "2026-04-29T00:00:00Z",
        },
        "output_root": "outputs/projects",
    }


def test_t9_cycle_progression_full_render():
    """T9: Cycle progression renderer produces a full document with every
    section: header, per-cycle table, K_d timeline SVG, plateau badge,
    synthesis priority list.
    """
    h = _make_cycle_history()
    out = render_cycle_progression(h)
    assert "<!DOCTYPE html>" in out
    assert out.count("<html") == out.count("</html>")
    # Header.
    assert "test_proj" in out
    assert "Cycle Progression" in out
    # Per-cycle table — rows for cycle 0/1/2.
    assert "<table" in out and "</table>" in out
    for cycle_n in ("0", "1", "2"):
        # Look for ``<td>0</td>`` / ``<td>1</td>`` / ``<td>2</td>`` in body.
        assert f"<td>{cycle_n}</td>" in out, (
            f"per-cycle row for cycle {cycle_n} not rendered"
        )
    # Per-cycle table includes best validated sequences.
    assert "ICVVQDW(MTR)H(MEH)HRCT" in out
    # K_d timeline SVG present.
    assert "kd-timeline-chart" in out
    assert "<svg" in out and "</svg>" in out
    # Plateau badge — recommendation == "continue".
    assert "plateau-badge" in out
    assert "Continue" in out
    # Synthesis priority list (latest cycle has 2 predicted ΔΔG → ordered list).
    assert "<ol>" in out and "</ol>" in out
    # Schema attribution.
    assert SCHEMA_VERSION in out


def test_t10_cycle_progression_pr_new_c_integration():
    """T10: PR-NEW-C integration — render consumes a CycleHistory.to_dict()
    output (mock) and produces all expected sections.
    """
    h = _make_cycle_history()

    class _FakeHistory:
        def __init__(self, d):
            self._d = d

        def to_dict(self):
            return self._d

    out = render_cycle_progression(_FakeHistory(h))
    assert "<!DOCTYPE html>" in out
    # Every renamed section class is present.
    for cls in ("cycle-overview", "cycle-table", "cycle-timeline",
                "cycle-recommendation", "cycle-priority"):
        assert cls in out, f"section class {cls!r} missing"
    # K_d formatting via the cycle_manager-aligned _fmt_kd_human.
    # 100 nM → '100 nM', 20 nM → '20 nM', 5 nM → '5 nM'.
    for kd_label in ("100 nM", "20 nM", "5 nM"):
        assert kd_label in out, f"K_d label {kd_label!r} not rendered"


def test_t10b_cycle_progression_plateau_recommendations():
    """T10b: Each plateau recommendation maps to a distinct badge class."""
    cases = [
        ("converged", "converged", "CONVERGED"),
        ("try_alternative", "try-alternative", "PLATEAU DETECTED"),
        ("continue", "continue", "Continue"),
    ]
    for rec, css_cls, badge_text in cases:
        h = _make_cycle_history()
        h["trajectory"]["recommendation"] = rec
        h["trajectory"]["plateau_detected"] = (rec != "continue")
        out = render_cycle_progression(h)
        assert f"plateau-badge {css_cls}" in out, (
            f"recommendation {rec!r} did not produce CSS class {css_cls!r}"
        )
        assert badge_text in out, (
            f"recommendation {rec!r} did not produce badge text {badge_text!r}"
        )


def test_t11_pr_new_d_html_format_wireup(tmp_path):
    """T11: scripts/updd_cli.py report --report-format html — Phase B
    integration test.

    Drives the production CLI through subprocess and verifies that html
    output is rendered via utils.html_report (rather than falling back to
    markdown).
    """
    r = _make_branched_ddg(ddg=-2.5, tier="X.B")
    in_json = tmp_path / "in.json"
    out_dir = tmp_path / "out"
    with open(in_json, "w", encoding="utf-8") as f:
        json.dump(r, f)
    cmd = [
        sys.executable, "scripts/updd_cli.py", "report",
        "--result", str(in_json),
        "--output", str(out_dir),
        "--report-format", "html",
        "--theme", "light",
    ]
    proc = subprocess.run(
        cmd, cwd=str(_REPO_ROOT),
        capture_output=True, text=True, timeout=30,
    )
    assert proc.returncode == 0, (
        f"updd_cli failed: stdout={proc.stdout!r} stderr={proc.stderr!r}"
    )
    out_html = out_dir / "report.html"
    assert out_html.is_file(), "report.html not written"
    text = out_html.read_text(encoding="utf-8")
    assert "<!DOCTYPE html>" in text
    assert SCHEMA_VERSION in text
    assert "X.B" in text
    # Output stdout mentions writing the HTML file.
    assert "report.html" in proc.stdout


def _make_branched_ddg_with_detachment(warning: str = "high_detachment",
                                       fraction: float = 0.65) -> dict:
    """Synthetic v0.2 branched_ddg result with a detachment_metric block."""
    r = _make_branched_ddg(var_label="MTR13", ddg=-2.5, tier="X.B")
    r["detachment_metric"] = {
        "per_seed_fraction": [0.6, 0.7, 0.65],
        "per_seed_detail": [
            {"snapshots_dir": "/tmp/s1", "n_frames": 25, "n_detached": 15,
             "fraction": 0.6},
            {"snapshots_dir": "/tmp/s2", "n_frames": 25, "n_detached": 17,
             "fraction": 0.68},
            {"snapshots_dir": "/tmp/s3", "n_frames": 25, "n_detached": 16,
             "fraction": 0.64},
        ],
        "aggregate_fraction": fraction,
        "threshold_angstrom": 5.0,
        "binder_chain": "B",
        "n_frames": 75,
        "n_detached": int(round(fraction * 75)),
        "n_skipped": 0,
        "warning": warning,
        "ncaa_resname": "MTR",
        "ncaa_resseq": 13,
    }
    return r


def test_t12_detachment_metric_banner_high():
    """T12: BranchedDDGResult with detachment_metric warning='high_detachment'
    produces a visible red banner; moderate produces a yellow banner; 'none'
    produces no banner.
    """
    r_high = _make_branched_ddg_with_detachment("high_detachment", 0.65)
    out_high = render_single_comparison(r_high)
    assert "detachment-banner high" in out_high
    assert "HIGH DETACHMENT" in out_high
    # Per-branch fraction details displayed.
    assert "0.65" in out_high
    assert "MTR" in out_high  # ncAA resname surfaced

    r_mod = _make_branched_ddg_with_detachment("moderate_detachment", 0.30)
    out_mod = render_single_comparison(r_mod)
    assert "detachment-banner moderate" in out_mod
    assert "MODERATE DETACHMENT" in out_mod

    r_none = _make_branched_ddg_with_detachment("none", 0.05)
    out_none = render_single_comparison(r_none)
    # No detachment banner when warning == "none".
    assert "detachment-banner high" not in out_none
    assert "detachment-banner moderate" not in out_none


def test_t12b_detachment_metric_in_multi_ranking():
    """T12b: Multi-ranking surfaces a detachment-summary section when at
    least one ranked result carries a non-'none' warning.
    """
    r1 = _make_branched_ddg_with_detachment("high_detachment", 0.65)
    r1["branches"]["variant"]["label"] = "VarA"
    r2 = _make_branched_ddg(var_label="VarB", ddg=-1.0, tier="X.B")
    out = render_multi_ranking([r1, r2])
    assert "detachment-summary" in out
    assert "Detachment Advisory" in out
    assert "VarA" in out


def test_t13_print_css_block_present():
    """T13: Output contains a comprehensive ``@media print`` block with
    page-break management, navigation/footer hidden, and black-on-white
    fallback.
    """
    r = _make_branched_ddg()
    out = render_single_comparison(r)
    assert "@media print" in out
    assert "page-break-before" in out or "page-break-inside" in out
    # Footer/nav hidden in print.
    assert "display: none" in out
    # Black-on-white fallback.
    assert "#ffffff" in out and "#000000" in out
    # Responsive design (mobile).
    assert "@media (max-width: 768px)" in out


def test_t14_backward_compat_v01_branched_ddg():
    """T14: A v0.1-shaped BranchedDDGResult (no detachment_metric, no
    cycle history) still renders cleanly under the v0.2 reader.
    """
    r = _make_branched_ddg()  # no detachment_metric in fixture
    out = render_single_comparison(r)
    assert "<!DOCTYPE html>" in out
    # No detachment banner since metric is absent.
    assert "detachment-banner high" not in out
    assert "detachment-banner moderate" not in out

    # And v0.1 cycle_manager input renders too (no trajectory block).
    h_v01 = {
        "schema": "cycle_manager/0.1",
        "project_id": "compat_proj",
        "target_pdb": "/tmp/T.pdb",
        "created_at": "2026-04-28T00:00:00Z",
        "cycles": [
            {"cycle_n": 0, "variants": ["X"], "predicted_ddg": [],
             "wet_lab_results": [], "synthesized": False},
        ],
    }
    out2 = render_cycle_progression(h_v01)
    assert "<!DOCTYPE html>" in out2
    assert "compat_proj" in out2
    # Without trajectory, plateau section renders the no-trajectory note.
    assert ("No trajectory data" in out2
            or "no trajectory" in out2.lower())


def test_t14b_v02_schema_in_output():
    """T14b: Generated documents declare html_report/0.2 schema."""
    r = _make_branched_ddg()
    out = render_single_comparison(r)
    assert "html_report/0.2" in out
