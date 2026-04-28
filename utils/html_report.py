#!/usr/bin/env python
"""utils/html_report.py — UPDD HTML Report Generator (industry-adoption polish).

Schema: ``html_report/0.1`` (PR-NEW-E Day 1, 2026-04-28)

Generates standalone HTML reports from PR-NEW-B :class:`BranchedDDGResult`
inputs (and, in a Day-2 stub, PR-NEW-C :class:`CycleHistory` inputs):

    * Single comparison       — WT vs single variant
    * Multi-ranking            — N variants ranked by ΔΔG (with pairwise matrix)
    * Cycle progression (stub) — Cycle 0 → N timeline, K_d trajectory

Design constraints (Day 1):
    * **Stdlib only.** No Jinja2, no matplotlib. Charts are inline ``<svg>``.
    * **Single-file output.** ``embed_assets=True`` (default and only Day-1 mode)
      writes one HTML document with embedded CSS + inline SVG. The
      ``embed_assets=False`` path is a Day-2 stub.
    * **Print-friendly.** Uses screen-only and ``@media print`` blocks so the
      same file can be PDF-exported via the browser.
    * **Light + dark themes.** ``theme="light"`` (default) or ``theme="dark"``.

Public API:
    * :class:`HTMLReportError` and subclasses
      (:class:`InputValidationError`, :class:`RenderError`)
    * :func:`render_single_comparison(result, *, theme, embed_assets) -> str`
    * :func:`render_multi_ranking(results, *, theme, embed_assets) -> str`
    * :func:`render_cycle_progression(history, *, theme) -> str` (Day-2 stub)
    * :func:`write_report(html, output_path) -> str`  (atomic write)

CLI::

    python -m utils.html_report \\
        --result <branched_ddg_result.json> \\
        --output report.html \\
        [--mode {single,multi,cycle}] \\
        [--theme {light,dark}]

Tier color palette (must match PR-NEW-B :func:`tier_classify` semantics):

    X.A — green   (CI95 excludes 0; sign + magnitude significant)
    X.B — amber   (|z_SE| ≥ 2.0, CI includes 0; sign-significant only)
    X.C — gray    (insufficient evidence)
    X.D — light-gray (insufficient sampling, n_seed < 3)
"""

from __future__ import annotations

import argparse
import html
import json
import math
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

# ---------------------------------------------------------------------------
# Constants / Schema
# ---------------------------------------------------------------------------

SCHEMA_VERSION = "html_report/0.1"
TOOL_NAME = "utils.html_report"

# Recognized tier set (mirror utils.branched_ddg.tier_classify).
_TIERS = ("X.A", "X.B", "X.C", "X.D")

# Tier label for accessibility / non-CSS readers.
_TIER_LABELS: Dict[str, str] = {
    "X.A": "Recommended for synthesis (sign + magnitude significant)",
    "X.B": "Recommended for synthesis (sign-significant only)",
    "X.C": "Insufficient evidence — collect more sampling",
    "X.D": "Insufficient sampling (n_seed < 3) — increase replicates",
}

# Tier short verdicts (used in the verdict line).
_TIER_VERDICTS: Dict[str, str] = {
    "X.A": "is recommended for synthesis",
    "X.B": "is conditionally recommended for synthesis",
    "X.C": "should not be prioritized without additional sampling",
    "X.D": "has insufficient sampling for a reliable verdict",
}


# ---------------------------------------------------------------------------
# Public exception classes
# ---------------------------------------------------------------------------

class HTMLReportError(Exception):
    """Base class for HTML report errors."""


class InputValidationError(HTMLReportError):
    """Raised when an input dict does not match the expected schema."""


class RenderError(HTMLReportError):
    """Raised when an internal render step fails."""


# ---------------------------------------------------------------------------
# Internal helpers — formatting
# ---------------------------------------------------------------------------

def _html_escape(s: Any) -> str:
    """Stdlib HTML escape wrapper."""
    if s is None:
        return ""
    return html.escape(str(s), quote=True)


def _format_float(x: Any, *, sig: int = 3, signed: bool = True) -> str:
    """Consistent number formatting for tables/charts.

    Returns an em-dash for ``None``/non-numeric. Uses fixed-precision
    (``sig`` digits) and an explicit ``+/-`` sign when ``signed=True``.
    """
    if x is None:
        return "—"
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return _html_escape(x)
    if not math.isfinite(xf):
        return "nan/inf"
    if signed:
        return f"{xf:+.{sig}f}"
    return f"{xf:.{sig}f}"


def _format_ci(ci: Any, *, sig: int = 3) -> str:
    """Format ``[lo, hi]`` confidence interval."""
    if not isinstance(ci, (list, tuple)) or len(ci) != 2:
        return "—"
    lo, hi = ci
    return f"[{_format_float(lo, sig=sig)}, {_format_float(hi, sig=sig)}]"


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ---------------------------------------------------------------------------
# Internal helpers — schema validation
# ---------------------------------------------------------------------------

def _validate_branched_ddg_input(d: Dict[str, Any]) -> None:
    """Validate that ``d`` matches schema ``branched_ddg/0.1`` (PR-NEW-B).

    Raises :class:`InputValidationError` on missing required fields.
    Does NOT enforce strict schema; tolerates extra keys.
    """
    if not isinstance(d, dict):
        raise InputValidationError(
            f"expected dict for branched_ddg input, got {type(d).__name__}"
        )
    schema = d.get("schema") or d.get("schema_version") or ""
    if not str(schema).startswith("branched_ddg/"):
        raise InputValidationError(
            f"schema mismatch: expected 'branched_ddg/*', got {schema!r}"
        )
    branches = d.get("branches")
    if not isinstance(branches, dict):
        raise InputValidationError("missing 'branches' object")
    for side in ("wt", "variant"):
        if side not in branches or not isinstance(branches[side], dict):
            raise InputValidationError(f"missing branch: {side}")
        b = branches[side]
        for f in ("label", "mean_of_seeds", "sigma_w_median",
                  "n_seed", "solvent_model"):
            if f not in b:
                raise InputValidationError(
                    f"branch {side!r} missing field {f!r}"
                )
    pw = d.get("pairwise")
    if not isinstance(pw, dict):
        raise InputValidationError("missing 'pairwise' object")
    for f in ("ddG", "SE", "tier"):
        if f not in pw:
            raise InputValidationError(f"pairwise missing field {f!r}")


def _validate_cycle_history_input(d: Dict[str, Any]) -> None:
    """Validate that ``d`` matches schema ``cycle_manager/0.1`` (PR-NEW-C)."""
    if not isinstance(d, dict):
        raise InputValidationError(
            f"expected dict for cycle_history input, got {type(d).__name__}"
        )
    schema = d.get("schema") or d.get("schema_version") or ""
    if not str(schema).startswith("cycle_manager/"):
        raise InputValidationError(
            f"schema mismatch: expected 'cycle_manager/*', got {schema!r}"
        )
    if "cycles" not in d or not isinstance(d["cycles"], list):
        raise InputValidationError("missing 'cycles' list")


def _load_input(path: str) -> Dict[str, Any]:
    """Load a JSON input file. Returns its top-level dict."""
    p = Path(path)
    if not p.is_file():
        raise InputValidationError(f"input file not found: {path}")
    try:
        with open(p, "r", encoding="utf-8") as f:
            return json.load(f)
    except json.JSONDecodeError as exc:
        raise InputValidationError(f"input is not valid JSON: {path}: {exc}")


# ---------------------------------------------------------------------------
# Internal helpers — CSS theme
# ---------------------------------------------------------------------------

def _inline_css(theme: str = "light") -> str:
    """Return an embedded ``<style>`` block. ``theme`` is ``light`` or ``dark``.

    Tier color tokens stay constant across themes (consistent semantic
    meaning); foreground/background/borders flip.
    """
    if theme not in ("light", "dark"):
        raise RenderError(f"unknown theme: {theme!r}")

    light = {
        "bg": "#ffffff",
        "fg": "#1a1a1a",
        "muted": "#5a5a5a",
        "panel": "#fafafa",
        "border": "#d0d0d0",
        "table_header_bg": "#eef2f7",
        "table_row_alt_bg": "#f7f7f9",
        "footer_fg": "#666666",
        "header_bg": "#f0f4fa",
        "axis": "#333333",
        "code_bg": "#f4f4f4",
    }
    dark = {
        "bg": "#1c1f24",
        "fg": "#e8e8e8",
        "muted": "#a8a8a8",
        "panel": "#262a31",
        "border": "#3a3f47",
        "table_header_bg": "#2d343d",
        "table_row_alt_bg": "#23272d",
        "footer_fg": "#9a9a9a",
        "header_bg": "#222730",
        "axis": "#cfcfcf",
        "code_bg": "#2a2f37",
    }
    palette = light if theme == "light" else dark

    return f"""<style>
:root {{
  --bg: {palette["bg"]};
  --fg: {palette["fg"]};
  --muted: {palette["muted"]};
  --panel: {palette["panel"]};
  --border: {palette["border"]};
  --th-bg: {palette["table_header_bg"]};
  --row-alt-bg: {palette["table_row_alt_bg"]};
  --footer-fg: {palette["footer_fg"]};
  --header-bg: {palette["header_bg"]};
  --axis: {palette["axis"]};
  --code-bg: {palette["code_bg"]};
  --tier-A-bg: #d6f5d6;
  --tier-A-fg: #155724;
  --tier-A-bd: #2d7a32;
  --tier-B-bg: #fff3cd;
  --tier-B-fg: #856404;
  --tier-B-bd: #b07d10;
  --tier-C-bg: #e2e3e5;
  --tier-C-fg: #4a4a4a;
  --tier-C-bd: #6c757d;
  --tier-D-bg: #f5f5f5;
  --tier-D-fg: #777777;
  --tier-D-bd: #aaaaaa;
  --bar-variant: #4f7fbf;
  --bar-wt: #b07d10;
}}
* {{ box-sizing: border-box; }}
html, body {{
  margin: 0;
  padding: 0;
  background: var(--bg);
  color: var(--fg);
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Helvetica Neue",
               Arial, sans-serif;
  font-size: 14px;
  line-height: 1.5;
}}
.report {{
  max-width: 920px;
  margin: 0 auto;
  padding: 24px 28px;
}}
header.report-header {{
  background: var(--header-bg);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 18px 22px;
  margin-bottom: 22px;
}}
header.report-header h1 {{
  margin: 0 0 6px 0;
  font-size: 22px;
  font-weight: 600;
}}
header.report-header .meta {{
  color: var(--muted);
  font-size: 12.5px;
}}
section {{
  margin-bottom: 22px;
  border: 1px solid var(--border);
  border-radius: 8px;
  background: var(--panel);
  padding: 16px 20px;
}}
section h2 {{
  margin: 0 0 12px 0;
  font-size: 17px;
  font-weight: 600;
  border-bottom: 1px solid var(--border);
  padding-bottom: 6px;
}}
table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 13px;
}}
th, td {{
  border: 1px solid var(--border);
  padding: 6px 10px;
  text-align: left;
  vertical-align: top;
}}
th {{
  background: var(--th-bg);
  font-weight: 600;
}}
tbody tr:nth-child(even) td {{
  background: var(--row-alt-bg);
}}
.summary-table td:nth-child(2) {{
  font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
  white-space: nowrap;
}}
.code, code {{
  font-family: ui-monospace, "SF Mono", Menlo, Consolas, monospace;
  background: var(--code-bg);
  padding: 1px 5px;
  border-radius: 3px;
  font-size: 12.5px;
}}
.tier-badge {{
  display: inline-block;
  padding: 2px 10px;
  border-radius: 12px;
  font-weight: 600;
  font-size: 12.5px;
  border-style: solid;
  border-width: 1px;
}}
.tier-A {{
  background: var(--tier-A-bg);
  color: var(--tier-A-fg);
  border-color: var(--tier-A-bd);
}}
.tier-B {{
  background: var(--tier-B-bg);
  color: var(--tier-B-fg);
  border-color: var(--tier-B-bd);
}}
.tier-C {{
  background: var(--tier-C-bg);
  color: var(--tier-C-fg);
  border-color: var(--tier-C-bd);
}}
.tier-D {{
  background: var(--tier-D-bg);
  color: var(--tier-D-fg);
  border-color: var(--tier-D-bd);
}}
section.verdict.tier-A {{ border-left: 5px solid var(--tier-A-bd); }}
section.verdict.tier-B {{ border-left: 5px solid var(--tier-B-bd); }}
section.verdict.tier-C {{ border-left: 5px solid var(--tier-C-bd); }}
section.verdict.tier-D {{ border-left: 5px solid var(--tier-D-bd); }}
section.verdict p {{
  margin: 6px 0 0 0;
  font-size: 14.5px;
}}
.sigma-chart {{
  display: block;
  margin: 6px auto 0 auto;
  max-width: 100%;
}}
.kv-grid {{
  display: grid;
  grid-template-columns: max-content auto;
  gap: 4px 14px;
  margin-bottom: 8px;
  font-size: 13px;
}}
.kv-grid .k {{ color: var(--muted); }}
.kv-grid .v {{ font-family: ui-monospace, monospace; }}
footer.report-footer {{
  margin-top: 28px;
  border-top: 1px solid var(--border);
  padding-top: 10px;
  color: var(--footer-fg);
  font-size: 11.5px;
  text-align: center;
}}
@media print {{
  body {{ background: white; color: black; }}
  section, header.report-header {{
    border-color: #888;
    background: white !important;
  }}
  .report {{ max-width: none; padding: 12px 16px; }}
}}
</style>"""


# ---------------------------------------------------------------------------
# Internal helpers — SVG charting
# ---------------------------------------------------------------------------

def _safe_get(d: Dict[str, Any], *keys: str, default: Any = None) -> Any:
    """Walk nested dicts safely. Returns ``default`` on any miss."""
    cur: Any = d
    for k in keys:
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def _branch_sigma_btwn(branch: Dict[str, Any]) -> Optional[float]:
    """Read ``sigma_btwn`` from a branch dict; ``None`` if degenerate."""
    s = branch.get("sigma_btwn")
    if s is None:
        return None
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _branch_sigma_w_median(branch: Dict[str, Any]) -> float:
    """Read ``sigma_w_median`` from a branch dict; 0.0 fallback."""
    try:
        return float(branch.get("sigma_w_median", 0.0) or 0.0)
    except (TypeError, ValueError):
        return 0.0


def _render_sigma_chart_svg(result: Dict[str, Any], *,
                            width: int = 520, height: int = 240) -> str:
    """Inline SVG bar chart for σ_btwn vs σ_w_median (variant + WT side-by-side).

    4 bars: variant σ_btwn, variant σ_w_median, wt σ_btwn, wt σ_w_median.
    Y-axis auto-scales to ``max(all σ) * 1.2``. Axis labels are SVG ``<text>``.
    """
    branches = result.get("branches") or {}
    var = branches.get("variant") or {}
    wt = branches.get("wt") or {}
    bars: List[Tuple[str, Optional[float], str]] = [
        ("var σ_btwn", _branch_sigma_btwn(var), "var(--bar-variant)"),
        ("var σ_w_med", _branch_sigma_w_median(var), "var(--bar-variant)"),
        ("wt σ_btwn", _branch_sigma_btwn(wt), "var(--bar-wt)"),
        ("wt σ_w_med", _branch_sigma_w_median(wt), "var(--bar-wt)"),
    ]
    # Substitute None → 0 for chart layout, but mark the bar.
    numeric = [(lbl, (0.0 if v is None else v), col, v is None)
               for lbl, v, col in bars]
    finite = [v for _, v, _, miss in numeric if not miss]
    if finite:
        ymax = max(finite) * 1.2
        if ymax <= 0:
            ymax = 1.0
    else:
        ymax = 1.0

    margin_l = 56
    margin_r = 16
    margin_t = 18
    margin_b = 50
    chart_w = width - margin_l - margin_r
    chart_h = height - margin_t - margin_b

    n = len(numeric)
    slot_w = chart_w / n
    bar_w = slot_w * 0.55

    parts: List[str] = []
    parts.append(
        f'<svg class="sigma-chart" viewBox="0 0 {width} {height}" '
        f'xmlns="http://www.w3.org/2000/svg" role="img" '
        f'aria-label="sigma decomposition bar chart">'
    )
    # Axes.
    parts.append(
        f'<line x1="{margin_l}" y1="{margin_t}" '
        f'x2="{margin_l}" y2="{margin_t + chart_h}" '
        f'stroke="var(--axis)" stroke-width="1"/>'
    )
    parts.append(
        f'<line x1="{margin_l}" y1="{margin_t + chart_h}" '
        f'x2="{margin_l + chart_w}" y2="{margin_t + chart_h}" '
        f'stroke="var(--axis)" stroke-width="1"/>'
    )
    # Y-axis label.
    parts.append(
        f'<text x="14" y="{margin_t + chart_h / 2}" '
        f'fill="var(--axis)" font-size="11" '
        f'transform="rotate(-90 14 {margin_t + chart_h / 2})" '
        f'text-anchor="middle">σ (kcal/mol)</text>'
    )
    # Y-axis ticks (0, ymax/2, ymax).
    for frac, val in ((0.0, 0.0), (0.5, ymax / 2.0), (1.0, ymax)):
        ty = margin_t + chart_h - frac * chart_h
        parts.append(
            f'<line x1="{margin_l - 4}" y1="{ty}" x2="{margin_l}" y2="{ty}" '
            f'stroke="var(--axis)" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{margin_l - 8}" y="{ty + 3}" fill="var(--axis)" '
            f'font-size="10" text-anchor="end">{val:.2f}</text>'
        )
    # Bars.
    for i, (lbl, val, color, missing) in enumerate(numeric):
        x_center = margin_l + slot_w * (i + 0.5)
        x_left = x_center - bar_w / 2
        h_px = (val / ymax) * chart_h if ymax > 0 else 0.0
        if h_px < 0:
            h_px = 0.0
        y_top = margin_t + chart_h - h_px
        if missing:
            parts.append(
                f'<rect x="{x_left:.1f}" y="{margin_t + chart_h - 4}" '
                f'width="{bar_w:.1f}" height="4" '
                f'fill="{color}" fill-opacity="0.25" '
                f'stroke="{color}" stroke-dasharray="3 2"/>'
            )
            parts.append(
                f'<text x="{x_center:.1f}" y="{margin_t + chart_h - 8}" '
                f'fill="var(--axis)" font-size="10" '
                f'text-anchor="middle">n/a</text>'
            )
        else:
            parts.append(
                f'<rect x="{x_left:.1f}" y="{y_top:.1f}" '
                f'width="{bar_w:.1f}" height="{h_px:.1f}" '
                f'fill="{color}"/>'
            )
            parts.append(
                f'<text x="{x_center:.1f}" y="{y_top - 4:.1f}" '
                f'fill="var(--axis)" font-size="10" '
                f'text-anchor="middle">{val:.2f}</text>'
            )
        # X-axis label.
        parts.append(
            f'<text x="{x_center:.1f}" y="{margin_t + chart_h + 16}" '
            f'fill="var(--axis)" font-size="10.5" '
            f'text-anchor="middle">{_html_escape(lbl)}</text>'
        )
    # Title (top).
    parts.append(
        f'<text x="{width / 2:.1f}" y="14" fill="var(--axis)" '
        f'font-size="11" text-anchor="middle">'
        f'σ decomposition (between-seed vs within-seed median)</text>'
    )
    parts.append('</svg>')
    return "".join(parts)


# ---------------------------------------------------------------------------
# Internal helpers — table renderers
# ---------------------------------------------------------------------------

def _tier_class(tier: Optional[str]) -> str:
    """Return CSS modifier class for a tier label, or ``tier-D`` fallback."""
    if not tier or not isinstance(tier, str):
        return "tier-D"
    suffix = tier.strip().split(".")[-1].upper()
    if suffix in {"A", "B", "C", "D"}:
        return f"tier-{suffix}"
    return "tier-D"


def _tier_badge(tier: Optional[str]) -> str:
    """Inline tier badge (``<span class='tier-badge tier-X'>``)."""
    cls = _tier_class(tier)
    return (
        f'<span class="tier-badge {cls}" '
        f'title="{_html_escape(_TIER_LABELS.get(tier or "", ""))}">'
        f"{_html_escape(tier or '—')}</span>"
    )


def _render_summary_table(result: Dict[str, Any]) -> str:
    """Summary stats table: ΔΔG, SE, df, CI95, z_SE, tier."""
    pw = result.get("pairwise") or {}
    ddg = pw.get("ddG")
    se = pw.get("SE")
    df = pw.get("df")
    t_crit = pw.get("t_crit")
    ci = pw.get("CI95")
    z = pw.get("z_SE")
    tier = pw.get("tier") or ""
    note = pw.get("note") or result.get("note")

    rows = [
        ("ΔΔG (kcal/mol)", _format_float(ddg, sig=3)),
        ("SE (kcal/mol)", _format_float(se, sig=3, signed=False)),
        ("df", _html_escape(df) if df is not None else "—"),
        ("t_crit (95%)",
         _format_float(t_crit, sig=3, signed=False) if t_crit is not None else "—"),
        ("CI95", _format_ci(ci)),
        ("z_SE", _format_float(z, sig=2)),
        ("Tier", _tier_badge(tier)),
    ]
    if note:
        rows.append(("Note", _html_escape(note)))

    body = "\n".join(
        f"<tr><th scope=\"row\">{_html_escape(k)}</th><td>{v}</td></tr>"
        for k, v in rows
    )
    return f'<table class="summary-table">\n<tbody>\n{body}\n</tbody>\n</table>'


def _render_branch_meta(branch: Dict[str, Any]) -> str:
    """Per-branch meta block (key/value grid)."""
    items = [
        ("label", _html_escape(branch.get("label", "—"))),
        ("n_seed", _html_escape(branch.get("n_seed", "—"))),
        ("solvent", _html_escape(branch.get("solvent_model", "—"))),
        ("mean_of_seeds",
         _format_float(branch.get("mean_of_seeds"), sig=3)),
        ("σ_btwn",
         _format_float(branch.get("sigma_btwn"), sig=3, signed=False)
         if branch.get("sigma_btwn") is not None else "—"),
        ("σ_w_median",
         _format_float(branch.get("sigma_w_median"), sig=3, signed=False)),
        ("SE", _format_float(branch.get("SE"), sig=3, signed=False)),
        ("charge_audit", _html_escape(branch.get("charge_audit", "—"))),
    ]
    parts = ['<div class="kv-grid">']
    for k, v in items:
        parts.append(
            f'<div class="k">{_html_escape(k)}</div><div class="v">{v}</div>'
        )
    parts.append("</div>")
    return "".join(parts)


def _render_verdict_line(result: Dict[str, Any]) -> str:
    """Natural-language verdict sentence keyed off tier."""
    pw = result.get("pairwise") or {}
    tier = pw.get("tier") or ""
    branches = result.get("branches") or {}
    var_label = _safe_get(branches, "variant", "label", default="variant")
    wt_label = _safe_get(branches, "wt", "label", default="WT")
    ddg = pw.get("ddG")
    ddg_str = _format_float(ddg, sig=2)
    verdict = _TIER_VERDICTS.get(tier, "has an unrecognised verdict tier")
    sign_label = (
        "stronger binding than WT" if isinstance(ddg, (int, float)) and ddg < 0
        else "weaker binding than WT" if isinstance(ddg, (int, float)) and ddg > 0
        else "indistinguishable from WT"
    )
    return (
        f"Variant <code>{_html_escape(var_label)}</code> "
        f"(vs WT <code>{_html_escape(wt_label)}</code>) "
        f"{verdict} (tier <strong>{_html_escape(tier or '—')}</strong>; "
        f"ΔΔG = {ddg_str} kcal/mol; predicted {sign_label})."
    )


# ---------------------------------------------------------------------------
# Internal helpers — multi-comparison
# ---------------------------------------------------------------------------

def _rank_results(results: List[Dict[str, Any]]) -> List[Tuple[int, Dict[str, Any]]]:
    """Stable-rank by ΔΔG (more-negative = higher rank, tier-aware tiebreak).

    Returns list of ``(rank, result)`` pairs. Rank starts at 1.
    """
    def _key(r: Dict[str, Any]) -> Tuple[int, float, str]:
        pw = r.get("pairwise") or {}
        # Tier order: A best, then B, C, D.
        tier_rank_map = {"X.A": 0, "X.B": 1, "X.C": 2, "X.D": 3}
        tier_rank = tier_rank_map.get(pw.get("tier") or "", 4)
        ddg = pw.get("ddG")
        try:
            ddg_f = float(ddg) if ddg is not None else math.inf
        except (TypeError, ValueError):
            ddg_f = math.inf
        label = _safe_get(r, "branches", "variant", "label", default="") or ""
        return (tier_rank, ddg_f, label)
    sorted_r = sorted(results, key=_key)
    return list(enumerate(sorted_r, start=1))


def _render_ranked_table(ranked: List[Tuple[int, Dict[str, Any]]]) -> str:
    """Top-level ranked variant table."""
    header = (
        "<thead><tr><th>Rank</th><th>Variant</th><th>WT</th>"
        "<th>ΔΔG (kcal/mol)</th><th>SE</th><th>CI95</th>"
        "<th>z_SE</th><th>Tier</th></tr></thead>"
    )
    rows: List[str] = []
    for rank, r in ranked:
        pw = r.get("pairwise") or {}
        var_label = _safe_get(r, "branches", "variant", "label", default="—")
        wt_label = _safe_get(r, "branches", "wt", "label", default="—")
        rows.append(
            "<tr>"
            f"<td>{rank}</td>"
            f"<td><code>{_html_escape(var_label)}</code></td>"
            f"<td><code>{_html_escape(wt_label)}</code></td>"
            f"<td>{_format_float(pw.get('ddG'), sig=3)}</td>"
            f"<td>{_format_float(pw.get('SE'), sig=3, signed=False)}</td>"
            f"<td>{_format_ci(pw.get('CI95'))}</td>"
            f"<td>{_format_float(pw.get('z_SE'), sig=2)}</td>"
            f"<td>{_tier_badge(pw.get('tier'))}</td>"
            "</tr>"
        )
    body = "<tbody>\n" + "\n".join(rows) + "\n</tbody>"
    return f"<table>\n{header}\n{body}\n</table>"


def _render_pairwise_matrix(ranked: List[Tuple[int, Dict[str, Any]]]) -> str:
    """N×N pairwise ΔΔG matrix.

    The matrix shows ``ΔΔG_ij = ΔΔG_i − ΔΔG_j`` (each result's variant-vs-WT
    ΔΔG, differenced). Diagonal is zero. Useful for industrial users to read
    "how much better is variant i than variant j (relative to the same WT)".
    """
    if not ranked:
        return '<p class="muted">No variants to compare.</p>'
    n = len(ranked)
    labels: List[str] = []
    ddgs: List[Optional[float]] = []
    for _, r in ranked:
        labels.append(_safe_get(r, "branches", "variant", "label", default="—"))
        try:
            ddgs.append(float(r.get("pairwise", {}).get("ddG")))
        except (TypeError, ValueError):
            ddgs.append(None)
    head_cells = "".join(
        f"<th>{_html_escape(lbl)}</th>" for lbl in labels
    )
    rows: List[str] = []
    for i, lbl_i in enumerate(labels):
        cells = [f"<th scope=\"row\"><code>{_html_escape(lbl_i)}</code></th>"]
        for j in range(n):
            if ddgs[i] is None or ddgs[j] is None:
                cells.append("<td>—</td>")
                continue
            diff = ddgs[i] - ddgs[j]
            cells.append(f"<td>{_format_float(diff, sig=2)}</td>")
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return (
        "<table>\n"
        f"<thead><tr><th>i \\ j</th>{head_cells}</tr></thead>\n"
        "<tbody>\n" + "\n".join(rows) + "\n</tbody>\n"
        "</table>\n"
        '<p style="font-size:12px;color:var(--muted);margin-top:6px;">'
        "Cell value = ΔΔG<sub>i</sub> − ΔΔG<sub>j</sub> (kcal/mol). "
        "More-negative = row variant binds more strongly than column variant."
        "</p>"
    )


def _render_synthesis_priority(ranked: List[Tuple[int, Dict[str, Any]]]) -> str:
    """Synthesis priority list — variants sorted by tier+ΔΔG."""
    if not ranked:
        return '<p class="muted">No variants in priority list.</p>'
    items: List[str] = []
    for rank, r in ranked:
        pw = r.get("pairwise") or {}
        var_label = _safe_get(r, "branches", "variant", "label", default="—")
        tier = pw.get("tier") or ""
        ddg = pw.get("ddG")
        verdict = _TIER_VERDICTS.get(tier, "")
        items.append(
            f'<li><strong>#{rank}</strong> '
            f'<code>{_html_escape(var_label)}</code> '
            f'{_tier_badge(tier)} — ΔΔG = {_format_float(ddg, sig=2)} kcal/mol '
            f'<em>({_html_escape(verdict)})</em></li>'
        )
    return "<ol>" + "".join(items) + "</ol>"


# ---------------------------------------------------------------------------
# Document scaffold
# ---------------------------------------------------------------------------

def _document_skeleton(title: str, theme: str, body: str) -> str:
    """Wrap a body string in the ``<!DOCTYPE>`` document scaffold."""
    css = _inline_css(theme)
    generated = _now_iso()
    return (
        "<!DOCTYPE html>\n"
        f'<html lang="en" data-theme="{_html_escape(theme)}">\n'
        "<head>\n"
        '  <meta charset="utf-8">\n'
        '  <meta name="viewport" content="width=device-width,initial-scale=1">\n'
        f'  <meta name="generator" content="{_html_escape(TOOL_NAME)} '
        f'({_html_escape(SCHEMA_VERSION)})">\n'
        f'  <meta name="generated" content="{_html_escape(generated)}">\n'
        f"  <title>{_html_escape(title)}</title>\n"
        f"  {css}\n"
        "</head>\n"
        '<body>\n<div class="report">\n'
        f"{body}\n"
        '<footer class="report-footer">'
        f"Generated by UPDD v0.7.x · "
        f"{_html_escape(TOOL_NAME)} {_html_escape(SCHEMA_VERSION)} · "
        f"{_html_escape(generated)}"
        "</footer>\n"
        "</div>\n</body>\n</html>\n"
    )


# ---------------------------------------------------------------------------
# Public renderers
# ---------------------------------------------------------------------------

def render_single_comparison(result: Union[Dict[str, Any], Any], *,
                             theme: str = "light",
                             embed_assets: bool = True) -> str:
    """Render a single WT-vs-variant comparison as a standalone HTML document.

    Args:
        result: A ``branched_ddg/0.1`` JSON dict (or any object with a
            ``to_dict()`` method that returns one).
        theme: ``"light"`` (default) or ``"dark"``.
        embed_assets: Day-1 only supports ``True`` (single-file output).
            ``False`` is reserved for Day 2 (separate ``report_assets/``).

    Returns:
        Full HTML document string.

    Raises:
        InputValidationError: on schema mismatch.
        RenderError: on internal render failure.
    """
    if not embed_assets:
        raise RenderError(
            "embed_assets=False is reserved for Day 2; only embedded mode "
            "is supported in html_report/0.1"
        )
    if hasattr(result, "to_dict") and callable(result.to_dict):  # type: ignore[attr-defined]
        result = result.to_dict()  # type: ignore[assignment]
    _validate_branched_ddg_input(result)  # type: ignore[arg-type]
    rd: Dict[str, Any] = result  # type: ignore[assignment]

    pw = rd.get("pairwise") or {}
    branches = rd.get("branches") or {}
    var_branch = branches.get("variant") or {}
    wt_branch = branches.get("wt") or {}
    var_label = var_branch.get("label", "variant")
    wt_label = wt_branch.get("label", "WT")
    target_id = rd.get("target_id") or "(target unset)"
    tier = pw.get("tier") or ""
    tier_cls = _tier_class(tier)
    sign_conv = (rd.get("sign_convention")
                 or pw.get("sign_convention") or "")

    body_parts: List[str] = []

    # Header.
    body_parts.append(
        '<header class="report-header">\n'
        f'  <h1>UPDD ΔΔG Report — {_html_escape(target_id)} '
        f'/ {_html_escape(var_label)} vs {_html_escape(wt_label)}</h1>\n'
        '  <div class="meta">\n'
        f'    Schema: <code>{_html_escape(rd.get("schema", "—"))}</code> '
        f'(<code>{_html_escape(rd.get("tool", "—"))}</code>) · '
        f'Generated: <code>{_html_escape(rd.get("generated_at", "—"))}</code> · '
        f'Solvent: <code>{_html_escape(rd.get("solvent_model", "—"))}</code> · '
        f'df strategy: <code>{_html_escape(rd.get("df_strategy", "—"))}</code>\n'
        '  </div>\n'
        '</header>'
    )

    # Summary section.
    body_parts.append(
        "<section class=\"summary\">\n"
        "  <h2>Summary</h2>\n"
        f"  {_render_summary_table(rd)}\n"
        f"  <p style=\"font-size:12px;color:var(--muted);margin:8px 0 0 0;\">"
        f"Sign convention: {_html_escape(sign_conv)}</p>\n"
        "</section>"
    )

    # σ decomposition.
    body_parts.append(
        '<section class="sigma">\n'
        '  <h2>σ Decomposition</h2>\n'
        f'  {_render_sigma_chart_svg(rd)}\n'
        '  <p style="font-size:12px;color:var(--muted);margin:8px 0 0 0;">'
        'Between-seed dispersion (<code>σ_btwn</code>) vs within-seed '
        'median dispersion (<code>σ_w_median</code>) for each branch. '
        'When <code>σ_btwn ≫ σ_w</code>, between-seed sampling dominates the '
        'standard error; when <code>σ_btwn ≪ σ_w</code>, more within-seed '
        'snapshots are unlikely to reduce SE.</p>\n'
        '</section>'
    )

    # Per-branch meta.
    body_parts.append(
        '<section class="branches">\n'
        '  <h2>Per-Branch Aggregates</h2>\n'
        f'  <h3 style="margin:10px 0 6px 0;font-size:14px;">'
        f'WT — <code>{_html_escape(wt_label)}</code></h3>\n'
        f'  {_render_branch_meta(wt_branch)}\n'
        f'  <h3 style="margin:14px 0 6px 0;font-size:14px;">'
        f'Variant — <code>{_html_escape(var_label)}</code></h3>\n'
        f'  {_render_branch_meta(var_branch)}\n'
        '</section>'
    )

    # Verdict.
    body_parts.append(
        f'<section class="verdict {tier_cls}">\n'
        '  <h2>Verdict</h2>\n'
        f'  <p>{_render_verdict_line(rd)}</p>\n'
        f'  <p style="font-size:12.5px;color:var(--muted);margin-top:6px;">'
        f'{_html_escape(_TIER_LABELS.get(tier, ""))}</p>\n'
        '</section>'
    )

    # Charge audit notice (if FAIL or WARN).
    audit_status = rd.get("charge_audit") or "PASS"
    if audit_status != "PASS":
        body_parts.append(
            '<section class="audit">\n'
            '  <h2>Charge Audit</h2>\n'
            f'  <p>Status: <strong>{_html_escape(audit_status)}</strong>. '
            'Inspect the JSON record for per-manifest details.</p>\n'
            '</section>'
        )

    title = f"UPDD ΔΔG — {target_id} {var_label}"
    return _document_skeleton(title, theme, "\n".join(body_parts))


def render_multi_ranking(results: List[Union[Dict[str, Any], Any]], *,
                         theme: str = "light",
                         embed_assets: bool = True) -> str:
    """Render a multi-variant ranked report.

    Args:
        results: List of ``branched_ddg/0.1`` JSON dicts (or dataclasses
            with ``to_dict()``).
        theme: ``"light"`` (default) or ``"dark"``.
        embed_assets: Day-1 only supports ``True``.

    Returns:
        Full HTML document string.
    """
    if not embed_assets:
        raise RenderError(
            "embed_assets=False is reserved for Day 2"
        )
    if not results:
        raise InputValidationError(
            "results must be a non-empty list"
        )
    norm: List[Dict[str, Any]] = []
    for r in results:
        if hasattr(r, "to_dict") and callable(r.to_dict):
            r = r.to_dict()
        _validate_branched_ddg_input(r)  # type: ignore[arg-type]
        norm.append(r)  # type: ignore[arg-type]

    ranked = _rank_results(norm)
    target_set = sorted({r.get("target_id") or "(target unset)" for r in norm})
    target_str = ", ".join(target_set)

    body: List[str] = []
    body.append(
        '<header class="report-header">\n'
        f'  <h1>UPDD ΔΔG Multi-Variant Ranking — '
        f'{_html_escape(target_str)}</h1>\n'
        '  <div class="meta">\n'
        f'    Schema: <code>{_html_escape(SCHEMA_VERSION)}</code> · '
        f'Generated: <code>{_html_escape(_now_iso())}</code> · '
        f'N variants: {len(norm)}\n'
        '  </div>\n'
        '</header>'
    )
    body.append(
        '<section class="ranking">\n'
        '  <h2>Ranked Variants</h2>\n'
        f'  {_render_ranked_table(ranked)}\n'
        '</section>'
    )
    body.append(
        '<section class="priority">\n'
        '  <h2>Synthesis Priority</h2>\n'
        f'  {_render_synthesis_priority(ranked)}\n'
        '</section>'
    )
    body.append(
        '<section class="pairwise-matrix">\n'
        '  <h2>Pairwise ΔΔG Matrix</h2>\n'
        f'  {_render_pairwise_matrix(ranked)}\n'
        '</section>'
    )
    title = f"UPDD ΔΔG Multi-Ranking — {target_str}"
    return _document_skeleton(title, theme, "\n".join(body))


def render_cycle_progression(history: Union[Dict[str, Any], Any], *,
                             theme: str = "light") -> str:
    """Render a cycle-progression report (Day-2 stub).

    Day-1 implementation: emits a placeholder document with the project ID,
    cycle count, and a TODO note. Full timeline + K_d trajectory chart is
    deferred to Day 2.
    """
    if hasattr(history, "to_dict") and callable(history.to_dict):
        history = history.to_dict()
    _validate_cycle_history_input(history)  # type: ignore[arg-type]
    h: Dict[str, Any] = history  # type: ignore[assignment]

    pid = h.get("project_id") or "(project unset)"
    n_cycles = len(h.get("cycles") or [])
    last_cycle = (h.get("cycles") or [{}])[-1] if n_cycles else {}
    last_n = last_cycle.get("cycle_n", "—") if last_cycle else "—"

    body: List[str] = []
    body.append(
        '<header class="report-header">\n'
        f'  <h1>UPDD Cycle Progression — {_html_escape(pid)}</h1>\n'
        '  <div class="meta">\n'
        f'    Schema: <code>{_html_escape(h.get("schema", "—"))}</code> · '
        f'Cycles: {n_cycles} · Last cycle: <code>{_html_escape(last_n)}</code> '
        f'· Generated: <code>{_html_escape(_now_iso())}</code>\n'
        '  </div>\n'
        '</header>'
    )
    body.append(
        '<section>\n'
        '  <h2>Project Overview</h2>\n'
        '  <div class="kv-grid">\n'
        f'    <div class="k">project_id</div><div class="v">'
        f'{_html_escape(pid)}</div>\n'
        f'    <div class="k">target</div><div class="v">'
        f'{_html_escape(h.get("target_pdb", "—"))}</div>\n'
        f'    <div class="k">created_at</div><div class="v">'
        f'{_html_escape(h.get("created_at", "—"))}</div>\n'
        f'    <div class="k">n_cycles</div><div class="v">{n_cycles}</div>\n'
        '  </div>\n'
        '</section>'
    )
    body.append(
        '<section>\n'
        '  <h2>Day-2 placeholder</h2>\n'
        '  <p>Full cycle timeline (Cycle 0 → N) with K<sub>d</sub> trajectory'
        ' and plateau indicator is scheduled for PR-NEW-E Day 2-3.</p>\n'
        '  <p style="font-size:12.5px;color:var(--muted);">'
        'Inspect the underlying <code>cycle_history.json</code> for per-cycle '
        'records and <code>utils.cycle_manager.compute_trajectory()</code> '
        'for trajectory metadata.</p>\n'
        '</section>'
    )
    title = f"UPDD Cycle Progression — {pid}"
    return _document_skeleton(title, theme, "\n".join(body))


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def write_report(html_text: str, output_path: str) -> str:
    """Atomic write of an HTML string. Returns the absolute path written.

    Uses ``os.replace()`` on a same-directory temp file to avoid partial
    writes if the process is interrupted.
    """
    if not isinstance(html_text, str):
        raise RenderError(
            f"html_text must be str, got {type(html_text).__name__}"
        )
    out = Path(output_path).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(
        prefix=".html_report_", suffix=".tmp",
        dir=str(out.parent),
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            f.write(html_text)
        os.replace(tmp_path, out)
    except Exception:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise
    return str(out)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="python -m utils.html_report",
        description=("UPDD HTML Report Generator (PR-NEW-E Day 1; "
                     f"schema {SCHEMA_VERSION})."),
    )
    ap.add_argument("--result", required=True,
                    help="Path to a branched_ddg/0.1 JSON file (single mode), "
                         "or a JSON file containing a list of such records "
                         "(multi mode), or a cycle_manager/0.1 JSON file "
                         "(cycle mode).")
    ap.add_argument("--output", required=True,
                    help="Output HTML file path (will be created or "
                         "overwritten atomically).")
    ap.add_argument("--mode", choices=("single", "multi", "cycle"),
                    default="single",
                    help="Render mode (default: single).")
    ap.add_argument("--theme", choices=("light", "dark"),
                    default="light",
                    help="CSS theme (default: light).")
    ap.add_argument("--embed-assets", action="store_true", default=True,
                    help="Embed CSS + SVG inline (only mode supported in "
                         "Day 1; flag retained for future compatibility).")
    return ap


def main(argv: Optional[List[str]] = None) -> int:
    ap = _build_parser()
    args = ap.parse_args(argv)
    try:
        payload = _load_input(args.result)
        if args.mode == "single":
            html_text = render_single_comparison(
                payload, theme=args.theme, embed_assets=True
            )
        elif args.mode == "multi":
            if not isinstance(payload, list):
                # Allow {"results": [...]} envelope.
                if isinstance(payload, dict) and isinstance(
                        payload.get("results"), list):
                    payload = payload["results"]
                else:
                    raise InputValidationError(
                        "multi mode expects a JSON list (or {\"results\": "
                        "[...]} envelope)"
                    )
            html_text = render_multi_ranking(
                payload, theme=args.theme, embed_assets=True
            )
        else:  # cycle
            html_text = render_cycle_progression(payload, theme=args.theme)
    except InputValidationError as exc:
        print(f"[html_report] InputValidationError: {exc}", file=sys.stderr)
        return 3
    except RenderError as exc:
        print(f"[html_report] RenderError: {exc}", file=sys.stderr)
        return 2
    out = write_report(html_text, args.output)
    print(f"[html_report] schema={SCHEMA_VERSION}")
    print(f"[html_report] mode={args.mode}  theme={args.theme}")
    print(f"[html_report] wrote: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
