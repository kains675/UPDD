#!/usr/bin/env python
"""scripts/updd_cli.py — UPDD User CLI (Capability Level 1 demo integration).

Schema: ``updd_cli/0.1`` (PR-NEW-D Day 1, 2026-04-28)

Unified subcommand interface that integrates the v0.7 Stage 1-5 user-facing
layer (per ``v07_user_verdict/v0.7 Section 4 finalized.md`` §4.1.1 — §4.1.3,
§7.0):

    updd check      Single sequence binder verification
    updd compare    Multi-sequence pairwise ranking
    updd refine     Iterative refinement (Cycle 1+, WT vs variants)
    updd discover   De novo discovery (Cycle 0, file-based input — Day 1)
    updd report     Generate report from existing branched_ddg results

Integrations:
    - :mod:`utils.sequence_parser`  (PR-NEW-A) — input layer
    - :mod:`utils.branched_ddg`     (PR-NEW-B) — ΔΔG aggregation engine

Day-1 scope (this revision):
    - argparse subcommand dispatch + CLIContext propagation
    - Sequence parsing (delegates to PR-NEW-A)
    - Branched ΔΔG fixture-driven `report` (delegates to PR-NEW-B)
    - **Mock pipeline** for `check` / `compare` / `refine` / `discover`
      (deterministic placeholder ΔG estimates; clearly labelled as
      ``"mock": true``). Real MD/PBSA hookup is Day 2-3.

Note on namespace: ``utils/updd_cli.py`` already exists as a separate legacy
TUI/menu module driving the existing ``UPDD.py`` orchestrator. **This file
deliberately lives at ``scripts/updd_cli.py``** and does not import or modify
the legacy module. There is no functional overlap.

CLI usage::

    python scripts/updd_cli.py check    --target /path/to/target.pdb \\
        --sequence "ICVVQDWGHHRCT" [--dry-run]
    python scripts/updd_cli.py compare  --target /path/to/target.pdb \\
        --sequences "ICVVQDWGHHRCT,ICVVQ(MTR)WGHHRCT"
    python scripts/updd_cli.py refine   --target /path/to/target.pdb \\
        --wt "ICVVQDWGHHRCT" \\
        --variants "ICVVQ(MTR)WGHHRCT,ICVVQ(NMA)WGHHRCT"
    python scripts/updd_cli.py discover --target /path/to/target.pdb \\
        --sequence-file candidates.txt
    python scripts/updd_cli.py report   --result branched_ddg.json
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# sys.path / package import shim — keep this script runnable directly via
# ``python scripts/updd_cli.py ...`` AND via ``python -m scripts.updd_cli``.
# ---------------------------------------------------------------------------

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from utils.sequence_parser import (  # noqa: E402
    NcAAValidationError,
    SequenceFormatError,
    SequenceParseResult,
    parse as seq_parse,
)
from utils import branched_ddg as bddg  # noqa: E402


# ---------------------------------------------------------------------------
# Schema / constants
# ---------------------------------------------------------------------------

SCHEMA_VERSION = "updd_cli/0.1"
TOOL_NAME = "scripts.updd_cli"
DEFAULT_OUTPUT_PREFIX = "outputs/cli_run"

# Subcommand name for diagnostic / output payload tagging.
_SUBCOMMANDS = ("check", "compare", "refine", "discover", "report")

# Tier classification thresholds when synthesising mock ΔΔG values
# (matches PR-NEW-B X.A / X.B / X.C / X.D semantics).
_MOCK_DDG_PER_AA = -0.5      # base per-residue contribution (ranking-only)
_MOCK_NCAA_BONUS = -1.5      # added for each ncAA position (Capability L1 hint)
_MOCK_BASELINE = -10.0       # WT-like reference baseline
_MOCK_SE = 1.5               # synthetic SE (placeholder)

# Verdict thresholds for `updd check` (Day-1 mock):
_CHECK_DG_BINDER_THRESHOLD = -8.0     # ΔG < -8 ⇒ "yes"
_CHECK_DG_NONBINDER_THRESHOLD = -5.0  # ΔG > -5 ⇒ "no"


# ---------------------------------------------------------------------------
# CLI context
# ---------------------------------------------------------------------------

@dataclass
class CLIContext:
    """Common options propagated to every subcommand handler.

    Constructed from the parsed ``argparse.Namespace`` via
    :meth:`from_args`. Subcommand handlers receive both the raw
    ``args`` and this typed view.
    """

    output_dir: Path
    solvent: str               # "auto" | "pbsa" | "gbsa"
    df_strategy: str           # "conservative" | "min"
    strict_charge_audit: bool
    verbose: bool
    quiet: bool
    dry_run: bool
    fmt: str                   # "json" | "markdown" | "both"
    target: Optional[Path]
    timestamp: str             # ISO-8601 UTC, used in filenames
    extras: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> "CLIContext":
        ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        out = getattr(args, "output", None)
        if out:
            output_dir = Path(out)
        else:
            output_dir = Path(f"{DEFAULT_OUTPUT_PREFIX}_{ts}")
        target = getattr(args, "target", None)
        return cls(
            output_dir=output_dir,
            solvent=getattr(args, "solvent", "auto"),
            df_strategy=getattr(args, "df_strategy", "conservative"),
            strict_charge_audit=getattr(args, "strict_charge_audit", False),
            verbose=getattr(args, "verbose", False),
            quiet=getattr(args, "quiet", False),
            dry_run=getattr(args, "dry_run", False),
            fmt=getattr(args, "format", "both"),
            target=Path(target) if target else None,
            timestamp=ts,
        )

    def setup_logging(self) -> logging.Logger:
        level = logging.INFO
        if self.verbose:
            level = logging.DEBUG
        elif self.quiet:
            level = logging.WARNING
        logger = logging.getLogger(TOOL_NAME)
        logger.setLevel(level)
        # Idempotent handler installation — re-runs in the same process
        # (e.g. test harness) must not multiply log lines.
        if not logger.handlers:
            h = logging.StreamHandler(sys.stderr)
            h.setLevel(level)
            h.setFormatter(logging.Formatter("[updd_cli] %(levelname)s: %(message)s"))
            logger.addHandler(h)
        else:
            for h in logger.handlers:
                h.setLevel(level)
        return logger


# ---------------------------------------------------------------------------
# Provenance payload helper
# ---------------------------------------------------------------------------

def _provenance_header(cmd: str, ctx: CLIContext) -> Dict[str, Any]:
    """Common JSON-output preamble: schema, tool, timestamp, subcommand."""
    return {
        "schema": SCHEMA_VERSION,
        "tool": TOOL_NAME,
        "subcommand": cmd,
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "options": {
            "target": str(ctx.target) if ctx.target else None,
            "solvent": ctx.solvent,
            "df_strategy": ctx.df_strategy,
            "strict_charge_audit": ctx.strict_charge_audit,
            "dry_run": ctx.dry_run,
            "format": ctx.fmt,
        },
    }


# ---------------------------------------------------------------------------
# Mock pipeline (Day 1 only) — deterministic placeholder
# ---------------------------------------------------------------------------

def _maybe_mock_pipeline(parsed: SequenceParseResult, ctx: CLIContext,
                         *, label: str = "candidate") -> Dict[str, Any]:
    """Day-1 mock ΔG estimate.

    Real MD/PBSA execution is not invoked from this CLI in the Day-1 scope
    (no GPU + no input PDB pipeline integration yet). The return value is
    a deterministic synthetic estimate purely for CLI flow validation.

    Determinism: ``ΔG_mock = -10.0 + (-0.5) × length + (-1.5) × n_ncaa``
    (matches the convention noted in the dispatch hint, with an explicit
    ncAA bonus making the ranking responsive to Capability Level 1's
    "ncAA differentiator" framing per §7.0).

    Returns a dict with keys::

        {
            "label": str,
            "dG_estimate": float,            # kcal/mol
            "dG_se": float,                  # placeholder SE
            "n_residues": int,
            "n_ncaa": int,
            "mock": True,
            "status": "mock" | "dry_run",
            "would_run": bool,
        }
    """
    n_res = parsed.length
    n_ncaa = len(parsed.ncaa_positions)
    if ctx.dry_run:
        return {
            "label": label,
            "status": "dry_run",
            "would_run": True,
            "n_residues": n_res,
            "n_ncaa": n_ncaa,
            "mock": True,
            "dG_estimate": None,
            "dG_se": None,
        }
    dG = _MOCK_BASELINE + _MOCK_DDG_PER_AA * n_res + _MOCK_NCAA_BONUS * n_ncaa
    return {
        "label": label,
        "status": "mock",
        "would_run": False,
        "n_residues": n_res,
        "n_ncaa": n_ncaa,
        "mock": True,
        "dG_estimate": float(dG),
        "dG_se": float(_MOCK_SE),
    }


def _classify_binder(dG: Optional[float]) -> Tuple[str, float]:
    """Return ``(verdict, confidence)`` for an estimated ΔG.

    Verdicts: ``"yes"`` / ``"no"`` / ``"undetermined"``. Confidence is a
    heuristic 0-1 score derived from the distance to the nearest decision
    boundary (Day-1 placeholder).
    """
    if dG is None:
        return "undetermined", 0.0
    if dG <= _CHECK_DG_BINDER_THRESHOLD:
        c = min(1.0, max(0.0, (_CHECK_DG_BINDER_THRESHOLD - dG) / 5.0 + 0.5))
        return "yes", round(c, 3)
    if dG >= _CHECK_DG_NONBINDER_THRESHOLD:
        c = min(1.0, max(0.0, (dG - _CHECK_DG_NONBINDER_THRESHOLD) / 5.0 + 0.5))
        return "no", round(c, 3)
    return "undetermined", 0.3


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _write_outputs(payload: Dict[str, Any], md_text: str,
                   *, ctx: CLIContext, basename: str
                   ) -> Tuple[Optional[Path], Optional[Path]]:
    """Write JSON + markdown per ``ctx.fmt``. Returns ``(json_path, md_path)``.

    Either path may be ``None`` if the format excludes it. Creates
    ``ctx.output_dir`` if needed.
    """
    json_path: Optional[Path] = None
    md_path: Optional[Path] = None
    if ctx.fmt in ("json", "both"):
        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        json_path = ctx.output_dir / f"{basename}.json"
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, default=_json_default)
    if ctx.fmt in ("markdown", "both"):
        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        md_path = ctx.output_dir / f"{basename}.md"
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(md_text)
    return json_path, md_path


def _json_default(obj: Any) -> Any:
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, tuple):
        return list(obj)
    if hasattr(obj, "to_dict"):
        return obj.to_dict()
    if hasattr(obj, "__dict__"):
        return obj.__dict__
    return str(obj)


def _safe_parse(seq: str, *, label: str, logger: logging.Logger
                ) -> SequenceParseResult:
    """Parse a sequence and re-raise with a user-friendly error chain.

    The parser raises :class:`NcAAValidationError` / :class:`SequenceFormatError`
    on bad input; we let them propagate so the dispatcher can format the
    error consistently across subcommands.
    """
    logger.debug(f"parsing sequence ({label!r}): {seq!r}")
    return seq_parse(seq)


def _check_target(ctx: CLIContext, *, required: bool, logger: logging.Logger) -> None:
    """Validate target path. Soft-fails under ``--dry-run`` (allowed).

    When ``required=True`` and the path is missing on disk, log a warning
    and (for non-dry-run) raise :class:`FileNotFoundError`.
    """
    if not ctx.target:
        if required:
            raise SystemExit(
                "[updd_cli] --target is required for this subcommand"
            )
        return
    if not ctx.target.exists():
        if ctx.dry_run:
            logger.warning(
                f"target path not found: {ctx.target}; permitted under --dry-run"
            )
        else:
            raise FileNotFoundError(f"target file not found: {ctx.target}")


# ---------------------------------------------------------------------------
# `updd check` — single sequence binder verdict
# ---------------------------------------------------------------------------

def cmd_check(args: argparse.Namespace, ctx: CLIContext) -> int:
    logger = ctx.setup_logging()
    _check_target(ctx, required=True, logger=logger)
    parsed = _safe_parse(args.sequence, label="check", logger=logger)
    pipe = _maybe_mock_pipeline(parsed, ctx, label="candidate")
    verdict, confidence = _classify_binder(pipe.get("dG_estimate"))

    payload: Dict[str, Any] = _provenance_header("check", ctx)
    payload["sequence"] = args.sequence
    payload["parsed"] = parsed.to_dict()
    payload["pipeline"] = pipe
    payload["dG_estimate"] = pipe.get("dG_estimate")
    payload["dG_se"] = pipe.get("dG_se")
    payload["binder_verdict"] = verdict
    payload["confidence"] = confidence

    md = _render_md_check(payload)
    json_path, md_path = _write_outputs(payload, md, ctx=ctx, basename="check")
    _print_summary_check(payload, json_path, md_path, ctx=ctx)
    return 0


def _render_md_check(payload: Dict[str, Any]) -> str:
    p = payload
    parsed = p["parsed"]
    pipe = p["pipeline"]
    lines: List[str] = []
    lines.append("# UPDD `check` Report")
    lines.append("")
    lines.append(f"- Schema: `{p['schema']}` ({p['tool']})")
    lines.append(f"- Generated: {p['generated_at']}")
    lines.append(f"- Subcommand: `check`")
    lines.append(f"- Target: `{p['options']['target']}`")
    lines.append(f"- Solvent: `{p['options']['solvent']}`  "
                 f"df_strategy: `{p['options']['df_strategy']}`")
    lines.append("")
    lines.append("## Input")
    lines.append("")
    lines.append(f"- Sequence: `{p['sequence']}`")
    lines.append(f"- Length: {parsed['length']}  "
                 f"ncAA count: {len(parsed['ncaa_positions'])}  "
                 f"q≈ {parsed['estimated_total_charge']:+d}")
    lines.append(f"- Format detected: `{parsed['format_detected']}`")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append("| Quantity | Value |")
    lines.append("|---|---|")
    dG = pipe.get("dG_estimate")
    dG_str = f"{dG:+.2f}" if isinstance(dG, (int, float)) else "—"
    lines.append(f"| ΔG estimate (kcal/mol) | {dG_str} |")
    se = pipe.get("dG_se")
    se_str = f"{se:.2f}" if isinstance(se, (int, float)) else "—"
    lines.append(f"| ΔG SE | {se_str} |")
    lines.append(f"| Binder verdict | **{p['binder_verdict']}** |")
    lines.append(f"| Confidence | {p['confidence']} |")
    lines.append(f"| Pipeline status | `{pipe.get('status')}` |")
    if pipe.get("mock"):
        lines.append("")
        lines.append("> **Note**: Day-1 mock pipeline. Real MD/PBSA "
                     "integration is Day 2-3.")
    return "\n".join(lines) + "\n"


def _print_summary_check(payload: Dict[str, Any], json_path: Optional[Path],
                         md_path: Optional[Path], *, ctx: CLIContext) -> None:
    if ctx.quiet:
        return
    pipe = payload["pipeline"]
    dG = pipe.get("dG_estimate")
    dG_str = f"{dG:+.3f}" if isinstance(dG, (int, float)) else "—"
    print(f"[updd_cli] check  schema={payload['schema']}  "
          f"verdict={payload['binder_verdict']}  "
          f"ΔG={dG_str}  status={pipe.get('status')}")
    if json_path:
        print(f"[updd_cli]   wrote: {json_path}")
    if md_path:
        print(f"[updd_cli]   wrote: {md_path}")


# ---------------------------------------------------------------------------
# `updd compare` — multi-sequence pairwise ranking
# ---------------------------------------------------------------------------

def cmd_compare(args: argparse.Namespace, ctx: CLIContext) -> int:
    logger = ctx.setup_logging()
    _check_target(ctx, required=True, logger=logger)
    seqs = _split_csv(args.sequences)
    if len(seqs) < 2:
        raise SystemExit(
            "[updd_cli] compare requires at least 2 sequences (got "
            f"{len(seqs)})"
        )
    entries: List[Dict[str, Any]] = []
    for i, s in enumerate(seqs, 1):
        parsed = _safe_parse(s, label=f"compare:{i}", logger=logger)
        pipe = _maybe_mock_pipeline(parsed, ctx, label=f"seq_{i}")
        entries.append({
            "index": i,
            "sequence": s,
            "parsed": parsed.to_dict(),
            "pipeline": pipe,
        })

    # Rank by ΔG ascending (more-negative is more favourable). Dry-run
    # entries (dG=None) are emitted in input order at the end.
    def _rank_key(e: Dict[str, Any]) -> Tuple[int, float]:
        dG = e["pipeline"].get("dG_estimate")
        return (0, dG) if isinstance(dG, (int, float)) else (1, 0.0)

    ranked = sorted(entries, key=_rank_key)
    for new_idx, e in enumerate(ranked, 1):
        e["rank"] = new_idx

    # Pairwise ΔΔG matrix vs the best (rank 1) candidate, plus tier proxy.
    pairwise: List[Dict[str, Any]] = []
    if ranked and ranked[0]["pipeline"].get("dG_estimate") is not None:
        ref = ranked[0]
        for e in ranked[1:]:
            ref_dG = ref["pipeline"]["dG_estimate"]
            this_dG = e["pipeline"].get("dG_estimate")
            if this_dG is None:
                continue
            ddG = this_dG - ref_dG
            se = (ref["pipeline"]["dG_se"] ** 2 + e["pipeline"]["dG_se"] ** 2) ** 0.5
            z = ddG / se if se > 0 else 0.0
            tier = _mock_tier_from_z(z)
            pairwise.append({
                "ref_index": ref["index"],
                "test_index": e["index"],
                "ddG": ddG,
                "SE_combined": se,
                "z_SE": z,
                "tier": tier,
            })

    payload: Dict[str, Any] = _provenance_header("compare", ctx)
    payload["sequences"] = seqs
    payload["ranked_entries"] = ranked
    payload["pairwise_vs_top"] = pairwise

    md = _render_md_compare(payload)
    json_path, md_path = _write_outputs(payload, md, ctx=ctx, basename="compare")
    _print_summary_compare(payload, json_path, md_path, ctx=ctx)
    return 0


def _mock_tier_from_z(z: float) -> str:
    """Day-1 mock tier from a synthetic z-score, mirroring PR-NEW-B
    semantics: |z|≥2 ⇒ X.B (sign-significant), else X.C (insufficient)."""
    if abs(z) >= 2.0:
        return "X.B"
    return "X.C"


def _render_md_compare(payload: Dict[str, Any]) -> str:
    p = payload
    lines: List[str] = []
    lines.append("# UPDD `compare` Report")
    lines.append("")
    lines.append(f"- Schema: `{p['schema']}` ({p['tool']})")
    lines.append(f"- Generated: {p['generated_at']}")
    lines.append(f"- Target: `{p['options']['target']}`")
    lines.append(f"- N sequences: {len(p['sequences'])}")
    lines.append("")
    lines.append("## Ranked Candidates")
    lines.append("")
    lines.append("| Rank | Index | Sequence | Length | nNcAA | ΔG (kcal/mol) | Status |")
    lines.append("|---|---|---|---|---|---|---|")
    for e in p["ranked_entries"]:
        pipe = e["pipeline"]
        dG = pipe.get("dG_estimate")
        dG_str = f"{dG:+.2f}" if isinstance(dG, (int, float)) else "—"
        lines.append(
            f"| {e['rank']} | {e['index']} | `{e['sequence']}` | "
            f"{e['parsed']['length']} | {len(e['parsed']['ncaa_positions'])} | "
            f"{dG_str} | {pipe.get('status')} |"
        )
    lines.append("")
    if p.get("pairwise_vs_top"):
        lines.append("## Pairwise ΔΔG vs Top")
        lines.append("")
        lines.append("| ref | test | ΔΔG | SE | z_SE | tier |")
        lines.append("|---|---|---|---|---|---|")
        for pw in p["pairwise_vs_top"]:
            lines.append(
                f"| seq_{pw['ref_index']} | seq_{pw['test_index']} | "
                f"{pw['ddG']:+.3f} | {pw['SE_combined']:.3f} | "
                f"{pw['z_SE']:+.2f} | {pw['tier']} |"
            )
        lines.append("")
    return "\n".join(lines) + "\n"


def _print_summary_compare(payload: Dict[str, Any], json_path: Optional[Path],
                           md_path: Optional[Path], *, ctx: CLIContext) -> None:
    if ctx.quiet:
        return
    print(f"[updd_cli] compare  schema={payload['schema']}  "
          f"n={len(payload['sequences'])}  "
          f"top_rank=seq_{payload['ranked_entries'][0]['index']}")
    for e in payload["ranked_entries"]:
        pipe = e["pipeline"]
        dG = pipe.get("dG_estimate")
        dG_str = f"{dG:+.3f}" if isinstance(dG, (int, float)) else "—"
        print(f"[updd_cli]   rank {e['rank']}: seq_{e['index']}  "
              f"ΔG={dG_str}  status={pipe.get('status')}")
    if json_path:
        print(f"[updd_cli]   wrote: {json_path}")
    if md_path:
        print(f"[updd_cli]   wrote: {md_path}")


# ---------------------------------------------------------------------------
# `updd refine` — Cycle 1+ iterative refinement (WT + variants)
# ---------------------------------------------------------------------------

def cmd_refine(args: argparse.Namespace, ctx: CLIContext) -> int:
    logger = ctx.setup_logging()
    _check_target(ctx, required=True, logger=logger)
    wt_seq = args.wt
    variants = _split_csv(args.variants)
    if not variants:
        raise SystemExit("[updd_cli] refine requires at least 1 variant")

    wt_parsed = _safe_parse(wt_seq, label="refine:wt", logger=logger)
    wt_pipe = _maybe_mock_pipeline(wt_parsed, ctx, label="WT")

    variant_entries: List[Dict[str, Any]] = []
    for i, v in enumerate(variants, 1):
        v_parsed = _safe_parse(v, label=f"refine:v{i}", logger=logger)
        v_pipe = _maybe_mock_pipeline(v_parsed, ctx, label=f"variant_{i}")
        ddG, se, z, tier, ci_lo, ci_hi = _ddg_from_two_branches(
            wt_pipe, v_pipe
        )
        variant_entries.append({
            "index": i,
            "sequence": v,
            "parsed": v_parsed.to_dict(),
            "pipeline": v_pipe,
            "ddG": ddG,
            "SE_combined": se,
            "z_SE": z,
            "tier": tier,
            "CI95": [ci_lo, ci_hi] if ddG is not None else None,
        })

    # Tier-then-ΔΔG ranking — mirrors PR-NEW-B X.A → X.B → X.C → X.D order.
    _tier_rank = {"X.A": 0, "X.B": 1, "X.C": 2, "X.D": 3, "—": 4}

    def _rank_key(e: Dict[str, Any]) -> Tuple[int, float]:
        return (
            _tier_rank.get(e.get("tier") or "—", 4),
            e["ddG"] if isinstance(e.get("ddG"), (int, float)) else 0.0,
        )

    ranked = sorted(variant_entries, key=_rank_key)
    for n, e in enumerate(ranked, 1):
        e["rank"] = n

    payload: Dict[str, Any] = _provenance_header("refine", ctx)
    payload["wt"] = {
        "sequence": wt_seq,
        "parsed": wt_parsed.to_dict(),
        "pipeline": wt_pipe,
    }
    payload["variants"] = ranked
    payload["compute_branched_ddg_compatible"] = True
    payload["sign_convention"] = "ddG = mean(variant) − mean(WT)"

    md = _render_md_refine(payload)
    json_path, md_path = _write_outputs(payload, md, ctx=ctx, basename="refine")
    _print_summary_refine(payload, json_path, md_path, ctx=ctx)
    return 0


def _ddg_from_two_branches(wt_pipe: Dict[str, Any], v_pipe: Dict[str, Any]
                           ) -> Tuple[Optional[float], Optional[float],
                                      Optional[float], Optional[str],
                                      Optional[float], Optional[float]]:
    """Compute mock ΔΔG / SE / z / tier / CI95 (Day-1 placeholder).

    For dry-run inputs (``dG_estimate is None``) returns all-None tuple.
    Otherwise:

        ddG  = v.dG − wt.dG
        SE   = sqrt(v.SE^2 + wt.SE^2)
        z_SE = ddG / SE  (signed)
        tier = mock_tier_from_z(z) when |z|<2 → 'X.C', else 'X.B'
               (CI95 inclusion uses 1.96·SE — Day-1 large-sample proxy)
    """
    wt_dG = wt_pipe.get("dG_estimate")
    v_dG = v_pipe.get("dG_estimate")
    if wt_dG is None or v_dG is None:
        return None, None, None, None, None, None
    ddG = v_dG - wt_dG
    se = (wt_pipe["dG_se"] ** 2 + v_pipe["dG_se"] ** 2) ** 0.5
    z = ddG / se if se > 0 else 0.0
    half = 1.96 * se
    ci_lo = ddG - half
    ci_hi = ddG + half
    excludes_zero = (ci_lo > 0) or (ci_hi < 0)
    if excludes_zero:
        tier = "X.A"
    elif abs(z) >= 2.0:
        tier = "X.B"
    else:
        tier = "X.C"
    return ddG, se, z, tier, ci_lo, ci_hi


def _render_md_refine(payload: Dict[str, Any]) -> str:
    p = payload
    wt = p["wt"]
    lines: List[str] = []
    lines.append("# UPDD `refine` Report (Cycle 1+)")
    lines.append("")
    lines.append(f"- Schema: `{p['schema']}` ({p['tool']})")
    lines.append(f"- Generated: {p['generated_at']}")
    lines.append(f"- Target: `{p['options']['target']}`")
    lines.append(f"- Sign convention: {p['sign_convention']}")
    lines.append("")
    lines.append("## WT (Reference)")
    lines.append("")
    lines.append(f"- Sequence: `{wt['sequence']}`")
    lines.append(f"- ΔG: "
                 + (f"{wt['pipeline']['dG_estimate']:+.2f}"
                    if isinstance(wt['pipeline'].get('dG_estimate'),
                                  (int, float)) else "—")
                 + f" ± {wt['pipeline'].get('dG_se') or '—'}")
    lines.append("")
    lines.append("## Variants (ranked)")
    lines.append("")
    lines.append("| Rank | Variant | ΔΔG | SE | z_SE | tier | CI95 |")
    lines.append("|---|---|---|---|---|---|---|")
    for v in p["variants"]:
        ddG = v.get("ddG")
        se = v.get("SE_combined")
        z = v.get("z_SE")
        tier = v.get("tier") or "—"
        ci = v.get("CI95")
        ddG_str = f"{ddG:+.3f}" if isinstance(ddG, (int, float)) else "—"
        se_str = f"{se:.3f}" if isinstance(se, (int, float)) else "—"
        z_str = f"{z:+.2f}" if isinstance(z, (int, float)) else "—"
        ci_str = (f"[{ci[0]:+.2f}, {ci[1]:+.2f}]"
                  if ci and isinstance(ci[0], (int, float)) else "—")
        lines.append(
            f"| {v['rank']} | `{v['sequence']}` | {ddG_str} | {se_str} | "
            f"{z_str} | {tier} | {ci_str} |"
        )
    lines.append("")
    lines.append("> Day-1 mock ΔΔG. Real `compute_branched_ddg` integration "
                 "(MD + PBSA per-seed summaries) is Day 2-3.")
    return "\n".join(lines) + "\n"


def _print_summary_refine(payload: Dict[str, Any], json_path: Optional[Path],
                          md_path: Optional[Path], *, ctx: CLIContext) -> None:
    if ctx.quiet:
        return
    n = len(payload["variants"])
    print(f"[updd_cli] refine  schema={payload['schema']}  "
          f"n_variants={n}  WT={payload['wt']['sequence']!r}")
    for v in payload["variants"]:
        ddG = v.get("ddG")
        ddG_str = f"{ddG:+.3f}" if isinstance(ddG, (int, float)) else "—"
        print(f"[updd_cli]   rank {v['rank']}: {v['sequence']!r}  "
              f"ΔΔG={ddG_str}  tier={v.get('tier') or '—'}")
    if json_path:
        print(f"[updd_cli]   wrote: {json_path}")
    if md_path:
        print(f"[updd_cli]   wrote: {md_path}")


# ---------------------------------------------------------------------------
# `updd discover` — Cycle 0 file-based input (Day-1 scope)
# ---------------------------------------------------------------------------

def cmd_discover(args: argparse.Namespace, ctx: CLIContext) -> int:
    logger = ctx.setup_logging()
    _check_target(ctx, required=True, logger=logger)
    src = Path(args.sequence_file)
    if not src.exists():
        raise FileNotFoundError(f"sequence file not found: {src}")
    raw_lines = [ln.strip() for ln in src.read_text(encoding="utf-8").splitlines()]
    seqs = [ln for ln in raw_lines if ln and not ln.startswith("#")]
    if not seqs:
        raise SystemExit(f"[updd_cli] sequence file is empty: {src}")

    candidates: List[Dict[str, Any]] = []
    for i, s in enumerate(seqs, 1):
        try:
            parsed = _safe_parse(s, label=f"discover:{i}", logger=logger)
        except (NcAAValidationError, SequenceFormatError) as exc:
            candidates.append({
                "index": i,
                "sequence": s,
                "status": "parse_error",
                "error": f"{type(exc).__name__}: {exc}",
            })
            continue
        pipe = _maybe_mock_pipeline(parsed, ctx, label=f"cand_{i}")
        candidates.append({
            "index": i,
            "sequence": s,
            "parsed": parsed.to_dict(),
            "pipeline": pipe,
            "status": pipe.get("status"),
        })

    # Sort by ΔG (best first); parse_error / dry_run rows go to bottom.
    def _rank_key(e: Dict[str, Any]) -> Tuple[int, float]:
        if e.get("status") in ("parse_error", "dry_run"):
            return (1, 0.0)
        dG = e.get("pipeline", {}).get("dG_estimate")
        if isinstance(dG, (int, float)):
            return (0, dG)
        return (1, 0.0)

    candidates.sort(key=_rank_key)
    for new_idx, c in enumerate(candidates, 1):
        c["rank"] = new_idx

    payload: Dict[str, Any] = _provenance_header("discover", ctx)
    payload["source"] = str(src.resolve())
    payload["n_input"] = len(seqs)
    payload["ncaa_budget"] = getattr(args, "ncaa_budget", None)
    payload["candidates"] = candidates
    payload["note"] = (
        "Day-1: file-based input only. Stage 1-2 (RFdiffusion + ProteinMPNN + AF2) "
        "integration is Day 2-3."
    )

    md = _render_md_discover(payload)
    json_path, md_path = _write_outputs(payload, md, ctx=ctx, basename="discover")
    _print_summary_discover(payload, json_path, md_path, ctx=ctx)
    return 0


def _render_md_discover(payload: Dict[str, Any]) -> str:
    p = payload
    lines: List[str] = []
    lines.append("# UPDD `discover` Report (Cycle 0, Day-1 scope)")
    lines.append("")
    lines.append(f"- Schema: `{p['schema']}` ({p['tool']})")
    lines.append(f"- Generated: {p['generated_at']}")
    lines.append(f"- Target: `{p['options']['target']}`")
    lines.append(f"- Source: `{p['source']}`  N input: {p['n_input']}")
    if p.get("ncaa_budget") is not None:
        lines.append(f"- ncAA budget: {p['ncaa_budget']}")
    lines.append("")
    lines.append("| Rank | Index | Sequence | Status | ΔG (kcal/mol) | Note |")
    lines.append("|---|---|---|---|---|---|")
    for c in p["candidates"]:
        if c.get("status") == "parse_error":
            lines.append(
                f"| {c['rank']} | {c['index']} | `{c['sequence']}` | "
                f"parse_error | — | {c.get('error', '')} |"
            )
            continue
        pipe = c.get("pipeline", {})
        dG = pipe.get("dG_estimate")
        dG_str = f"{dG:+.2f}" if isinstance(dG, (int, float)) else "—"
        lines.append(
            f"| {c['rank']} | {c['index']} | `{c['sequence']}` | "
            f"{c.get('status')} | {dG_str} |  |"
        )
    lines.append("")
    lines.append(f"> {p['note']}")
    return "\n".join(lines) + "\n"


def _print_summary_discover(payload: Dict[str, Any], json_path: Optional[Path],
                            md_path: Optional[Path], *, ctx: CLIContext) -> None:
    if ctx.quiet:
        return
    n_ok = sum(1 for c in payload["candidates"]
               if c.get("status") not in ("parse_error",))
    n_err = sum(1 for c in payload["candidates"]
                if c.get("status") == "parse_error")
    print(f"[updd_cli] discover  schema={payload['schema']}  "
          f"n_input={payload['n_input']}  parsed_ok={n_ok}  parse_errors={n_err}")
    if json_path:
        print(f"[updd_cli]   wrote: {json_path}")
    if md_path:
        print(f"[updd_cli]   wrote: {md_path}")


# ---------------------------------------------------------------------------
# `updd report` — render existing branched_ddg result
# ---------------------------------------------------------------------------

def cmd_report(args: argparse.Namespace, ctx: CLIContext) -> int:
    logger = ctx.setup_logging()
    src = Path(args.result)
    if not src.exists():
        raise FileNotFoundError(f"result file not found: {src}")
    with open(src, "r", encoding="utf-8") as f:
        result_json = json.load(f)

    fmt = getattr(args, "format", None) or ctx.fmt
    # Normalise; `report` allows html as a future hook.
    if fmt == "html":
        logger.warning("HTML output is reserved for PR-NEW-E; falling back to markdown")
        fmt = "markdown"

    md = _render_md_report(result_json)
    payload: Dict[str, Any] = _provenance_header("report", ctx)
    payload["source"] = str(src.resolve())
    payload["result"] = result_json
    payload["rendered_format"] = fmt

    json_path: Optional[Path] = None
    md_path: Optional[Path] = None
    if fmt in ("json", "both"):
        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        json_path = ctx.output_dir / "report.json"
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2, default=_json_default)
    if fmt in ("markdown", "both"):
        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        md_path = ctx.output_dir / "report.md"
        with open(md_path, "w", encoding="utf-8") as f:
            f.write(md)

    if not ctx.quiet:
        schema = result_json.get("schema", "unknown")
        ddg = result_json.get("pairwise", {}).get("ddG")
        tier = result_json.get("pairwise", {}).get("tier")
        ddg_str = f"{ddg:+.3f}" if isinstance(ddg, (int, float)) else "—"
        print(f"[updd_cli] report  schema={SCHEMA_VERSION}  "
              f"input_schema={schema}  ΔΔG={ddg_str}  tier={tier}")
        if json_path:
            print(f"[updd_cli]   wrote: {json_path}")
        if md_path:
            print(f"[updd_cli]   wrote: {md_path}")
    return 0


def _render_md_report(rj: Dict[str, Any]) -> str:
    """Render a markdown summary from a ``branched_ddg/0.1``-shaped JSON.

    Tolerant of partial / future-schema results (missing keys → "—").
    """
    schema = rj.get("schema", "unknown")
    target = rj.get("target_id") or "(target unset)"
    pw = rj.get("pairwise", {})
    branches = rj.get("branches", {})
    wt = branches.get("wt", {})
    var = branches.get("variant", {})
    lines: List[str] = []
    lines.append(f"# UPDD `report` — Branched ΔΔG Summary")
    lines.append("")
    lines.append(f"- Source schema: `{schema}`")
    lines.append(f"- Target: {target}")
    lines.append(f"- Solvent model: {rj.get('solvent_model', '—')}")
    lines.append(f"- df strategy: {rj.get('df_strategy', '—')}")
    sc = rj.get("sign_convention") or pw.get("sign_convention") or "—"
    lines.append(f"- Sign convention: {sc}")
    lines.append("")
    lines.append("## Pairwise")
    lines.append("")
    lines.append("| Quantity | Value |")
    lines.append("|---|---|")
    ddg = pw.get("ddG")
    se = pw.get("SE")
    z = pw.get("z_SE")
    ci = pw.get("CI95")
    df = pw.get("df")
    tier = pw.get("tier")

    def _f(x, d=3):
        return f"{x:+.{d}f}" if isinstance(x, (int, float)) else "—"

    lines.append(f"| ΔΔG (kcal/mol) | {_f(ddg)} |")
    lines.append(f"| SE | {se if se is not None else '—'} |")
    lines.append(f"| df | {df if df is not None else '—'} |")
    lines.append(f"| CI95 | {ci if ci is not None else '—'} |")
    lines.append(f"| z_SE | {_f(z, 2)} |")
    lines.append(f"| tier | **{tier or '—'}** |")
    note = pw.get("note") or rj.get("note")
    if note:
        lines.append(f"| note | {note} |")
    lines.append("")
    lines.append("## Branches")
    lines.append("")
    for name, b in (("WT", wt), ("Variant", var)):
        if not b:
            continue
        lines.append(f"### {name} — `{b.get('label', '?')}`")
        lines.append("")
        lines.append(f"- branch_dir: `{b.get('branch_dir', '?')}`")
        lines.append(f"- n_seed: {b.get('n_seed', '—')}  "
                     f"degenerate: {b.get('degenerate', '—')}")
        lines.append(f"- mean_of_seeds: {_f(b.get('mean_of_seeds'))}  "
                     f"σ_btwn: {b.get('sigma_btwn', '—')}  "
                     f"σ_w_median: {b.get('sigma_w_median', '—')}")
        lines.append(f"- SE: {b.get('SE', '—')}  "
                     f"solvent: {b.get('solvent_model', '—')}  "
                     f"charge_audit: {b.get('charge_audit', '—')}")
        lines.append("")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# argparse subparser builders
# ---------------------------------------------------------------------------

def _add_common_options(p: argparse.ArgumentParser, *, target_required: bool
                        ) -> None:
    """Common options for every subcommand. ``--target`` is required for
    check / compare / refine / discover; optional for report.
    """
    p.add_argument("--target", required=target_required,
                   help="Target PDB file (required for check/compare/refine/discover; "
                        "optional for report).")
    p.add_argument("--output", default=None,
                   help="Output directory (default: outputs/cli_run_<timestamp>/).")
    p.add_argument("--solvent", choices=("auto", "pbsa", "gbsa"), default="auto",
                   help="Solvent model preference (default: auto).")
    p.add_argument("--df-strategy", choices=("conservative", "min"),
                   default="conservative",
                   help="Combined-df rule for ΔΔG (default: conservative).")
    p.add_argument("--strict-charge-audit", action="store_true",
                   help="Enable PR-NEW-B Q4 strict charge audit mode.")
    p.add_argument("--verbose", action="store_true",
                   help="DEBUG-level logging.")
    p.add_argument("--quiet", action="store_true",
                   help="Suppress stdout summary lines.")
    p.add_argument("--dry-run", action="store_true",
                   help="Validate inputs only; emit a placeholder result.")
    p.add_argument("--format", choices=("json", "markdown", "both"),
                   default="both",
                   help="Output format (default: both).")


def _add_check(sub: Any) -> None:
    p = sub.add_parser("check",
                       help="Single sequence binder verification.",
                       description="Capability Level 1: single sequence verdict.")
    _add_common_options(p, target_required=True)
    p.add_argument("--sequence", required=True,
                   help="Sequence (canonical 1-letter or PR-NEW-A formats).")


def _add_compare(sub: Any) -> None:
    p = sub.add_parser("compare",
                       help="Multi-sequence pairwise ranking.",
                       description="Capability Level 1: rank N sequences by ΔG.")
    _add_common_options(p, target_required=True)
    p.add_argument("--sequences", required=True,
                   help="Comma-separated sequences.")


def _add_refine(sub: Any) -> None:
    p = sub.add_parser("refine",
                       help="Iterative refinement (Cycle 1+, WT vs variants).",
                       description=(
                           "Capability Level 1: WT branch + per-variant branches "
                           "with ΔΔG via PR-NEW-B compute_branched_ddg semantics."))
    _add_common_options(p, target_required=True)
    p.add_argument("--wt", required=True, help="WT (reference) sequence.")
    p.add_argument("--variants", required=True,
                   help="Comma-separated variant sequences.")


def _add_discover(sub: Any) -> None:
    p = sub.add_parser("discover",
                       help="De novo discovery (Cycle 0, file-based — Day 1).",
                       description=(
                           "Day-1 scope: file-based candidate list. Real Stage 1-2 "
                           "(RFdiffusion + ProteinMPNN + AF2) hookup is Day 2-3."))
    _add_common_options(p, target_required=True)
    p.add_argument("--sequence-file", required=True,
                   help="Plain-text file with one sequence per line ('#' = comment).")
    p.add_argument("--ncaa-budget", type=int, default=None,
                   help="Advisory ncAA count budget (Day-1: not enforced).")


def _add_report(sub: Any) -> None:
    p = sub.add_parser("report",
                       help="Generate human-readable report from branched_ddg.json.",
                       description=(
                           "Reads a branched_ddg/0.1 result and renders a markdown / "
                           "JSON summary. HTML is reserved for PR-NEW-E."))
    _add_common_options(p, target_required=False)
    p.add_argument("--result", required=True,
                   help="Path to a branched_ddg JSON result.")
    p.add_argument("--report-format", dest="format",
                   choices=("json", "markdown", "html", "both"),
                   default="markdown",
                   help="Render format (default: markdown). 'html' is "
                        "reserved for PR-NEW-E and currently falls back to markdown.")


# ---------------------------------------------------------------------------
# Misc helpers
# ---------------------------------------------------------------------------

def _split_csv(s: str) -> List[str]:
    """Split a comma-separated CLI string, trimming whitespace and empties."""
    if not s:
        return []
    return [tok.strip() for tok in s.split(",") if tok.strip()]


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="updd",
        description=(
            f"UPDD CLI {SCHEMA_VERSION} — Capability Level 1 user-facing layer. "
            "Integrates PR-NEW-A (sequence parser) + PR-NEW-B (branched ΔΔG)."
        ),
    )
    parser.add_argument(
        "--version", action="version",
        version=f"{SCHEMA_VERSION} ({TOOL_NAME})",
    )
    sub = parser.add_subparsers(
        dest="cmd", required=True,
        title="subcommands", metavar="{check,compare,refine,discover,report}",
    )
    _add_check(sub)
    _add_compare(sub)
    _add_refine(sub)
    _add_discover(sub)
    _add_report(sub)

    args = parser.parse_args(argv)
    ctx = CLIContext.from_args(args)
    handlers = {
        "check": cmd_check,
        "compare": cmd_compare,
        "refine": cmd_refine,
        "discover": cmd_discover,
        "report": cmd_report,
    }
    try:
        return handlers[args.cmd](args, ctx)
    except NcAAValidationError as exc:
        print(f"[updd_cli] NcAAValidationError: {exc}", file=sys.stderr)
        if getattr(exc, "suggestions", None):
            print(f"[updd_cli]   suggestions: {', '.join(exc.suggestions)}",
                  file=sys.stderr)
        return 4
    except SequenceFormatError as exc:
        print(f"[updd_cli] SequenceFormatError: {exc}", file=sys.stderr)
        return 3
    except FileNotFoundError as exc:
        print(f"[updd_cli] FileNotFoundError: {exc}", file=sys.stderr)
        return 5
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001
        print(f"[updd_cli] {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
