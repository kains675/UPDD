#!/usr/bin/env python
"""utils/branched_ddg.py — UPDD Branched ΔΔG Engine.

Schema: ``branched_ddg/0.1`` (PR-NEW-B Option A, 2026-04-28)

Aggregates per-seed ``mmpbsa_summary.json`` files (or ``mmgbsa_summary.json``)
into a same-protocol pairwise WT vs Variant ΔΔG result with full statistical
decomposition (σ_btwn, σ_w_median, SE, df, t_crit, CI95, z_SE, tier).

Capabilities:
    1. Per-seed JSON aggregation (auto-detect PBSA / GBSA solvent)
    2. WT branch + Variant branch processing
    3. Tier classification (X.A / X.B / X.C / X.D)
    4. WT control bit-identical verification (PR-12)
    5. Hybrid charge audit (default advisory; ``--strict-charge-audit`` = gate)

CLI::

    python -m utils.branched_ddg \\
        --wt outputs/1EBP_WT_calib_*/mmpbsa_results_fix4_postpatch \\
        --variant outputs/1EBP_MTR13_calib_*/mmpbsa_results_fix4_postpatch \\
        --output /tmp/demo_1ebp/

Tier rules (formalized §7 of capability analysis):

    X.A — CI95 excludes 0 (sign + magnitude significant)
    X.B — |z_SE| ≥ 2.0 AND CI95 includes 0 (sign-significant only)
    X.C — |z_SE| < 2.0 AND CI95 includes 0 (insufficient evidence)
    X.D — n_seed < 3 (insufficient sampling, conservative SE proxy used)

Statistical conventions (matching Phase α §3.4):

    per-seed mean     = mean(per-snap ΔG)              # arithmetic
    per-seed σ_w      = stdev(per-snap ΔG, ddof=1)     # sample std
    mean_of_seeds     = mean(per-seed means)           # arithmetic
    σ_btwn            = stdev(per-seed means, ddof=1)  # sample std
    σ_w_median        = median(per-seed σ_w)
    SE_branch         = σ_btwn / sqrt(n_seed)          # n_seed >= 2
    SE_branch (degen.) = σ_w / sqrt(n_snaps)           # n_seed == 1
    ΔΔG               = mean_of_seeds(variant) − mean_of_seeds(wt)
    SE_combined       = sqrt(SE_variant² + SE_wt²)
    df                = (n_variant − 1) + (n_wt − 1)   # combined; min(n−1) when degenerate
    CI95              = ΔΔG ± t_crit(df) × SE_combined
    z_SE              = ΔΔG / SE_combined              # signed (sign convention)
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import math
import os
import re
import statistics
import sys
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Constants / Schema
# ---------------------------------------------------------------------------

SCHEMA_VERSION = "branched_ddg/0.1"
TOOL_NAME = "utils.branched_ddg"

# Two-tail Student's t critical at α=0.05 (95% CI). Hard-coded LUT for
# small df. Beyond df=10 we fall back to z=1.96 (large-sample limit).
_T_CRIT_95: Dict[int, float] = {
    1: 12.706,
    2: 4.303,
    3: 3.182,
    4: 2.776,
    5: 2.571,
    6: 2.447,
    7: 2.365,
    8: 2.306,
    9: 2.262,
    10: 2.228,
}
_Z_LARGE = 1.960


# ---------------------------------------------------------------------------
# Public exception classes
# ---------------------------------------------------------------------------

class BranchedDDGError(Exception):
    """Base class for branched-ΔΔG errors."""


class ChargeAuditError(BranchedDDGError):
    """Raised in ``--strict-charge-audit`` mode when any branch fails Σq verify."""


class InsufficientSamplingError(BranchedDDGError):
    """Raised when a branch has zero usable seeds.

    Note: ``n_seed < 3`` does NOT raise — it is reported as tier ``X.D``.
    Only fully-empty branches escalate to this error.
    """


# ---------------------------------------------------------------------------
# Result containers
# ---------------------------------------------------------------------------

@dataclass
class SeedRecord:
    seed_label: str
    mean_dg: float
    sigma_w: float
    n_snap: int
    dgs: List[float]
    source_json: str
    schema_version: str
    solvent_model: str

    def to_dict(self) -> Dict[str, Any]:
        return {
            "seed_label": self.seed_label,
            "mean_dG": self.mean_dg,
            "sigma_w": self.sigma_w,
            "n_snap": self.n_snap,
            "dgs": self.dgs,
            "source_json": self.source_json,
            "schema_version": self.schema_version,
            "solvent_model": self.solvent_model,
        }


@dataclass
class BranchAggregate:
    label: str
    branch_dir: str
    per_seed: Dict[str, SeedRecord]
    mean_of_seeds: float
    sigma_btwn: Optional[float]
    sigma_w_median: float
    sigma_w_max: float
    SE: float
    n_seed: int
    n_snap_per_seed: int
    solvent_model: str
    degenerate: bool
    charge_audit: str        # "PASS" | "WARN" | "FAIL"
    charge_audit_detail: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "label": self.label,
            "branch_dir": self.branch_dir,
            "per_seed": {k: v.to_dict() for k, v in self.per_seed.items()},
            "mean_of_seeds": self.mean_of_seeds,
            "sigma_btwn": self.sigma_btwn,
            "sigma_w_median": self.sigma_w_median,
            "sigma_w_max": self.sigma_w_max,
            "SE": self.SE,
            "n_seed": self.n_seed,
            "n_snap_per_seed": self.n_snap_per_seed,
            "solvent_model": self.solvent_model,
            "degenerate": self.degenerate,
            "charge_audit": self.charge_audit,
            "charge_audit_detail": self.charge_audit_detail,
        }


@dataclass
class BranchedDDGResult:
    """Schema-typed branched ΔΔG result.

    This is the public dataclass returned by :func:`compute_branched_ddg`.
    Use :meth:`to_dict` to obtain a JSON-serializable representation that
    matches schema ``branched_ddg/0.1``.
    """

    schema: str
    tool: str
    generated_at: str
    target_id: str
    solvent_model: str
    df_strategy: str
    wt: BranchAggregate
    variant: BranchAggregate
    ddg: float
    SE_combined: float
    df: int
    t_crit: float
    CI95: Tuple[float, float]
    z_SE: float
    tier: str
    sign_convention: str
    note: Optional[str]
    charge_audit: str          # rolled-up across both branches
    charge_audit_strict: bool
    wt_bit_identical: Optional[bool]
    wt_bit_identical_detail: Optional[Dict[str, Any]]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema": self.schema,
            "tool": self.tool,
            "generated_at": self.generated_at,
            "target_id": self.target_id,
            "solvent_model": self.solvent_model,
            "df_strategy": self.df_strategy,
            "branches": {
                "wt": self.wt.to_dict(),
                "variant": self.variant.to_dict(),
            },
            "pairwise": {
                "ddG": self.ddg,
                "SE": self.SE_combined,
                "df": self.df,
                "t_crit": self.t_crit,
                "CI95": list(self.CI95),
                "z_SE": self.z_SE,
                "tier": self.tier,
                "sign_convention": self.sign_convention,
                "note": self.note,
            },
            "charge_audit": self.charge_audit,
            "charge_audit_strict": self.charge_audit_strict,
            "wt_bit_identical": self.wt_bit_identical,
            "wt_bit_identical_detail": self.wt_bit_identical_detail,
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Match ``..._calib_s<NUM>`` or ``..._calib_s<NUM>_anything`` (e.g.
# ``_reseed55``). Production layout. Falls back to bare ``_s<NUM>``
# (synthetic / test layouts).
_SEED_RE = re.compile(r"_calib_(s\d+(?:_[A-Za-z0-9]+)?)/?")
_SEED_BARE_RE = re.compile(r"_(s\d+(?:_[A-Za-z0-9]+)?)/?")


def _t_crit(df: int) -> float:
    """Two-tail Student's t critical at α=0.05 for the given df.

    Hard-coded LUT for ``df ∈ {1..10}``; falls back to ``z=1.96`` for ``df > 10``
    (large-sample normal approximation). ``df < 1`` raises.
    """
    if df < 1:
        raise ValueError(f"df must be >= 1, got {df}")
    if df in _T_CRIT_95:
        return _T_CRIT_95[df]
    return _Z_LARGE


def _md5sum_file(path: str) -> str:
    """Stream-hash a file. Returns hex digest. Raises if file missing."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _detect_solvent(summary: Dict[str, Any]) -> str:
    """Read solvent model from a summary JSON (PBSA / GBSA).

    Returns one of ``"pbsa"``, ``"gbsa"``, or ``"unknown"``.

    Phase α convention (schema 0.6.9):
        - run_mmpbsa.py writes ``solvent_model: "pbsa"``.
        - run_mmgbsa.py writes ``solvent_model: "gbsa"`` (or absent on older
          versions; fall back to ``gb_model`` field if so).
    """
    sm = summary.get("solvent_model")
    if isinstance(sm, str) and sm:
        s = sm.lower().strip()
        if s in {"pbsa", "gbsa"}:
            return s
        # Older mmgbsa schemas use ``gb_model`` ("n2", "obc", "obc2").
        # Treat as GBSA umbrella.
        if s in {"n2", "obc", "obc2", "gbn2", "obc2_gbn", "gb"}:
            return "gbsa"
        return s
    gb = summary.get("gb_model")
    if isinstance(gb, str) and gb:
        return "gbsa"
    return "unknown"


def _seed_label_from_dir(branch_dir: str, fallback_idx: int = 0) -> str:
    """Extract canonical seed label (e.g. ``s7``, ``s19_reseed55``) from a dir.

    Walks parent dirs (the seed token is normally on the parent of the
    ``mmpbsa_results_*`` dir). Falls back to ``unk{idx}`` if not parseable.
    """
    p = Path(branch_dir).resolve()
    cands = (p, p.parent, p.parent.parent)
    # Pass 1: production pattern (`_calib_s*`).
    for cand in cands:
        m = _SEED_RE.search(str(cand) + "/")
        if m:
            return m.group(1)
    # Pass 2: bare `_s<N>` token (synthetic test layouts).
    for cand in cands:
        m = _SEED_BARE_RE.search(str(cand) + "/")
        if m:
            return m.group(1)
    return f"unk{fallback_idx}"


def _load_seed_summary(summary_json: str, *, seed_label: Optional[str] = None) -> SeedRecord:
    """Load one ``mmpbsa_summary.json`` / ``mmgbsa_summary.json`` into SeedRecord.

    Reads ``schema_version``, ``solvent_model`` / ``gb_model`` for solvent
    auto-detect, and ``results[i].delta_g_kcal`` for per-snap ΔG list.
    Computes mean and σ_w (ddof=1) from the per-snap array, NOT from the
    summary's ``mean_dg`` field — this is intentional, so that seeds that
    were re-aggregated externally are still consistent.
    """
    with open(summary_json, "r", encoding="utf-8") as f:
        summary = json.load(f)
    schema = str(summary.get("schema_version", "unknown"))
    solvent = _detect_solvent(summary)
    rows = summary.get("results") or []
    dgs: List[float] = []
    for r in rows:
        if not isinstance(r, dict):
            continue
        try:
            dgs.append(float(r["delta_g_kcal"]))
        except (KeyError, TypeError, ValueError):
            continue
    if not dgs:
        raise BranchedDDGError(
            f"no per-snap ΔG values parseable from {summary_json}"
        )
    n = len(dgs)
    mean = statistics.fmean(dgs)
    if n >= 2:
        sigma_w = statistics.stdev(dgs)
    else:
        sigma_w = 0.0
    if seed_label is None:
        seed_label = _seed_label_from_dir(os.path.dirname(summary_json))
    return SeedRecord(
        seed_label=seed_label,
        mean_dg=mean,
        sigma_w=sigma_w,
        n_snap=n,
        dgs=dgs,
        source_json=summary_json,
        schema_version=schema,
        solvent_model=solvent,
    )


def _find_seed_summaries(branch_dir: str) -> List[Tuple[str, str]]:
    """Locate ``*summary.json`` files under ``branch_dir``.

    The branch_dir spec accepts three forms:
        1. A glob pattern that matches one or more ``mmpbsa_results_*`` dirs
           directly (e.g. ``outputs/1EBP_WT_calib_s*/mmpbsa_results_fix4_postpatch``).
        2. A single concrete ``mmpbsa_results_*`` dir (returns just that).
        3. A parent dir under which one or more ``*_calib_s*/mmpbsa_results_*``
           subdirs live (form 3 walks up to two levels deep).

    Returns a sorted list of ``(seed_label, summary_json_path)``.
    """
    candidates: List[str] = []
    expanded = sorted(glob.glob(branch_dir))
    if expanded:
        candidates = expanded
    else:
        candidates = [branch_dir]
    summaries: List[Tuple[str, str]] = []
    seen_paths: set = set()

    def _claim(seed_dir: str, summary_path: str) -> None:
        if summary_path in seen_paths:
            return
        seen_paths.add(summary_path)
        seed = _seed_label_from_dir(seed_dir, fallback_idx=len(summaries))
        summaries.append((seed, summary_path))

    for c in candidates:
        if not os.path.isdir(c):
            continue
        # Form 1/2: mmpbsa_summary.json directly under c.
        direct_hit = False
        for nm in ("mmpbsa_summary.json", "mmgbsa_summary.json"):
            p = os.path.join(c, nm)
            if os.path.isfile(p):
                _claim(c, p)
                direct_hit = True
                break
        if direct_hit:
            continue
        # Form 3a: c contains *_calib_s* dirs whose mmpbsa_results_*/ holds the
        # summary; or c contains mmpbsa_results_* dirs directly.
        for sub in sorted(os.listdir(c)):
            subp = os.path.join(c, sub)
            if not os.path.isdir(subp):
                continue
            # Direct child has the summary?
            sub_hit = False
            for nm in ("mmpbsa_summary.json", "mmgbsa_summary.json"):
                p = os.path.join(subp, nm)
                if os.path.isfile(p):
                    _claim(subp, p)
                    sub_hit = True
                    break
            if sub_hit:
                continue
            # Form 3b: subp is a seed dir; look inside its mmpbsa_results_*/.
            for sub2 in sorted(os.listdir(subp)):
                sub2p = os.path.join(subp, sub2)
                if not os.path.isdir(sub2p):
                    continue
                for nm in ("mmpbsa_summary.json", "mmgbsa_summary.json"):
                    p = os.path.join(sub2p, nm)
                    if os.path.isfile(p):
                        _claim(sub2p, p)
                        break
    # Stable order: (seed_label, path)
    summaries.sort()
    return summaries


def _charge_audit_branch(branch_dir: str, *, strict: bool = False) -> Tuple[str, Dict[str, Any]]:
    """Run a Σq audit on every ``params/*_params_manifest.json`` under the
    seed dirs that compose the branch.

    Best-effort: if no ``params/`` subdir is present (e.g. a WT branch with
    no ncAA), returns ``("PASS", {note: "no params manifests — nothing to audit"})``.

    Tolerance: 1e-6 absolute Σq. Method: reuses
    ``utils_draft.sigma_q_verify.verify_one`` for the 4-method consensus.

    Mode:
        ``strict=False`` (default) — non-zero Σq (any seed/manifest) → return
        ``("WARN", details)``. Pipeline continues.
        ``strict=True`` — non-zero Σq raises :class:`ChargeAuditError`.
    """
    detail: Dict[str, Any] = {"manifests": [], "max_abs_sigma_q": 0.0}
    try:
        # Lazy import — sigma_q_verify is in utils_draft which may not be on path.
        repo_root = Path(__file__).resolve().parents[1]
        sys.path.insert(0, str(repo_root))
        from utils_draft.sigma_q_verify import verify_one  # type: ignore
    except Exception as exc:  # noqa: BLE001
        detail["error"] = f"sigma_q_verify import failed: {type(exc).__name__}: {exc}"
        return "WARN", detail

    seed_dirs = _branch_seed_dirs(branch_dir)
    if not seed_dirs:
        detail["note"] = "branch_dir resolved to no seeds — nothing to audit"
        return "PASS", detail

    any_warn = False
    any_fail = False
    for sd in seed_dirs:
        params_dir = os.path.join(sd, "params")
        if not os.path.isdir(params_dir):
            continue
        for fn in sorted(os.listdir(params_dir)):
            if not fn.endswith("_params_manifest.json"):
                continue
            mpath = os.path.join(params_dir, fn)
            try:
                with open(mpath, "r", encoding="utf-8") as f:
                    manifest = json.load(f)
            except Exception as exc:  # noqa: BLE001
                detail["manifests"].append(
                    {"manifest": mpath, "status": "PARSE_ERROR",
                     "error": f"{type(exc).__name__}: {exc}"}
                )
                any_warn = True
                continue
            xml_path = manifest.get("xml_path")
            if not xml_path or not os.path.isfile(xml_path):
                # Best-effort: try sibling ``<code>_gaff2.xml`` next to manifest.
                code = manifest.get("ncaa_code") or fn.split("_params_manifest.json")[0]
                cand = os.path.join(params_dir, f"{code}_gaff2.xml")
                if os.path.isfile(cand):
                    xml_path = cand
                else:
                    detail["manifests"].append(
                        {"manifest": mpath, "status": "XML_NOT_FOUND",
                         "xml_path": manifest.get("xml_path")}
                    )
                    any_warn = True
                    continue
            try:
                vr = verify_one(Path(xml_path), None, tol=1.0e-6)
                sigma = vr.consensus_sigma_q
                if sigma is None:
                    detail["manifests"].append(
                        {"manifest": mpath, "xml": xml_path,
                         "status": "NO_METHOD_RAN"}
                    )
                    any_warn = True
                    continue
                abs_dev = abs(sigma)
                detail["max_abs_sigma_q"] = max(detail["max_abs_sigma_q"], abs_dev)
                status = "PASS" if abs_dev <= 1.0e-6 else "FAIL"
                detail["manifests"].append({
                    "manifest": mpath,
                    "xml": xml_path,
                    "sigma_q": sigma,
                    "abs_dev_from_neutral": abs_dev,
                    "status": status,
                })
                if status == "FAIL":
                    any_fail = True
            except Exception as exc:  # noqa: BLE001
                detail["manifests"].append(
                    {"manifest": mpath, "xml": xml_path,
                     "status": "VERIFY_EXC",
                     "error": f"{type(exc).__name__}: {exc}"}
                )
                any_warn = True

    if any_fail:
        if strict:
            raise ChargeAuditError(
                f"strict charge audit FAILED for branch {branch_dir!r}: "
                f"max |Σq| = {detail['max_abs_sigma_q']:.3e}; details: {detail}"
            )
        return "WARN", detail
    if any_warn:
        return "WARN", detail
    return "PASS", detail


def _branch_seed_dirs(branch_dir: str) -> List[str]:
    """Resolve ``branch_dir`` to a list of seed dirs (parent of mmpbsa_results_*).

    Same dispatch logic as :func:`_find_seed_summaries`, but returns the
    parent directory of each summary JSON (the ``*_calib_s*`` dir).
    """
    out: List[str] = []
    for _seed, sj in _find_seed_summaries(branch_dir):
        # sj = .../X_calib_sN/mmpbsa_results_*/mmpbsa_summary.json
        parent = os.path.dirname(os.path.dirname(sj))
        out.append(parent)
    # Dedup, preserve order.
    seen: set = set()
    dedup: List[str] = []
    for p in out:
        if p not in seen:
            seen.add(p)
            dedup.append(p)
    return dedup


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------

def aggregate_branch(branch_dir: str, *,
                     label: Optional[str] = None,
                     expected_solvent: Optional[str] = None,
                     charge_audit: bool = False,
                     strict_charge_audit: bool = False) -> BranchAggregate:
    """Aggregate all per-seed summaries under ``branch_dir`` into a BranchAggregate.

    Args:
        branch_dir: Glob pattern, single ``mmpbsa_results_*`` dir, or parent dir.
        label: Free-form label written to the result. Defaults to the basename.
        expected_solvent: If set, raise on mismatch with detected solvent.
        charge_audit: Whether to run the Σq audit.
        strict_charge_audit: If True, FAIL → :class:`ChargeAuditError`.

    Returns:
        BranchAggregate with per-seed records and aggregate stats.

    Raises:
        InsufficientSamplingError: if zero seeds resolvable.
        ValueError: on solvent mismatch with ``expected_solvent``.
        ChargeAuditError: in strict mode if any manifest fails.
    """
    summaries = _find_seed_summaries(branch_dir)
    if not summaries:
        raise InsufficientSamplingError(
            f"branch_dir {branch_dir!r} resolved to zero summary JSONs"
        )

    per_seed: Dict[str, SeedRecord] = {}
    solvents: List[str] = []
    n_snap_seen: List[int] = []
    for seed_label, sj in summaries:
        rec = _load_seed_summary(sj, seed_label=seed_label)
        # Avoid clobbering on duplicate seed_label from messy globs.
        if seed_label in per_seed:
            seed_label = f"{seed_label}_dup{len(per_seed)}"
            rec.seed_label = seed_label
        per_seed[seed_label] = rec
        solvents.append(rec.solvent_model)
        n_snap_seen.append(rec.n_snap)

    # Solvent consistency within branch.
    distinct = sorted(set(s for s in solvents if s and s != "unknown"))
    if len(distinct) > 1:
        raise ValueError(
            f"solvent mismatch within branch {branch_dir!r}: "
            f"per-seed solvent_model = {sorted(set(solvents))}"
        )
    branch_solvent = distinct[0] if distinct else (solvents[0] if solvents else "unknown")
    if expected_solvent is not None and branch_solvent != expected_solvent:
        raise ValueError(
            f"solvent mismatch — branch {branch_dir!r} is {branch_solvent}, "
            f"expected {expected_solvent}"
        )

    seed_means = [r.mean_dg for r in per_seed.values()]
    seed_sigmas = [r.sigma_w for r in per_seed.values()]
    n_seed = len(seed_means)

    mean_of_seeds = statistics.fmean(seed_means)
    if n_seed >= 2:
        sigma_btwn: Optional[float] = statistics.stdev(seed_means)
        SE = sigma_btwn / math.sqrt(n_seed)
        degenerate = False
    else:
        sigma_btwn = None
        # Degenerate single-seed branch: SE proxy = σ_w / sqrt(n_snap).
        sole = next(iter(per_seed.values()))
        SE = sole.sigma_w / math.sqrt(max(sole.n_snap, 1))
        degenerate = True

    sigma_w_median = statistics.median(seed_sigmas) if seed_sigmas else 0.0
    sigma_w_max = max(seed_sigmas) if seed_sigmas else 0.0
    n_snap_per_seed = n_snap_seen[0] if n_snap_seen and len(set(n_snap_seen)) == 1 else (
        n_snap_seen[0] if n_snap_seen else 0
    )

    # Charge audit (optional).
    audit_status = "PASS"
    audit_detail: Dict[str, Any] = {"requested": charge_audit, "strict": strict_charge_audit}
    if charge_audit:
        audit_status, more = _charge_audit_branch(branch_dir, strict=strict_charge_audit)
        audit_detail.update(more)

    if label is None:
        label = os.path.basename(os.path.normpath(branch_dir).rstrip("/")) or "branch"

    return BranchAggregate(
        label=label,
        branch_dir=branch_dir,
        per_seed=per_seed,
        mean_of_seeds=mean_of_seeds,
        sigma_btwn=sigma_btwn,
        sigma_w_median=sigma_w_median,
        sigma_w_max=sigma_w_max,
        SE=SE,
        n_seed=n_seed,
        n_snap_per_seed=n_snap_per_seed,
        solvent_model=branch_solvent,
        degenerate=degenerate,
        charge_audit=audit_status,
        charge_audit_detail=audit_detail,
    )


def tier_classify(ddg: float, ci_lower: float, ci_upper: float,
                  z_se: float, n_seed: int) -> str:
    """Classify ΔΔG into X.A / X.B / X.C / X.D.

    Rules (formalized §7 of capability analysis, 2026-04-28):

        X.D — n_seed < 3 (insufficient sampling, conservative SE proxy used)
        X.A — CI95 excludes 0 (sign + magnitude significant)
        X.B — |z_SE| ≥ 2.0 AND CI95 includes 0 (sign-significant only)
        X.C — |z_SE| < 2.0 AND CI95 includes 0 (insufficient evidence)

    The X.D rule fires *first* — even a "high z_SE" with n_seed=1 is X.D
    because the SE estimator is a proxy (σ_w/√n_snap) rather than σ_btwn,
    and the user explicitly demanded conservative reporting in this regime.
    """
    if n_seed < 3:
        return "X.D"
    excludes_zero = (ci_lower > 0.0) or (ci_upper < 0.0)
    if excludes_zero:
        return "X.A"
    if abs(z_se) >= 2.0:
        return "X.B"
    return "X.C"


def compute_branched_ddg(wt_dir: str, variant_dir: str, *,
                         wt_label: Optional[str] = None,
                         variant_label: Optional[str] = None,
                         target_id: str = "",
                         df_strategy: str = "conservative",
                         expected_solvent: Optional[str] = None,
                         charge_audit: bool = False,
                         strict_charge_audit: bool = False,
                         wt_pre_dir: Optional[str] = None) -> BranchedDDGResult:
    """Compute branched ΔΔG = ⟨⟨ΔG_variant⟩⟩ − ⟨⟨ΔG_wt⟩⟩ end-to-end.

    Args:
        wt_dir, variant_dir: Branch directories (see :func:`aggregate_branch`).
        wt_label, variant_label: Free-form labels written to result.
        target_id: Optional target identifier (e.g. ``"1EBP"``); propagated
            to the result.
        df_strategy: ``"conservative"`` (default) uses
            ``df = (n_v − 1) + (n_wt − 1)`` (combined). ``"min"`` uses
            ``df = min(n_v − 1, n_wt − 1)``.
        expected_solvent: If set (``"pbsa"`` / ``"gbsa"``), raise on
            cross-branch mismatch.
        charge_audit: Run Σq audit on each branch.
        strict_charge_audit: Failures raise :class:`ChargeAuditError`.
        wt_pre_dir: If provided, run the PR-12 WT bit-identical check by
            comparing every WT seed's summary JSON md5 against the same-named
            file under ``wt_pre_dir`` (mirroring the per-seed structure).

    Returns:
        BranchedDDGResult.

    Raises:
        ValueError: on solvent mismatch between WT and variant.
        InsufficientSamplingError: on zero-seed branch.
        ChargeAuditError: under strict charge audit.
    """
    wt = aggregate_branch(
        wt_dir, label=wt_label or "wt",
        expected_solvent=expected_solvent,
        charge_audit=charge_audit,
        strict_charge_audit=strict_charge_audit,
    )
    variant = aggregate_branch(
        variant_dir, label=variant_label or "variant",
        expected_solvent=expected_solvent,
        charge_audit=charge_audit,
        strict_charge_audit=strict_charge_audit,
    )

    # Same-protocol invariant: solvent must match across branches.
    if wt.solvent_model != "unknown" and variant.solvent_model != "unknown" \
            and wt.solvent_model != variant.solvent_model:
        raise ValueError(
            f"solvent mismatch — WT={wt.solvent_model}, "
            f"Variant={variant.solvent_model}; same-protocol invariant violated"
        )
    solvent = wt.solvent_model if wt.solvent_model != "unknown" else variant.solvent_model

    # ΔΔG and combined SE.
    ddg = variant.mean_of_seeds - wt.mean_of_seeds
    SE_combined = math.sqrt(variant.SE ** 2 + wt.SE ** 2)

    # df: conservative (sum) vs min.
    if df_strategy == "min":
        df_v = max(1, variant.n_seed - 1)
        df_w = max(1, wt.n_seed - 1)
        df = min(df_v, df_w)
    else:  # "conservative" — combined
        df = max(1, (variant.n_seed - 1) + (wt.n_seed - 1))

    t_crit = _t_crit(df)
    half = t_crit * SE_combined
    ci_lo = ddg - half
    ci_hi = ddg + half
    z_se = ddg / SE_combined if SE_combined > 0 else 0.0

    n_min = min(variant.n_seed, wt.n_seed)
    tier = tier_classify(ddg, ci_lo, ci_hi, z_se, n_min)
    note = None
    if tier == "X.D":
        if variant.degenerate or wt.degenerate:
            note = (
                f"insufficient sampling: variant n_seed={variant.n_seed}, "
                f"wt n_seed={wt.n_seed}; SE proxy used for degenerate branch"
            )
        else:
            note = (
                f"insufficient sampling: min(n_seed)={n_min} < 3"
            )

    # Roll up branch-level audit into a single status (FAIL > WARN > PASS).
    rolled = "PASS"
    statuses = [wt.charge_audit, variant.charge_audit]
    if "FAIL" in statuses:
        rolled = "FAIL"
    elif "WARN" in statuses:
        rolled = "WARN"

    # PR-12 WT bit-identical (optional).
    wt_bit_id: Optional[bool] = None
    wt_bit_id_detail: Optional[Dict[str, Any]] = None
    if wt_pre_dir is not None:
        bit_res = verify_wt_bit_identical(wt_pre_dir, wt_dir)
        wt_bit_id = bit_res["all_match"]
        wt_bit_id_detail = bit_res

    return BranchedDDGResult(
        schema=SCHEMA_VERSION,
        tool=TOOL_NAME,
        generated_at=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        target_id=target_id,
        solvent_model=solvent,
        df_strategy=df_strategy,
        wt=wt,
        variant=variant,
        ddg=ddg,
        SE_combined=SE_combined,
        df=df,
        t_crit=t_crit,
        CI95=(ci_lo, ci_hi),
        z_SE=z_se,
        tier=tier,
        sign_convention=f"ddG = mean(variant={variant.label}) − mean(wt={wt.label})",
        note=note,
        charge_audit=rolled,
        charge_audit_strict=strict_charge_audit,
        wt_bit_identical=wt_bit_id,
        wt_bit_identical_detail=wt_bit_id_detail,
    )


def verify_wt_bit_identical(wt_pre_dir: str, wt_post_dir: str) -> Dict[str, Any]:
    """PR-12 helper: md5-compare per-seed summary JSONs between two WT branches.

    For each seed found under ``wt_post_dir``, look up the same seed_label
    under ``wt_pre_dir`` and md5sum-compare the corresponding ``mmpbsa_summary.json``
    (or ``mmgbsa_summary.json``) file. Returns:

        {
            "all_match": bool,
            "n_compared": int,
            "matches": int,
            "diffs": [
                {"seed": "s7", "pre_md5": "...", "post_md5": "...",
                 "pre_path": "...", "post_path": "..."},
                ...
            ],
            "missing_in_pre": ["sN", ...],
            "missing_in_post": ["sM", ...],
        }
    """
    pre_summaries = dict(_find_seed_summaries(wt_pre_dir))
    post_summaries = dict(_find_seed_summaries(wt_post_dir))

    # Some "missing" seeds may exist on one side only.
    pre_set = set(pre_summaries)
    post_set = set(post_summaries)
    missing_in_pre = sorted(post_set - pre_set)
    missing_in_post = sorted(pre_set - post_set)
    common = sorted(pre_set & post_set)

    diffs: List[Dict[str, Any]] = []
    matches = 0
    for seed in common:
        pre_p = pre_summaries[seed]
        post_p = post_summaries[seed]
        try:
            pre_md5 = _md5sum_file(pre_p)
            post_md5 = _md5sum_file(post_p)
        except OSError as exc:
            diffs.append({
                "seed": seed, "pre_path": pre_p, "post_path": post_p,
                "error": f"{type(exc).__name__}: {exc}",
            })
            continue
        if pre_md5 == post_md5:
            matches += 1
        else:
            diffs.append({
                "seed": seed, "pre_md5": pre_md5, "post_md5": post_md5,
                "pre_path": pre_p, "post_path": post_p,
            })
    n_compared = len(common)
    all_match = (matches == n_compared) and not missing_in_pre and not missing_in_post
    return {
        "all_match": all_match,
        "n_compared": n_compared,
        "matches": matches,
        "diffs": diffs,
        "missing_in_pre": missing_in_pre,
        "missing_in_post": missing_in_post,
    }


def write_report(result: BranchedDDGResult, output_dir: str) -> Tuple[str, str]:
    """Write JSON + Markdown report. Returns (json_path, md_path).

    Files:
        ``branched_ddg.json`` — schema ``branched_ddg/0.1`` payload.
        ``branched_ddg.md``   — human-readable summary table.
    """
    os.makedirs(output_dir, exist_ok=True)
    json_path = os.path.join(output_dir, "branched_ddg.json")
    md_path = os.path.join(output_dir, "branched_ddg.md")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(result.to_dict(), f, indent=2, default=_json_default)
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(_render_md(result))
    return json_path, md_path


def _json_default(obj: Any) -> Any:
    if isinstance(obj, tuple):
        return list(obj)
    if hasattr(obj, "to_dict"):
        return obj.to_dict()
    if hasattr(obj, "__dict__"):
        return obj.__dict__
    return str(obj)


def _fmt(x: Optional[float], digits: int = 3) -> str:
    if x is None:
        return "—"
    if isinstance(x, float) and not math.isfinite(x):
        return "nan/inf"
    return f"{x:+.{digits}f}" if isinstance(x, (int, float)) else str(x)


def _render_md(r: BranchedDDGResult) -> str:
    """Render a Phase α track_b_*.md-flavored Markdown report."""
    lines: List[str] = []
    lines.append(f"# Branched ΔΔG — {r.target_id or '(target unset)'}")
    lines.append("")
    lines.append(f"- **Schema**: `{r.schema}` (`{r.tool}`)")
    lines.append(f"- **Generated**: {r.generated_at}")
    lines.append(f"- **Solvent model**: {r.solvent_model}")
    lines.append(f"- **df strategy**: {r.df_strategy}  (df={r.df}, t_crit={r.t_crit:.3f})")
    lines.append(f"- **Sign convention**: {r.sign_convention}")
    lines.append(f"- **Charge audit**: {r.charge_audit}"
                 + (" (strict)" if r.charge_audit_strict else ""))
    if r.wt_bit_identical is not None:
        lines.append(f"- **WT bit-identical (PR-12)**: {r.wt_bit_identical}")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append("")
    lines.append("| Quantity | Value |")
    lines.append("|---|---|")
    lines.append(f"| ΔΔG (kcal/mol) | {_fmt(r.ddg)} |")
    lines.append(f"| SE_combined | {r.SE_combined:.3f} |")
    lines.append(f"| CI95 | [{_fmt(r.CI95[0])}, {_fmt(r.CI95[1])}] |")
    lines.append(f"| z_SE | {_fmt(r.z_SE, 2)} |")
    lines.append(f"| Tier | **{r.tier}** |")
    if r.note:
        lines.append(f"| Note | {r.note} |")
    lines.append("")
    lines.append("## Per-Branch Aggregates")
    lines.append("")
    for tag, b in (("WT", r.wt), ("Variant", r.variant)):
        lines.append(f"### {tag} — `{b.label}` ({b.branch_dir})")
        lines.append("")
        lines.append(f"- n_seed = {b.n_seed} (degenerate: {b.degenerate})")
        sb_str = "—" if b.sigma_btwn is None else f"{b.sigma_btwn:.3f}"
        lines.append(f"- mean_of_seeds = {b.mean_of_seeds:+.3f}, "
                     f"σ_btwn = {sb_str}, "
                     f"σ_w_median = {b.sigma_w_median:.3f}, "
                     f"σ_w_max = {b.sigma_w_max:.3f}")
        lines.append(f"- SE = {b.SE:.3f}, "
                     f"solvent = {b.solvent_model}, "
                     f"charge_audit = {b.charge_audit}")
        lines.append("")
        lines.append("| Seed | mean_dG | σ_w | n_snap |")
        lines.append("|---|---|---|---|")
        for sl, rec in sorted(b.per_seed.items()):
            lines.append(f"| {sl} | {rec.mean_dg:+.3f} | "
                         f"{rec.sigma_w:.3f} | {rec.n_snap} |")
        lines.append("")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="python -m utils.branched_ddg",
        description="UPDD Branched ΔΔG Engine (PR-NEW-B Option A; schema "
                    f"{SCHEMA_VERSION}).",
    )
    ap.add_argument("--wt", required=True,
                    help="WT branch dir (glob, single mmpbsa_results_* dir, "
                         "or parent dir).")
    ap.add_argument("--variant", required=True,
                    help="Variant branch dir (same forms as --wt).")
    ap.add_argument("--output", required=True,
                    help="Output directory; will receive branched_ddg.{json,md}.")
    ap.add_argument("--wt-label", default=None,
                    help="Free-form WT label for the report (default: dir basename).")
    ap.add_argument("--variant-label", default=None,
                    help="Free-form variant label.")
    ap.add_argument("--target-id", default="",
                    help="Target identifier (e.g. 1EBP). Propagated to JSON.")
    ap.add_argument("--df-strategy", choices=("conservative", "min"),
                    default="conservative",
                    help="Combined-df rule. conservative = (n_v-1)+(n_wt-1), "
                         "min = min(n_v-1, n_wt-1).")
    ap.add_argument("--expected-solvent", choices=("pbsa", "gbsa"), default=None,
                    help="If set, raise on solvent mismatch.")
    ap.add_argument("--charge-audit", action="store_true",
                    help="Run hybrid Σq charge audit on each branch (advisory).")
    ap.add_argument("--strict-charge-audit", action="store_true",
                    help="With --charge-audit: any FAIL → ChargeAuditError "
                         "(non-zero exit).")
    ap.add_argument("--check-wt-bit-identical", default=None,
                    help="PR-12 helper: WT pre-patch dir to md5-compare against.")
    return ap


def main(argv: Optional[List[str]] = None) -> int:
    ap = _build_parser()
    args = ap.parse_args(argv)
    if args.strict_charge_audit and not args.charge_audit:
        # Promote: strict implies enable.
        args.charge_audit = True
    try:
        result = compute_branched_ddg(
            wt_dir=args.wt,
            variant_dir=args.variant,
            wt_label=args.wt_label,
            variant_label=args.variant_label,
            target_id=args.target_id,
            df_strategy=args.df_strategy,
            expected_solvent=args.expected_solvent,
            charge_audit=args.charge_audit,
            strict_charge_audit=args.strict_charge_audit,
            wt_pre_dir=args.check_wt_bit_identical,
        )
    except ChargeAuditError as exc:
        print(f"[branched_ddg] ChargeAuditError: {exc}", file=sys.stderr)
        return 2
    except (BranchedDDGError, ValueError) as exc:
        print(f"[branched_ddg] ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 3
    json_path, md_path = write_report(result, args.output)
    print(f"[branched_ddg] schema={result.schema}")
    print(f"[branched_ddg] ddG={result.ddg:+.4f}  "
          f"SE={result.SE_combined:.4f}  "
          f"z_SE={result.z_SE:+.3f}  "
          f"tier={result.tier}")
    print(f"[branched_ddg] CI95=[{result.CI95[0]:+.4f}, {result.CI95[1]:+.4f}] "
          f"(df={result.df}, t_crit={result.t_crit:.3f})")
    print(f"[branched_ddg] charge_audit={result.charge_audit}"
          + (f" (strict)" if result.charge_audit_strict else ""))
    if result.wt_bit_identical is not None:
        print(f"[branched_ddg] wt_bit_identical={result.wt_bit_identical}")
    print(f"[branched_ddg] wrote: {json_path}")
    print(f"[branched_ddg] wrote: {md_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
