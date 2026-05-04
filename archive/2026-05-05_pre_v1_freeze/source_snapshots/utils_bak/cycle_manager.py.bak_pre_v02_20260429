#!/usr/bin/env python
"""utils/cycle_manager.py — UPDD Iterative Cycle Manager.

Schema: ``cycle_manager/0.1`` (PR-NEW-C Day 1, 2026-04-28)

Capability Level 1 project tracking + cycle history layer for the
iterative peptide drug discovery workflow described in
``v07_user_verdict/v0.7 Section 4 finalized.md`` §4.1.3.

Manages multi-cycle peptide drug discovery projects:
    - Cycle 0 (de novo) → wet lab feedback → Cycle 1 WT
    - Cycle 1+ (refinement) → wet lab feedback → Cycle 2 WT
    - ... continuing until convergence or plateau

State persisted as JSON at
``<output_root>/<project_id>/cycle_history.json`` (default
``output_root='outputs/projects'``). Integrates with PR-NEW-A
(:mod:`utils.sequence_parser` — sequence input) and PR-NEW-B
(:mod:`utils.branched_ddg` — ΔΔG predictions).

Public API:
    - :class:`CycleHistory`       — full project state (dataclass).
    - :class:`CycleEntry`         — per-cycle record.
    - :class:`WetLabResult`       — assay readout (ITC/SPR/FP/BLI).
    - :class:`PredictedDDG`       — predicted ΔΔG record (PR-NEW-B linkage).
    - :class:`CycleManagerError`  + subclasses
      (:class:`ProjectNotFoundError`, :class:`InvalidCycleError`,
      :class:`KdParseError`).
    - :func:`init_project`        — create new project.
    - :func:`add_cycle`           — register Cycle N entry.
    - :func:`add_result`          — record wet lab K_d / IC50.
    - :func:`set_next_wt`         — update next-cycle WT (best validated).
    - :func:`compute_trajectory`  — ΔK_d trajectory + plateau detection.
    - :func:`detect_plateau`      — plateau classifier.
    - :func:`generate_report`     — markdown summary.

CLI::

    python -m utils.cycle_manager init       --project-id X --target T.pdb
    python -m utils.cycle_manager add-cycle  --project-id X --cycle 0 \\
        --variants "ICVVQDWGHHRCT,ICVVQDW(MTR)HHRCT"
    python -m utils.cycle_manager add-result --project-id X --cycle 0 \\
        --sequence "ICVVQDW(MTR)HHRCT" --kd 80nM --assay ITC
    python -m utils.cycle_manager set-wt     --project-id X --source-cycle 0 \\
        --best auto
    python -m utils.cycle_manager trajectory --project-id X
    python -m utils.cycle_manager report     --project-id X

K_d → ΔG conversion uses the canonical biochemistry relation
``ΔG_bind = R T ln(K_d / 1 M)`` (equivalently ``-R T ln(K_a)``) with
``R = 1.987204e-3 kcal mol⁻¹ K⁻¹`` and ``T = 298.15 K``. Sign convention:
ΔG_bind < 0 for tight binders (K_d ≪ 1 M); 80 nM at 298.15 K maps to
ΔG ≈ -9.68 kcal/mol.

Plateau detection (default): the last ``lookback=3`` cycle-over-cycle
``|ΔΔG_wet|`` values are all below ``threshold_kcal=0.5`` (≈ 2.3× K_d
ratio) → ``plateau_detected = True``.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import tempfile
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union


# ---------------------------------------------------------------------------
# Constants / Schema
# ---------------------------------------------------------------------------

SCHEMA_VERSION = "cycle_manager/0.1"
TOOL_NAME = "utils.cycle_manager"

# Gas constant in kcal mol⁻¹ K⁻¹ (CODATA 2018: R = 8.314462618 J mol⁻¹ K⁻¹).
R_KCAL_PER_MOL_K = 1.987204258640832e-3
DEFAULT_TEMPERATURE_K = 298.15

# K_d unit-suffix table → multiplier to convert to molar (M).
_KD_UNIT_TABLE: Dict[str, float] = {
    "":   1.0,        # no suffix → molar
    "m":  1.0,
    "mm": 1.0e-3,     # millimolar
    "um": 1.0e-6,     # micromolar
    "μm": 1.0e-6,     # micromolar (Greek mu)
    "nm": 1.0e-9,     # nanomolar
    "pm": 1.0e-12,    # picomolar
    "fm": 1.0e-15,    # femtomolar
}

# K_d parser regex: ``80``, ``80nM``, ``80 nM``, ``8e-8``, ``0.08uM``,
# ``80,000 pM`` (commas stripped before regex).
_KD_RE = re.compile(
    r"^\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s*([a-zA-Zμ]{0,2})?\s*$"
)


# ---------------------------------------------------------------------------
# Public exception classes
# ---------------------------------------------------------------------------

class CycleManagerError(Exception):
    """Base class for cycle-manager errors."""


class ProjectNotFoundError(CycleManagerError):
    """Raised when ``cycle_history.json`` for the requested project is missing."""


class InvalidCycleError(CycleManagerError):
    """Raised when a cycle operation is structurally invalid

    (e.g. duplicate cycle_n, negative cycle_n, missing prerequisite cycle).
    """


class KdParseError(CycleManagerError):
    """Raised when a K_d string cannot be parsed into a float in molar.

    Also used for downstream numerical pathologies such as zero-confidence
    multi-assay aggregation (NaN result).
    """


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class WetLabResult:
    """One wet-lab assay readout for a specific sequence in a cycle.

    Attributes:
        sequence: Variant sequence (string form, may contain ncAA tags).
        kd_molar: Measured K_d in molar (M). Always > 0.
        assay:    Assay type label (free-form: ``ITC``, ``SPR``, ``FP``,
                  ``BLI``, ``ELISA``, etc.).
        ddg_calc: ΔG = -RT ln(K_d) computed at ``temperature_k`` (kcal/mol).
        confidence: Optional weight (≥ 0). Used for inverse-variance-style
                  weighted aggregation when multiple assays are present.
                  ``None`` (default) is treated as 1.0 for simple averaging.
        temperature_k: Temperature at which K_d was measured.
        notes:    Free-form text (e.g. buffer conditions, replicate count).
    """

    sequence: str
    kd_molar: float
    assay: str = "unknown"
    ddg_calc: Optional[float] = None
    confidence: Optional[float] = None
    temperature_k: float = DEFAULT_TEMPERATURE_K
    notes: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "WetLabResult":
        return cls(
            sequence=str(d["sequence"]),
            kd_molar=float(d["kd_molar"]),
            assay=str(d.get("assay", "unknown")),
            ddg_calc=(float(d["ddg_calc"]) if d.get("ddg_calc") is not None else None),
            confidence=(float(d["confidence"]) if d.get("confidence") is not None else None),
            temperature_k=float(d.get("temperature_k", DEFAULT_TEMPERATURE_K)),
            notes=(str(d["notes"]) if d.get("notes") is not None else None),
        )


@dataclass
class PredictedDDG:
    """One predicted ΔΔG record (typically sourced from PR-NEW-B output).

    Attributes:
        sequence: Variant sequence.
        ddg:      Predicted ΔΔG (kcal/mol; sign convention =
                  ``mean(variant) − mean(WT)``).
        se:       Standard error of the prediction (kcal/mol).
        tier:     Branched-ΔΔG tier classification (``X.A``/``X.B``/
                  ``X.C``/``X.D``). Optional.
        z_se:     ``ddg / se`` (signed). Optional.
        ci95:     2-tuple lower/upper 95% confidence interval (kcal/mol).
                  Optional; stored as a list when serialised.
        source:   Provenance hint (e.g. branched_ddg JSON path).
    """

    sequence: str
    ddg: float
    se: Optional[float] = None
    tier: Optional[str] = None
    z_se: Optional[float] = None
    ci95: Optional[Tuple[float, float]] = None
    source: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "sequence": self.sequence,
            "ddg": self.ddg,
            "se": self.se,
            "tier": self.tier,
            "z_se": self.z_se,
            "ci95": (list(self.ci95) if self.ci95 is not None else None),
            "source": self.source,
        }
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "PredictedDDG":
        ci95 = d.get("ci95")
        if ci95 is not None and not isinstance(ci95, tuple):
            ci95 = tuple(ci95)  # type: ignore[assignment]
        return cls(
            sequence=str(d["sequence"]),
            ddg=float(d["ddg"]),
            se=(float(d["se"]) if d.get("se") is not None else None),
            tier=(str(d["tier"]) if d.get("tier") is not None else None),
            z_se=(float(d["z_se"]) if d.get("z_se") is not None else None),
            ci95=ci95,  # type: ignore[arg-type]
            source=(str(d["source"]) if d.get("source") is not None else None),
        )


@dataclass
class CycleEntry:
    """One cycle within an iterative peptide-design project.

    Attributes:
        cycle_n: 0-indexed cycle number. ``cycle_n == 0`` is de novo
                 (no WT — ``wt_sequence`` may be ``None``); ``cycle_n ≥ 1``
                 is refinement (``wt_sequence`` should be set).
        wt_sequence: WT reference for this cycle (None for Cycle 0 de novo).
        variants: List of variant sequences submitted for evaluation.
        predicted_ddg: List of :class:`PredictedDDG` (PR-NEW-B output).
        wet_lab_results: List of :class:`WetLabResult` ingested for this
            cycle (one or more per sequence; can be empty until wet lab
            data return).
        best_validated_sequence: Sequence picked from ``wet_lab_results``
            for promotion to the next cycle's WT (set by
            :func:`set_next_wt`). ``None`` until selection.
        synthesized: Whether the variants were synthesised (informational).
        notes: Free-form per-cycle annotation.
        created_at: ISO-8601 timestamp.
    """

    cycle_n: int
    wt_sequence: Optional[str] = None
    variants: List[str] = field(default_factory=list)
    predicted_ddg: List[PredictedDDG] = field(default_factory=list)
    wet_lab_results: List[WetLabResult] = field(default_factory=list)
    best_validated_sequence: Optional[str] = None
    synthesized: bool = False
    notes: Optional[str] = None
    created_at: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "cycle_n": self.cycle_n,
            "wt_sequence": self.wt_sequence,
            "variants": list(self.variants),
            "predicted_ddg": [p.to_dict() for p in self.predicted_ddg],
            "wet_lab_results": [w.to_dict() for w in self.wet_lab_results],
            "best_validated_sequence": self.best_validated_sequence,
            "synthesized": self.synthesized,
            "notes": self.notes,
            "created_at": self.created_at,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "CycleEntry":
        return cls(
            cycle_n=int(d["cycle_n"]),
            wt_sequence=(str(d["wt_sequence"])
                         if d.get("wt_sequence") is not None else None),
            variants=list(d.get("variants") or []),
            predicted_ddg=[PredictedDDG.from_dict(p)
                           for p in (d.get("predicted_ddg") or [])],
            wet_lab_results=[WetLabResult.from_dict(w)
                             for w in (d.get("wet_lab_results") or [])],
            best_validated_sequence=(
                str(d["best_validated_sequence"])
                if d.get("best_validated_sequence") is not None else None
            ),
            synthesized=bool(d.get("synthesized", False)),
            notes=(str(d["notes"]) if d.get("notes") is not None else None),
            created_at=str(d.get("created_at", "")),
        )


@dataclass
class CycleHistory:
    """Full persisted state for one iterative peptide-design project."""

    schema: str
    project_id: str
    target_pdb: str
    created_at: str
    metadata: Dict[str, Any] = field(default_factory=dict)
    cycles: List[CycleEntry] = field(default_factory=list)
    trajectory: Dict[str, Any] = field(default_factory=dict)
    output_root: str = "outputs/projects"

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema": self.schema,
            "project_id": self.project_id,
            "target_pdb": self.target_pdb,
            "created_at": self.created_at,
            "metadata": dict(self.metadata),
            "cycles": [c.to_dict() for c in self.cycles],
            "trajectory": dict(self.trajectory),
            "output_root": self.output_root,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "CycleHistory":
        return cls(
            schema=str(d.get("schema", SCHEMA_VERSION)),
            project_id=str(d["project_id"]),
            target_pdb=str(d["target_pdb"]),
            created_at=str(d.get("created_at", "")),
            metadata=dict(d.get("metadata") or {}),
            cycles=[CycleEntry.from_dict(c) for c in (d.get("cycles") or [])],
            trajectory=dict(d.get("trajectory") or {}),
            output_root=str(d.get("output_root", "outputs/projects")),
        )


# ---------------------------------------------------------------------------
# Internal helpers — paths, persistence, K_d parsing, ΔG conversion
# ---------------------------------------------------------------------------

def _project_dir(project_id: str, output_root: str) -> Path:
    """Return the canonical directory for a project's persisted state."""
    if not project_id or not project_id.strip():
        raise CycleManagerError("project_id must be non-empty")
    return Path(output_root) / project_id


def _history_path(project_id: str, output_root: str) -> Path:
    """Return the ``cycle_history.json`` path for ``project_id``."""
    return _project_dir(project_id, output_root) / "cycle_history.json"


def _now_iso() -> str:
    """ISO-8601 UTC timestamp (Z suffix)."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _save_history(history: CycleHistory) -> None:
    """Persist the history to ``cycle_history.json`` (atomic write).

    Uses NamedTemporaryFile in the same directory + ``os.replace`` for
    atomicity: a partial write cannot leave the file truncated.
    """
    target_dir = _project_dir(history.project_id, history.output_root)
    target_dir.mkdir(parents=True, exist_ok=True)
    target_path = target_dir / "cycle_history.json"
    payload = history.to_dict()
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=str(target_dir), delete=False, suffix=".tmp"
    ) as tf:
        json.dump(payload, tf, indent=2, default=str)
        tmp_name = tf.name
    os.replace(tmp_name, target_path)


def _load_history(project_id: str, *,
                  output_root: str = "outputs/projects") -> CycleHistory:
    """Load a project's ``cycle_history.json`` from disk.

    Raises:
        ProjectNotFoundError: if the JSON file is missing.
        CycleManagerError: on schema-decode errors.
    """
    path = _history_path(project_id, output_root)
    if not path.is_file():
        raise ProjectNotFoundError(
            f"project {project_id!r} not found at {path}"
        )
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
    except json.JSONDecodeError as exc:
        raise CycleManagerError(
            f"failed to decode {path}: {exc}"
        ) from exc
    history = CycleHistory.from_dict(data)
    # Re-anchor to the caller's output_root so subsequent saves write where
    # the caller expects (fixtures often pass tmp_path here even if the JSON
    # was written with a different on-disk root).
    history.output_root = output_root
    return history


def _parse_kd(raw: Union[str, float, int]) -> float:
    """Parse a K_d string (``80nM``, ``8e-8``, ``80 nM``, ``0.08uM``,
    ``80,000 pM``, ``5 fM``, etc.) into a float in molar (M).

    Numeric inputs (float / int) are accepted as already-molar values.

    Raises:
        KdParseError: on unparseable input, non-positive values, or unknown
            unit suffixes.
    """
    if isinstance(raw, (int, float)) and not isinstance(raw, bool):
        if raw <= 0 or not math.isfinite(float(raw)):
            raise KdParseError(f"K_d must be > 0 and finite; got {raw!r}")
        return float(raw)
    if not isinstance(raw, str):
        raise KdParseError(
            f"K_d must be a string or numeric; got {type(raw).__name__}"
        )
    s = raw.strip().replace(",", "")
    if not s:
        raise KdParseError("K_d string is empty")
    m = _KD_RE.match(s)
    if not m:
        raise KdParseError(f"could not parse K_d string {raw!r}")
    num_s, unit_raw = m.group(1), (m.group(2) or "")
    try:
        num = float(num_s)
    except ValueError as exc:
        raise KdParseError(
            f"K_d numeric component {num_s!r} not a float (input {raw!r})"
        ) from exc
    if num <= 0 or not math.isfinite(num):
        raise KdParseError(
            f"K_d must be > 0 and finite; got {num!r} (input {raw!r})"
        )
    unit_key = unit_raw.lower().strip()
    # Normalise Greek mu micro (μ) to ``u``.
    unit_key = unit_key.replace("μ", "u")
    if unit_key not in _KD_UNIT_TABLE:
        raise KdParseError(
            f"unknown K_d unit {unit_raw!r} in {raw!r}; "
            f"supported: {sorted(_KD_UNIT_TABLE)}"
        )
    return num * _KD_UNIT_TABLE[unit_key]


def _kd_to_ddg(kd_molar: float, *, T: float = DEFAULT_TEMPERATURE_K) -> float:
    """Convert K_d (in molar) to ΔG_bind (kcal/mol) at temperature ``T``.

    Canonical biochemistry relation::

        K_a = 1 / K_d              (association constant)
        ΔG_bind = -R T ln(K_a) = +R T ln(K_d / 1 M)

    Sign convention: ΔG_bind < 0 for tight binders (K_d ≪ 1 M, the
    common regime for peptide/protein binding). Example: K_d = 80 nM at
    298.15 K → ΔG_bind ≈ -9.68 kcal/mol.

    Notes:
        - ``R = 1.987204e-3 kcal mol⁻¹ K⁻¹`` (CODATA 2018-derived).
        - K_d must be in absolute molar; do **not** pre-divide.
        - Equivalent form: ``ΔG = -R T ln(1 / K_d)`` (literal binding
          equilibrium with K_a as the equilibrium constant).
    """
    if kd_molar <= 0 or not math.isfinite(kd_molar):
        raise KdParseError(f"K_d must be > 0 and finite; got {kd_molar!r}")
    if T <= 0 or not math.isfinite(T):
        raise CycleManagerError(f"temperature must be > 0 K; got {T!r}")
    return R_KCAL_PER_MOL_K * T * math.log(kd_molar)


def _aggregate_assays(results: List[WetLabResult]) -> WetLabResult:
    """Aggregate multiple wet-lab readouts for the same sequence.

    Algorithm:
        - If every record has ``confidence is None`` (or all-zero) →
          arithmetic mean of K_d values (and the corresponding ΔG).
        - If any record has positive confidence → confidence-weighted mean.
        - If the total confidence weight is zero (all confidences provided
          but all zero) → :class:`KdParseError` (caller misconfiguration).

    The aggregate's ``assay`` field becomes a comma-joined sorted unique
    list (e.g. ``"ITC+SPR"``).

    Returns a fresh :class:`WetLabResult` with ``confidence`` carrying the
    summed weight (informational; not re-applied if aggregated again).
    """
    if not results:
        raise KdParseError("cannot aggregate empty assay list")
    seq_set = {r.sequence for r in results}
    if len(seq_set) > 1:
        raise KdParseError(
            f"multi-assay aggregation requires identical sequences; "
            f"got {sorted(seq_set)}"
        )
    confs: List[Optional[float]] = [
        (None if r.confidence is None else float(r.confidence))
        for r in results
    ]
    has_explicit = any(c is not None for c in confs)
    has_positive = any((c is not None and c > 0) for c in confs)
    if has_explicit and not has_positive:
        # All confidences provided but all zero → user misconfiguration.
        raise KdParseError(
            "multi-assay aggregation: all confidences are zero (or non-finite); "
            "remove confidence values to fall back to simple mean, or supply "
            "at least one positive weight"
        )
    if has_positive:
        # Confidence-weighted mean. Records with confidence=None are excluded
        # from the weighted average (only positively-weighted records count).
        total_w = sum(c for c in confs
                      if c is not None and c > 0 and math.isfinite(c))
        if total_w <= 0 or not math.isfinite(total_w):
            raise KdParseError(
                "multi-assay aggregation: total confidence weight is zero "
                "(all positive weights summed to zero or non-finite)"
            )
        kd = sum(
            r.kd_molar * (c if (c is not None and c > 0
                                and math.isfinite(c)) else 0.0)
            for r, c in zip(results, confs)
        ) / total_w
    else:
        # Simple unweighted mean.
        kd = sum(r.kd_molar for r in results) / len(results)
        total_w = float(len(results))
    if kd <= 0 or not math.isfinite(kd):
        raise KdParseError(
            f"multi-assay aggregation produced non-finite/non-positive K_d: {kd!r}"
        )
    T_avg = sum(r.temperature_k for r in results) / len(results)
    ddg = _kd_to_ddg(kd, T=T_avg)
    assay_join = "+".join(sorted({r.assay for r in results if r.assay}))
    note_join = "; ".join(
        f"{r.assay}: K_d={r.kd_molar:.3e}M (conf={r.confidence})"
        for r in results
    )
    return WetLabResult(
        sequence=results[0].sequence,
        kd_molar=kd,
        assay=assay_join or "aggregate",
        ddg_calc=ddg,
        confidence=total_w,
        temperature_k=T_avg,
        notes=f"aggregate of {len(results)} assays — {note_join}",
    )


def _pick_best_validated(cycle: CycleEntry,
                         criterion: str = "lowest_kd") -> Optional[str]:
    """Pick the best-validated sequence from a cycle's wet-lab results.

    Criterion:
        ``"lowest_kd"`` (default) — minimum K_d_molar (= tightest binder).

    Returns the sequence string, or ``None`` if no wet-lab data exist.
    Aggregates multiple assays per sequence first via
    :func:`_aggregate_assays` so that the comparison is on a single K_d
    per sequence.
    """
    if not cycle.wet_lab_results:
        return None
    grouped: Dict[str, List[WetLabResult]] = {}
    for r in cycle.wet_lab_results:
        grouped.setdefault(r.sequence, []).append(r)
    best_seq: Optional[str] = None
    best_kd = math.inf
    for seq, recs in grouped.items():
        if len(recs) == 1:
            kd = recs[0].kd_molar
        else:
            try:
                kd = _aggregate_assays(recs).kd_molar
            except KdParseError:
                # Fall back to arithmetic mean of K_d values.
                kd = sum(r.kd_molar for r in recs) / len(recs)
        if criterion == "lowest_kd":
            if kd < best_kd:
                best_kd = kd
                best_seq = seq
        else:
            raise CycleManagerError(
                f"unknown best-validated criterion {criterion!r}; "
                "supported: lowest_kd"
            )
    return best_seq


# ---------------------------------------------------------------------------
# Public API — project lifecycle
# ---------------------------------------------------------------------------

def init_project(project_id: str, target_pdb: str, *,
                 output_root: str = "outputs/projects",
                 overwrite: bool = False,
                 **metadata: Any) -> CycleHistory:
    """Create a new project under ``<output_root>/<project_id>/``.

    Args:
        project_id: Free-form project identifier (used as directory name).
        target_pdb: Path / identifier of the target PDB (informational).
        output_root: Root directory under which projects are persisted.
        overwrite: If True, replace any existing project state. Default
            False raises :class:`InvalidCycleError` on conflict.
        **metadata: Free-form key/value pairs stored in
            ``CycleHistory.metadata`` (e.g. ``ncaa_budget=5``,
            ``operator='user'``, ``notes='1EBP design 2026Q2'``).

    Returns:
        The newly created :class:`CycleHistory`.
    """
    path = _history_path(project_id, output_root)
    if path.is_file() and not overwrite:
        raise InvalidCycleError(
            f"project {project_id!r} already exists at {path}; "
            "pass overwrite=True to replace"
        )
    history = CycleHistory(
        schema=SCHEMA_VERSION,
        project_id=project_id,
        target_pdb=target_pdb,
        created_at=_now_iso(),
        metadata=dict(metadata),
        cycles=[],
        trajectory={},
        output_root=output_root,
    )
    _save_history(history)
    return history


def add_cycle(project_id: str, cycle_n: int, *,
              wt_sequence: Optional[str] = None,
              variant_sequences: Optional[List[str]] = None,
              predicted_ddg: Optional[List[Union[PredictedDDG, Dict[str, Any]]]] = None,
              notes: Optional[str] = None,
              synthesized: bool = False,
              output_root: str = "outputs/projects") -> CycleEntry:
    """Register a new cycle on the project.

    Args:
        project_id: Project to mutate.
        cycle_n: 0-indexed cycle number. Must be unique within the project.
        wt_sequence: WT reference (``None`` for Cycle 0 de novo).
        variant_sequences: Sequences to be evaluated this cycle. Free-form
            strings; ncAA round-trip via PR-NEW-A is the consumer's
            responsibility (this module stores raw strings).
        predicted_ddg: Optional list of :class:`PredictedDDG` records (or
            dicts; will be normalised). Typically populated after
            PR-NEW-B :func:`utils.branched_ddg.compute_branched_ddg`
            completes.
        notes: Free-form annotation.
        synthesized: Whether variants were physically synthesised.
        output_root: Persistence root.

    Raises:
        InvalidCycleError: on duplicate ``cycle_n`` or negative ``cycle_n``.
        ProjectNotFoundError: on missing project.
    """
    if cycle_n < 0:
        raise InvalidCycleError(f"cycle_n must be >= 0; got {cycle_n}")
    history = _load_history(project_id, output_root=output_root)
    for c in history.cycles:
        if c.cycle_n == cycle_n:
            raise InvalidCycleError(
                f"cycle {cycle_n} already exists in project {project_id!r}"
            )
    pdds: List[PredictedDDG] = []
    for p in (predicted_ddg or []):
        if isinstance(p, PredictedDDG):
            pdds.append(p)
        elif isinstance(p, dict):
            pdds.append(PredictedDDG.from_dict(p))
        else:
            raise CycleManagerError(
                f"predicted_ddg entry must be PredictedDDG or dict; "
                f"got {type(p).__name__}"
            )
    entry = CycleEntry(
        cycle_n=cycle_n,
        wt_sequence=wt_sequence,
        variants=list(variant_sequences or []),
        predicted_ddg=pdds,
        wet_lab_results=[],
        best_validated_sequence=None,
        synthesized=bool(synthesized),
        notes=notes,
        created_at=_now_iso(),
    )
    history.cycles.append(entry)
    history.cycles.sort(key=lambda c: c.cycle_n)
    _save_history(history)
    return entry


def add_result(project_id: str, cycle_n: int, sequence: str,
               wet_lab: Optional[WetLabResult] = None, *,
               kd: Union[str, float, None] = None,
               assay: str = "unknown",
               confidence: Optional[float] = None,
               temperature_k: float = DEFAULT_TEMPERATURE_K,
               notes: Optional[str] = None,
               output_root: str = "outputs/projects") -> WetLabResult:
    """Record a wet-lab assay readout for a sequence in a specific cycle.

    Either pass a fully-built :class:`WetLabResult` via ``wet_lab`` OR pass
    individual fields (``kd`` accepts the same string forms as
    :func:`_parse_kd`).

    Side-effect: ``ddg_calc`` is computed if the caller did not set it.

    Raises:
        ProjectNotFoundError: on missing project.
        InvalidCycleError: if ``cycle_n`` not registered.
        KdParseError: on K_d parse failure.
    """
    history = _load_history(project_id, output_root=output_root)
    target_cycle: Optional[CycleEntry] = None
    for c in history.cycles:
        if c.cycle_n == cycle_n:
            target_cycle = c
            break
    if target_cycle is None:
        raise InvalidCycleError(
            f"cycle {cycle_n} not registered in project {project_id!r}; "
            "call add_cycle first"
        )
    if wet_lab is None:
        if kd is None:
            raise CycleManagerError(
                "add_result requires either wet_lab=WetLabResult or kd=..."
            )
        kd_molar = _parse_kd(kd)
        wet_lab = WetLabResult(
            sequence=sequence,
            kd_molar=kd_molar,
            assay=assay,
            confidence=confidence,
            temperature_k=temperature_k,
            notes=notes,
        )
    if wet_lab.ddg_calc is None:
        wet_lab.ddg_calc = _kd_to_ddg(wet_lab.kd_molar, T=wet_lab.temperature_k)
    target_cycle.wet_lab_results.append(wet_lab)
    _save_history(history)
    return wet_lab


def set_next_wt(project_id: str, source_cycle: int, *,
                best: Union[str, None] = "auto",
                output_root: str = "outputs/projects") -> str:
    """Mark the best-validated sequence in ``source_cycle`` and return it.

    The sequence is stored on the source cycle's
    ``best_validated_sequence`` field. To use it as the next cycle's WT,
    the caller invokes :func:`add_cycle` with
    ``wt_sequence=set_next_wt(...)``.

    Args:
        project_id: Project.
        source_cycle: The cycle from which to pick the best validated
            sequence.
        best: ``"auto"`` (default) → use :func:`_pick_best_validated` with
            criterion ``lowest_kd``. Any other string is taken verbatim
            as the chosen sequence.
        output_root: Persistence root.

    Returns:
        The chosen sequence.

    Raises:
        ProjectNotFoundError: on missing project.
        InvalidCycleError: if ``source_cycle`` is unknown or has no
            wet-lab data when ``best=='auto'``.
    """
    history = _load_history(project_id, output_root=output_root)
    target_cycle: Optional[CycleEntry] = None
    for c in history.cycles:
        if c.cycle_n == source_cycle:
            target_cycle = c
            break
    if target_cycle is None:
        raise InvalidCycleError(
            f"source cycle {source_cycle} not registered in project "
            f"{project_id!r}"
        )
    if best == "auto" or best is None:
        chosen = _pick_best_validated(target_cycle, criterion="lowest_kd")
        if chosen is None:
            raise InvalidCycleError(
                f"source cycle {source_cycle} has no wet-lab results — "
                "cannot auto-pick best validated sequence"
            )
    else:
        chosen = str(best)
    target_cycle.best_validated_sequence = chosen
    _save_history(history)
    return chosen


# ---------------------------------------------------------------------------
# Public API — trajectory / plateau / report
# ---------------------------------------------------------------------------

def compute_trajectory(project_id: str, *,
                       output_root: str = "outputs/projects",
                       lookback: int = 3,
                       threshold_kcal: float = 0.5) -> Dict[str, Any]:
    """Compute the project's wet-lab K_d trajectory and detect plateau.

    For each cycle that has at least one wet-lab record, take the
    ``best_validated_sequence`` (or ``_pick_best_validated`` fallback) and
    its aggregated ΔG. The trajectory is the cycle-over-cycle differences:

        ``delta_kd_kcal[i] = ΔG_cycle_{i+1} − ΔG_cycle_i``  (i ≥ 0)
        ``ratio_per_cycle[i] = K_d_{i+1} / K_d_i``

    Plateau: ``detect_plateau(history, lookback, threshold_kcal)``.

    Recommendation:
        ``"continue"``   — clear ongoing improvement, no plateau.
        ``"converged"``  — plateau **and** final ΔG ≤ a stringent target
                           (``-12 kcal/mol`` ≈ K_d ~ 1 nM); informative
                           default.
        ``"try_alternative"`` — plateau but final ΔG worse than ``-9 kcal/mol``
                           (K_d ~ 250 nM): plateaued at a weak binder, the
                           ncAA strategy may need rethinking.

    Returns the trajectory dict and persists it to
    ``history.trajectory``.
    """
    history = _load_history(project_id, output_root=output_root)
    sorted_cycles = sorted(history.cycles, key=lambda c: c.cycle_n)
    per_cycle: List[Dict[str, Any]] = []
    for c in sorted_cycles:
        seq = c.best_validated_sequence
        recs = c.wet_lab_results
        if not recs:
            per_cycle.append({"cycle_n": c.cycle_n, "kd_molar": None,
                              "ddg": None, "sequence": seq})
            continue
        if seq is None:
            seq = _pick_best_validated(c)
        # Aggregate that sequence's wet-lab data (if multiple assays).
        seq_recs = [r for r in recs if r.sequence == seq] if seq else recs
        if len(seq_recs) == 1:
            kd = seq_recs[0].kd_molar
            T = seq_recs[0].temperature_k
        elif len(seq_recs) > 1:
            agg = _aggregate_assays(seq_recs)
            kd = agg.kd_molar
            T = agg.temperature_k
        else:
            kd = None
            T = DEFAULT_TEMPERATURE_K
        ddg = _kd_to_ddg(kd, T=T) if kd is not None else None
        per_cycle.append({
            "cycle_n": c.cycle_n,
            "kd_molar": kd,
            "ddg": ddg,
            "sequence": seq,
        })

    delta_kd_kcal: List[float] = []
    ratio_per_cycle: List[float] = []
    prev: Optional[Dict[str, Any]] = None
    for cur in per_cycle:
        if prev is not None and prev["ddg"] is not None and cur["ddg"] is not None:
            delta_kd_kcal.append(cur["ddg"] - prev["ddg"])
            ratio_per_cycle.append(cur["kd_molar"] / prev["kd_molar"])
        prev = cur

    plateau = detect_plateau(history, lookback=lookback,
                             threshold_kcal=threshold_kcal,
                             _per_cycle_cache=per_cycle)
    final_ddg = next((p["ddg"] for p in reversed(per_cycle)
                      if p["ddg"] is not None), None)
    if not plateau:
        recommendation = "continue"
    elif final_ddg is not None and final_ddg <= -12.0:
        recommendation = "converged"
    elif final_ddg is not None and final_ddg > -9.0:
        recommendation = "try_alternative"
    else:
        recommendation = "converged"

    trajectory = {
        "per_cycle": per_cycle,
        "delta_kd_kcal": delta_kd_kcal,
        "ratio_per_cycle": ratio_per_cycle,
        "plateau_detected": plateau,
        "lookback": lookback,
        "threshold_kcal": threshold_kcal,
        "recommendation": recommendation,
        "final_ddg_kcal": final_ddg,
        "computed_at": _now_iso(),
    }
    history.trajectory = trajectory
    _save_history(history)
    return trajectory


def detect_plateau(history: CycleHistory, *,
                   lookback: int = 3,
                   threshold_kcal: float = 0.5,
                   _per_cycle_cache: Optional[List[Dict[str, Any]]] = None
                   ) -> bool:
    """Return True iff the last ``lookback`` cycle-over-cycle |ΔΔG_wet|
    values are all below ``threshold_kcal``.

    Behaviour:
        - Requires at least ``lookback + 1`` wet-lab-bearing cycles to
          fire (i.e. ``lookback`` deltas). Below that → False.
        - Uses the same per-cycle aggregation logic as
          :func:`compute_trajectory`.
    """
    if _per_cycle_cache is not None:
        per_cycle = _per_cycle_cache
    else:
        per_cycle = []
        for c in sorted(history.cycles, key=lambda c: c.cycle_n):
            recs = c.wet_lab_results
            if not recs:
                per_cycle.append({"cycle_n": c.cycle_n, "ddg": None,
                                  "kd_molar": None})
                continue
            seq = c.best_validated_sequence or _pick_best_validated(c)
            seq_recs = [r for r in recs if r.sequence == seq] if seq else recs
            if len(seq_recs) == 1:
                kd = seq_recs[0].kd_molar
                T = seq_recs[0].temperature_k
            elif len(seq_recs) > 1:
                agg = _aggregate_assays(seq_recs)
                kd = agg.kd_molar
                T = agg.temperature_k
            else:
                kd = None
                T = DEFAULT_TEMPERATURE_K
            per_cycle.append({
                "cycle_n": c.cycle_n,
                "kd_molar": kd,
                "ddg": (_kd_to_ddg(kd, T=T) if kd is not None else None),
            })
    valid = [p for p in per_cycle if p.get("ddg") is not None]
    if len(valid) < lookback + 1:
        return False
    deltas = [valid[i + 1]["ddg"] - valid[i]["ddg"]
              for i in range(len(valid) - 1)]
    last = deltas[-lookback:]
    return all(abs(d) < threshold_kcal for d in last)


# ---------------------------------------------------------------------------
# Public API — markdown report
# ---------------------------------------------------------------------------

def generate_report(project_id: str, *,
                    format: str = "markdown",
                    output_root: str = "outputs/projects",
                    output_path: Optional[str] = None) -> str:
    """Render a markdown summary of the project's cycle history.

    Args:
        project_id: Project to report on.
        format: Currently only ``"markdown"`` (Day-1 scope; HTML/PDF
            deferred).
        output_root: Persistence root.
        output_path: If set, also write the rendered text to this path.

    Returns:
        The rendered text.
    """
    if format != "markdown":
        raise CycleManagerError(
            f"format {format!r} not supported (Day-1: only 'markdown')"
        )
    history = _load_history(project_id, output_root=output_root)
    if not history.trajectory:
        try:
            compute_trajectory(project_id, output_root=output_root)
            history = _load_history(project_id, output_root=output_root)
        except Exception:  # noqa: BLE001
            pass
    text = _render_markdown(history)
    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(text)
    return text


def _render_markdown(history: CycleHistory) -> str:
    """Internal: render :class:`CycleHistory` as a Phase α-flavoured
    markdown report.
    """
    lines: List[str] = []
    lines.append(f"# UPDD Cycle Manager Report — {history.project_id}")
    lines.append("")
    lines.append(f"- **Schema**: `{history.schema}` (`{TOOL_NAME}`)")
    lines.append(f"- **Target PDB**: `{history.target_pdb}`")
    lines.append(f"- **Created**: {history.created_at}")
    if history.metadata:
        lines.append(f"- **Metadata**: {json.dumps(history.metadata, sort_keys=True)}")
    lines.append("")
    lines.append("## Cycle Summary")
    lines.append("")
    if not history.cycles:
        lines.append("*(no cycles registered)*")
        lines.append("")
    else:
        for c in sorted(history.cycles, key=lambda x: x.cycle_n):
            lines.append(f"## Cycle {c.cycle_n}")
            lines.append("")
            wt = c.wt_sequence or "*(de novo — no WT)*"
            lines.append(f"- WT sequence: `{wt}`")
            lines.append(f"- Variants ({len(c.variants)}): "
                         + (", ".join(f"`{v}`" for v in c.variants)
                            if c.variants else "*(none)*"))
            if c.predicted_ddg:
                lines.append("")
                lines.append("### Predicted ΔΔG (PR-NEW-B)")
                lines.append("")
                lines.append("| Sequence | ΔΔG (kcal/mol) | SE | Tier |")
                lines.append("|---|---|---|---|")
                for p in c.predicted_ddg:
                    se_s = (f"{p.se:.3f}" if p.se is not None else "—")
                    lines.append(
                        f"| `{p.sequence}` | {p.ddg:+.3f} | "
                        f"{se_s} | {p.tier or '—'} |"
                    )
            if c.wet_lab_results:
                lines.append("")
                lines.append("### Wet Lab Results")
                lines.append("")
                lines.append("| Sequence | K_d | Assay | ΔG (kcal/mol) | T (K) |")
                lines.append("|---|---|---|---|---|")
                for w in c.wet_lab_results:
                    ddg_s = (f"{w.ddg_calc:+.3f}"
                             if w.ddg_calc is not None else "—")
                    lines.append(
                        f"| `{w.sequence}` | {_fmt_kd(w.kd_molar)} | "
                        f"{w.assay} | {ddg_s} | {w.temperature_k:.2f} |"
                    )
            if c.best_validated_sequence:
                lines.append("")
                lines.append(
                    f"- **Best validated**: `{c.best_validated_sequence}`"
                )
            if c.notes:
                lines.append(f"- Notes: {c.notes}")
            lines.append("")
    if history.trajectory:
        lines.append("## Trajectory")
        lines.append("")
        per = history.trajectory.get("per_cycle") or []
        if per:
            lines.append("| Cycle | Sequence | K_d | ΔG (kcal/mol) |")
            lines.append("|---|---|---|---|")
            for p in per:
                seq = p.get("sequence") or "—"
                kd = p.get("kd_molar")
                ddg = p.get("ddg")
                lines.append(
                    f"| {p['cycle_n']} | `{seq}` | "
                    f"{_fmt_kd(kd) if kd is not None else '—'} | "
                    f"{ddg:+.3f}".rstrip("—") + (
                        " |" if ddg is not None else "— |"
                    )
                )
            lines.append("")
        deltas = history.trajectory.get("delta_kd_kcal") or []
        ratios = history.trajectory.get("ratio_per_cycle") or []
        if deltas:
            lines.append("| Step | ΔΔG_wet (kcal/mol) | K_d ratio |")
            lines.append("|---|---|---|")
            for i, (d, r) in enumerate(zip(deltas, ratios)):
                lines.append(
                    f"| {i}→{i + 1} | {d:+.3f} | {r:.3e} |"
                )
            lines.append("")
        plateau = history.trajectory.get("plateau_detected")
        rec = history.trajectory.get("recommendation")
        lines.append(f"- **Plateau detected**: {plateau}")
        lines.append(f"- **Recommendation**: `{rec}`")
        lines.append("")
    return "\n".join(lines) + "\n"


def _fmt_kd(kd_molar: Optional[float]) -> str:
    """Pretty-print a K_d in molar as the most appropriate sub-unit."""
    if kd_molar is None or not math.isfinite(kd_molar) or kd_molar <= 0:
        return "—"
    # Pick the unit whose mantissa is in [1, 1000).
    units = [("M", 1.0), ("mM", 1.0e-3), ("uM", 1.0e-6),
             ("nM", 1.0e-9), ("pM", 1.0e-12), ("fM", 1.0e-15)]
    for label, mult in units:
        val = kd_molar / mult
        if 1.0 <= val < 1000.0:
            return f"{val:.3g} {label}"
    return f"{kd_molar:.3e} M"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="python -m utils.cycle_manager",
        description=(
            f"UPDD Iterative Cycle Manager (PR-NEW-C; schema {SCHEMA_VERSION})."
        ),
    )
    ap.add_argument("--output-root", default="outputs/projects",
                    help="Root directory for project state. Default "
                         "'outputs/projects'. May be repeated on the "
                         "subcommand for ergonomic flexibility.")

    def _add_output_root(p: argparse.ArgumentParser) -> None:
        # Allow ``--output-root`` after the subcommand too. Stored under a
        # distinct dest so the parent parser's value is preserved when the
        # subcommand flag is absent.
        p.add_argument("--output-root", default=None, dest="output_root_sub",
                       help="Override parent parser's --output-root.")

    sub = ap.add_subparsers(dest="cmd", required=True)

    p_init = sub.add_parser("init", help="Create a new project.")
    p_init.add_argument("--project-id", required=True)
    p_init.add_argument("--target", "--target-pdb", dest="target_pdb",
                        required=True)
    p_init.add_argument("--operator", default=None)
    p_init.add_argument("--ncaa-budget", type=int, default=None)
    p_init.add_argument("--notes", default=None)
    p_init.add_argument("--overwrite", action="store_true")
    _add_output_root(p_init)

    p_add = sub.add_parser("add-cycle",
                           help="Register a new cycle on the project.")
    p_add.add_argument("--project-id", required=True)
    p_add.add_argument("--cycle", type=int, required=True, dest="cycle_n")
    p_add.add_argument("--wt", default=None, dest="wt_sequence",
                       help="WT sequence (omit for Cycle 0 de novo).")
    p_add.add_argument("--variants", default="",
                       help="Comma-separated variant sequences.")
    p_add.add_argument("--notes", default=None)
    p_add.add_argument("--synthesized", action="store_true")
    _add_output_root(p_add)

    p_res = sub.add_parser("add-result",
                           help="Record a wet-lab K_d for a sequence.")
    p_res.add_argument("--project-id", required=True)
    p_res.add_argument("--cycle", type=int, required=True, dest="cycle_n")
    p_res.add_argument("--sequence", required=True)
    p_res.add_argument("--kd", required=True,
                       help="K_d value (e.g. '80nM', '8e-8', '0.08uM').")
    p_res.add_argument("--assay", default="unknown")
    p_res.add_argument("--confidence", type=float, default=None)
    p_res.add_argument("--temperature", type=float,
                       default=DEFAULT_TEMPERATURE_K, dest="temperature_k")
    p_res.add_argument("--notes", default=None)
    _add_output_root(p_res)

    p_wt = sub.add_parser("set-wt",
                          help="Pick best-validated sequence from a cycle.")
    p_wt.add_argument("--project-id", required=True)
    p_wt.add_argument("--source-cycle", type=int, required=True)
    p_wt.add_argument("--best", default="auto",
                      help="'auto' (default) or a literal sequence string.")
    _add_output_root(p_wt)

    p_traj = sub.add_parser("trajectory",
                            help="Compute K_d trajectory + plateau detection.")
    p_traj.add_argument("--project-id", required=True)
    p_traj.add_argument("--lookback", type=int, default=3)
    p_traj.add_argument("--threshold-kcal", type=float, default=0.5)
    _add_output_root(p_traj)

    p_rep = sub.add_parser("report", help="Generate markdown report.")
    p_rep.add_argument("--project-id", required=True)
    p_rep.add_argument("--output", default=None,
                       help="Optional path to write the rendered markdown.")
    _add_output_root(p_rep)

    return ap


def main(argv: Optional[List[str]] = None) -> int:
    ap = _build_parser()
    args = ap.parse_args(argv)
    # If the user supplied --output-root after the subcommand, prefer that.
    sub_root = getattr(args, "output_root_sub", None)
    if sub_root is not None:
        args.output_root = sub_root
    try:
        if args.cmd == "init":
            md: Dict[str, Any] = {}
            if args.operator:
                md["operator"] = args.operator
            if args.ncaa_budget is not None:
                md["ncaa_budget"] = args.ncaa_budget
            if args.notes:
                md["notes"] = args.notes
            history = init_project(
                project_id=args.project_id,
                target_pdb=args.target_pdb,
                output_root=args.output_root,
                overwrite=args.overwrite,
                **md,
            )
            print(f"[cycle_manager] init OK — schema={history.schema} "
                  f"project_id={history.project_id} "
                  f"path={_history_path(args.project_id, args.output_root)}")
            return 0
        if args.cmd == "add-cycle":
            variants = [v.strip() for v in (args.variants or "").split(",")
                        if v.strip()]
            entry = add_cycle(
                project_id=args.project_id,
                cycle_n=args.cycle_n,
                wt_sequence=args.wt_sequence,
                variant_sequences=variants,
                notes=args.notes,
                synthesized=args.synthesized,
                output_root=args.output_root,
            )
            print(f"[cycle_manager] add-cycle OK — cycle_n={entry.cycle_n} "
                  f"variants={len(entry.variants)} "
                  f"wt={entry.wt_sequence!r}")
            return 0
        if args.cmd == "add-result":
            wlr = add_result(
                project_id=args.project_id,
                cycle_n=args.cycle_n,
                sequence=args.sequence,
                kd=args.kd,
                assay=args.assay,
                confidence=args.confidence,
                temperature_k=args.temperature_k,
                notes=args.notes,
                output_root=args.output_root,
            )
            print(f"[cycle_manager] add-result OK — sequence={wlr.sequence} "
                  f"K_d={_fmt_kd(wlr.kd_molar)} ΔG={wlr.ddg_calc:+.3f} "
                  f"assay={wlr.assay}")
            return 0
        if args.cmd == "set-wt":
            chosen = set_next_wt(
                project_id=args.project_id,
                source_cycle=args.source_cycle,
                best=args.best,
                output_root=args.output_root,
            )
            print(f"[cycle_manager] set-wt OK — chosen={chosen!r}")
            return 0
        if args.cmd == "trajectory":
            traj = compute_trajectory(
                project_id=args.project_id,
                output_root=args.output_root,
                lookback=args.lookback,
                threshold_kcal=args.threshold_kcal,
            )
            print(f"[cycle_manager] trajectory OK — "
                  f"plateau={traj['plateau_detected']} "
                  f"recommendation={traj['recommendation']!r}")
            print(f"[cycle_manager]   delta_kd_kcal={traj['delta_kd_kcal']}")
            print(f"[cycle_manager]   ratio_per_cycle={traj['ratio_per_cycle']}")
            return 0
        if args.cmd == "report":
            text = generate_report(
                project_id=args.project_id,
                output_root=args.output_root,
                output_path=args.output,
            )
            print(text)
            return 0
    except ProjectNotFoundError as exc:
        print(f"[cycle_manager] ProjectNotFoundError: {exc}", file=sys.stderr)
        return 4
    except InvalidCycleError as exc:
        print(f"[cycle_manager] InvalidCycleError: {exc}", file=sys.stderr)
        return 5
    except KdParseError as exc:
        print(f"[cycle_manager] KdParseError: {exc}", file=sys.stderr)
        return 6
    except CycleManagerError as exc:
        print(f"[cycle_manager] ERROR: {type(exc).__name__}: {exc}",
              file=sys.stderr)
        return 2
    print(f"[cycle_manager] unknown subcommand {args.cmd!r}", file=sys.stderr)
    return 3


if __name__ == "__main__":
    sys.exit(main())
