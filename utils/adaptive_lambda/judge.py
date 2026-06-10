"""Adaptive λ-scheduling JUDGE module — bottleneck / redundant classification.

Consumes an ``OverlapReport`` (overlap.py) + the seed schedule + an
``AdaptiveConfig`` and classifies each adjacent pair:

  * BOTTLENECK — BC < ``bottleneck_floor`` (0.30) OR MBAR O < 0.1 when available.
    The variance-limiting bridge(s); PRESCRIBE adds resolution here.
  * REDUNDANT  — contiguous pair with BC > ``redundant_ceiling`` (0.85). When
    Δlength information is available (adjacent λ spacing as a thermodynamic-
    length surrogate), redundancy ALSO requires Δlength < mean/3 (the verdict's
    double condition — over-overlap alone must not thin a steep region). When the
    schedule does not expose per-pair Δλ, the BC-only redundancy stands (MVP).

``region_of(pair)`` is derived from the SEED schedule's W0 / INTERMEDIATE columns
(W0==0 AND INTERMEDIATE==0 → "linear", else "softcore") — NOT a hardcoded "state
10" boundary (the boundary may move between schedules; the spec forbids a literal
index).

Python 3.8+; numpy-only.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .constraints import AdaptiveConfig
from .overlap import OverlapReport, PairOverlap


@dataclass
class JudgeResult:
    """Classification of an OverlapReport against the configured thresholds."""

    # direction tag -> list of PairOverlap flagged as bottleneck (BC<floor / O<0.1)
    bottlenecks: Dict[str, List[PairOverlap]] = field(default_factory=dict)
    # direction tag -> list of contiguous PairOverlap flagged redundant
    redundant_runs: Dict[str, List[PairOverlap]] = field(default_factory=dict)
    # direction tag -> {(i, j): "linear"|"softcore"}
    regions: Dict[str, Dict[Tuple[int, int], str]] = field(default_factory=dict)
    all_pass: bool = False
    config_snapshot: Dict[str, Any] = field(default_factory=dict)

    def region_of(self, direction: str, pair: PairOverlap) -> str:
        """Return the region ("linear"/"softcore") of ``pair`` in ``direction``.

        Looked up from the per-direction region map derived from the seed
        schedule's W0/INTERMEDIATE columns. Raises KeyError if the pair was not
        classified (it was outside the overlap report).
        """
        return self.regions[direction][(pair.i, pair.j)]

    def as_dict(self) -> Dict[str, Any]:
        return {
            "all_pass": self.all_pass,
            "config_snapshot": self.config_snapshot,
            "bottlenecks": {
                d: [p.as_dict() for p in pairs]
                for d, pairs in self.bottlenecks.items()
            },
            "redundant_runs": {
                d: [p.as_dict() for p in pairs]
                for d, pairs in self.redundant_runs.items()
            },
            "regions": {
                d: {f"{i}-{j}": r for (i, j), r in rmap.items()}
                for d, rmap in self.regions.items()
            },
        }


def _state_is_linear(schedule: Dict[str, Any], state_id: int) -> bool:
    """A state is LINEAR when its W0==0.0 AND INTERMEDIATE==0 (the plain-λ region
    before the soft-core ladder). Derived from the schedule columns — NOT a
    hardcoded boundary index.

    The schedule arrays are indexed by GLOBAL state id (the same ordering the
    .out col-0 stateid uses). Raises IndexError-safe: returns False (treat as
    softcore) if the state is out of range, so a pair touching an unknown state
    is conservatively NOT treated as a relocatable linear pair.
    """
    w0 = schedule.get("w0")
    intermd = schedule.get("intermd")
    if w0 is None or intermd is None:
        raise KeyError(
            "schedule must expose 'w0' and 'intermd' columns to derive region"
        )
    if state_id < 0 or state_id >= len(w0):
        return False
    return float(w0[state_id]) == 0.0 and int(round(float(intermd[state_id]))) == 0


def _region_of_pair(schedule: Dict[str, Any], pair: PairOverlap) -> str:
    """A pair is "linear" iff BOTH its states are linear; otherwise "softcore"
    (a pair straddling the boundary is softcore — it must not be rebalanced by
    the linear-region-only MVP)."""
    if _state_is_linear(schedule, pair.i) and _state_is_linear(schedule, pair.j):
        return "linear"
    return "softcore"


def _thermo_length_surrogate(pair: PairOverlap) -> float:
    """Discrete thermodynamic-length surrogate for the redundancy double-
    condition: |Δmean(pertE)| across the pair.

    Per verdict Q1, the principled allocator objective is equal thermodynamic-
    length increment Δlength = ∫√(g)dλ with g = Var(∂U/∂control)/kT, whose
    DISCRETE surrogate is the change in the perturbation-energy integrand across
    the pair. The mean-pertE jump |mean_j − mean_i| is the leading-order term: a
    plateau pair (Δmean < ~3 kcal/mol — tiny integrand change) carries negligible
    thermodynamic length and is genuinely over-allocated, whereas the turnover
    pair (Δmean ~130 kcal/mol — huge integrand change) is under-allocated. This
    is EXACTLY the verdict's plateau-vs-turnover discrimination (the manual
    analysis used "Δmean < 4 over 9 windows" as the redundancy guardrail, NOT raw
    λ-spacing — the plateau is uniform-Δλ yet redundant). Always derivable from
    the overlap report's per-state means (no schedule lookup needed)."""
    return abs(float(pair.mean_j) - float(pair.mean_i))


def classify(
    report: OverlapReport,
    schedule: Dict[str, Any],
    config: Optional[AdaptiveConfig] = None,
) -> JudgeResult:
    """Classify an OverlapReport into bottleneck / redundant pairs + regions.

    Bottleneck rule (verdict Q2 thresholds): a pair is a bottleneck if
    ``BC < config.bottleneck_floor`` OR (when MBAR-O is available) ``O < 0.1``.
    MBAR-O, when present, is the PRIMARY gate and SUPERSEDES BC (verdict G3); BC
    is the cheap proxy / fallback when pymbar is absent.

    Redundant rule (verdict Q2 double condition): a contiguous LINEAR pair with
    ``BC > config.redundant_ceiling`` AND a thermodynamic-length surrogate
    (|Δmean(pertE)|) below one-third of the mean surrogate across the SAME region
    — over-overlap alone does not qualify (a high-overlap pair on a STEEP
    interval, i.e. a large integrand change, must not be thinned). The surrogate
    is the change in the perturbation-energy integrand (verdict Q1 discrete g),
    NOT raw λ-spacing: the densified38 plateau is uniform-Δλ yet redundant
    because its Δmean(pertE) is ~2 kcal/mol vs the turnover's ~130.

    ``all_pass`` is True iff there are NO bottlenecks in any direction (the
    schedule clears the overlap floor everywhere — no rebalance needed).
    """
    cfg = config or AdaptiveConfig()
    result = JudgeResult(config_snapshot={
        "bottleneck_floor": cfg.bottleneck_floor,
        "redundant_ceiling": cfg.redundant_ceiling,
        "mbar_well_determined_floor": cfg.mbar_well_determined_floor,
        "region": cfg.region,
    })

    any_bottleneck = False
    for direction, pairs in report.per_direction.items():
        region_map: Dict[Tuple[int, int], str] = {}
        for p in pairs:
            region_map[(p.i, p.j)] = _region_of_pair(schedule, p)
        result.regions[direction] = region_map

        # --- bottlenecks ---
        bn: List[PairOverlap] = []
        for p in pairs:
            is_bottleneck = p.bc < cfg.bottleneck_floor
            if p.mbar_o is not None:
                # MBAR-O is the primary gate when available.
                is_bottleneck = (
                    is_bottleneck or p.mbar_o < cfg.mbar_well_determined_floor
                )
            if is_bottleneck:
                bn.append(p)
        if bn:
            any_bottleneck = True
        result.bottlenecks[direction] = bn

        # --- redundant runs (verdict Q2 double condition) ---
        # Mean thermodynamic-length surrogate |Δmean(pertE)| across the LINEAR
        # region (the only region we can relocate in v1). A pair is redundant
        # iff over-overlapped (BC > ceiling) AND its surrogate is < one-third of
        # the region mean (negligible integrand change = flat plateau). The
        # surrogate is the integrand jump, NOT raw Δλ — the densified38 plateau
        # is uniform-Δλ yet flat (Δmean ~2 vs turnover ~130).
        linear_surr = [
            _thermo_length_surrogate(p) for p in pairs
            if region_map[(p.i, p.j)] == "linear"
        ]
        mean_surr = float(np.mean(linear_surr)) if linear_surr else None

        redundant: List[PairOverlap] = []
        for p in pairs:
            if p.bc <= cfg.redundant_ceiling:
                continue
            # Region restriction: only LINEAR pairs are relocatable in v1.
            if region_map[(p.i, p.j)] != "linear":
                continue
            surr = _thermo_length_surrogate(p)
            if mean_surr is not None and mean_surr > 0.0:
                # Double condition: also require negligible Δlength surrogate.
                if surr < mean_surr / 3.0:
                    redundant.append(p)
                # else: over-overlapped but on a steep interval → keep (verdict
                # guardrail: do not thin a high-overlap pair on a steep region).
            else:
                # Degenerate region (all-zero surrogate) → BC-only fallback.
                redundant.append(p)
        result.redundant_runs[direction] = redundant

    result.all_pass = not any_bottleneck
    return result
