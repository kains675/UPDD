"""Track B adaptive λ-scheduling engine (MVP — advisory, one-shot, numpy-only).

This package implements the **MVP** of the automated adaptive
λ-scheduling design for Track B ATM ABFE/RBFE — a GATED PROPOSER / DIAGNOSTIC,
NOT an autonomous self-deployer. It DIAGNOSEs per-adjacent-pair overlap from an
existing per-direction pilot's ``.out`` files, JUDGEs bottleneck / redundant
intervals, and PRESCRIBEs a count-fixed, endpoint-pinned, linear-region-only
rebalance — the exact generalization of the manual densified38 → densified38v2
move.

Scope (MVP, this package):
  * ``overlap.py``    — DIAGNOSE (Bhattacharyya gate proxy + optional MBAR).
  * ``judge.py``      — JUDGE (bottleneck / redundant / region classification).
  * ``prescribe.py``  — PRESCRIBE (count-fixed linear-region rebalance,
                        delegates ladder math to the EXISTING
                        ``_build_densified38_arrays``).
  * ``constraints.py``— GUARDS (load-bearing fail-loud invariance asserts).
  * ``schedule_io.py``— registry-dict round-trip + validated JSON emit.

Explicitly OUT OF SCOPE for the MVP (Stage 2+, separate work):
  * ``engine.py``     — iterate-to-convergence state machine, MBAR gate,
                        bootstrap CI, damping, additive count-growth fallback.
  * any auto-launch / subprocess pilot loop.
  * cross-region (linear ↔ soft-core) λ transfer.

Regime: ranking-only. The redistribution is ΔG-unbiased BY CONSTRUCTION
**iff** the constraints in ``constraints.py`` are code-enforced — this is the
load-bearing correctness condition for the whole package.

Hard deps: numpy only. pymbar is OPTIONAL (auto-detected); when absent the
Bhattacharyya coefficient remains the gate (MVP must NOT hard-depend on pymbar).
"""

from .constraints import (
    AdaptiveConfig,
    AdaptiveConstraintError,
    EndpointInvariantError,
    RegionLockedError,
    SymmetryMirrorError,
    CountPolicyError,
    endpoint_invariant,
    region_locked,
    symmetry_mirror,
    count_policy,
)
from .overlap import (
    PairOverlap,
    OverlapReport,
    extract_pertE_by_state,
    extract_samples_by_state,
    bhattacharyya_coefficient,
    histogram_intersection_advisory,
    adjacent_overlaps,
    build_softcore_neg_pot,
    mbar_overlap_matrix,
    estimate_overlap,
)
from .judge import JudgeResult, classify
from .prescribe import propose_rebalance
from .schedule_io import (
    to_registry_dict,
    from_registry_dict,
    emit_validated_schedule,
    load_schedule_dict,
)

__all__ = [
    # constraints
    "AdaptiveConfig",
    "AdaptiveConstraintError",
    "EndpointInvariantError",
    "RegionLockedError",
    "SymmetryMirrorError",
    "CountPolicyError",
    "endpoint_invariant",
    "region_locked",
    "symmetry_mirror",
    "count_policy",
    # overlap
    "PairOverlap",
    "OverlapReport",
    "extract_pertE_by_state",
    "extract_samples_by_state",
    "bhattacharyya_coefficient",
    "histogram_intersection_advisory",
    "adjacent_overlaps",
    "build_softcore_neg_pot",
    "mbar_overlap_matrix",
    "estimate_overlap",
    # judge
    "JudgeResult",
    "classify",
    # prescribe
    "propose_rebalance",
    # schedule_io
    "to_registry_dict",
    "from_registry_dict",
    "emit_validated_schedule",
    "load_schedule_dict",
]
