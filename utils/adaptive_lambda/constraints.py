"""Adaptive λ-scheduling GUARDS — LOAD-BEARING, fail-loud.

These assertions are the single most important correctness requirement of
this package. The ΔΔG-unbias guarantee (ranking-safe) of an interior-window
rebalance is valid **only if** endpoint-invariance + within-region-locking +
reverse-symmetry + count-policy are CODE-ENFORCED. If any of these can be
silently violated, the rebalance is not provably unbiased.

Every guard RAISES a typed ``AdaptiveConstraintError`` subclass on violation —
never returns a soft False, never warns-and-continues. A proposed schedule that
trips any guard must be discarded by the caller (prescribe.py asserts all four
before returning a proposal).

Python 3.8+ compatible (Optional/List/Dict/Tuple from typing). numpy-only.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Exception hierarchy
# ---------------------------------------------------------------------------
class AdaptiveConstraintError(Exception):
    """Base class for a load-bearing constraint violation."""


class EndpointInvariantError(AdaptiveConstraintError):
    """Raised when a physical endpoint (λ=0 coupled OR apex W0=1/λ=0.5) moved."""


class RegionLockedError(AdaptiveConstraintError):
    """Raised when a non-linear (soft-core ladder) state was mutated."""


class SymmetryMirrorError(AdaptiveConstraintError):
    """Raised when the backward half is not the whole-tuple reverse of forward."""


class CountPolicyError(AdaptiveConstraintError):
    """Raised when the state count changed under a ``fixed`` count policy."""


# ---------------------------------------------------------------------------
# Config — decision thresholds the gate depends on. Do NOT edit these without
# a documented scientific review.
# ---------------------------------------------------------------------------
@dataclass
class AdaptiveConfig:
    """Adaptive λ-scheduling parameters (MVP subset + Stage-2 keys).

    The MVP (one-shot, advisory) uses ``bottleneck_floor`` / ``redundant_ceiling``
    / ``count_policy`` / ``damping_max_moves`` / ``region``. The remaining keys
    (``target_margin`` / ``max_iterations`` / ``pilot_max_samples``) are recorded
    here so Stage-2's engine.py inherits the SAME source-of-truth values without
    re-deriving them.
    """

    # Bhattacharyya / overlap GATE thresholds (the cheap-proxy gate; MBAR-O is
    # the primary gate when computable but BC stays the gate when pymbar absent).
    bottleneck_floor: float = 0.30      # BC < this (or O < 0.1) → bottleneck
    redundant_ceiling: float = 0.85     # BC > this (and Δlength<mean/3) → redundant
    target_margin: float = 0.40         # post-fix turnover target BC (margin)

    # Engine iteration (Stage 2 — recorded for SSOT, unused in MVP one-shot).
    max_iterations: int = 4

    # Rebalance policy.
    count_policy: str = "fixed"         # "fixed" (MVP) | "additive" (Stage-2 gated)
    damping_max_moves: int = 2          # max linear states relocated per proposal

    # Pilot sizing (Stage 2 — recorded for SSOT).
    pilot_max_samples: int = 25

    # Region the MVP rebalances within (linear λ segment only; cross-region is a
    # potential-shape re-parameterization needing a separate verdict).
    region: str = "linear"

    # Hard estimability floor for MBAR overlap (advisory in MVP; Stage-2 gate).
    mbar_estimability_floor: float = 0.03
    mbar_well_determined_floor: float = 0.10

    # Numerical tolerance for endpoint / apex float comparisons.
    endpoint_atol: float = 1e-9

    def __post_init__(self) -> None:
        if self.count_policy not in ("fixed", "additive"):
            raise ValueError(
                f"count_policy must be 'fixed' or 'additive', got "
                f"{self.count_policy!r}"
            )
        if self.region not in ("linear", "softcore"):
            raise ValueError(
                f"region must be 'linear' or 'softcore', got {self.region!r}"
            )
        if not (0.0 < self.bottleneck_floor < self.redundant_ceiling < 1.0):
            raise ValueError(
                "thresholds must satisfy 0 < bottleneck_floor < "
                "redundant_ceiling < 1"
            )


# ---------------------------------------------------------------------------
# Internal helpers — extract the per-state arrays from a schedule dict in the
# canonical _schedule_dict shape (lambdas / lambdas_1 / lambdas_2 / directions /
# intermd / w0 / alpha / u0 / n_states).
# ---------------------------------------------------------------------------
_REQUIRED_KEYS = (
    "lambdas_1", "lambdas_2", "directions", "intermd", "w0", "alpha", "u0",
)


def _require_schedule_shape(schedule: Dict[str, Any], label: str) -> None:
    """Fail loud if the schedule dict is missing required per-state arrays or
    has inconsistent array lengths."""
    missing = [k for k in _REQUIRED_KEYS if k not in schedule]
    if missing:
        raise AdaptiveConstraintError(
            f"{label}: schedule dict missing keys {missing} "
            f"(have {sorted(schedule)})"
        )
    n = len(schedule["directions"])
    for k in _REQUIRED_KEYS:
        if len(schedule[k]) != n:
            raise AdaptiveConstraintError(
                f"{label}: array {k!r} length {len(schedule[k])} != "
                f"n_states {n}"
            )
    if "lambdas" in schedule and len(schedule["lambdas"]) != n:
        raise AdaptiveConstraintError(
            f"{label}: 'lambdas' length {len(schedule['lambdas'])} != "
            f"n_states {n}"
        )


def _forward_count(directions: List[int]) -> int:
    """Number of leading DIRECTION>=0 (forward) states. Mirrors
    ``_derive_state_counts_from_directions`` semantics (contiguous +1 block
    then -1 block) but returns only the forward count.
    """
    fwd = sum(1 for d in directions if float(d) >= 0)
    return fwd


def _is_linear_state(w0: float, intermd: Any) -> bool:
    """A LINEAR-segment state has W0==0.0 AND INTERMEDIATE==0 (the plain-λ region
    before the soft-core ladder). Anything with W0!=0 or INTERMEDIATE!=0 is a
    soft-core ladder state (verdict Q4; densified38 ladder starts at state 10).
    """
    return float(w0) == 0.0 and int(round(float(intermd))) == 0


def _apex_signature(schedule: Dict[str, Any], idx: int) -> Tuple[float, ...]:
    """The byte-significant tuple identifying the apex / endpoint state.

    Includes every alchemical control parameter so a change to ANY of them
    (λ1, λ2, W0, ALPHA, U0, INTERMEDIATE) trips the guard — the apex is the
    fully-decoupled physical state and ΔG depends on it.
    """
    return (
        float(schedule["lambdas_1"][idx]),
        float(schedule["lambdas_2"][idx]),
        float(schedule["w0"][idx]),
        float(schedule["alpha"][idx]),
        float(schedule["u0"][idx]),
        float(schedule["intermd"][idx]),
    )


# ---------------------------------------------------------------------------
# GUARD 1 — endpoint invariant (λ=0 coupled + apex W0=1/λ=0.5 byte-identical)
# ---------------------------------------------------------------------------
def endpoint_invariant(
    seed: Dict[str, Any],
    proposed: Dict[str, Any],
    config: Optional[AdaptiveConfig] = None,
) -> bool:
    """Assert both physical endpoints are invariant between ``seed`` and
    ``proposed``. Raises ``EndpointInvariantError`` on any change.

    Endpoints (per direction):
      * λ=0 coupled — the FIRST forward state's λ1 must be 0.0 and unchanged.
      * apex — the W0==1.0 / λ2==0.5 fully-decoupled state must be byte-identical
        in EVERY control parameter (the apex_signature). The apex is located by
        ``w0 == 1.0`` (not by a hardcoded index) so the guard is schedule-shape
        agnostic.

    ΔG = G(apex) − G(λ=0) is a state function (Zwanzig 1954; Kirkwood 1935); it
    depends ONLY on these endpoints. Holding them invariant is what makes the
    interior rebalance ΔΔG-unbiased (ranking-safe). This guard is the empirical
    enforcement of that proof.
    """
    cfg = config or AdaptiveConfig()
    atol = cfg.endpoint_atol
    _require_schedule_shape(seed, "endpoint_invariant.seed")
    _require_schedule_shape(proposed, "endpoint_invariant.proposed")

    # --- λ=0 coupled endpoint (forward first state) ---
    if not np.isclose(float(seed["lambdas_1"][0]), 0.0, atol=atol):
        raise EndpointInvariantError(
            f"seed forward first state λ1={seed['lambdas_1'][0]} is not the "
            f"coupled endpoint λ=0.0 — schedule does not start coupled"
        )
    if not np.isclose(
        float(proposed["lambdas_1"][0]), float(seed["lambdas_1"][0]), atol=atol
    ):
        raise EndpointInvariantError(
            f"coupled endpoint λ=0 MOVED: seed λ1[0]={seed['lambdas_1'][0]} → "
            f"proposed λ1[0]={proposed['lambdas_1'][0]}"
        )

    # --- apex endpoint (located by W0==1.0, byte-identical in every param) ---
    seed_apex = [i for i, w in enumerate(seed["w0"]) if float(w) == 1.0]
    prop_apex = [i for i, w in enumerate(proposed["w0"]) if float(w) == 1.0]
    if not seed_apex:
        raise EndpointInvariantError(
            "seed has no apex state (no W0==1.0) — cannot verify decoupled "
            "endpoint"
        )
    if len(seed_apex) != len(prop_apex):
        raise EndpointInvariantError(
            f"apex (W0==1.0) state COUNT changed: seed has {len(seed_apex)} "
            f"apex states, proposed has {len(prop_apex)}"
        )
    for si, pi in zip(seed_apex, prop_apex):
        ssig = _apex_signature(seed, si)
        psig = _apex_signature(proposed, pi)
        if not np.allclose(ssig, psig, atol=atol):
            raise EndpointInvariantError(
                f"apex endpoint MUTATED: seed state {si} {ssig} != "
                f"proposed state {pi} {psig} (λ1,λ2,W0,ALPHA,U0,INTERMEDIATE)"
            )
    return True


# ---------------------------------------------------------------------------
# GUARD 2 — region locked (only linear-segment states may change)
# ---------------------------------------------------------------------------
def region_locked(
    seed: Dict[str, Any],
    proposed: Dict[str, Any],
    config: Optional[AdaptiveConfig] = None,
) -> bool:
    """Assert that ONLY linear-segment states (W0==0 AND INTERMEDIATE==0)
    changed; every soft-core ladder state (W0!=0 or INTERMEDIATE!=0) must be
    byte-identical between ``seed`` and ``proposed``. Raises ``RegionLockedError``
    on any ladder mutation.

    The linear ↔ soft-core boundary must not be crossed in v1: the soft-core
    ladder's (λ1,λ2,W0,α,U0) tuple is a designed annealing path; moving windows
    across the boundary is a potential-shape re-parameterization, not a spacing
    move (verdict Q4-ii). The MVP rebalances WITHIN the linear region only.

    Comparison is done on the ladder slice. Because the linear segment may be
    re-spaced (different λ values, possibly different per-direction linear
    COUNT under a future additive policy), we identify ladder states by their
    own (W0,INTERMEDIATE) signature in EACH schedule and compare the ladder
    tuples directionally.
    """
    cfg = config or AdaptiveConfig()
    atol = cfg.endpoint_atol
    _require_schedule_shape(seed, "region_locked.seed")
    _require_schedule_shape(proposed, "region_locked.proposed")

    seed_ladder = _ladder_tuples(seed)
    prop_ladder = _ladder_tuples(proposed)
    if len(seed_ladder) != len(prop_ladder):
        raise RegionLockedError(
            f"soft-core ladder window COUNT changed: seed has "
            f"{len(seed_ladder)} ladder states, proposed has "
            f"{len(prop_ladder)} — cross-region transfer / ladder mutation is "
            f"forbidden in v1 (region={cfg.region!r})"
        )
    for k, (s_tup, p_tup) in enumerate(zip(seed_ladder, prop_ladder)):
        if not np.allclose(s_tup, p_tup, atol=atol):
            raise RegionLockedError(
                f"soft-core ladder window {k} MUTATED: seed {s_tup} != "
                f"proposed {p_tup} (λ1,λ2,W0,ALPHA,U0). The ladder must be "
                f"byte-identical; v1 rebalances the LINEAR segment only."
            )
    return True


def _ladder_tuples(schedule: Dict[str, Any]) -> List[Tuple[float, ...]]:
    """Return the ordered list of soft-core ladder per-state tuples
    (λ1,λ2,W0,ALPHA,U0) for the FORWARD half. A ladder state is W0!=0 or
    INTERMEDIATE!=0. Forward half = leading DIRECTION>=0 block.
    """
    directions = [float(d) for d in schedule["directions"]]
    fwd = _forward_count(directions)
    tuples: List[Tuple[float, ...]] = []
    for i in range(fwd):
        if not _is_linear_state(schedule["w0"][i], schedule["intermd"][i]):
            tuples.append((
                float(schedule["lambdas_1"][i]),
                float(schedule["lambdas_2"][i]),
                float(schedule["w0"][i]),
                float(schedule["alpha"][i]),
                float(schedule["u0"][i]),
            ))
    return tuples


# ---------------------------------------------------------------------------
# GUARD 3 — symmetry mirror (backward == whole-tuple reverse of forward)
# ---------------------------------------------------------------------------
def symmetry_mirror(
    schedule: Dict[str, Any],
    config: Optional[AdaptiveConfig] = None,
) -> bool:
    """Assert the backward half is the whole-tuple reverse of the forward half
    for every per-state array EXCEPT ``directions`` (forward=+1, backward=-1).
    Raises ``SymmetryMirrorError`` on desync.

    The backward (dminus) half is built as ``list(reversed(forward))`` in
    ``_build_densified38_arrays``; the rebalancer edits the forward array ONLY
    and the backward half inherits the mirror automatically (verdict Q4-iv). This
    guard verifies that property held — a desynced backward half means a bug that
    could bias ΔG.

    The forward / backward split is derived from the DIRECTION column (leading
    +1 block, trailing -1 block); the two halves must be equal length.
    """
    cfg = config or AdaptiveConfig()
    atol = cfg.endpoint_atol
    _require_schedule_shape(schedule, "symmetry_mirror")

    directions = [float(d) for d in schedule["directions"]]
    n = len(directions)
    fwd = _forward_count(directions)
    bwd = n - fwd
    if fwd == 0 or bwd == 0:
        raise SymmetryMirrorError(
            f"schedule must have both forward (+1) and backward (-1) states; "
            f"got fwd={fwd} bwd={bwd}"
        )
    if fwd != bwd:
        raise SymmetryMirrorError(
            f"forward / backward halves must be equal length for the mirror "
            f"symmetry; got fwd={fwd} bwd={bwd}"
        )
    # Contiguity: leading block +1, trailing block -1.
    if not (all(d >= 0 for d in directions[:fwd])
            and all(d < 0 for d in directions[fwd:])):
        raise SymmetryMirrorError(
            f"DIRECTION must be a contiguous forward(+1) then backward(-1) "
            f"block; got {schedule['directions']!r}"
        )

    mirror_keys = ("lambdas_1", "lambdas_2", "w0", "alpha", "u0", "intermd")
    if "lambdas" in schedule:
        mirror_keys = mirror_keys + ("lambdas",)
    for key in mirror_keys:
        arr = [float(v) for v in schedule[key]]
        fwd_half = arr[:fwd]
        bwd_half = arr[fwd:]
        if not np.allclose(bwd_half, list(reversed(fwd_half)), atol=atol):
            raise SymmetryMirrorError(
                f"backward half of {key!r} is NOT the whole-tuple reverse of "
                f"the forward half: forward={fwd_half} backward={bwd_half} "
                f"(expected {list(reversed(fwd_half))})"
            )
    # DIRECTION column itself: forward all +1, backward all -1.
    if not (all(d == 1.0 for d in directions[:fwd])
            and all(d == -1.0 for d in directions[fwd:])):
        raise SymmetryMirrorError(
            f"DIRECTION column must be exactly +1 (forward) / -1 (backward); "
            f"got {schedule['directions']!r}"
        )
    return True


# ---------------------------------------------------------------------------
# GUARD 4 — count policy (n_states preserved under "fixed")
# ---------------------------------------------------------------------------
def count_policy(
    seed: Dict[str, Any],
    proposed: Dict[str, Any],
    policy: str = "fixed",
) -> bool:
    """Assert the proposed schedule honours the count policy.

    ``policy="fixed"`` (MVP default): ``n_states`` MUST be preserved — a
    count-neutral rebalance. Raises ``CountPolicyError`` if the count changed.

    ``policy="additive"`` (Stage-2 gated fallback): the count may only GROW
    (never shrink); a decrease still raises. Count growth requires explicit
    explicit user approval — that approval is the CALLER's
    responsibility; this guard only enforces the monotonicity.
    """
    if policy not in ("fixed", "additive"):
        raise ValueError(
            f"policy must be 'fixed' or 'additive', got {policy!r}"
        )
    n_seed = len(seed["directions"])
    n_prop = len(proposed["directions"])
    if policy == "fixed":
        if n_prop != n_seed:
            raise CountPolicyError(
                f"count_policy='fixed' but n_states changed: seed {n_seed} → "
                f"proposed {n_prop}"
            )
    else:  # additive
        if n_prop < n_seed:
            raise CountPolicyError(
                f"count_policy='additive' but n_states SHRANK: seed {n_seed} "
                f"→ proposed {n_prop} (additive may only grow)"
            )
    return True


def assert_all_constraints(
    seed: Dict[str, Any],
    proposed: Dict[str, Any],
    config: Optional[AdaptiveConfig] = None,
) -> bool:
    """Run all four load-bearing guards on a (seed, proposed) pair.

    Raises the first violated guard's typed error. Returns True iff ALL pass.
    This is the single entry point prescribe.py calls before returning a
    proposal so a malformed rebalance can never escape the module.
    """
    cfg = config or AdaptiveConfig()
    count_policy(seed, proposed, policy=cfg.count_policy)
    endpoint_invariant(seed, proposed, config=cfg)
    region_locked(seed, proposed, config=cfg)
    symmetry_mirror(proposed, config=cfg)
    return True
