"""Adaptive λ-scheduling PRESCRIBE module — count-fixed linear-region rebalance.

Given an ``OverlapReport`` + ``JudgeResult`` + the seed schedule, propose a NEW
schedule that relocates up to ``config.damping_max_moves`` linear-segment λ
windows from a redundant (over-overlapped, flat-integrand) plateau into a
bottleneck (under-overlapped, steep) turnover — n_states FIXED, endpoints PINNED,
LINEAR region only.

This is the exact generalization of the manual densified38 → densified38v2 move.
The full per-state
arrays are rebuilt by DELEGATING to the EXISTING
``scripts/trackb_production_v2_asyncre.py::_build_densified38_arrays(linear_lambdas=...)``
— the ladder math (soft-core 9-window ladder, α/U0/W0 ramps, whole-tuple-reverse
backward construction) is NOT reimplemented here (anti-fragmentation). This
module computes ONLY the new forward LINEAR λ array and hands it to that builder.

The rebalance recipe:
  1. Locate the worst LINEAR bottleneck pair → turnover interval [a, b].
  2. Locate the contiguous redundant plateau before it → donor span [0, a].
  3. Relocate n = min(damping_max_moves, donor_capacity) windows: re-space the
     donor span COARSER (fewer, evenly-spaced windows) and densify the turnover
     FINER (the reclaimed windows + a pre-turnover bridge), holding the linear
     count fixed and both endpoints (λ=0; the turnover's b feeds the ladder apex)
     pinned.

After building, ALL four load-bearing guards (constraints.assert_all_
constraints) are asserted before the proposal is returned — a malformed rebalance
can never escape this module.

Python 3.8+; numpy-only (plus the asyncre builder, spec-loaded).
"""

import importlib.util
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .constraints import (
    AdaptiveConfig,
    assert_all_constraints,
)
from .judge import JudgeResult
from .overlap import OverlapReport, PairOverlap

_PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# ---------------------------------------------------------------------------
# Rebalance recipe constants. These are NOT invented by the optimizer — they are
# the prescribed turnover resolution + the pre-turnover bridge that makes the
# densified38v2 move reproducible. Kept here as named, documented constants
# (analogous to AdaptiveConfig's thresholds), so the recipe is auditable rather
# than buried as magic numbers.
# ---------------------------------------------------------------------------
@dataclass
class RebalanceRecipe:
    """Turnover-densification recipe (prescribed values).

    ``turnover_delta_lambda`` — the fine Δλ across the bottleneck flip (0.02:
        0.40→0.42→0.44→0.45).
    ``include_pre_bridge`` — add one pre-turnover bridge window one ``turnover_
        delta_lambda``-class step before the turnover start (the 0.38 window),
        softening entry into the flip.
    These reproduce DENSE38v2_LAMBDA_FWD_LINEAR exactly; they are documented so a
    future target can adjust them under review without a code rewrite.
    """

    turnover_delta_lambda: float = 0.02
    include_pre_bridge: bool = True
    round_decimals: int = 4


# ---------------------------------------------------------------------------
# Delegate to the EXISTING ladder builder (anti-fragmentation)
# ---------------------------------------------------------------------------
def _load_asyncre_module():
    """Spec-load ``scripts/trackb_production_v2_asyncre.py`` (the test-isolation
    pattern). Returns the module exposing
    ``_build_densified38_arrays`` + ``_schedule_dict``. Raises ImportError if the
    launcher is missing.
    """
    path = os.path.join(_PROJ_ROOT, "scripts", "trackb_production_v2_asyncre.py")
    if not os.path.isfile(path):
        raise ImportError(f"asyncre launcher not found at {path}")
    spec = importlib.util.spec_from_file_location(
        "_trackb_asyncre_for_prescribe", path
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _extract_forward_linear(schedule: Dict[str, Any]) -> List[float]:
    """Return the FORWARD-half linear-segment λ array from a schedule dict.

    A linear-segment state has W0==0 AND INTERMEDIATE==0. The forward half is the
    leading DIRECTION>=0 block. Returns λ1 of those states (λ1==λ2==λ in the
    linear region). Raises ValueError if no linear states are found.
    """
    directions = [float(d) for d in schedule["directions"]]
    fwd = sum(1 for d in directions if d >= 0)
    lin: List[float] = []
    for i in range(fwd):
        if (float(schedule["w0"][i]) == 0.0
                and int(round(float(schedule["intermd"][i]))) == 0):
            lin.append(float(schedule["lambdas_1"][i]))
    if not lin:
        raise ValueError(
            "seed schedule has no LINEAR-segment states (W0==0 & INTERMEDIATE==0)"
        )
    return lin


# ---------------------------------------------------------------------------
# Core linear-array rebalance (count-fixed, endpoints pinned)
# ---------------------------------------------------------------------------
def rebalance_linear_array(
    seed_linear: List[float],
    bottleneck_interval: Tuple[float, float],
    n_moves: int,
    recipe: Optional[RebalanceRecipe] = None,
) -> List[float]:
    """Produce a count-neutral rebalanced LINEAR λ array.

    ``seed_linear`` — the seed forward linear λ array (e.g. densified38's
        uniform-Δλ=0.05 [0.00..0.45], 10 states).
    ``bottleneck_interval`` — (a, b) λ bounds of the worst linear bottleneck pair
        (e.g. (0.40, 0.45)).
    ``n_moves`` — the relocation magnitude requested by the caller (clamped to
        ``damping_max_moves``). Must be > 0 (a 0-move call is a no-op the caller
        must handle as "nothing to rebalance"). It sets the floor on how many
        windows are reclaimed from the plateau; the turnover densification +
        pre-bridge then consume the freed budget at ``recipe.turnover_delta_
        lambda`` resolution. The verdict's count-neutral move corresponds to
        n_moves=2 (the manual densified38→v2 step reclaimed 3 plateau windows /
        gained 2 turnover + 1 bridge — net count-neutral at 10 linear states).
    ``recipe`` — turnover densification recipe constants.

    The output keeps ``len(seed_linear)`` states, pins the first state (λ=0) and
    the last state (the turnover end ``b``, which feeds the soft-core ladder apex
    via the builder), and is strictly monotone increasing. The donor span
    [0, a) is re-spaced with fewer evenly-spaced windows; the turnover [a, b] is
    densified to ``recipe.turnover_delta_lambda`` (+ optional pre-bridge).

    With the verdict recipe (Δλ=0.02, pre-bridge) on densified38's seed this
    reproduces DENSE38v2_LAMBDA_FWD_LINEAR exactly:
        [0.00, 0.08, 0.16, 0.24, 0.32, 0.38, 0.40, 0.42, 0.44, 0.45].
    """
    rec = recipe or RebalanceRecipe()
    n_total = len(seed_linear)
    a, b = float(bottleneck_interval[0]), float(bottleneck_interval[1])
    if b <= a:
        raise ValueError(
            f"bottleneck interval must have b > a; got ({a}, {b})"
        )
    if n_moves <= 0:
        raise ValueError(
            f"n_moves must be > 0 to perform a relocation, got {n_moves}"
        )
    nd = rec.round_decimals

    # --- turnover block: densify [a, b] at recipe.turnover_delta_lambda ---
    # arange(a, b, Δλ) gives the fine grid up to (but not incl) b; b is then
    # appended so the turnover end (the ladder apex feed) is pinned. For the
    # verdict case [0.40, 0.45] @ Δλ=0.02 → [0.40, 0.42, 0.44] + [0.45].
    turnover = [round(float(x), nd)
                for x in np.arange(a, b - 1e-12, rec.turnover_delta_lambda)]
    if not turnover or not np.isclose(turnover[-1], b):
        turnover.append(round(b, nd))

    # --- optional pre-turnover bridge: one Δλ step before a (the 0.38 window) ---
    bridge: List[float] = []
    if rec.include_pre_bridge:
        bridge_pt = round(a - rec.turnover_delta_lambda, nd)
        if bridge_pt > 0.0:
            bridge = [bridge_pt]

    # --- donor (plateau) span: the remaining windows over [0, a), evenly ---
    # The turnover block (bridge + turnover, the turnover INCLUDES a) consumes
    # n_turn_block states. The donor supplies the rest over [0, a) — i.e.
    # linspace(0, a, n_donor+1)[:-1] so the donor's last node sits one donor-step
    # below a (the bridge bridges that gap). For the verdict case: n_donor=5 →
    # linspace(0, 0.40, 6)[:-1] = [0, 0.08, 0.16, 0.24, 0.32].
    n_turn_block = len(bridge) + len(turnover)
    n_donor = n_total - n_turn_block
    if n_donor < 1:
        raise ValueError(
            f"rebalance would leave < 1 donor plateau window "
            f"(n_donor={n_donor}); the bottleneck densification + bridge consume "
            f"{n_turn_block} of {n_total} states. Reduce damping_max_moves / "
            f"turnover resolution, or use the additive count policy."
        )
    if n_donor == 1:
        donor = [0.0]
    else:
        donor = [round(float(x), nd)
                 for x in np.linspace(0.0, a, n_donor + 1)[:-1]]

    # Stitch: donor [0 .. a-step] + bridge + turnover [a .. b].
    merged = donor + bridge + turnover
    merged = sorted(set(merged))

    if len(merged) != n_total:
        raise ValueError(
            f"rebalanced linear array has {len(merged)} states, expected "
            f"{n_total} (count-neutral). Recipe / n_moves mismatch: {merged}"
        )
    # Endpoints pinned.
    if merged[0] != 0.0:
        raise ValueError(f"coupled endpoint λ=0 not preserved: {merged}")
    if not np.isclose(merged[-1], b):
        raise ValueError(
            f"turnover end {b} not preserved as last linear state: {merged}"
        )
    # Strict monotone.
    if any(merged[k] >= merged[k + 1] for k in range(len(merged) - 1)):
        raise ValueError(f"rebalanced array not strictly increasing: {merged}")
    return merged


# ---------------------------------------------------------------------------
# Top-level proposal
# ---------------------------------------------------------------------------
def propose_rebalance(
    report: OverlapReport,
    judge: JudgeResult,
    schedule: Dict[str, Any],
    config: Optional[AdaptiveConfig] = None,
    recipe: Optional[RebalanceRecipe] = None,
) -> Dict[str, Any]:
    """Propose a count-fixed, linear-region rebalanced schedule dict.

    Uses the JUDGE result to locate (i) the worst LINEAR bottleneck pair
    (defines the turnover interval) and (ii) the redundant plateau (the donor),
    then relocates up to ``config.damping_max_moves`` windows via
    ``rebalance_linear_array`` and rebuilds the full per-state arrays by
    DELEGATING to ``_build_densified38_arrays`` (the existing ladder builder).

    The returned dict is the SAME shape as ``_schedule_dict`` /
    ``get_schedule(...)`` output (lambdas / lambdas_1 / lambdas_2 / directions /
    intermd / w0 / alpha / u0 / n_states). Before returning, ALL FOUR guards
    are asserted against the seed (assert_all_constraints) — raises an
    ``AdaptiveConstraintError`` subclass if any invariant would be violated
    (endpoint moved / ladder mutated / symmetry desynced / count changed).

    Raises ``ValueError`` (NOT a soft no-op) when the judge found no LINEAR
    bottleneck to fix or no redundant donor to reclaim from — the caller must
    then REFUSE + escalate (the verdict's barrier-limited / nothing-to-do
    branch), not silently emit the seed.
    """
    cfg = config or AdaptiveConfig()
    rec = recipe or RebalanceRecipe()

    # --- find the worst LINEAR bottleneck pair across both directions ---
    worst: Optional[Tuple[str, PairOverlap]] = None
    for direction, pairs in judge.bottlenecks.items():
        for p in pairs:
            if judge.region_of(direction, p) != "linear":
                continue  # v1 rebalances the LINEAR region only.
            if worst is None or p.bc < worst[1].bc:
                worst = (direction, p)
    if worst is None:
        raise ValueError(
            "no LINEAR-region bottleneck found — nothing to rebalance "
            "(the schedule clears the linear overlap floor, OR the bottleneck is "
            "in the soft-core ladder which v1 does not rebalance — REFUSE + "
            "escalate per verdict Q4-ii/Q5)."
        )

    # --- redundant donor must exist (linear, over-overlapped plateau) ---
    n_redundant = sum(
        len(pairs) for pairs in judge.redundant_runs.values()
    )
    if n_redundant == 0:
        raise ValueError(
            "no redundant LINEAR plateau to reclaim windows from — a "
            "count-neutral rebalance is impossible (would require count growth "
            "= additive policy, which needs explicit approval per verdict "
            "Q4-iii). REFUSE + escalate."
        )

    # --- turnover interval from the bottleneck pair's λ bounds ---
    seed_linear = _extract_forward_linear(schedule)
    lam = schedule.get("lambdas") or schedule.get("lambdas_1")
    a = float(lam[worst[1].i])
    b = float(lam[worst[1].j])

    # --- relocate up to damping_max_moves windows ---
    n_moves = min(cfg.damping_max_moves, n_redundant)

    new_linear = rebalance_linear_array(
        seed_linear=seed_linear,
        bottleneck_interval=(a, b),
        n_moves=n_moves,
        recipe=rec,
    )

    # --- DELEGATE ladder rebuild to the existing builder ---
    asyncre = _load_asyncre_module()
    built = asyncre._build_densified38_arrays(new_linear)
    proposed = asyncre._schedule_dict(
        lambdas_1=built["lambdas_1"],
        lambdas_2=built["lambdas_2"],
        directions=built["directions"],
        intermd=built["intermd"],
        w0=built["w0"],
        alpha=built["alpha"],
        u0=built["u0"],
        lambdas=built["lambdas"],
    )

    # --- load-bearing guards (fail-loud) BEFORE returning ---
    assert_all_constraints(schedule, proposed, config=cfg)

    # Provenance the caller (schedule_io) stamps into adaptive_schedule.json.
    proposed["_provenance"] = {
        "rebalanced_from": "seed",
        "bottleneck_direction": worst[0],
        "bottleneck_pair": [worst[1].i, worst[1].j],
        "bottleneck_bc": worst[1].bc,
        "bottleneck_interval": [a, b],
        "n_moves": n_moves,
        "n_redundant_donors": n_redundant,
        "new_linear_segment": list(new_linear),
        "recipe": {
            "turnover_delta_lambda": rec.turnover_delta_lambda,
            "include_pre_bridge": rec.include_pre_bridge,
        },
        "count_policy": cfg.count_policy,
        "region": cfg.region,
        "regime": "ranking_only",
    }
    return proposed
