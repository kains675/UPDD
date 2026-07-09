# -*- coding: utf-8 -*-
"""Track B in-place residue-4 RBFE — System serializer + λ-ladder + asyncre bridge.

This module is the BRIDGE from the validated in-place fused build
(``atm_trackB_setup.build_inplace_res4_fused_system``, genuine HE1<->methyl swap,
S0-harmonized RBFE XML, both legs Tier-1 + Tier-2 PASS) to a runnable
Hamiltonian replica-exchange λ-LADDER. A single λ=0.5 frame proves the box is
finite + the swap is wired; an actual ΔΔG_bind needs the ladder run + UWHAM.

ARCHITECTURE / INTEGRATION FINDING (R-18, honest):
    The validated per-direction asyncre DRIVER
    (``scripts/trackb_per_direction_production.py``) is hard-wired to the
    upstream ABFE engine: it spawns ``atom_openmm.abfe_production`` whose
    ``OMMSystemABFE.create_system`` (a) deserializes a PLAIN System XML (NO
    ATMForce) and (b) builds its OWN ATMForce at runtime from the cntl
    ``LIGAND_ATOMS`` + ``DISPLACEMENT`` (whole-binder decoupling). The upstream
    RBFE class ``OMMSystemRBFE`` likewise rebuilds its ATMForce from a plain
    System + a TWO-attach-atom var-region protocol
    (``add_common_var_atoms_to_atmforce``, ommsystem.py:773-776:
    ``ParticleOffsetDisplacement(lig2_attach, lig1_attach)`` for lig1's var
    atoms and the mirror for lig2's). That protocol assumes the upstream
    TWO-COPY overlay where ``lig1_attach != lig2_attach``.

    Our in-place fused System is a SINGLE-SHARED-CORE box: both var groups
    ({CM,HM1-3} and {HE1}) attach the SAME physical NE1, so feeding it through
    the upstream var-region protocol yields offset 0 for BOTH groups => u1 == u0
    (the documented single-shared-core NULL op). The genuine HE1 decoupling
    (massless dummy NE1 reference + internal-angle zeroing) that makes u1 != u0
    lives ONLY in ``atm_trackB_setup.attach_inplace_swap_atmforce`` — it is NOT
    reproducible from a plain System XML + upstream cntl keywords.

    => The per-direction ABFE driver CANNOT consume the in-place System as-is
       without its displacement/LIGAND_ATOMS assumptions (a real blocker, NOT
       forced). The faithful bridge is to serialize the System WITH the genuine
       ATMForce already attached, and run the λ-ladder through a THIN in-process
       asyncre adapter that (i) NEVER rebuilds the ATMForce, only sets its
       per-state global parameters, and (ii) REUSES the driver's mixing gate
       VERBATIM by emitting a driver log it already knows how to parse
       (``Replica N new state M``). The ABFE displacement path is untouched.

This module is PURELY ADDITIVE (off-campaign C8): it imports the validated
build verbatim (C1), introduces no change to the ABFE schedules / driver path /
S0 charges / shared XML / densified arrays.

Ranking-only (R-11); this is a v0.8 PREDICTION-test bridge (R-18), NOT a
converged ΔΔG_bind. The initial window count is a PARAMETER, validated
empirically by the short asyncre smoke (mixing / round-trips), NOT pre-asserted.

DOI references:
  - Gallicchio 2021 J Chem Theory Comput, DOI 10.1021/acs.jctc.1c00753
  - Azimi et al. 2022 J Chem Inf Model 62(2):309, DOI 10.1021/acs.jcim.1c01129
  - Mey et al. 2020 LiveCoMS, DOI 10.33011/livecoms.2.1.18378
"""

from __future__ import annotations

import gc
import math
import os
import random as _random
import sys
from typing import Any, Dict, List, Optional, Tuple

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import atm_trackB_setup as ats  # noqa: E402


# ---------------------------------------------------------------------------
# RBFE λ-ladder (Task 2).
#
# The ATM cntl STRUCTURE is identical to the ABFE ladder (per-state LAMBDA1 /
# LAMBDA2 / ALPHA / U0 / W0COEFF / DIRECTION / INTERMEDIATE), but the in-place
# HE1<->methyl edit is a ~4-5-atom perturbation with NO whole-binder
# displacement cliff to bridge — so the soft-core perturbation magnitude is
# ones-to-tens of kcal/mol, not the ~191 kcal/mol ABFE decoupling collapse.
# Therefore FAR fewer windows are needed than the 30-40-state ABFE ladders.
#
# Following the upstream temoa-g1 RBFE example pattern (a single symmetric λ
# 0->0.5->0 ramp with the λ2 partner offset in the soft-core anneal band and a
# W0 ramp), we build a forward half (λ1 0->0.5) and a backward half (the
# whole-tuple reverse with DIRECTION=-1) so the ladder is symmetric. The
# INTERMEDIATE (soft-core) band sits at the λ=0.5 apex where λ1 < λ2 so the
# ilogistic softplus is active; the linear segment uses λ1 == λ2, W0 = 0.
#
# The window COUNT is a PARAMETER (n_windows_half). Do NOT hard-assume the
# converged count — it is validated empirically (Task 4 short smoke measures
# adjacent replica-exchange mixing; add windows where overlap is thin).
# ---------------------------------------------------------------------------

# Soft-core canon (NOT re-tuned; mirrors atm_trackB_setup ATS_* + upstream).
RBFE_ALPHA_DEFAULT = 0.10
RBFE_U0_DEFAULT = 30.0      # Uh knee for a ~tens-kcal/mol perturbation (NOT the
                            # 110 ABFE value tuned for the ~191 kcal/mol cliff).
RBFE_W0_APEX = 1.0          # W0 ramp apex at the λ=0.5 soft-core midpoint.
RBFE_UMAX_KCAL = ats.ATS_UMAX_KCAL      # 200.0  (global, fixed)
RBFE_UBCORE_KCAL = ats.ATS_UBCORE_KCAL  # 100.0  (global, fixed)
RBFE_ACORE = ats.ATS_ACORE              # 0.0625 (global, fixed)
RBFE_TEMP_K = 300.0

# Staged-minimization (OPT-IN, W4A large-box stability). The W4A 9-heavy fused-
# indole two-copy box (~292k atoms) under-minimizes at the standard small budget:
# a fresh PME water shell carries close contacts a few-hundred-iter minimize cannot
# relieve, so the first integration step detonates (cycle-0 NaN). The validated fix
# (W4A/w4a_bound_smoke.py) is a STAGED relax: minimize FIRST at the reference state
# (soft-core OFF) with the project production-floor iteration budget, then a short
# polish at the assigned state, then a brief MD warmup. These constants are the
# staged-path FLOOR + warmup length; they are used ONLY when staged_min=True (the
# default-off path never references them), so the non-staged minimize budget is
# untouched. P1 (2026-04-22 cycle-0 NaN incident): 5000-iter minimum.
STAGED_MIN_ITERS_FLOOR = 5000
STAGED_WARMUP_STEPS = 200

# Canonical ATS two-copy U0 (Uh) knee. The two-copy box's swap is a REAL
# coordinate transfer (copy-2 displaced into bulk -> swapped to the site), so its
# soft-core perturbation magnitude is in the ABFE-cliff band, not the single-
# shared-core null-op tens-of-kcal. The standard ATS/ABFE Uh = 110 kcal/mol is
# the matching knee (NOT the RBFE_U0_DEFAULT=30 tuned for the single-shared-core
# null op). umax/ubcore/acore are the same global canon (NOT re-tuned).
ATS_TWOCOPY_U0_DEFAULT = 110.0


def build_rbfe_ladder(
    n_windows_half: int = 8,
    softcore_band: int = 2,
    alpha: float = RBFE_ALPHA_DEFAULT,
    u0_kcal: float = RBFE_U0_DEFAULT,
    w0_apex_kcal: float = RBFE_W0_APEX,
    n_apex_bridge: int = 0,
    apex_band: float = 0.5,
    single_direction: Optional[str] = None,
) -> Dict[str, Any]:
    """Build the symmetric in-place RBFE λ-ladder (forward + backward halves).

    ``n_windows_half`` is the number of states in the forward (DIRECTION=+1)
    half, λ1 sweeping 0 -> 0.5 inclusive. The backward half is the whole-tuple
    reverse with DIRECTION=-1. Total states = 2 * n_windows_half (+ the apex
    bridge, see ``n_apex_bridge``).

    ``single_direction`` (DEFAULT None -> the full symmetric two-direction
    ladder, byte-identical to the pre-existing behaviour) selects a STANDALONE
    single-direction ladder — the DECISIVE fork test for whether the per-
    direction-separate production estimator is viable on the single-shared-core
    box:

      * ``"forward"`` -> ONLY the forward (DIRECTION=+1) half: ``n_windows_half``
        (+ ``n_apex_bridge``) states, λ1 0 -> 0.5, state 0 = the λ=0 endpoint,
        the LAST state = the λ=0.5 apex. NO backward states, NO direction-flip
        apex BOUNDARY. This is exactly the dplus standalone REXEE ladder the ABFE
        per-direction driver runs (a separate asyncre run per direction, merged
        post-hoc by UWHAM at the shared λ=0.5 apex λ-state — the two directions
        never exchange with each other).
      * ``"backward"`` -> ONLY the backward (DIRECTION=-1) half (the whole-tuple
        reverse, state 0 = apex side, last state = λ=0 endpoint). The dminus
        standalone ladder, for symmetry.

    For a single-direction ladder the mixing gate's two ends (state 0 / state
    K-1) are the TWO PHYSICAL λ-endpoints of THAT direction (λ=0 and λ=0.5
    apex), so ``round_trips`` measures the INTERNAL 0<->apex traversal of one
    standalone direction — NOT the combined-ladder apex Dir-flip handoff (which
    the just-completed finding proved is a STRUCTURAL wall in the single-shared-
    core box, see the IMPORTANT note below). If a single-direction standalone
    ladder round-trips (round_trips > 0) the per-direction-separate estimator is
    viable without the two-copy overlay; if it does NOT, the single-shared-core
    box genuinely cannot sample a converged direction and the two-copy overlay
    is required (R-18 — the smoke reports the numbers, this builder does not
    pre-assert the verdict).

    ``softcore_band`` is the number of states at the λ=0.5 APEX of each half
    that carry the INTERMEDIATE soft-core anneal (λ1 < λ2, W0 ramp, INTERMEDIATE
    flag = 1). The remaining (linear) states use λ1 == λ2, W0 = 0,
    INTERMEDIATE = 0 — the same convention as the ABFE forward leg.

    ``n_apex_bridge`` (DEFAULT 0 — byte-identical to the pre-bridge ladder) adds
    that many EXTRA soft-core windows on EACH half just below the λ=0.5 apex, in
    a CONCAVE W0 grid that compresses spacing as W0 -> 1.0 (the same "finer W0
    near saturation" lever the densified ABFE ladder used at its backward
    W0-peak; see ``scripts/trackb_production_v2_asyncre.py`` DENSE38 escalation).
    ``apex_band`` is the fraction of the soft-core band (measured from λ=0.5
    inward) the bridge windows occupy (default 0.5 = the half nearest the apex).

    IMPORTANT (R-18, empirical finding): for the in-place SINGLE-SHARED-CORE box
    the apex handoff (forward state Dir=+1 base=u0 <-> backward state Dir=-1
    base=u1) wall is INVARIANT to W0/λ/α because the ATM hybrid potential bases
    on the UN-softened ``select(step(Direction), u0, u1)`` term (ommsystem.py:817)
    and the soft-core caps only the PERTURBATION ``usc``, not the base. With the
    in-place box's ~56,000-kcal/mol ``u1-u0`` clash (both partners attach the
    SAME NE1, so the decoupled partner overlaps the coupled one), this apex-base
    barrier (~233,000 kJ/mol exchange Δ) is NOT a soft-core-bridgeable W0-peak.
    The bridge below is correct infrastructure (and would close a genuine
    moderate-perturbation W0-peak wall) but does NOT close the single-shared-core
    apex — that needs the two-copy overlay where u0/u1 each physically displace
    their partner. Window/bridge counts are PARAMETERS validated empirically by
    the asyncre mixing smoke (R-18 — NOT pre-assumed converged).

    Returns the schedule dict in the SAME shape ``write_cntl_file`` /
    ``schedule_io`` consume: ``lambdas`` / ``lambdas_1`` / ``lambdas_2`` /
    ``directions`` / ``intermd`` / ``alpha`` / ``u0`` / ``w0`` plus the GLOBAL
    soft-core canon (``umax`` / ``ubcore`` / ``acore``) and metadata.
    """
    if n_windows_half < 2:
        raise ValueError(
            "build_rbfe_ladder: n_windows_half must be >= 2 (need at least the "
            "λ=0 and λ=0.5 end states), got %d" % (n_windows_half,))
    if softcore_band < 1 or softcore_band > n_windows_half:
        raise ValueError(
            "build_rbfe_ladder: softcore_band must be in [1, n_windows_half], "
            "got %d (n_windows_half=%d)" % (softcore_band, n_windows_half))
    if n_apex_bridge < 0:
        raise ValueError(
            "build_rbfe_ladder: n_apex_bridge must be >= 0, got %d"
            % (n_apex_bridge,))
    if not (0.0 < apex_band <= 1.0):
        raise ValueError(
            "build_rbfe_ladder: apex_band must be in (0, 1], got %r"
            % (apex_band,))
    if single_direction not in (None, "forward", "backward"):
        raise ValueError(
            "build_rbfe_ladder: single_direction must be one of "
            "{None, 'forward', 'backward'}, got %r" % (single_direction,))

    # Forward half: λ1 linearly 0 -> 0.5 over n_windows_half states.
    lam1_fwd: List[float] = []
    lam2_fwd: List[float] = []
    w0_fwd: List[float] = []
    inter_fwd: List[int] = []
    alpha_fwd: List[float] = []
    u0_fwd: List[float] = []

    n = n_windows_half
    # The last ``softcore_band`` states (closest to λ=0.5) are the anneal band.
    band_start = n - softcore_band
    for i in range(n):
        lam1 = round(0.5 * i / (n - 1), 6)
        if i >= band_start:
            # Soft-core anneal: λ2 leads λ1 by one ladder step so the ilogistic
            # softplus is active; W0 ramps 0 -> w0_apex across the band; the
            # state is flagged INTERMEDIATE.
            step = 0.5 / (n - 1)
            lam2 = round(min(0.5, lam1 + step), 6)
            band_pos = i - band_start + 1          # 1..softcore_band
            w0 = round(w0_apex_kcal * band_pos / softcore_band, 6)
            inter = 1
        else:
            lam2 = lam1
            w0 = 0.0
            inter = 0
        lam1_fwd.append(lam1)
        lam2_fwd.append(lam2)
        w0_fwd.append(w0)
        inter_fwd.append(inter)
        alpha_fwd.append(alpha)
        u0_fwd.append(u0_kcal)

    # Apex bridge: insert n_apex_bridge EXTRA soft-core windows on the forward
    # half just below the λ=0.5 apex (the W0 -> 1.0 saturation edge). They are
    # placed strictly BETWEEN the penultimate band state and the apex state, with
    # a CONCAVE W0 grid (compressed as W0 -> 1.0), λ1==λ2 fixed at the apex λ
    # (0.5) so the bridge lives entirely in the W0 anneal (no new λ insertion).
    # This MUTATES only the apex-saturation approach; default n_apex_bridge=0
    # leaves the arrays byte-identical to the pre-bridge ladder.
    if n_apex_bridge > 0:
        apex_lam = lam1_fwd[-1]                       # 0.5
        w0_pre = w0_fwd[-2] if n >= 2 else 0.0        # W0 of the penultimate band
        w0_apex = w0_fwd[-1]                          # the apex W0 (== w0_apex_kcal)
        a_pre = alpha_fwd[-1]
        u_pre = u0_fwd[-1]
        # Concave interpolation: fraction f^p with p>1 compresses near f=1
        # (the apex). The bridge spans the apex_band fraction of [w0_pre, w0_apex]
        # nearest the apex, so the EXISTING penultimate->apex gap is subdivided
        # with finer steps approaching saturation.
        lo = w0_apex - apex_band * (w0_apex - w0_pre)
        bridge = []
        for j in range(1, n_apex_bridge + 1):
            f = j / float(n_apex_bridge + 1)          # 0 < f < 1
            fc = f ** 2                                # concave (denser near apex)
            w0_b = round(lo + (w0_apex - lo) * fc, 6)
            bridge.append((apex_lam, apex_lam, w0_b, a_pre, u_pre))
        # Splice the bridge windows in just before the apex state.
        insert_at = len(lam1_fwd) - 1
        for (l1b, l2b, w0b, ab, ub) in bridge:
            lam1_fwd.insert(insert_at, l1b)
            lam2_fwd.insert(insert_at, l2b)
            w0_fwd.insert(insert_at, w0b)
            inter_fwd.insert(insert_at, 1)
            alpha_fwd.insert(insert_at, ab)
            u0_fwd.insert(insert_at, ub)
            insert_at += 1

    # STANDALONE single-direction ladder (the DECISIVE fork test). Build ONLY the
    # forward (or reversed-forward = backward) half as a self-contained REXEE
    # ladder with ONE DIRECTION value — NO backward states, NO direction-flip apex
    # BOUNDARY. The mixing gate's two ends (state 0 / state K-1) are the two
    # PHYSICAL λ-endpoints of THIS direction (λ=0 endpoint and λ=0.5 apex), so its
    # round_trips measure the INTERNAL 0<->apex traversal of one standalone
    # direction — exactly the per-direction-separate production geometry (each of
    # dplus/dminus runs as its own asyncre ladder, merged post-hoc by UWHAM at the
    # shared apex; the two directions never exchange). This is the fully reusable
    # forward-half construction above (C1) wired into a single-direction ladder; it
    # does NOT touch the symmetric two-direction return path below.
    if single_direction is not None:
        if single_direction == "forward":
            sd_lam1 = list(lam1_fwd)
            sd_lam2 = list(lam2_fwd)
            sd_w0 = list(w0_fwd)
            sd_inter = list(inter_fwd)
            sd_alpha = list(alpha_fwd)
            sd_u0 = list(u0_fwd)
            sd_dir = 1
        else:  # "backward" — the whole-tuple reverse, DIRECTION = -1.
            sd_lam1 = list(reversed(lam1_fwd))
            sd_lam2 = list(reversed(lam2_fwd))
            sd_w0 = list(reversed(w0_fwd))
            sd_inter = list(reversed(inter_fwd))
            sd_alpha = list(reversed(alpha_fwd))
            sd_u0 = list(reversed(u0_fwd))
            sd_dir = -1
        sd_directions = [sd_dir] * len(sd_lam1)
        sd_total = len(sd_lam1)
        return {
            "lambdas": list(sd_lam1),
            "lambdas_1": sd_lam1,
            "lambdas_2": sd_lam2,
            "directions": sd_directions,
            "intermd": sd_inter,
            "alpha": sd_alpha,
            "u0": sd_u0,
            "w0": sd_w0,
            "umax": RBFE_UMAX_KCAL,
            "ubcore": RBFE_UBCORE_KCAL,
            "acore": RBFE_ACORE,
            "n_states": sd_total,
            "n_windows_half": sd_total,
            "n_windows_half_linear": n,
            "softcore_band": softcore_band,
            "n_apex_bridge": n_apex_bridge,
            "apex_band": apex_band,
            "single_direction": single_direction,
            "temperature_K": RBFE_TEMP_K,
            "schedule_name": "inplace_rbfe_%s_%dw" % (single_direction, sd_total),
            "regime": "ranking_only",
            "note": ("In-place residue-4 RBFE STANDALONE %s λ-ladder "
                     "(HE1<->methyl, ~4-5 atom perturbation, DIRECTION=%d only, "
                     "NO apex direction-flip boundary). The DECISIVE fork test: "
                     "does ONE direction round-trip INTERNALLY (state 0 <-> apex)? "
                     "round_trips > 0 => the per-direction-separate estimator is "
                     "viable on the single-shared-core box (apex wall was a "
                     "combined-ladder artifact, NO two-copy overlay needed). "
                     "Window count is a PARAMETER validated empirically by the "
                     "asyncre mixing smoke (R-18)." % (single_direction, sd_dir)),
        }

    # Backward half = whole-tuple reverse of the forward half, DIRECTION = -1
    # (λ1 stays λ1, λ2 stays λ2 — NOT a λ1<->λ2 swap; symmetric ladder).
    lam1_bwd = list(reversed(lam1_fwd))
    lam2_bwd = list(reversed(lam2_fwd))
    w0_bwd = list(reversed(w0_fwd))
    inter_bwd = list(reversed(inter_fwd))
    alpha_bwd = list(reversed(alpha_fwd))
    u0_bwd = list(reversed(u0_fwd))

    lambdas_1 = lam1_fwd + lam1_bwd
    lambdas_2 = lam2_fwd + lam2_bwd
    w0 = w0_fwd + w0_bwd
    intermd = inter_fwd + inter_bwd
    alpha_arr = alpha_fwd + alpha_bwd
    u0_arr = u0_fwd + u0_bwd
    # The forward half now has (n + n_apex_bridge) states; the backward half
    # mirrors it. DIRECTION splits at the per-half boundary (NOT at n).
    half = len(lam1_fwd)
    directions = [1] * half + [-1] * half
    # LAMBDAS column == λ1 for the symmetric ladder (matches the ABFE default).
    lambdas = list(lambdas_1)

    total = 2 * half
    return {
        "lambdas": lambdas,
        "lambdas_1": lambdas_1,
        "lambdas_2": lambdas_2,
        "directions": directions,
        "intermd": intermd,
        "alpha": alpha_arr,
        "u0": u0_arr,
        "w0": w0,
        "umax": RBFE_UMAX_KCAL,
        "ubcore": RBFE_UBCORE_KCAL,
        "acore": RBFE_ACORE,
        "n_states": total,
        "n_windows_half": half,
        "n_windows_half_linear": n,
        "softcore_band": softcore_band,
        "n_apex_bridge": n_apex_bridge,
        "apex_band": apex_band,
        "single_direction": None,
        "temperature_K": RBFE_TEMP_K,
        "schedule_name": "inplace_rbfe_%dw" % (total,),
        "regime": "ranking_only",
        "note": ("In-place residue-4 RBFE λ-ladder (HE1<->methyl, ~4-5 atom "
                 "perturbation, NO displacement cliff). Window/apex-bridge counts "
                 "are PARAMETERS validated empirically by the asyncre mixing smoke "
                 "(R-18) — NOT a pre-assumed converged ladder. NOTE: the apex "
                 "direction-flip wall in the single-shared-core box is NOT "
                 "soft-core-bridgeable (the ATM base term u0/u1 is un-softened); "
                 "it requires the two-copy overlay."),
    }


def build_ats_standard_ladder(
    n_windows_half: int = 6,
    single_direction: str = "forward",
    alpha: float = RBFE_ALPHA_DEFAULT,
    u0_kcal: float = ATS_TWOCOPY_U0_DEFAULT,
    w0_apex_kcal: float = 0.0,
    *,
    lambda1_rampdown: Optional[List[float]] = None,
    lambda2_rampup: Optional[List[float]] = None,
) -> Dict[str, Any]:
    """Build the CANONICAL single-direction ATS λ-schedule (spec C).

    This is the standard AToM (ATM/ATS) per-leg ladder used with the CANONICAL
    TWO-COPY box, distinct from :func:`build_rbfe_ladder` (which targets the
    single-shared-core box with a λ1==λ2 linear ramp + a small apex anneal band).
    For 11 states (``n_windows_half=6``) it reproduces the the specified
    schedule exactly:

      λ1 = [0, 0, 0, 0, 0, 0,   0.1, 0.2, 0.3, 0.4, 0.5]
      λ2 = [0, 0.1, 0.2, 0.3, 0.4, 0.5,   0.5, 0.5, 0.5, 0.5, 0.5]

    i.e. the standard two-phase ATM ramp: the FIRST ``n_windows_half`` states hold
    λ1=0 while λ2 climbs 0 -> 0.5 (the "leg up" — the soft-core alchemical region
    where λ2 != λ1 so the ilogistic softplus is active), then the LAST
    ``n_windows_half - 1`` states hold λ2=0.5 while λ1 climbs 0.1 -> 0.5 to meet it
    (the "leg down"). The two overlapping λ1=0/λ2=0.5 ... λ1=0.5/λ2=0.5 states are
    the symmetric λ=0.5 apex of the ATM hybrid. Total states = 2*n_windows_half-1.

    Every state of a single-leg ATS ladder carries the SAME DIRECTION (+1 for the
    forward/"leg up→down" ladder, -1 for the reverse) — there is NO direction-flip
    boundary (that lives only in the COMBINED two-direction ladder, which is the
    single-shared-core artifact this canonical schedule replaces). The two
    directions are run as SEPARATE asyncre runs (per-direction-separate estimator)
    and UWHAM-merged at the shared apex, exactly the ABFE geometry.

    ``single_direction`` ("forward" default, or "backward" = the whole-tuple
    reverse with DIRECTION=-1) selects which standalone direction this ladder is.

    ``lambda1_rampdown`` (keyword-only, default ``None``) lets the caller densify
    the leg-down (λ2=0.5, λ1 climbing) phase at the leg-switch boundary WITHOUT
    touching the soft-core canon (C7) or the leg-up λ2 ramp. ``None`` keeps the
    legacy uniform λ1 ramp (``round(0.5*i/(n-1))`` over the last n-1 states) — so
    ``n_windows_half=6, lambda1_rampdown=None`` is byte-identical to the historical
    11-state schedule. When provided it is the EXPLICIT list of leg-down λ1 knots
    (the apex λ1=0/λ2=0.5 state is already placed by the leg-up phase, so the list
    starts ABOVE 0): every value must be in (0, 0.5], strictly increasing, and end
    at exactly 0.5. The leg-down then has ``len(lambda1_rampdown)`` states and the
    total becomes ``n_windows_half + len(lambda1_rampdown)``. This is an
    interior-λ reshaping of the leg-switch handoff — ΔG-unbiased (Kirkwood path
    independence): the endpoints (decoupled / coupled apex) and the soft-core
    constants are unchanged, so it is a robustness (overlap) lever, not a free
    energy change (ranking-only, R-11). Example: ``[0.05,0.1,0.2,0.3,0.4,0.5]``
    inserts a λ1=0.05 bridge window right after the leg-switch apex (12 states for
    ``n_windows_half=6``).

    ``lambda2_rampup`` (keyword-only, default ``None``) is the analogous lever on
    the OTHER axis: it densifies the leg-up (λ1=0, λ2 climbing) soft-core phase at
    the deep-λ2 decouple tail WITHOUT touching the soft-core canon (C7) or the
    leg-down λ1 ramp. ``None`` keeps the legacy uniform λ2 ramp
    (``round(0.5*i/(n-1))`` over the n leg-up states) — so
    ``n_windows_half=6, lambda2_rampup=None`` is byte-identical to the historical
    11-state schedule. When provided it is the EXPLICIT list of leg-up λ2 knots:
    unlike λ1 (whose decoupled λ1=0 endpoint is placed by the leg-up phase), the
    leg-up phase OWNS the genuine decoupled endpoint (λ2=λ1=0), so this list MUST
    start at exactly 0.0 and MUST end at exactly 0.5 (the apex λ2=0.5 where the
    leg-down begins). Every value must be in [0, 0.5] and strictly increasing. The
    leg-up then has ``len(lambda2_rampup)`` states. This is the same interior-λ
    reshaping (ΔG-unbiased; ranking-only, R-11) applied to the leg-up tail — the
    overlap lever for the deep-λ2 decouple bonds (λ2 0.2↔0.1↔0.0) that the
    leg-down ``lambda1_rampdown`` cannot reach. Example:
    ``[0.0,0.05,0.1,0.15,0.2,0.3,0.4,0.5]`` inserts a λ2=0.05 bridge between
    0/0.1 and a λ2=0.15 bridge between 0.1/0.2 (8 leg-up states for
    ``n_windows_half=6``). ``lambda1_rampdown`` (leg-down) and ``lambda2_rampup``
    (leg-up) are INDEPENDENT, composable axes — when both are given, the leg-up
    has ``len(lambda2_rampup)`` states and the leg-down ``len(lambda1_rampdown)``
    states, sharing the single λ1=0/λ2=0.5 apex (placed once, by the leg-up), so
    total = ``len(lambda2_rampup) + len(lambda1_rampdown)``.

    ``u0_kcal`` (Uh) defaults to the ATS canon 110 kcal/mol (the two-copy swap is
    a real ABFE-band transfer, NOT the tens-of-kcal single-shared-core null op).
    ``alpha`` = 0.10, ``w0_apex_kcal`` = 0.0 (the symmetric two-copy ΔG-unbiased
    apex — W0=0 means the two endpoints are NOT biased toward each other; the
    free energy comes out of the unbiased overlap). umax/ubcore/acore are the
    global canon, NOT re-tuned (C7).

    Returns the schedule dict in the SAME shape ``write_cntl_file`` /
    ``schedule_io`` / :class:`InplaceRbfeLadder` consume.
    """
    if n_windows_half < 2:
        raise ValueError(
            "build_ats_standard_ladder: n_windows_half must be >= 2 (need at "
            "least the λ=0 and λ=0.5 apex states), got %d" % (n_windows_half,))
    if single_direction not in ("forward", "backward"):
        raise ValueError(
            "build_ats_standard_ladder: single_direction must be 'forward' or "
            "'backward', got %r" % (single_direction,))

    n = n_windows_half
    lam1: List[float] = []
    lam2: List[float] = []
    inter: List[int] = []
    # Phase 1 ("leg up"): λ1 = 0, λ2 climbs 0 -> 0.5.
    if lambda2_rampup is None:
        # Legacy uniform ramp: λ2 = 0 -> 0.5 over n states. This branch is
        # byte-identical to the historical 11-state schedule for n=6.
        for i in range(n):
            l2 = round(0.5 * i / (n - 1), 6)
            lam1.append(0.0)
            lam2.append(l2)
            # λ2 != λ1 (except the i=0 endpoint) -> the soft-core alchemical region.
            inter.append(0 if i == 0 else 1)
    else:
        # Explicit leg-up λ2 knots (densify the deep-λ2 decouple tail). The leg-up
        # phase OWNS the genuine decoupled endpoint (λ2=λ1=0), so the knots MUST
        # start at 0.0; the last must be 0.5 (the apex λ1=0/λ2=0.5 where the
        # leg-down begins).
        knots2 = [float(x) for x in lambda2_rampup]
        if len(knots2) < 2:
            raise ValueError(
                "build_ats_standard_ladder: lambda2_rampup must have at least "
                "two knots (the λ2=0.0 decoupled endpoint and the λ2=0.5 apex), "
                "got %r" % (knots2,))
        for x in knots2:
            if not (0.0 <= x <= 0.5):
                raise ValueError(
                    "build_ats_standard_ladder: lambda2_rampup values must be "
                    "in [0, 0.5] (the leg-up climbs λ2 from the decoupled "
                    "endpoint 0.0 to the apex 0.5), got %r" % (x,))
        if knots2[0] != 0.0:
            raise ValueError(
                "build_ats_standard_ladder: lambda2_rampup must start at exactly "
                "0.0 (the genuine decoupled λ2=λ1=0 endpoint is placed by the "
                "leg-up phase), got %r" % (knots2[0],))
        if knots2[-1] != 0.5:
            raise ValueError(
                "build_ats_standard_ladder: lambda2_rampup must end at exactly "
                "0.5 (the apex λ1=0/λ2=0.5 where the leg-down begins), got %r"
                % (knots2[-1],))
        for a, b in zip(knots2, knots2[1:]):
            if not (b > a):
                raise ValueError(
                    "build_ats_standard_ladder: lambda2_rampup must be strictly "
                    "increasing, got %r" % (knots2,))
        for i, x in enumerate(knots2):
            l2 = round(x, 6)
            lam1.append(0.0)
            lam2.append(l2)
            # λ2 != λ1 (except the i=0 endpoint) -> the soft-core alchemical region.
            inter.append(0 if i == 0 else 1)
    # Phase 2 ("leg down"): λ2 = 0.5, λ1 climbs to 0.5.
    if lambda1_rampdown is None:
        # Legacy uniform ramp: λ1 = 0.1 -> 0.5 over the next n-1 states (skip i=0,
        # which is the apex λ1=0/λ2=0.5 state already placed above). This branch is
        # byte-identical to the historical 11-state schedule for n=6.
        for i in range(1, n):
            l1 = round(0.5 * i / (n - 1), 6)
            lam1.append(l1)
            lam2.append(0.5)
            # λ1 != λ2 until the final λ1=λ2=0.5 state (the symmetric apex endpoint).
            inter.append(0 if i == n - 1 else 1)
    else:
        # Explicit leg-down λ1 knots (densify the leg-switch handoff). The apex
        # λ1=0/λ2=0.5 state is already placed by the leg-up phase, so the knots
        # start ABOVE 0; the last must be 0.5 (the symmetric λ1=λ2=0.5 endpoint).
        knots = [float(x) for x in lambda1_rampdown]
        if len(knots) < 1:
            raise ValueError(
                "build_ats_standard_ladder: lambda1_rampdown must have at least "
                "one knot (the λ1=0.5 endpoint), got an empty list")
        for x in knots:
            if not (0.0 < x <= 0.5):
                raise ValueError(
                    "build_ats_standard_ladder: lambda1_rampdown values must be "
                    "in (0, 0.5] (the apex λ1=0 state is placed by the leg-up "
                    "phase; the final value must be 0.5), got %r" % (x,))
        for a, b in zip(knots, knots[1:]):
            if not (b > a):
                raise ValueError(
                    "build_ats_standard_ladder: lambda1_rampdown must be strictly "
                    "increasing, got %r" % (knots,))
        if knots[-1] != 0.5:
            raise ValueError(
                "build_ats_standard_ladder: lambda1_rampdown must end at exactly "
                "0.5 (the symmetric λ1=λ2=0.5 apex endpoint), got %r"
                % (knots[-1],))
        for x in knots:
            l1 = round(x, 6)
            lam1.append(l1)
            lam2.append(0.5)
            # λ1 != λ2 until the final λ1=λ2=0.5 state (the symmetric apex endpoint).
            inter.append(0 if l1 == 0.5 else 1)

    # Leg-up window count = n for the legacy uniform λ2 ramp; with an explicit
    # lambda2_rampup it is len(lambda2_rampup). The leg-down count is independent
    # (n-1 legacy or len(lambda1_rampdown)).
    n_legup = n if lambda2_rampup is None else len(lambda2_rampup)
    # total == 2*n - 1 for the all-legacy schedule; with explicit knot lists it is
    # n_legup (leg-up) + n_legdown (leg-down), sharing the single apex placed by
    # the leg-up phase (no double-count).
    total = len(lam1)
    w0 = [w0_apex_kcal] * total
    alpha_arr = [alpha] * total
    u0_arr = [u0_kcal] * total

    if single_direction == "forward":
        sd_dir = 1
        sd_lam1, sd_lam2, sd_w0 = lam1, lam2, w0
        sd_inter, sd_alpha, sd_u0 = inter, alpha_arr, u0_arr
    else:
        sd_dir = -1
        sd_lam1 = list(reversed(lam1))
        sd_lam2 = list(reversed(lam2))
        sd_w0 = list(reversed(w0))
        sd_inter = list(reversed(inter))
        sd_alpha = list(reversed(alpha_arr))
        sd_u0 = list(reversed(u0_arr))
    directions = [sd_dir] * total

    return {
        "lambdas": list(sd_lam1),
        "lambdas_1": sd_lam1,
        "lambdas_2": sd_lam2,
        "directions": directions,
        "intermd": sd_inter,
        "alpha": sd_alpha,
        "u0": sd_u0,
        "w0": sd_w0,
        "umax": RBFE_UMAX_KCAL,
        "ubcore": RBFE_UBCORE_KCAL,
        "acore": RBFE_ACORE,
        "n_states": total,
        "n_windows_half": total,
        "n_windows_half_linear": n_legup,
        "softcore_band": n_legup,     # the whole leg-up phase is soft-core
        "n_apex_bridge": 0,
        "apex_band": 1.0,
        "single_direction": single_direction,
        "schedule_kind": "ats_standard",
        "lambda1_rampdown": (list(lambda1_rampdown)
                             if lambda1_rampdown is not None else None),
        "lambda2_rampup": (list(lambda2_rampup)
                           if lambda2_rampup is not None else None),
        "temperature_K": RBFE_TEMP_K,
        "schedule_name": (
            "ats_standard_%s_%dw" % (single_direction, total)
            if (lambda1_rampdown is None and lambda2_rampup is None)
            else "ats_standard_%s_%dw_bridge%d"
                 % (single_direction, total, total)),
        "regime": "ranking_only",
        "note": ("CANONICAL ATS single-leg %s λ-schedule (spec C): λ1=0 / "
                 "λ2 climbs 0->0.5 (leg up) then λ2=0.5 / λ1 climbs 0->0.5 (leg "
                 "down), single DIRECTION=%d, NO direction-flip boundary. Used "
                 "with the CANONICAL TWO-COPY box (real coordinate transfer); the "
                 "two directions run SEPARATELY + UWHAM-merge at the shared apex. "
                 "Soft-core canon umax/ubcore/acore NOT re-tuned (C7); Uh=%g "
                 "kcal/mol (ATS canon). Window count is a PARAMETER validated by "
                 "the pilot (R-18)." % (single_direction, sd_dir, u0_kcal)),
    }


# ---------------------------------------------------------------------------
# System-XML serializer (Task 1).
#
# Serializes the in-place fused System WITH the genuine ATMForce already
# attached + the matching topology PDB the asyncre adapter loads. The upstream
# OMMSystem{ABFE,RBFE} classes would REBUILD the ATMForce from cntl keywords
# (incompatible with the single-shared-core in-place box — see module docstring),
# so the serialized System is consumed by the in-process adapter below, NOT by
# OMMSystemABFE.create_system.
# ---------------------------------------------------------------------------
def _write_system_and_pdb(system, modeller, sys_xml, pdb_path):
    """Serialize ``system`` to ``sys_xml`` + write the matching topology PDB.

    Shared by the single-core and two-copy serialize paths so both produce the
    SAME on-disk contract the loader + ladder consume. Handles the
    topology/position lockstep: if the System has MORE particles than the topology
    (the single-core genuine swap adds a MASSLESS dummy NE1 reference particle to
    the System + positions but not the topology), a matching topology atom is
    materialised per extra particle in a dedicated "X" chain so PDBFile.writeFile
    (which requires topology and positions to match) round-trips. This is
    serialization-only bookkeeping — it does NOT touch the validated build's
    System/energy (R-7). The two-copy box has NO such extra particle (its
    topology == System particle count), so n_extra == 0 and no "X" chain is
    added. Returns the final topology atom count.
    """
    import openmm as mm
    import openmm.unit as unit
    from openmm.app import PDBFile, element as _app_element

    with open(sys_xml, "w") as fh:
        fh.write(mm.XmlSerializer.serialize(system))

    topology = modeller.topology
    positions = modeller.positions   # a unit-wrapped position list (Quantity)
    n_extra = system.getNumParticles() - topology.getNumAtoms()
    if n_extra < 0:
        raise RuntimeError(
            "_write_system_and_pdb: topology (%d atoms) exceeds System particles "
            "(%d) — unexpected build state."
            % (topology.getNumAtoms(), system.getNumParticles()))
    if n_extra > 0:
        ref_chain = topology.addChain(id="X")
        ref_res = topology.addResidue("REFX", ref_chain, id="900")
        for _ in range(n_extra):
            topology.addAtom("DUM", _app_element.hydrogen, ref_res)

    # PDBFile.writeFile wants a unit-wrapped position sequence; strip + re-wrap in
    # nm so the (possibly Python-list) positions become a clean Quantity vector.
    pos_nm = [p.value_in_unit(unit.nanometer) for p in positions]
    pos_q = pos_nm * unit.nanometer
    with open(pdb_path, "w") as fh:
        PDBFile.writeFile(topology, pos_q, fh, keepIds=True)
    return topology.getNumAtoms()


def _serialize_twocopy_system(
    *,
    leg: str,
    out_dir: str,
    tag: str,
    seed: str,
    binder_chain: str,
    solvate: bool,
    harmonize_common_charges: bool,
    displacement_nm: float,
    mtr_ncaa_xml: Optional[str],
    constraints: Any,
    mutation_spec: Optional[Any] = None,
    auto_search_displacement: bool = False,
    accept_sep_nm: float = ats.ATS_TWOCOPY_ACCEPT_SEP_NM,
    leg_inputs: Optional[Dict[str, Any]] = None,
    appearing_h_retry_k: Optional[int] = None,
    carve_void_waters: bool = False,
    carve_cutoff_nm: float = ats.ATS_CARVE_VOID_CUTOFF_NM,
) -> Dict[str, Any]:
    """Build + serialize the CANONICAL ATS TWO-COPY box for one leg.

    The two-copy builder (``ats.build_inplace_res4_twocopy_system``) already
    attaches its own swap ATMForce (distinct attach atoms, real coordinate
    transfer), so this only serializes the merged System + topology onto the SAME
    on-disk contract the single-core path produces (so the loader + ladder are
    construction-agnostic). The build dict shape differs from single-core, so the
    return is mapped explicitly:

      - ``outcome`` must be ``"twocopy_attached"`` (the MC1 charge-discontinuity
        outcome RAISES — a real finding, resolve via the harmonized RBFE XML);
      - ``genuine_decouple_dir`` is the UNIT of the copy-2 displacement vector
        (the bulk copy IS the decoupled partner — its outward direction is the
        decouple direction the C8 bound-leg gate validates); for the free leg the
        same finite vector is returned (any bulk direction is valid free).
    """
    build_kwargs = dict(
        leg=leg, seed=seed, binder_chain=binder_chain, solvate=solvate,
        harmonize_common_charges=harmonize_common_charges,
        displacement_nm=displacement_nm, mtr_ncaa_xml=mtr_ncaa_xml,
        constraints=constraints, spec=mutation_spec,
        auto_search_displacement=auto_search_displacement,
        accept_sep_nm=accept_sep_nm, leg_inputs=leg_inputs,
        # Opt-in void-water carve (default OFF -> byte-identical): delete the whole
        # bulk waters that penetrate the swap-displaced disappearing-heavy volume
        # (backward-endpoint NaN-crash fix). Forwarded verbatim to the builder.
        carve_void_waters=carve_void_waters, carve_cutoff_nm=carve_cutoff_nm,
    )
    # P3-#116 FIX2 (opt-in, default None -> byte-identical legacy build): the bounded
    # R2-retry + deterministic per-unit appearing-H placement. When appearing_h_retry_k
    # is set the build routes through build_inplace_res4_twocopy_system_r2_retry, which
    # seeds the addHydrogens jitter deterministically per unit and re-places on the rare
    # R2 seed-clash FAIL (up to K attempts, K exhaustion is fail-loud). The retry is a
    # PURE placement fix — the forwarded build_kwargs (soft-core canon, λ-schedule,
    # templates/charges, box/PME, frozen FE core) are untouched. Default None keeps
    # every existing MTR/V3I/A9G/W4A caller on the exact legacy build.
    if appearing_h_retry_k is not None:
        # Per-unit identity for the deterministic seed: velocity-seed label + leg +
        # mutation name (direction is applied downstream in the ladder, so the build
        # is direction-agnostic). Stable across processes / build order.
        unit_key = "%s|%s|%s" % (
            seed, leg, getattr(mutation_spec, "name", None) or "default")
        build = ats.build_inplace_res4_twocopy_system_r2_retry(
            unit_key=unit_key, retry_k=appearing_h_retry_k, **build_kwargs)
    else:
        build = ats.build_inplace_res4_twocopy_system(**build_kwargs)
    if build.get("outcome") != "twocopy_attached":
        raise RuntimeError(
            "serialize_inplace_rbfe_system(construction='twocopy'): the two-copy "
            "build did NOT attach (outcome=%r). This is a real finding (MC1 "
            "charge discontinuity between copy-1 MTR and copy-2 WT common cores, "
            "or a failed seed), not a serialization step — resolve it before "
            "serializing. Supply an S0-harmonized RBFE XML (Σ|Δq|=0 common core) "
            "or pass harmonize_common_charges=True for the diagnostic path."
            % (build.get("outcome"),))

    fused = build["fused_build"]
    system = fused["system"]
    modeller = fused["modeller"]
    atm_index = build["swap"]["atm_force_index"]

    sys_xml = os.path.join(out_dir, "inplace_rbfe_%s_sys.xml" % (tag,))
    pdb_path = os.path.join(out_dir, "inplace_rbfe_%s.pdb" % (tag,))
    n_atoms = _write_system_and_pdb(system, modeller, sys_xml, pdb_path)

    # C8 decouple direction = the UNIT of the copy-2 bulk displacement vector. The
    # displaced bulk copy IS the decoupled partner, so the vector that carries it
    # into bulk (NE1-local-outward, finite + non-degenerate by construction) is the
    # decouple direction. Normalize (the build stores the d-SCALED vector).
    dvec = fused.get("displacement_vector_nm") or build.get("displacement_vector_nm")
    decouple_dir = None
    if dvec is not None:
        mag = math.sqrt(sum(float(c) ** 2 for c in dvec))
        if mag > 1e-9:
            decouple_dir = tuple(float(c) / mag for c in dvec)

    return {
        "leg": leg,
        "tag": tag,
        "seed": seed,
        "construction": "twocopy",
        "sys_xml_path": sys_xml,
        "pdb_path": pdb_path,
        "atmforce_index": atm_index,
        "n_atoms": n_atoms,
        "n_copy1": fused.get("n_copy1"),
        "solvated": solvate,
        "swap_mode": "twocopy",
        "displacement_vector_nm": [float(c) for c in dvec] if dvec else None,
        # task #100/#6: how d was chosen (fixed_direction vs auto_search) + the
        # per-build search trail (selected dir/magnitude/min-image sep) so the
        # run_manifest records it for the integrity audit + post-run decoupling check.
        "displacement_mode": build.get("displacement_mode"),
        "displacement_log": build.get("displacement_log"),
        # Opt-in void-water carve report (None when carve_void_waters=False) —
        # surfaced so the launcher can log per-leg (free vs bound) carve counts.
        "carve_report": fused.get("carve_report"),
        "common_charges_harmonized": build.get("common_charges_harmonized"),
        "mtr_ncaa_xml": build.get("mtr_ncaa_xml"),
        # C8 SIGN-critical: for the two-copy box the decouple direction is the
        # copy-2 bulk displacement unit vector (finite for BOTH legs — the bulk
        # copy is always displaced). The launcher's C8 gate asserts non-None +
        # finite + unit-magnitude for the bound leg.
        "genuine_decouple_dir": decouple_dir,
        "mc1_passed": True,
        "alchemical_atoms": {
            "common_attach_ne1": fused["alchemical_atoms"]["copy1_ne1"],
            "copy1_ne1": fused["alchemical_atoms"]["copy1_ne1"],
            "copy2_ne1": fused["alchemical_atoms"]["copy2_ne1"],
            "mtr_var": list(fused["alchemical_atoms"]["mtr_only"]),
            "wt_var": list(fused["alchemical_atoms"]["wt_only"]),
        },
        "separation": build.get("separation"),
        "swap": build.get("swap"),
        "regime": "ranking_only",
        # The build dict is kept so an in-process caller (the smoke) can run the
        # ladder WITHOUT re-deserializing (deserialize is verified separately).
        "_build": build,
    }


def serialize_inplace_rbfe_system(
    leg: str = "free",
    out_dir: str = ".",
    seed: str = "s7",
    binder_chain: str = "B",
    solvate: bool = True,
    harmonize_common_charges: bool = False,
    swap_mode: str = "genuine",
    genuine_decouple_nm: float = 1.2,
    mtr_ncaa_xml: Optional[str] = None,
    constraints: Any = None,
    tag: Optional[str] = None,
    construction: str = "single_core",
    displacement_nm: float = ats.ATS_TWOCOPY_DISPLACEMENT_NM,
    mutation_spec: Optional[Any] = None,
    auto_search_displacement: bool = False,
    accept_sep_nm: float = ats.ATS_TWOCOPY_ACCEPT_SEP_NM,
    leg_inputs: Optional[Dict[str, Any]] = None,
    appearing_h_retry_k: Optional[int] = None,
    carve_void_waters: bool = False,
    carve_cutoff_nm: float = ats.ATS_CARVE_VOID_CUTOFF_NM,
) -> Dict[str, Any]:
    """Build + serialize the in-place fused RBFE System for one leg.

    Writes ``<out_dir>/inplace_rbfe_<tag>_sys.xml`` (the OpenMM System,
    XmlSerializer, WITH the genuine ATMForce attached) and
    ``<out_dir>/inplace_rbfe_<tag>.pdb`` (the matching topology + coordinates +
    periodic box). ``tag`` defaults to ``<leg>``.

    ``construction`` selects which in-place box to serialize (DEFAULT
    ``"single_core"`` — byte-identical to the pre-existing legacy behaviour, the
    path densify_pilot / the existing ladder use):

      * ``"single_core"`` (default): the SINGLE-SHARED-CORE in-place fused box
        (``ats.build_inplace_res4_fused_system``) — both var groups attach the
        SAME NE1. Validated for the λ=0.5 finite frame + per-direction standalone
        ladders, but its apex Dir-flip is a structural wall and its swap is a
        soft-core-capped null op (pertE saturated ~150, overlap ~0). UNTOUCHED.
      * ``"twocopy"`` (opt-in): the CANONICAL ATS TWO-COPY box
        (``ats.build_inplace_res4_twocopy_system`` already attaches its own swap
        ATMForce) — copy-1 (MTR) at the site + copy-2 (WT) displaced by
        ``displacement_nm`` (~40 Å) into bulk, common coordinates swapped with
        DISTINCT attach atoms, clash avoided by spatial separation (NO inter-copy
        exclusions). The swap is a REAL coordinate transfer (one-frame |u1-u0|
        finite + non-saturated, ~few kcal/mol, validated). This is the path that
        can yield a converged ΔΔG (pilot needed to confirm; R-18).

    ``auto_search_displacement`` (two-copy ONLY, default False -> byte-identical
    fixed-direction legacy path): when True the copy-2 bulk displacement is chosen
    by the builder's direction-aware cone search (maximises the copy1<->copy2 +
    periodic-image min heavy-atom distance, escalates the magnitude only if no
    direction clears ``accept_sep_nm``). This recovers the bound-leg box build
    where a fixed-direction d drives copy-2's binder through copy-1's receptor
    body. d-/direction-NEUTRAL (the swap is partner-offset based; u1-u0 is
    d-invariant given full decoupling), so ranking-safe. ``accept_sep_nm`` is
    consulted only by the auto-search; the post-solvate C6 separation assert
    always enforces the 1.0 nm clash floor + the periodic-image gate.

    ``constraints=None`` (the DEFAULT here) matches the Tier-2 R3 requirement
    that the appearing/disappearing alch H carry NO SHAKE (a 1 fs unconstrained
    integrator is the validated stable choice; 2 fs without constraints is the
    documented instability signature). Pass ``constraints=HBonds`` (from
    ``openmm.app``) only for a 2 fs ladder.

    Returns a dict with the written paths + the build's ATMForce index +
    alchemical bookkeeping. Raises if the build surfaces the MC1 charge
    discontinuity (the System was NOT attached — a real finding, not a pass);
    pass ``harmonize_common_charges=True`` for the mechanical path or supply an
    S0-harmonized RBFE XML so MC1 passes on-disk.
    """
    if construction not in ("single_core", "twocopy"):
        raise ValueError(
            "serialize_inplace_rbfe_system: construction must be 'single_core' "
            "or 'twocopy', got %r" % (construction,))

    if tag is None:
        tag = leg
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # --- TWO-COPY (opt-in) ------------------------------------------------
    # The canonical ATS two-copy box already attaches its own swap ATMForce in
    # build_inplace_res4_twocopy_system; we only serialize the merged System +
    # topology. The build dict shape differs from single-core (outcome
    # "twocopy_attached", different alchemical_atoms keys, no genuine_decouple_dir
    # — the bulk decouple IS the copy-2 displacement vector), so this branch maps
    # it onto the SAME serialized contract the loader + ladder consume.
    if construction == "twocopy":
        return _serialize_twocopy_system(
            leg=leg, out_dir=out_dir, tag=tag, seed=seed,
            binder_chain=binder_chain, solvate=solvate,
            harmonize_common_charges=harmonize_common_charges,
            displacement_nm=displacement_nm, mtr_ncaa_xml=mtr_ncaa_xml,
            constraints=constraints, mutation_spec=mutation_spec,
            auto_search_displacement=auto_search_displacement,
            accept_sep_nm=accept_sep_nm, leg_inputs=leg_inputs,
            appearing_h_retry_k=appearing_h_retry_k,
            carve_void_waters=carve_void_waters, carve_cutoff_nm=carve_cutoff_nm)

    # --- SINGLE-CORE (legacy default; byte-identical) ---------------------
    # mutation_spec is two-copy-only (single-core is the MTR<->Trp single-shared-
    # core path); a non-None spec on single_core is a wiring error, not silently
    # ignored.
    if mutation_spec is not None:
        raise ValueError(
            "serialize_inplace_rbfe_system: mutation_spec is only supported with "
            "construction='twocopy' (the single_core path is the MTR<->Trp "
            "single-shared-core build).")
    # leg_inputs (the single-scaffold folding-thermocycle input override) is a
    # two-copy-only knob; a non-None override on single_core is a wiring error
    # (not silently ignored — the single_core build uses the 2QKI resolver).
    if leg_inputs is not None:
        raise ValueError(
            "serialize_inplace_rbfe_system: leg_inputs is only supported with "
            "construction='twocopy' (the single_core path resolves 2QKI inputs "
            "directly).")
    # auto_search_displacement is a two-copy-only displacement-construction knob
    # (single_core has no copy-2 bulk displacement to search); a True flag on
    # single_core is a wiring error, not silently ignored. accept_sep_nm is only
    # consulted by the auto-search, so its default is harmless on single_core.
    if auto_search_displacement:
        raise ValueError(
            "serialize_inplace_rbfe_system: auto_search_displacement is only "
            "supported with construction='twocopy' (the single_core path has no "
            "copy-2 bulk displacement to search).")
    # carve_void_waters is a two-copy-only build knob (the void arises from the
    # displaced bulk copy's disappearing-heavy atoms swapping into the partner
    # copy's site — single_core has no displaced bulk copy). A True flag on
    # single_core is a wiring error, not silently ignored.
    if carve_void_waters:
        raise ValueError(
            "serialize_inplace_rbfe_system: carve_void_waters is only supported "
            "with construction='twocopy' (the single_core path has no displaced "
            "bulk copy whose disappearing-heavy atoms swap into a water-filled "
            "void).")
    # appearing_h_retry_k (the P3-#116 FIX2 deterministic + bounded-retry appearing-H
    # placement) is a two-copy-only knob: the single_core MTR<->Trp path reads both
    # endpoints from their own final.pdb and never runs the canonical mutated-copy
    # addHydrogens placement, so a non-None value is a wiring error (not silently
    # ignored).
    if appearing_h_retry_k is not None:
        raise ValueError(
            "serialize_inplace_rbfe_system: appearing_h_retry_k is only supported "
            "with construction='twocopy' (the single_core path does not place the "
            "canonical mutated-copy appearing hydrogens).")
    build = ats.build_inplace_res4_fused_system(
        leg=leg, seed=seed, binder_chain=binder_chain, solvate=solvate,
        harmonize_common_charges=harmonize_common_charges, swap_mode=swap_mode,
        genuine_decouple_nm=genuine_decouple_nm, mtr_ncaa_xml=mtr_ncaa_xml,
        constraints=constraints,
    )
    if build.get("outcome") != "fused_attached":
        raise RuntimeError(
            "serialize_inplace_rbfe_system: the fused build did NOT attach "
            "(outcome=%r). This is a real finding (MC1 charge discontinuity or "
            "a failed seed), not a serialization step — resolve it before "
            "serializing. Pass harmonize_common_charges=True for the mechanical "
            "path or supply an S0-harmonized RBFE XML." % (build.get("outcome"),))

    fused = build["fused_build"]
    system = fused["system"]
    modeller = fused["modeller"]
    atm_index = build["swap"]["atm_force_index"]

    sys_xml = os.path.join(out_dir, "inplace_rbfe_%s_sys.xml" % (tag,))
    pdb_path = os.path.join(out_dir, "inplace_rbfe_%s.pdb" % (tag,))

    # Serialize the System + write the matching topology PDB (handles the
    # single-core genuine swap's MASSLESS dummy NE1 reference particle via the
    # shared lockstep helper — System particle count == topology atoms + 1, so
    # one "X"-chain placeholder atom is materialised; serialization-only, the
    # build's System/energy is untouched, R-7).
    n_atoms = _write_system_and_pdb(system, modeller, sys_xml, pdb_path)
    return {
        "leg": leg,
        "tag": tag,
        "seed": seed,
        "construction": "single_core",
        "sys_xml_path": sys_xml,
        "pdb_path": pdb_path,
        "atmforce_index": atm_index,
        "n_atoms": n_atoms,
        "solvated": solvate,
        "swap_mode": swap_mode,
        "common_charges_harmonized": build.get("common_charges_harmonized"),
        "mtr_ncaa_xml": build.get("mtr_ncaa_xml"),
        # C8 SIGN-critical: the genuine HE1 decouple direction. For the BOUND leg
        # this MUST be a finite outward (low-density) unit vector — the vector
        # ``compute_decouple_direction`` derives as NE1 - centroid(local heavy
        # atoms), which clears BOTH the receptor and the binder's own fold so the
        # decoupled HE1 lands in bulk (a naive away-from-receptor-centroid vector
        # is WRONG for this pose — it points back through the peptide ring). For
        # the FREE leg it is None by design (any direction is bulk for a free
        # peptide; the legacy fixed +Z is used). The launcher's C8 gate asserts
        # non-None + finite + unit-magnitude for the bound leg, fail-loud.
        "genuine_decouple_dir": build.get("genuine_decouple_dir"),
        "mc1_passed": True,
        "alchemical_atoms": {
            "common_attach_ne1": fused["alchemical_atoms"]["common"][0],
            "mtr_var": list(fused["alchemical_atoms"]["mtr_only"]),
            "wt_var_fused": list(fused["alchemical_atoms"]["wt_only"]),
        },
        "regime": "ranking_only",
        # The build dict is kept so an in-process caller (the smoke) can run the
        # ladder WITHOUT re-deserializing (deserialize is verified separately).
        "_build": build,
    }


def load_serialized_system(sys_xml_path: str, pdb_path: str) -> Dict[str, Any]:
    """Deserialize a serialized in-place RBFE System + topology (load-only check).

    This is the per-direction-driver-side load contract: deserialize the System
    XML, load the topology PDB, and confirm the ATMForce survived round-trip.
    Returns the system / topology / positions + the discovered ATMForce index.
    Raises if no ATMForce is present (the in-place box MUST carry its own
    ATMForce — unlike an ABFE plain System the driver would wrap at runtime).
    """
    import openmm as mm
    from openmm.app import PDBFile

    with open(sys_xml_path) as fh:
        system = mm.XmlSerializer.deserialize(fh.read())
    pdb = PDBFile(pdb_path)

    atm_index = None
    for i in range(system.getNumForces()):
        if isinstance(system.getForce(i), mm.ATMForce):
            atm_index = i
            break
    if atm_index is None:
        raise RuntimeError(
            "load_serialized_system: the deserialized System carries NO "
            "ATMForce. The in-place RBFE System must be serialized WITH its "
            "genuine ATMForce attached (the upstream driver would otherwise "
            "rebuild it from cntl keywords, which is incompatible with the "
            "single-shared-core box).")
    return {
        "system": system,
        "topology": pdb.topology,
        "positions": pdb.positions,
        "atmforce_index": atm_index,
        "n_atoms": pdb.topology.getNumAtoms(),
    }


class _LoadedBuildAdapter(object):
    """A minimal build-shaped view over a DESERIALIZED in-place box.

    ``ats.compute_decouple_direction`` consumes a ``build`` dict that exposes
    ``build["modeller"].topology`` / ``build["modeller"].positions`` and
    ``build["system"]``. A reused box (from :func:`load_serialized_system`) gives
    those three objects directly but is not a build dict, so this adapter wraps
    them with the SAME attribute / key contract — no System rebuild, no
    re-solvation. It is read-only and dict-like only for the keys
    ``compute_decouple_direction`` reads (``system``).
    """

    def __init__(self, loaded: Dict[str, Any]) -> None:
        self._system = loaded["system"]
        self.modeller = _ModellerView(loaded["topology"], loaded["positions"])

    def _resolve(self, key):
        # The keys ``compute_decouple_direction`` reads off a build dict:
        # ``build["modeller"]`` (subscript) and ``build.get("system")``.
        if key == "system":
            return self._system
        if key == "modeller":
            return self.modeller
        raise KeyError(key)

    def get(self, key, default=None):
        try:
            return self._resolve(key)
        except KeyError:
            return default

    def __getitem__(self, key):
        return self._resolve(key)


class _ModellerView(object):
    """A read-only ``modeller``-shaped view (topology + positions only)."""

    def __init__(self, topology, positions) -> None:
        self.topology = topology
        self.positions = positions


def recompute_decouple_direction_from_loaded(
    loaded: Dict[str, Any],
    binder_chain: str = "B",
    resnum: int = ats.ALCH_RESNUM,
    shell_nm: float = 1.0,
    decouple_nm: float = 1.2,
) -> Optional[Tuple[float, float, float]]:
    """Re-derive the C8 genuine decouple direction from a REUSED (deserialized) box.

    This is the box-reuse counterpart of the ``genuine_decouple_dir`` that
    ``serialize_inplace_rbfe_system`` returns when it BUILDS a fresh box. When a
    saved box is reused (``inplace_rbfe_<tag>_sys.xml`` + ``.pdb`` already on
    disk), the System is not rebuilt, so the decouple direction must be recovered.

    It is recomputed from the saved box's OWN geometry (``NE1 - centroid(local
    heavy atoms)`` via :func:`ats.compute_decouple_direction`) rather than trusted
    from a stale value — this RE-VALIDATES C8 against the actual box that will be
    sampled. The saved ``.pdb`` carries the addSolvent-output coordinates of the
    original build (the producing direction does not evolve the serialized box's
    stored coordinates), so the recomputed vector is consistent with that box. The
    massless dummy reference particle is excluded by the ``mass <= 1.5`` filter
    inside ``compute_decouple_direction`` (same as the fresh-build path), so the
    direction is computed over the real heavy-atom neighbourhood only.

    Returns the unit outward vector for the bound box (a finite 3-tuple), or
    ``None`` for the free box / when NE1 or its local shell cannot be resolved
    (the caller's C8 gate accepts ``None`` only for the free leg).
    """
    return ats.compute_decouple_direction(
        _LoadedBuildAdapter(loaded), binder_chain=binder_chain, resnum=resnum,
        shell_nm=shell_nm, decouple_nm=decouple_nm)


# ---------------------------------------------------------------------------
# Thin in-process asyncre adapter (Task 3).
#
# Runs the λ-ladder as a Hamiltonian replica-exchange on the in-place fused
# System WITHOUT rebuilding the ATMForce: each replica owns a Context; the
# per-state λ-tuple is applied via context.setParameter on the ATMForce
# globals; cross-state energies for the exchange acceptance come from
# ATMForce.getPerturbationEnergy(context) (u0/u1/pert per replica) recombined
# with the target state's (λ1,λ2,α,Uh,W0,Direction) — the SAME ATM hybrid
# potential the upstream engine evaluates, so this is the engine's math, not a
# hand-rolled surrogate.
#
# Crucially it EMITS a driver log line per replica per cycle in the
# ``Replica N new state M`` form the per-direction driver's mixing gate
# (parse_state_transitions_from_log / check_atm_mixing) already parses — so the
# window-count validator is REUSED VERBATIM, not reimplemented.
# ---------------------------------------------------------------------------
_KB_KCAL = 0.0019872041     # Boltzmann constant kcal/mol/K (kT at 300 K).


def _atm_softcore_components_kj(
    pert_kj: float, u0_kj: float, direction: float,
    umax_kj: float, ubcore_kj: float, acore: float,
) -> Tuple[float, float]:
    """Return (base_kj, usc_kj) — the ATM hybrid base + the SOFT-CORE-CAPPED
    perturbation, in kJ/mol.

    ``base`` = u0 (Direction>=0) or u1 (Direction<0) — the λ-INDEPENDENT
    reference. ``usc`` = sign(Direction)*softcore(u1-u0) — the soft-core-capped
    perturbation the ilogistic hybrid term operates on (NOT the raw u1-u0).
    This is the value the UWHAM estimator expects in its ``pertE`` column:
    upstream ``uwham._bias_fcn`` operates DIRECTLY on ``pertE`` with no soft-core
    cap of its own, so ``pertE`` MUST already be ``usc`` (a raw u1-u0 of tens of
    thousands of kcal/mol would overflow the WHAM solve).
    """
    u1_kj = u0_kj + pert_kj
    base = u0_kj if direction >= 0 else u1_kj
    sign = 1.0 if direction >= 0 else -1.0
    u = sign * (u1_kj - u0_kj)        # UOffset = 0

    # Soft-core cap usc(u) (ommsystem softCoreExpression).
    if acore <= 0.0 or u <= ubcore_kj:
        usc = u
    else:
        y = (u - ubcore_kj) / (umax_kj - ubcore_kj)
        z = 1.0 + 2.0 * (y / acore) + 2.0 * (y / acore) ** 2
        fsc = (z ** acore - 1.0) / (z ** acore + 1.0)
        usc = (umax_kj - ubcore_kj) * fsc + ubcore_kj
    return base, usc


def _atm_state_energy_kj(
    pert_kj: float, u0_kj: float,
    lambda1: float, lambda2: float, alpha_per_kj: float,
    uh_kj: float, w0_kj: float, direction: float,
    umax_kj: float, ubcore_kj: float, acore: float,
) -> float:
    """Recompute the ATM hybrid potential U_state for a sampled configuration.

    Reproduces the upstream ATMForce energy expression (ommsystem.py:817-823 /
    atm_trackB_setup._ATS_*) in Python from the per-replica raw perturbation
    ``pert = u1 - u0`` and the reference base ``u0`` (both from
    ``getPerturbationEnergy``), so the replica-exchange acceptance uses the
    engine's own potential — NOT a linear-λ surrogate. All energies in kJ/mol.

    U = base + hybrid(usc), where base = u0 (Direction>=0) or u1 (Direction<0),
    usc = softcore-capped u, u = sign(Direction)*(u1 - (u0 + UOffset)),
    UOffset = 0 here.
    """
    base, usc = _atm_softcore_components_kj(
        pert_kj, u0_kj, direction, umax_kj, ubcore_kj, acore)

    # Alchemical hybrid term (ilogistic softplus when λ2 != λ1).
    if lambda2 != lambda1:
        if alpha_per_kj > 0.0:
            hybrid = (((lambda2 - lambda1) / alpha_per_kj)
                      * math.log(1.0 + math.exp(-alpha_per_kj * (usc - uh_kj)))
                      + lambda2 * usc + w0_kj)
        else:
            hybrid = lambda2 * usc + w0_kj
    else:
        hybrid = lambda2 * usc + w0_kj
    return base + hybrid


# ---------------------------------------------------------------------------
# FIX-A: re-seeding / non-identity ladder initialisation (OPT-IN). Addresses the
# bound-leg ladder non-mixing pathology (decisive bond seal / fragmented ladder),
# validated across 3 wall seeds. FE-UNBIASED — an INITIAL-CONDITION change only:
# replica-exchange equilibrium is initial-condition-independent (Sugita-Okamoto
# 1999; Chodera-Shirts 2011 DOI 10.1063/1.3660669). DEFAULT OFF reproduces the
# current identity init byte-for-byte (V3I/A9G/MTR unaffected). C7 soft-core (U0/Uh=110,
# alpha=0.10, UMAX=200, UBCORE=100, ACORE=0.0625), the λ-schedule, every energy
# and the Metropolis criterion are FROZEN — this only changes WHERE each walker
# starts at t=0.
# ---------------------------------------------------------------------------
def make_reseed_permutation(n_states: int,
                            reseed_perm_seed: Optional[int] = None) -> List[int]:
    """Return the replica->state assignment used to initialise the ladder.

    ``perm[r]`` = the state replica ``r`` starts in at t=0.

    * ``reseed_perm_seed is None`` (DEFAULT) -> identity ``list(range(n_states))``,
      EXACTLY the legacy ``self.replica_state = list(range(self.n_states))`` init.
      This is the byte/behaviour-identical default-OFF path.
    * ``reseed_perm_seed`` is an int -> a RANDOMIZED bijection produced by a LOCAL
      ``random.Random(reseed_perm_seed)`` (NOT the global RNG, NOT date-based), so
      the permutation is fully reproducible from the logged integer and is recorded
      to the run_manifest. The result is always a valid bijection (a permutation of
      ``range(n_states)``), which is the only correctness requirement the exchange
      logic places on the init assignment (``run_cycle`` rebuilds
      ``state_to_replica`` from ``replica_state`` every cycle, so ANY bijection is
      a valid starting assignment and the stationary ensemble is unchanged).

    Reproducibility: ``make_reseed_permutation(n, k)`` returns the SAME list for
    the same ``(n, k)`` regardless of process / wall-clock / global RNG state.
    """
    n = int(n_states)
    if n <= 0:
        raise ValueError("make_reseed_permutation: n_states must be > 0")
    if reseed_perm_seed is None:
        return list(range(n))
    rng = _random.Random(int(reseed_perm_seed))
    perm = list(range(n))
    rng.shuffle(perm)
    return perm


def decoupled_endpoint_state(schedule: Dict[str, Any]) -> int:
    """Index of the genuine decoupled endpoint = the state with minimal λ2
    (ties -> minimal λ1) — the fully soft-core-decoupled basin (λ1=λ2=0) where
    the RCA's basin-locked walker parks.

    DIRECTION-CORRECT for BOTH legs (this is the false-green guard the smoke
    fixed): the decoupled endpoint is defined by its λ-VALUE (λ2=0), not by a
    fixed index. For the forward (dplus) leg that is state 0; for the backward
    (dminus) leg that is the LAST state. A max-λ1 rule (the prototype's forward-
    only band) would pick the COUPLED apex on a backward schedule and re-seed the
    wrong end. With ``lambda1_rampdown`` (densified leg-down) the endpoint index
    shifts but the λ2=0 rule still resolves it correctly.
    """
    best = None
    best_key = None
    for k in range(int(schedule["n_states"])):
        key = (schedule["lambdas_2"][k], schedule["lambdas_1"][k])
        if best is None or key < best_key:
            best = k
            best_key = key
    return best


def decoupled_band_states(schedule: Dict[str, Any],
                          band_lambda2_max: float = 0.25) -> List[int]:
    """States in the decoupled-endpoint band = the low-λ2 tail
    (``λ2 <= band_lambda2_max``) — the soft-core-saturated tail that holds the
    basin-locked walker. A2 re-seeds ONLY these contexts with the pre-relaxed
    decoupled config; everything with higher λ2 keeps the default coupled seed (no
    cross-contamination).

    DIRECTION-CORRECT by λ-VALUE: e.g. baseline 11-state dminus -> {8,9,10};
    densified 12-state dminus -> {9,10,11}; forward -> {0,1,2}. The band is the
    same physical set of states in either direction (just at different indices).
    """
    out = []
    for k in range(int(schedule["n_states"])):
        if schedule["lambdas_2"][k] <= float(band_lambda2_max):
            out.append(k)
    return out


def _subset_topology(topology, atom_indices: List[int]):
    """Return a new ``openmm.app.Topology`` holding only ``atom_indices``.

    DCDFile records the FULL position array unless given a subset topology with
    a matching subset of coordinates. The subset preserves chain/residue/atom
    grouping (so mdtraj/MDAnalysis can resolve res-4 χ / pucker / φ/ψ from the
    written frames) and copies the source periodic box vectors. Used only by the
    OPT-IN DCD path; the default (dcd_dir=None) run never calls it. Bonds are
    intentionally NOT copied (DCD stores coordinates only; the analysis tool
    re-derives connectivity from the residue templates).
    """
    import openmm.app as app

    sub = app.Topology()
    keep = set(int(a) for a in atom_indices)
    new_chain = {}
    new_res = {}
    for atom in topology.atoms():
        if atom.index not in keep:
            continue
        ch = atom.residue.chain
        if ch.index not in new_chain:
            new_chain[ch.index] = sub.addChain(ch.id)
        res = atom.residue
        if res.index not in new_res:
            new_res[res.index] = sub.addResidue(
                res.name, new_chain[ch.index], res.id)
        sub.addAtom(atom.name, atom.element, new_res[res.index])
    box = topology.getPeriodicBoxVectors()
    if box is not None:
        sub.setPeriodicBoxVectors(box)
    return sub


class InplaceRbfeLadder(object):
    """In-process Hamiltonian replica-exchange driver for the in-place RBFE box.

    One replica (= one Context) per ladder state. Each cycle: propagate every
    replica a few MD steps at its CURRENT state's λ-tuple, then attempt nearest-
    neighbor state exchanges (Metropolis on the ATM hybrid potential). Writes a
    ``Replica N new state M`` log line per replica per cycle so the per-direction
    driver's mixing gate parses it unchanged.

    This adapter REUSES the validated build's ATMForce verbatim (it only sets
    global parameters); it does NOT rebuild forces and does NOT touch the ABFE
    displacement path.
    """

    def __init__(self, system, positions, schedule, platform_name="CUDA",
                 temperature_K=RBFE_TEMP_K, timestep_fs=1.0,
                 friction_per_ps=1.0, log_path=None, seed=None,
                 minimize_iters=500, backward_equil_steps=500,
                 out_dir=None, out_basename="trackb",
                 staged_min=False, staged_min_iters=STAGED_MIN_ITERS_FLOOR,
                 staged_warmup_steps=STAGED_WARMUP_STEPS,
                 reseed_perm_seed=None, reseed_endpoint=False,
                 reseed_endpoint_band_lambda2_max=0.25,
                 reseed_endpoint_equil_steps=2000,
                 reseed_endpoint_minimize_iters=500,
                 dcd_dir=None, dcd_topology=None, dcd_atom_indices=None,
                 dcd_stride_cycles=1):
        import openmm as mm
        import openmm.app  # noqa: F401 — registers mm.app.DCDFile for the opt-in DCD path
        import openmm.unit as unit

        self.mm = mm
        self.unit = unit
        self.schedule = schedule
        # OPT-IN staged minimization (W4A 9-heavy fused-indole large-box stability;
        # validated in W4A/w4a_bound_smoke.py). DEFAULT OFF => the existing single-
        # stage minimize path runs unchanged (V3I/MTR/A9G/V3A byte-identical). When
        # ON, each replica is staged-relaxed: (a) >=5000-iter minimize at the
        # reference state (Lambda1=Lambda2=0, soft-core OFF) to relieve the fresh
        # PME water shell, (b) a short polish minimize at the assigned state, then
        # (c) setVelocitiesToTemperature + a short MD warmup. The standard
        # ``minimize_iters`` default is NOT changed (staged uses its OWN floor so
        # the non-staged path's iteration count is untouched). See the per-replica
        # loop below for the staged branch.
        self.staged_min = bool(staged_min)
        self.staged_min_iters = max(int(staged_min_iters),
                                    STAGED_MIN_ITERS_FLOOR) if staged_min else \
            int(staged_min_iters)
        self.staged_warmup_steps = int(staged_warmup_steps)
        # OPT-IN FIX-A re-seeding (W4A bound-leg ladder mixing). DEFAULT OFF =>
        # the legacy identity init + default coordinate seed run unchanged
        # (V3I/A9G/MTR byte-identical). reseed_perm_seed=None keeps the identity
        # replica->state map; reseed_endpoint=False keeps every context seeded
        # from the same coupled ``positions``. When ON, A1 sets a logged-seed
        # permutation as the t=0 assignment and A2 seeds the decoupled-band
        # contexts from a standalone equilibration at the genuine decoupled
        # endpoint's own λ-tuple. FE-unbiased (init-condition only); see the
        # per-state context loop below and ``make_reseed_permutation`` /
        # ``decoupled_endpoint_state`` / ``decoupled_band_states``.
        self.reseed_perm_seed = (None if reseed_perm_seed is None
                                 else int(reseed_perm_seed))
        self.reseed_endpoint = bool(reseed_endpoint)
        self.reseed_endpoint_band_lambda2_max = \
            float(reseed_endpoint_band_lambda2_max)
        self.reseed_endpoint_equil_steps = int(reseed_endpoint_equil_steps)
        self.reseed_endpoint_minimize_iters = \
            int(reseed_endpoint_minimize_iters)
        self.n_states = schedule["n_states"]
        self.temperature_K = temperature_K
        self.kT_kj = (unit.MOLAR_GAS_CONSTANT_R * temperature_K * unit.kelvin
                      ).value_in_unit(unit.kilojoule_per_mole)
        self.log_path = log_path
        self._log_fh = open(log_path, "w") if log_path else None
        self._cycle = 0
        import datetime as _dt
        self._log_base_ts = _dt.datetime(2000, 1, 1, 0, 0, 0)

        # UWHAM-consumable per-WALKER ``.out`` writers (Task 1, production). When
        # ``out_dir`` is set the adapter writes one ``<out_dir>/r{r}/<basename>.out``
        # per walker, appending the abfe_production-format row each cycle:
        #   stateid temperature direction lambda1 lambda2 alpha u0 w0 potE pertE 0
        # (col 0 = stateid, col 8 = potE, col 9 = pertE — exactly the columns the
        # UWHAM estimator ``_calculate_uwham_multi_intermediate`` reads). potE is
        # the FULL ATM hybrid potential at the walker's CURRENT state and pertE is
        # the raw perturbation u1-u0 (both in kcal/mol, matching the upstream
        # ommreplica.save_out convention). This is the same per-walker layout the
        # asyncre engine emits, so the existing merge + UWHAM path consumes it
        # unchanged — the SHORT smoke leaves out_dir=None so it is purely additive.
        self.out_dir = out_dir
        self.out_basename = out_basename
        self._out_fhs = None
        if out_dir is not None:
            self._out_fhs = []
            for r in range(self.n_states):
                rdir = os.path.join(out_dir, "r%d" % (r,))
                os.makedirs(rdir, exist_ok=True)
                self._out_fhs.append(
                    open(os.path.join(rdir, out_basename + ".out"), "w"))

        # OPT-IN per-walker DCD trajectory (probe diagnostic). DEFAULT
        # ``dcd_dir=None`` => NO trajectory is written and EVERY path below is
        # byte-identical to the legacy run (the production legs do not pass it).
        # When set, one ``<dcd_dir>/r{r}/<basename>.dcd`` is written per walker
        # Context: each cycle (every ``dcd_stride_cycles`` cycles) one frame of
        # the ``dcd_atom_indices`` subset is appended. The subset is the
        # alchemical-region + receptor-context atoms (NOT the ~290k PME waters),
        # so the frames resolve the res-4 sidechain χ1/χ2 swap, ring pucker, and
        # res-4 backbone φ/ψ for the under-sampling-vs-real-basin judgment. This
        # is OBSERVATION ONLY: it reads each Context's positions AFTER all
        # estimator logic (energy / exchange / .out row) is complete, so it
        # cannot perturb the dgbind1 / UWHAM result. A walker frame is labelled
        # by the state it occupies via the per-cycle ``.out`` stateid column +
        # the frame index (a frame and an .out row are written in lockstep).
        self.dcd_dir = dcd_dir
        self.dcd_stride_cycles = max(1, int(dcd_stride_cycles))
        self._dcd_files = None
        self._dcd_atom_indices = None
        if dcd_dir is not None:
            if dcd_topology is None:
                raise ValueError(
                    "InplaceRbfeLadder: dcd_dir set but dcd_topology is None "
                    "(the DCDFile header needs the topology atom set). Pass the "
                    "loaded box topology.")
            n_top = dcd_topology.getNumAtoms()
            if dcd_atom_indices is None:
                # Full box (every atom) — large; the launcher normally supplies
                # a binder+receptor subset to keep the probe DCD ~MB not ~GB.
                self._dcd_atom_indices = list(range(n_top))
            else:
                idx = [int(a) for a in dcd_atom_indices]
                if not idx:
                    raise ValueError(
                        "InplaceRbfeLadder: dcd_atom_indices is empty (no atoms "
                        "to record).")
                for a in idx:
                    if a < 0 or a >= n_top:
                        raise ValueError(
                            "InplaceRbfeLadder: dcd_atom_indices entry %d out of "
                            "range [0, %d)." % (a, n_top))
                self._dcd_atom_indices = idx
            self._dcd_sub_topology = _subset_topology(
                dcd_topology, self._dcd_atom_indices)
            self._dcd_files = []
            for r in range(self.n_states):
                rdir = os.path.join(dcd_dir, "r%d" % (r,))
                os.makedirs(rdir, exist_ok=True)
                fh = open(os.path.join(rdir, out_basename + ".dcd"), "wb")
                dcdf = mm.app.DCDFile(
                    fh, self._dcd_sub_topology,
                    dt=timestep_fs * unit.femtoseconds)
                self._dcd_files.append((fh, dcdf))

        # Resolve the ATMForce on the system.
        self.atm_index = None
        for i in range(system.getNumForces()):
            if isinstance(system.getForce(i), mm.ATMForce):
                self.atm_index = i
                break
        if self.atm_index is None:
            raise RuntimeError(
                "InplaceRbfeLadder: system carries no ATMForce.")
        self.atmforce = system.getForce(self.atm_index)

        # Soft-core globals (kJ/mol) — fixed across states.
        self.umax_kj = ats._kcal_to_kj(schedule["umax"])
        self.ubcore_kj = ats._kcal_to_kj(schedule["ubcore"])
        self.acore = schedule["acore"]

        # One Context per state. Each owns its own integrator (Langevin) so MD
        # propagation at the per-state λ runs independently.
        try:
            self.platform = mm.Platform.getPlatformByName(platform_name)
        except Exception:
            self.platform = mm.Platform.getPlatformByName("Reference")
            platform_name = "Reference"
        self.platform_name = platform_name

        # replica r currently occupies state self.replica_state[r]. FIX-A (A1,
        # OPT-IN): with ``reseed_perm_seed`` set this is a logged-seed permutation
        # of the states rather than the identity map (DEFAULT None => identity =
        # ``list(range(self.n_states))`` byte-identical). The permutation is set
        # BEFORE the per-state context loop below, so each context is minimized at
        # the PERMUTED state it will occupy (``state_k = self.replica_state[k]``),
        # not at an identity state — the in-__init__ integration is exact (no
        # post-construction relabel). The exchange logic rebuilds state_to_replica
        # from replica_state every cycle, so any bijection is a valid t=0
        # assignment; the stationary ensemble is unchanged (FE-unbiased).
        self.replica_state = make_reseed_permutation(
            self.n_states, self.reseed_perm_seed)

        # Identify the backward (Direction<0) ladder endpoint state (λ minimal on
        # the backward half) — the u1 = HE1-decoupled basin the backward replicas
        # must relax into BEFORE the soft-core anneal-edge.
        self._backward_endpoint = self._find_backward_endpoint()

        # FIX-A (A2, OPT-IN): pre-relax the genuine decoupled endpoint (min-λ2
        # corner, direction-correct for BOTH dplus and dminus) and capture its
        # relaxed configuration so the decoupled-BAND contexts can be seeded from
        # the decoupled basin instead of the coupled ``positions`` (which leaves
        # the high-usc walker basin-locked behind the soft-core wall). DEFAULT
        # reseed_endpoint=False => this is skipped and every context keeps the
        # default coupled seed (byte-identical). The relaxed config comes from a
        # STANDALONE equilibration AT the decoupled state's OWN λ-tuple (no
        # artificial/biased config) — it changes KINETICS (starting
        # basin), not the EQUILIBRIUM the estimator averages over (FE-unbiased).
        self._reseed_band = []
        self._reseed_relaxed_positions = None
        self._reseed_endpoint_state = None
        if self.reseed_endpoint:
            self._reseed_band = decoupled_band_states(
                self.schedule, self.reseed_endpoint_band_lambda2_max)
            if self._reseed_band:
                self._reseed_endpoint_state = decoupled_endpoint_state(
                    self.schedule)
                self._reseed_relaxed_positions = self._equilibrate_endpoint(
                    system, positions, self._reseed_endpoint_state,
                    seed=seed,
                    minimize_iters=self.reseed_endpoint_minimize_iters,
                    equil_steps=self.reseed_endpoint_equil_steps)

        # One Context per state. Each replica is initialised at its OWN state's
        # λ-tuple BEFORE minimization (the per-state ATMForce globals change the
        # potential, so minimization must run at the state the replica occupies).
        # A light minimization at the assigned state is REQUIRED on the fresh
        # PME-solvated box: the genuine swap's dummy-NE1 reference + the
        # re-solvated water placement leave local strain that detonates the FIRST
        # integration step otherwise (the validated Tier-2 path minimizes 500
        # iters before integrating — same requirement here).
        #
        # PER-DIRECTION SEEDING (the load-bearing fix for the backward soft-core
        # anneal-edge instability): the backward half's base potential is u1 (the
        # HE1-decoupled state), whose equilibrium geometry differs from the u0
        # (WT-coupled) starting positions. Minimizing a backward anneal-edge
        # state from the u0 geometry lands on a singular (PE=inf) configuration
        # that NaNs on the first step. So backward (Direction<0) replicas are
        # FIRST equilibrated at the backward ENDPOINT (λ=0, base=u1) to relax into
        # the u1 basin, THEN re-set to their assigned state. This mirrors the
        # per-direction structprep lesson (the ABFE dminus leg equilibrates from a
        # dminus base, not the dplus base).
        self.contexts = []
        self.integrators = []
        for k in range(self.n_states):
            integ = mm.LangevinMiddleIntegrator(
                temperature_K * unit.kelvin,
                friction_per_ps / unit.picosecond,
                timestep_fs * unit.femtoseconds)
            if seed is not None:
                integ.setRandomNumberSeed(int(seed) + k)
            ctx = mm.Context(system, integ, self.platform)
            state_k = self.replica_state[k]
            # FIX-A (A2): a context whose assigned state is in the decoupled band
            # is seeded from the PRE-RELAXED decoupled-basin config (so it starts
            # inside the decoupled basin, not behind the soft-core wall); every
            # other context keeps the default coupled ``positions``. DEFAULT
            # reseed_endpoint=False => _reseed_relaxed_positions is None and ALL
            # contexts get the coupled positions (byte-identical).
            if self._reseed_relaxed_positions is not None \
                    and state_k in self._reseed_band:
                ctx.setPositions(self._reseed_relaxed_positions)
            else:
                ctx.setPositions(positions)
            is_backward = (self.schedule["directions"][state_k] < 0)

            vel_seed = (int(seed) + k) if seed is not None else 0

            if self.staged_min:
                # OPT-IN staged relax (W4A large-box stability). The backward
                # u1-basin pre-equilibration is preserved (the backward anneal-edge
                # still relaxes into its endpoint basin first), but the minimization
                # at EACH state is staged: reference-state (soft-core OFF) >=5000-iter
                # relax of the fresh PME shell, then a short polish minimize at the
                # state the context will occupy. This path is reached ONLY when
                # staged_min=True; default-off leaves the single-stage path below
                # byte-identical.
                if is_backward and self._backward_endpoint is not None \
                        and backward_equil_steps > 0:
                    self._staged_minimize_replica(ctx, self._backward_endpoint)
                    ctx.setVelocitiesToTemperature(
                        temperature_K * unit.kelvin, vel_seed)
                    integ.step(int(backward_equil_steps))
                    # Switch to the replica's assigned anneal-edge state and stage-
                    # relax there too (the assigned state's soft-core is active).
                    self._staged_minimize_replica(ctx, state_k)
                else:
                    self._staged_minimize_replica(ctx, state_k)
                # Re-seed velocities + a short warmup at the assigned state (a fresh
                # hot burst on a just-minimized large box is the detonation source,
                # not the box itself).
                ctx.setVelocitiesToTemperature(
                    temperature_K * unit.kelvin, vel_seed)
                if self.staged_warmup_steps > 0:
                    integ.step(int(self.staged_warmup_steps))
            elif is_backward and self._backward_endpoint is not None \
                    and backward_equil_steps > 0:
                # Relax into the u1 basin at the backward endpoint first.
                self._set_state(ctx, self._backward_endpoint)
                if minimize_iters and minimize_iters > 0:
                    mm.LocalEnergyMinimizer.minimize(
                        ctx, maxIterations=int(minimize_iters))
                ctx.setVelocitiesToTemperature(
                    temperature_K * unit.kelvin, vel_seed)
                integ.step(int(backward_equil_steps))
                # Now switch to the replica's assigned (anneal-edge) state.
                self._set_state(ctx, state_k)
            else:
                self._set_state(ctx, state_k)
                if minimize_iters and minimize_iters > 0:
                    mm.LocalEnergyMinimizer.minimize(
                        ctx, maxIterations=int(minimize_iters))
            if not self.staged_min:
                ctx.setVelocitiesToTemperature(
                    temperature_K * unit.kelvin, vel_seed)
            self.integrators.append(integ)
            self.contexts.append(ctx)

        self._apply_all_states()

    def _find_backward_endpoint(self):
        """Index of the backward (Direction<0) state with minimal λ1 (the
        backward ladder endpoint = the u1-decoupled basin). Returns None if the
        ladder has no backward states."""
        s = self.schedule
        best = None
        best_lam = None
        for k in range(self.n_states):
            if s["directions"][k] < 0:
                if best is None or s["lambdas_1"][k] < best_lam:
                    best = k
                    best_lam = s["lambdas_1"][k]
        return best

    # -- per-state ATMForce global parameters ------------------------------
    def _set_state(self, ctx, state_idx):
        s = self.schedule
        alpha_per_kj = (s["alpha"][state_idx] / ats._kcal_to_kj(1.0)
                        if s["alpha"][state_idx] else 0.0)
        ctx.setParameter(self.atmforce.Lambda1(), s["lambdas_1"][state_idx])
        ctx.setParameter(self.atmforce.Lambda2(), s["lambdas_2"][state_idx])
        ctx.setParameter(self.atmforce.Alpha(), alpha_per_kj)
        ctx.setParameter(self.atmforce.Uh(), ats._kcal_to_kj(s["u0"][state_idx]))
        ctx.setParameter(self.atmforce.W0(), ats._kcal_to_kj(s["w0"][state_idx]))
        ctx.setParameter(self.atmforce.Umax(), self.umax_kj)
        ctx.setParameter(self.atmforce.Ubcore(), self.ubcore_kj)
        ctx.setParameter(self.atmforce.Acore(), self.acore)
        ctx.setParameter(self.atmforce.Direction(),
                         float(s["directions"][state_idx]))

    def _apply_all_states(self):
        for r in range(self.n_states):
            self._set_state(self.contexts[r], self.replica_state[r])

    # -- FIX-A (A2): standalone decoupled-endpoint equilibration (OPT-IN) ---
    def _equilibrate_endpoint(self, system, positions, endpoint_state, *,
                              seed=None, minimize_iters=500, equil_steps=2000):
        """Relax a STANDALONE context at ``endpoint_state``'s own λ-tuple and
        return the relaxed positions (the pre-relaxed decoupled-basin config A2
        seeds the band contexts with). Mirrors the validated W4A endpoint re-seed
        (W4A/w4a_reseed_proto.py / w4a_mixing_smoke.apply_fixA_dminus_endpoint_reseed):

          (a) a fresh context at the decoupled endpoint state's full ATMForce
              globals (NOT an artificial / biased config — a genuine short MD
              equilibration AT that state);
          (b) a light minimize (the same budget as the per-state seed minimize)
              to relieve the fresh-shell strain;
          (c) setVelocitiesToTemperature + a short MD equilibration to settle the
              walker into the decoupled basin.

        It uses its OWN throwaway integrator/context (offset RNG seed so it does
        not correlate with any ladder replica) and returns ONLY positions — the
        ladder's own contexts adopt them, the standalone context is discarded.
        Used ONLY when ``reseed_endpoint`` is True (default path never calls it).
        """
        mm = self.mm
        unit = self.unit
        integ = mm.LangevinMiddleIntegrator(
            self.temperature_K * unit.kelvin, 1.0 / unit.picosecond,
            1.0 * unit.femtoseconds)
        if seed is not None:
            integ.setRandomNumberSeed(int(seed) + 9973)
        relax_ctx = mm.Context(system, integ, self.platform)
        relax_ctx.setPositions(positions)
        self._set_state(relax_ctx, endpoint_state)
        if minimize_iters and minimize_iters > 0:
            mm.LocalEnergyMinimizer.minimize(
                relax_ctx, maxIterations=int(minimize_iters))
        relax_ctx.setVelocitiesToTemperature(
            self.temperature_K * unit.kelvin,
            (int(seed) + 9973) if seed is not None else 0)
        if equil_steps and equil_steps > 0:
            integ.step(int(equil_steps))
        relaxed = relax_ctx.getState(getPositions=True).getPositions()
        del relax_ctx, integ
        return relaxed

    # -- staged minimization (OPT-IN, W4A large-box stability) -------------
    def _staged_minimize_replica(self, ctx, state_idx):
        """Stage-relax one replica context for the assigned ``state_idx`` (OPT-IN).

        Reproduces the validated W4A/w4a_bound_smoke.py staged-minimization:

          (a) set the FULL per-state ATMForce globals for ``state_idx`` (so Uh / W0
              / Alpha / Direction / Umax / Ubcore / Acore match the assigned state),
              then OVERRIDE Lambda1=Lambda2=0 to put the context at the REFERENCE
              state where the ATM potential is the plain physical energy of the two
              resident copies (soft-core hybrid OFF, well-conditioned), and minimize
              there with the staged floor (>=5000 iters, tolerance->0) to relieve the
              fresh PME water-shell contacts;
          (b) restore the assigned state's Lambda1/Lambda2 (the soft-core hybrid is
              now active) and do a short polish minimize (max(500, floor//5) iters)
              to settle the alchemical region WITHOUT re-introducing the large-box
              fresh-shell strain already relieved in (a).

        The soft-core canon (Umax/Ubcore/Acore) and every other ATM global are the
        per-state schedule values — this is a MINIMIZATION-QUALITY relax only, it
        does NOT change the soft-core constants or the decouple physics. Only called
        when ``self.staged_min`` is True (the default-off ladder never invokes it).
        """
        mm = self.mm
        s = self.schedule
        # (a) Reference-state relax: assigned per-state globals, but Lambda1=Lambda2=0.
        self._set_state(ctx, state_idx)
        ctx.setParameter(self.atmforce.Lambda1(), 0.0)
        ctx.setParameter(self.atmforce.Lambda2(), 0.0)
        mm.LocalEnergyMinimizer.minimize(
            ctx, 0.0, int(self.staged_min_iters))
        # (b) Polish at the assigned state (restore the state's λ-tuple).
        ctx.setParameter(self.atmforce.Lambda1(), s["lambdas_1"][state_idx])
        ctx.setParameter(self.atmforce.Lambda2(), s["lambdas_2"][state_idx])
        mm.LocalEnergyMinimizer.minimize(
            ctx, 0.0, max(500, int(self.staged_min_iters) // 5))

    # -- per-replica raw perturbation (u0, u1-u0) --------------------------
    def _raw_pert(self, ctx):
        """Return (u0_kj, pert_kj) for the replica's CURRENT configuration."""
        pe = self.atmforce.getPerturbationEnergy(ctx)
        # OpenMM 8.4 ATMForce.getPerturbationEnergy -> (u1, u0, energy) tuple of
        # Quantity. The perturbation u = u1 - u0.
        u1 = pe[0].value_in_unit(self.unit.kilojoule_per_mole)
        u0 = pe[1].value_in_unit(self.unit.kilojoule_per_mole)
        return u0, (u1 - u0)

    def _state_energy_kj(self, u0_kj, pert_kj, state_idx):
        s = self.schedule
        alpha_per_kj = (s["alpha"][state_idx] / ats._kcal_to_kj(1.0)
                        if s["alpha"][state_idx] else 0.0)
        return _atm_state_energy_kj(
            pert_kj=pert_kj, u0_kj=u0_kj,
            lambda1=s["lambdas_1"][state_idx], lambda2=s["lambdas_2"][state_idx],
            alpha_per_kj=alpha_per_kj,
            uh_kj=ats._kcal_to_kj(s["u0"][state_idx]),
            w0_kj=ats._kcal_to_kj(s["w0"][state_idx]),
            direction=float(s["directions"][state_idx]),
            umax_kj=self.umax_kj, ubcore_kj=self.ubcore_kj, acore=self.acore)

    # -- one RE cycle ------------------------------------------------------
    def run_cycle(self, md_steps=250, exchange_attempts=None):
        """Propagate every replica ``md_steps``, then attempt NN exchanges.

        Returns a per-cycle dict (any NaN flagged). Emits one
        ``Replica r new state s`` log line per replica AFTER the exchange sweep
        (the per-direction driver's parser keys on these lines).
        """
        import random
        # 1) Propagate every replica at its current state. A per-replica NaN
        #    (OpenMM raises "Particle coordinate is NaN" on step) is CAUGHT and
        #    FLAGGED per replica rather than aborting the whole run — so the
        #    smoke reports WHICH ladder state(s) are unstable (an honest
        #    empirical signal: an unstable window needs a schedule change, not a
        #    hidden crash). A replica that NaNs is marked dead for the cycle and
        #    excluded from the exchange sweep (its energy is unusable).
        nan_seen = False
        nan_states = []
        dead = [False] * self.n_states
        for r in range(self.n_states):
            try:
                self.integrators[r].step(md_steps)
                st = self.contexts[r].getState(getEnergy=True)
                pe = st.getPotentialEnergy().value_in_unit(
                    self.unit.kilojoule_per_mole)
                if math.isnan(pe) or math.isinf(pe):
                    dead[r] = True
            except Exception:                       # OpenMM NaN on step()
                dead[r] = True
            if dead[r]:
                nan_seen = True
                nan_states.append(self.replica_state[r])

        # 2) Cache each replica's raw (u0, pert) at its current configuration.
        #    A dead replica has no usable energy (None sentinel).
        raw = []
        for r in range(self.n_states):
            if dead[r]:
                raw.append(None)
            else:
                raw.append(self._raw_pert(self.contexts[r]))

        # 3) Nearest-neighbor state exchanges (Metropolis). Alternate even/odd
        #    adjacent state pairs each cycle to avoid systematic bias.
        if exchange_attempts is None:
            exchange_attempts = self.n_states
        # Map state -> replica currently in it.
        state_to_replica = {self.replica_state[r]: r
                            for r in range(self.n_states)}
        n_accepted = 0
        parity = self._cycle % 2
        for s_lo in range(parity, self.n_states - 1, 2):
            s_hi = s_lo + 1
            r_lo = state_to_replica[s_lo]
            r_hi = state_to_replica[s_hi]
            # Skip the exchange if either replica NaN'd this cycle (no usable
            # energy) — its state can still be visited again next cycle.
            if raw[r_lo] is None or raw[r_hi] is None:
                continue
            u0_lo, pert_lo = raw[r_lo]
            u0_hi, pert_hi = raw[r_hi]
            # ΔΔ = [U_shi(x_lo) + U_slo(x_hi)] - [U_slo(x_lo) + U_shi(x_hi)]
            e_lo_at_lo = self._state_energy_kj(u0_lo, pert_lo, s_lo)
            e_lo_at_hi = self._state_energy_kj(u0_lo, pert_lo, s_hi)
            e_hi_at_hi = self._state_energy_kj(u0_hi, pert_hi, s_hi)
            e_hi_at_lo = self._state_energy_kj(u0_hi, pert_hi, s_lo)
            delta = ((e_lo_at_hi + e_hi_at_lo) - (e_lo_at_lo + e_hi_at_hi))
            if (not math.isnan(delta)) and (
                    delta <= 0.0 or random.random() < math.exp(-delta / self.kT_kj)):
                # Swap the two replicas' states.
                self.replica_state[r_lo] = s_hi
                self.replica_state[r_hi] = s_lo
                state_to_replica[s_lo] = r_hi
                state_to_replica[s_hi] = r_lo
                n_accepted += 1

        # 4) Re-apply states (a replica that changed state needs its λ updated).
        self._apply_all_states()

        # 5) Emit the driver-compatible log lines (one per replica). The
        #    per-direction driver's parser groups samples into cycles by DISTINCT
        #    HH:MM:SS timestamp (its warmup boundary), so we synthesise a
        #    monotonically-increasing per-cycle timestamp (base + cycle seconds)
        #    rather than the wall clock — a fast (GPU/CPU) run can emit several
        #    cycles within the same wall-clock second, which would otherwise
        #    collapse them into ONE cycle and defeat the warmup/transition logic.
        self._cycle += 1
        if self._log_fh is not None:
            import datetime
            ts = (self._log_base_ts
                  + datetime.timedelta(seconds=self._cycle)
                  ).strftime("%Y-%m-%d %H:%M:%S")
            for r in range(self.n_states):
                self._log_fh.write(
                    "%s Replica %d new state %d\n"
                    % (ts, r, self.replica_state[r]))
            self._log_fh.flush()

        # 6) UWHAM-consumable per-walker .out row (production). Each walker r
        #    records the sample it just propagated, labelled with the state it
        #    NOW occupies (after the exchange sweep — the abfe_production
        #    convention: the stateid column is the walker's post-swap state, the
        #    energies are this cycle's configuration). A walker that NaN'd this
        #    cycle has no usable energy (raw is None) so it writes NO row — UWHAM
        #    simply sees one fewer sample for that walker; the empirical NaN
        #    signal is already surfaced via nan_states (NOT hidden). potE / pertE
        #    are written in kcal/mol (the UWHAM estimator's expected units).
        if self._out_fhs is not None:
            s = self.schedule
            kcal = ats._kcal_to_kj(1.0)
            for r in range(self.n_states):
                if raw[r] is None:
                    continue
                u0_kj, pert_kj = raw[r]
                st = self.replica_state[r]
                # potE = the FULL ATM hybrid potential at the walker's current
                # state; pertE = the SOFT-CORE-CAPPED perturbation usc (NOT the
                # raw u1-u0). UWHAM's _bias_fcn operates DIRECTLY on the pertE
                # column with no cap of its own, so it MUST receive usc — a raw
                # u1-u0 of tens of thousands of kcal/mol would overflow the WHAM
                # solve (and would not reconstruct the engine's capped hybrid).
                # e0 = potE - bias_fcn(usc) then recovers the λ-independent base,
                # exactly as the upstream ommreplica.save_out convention.
                direction = float(s["directions"][st])
                _base_kj, usc_kj = _atm_softcore_components_kj(
                    pert_kj, u0_kj, direction,
                    self.umax_kj, self.ubcore_kj, self.acore)
                pot_kcal = self._state_energy_kj(u0_kj, pert_kj, st) / kcal
                usc_kcal = usc_kj / kcal
                self._out_fhs[r].write(
                    "%d %.1f %d %.6f %.6f %.6f %.6f %.6f %.10f %.10f 0\n"
                    % (st, self.temperature_K, int(direction),
                       s["lambdas_1"][st], s["lambdas_2"][st],
                       s["alpha"][st], s["u0"][st], s["w0"][st],
                       pot_kcal, usc_kcal))
            for fh in self._out_fhs:
                fh.flush()

        # 7) OPT-IN per-walker DCD frame (probe diagnostic). DEFAULT
        #    dcd_dir=None => self._dcd_files is None and this is skipped entirely
        #    (byte-identical). When on, append ONE frame of the recorded atom
        #    subset per walker every dcd_stride_cycles cycles, AFTER the .out row
        #    is written so frame index k maps to .out row k for that walker (the
        #    state label is the stateid column of the matching .out row). A walker
        #    with no usable energy this cycle (raw is None) wrote no .out row, so
        #    it writes no frame either — the two streams stay in lockstep. Reading
        #    getPositions here is observation only; it does not advance dynamics or
        #    alter the estimator.
        if self._dcd_files is not None \
                and (self._cycle % self.dcd_stride_cycles == 0):
            idx = self._dcd_atom_indices
            for r in range(self.n_states):
                if raw[r] is None:
                    continue
                fh, dcdf = self._dcd_files[r]
                st_full = self.contexts[r].getState(
                    getPositions=True, enforcePeriodicBox=True)
                all_pos = st_full.getPositions(asNumpy=True)
                box = st_full.getPeriodicBoxVectors()
                dcdf.writeModel(all_pos[idx], periodicBoxVectors=box)
            for fh, _dcdf in self._dcd_files:
                fh.flush()

        return {
            "cycle": self._cycle,
            "n_accepted": n_accepted,
            "n_pairs": (self.n_states - 1 + 1 - parity) // 2,
            "nan_seen": nan_seen,
            "nan_states": sorted(set(nan_states)),
            "replica_state": list(self.replica_state),
        }

    def close(self):
        # Order matters (data integrity before resource teardown): flush + close
        # the file handles FIRST so the last per-cycle row is durably written,
        # THEN release the OpenMM Contexts/Integrators. Each ladder holds
        # n_states live Contexts (~13 GiB/worker at peak); without explicit
        # release the dplus ladder's Contexts stay resident while the dminus
        # ladder is built, and a 4-concurrent pool over-commits host RAM. The
        # method is idempotent: it is safe to call twice, and safe when called
        # on a partially-constructed ladder (an __init__ that raised before the
        # Contexts were created), so both the fail-fast and the normal paths can
        # funnel through one teardown.
        log_fh = getattr(self, "_log_fh", None)
        if log_fh is not None:
            try:
                log_fh.flush()
            except Exception:
                pass
            log_fh.close()
            self._log_fh = None
        out_fhs = getattr(self, "_out_fhs", None)
        if out_fhs is not None:
            for fh in out_fhs:
                try:
                    fh.flush()
                except Exception:
                    pass
                try:
                    fh.close()
                except Exception:
                    pass
            self._out_fhs = None
        # OPT-IN DCD handles (probe diagnostic): flush + close the trajectory
        # byte streams so the last frame is durable. Each entry is an
        # (open file handle, DCDFile) tuple; only the file handle is closed (the
        # DCDFile is a thin writer over it). Idempotent + exception-swallowing,
        # matching the .out teardown above (data integrity before resource free).
        dcd_files = getattr(self, "_dcd_files", None)
        if dcd_files is not None:
            for fh, _dcdf in dcd_files:
                try:
                    fh.flush()
                except Exception:
                    pass
                try:
                    fh.close()
                except Exception:
                    pass
            self._dcd_files = None

        # Resource teardown: SWIG-backed OpenMM Context/Integrator objects are
        # not reliably reclaimed by gc.collect() alone (reference cycles through
        # the SWIG proxy can keep the underlying C++ object alive), so drop our
        # own references with explicit ``del`` per object before collecting.
        contexts = getattr(self, "contexts", None)
        if contexts:
            for _ctx in contexts:
                try:
                    del _ctx
                except Exception:
                    pass
        self.contexts = []
        integrators = getattr(self, "integrators", None)
        if integrators:
            for _integ in integrators:
                try:
                    del _integ
                except Exception:
                    pass
        self.integrators = []
        gc.collect()
