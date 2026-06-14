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
) -> Dict[str, Any]:
    """Build + serialize the in-place fused RBFE System for one leg.

    Writes ``<out_dir>/inplace_rbfe_<tag>_sys.xml`` (the OpenMM System,
    XmlSerializer, WITH the genuine ATMForce attached) and
    ``<out_dir>/inplace_rbfe_<tag>.pdb`` (the matching topology + coordinates +
    periodic box). ``tag`` defaults to ``<leg>``.

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
    import openmm as mm

    if tag is None:
        tag = leg
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

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

    with open(sys_xml, "w") as fh:
        fh.write(mm.XmlSerializer.serialize(system))

    # Topology/position lockstep: the genuine swap adds a MASSLESS dummy NE1
    # reference particle to the System + NonbondedForce + positions (the minimal
    # in-place stand-in for the two-copy overlay's distinct partner-attach), but
    # NOT to the topology (it has no chemical identity). So System particle count
    # == len(positions) == topology atoms + 1. PDBFile.writeFile requires
    # topology and positions to match, so we materialise a matching topology atom
    # for each extra (massless) System particle, in a dedicated "X" chain (always
    # contiguity-legal). This is serialization-only bookkeeping — it does NOT
    # touch the validated build's System/energy (R-7); the dummy already exists in
    # the System, we only give it a topology slot so the PDB round-trips.
    import openmm.unit as unit
    from openmm.app import PDBFile, element as _app_element
    topology = modeller.topology
    positions = modeller.positions   # a unit-wrapped position list (Quantity)
    n_extra = system.getNumParticles() - topology.getNumAtoms()
    if n_extra < 0:
        raise RuntimeError(
            "serialize_inplace_rbfe_system: topology (%d atoms) exceeds System "
            "particles (%d) — unexpected build state."
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

    n_atoms = topology.getNumAtoms()
    return {
        "leg": leg,
        "tag": tag,
        "seed": seed,
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
                 out_dir=None, out_basename="trackb"):
        import openmm as mm
        import openmm.unit as unit

        self.mm = mm
        self.unit = unit
        self.schedule = schedule
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

        # replica r currently occupies state self.replica_state[r].
        self.replica_state = list(range(self.n_states))

        # Identify the backward (Direction<0) ladder endpoint state (λ minimal on
        # the backward half) — the u1 = HE1-decoupled basin the backward replicas
        # must relax into BEFORE the soft-core anneal-edge.
        self._backward_endpoint = self._find_backward_endpoint()

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
            ctx.setPositions(positions)
            state_k = self.replica_state[k]
            is_backward = (self.schedule["directions"][state_k] < 0)

            if is_backward and self._backward_endpoint is not None \
                    and backward_equil_steps > 0:
                # Relax into the u1 basin at the backward endpoint first.
                self._set_state(ctx, self._backward_endpoint)
                if minimize_iters and minimize_iters > 0:
                    mm.LocalEnergyMinimizer.minimize(
                        ctx, maxIterations=int(minimize_iters))
                ctx.setVelocitiesToTemperature(
                    temperature_K * unit.kelvin,
                    (int(seed) + k) if seed is not None else 0)
                integ.step(int(backward_equil_steps))
                # Now switch to the replica's assigned (anneal-edge) state.
                self._set_state(ctx, state_k)
            else:
                self._set_state(ctx, state_k)
                if minimize_iters and minimize_iters > 0:
                    mm.LocalEnergyMinimizer.minimize(
                        ctx, maxIterations=int(minimize_iters))
            ctx.setVelocitiesToTemperature(
                temperature_K * unit.kelvin,
                (int(seed) + k) if seed is not None else 0)
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
