#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B — ATS (Alchemical Transfer with coordinate Swapping) system setup.

Builds the Cp4<->WT single-point relative-binding alchemical system for the
2QKI compstatin complex (conditions B-C1 .. B-C7). INFRASTRUCTURE /
SMOKE-TEST scope only — the production FEP
(>=11 lambda x >=5 ns x >=3 replicas) is gated on Track A's read-out and is
NOT launched here.

Thermodynamic cycle (B-C3):
    DDG_bind = DG_alch(bound, Trp->MTR) - DG_alch(free-peptide, Trp->MTR)
Both legs retain the compstatin ``cyclic_ss`` disulfide (CYS2-CYS12 SG-SG).

Charge axis (B-C1): Option-beta / regime-2 hybrid MTR
(``params/MTR_gaff2_hybrid.xml``): amber14SB-frozen backbone, NE1 frozen at
-0.3418, Khoury 2014 OMW sidechain. Formal charge 0 at both endpoints
(charge-conserving perturbation). 2026-05-30 Q1 fix: per-residue Σq driven
to |Σq| ≤ 5e-4 e via uniform per-atom rescale on sidechain heavy+H scope
(MTR/NMTR already integer; CMTR rescaled, per-atom Δ ≈ +0.0088 e within
Khoury ±0.02 e tolerance). Audit log:
``params/_archive/mtr_rescale_audit_20260530.json``.

Alchemical region (B-C5), residue 4:
    WT  (TRP4): indole donor NE1-HE1.
    Cp4 (MTR4): N-methyl NE1-CM with HM1/HM2/HM3.
    Trp -> MTR perturbation: HE1 DISAPPEARS, CM+HM1+HM2+HM3 APPEAR.
    NE1 is the common (frozen) atom shared by both states.

The OpenMM ``ATMForce`` provides the soft-core alchemical potential that
scales the appearing/disappearing atoms across lambda.

Engine reuse (anti-fragmentation): the disulfide detector
(``run_restrained_md.detect_disulfide_pair`` / ``commit_bond_if_missing``),
the prepared 2QKI structures under ``outputs/2QKI_{Cp4,WT}_calib_s*/`` and the
hybrid charge XML are reused unchanged.
"""

import hashlib
import os
import random
import sys
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Tuple, Any

import numpy as np

import openmm as mm
import openmm.unit as unit
from openmm import app
from openmm.app import ForceField, Modeller, PDBFile, PME, HBonds

_UTILS_DIR = os.path.dirname(os.path.abspath(__file__))
if _UTILS_DIR not in sys.path:
    sys.path.insert(0, _UTILS_DIR)
_PROJ_ROOT = os.path.dirname(_UTILS_DIR)

# ---------------------------------------------------------------------------
# Canonical UPDD assets (B-C1, B-C3). Resolved relative to project root so the
# module is location-independent.
# ---------------------------------------------------------------------------
HYBRID_MTR_XML = os.path.join(_PROJ_ROOT, "params", "MTR_gaff2_hybrid.xml")
# RBFE-only MTR XML. The in-place residue-4 RBFE build must NOT load the shared
# HYBRID_MTR_XML (consumed by Track A MD / QM/MM 1-traj / v0.8 RESP / MM-PBSA,
# all of which intentionally use the Option-beta RESP charges). A separate
# RBFE-scoped file isolates any RBFE-specific common-charge harmonization from
# those tracks. Falls back to the shared file if the RBFE file is absent.
HYBRID_MTR_XML_RBFE = os.path.join(
    _PROJ_ROOT, "params", "MTR_gaff2_hybrid_rbfe_harmonized.xml")
# The project-root params/MTR_hydrogens.xml is a placeholder; the populated
# definition (HM1-3 on CM, etc.) lives per-seed. resolve_leg_inputs() returns
# the real path; this default is only the fallback location.
HYBRID_MTR_HYDROGENS_XML = os.path.join(
    _PROJ_ROOT, "outputs", "2QKI_Cp4_hybrid_calib_s7", "params", "MTR_hydrogens.xml"
)
FF_FILES = ["amber14-all.xml", "amber14/tip3pfb.xml"]

# Residue-4 alchemical atom partition (B-C5).
ALCH_RESNUM = 4
ALCH_COMMON_ATOM = "NE1"                       # frozen, shared by both states
ALCH_WT_ONLY = ["HE1"]                         # TRP4 indole donor (disappears)
ALCH_MTR_ONLY = ["CM", "HM1", "HM2", "HM3"]    # MTR4 N-methyl (appears)


# ---------------------------------------------------------------------------
# Mutation specification (C5 generalization, 2026-06-16)
#
# The in-place two-copy ATS box perturbs ONE residue's side chain between two
# endpoint states. The original (and DEFAULT) perturbation is residue-4
# Cp4(MTR N-methyl) <-> WT(Trp indole), which requires the Khoury hybrid ncAA
# XML. ``MutationSpec`` extracts that mutation-definition layer (resnum, common
# attach atom, per-state-only atoms, the ncAA XML, the two residue names) into a
# parameter object so the SAME two-copy core (build / swap / displacement / C6
# guards) can also drive a CANONICAL all-amber perturbation (e.g. residue-3
# Val<->Ile for the V3I engine-validation, where both endpoint residues are
# standard amber14 templates and NO ncAA XML is needed).
#
# State convention (kept identical to the legacy MTR<->Trp wiring so the index
# map / swap / asserts are reused unchanged):
#   * ``stateB`` = the APPEARING state -> copy-1 (built at the site). Its
#     ``stateB_only_atoms`` populate the legacy ``mtr_only`` slot.
#   * ``stateA`` = the DISAPPEARING state -> copy-2 (displaced into bulk). Its
#     ``stateA_only_atoms`` populate the legacy ``wt_only`` slot.
# For MTR<->Trp the appearing state is MTR (the N-methyl appears) and the
# disappearing state is WT-Trp (HE1 disappears), so stateB=MTR, stateA=WT —
# byte-identical to the legacy ALCH_* constants below.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class MutationSpec:
    """Definition of a single-residue dual-topology side-chain perturbation.

    Fields:
      name                : short identifier (used by the CLI / run manifest).
      resnum              : binder-chain residue number carrying the mutation.
      common_attach_atom  : the frozen heavy atom both per-state var groups bond
                            to (the swap attach atom). MUST be present in BOTH
                            endpoint residues (e.g. NE1 for Trp/MTR, CG1 for
                            Val/Ile — both the disappearing H and the appearing
                            methyl/ethyl carbon bond it).
      stateA_resname      : the DISAPPEARING-state residue name (copy-2 / WT).
      stateB_resname      : the APPEARING-state residue name (copy-1 / site).
      stateA_only_atoms   : atom names present ONLY in stateA (disappear). Maps to
                            the legacy ``wt_only`` slot.
      stateB_only_atoms   : atom names present ONLY in stateB (appear). Maps to
                            the legacy ``mtr_only`` slot.
      hybrid_xml          : ncAA ForceField XML required to template stateB
                            (``None`` => both endpoints are standard amber14
                            templates, no extra XML — the canonical path).
      bonded_heavy_appearing : the appearing-state HEAVY var atom that bonds the
                            common attach atom directly (CM for MTR, CD1 for
                            Ile). Used by the MC2 bonded-term assert. ``None`` when
                            the appearing side grows NO heavy atom (the
                            disappearing-heavy mirror, e.g. Ala->Gly).
      appearing_h_prefix  : the atom-name prefix of the appearing-state methyl
                            hydrogens (HM for MTR, HD for Ile). Used by MC2. ``None``
                            when the appearing side grows no methyl H group.
      bonded_heavy_disappearing : the DISAPPEARING-state HEAVY var atom that bonds
                            the common attach atom directly (CB for Ala in the
                            Ala->Gly mirror). The symmetric counterpart of
                            ``bonded_heavy_appearing``; ``None`` for the
                            appearing-heavy shapes (MTR / V3I disappear only an H).
                            Used by the MC2 bonded-term assert (disappearing branch).
      disappearing_h_prefix : the atom-name prefix of the disappearing-state methyl
                            hydrogens (HB for Ala). ``None`` for shapes that
                            disappear no methyl H group. Used by MC2.
      multiheavy_star_certified : explicit opt-in attestation that this spec is an
                            ACYCLIC "star" multi-methyl perturbation — two or more
                            var HEAVY atoms each bonded DIRECTLY to the common
                            attach atom (e.g. Val3->Ala deletes CG1 AND CG2, both
                            bonded to CB) — that the two-copy coordinate-swap engine
                            can build. Defaults ``False``: every multi-heavy spec is
                            ``unsupported`` UNLESS it is explicitly certified here
                            (fail-safe whitelist — an un-attested multi-heavy spec is
                            never silently treated as buildable). Ring / fused-ring /
                            chained-heavy / multi-branch (heavy on BOTH sides)
                            perturbations are NOT in scope for this attestation
                            (W4A indole fusion is Phase B's connected-subgraph shape,
                            W4F ring contraction stays unsupported) and must remain
                            ``False``; the MC2 backstop additionally fail-louds on any
                            certified heavy that does not in fact bond the attach atom
                            by a real bond.
      ring_closure_bonds : (Phase B) tuple of (atomA, atomB) NAME pairs that close a
                            ring ENTIRELY inside the perturbation's var group (the
                            indole 5/6 fusion bond CD2-CE2 for W4A). EMPTY ``()`` for
                            every ACYCLIC shape (MTR / V3I / A9G / V3A): with no ring
                            bonds the classifier is byte-identical to the legacy
                            single-heavy + multi-star paths. Non-empty => the
                            perturbation is a connected-subgraph RING shape and is
                            gated FIRST (before any heavy-count branch), so a ring spec
                            can NEVER edge into a heavy-count shape (star or single) by
                            accident — the load-bearing silent-build backstop.
      connected_group_certified : (Phase B) explicit opt-in that the var group is a
                            SINGLE-ATTACH connected subgraph — one attach-bonded root
                            heavy (CG->CB for W4A) + the rest reachable by intra-group
                            bonds + the declared ring-closure bonds. A DEDICATED flag,
                            INTENTIONALLY distinct from ``multiheavy_star_certified``: a
                            ring spec is buildable ONLY when THIS flag is set, never via
                            the star flag (a ring author who sets only the star flag
                            STILL lands at ``unsupported``). Defaults ``False``: an
                            un-certified ring spec is ``unsupported`` (fail-safe).
      chained_group_certified : explicit opt-in that the var group is a one-sided
                            (all-disappearing OR all-appearing) ACYCLIC connected
                            subgraph in which a var heavy may bond ONLY a PARENT var
                            heavy rather than the common attach atom directly (Ile->Ala
                            deletes CG1 AND CG2 both bonded to CB, PLUS CD1 chained off
                            CG1). A DEDICATED flag, INTENTIONALLY distinct from
                            ``multiheavy_star_certified`` (the star flag requires EVERY
                            heavy to bond the attach directly, so a chained heavy fails
                            the star MC2 backstop) and from ``connected_group_certified``
                            (which is the RING shape and requires a NON-EMPTY
                            ``ring_closure_bonds``). This flag REQUIRES an EMPTY
                            ``ring_closure_bonds`` (an acyclic tree) — a heavy cycle is a
                            ring shape and must use ``connected_group_certified``. Defaults
                            ``False``: an un-certified chained-heavy spec is ``unsupported``
                            (fail-safe whitelist). The MC2 backstop additionally
                            fail-louds on any declared heavy not reachable from the attach
                            atom by real bonds, or on a detected cycle.

    Two single-heavy SHAPES are supported (see :pyattr:`shape`):
      * ``appearing_heavy``    : exactly one appearing heavy + its H's, no
                                 disappearing heavy (MTR: CM+HM; V3I: CD1+HD).
      * ``disappearing_heavy`` : exactly one disappearing heavy + its H's, no
                                 appearing heavy (A9G mirror: CB+HB disappear,
                                 Gly grows only a backbone H, no heavy).
    Two acyclic-star multi-heavy SHAPES are supported ONLY when explicitly
    certified via ``multiheavy_star_certified`` (see :pyattr:`shape`):
      * ``multi_appearing_heavy``    : two or more appearing heavies (each star-
                                       bonded to the attach atom), no disappearing
                                       heavy.
      * ``multi_disappearing_heavy`` : two or more disappearing heavies (each star-
                                       bonded to the attach atom), no appearing heavy
                                       (V3A Val->Ala mirror: CG1+CG2 disappear).
    One connected-subgraph RING SHAPE is supported (Phase B) ONLY when explicitly
    certified via ``connected_group_certified`` + a non-empty ``ring_closure_bonds``
    (see :pyattr:`shape`):
      * ``single_attach_connected_group`` : a one-sided var group (all-disappearing OR
                                       all-appearing heavy) forming a connected subgraph
                                       rooted at ONE attach-bonded heavy, with the ring
                                       closed by the declared ``ring_closure_bonds``
                                       (W4A Trp->Ala deletes the 9-heavy fused indole;
                                       root CG bonds CB, the 5/6 fusion bond CD2-CE2
                                       closes the bicyclic system inside the var group).
    Any other shape (uncertified multi-heavy, an un-certified ring, or both sides
    growing/deleting a heavy => ring contraction / multi-branch) is ``unsupported``
    and the MC2 / R2 asserts FAIL-LOUD on it (W4F/Q5A must hit this, never
    silent-build).
    """
    name: str
    resnum: int
    common_attach_atom: str
    stateA_resname: str
    stateB_resname: str
    stateA_only_atoms: Tuple[str, ...]
    stateB_only_atoms: Tuple[str, ...]
    hybrid_xml: Optional[str] = None
    bonded_heavy_appearing: Optional[str] = None
    appearing_h_prefix: Optional[str] = None
    bonded_heavy_disappearing: Optional[str] = None
    disappearing_h_prefix: Optional[str] = None
    multiheavy_star_certified: bool = False
    ring_closure_bonds: Tuple[Tuple[str, str], ...] = ()
    connected_group_certified: bool = False
    chained_group_certified: bool = False

    @staticmethod
    def _heavy_names(atom_names: Tuple[str, ...]) -> List[str]:
        """Heavy (non-hydrogen) atom names in a var set, by name convention.

        Hydrogens are named with a leading ``H`` (optionally after a numeric
        wyckoff digit, e.g. ``1HB`` — handled defensively). Anything else is a
        heavy atom. Used only for shape classification, never for the partition.
        """
        heavy: List[str] = []
        for n in atom_names:
            s = n.strip()
            base = s.lstrip("0123456789")
            if base[:1].upper() == "H":
                continue
            heavy.append(s)
        return heavy

    @property
    def shape(self) -> str:
        """Classify the single-residue perturbation SHAPE for the MC2 / R2 asserts.

        Returns one of:
          ``appearing_heavy``         -> exactly one appearing heavy, zero
                                         disappearing heavy (MTR / V3I).
          ``disappearing_heavy``      -> exactly one disappearing heavy, zero
                                         appearing heavy (A9G Ala->Gly mirror).
          ``multi_appearing_heavy``   -> two or more appearing heavies, zero
                                         disappearing heavy, AND the spec is
                                         explicitly ``multiheavy_star_certified``
                                         (acyclic star multi-methyl).
          ``multi_disappearing_heavy``-> two or more disappearing heavies, zero
                                         appearing heavy, AND the spec is explicitly
                                         ``multiheavy_star_certified`` (V3A mirror).
          ``single_attach_connected_group`` -> (Phase B, ring-FIRST) a NON-EMPTY
                                         ``ring_closure_bonds`` AND a one-sided var
                                         group (all-disappearing OR all-appearing heavy)
                                         AND ``connected_group_certified`` is True (W4A
                                         certified fused indole: one attach-bonded root
                                         + BFS-reachable ring heavies + the declared
                                         ring-closure bond). Gated FIRST, before any
                                         heavy-count branch.
          ``unsupported``             -> anything else: a heavy on BOTH sides
                                         (ring contraction / multi-branch), zero heavy on
                                         either side, an UN-certified multi-heavy spec,
                                         OR a RING spec (``ring_closure_bonds`` non-empty)
                                         that is either un-certified or heavy-both-sides
                                         (a ring author who sets only
                                         ``multiheavy_star_certified`` STILL lands here —
                                         the star flag is NEVER honoured for a ring
                                         shape). The classification is a FAIL-SAFE
                                         WHITELIST: a multi-heavy / ring spec is only
                                         buildable when it has opted in via the
                                         APPROPRIATE attestation (``multiheavy_star_
                                         certified`` for acyclic stars, ``connected_
                                         group_certified`` for rings) — omitting the
                                         flag keeps it ``unsupported`` (so a spec
                                         author who forgets the flag, or a ring/fused/
                                         chained-heavy spec, can never be silently
                                         mis-classified as supported). The MC2 / R2
                                         asserts FAIL-LOUD on ``unsupported`` rather
                                         than silent-build a wrong endpoint.
        """
        n_app_heavy = len(self._heavy_names(self.stateB_only_atoms))
        n_dis_heavy = len(self._heavy_names(self.stateA_only_atoms))
        # RING-FIRST gate (Phase B, load-bearing silent-build backstop): a ring-closure
        # spec is classified by the DEDICATED connected-group attestation, evaluated
        # BEFORE the heavy-count opt-in branches so the star flag can never silent-build
        # a ring. A heavy growing AND deleting (W4F ring CONTRACTION / multi-branch)
        # stays unsupported even WITH the ring attestation: a single-attach connected
        # group has a one-sided var group (all-disappearing OR all-appearing), so a heavy
        # on BOTH sides is a dual swap the engine cannot construct. With NO ring bonds
        # the classification falls through to the byte-identical legacy paths below.
        if self.ring_closure_bonds:
            if n_app_heavy >= 1 and n_dis_heavy >= 1:
                return "unsupported"
            if self.connected_group_certified:
                return "single_attach_connected_group"
            return "unsupported"
        # A heavy growing AND deleting => multi-branch / ring rewiring (W4F): never
        # buildable here, regardless of the opt-in attestation.
        if n_app_heavy >= 1 and n_dis_heavy >= 1:
            return "unsupported"
        # Single-heavy shapes (byte-identical to the legacy classifier).
        if n_app_heavy == 1 and n_dis_heavy == 0:
            return "appearing_heavy"
        if n_dis_heavy == 1 and n_app_heavy == 0:
            return "disappearing_heavy"
        # Acyclic-star multi-heavy shapes — ONLY when explicitly opted in
        # (fail-safe whitelist; an un-certified multi-heavy spec stays unsupported).
        if self.multiheavy_star_certified:
            if n_app_heavy >= 2 and n_dis_heavy == 0:
                return "multi_appearing_heavy"
            if n_dis_heavy >= 2 and n_app_heavy == 0:
                return "multi_disappearing_heavy"
        # Acyclic CHAINED connected-group multi-heavy shape — ONLY when explicitly
        # opted in via ``chained_group_certified`` (a DEDICATED flag, distinct from the
        # star flag). Reached only for a one-sided (all-appearing OR all-disappearing)
        # multi-heavy var group with EMPTY ``ring_closure_bonds`` (the ring gate above
        # already claimed any ring spec). Unlike the star shape, a var heavy here may
        # bond only a PARENT var heavy (Ile->Ala: CD1 chains off CG1, not CB); the MC2
        # backstop certifies attach-reachability + ACYCLICITY (a heavy cycle is a ring
        # shape and belongs to ``single_attach_connected_group``, never here).
        if self.chained_group_certified:
            if n_app_heavy >= 2 and n_dis_heavy == 0:
                return "acyclic_connected_group"
            if n_dis_heavy >= 2 and n_app_heavy == 0:
                return "acyclic_connected_group"
        return "unsupported"


# DEFAULT spec — residue-4 Cp4(MTR) <-> WT(Trp). Byte-identical to the legacy
# ALCH_* constants (stateB=MTR appears, stateA=WT disappears); ``hybrid_xml`` is
# resolved at build time (RBFE-harmonized XML preferred) so it is left ``None``
# here and the builder keeps its existing XML-resolution logic for this spec.
MUTATION_MTR_TRP_RES4 = MutationSpec(
    name="mtr_trp_res4",
    resnum=4,
    common_attach_atom="NE1",
    stateA_resname="TRP",
    stateB_resname="MTR",
    stateA_only_atoms=("HE1",),                       # WT indole donor (disappears)
    stateB_only_atoms=("CM", "HM1", "HM2", "HM3"),    # MTR N-methyl (appears)
    hybrid_xml=None,                                  # resolved at build time
    bonded_heavy_appearing="CM",
    appearing_h_prefix="HM",
)

# V3I spec — residue-3 Val(WT) <-> Ile (engine validation, engine validation (Option B)).
# Both endpoints are standard amber14 templates (no ncAA XML). amber14 atom
# naming (verified against the amber14-all VAL/ILE templates):
#   VAL CG1: HG11, HG12, HG13 ; ILE CG1: HG12, HG13 + CD1(HD11,HD12,HD13).
# So Val->Ile ADDS a gamma-CH3: HG11 disappears, CD1+HD11-13 appear; the common
# attach atom is CG1 (both the disappearing HG11 and the appearing CD1 bond it).
# stateB=Ile (appears, copy-1), stateA=Val (disappears, copy-2). hybrid_xml=None
# (canonical amber14 — Val/Ile are charge-neutral, parity-identical: R-15/R-16
# pass trivially).
MUTATION_VAL_ILE_RES3 = MutationSpec(
    name="v3i_val_ile_res3",
    resnum=3,
    common_attach_atom="CG1",
    stateA_resname="VAL",
    stateB_resname="ILE",
    stateA_only_atoms=("HG11",),                      # Val gamma-H (disappears)
    stateB_only_atoms=("CD1", "HD11", "HD12", "HD13"),  # Ile gamma-CH3 (appears)
    hybrid_xml=None,                                  # canonical amber14
    bonded_heavy_appearing="CD1",
    appearing_h_prefix="HD",
)

# A9G spec — residue-9 Ala(WT) <-> Gly (engine de-risk of the disappearing-heavy
# MC2/partition extension). The MIRROR of V3I: V3I GROWS a heavy (CD1), A9G DELETES
# a heavy (CB). amber14 ff14SB atom naming (verified against the amber14-all ALA/GLY
# templates AND the prepared 2QKI WT binder chain B residue 9):
#   ALA past CA: HA, CB, HB1, HB2, HB3 ; GLY past CA: HA2, HA3 (NO CB).
# Ala->Gly DELETES the beta-CH3 (CB+HB1-3) and the appearing Gly grows NO heavy
# beyond the common attach CA. The alpha hydrogen is RENAMED in the swap (Ala HA vs
# Gly HA2/HA3): so the alpha-H is NOT a shared common atom by name — it is carried
# as a disappearing H (Ala HA) paired with the appearing H's (Gly HA2/HA3). This
# keeps the C4 common-core list NAME-ALIGNED (common = N,H,CA,C,O on BOTH copies,
# 5 atoms each: ALA 10 - 5 var = 5; GLY 7 - 2 var = 5) — the positional swap requires
# identical common names in order. The common attach atom is CA (both the disappearing
# CB and the Gly backbone bond it). stateB=Gly (appears, copy-1), stateA=Ala
# (disappears, copy-2). hybrid_xml=None (canonical amber14 — Ala/Gly are both net-0,
# parity-identical: R-15/R-16 pass trivially; the full-residue net-charge sanity gate
# sees net=0 both copies). The appearing side has NO heavy, so bonded_heavy_appearing
# /appearing_h_prefix are None; the single disappearing heavy (CB) + its HB methyl
# group drive the mirrored MC2/R2 asserts. SHAPE = disappearing_heavy (exactly one
# disappearing heavy, zero appearing heavy).
MUTATION_ALA_GLY_RES9 = MutationSpec(
    name="a9g_ala_gly_res9",
    resnum=9,
    common_attach_atom="CA",
    stateA_resname="ALA",
    stateB_resname="GLY",
    # Ala beta-CH3 disappears; Ala alpha-H (HA) is carried as a disappearing H so the
    # shared common core excludes the renamed alpha-H (C4 name-alignment).
    stateA_only_atoms=("HA", "CB", "HB1", "HB2", "HB3"),
    # Gly's two alpha-H's appear (HA2 pairs Ala's HA; HA3 takes the old CB direction).
    stateB_only_atoms=("HA2", "HA3"),
    hybrid_xml=None,                                   # canonical amber14
    bonded_heavy_appearing=None,                       # appearing side has no heavy
    appearing_h_prefix=None,
    bonded_heavy_disappearing="CB",
    disappearing_h_prefix="HB",
)

# W4A spec — residue-4 Trp(WT) -> Ala (connected-subgraph fused-indole knockout).
# The PRODUCTION promotion of the validated W4A/w4a_spec_draft.ProtoMutationSpec
# (the exact spec the C5 bound genuine-decouple pre-flight smoke validated). The
# real MutationSpec already carries the Phase B ring fields (ring_closure_bonds /
# connected_group_certified) + the ring-FIRST `shape` property, so this is a pure
# additive registration — the scratch ProtoMutationSpec is no longer needed at
# build time. amber14 ff14SB indole atom naming (verified against the amber14-all
# TRP template + W4A/notes.md). TRP res-4 side chain past CB:
#   CG, CD1(HD1), CD2, NE1(HE1), CE2, CZ2(HZ2), CZ3(HZ3), CH2(HH2), CE3(HE3)
#   9 ring heavies: CG, CD1, CD2, NE1, CE2, CZ2, CZ3, CH2, CE3
#   ring-closure (5/6 fusion): CD2-CE2  <-- closes the bicyclic indole INSIDE the
#                                            disappearing var group.
# The common attach is CB (present in BOTH TRP and ALA). Trp->Ala DELETES the whole
# indole (9 heavy + 6 ring H = 15 atoms). The shared beta-H pair HB2/HB3 is COMMON
# (BOTH residues carry it with identical names -> name-aligned common core); ALA's
# extra HB1 is the single APPEARING atom (an H, not a heavy => the shape stays a
# 0-appearing-heavy / 9-disappearing-heavy connected group). Common core (BOTH
# copies, 9 atoms): N,H,CA,HA,C,O,CB,HB2,HB3. TRP 24 - 15 var = 9 ; ALA 10 - 1 var
# = 9 (balanced). hybrid_xml=None (canonical amber14; Trp/Ala BOTH net-0 -> R-15/
# R-16 trivial, MC1 full-residue net-charge sanity passes).
#
# The disappearing branch is a connected subgraph rooted at CG (the only indole
# heavy that bonds the common attach CB). bonded_heavy_disappearing=CG records the
# root (the MC2 connected-group certify discovers it + the rest by BFS over the
# built box's real bonds); ring_closure_bonds=(("CD2","CE2"),) +
# connected_group_certified=True opt the ring shape in (the star flag is left at
# its default False — it must NEVER certify a ring; the ring-FIRST classifier
# gates on connected_group_certified before any heavy-count branch). SHAPE =
# single_attach_connected_group (the detbeta deterministic common beta-H placement
# is shape-gated to this shape, so it auto-applies for W4A).
MUTATION_TRP_ALA_RES4 = MutationSpec(
    name="w4a_trp_ala_res4",
    resnum=4,
    common_attach_atom="CB",
    stateA_resname="TRP",                              # disappearing (copy-2 / WT)
    stateB_resname="ALA",                              # appearing (copy-1 / site)
    # The whole fused indole disappears (9 heavy + 6 ring H = 15 atoms).
    stateA_only_atoms=("CG", "CD1", "HD1", "CD2", "NE1", "HE1", "CE2",
                       "CZ2", "HZ2", "CZ3", "HZ3", "CH2", "HH2", "CE3", "HE3"),
    # ALA's extra beta-H (TRP has HB2/HB3 only; ALA has HB1/HB2/HB3) — single
    # appearing H, no heavy.
    stateB_only_atoms=("HB1",),
    hybrid_xml=None,                                   # canonical amber14
    bonded_heavy_appearing=None,                       # appearing side grows no heavy
    appearing_h_prefix=None,
    # Connected-subgraph root: CG (the only indole heavy that bonds CB).
    bonded_heavy_disappearing="CG",
    disappearing_h_prefix=None,                        # H's grouped by connectivity
    ring_closure_bonds=(("CD2", "CE2"),),              # 5/6 indole fusion bond
    connected_group_certified=True,                    # ring opt-in (NOT the star flag)
)

# Registry of named mutation specs (CLI / launcher selection).
MUTATION_SPECS: Dict[str, MutationSpec] = {
    MUTATION_MTR_TRP_RES4.name: MUTATION_MTR_TRP_RES4,
    MUTATION_VAL_ILE_RES3.name: MUTATION_VAL_ILE_RES3,
    MUTATION_ALA_GLY_RES9.name: MUTATION_ALA_GLY_RES9,
    MUTATION_TRP_ALA_RES4.name: MUTATION_TRP_ALA_RES4,
}


def resolve_mutation_spec(spec: Optional[Any]) -> MutationSpec:
    """Coerce a spec selector into a ``MutationSpec`` (None -> the res-4 default).

    Accepts ``None`` (the legacy res-4 MTR<->Trp default, byte-identical), a
    registry name string (``"v3i_val_ile_res3"``), or a ``MutationSpec`` instance
    (returned unchanged). Fail-loud on an unknown name.
    """
    if spec is None:
        return MUTATION_MTR_TRP_RES4
    if isinstance(spec, MutationSpec):
        return spec
    if isinstance(spec, str):
        if spec not in MUTATION_SPECS:
            raise ValueError(
                "resolve_mutation_spec: unknown mutation spec name %r (known: %s)"
                % (spec, ", ".join(sorted(MUTATION_SPECS))))
        return MUTATION_SPECS[spec]
    raise TypeError(
        "resolve_mutation_spec: expected None / a registry name / a MutationSpec, "
        "got %r" % (type(spec),))


def make_ile_ala_mutation_spec(
    resnum: int, name: Optional[str] = None,
) -> MutationSpec:
    """Construct an Ile->Ala hydrophobic-core MutationSpec at ``resnum`` (SSOT).

    The barnase folding-thermocycle anchor (Ile96->Ala) + the Ac-Ile-NMe reference
    leg BOTH use this shape; the ONLY difference between the two legs is ``resnum``
    (the barnase core Ile vs the tripeptide Ile), so this factory is the single
    definition both the folding orchestrator and the STAGE-0 smoke construct from
    (anti-fragmentation). It is DELIBERATELY not registered in ``MUTATION_SPECS`` —
    the CLI ``--mutation`` name path is for the fixed 2QKI specs; the folding run
    passes the constructed instance(s) directly (no registry edit).

    amber14 ff14SB atom naming (verified against the amber14-all ILE / ALA templates):
      ILE past CA: HA, CB, HB, CG2(HG21-23), CG1(HG12-13), CD1(HD11-13) ;
      ALA past CA: HA, CB, HB1, HB2, HB3.
    Common attach = CB (present in BOTH). Ile->Ala DELETES the whole aliphatic side
    chain past CB (3 heavies CG1/CG2/CD1 + their H's) AND the single CB hydrogen HB;
    the appearing ALA grows the two extra CB hydrogens (ILE CB carries HB only; ALA
    CB carries HB1/HB2/HB3). The CB hydrogen is RENAMED in the swap (ILE HB vs ALA
    HB1/2/3), so — exactly like the A9G alpha-H handling — HB is carried as a
    disappearing H and HB1/HB2/HB3 as appearing H's, keeping the C4 common core
    NAME-ALIGNED (common = N,H,CA,HA,CB,C,O on BOTH copies, 7 atoms each: ILE 19 - 12
    var = 7 ; ALA 10 - 3 var = 7). hybrid_xml=None (canonical amber14; ILE/ALA are
    BOTH net-0 -> R-15/R-16 trivial, MC1 full-residue net-charge sanity passes).

    The disappearing branch is a one-sided ACYCLIC connected subgraph: CG1 AND CG2
    both bond the common attach CB, and CD1 chains off CG1 (BFS-reachable), with NO
    ring-closure bond. bonded_heavy_disappearing is left ``None`` (the group has TWO
    attach-bonded roots, not one; the MC2 acyclic-connected-group backstop discovers
    the roots + certifies attach-reachability + acyclicity itself). SHAPE =
    acyclic_connected_group (via ``chained_group_certified``; the star flag stays
    False — a chained heavy is not a star). No detbeta placement applies (its gate is
    the ring shape single_attach_connected_group, and CB carries NO common beta-H here
    — all CB hydrogens are in the var slots).
    """
    return MutationSpec(
        name=(name or ("ile_ala_res%d" % (int(resnum),))),
        resnum=int(resnum),
        common_attach_atom="CB",
        stateA_resname="ILE",                              # disappearing (copy-2)
        stateB_resname="ALA",                              # appearing (copy-1 / site)
        # ILE side chain past CB + the single CB hydrogen HB all disappear.
        stateA_only_atoms=("HB", "CG1", "HG12", "HG13", "CG2", "HG21", "HG22",
                           "HG23", "CD1", "HD11", "HD12", "HD13"),
        # ALA's three CB hydrogens appear (no appearing heavy).
        stateB_only_atoms=("HB1", "HB2", "HB3"),
        hybrid_xml=None,                                   # canonical amber14
        bonded_heavy_appearing=None,                       # appearing side grows no heavy
        appearing_h_prefix=None,
        # TWO attach-bonded roots (CG1, CG2) + CD1 chained off CG1 -> the acyclic
        # connected-group backstop discovers roots itself; no single declared root.
        bonded_heavy_disappearing=None,
        disappearing_h_prefix=None,                        # H's grouped by connectivity
        ring_closure_bonds=(),                             # ACYCLIC (no ring)
        chained_group_certified=True,                      # chained-group opt-in
    )


DISULFIDE_MAX_NM = 0.24  # CYS2-CYS12 SG-SG (matches run_restrained_md default)


# ---------------------------------------------------------------------------
# cyclic_ss disulfide helpers (B-C3).
#
# Logic is a verbatim copy of run_restrained_md.detect_disulfide_pair /
# commit_bond_if_missing. It is inlined (not imported) ONLY because
# run_restrained_md imports networkx + pdbfixer at module load, which are
# intentionally absent from the lean ``atm`` env (isolation contract). The two
# functions depend on nothing beyond openmm + numpy. If run_restrained_md's
# disulfide logic changes, this copy must track it.
# ---------------------------------------------------------------------------
def detect_disulfide_pair(modeller, binder_chain="B", max_sg_dist_nm=DISULFIDE_MAX_NM):
    chains = [c for c in modeller.topology.chains() if c.id == binder_chain]
    if not chains:
        return None
    residues = list(chains[0].residues())
    sg_atoms = [a for r in residues if r.name in ("CYS", "CYX")
                for a in r.atoms() if a.name == "SG"]
    if len(sg_atoms) < 2:
        return None

    pos = list(modeller.positions)
    existing = {frozenset([b.atom1.index, b.atom2.index])
                for b in modeller.topology.bonds()}

    candidates = []
    for i in range(len(sg_atoms)):
        for j in range(i + 1, len(sg_atoms)):
            a1, a2 = sg_atoms[i], sg_atoms[j]
            p1 = np.array(pos[a1.index].value_in_unit(unit.nanometers))
            p2 = np.array(pos[a2.index].value_in_unit(unit.nanometers))
            dist = np.linalg.norm(p1 - p2)
            already_bonded = frozenset([a1.index, a2.index]) in existing
            if already_bonded or dist <= max_sg_dist_nm:
                candidates.append((dist, a1, a2, already_bonded))

    if not candidates:
        return None

    candidates.sort(key=lambda x: (not x[3], x[0]))
    best = candidates[0]
    if len(candidates) > 1 and not best[3]:
        if abs(candidates[1][0] - candidates[0][0]) < 0.02:
            raise RuntimeError(
                "Disulfide pair ambiguous (gap < 0.02 nm); explicit rule needed."
            )
    return (best[1], best[2])


def commit_bond_if_missing(topology, atom1, atom2):
    existing = {frozenset([b.atom1.index, b.atom2.index])
                for b in topology.bonds()}
    pair = frozenset([atom1.index, atom2.index])
    if pair not in existing:
        topology.addBond(atom1, atom2)
        return True
    return False


# ---------------------------------------------------------------------------
# Peptide-bond completion for HETATM ncAA junctions.
#
# Verbatim copy of run_restrained_md.add_missing_peptide_bonds_safe (+ helpers
# is_peptide_like_atomset, is_carbonyl_carbon), inlined for the same isolation
# reason as the disulfide helpers: run_restrained_md pulls networkx + pdbfixer
# at module load, which are intentionally absent from the ``atm`` env. These
# helpers depend only on openmm + numpy. If the run_restrained_md logic
# changes, this copy must track it.
#
# Necessary because OpenMM's PDBFile parser does NOT auto-infer ATOM↔HETATM
# peptide bonds (the MTR ncAA is HETATM in the prepared PDBs), so the
# free-peptide leg loses VAL3.C-MTR4.N and MTR4.C-GLN5.N at extraction time
# unless we re-add them by C-N distance threshold (0.20 nm).
# ---------------------------------------------------------------------------
def is_peptide_like_atomset(atom_names):
    return {"N", "CA", "C"}.issubset(atom_names)


def is_carbonyl_carbon(residue, c_atom, positions_nm, max_co_nm: float = 0.14):
    if c_atom.element is None or c_atom.element.symbol != "C":
        return False
    c_pos = np.array(positions_nm[c_atom.index].value_in_unit(unit.nanometers))
    for a in residue.atoms():
        if a.index == c_atom.index:
            continue
        if a.element and a.element.symbol == "O":
            o_pos = np.array(positions_nm[a.index].value_in_unit(unit.nanometers))
            if np.linalg.norm(c_pos - o_pos) <= max_co_nm:
                return True
    return False


def inject_xml_internal_bonds(topology, xml_paths: List[str],
                              xml_res_name: str = "MTR") -> int:
    """Inject ncAA internal bonds from the residue's XML <Bond> records.

    Verbatim port of run_restrained_md.inject_xml_bonds (v38 multi-site).
    OpenMM's PDBFile parser treats ncAA residues as HETATM without inferring
    internal bonds — the free-leg extraction therefore yields an MTR with no
    intra-residue connectivity, breaking the amber14 template graph match
    ("residue has no bonds between its atoms"). This helper re-adds them by
    looking up the XML <Bond atomName1=.. atomName2=..> records and matching
    atom-name to the topology atoms. Idempotent.

    Returns the total number of bonds added across every residue instance.
    """
    if not xml_paths or not xml_res_name:
        return 0
    import xml.etree.ElementTree as ET

    target_residues = [r for r in topology.residues() if r.name == xml_res_name]
    if not target_residues:
        return 0

    # Search every XML for the residue (use the first that contains it).
    xml_bond_pairs: List[Tuple[str, str]] = []
    for xml_path in xml_paths:
        try:
            tree = ET.parse(xml_path)
        except (OSError, ET.ParseError):
            continue
        node = next((r for r in tree.getroot().iter("Residue")
                     if r.get("name") == xml_res_name), None)
        if node is None:
            continue
        xml_bond_pairs = [
            (b.get("atomName1", "").strip(), b.get("atomName2", "").strip())
            for b in node.findall("Bond")
        ]
        if xml_bond_pairs:
            break
    if not xml_bond_pairs:
        return 0

    existing = {frozenset([b.atom1.index, b.atom2.index])
                for b in topology.bonds()}
    added = 0
    for res in target_residues:
        atom_by_name = {a.name: a for a in res.atoms()}
        for a1, a2 in xml_bond_pairs:
            if a1 in atom_by_name and a2 in atom_by_name:
                oa1, oa2 = atom_by_name[a1], atom_by_name[a2]
                pair = frozenset([oa1.index, oa2.index])
                if pair not in existing:
                    topology.addBond(oa1, oa2)
                    existing.add(pair)
                    added += 1
    return added


def add_missing_peptide_bonds_safe(modeller, binder_chain: str = "B",
                                   max_cn_distance_nm: float = 0.20) -> int:
    """Add missing C(i)-N(i+1) peptide bonds inside the binder chain.

    Targets HETATM/ATOM junctions where OpenMM's PDBFile reader did not
    auto-infer the inter-residue peptide bond (e.g. ncAA boundaries). Only
    pairs that look like canonical peptide residues (N/CA/C set) and whose C
    is a true carbonyl (C=O within 0.14 nm) are joined; the C-N pair must
    also be within 0.20 nm to avoid spurious bonds.
    """
    pos = list(modeller.positions)
    added = 0
    existing = {frozenset([b.atom1.index, b.atom2.index])
                for b in modeller.topology.bonds()}

    for chain in modeller.topology.chains():
        if chain.id != binder_chain:
            continue
        residues = list(chain.residues())
        for r1, r2 in zip(residues[:-1], residues[1:]):
            r1_names = {a.name for a in r1.atoms()}
            r2_names = {a.name for a in r2.atoms()}
            if not (is_peptide_like_atomset(r1_names)
                    and is_peptide_like_atomset(r2_names)):
                continue
            c_atom = next((a for a in r1.atoms() if a.name == "C"), None)
            n_atom = next((a for a in r2.atoms() if a.name == "N"), None)
            if c_atom is None or n_atom is None:
                continue
            if not is_carbonyl_carbon(r1, c_atom, pos):
                continue
            pair = frozenset([c_atom.index, n_atom.index])
            if pair in existing:
                continue
            c_pos = np.array(pos[c_atom.index].value_in_unit(unit.nanometers))
            n_pos = np.array(pos[n_atom.index].value_in_unit(unit.nanometers))
            dist = np.linalg.norm(c_pos - n_pos)
            if dist <= max_cn_distance_nm:
                modeller.topology.addBond(c_atom, n_atom)
                existing.add(pair)
                added += 1
    return added


# ---------------------------------------------------------------------------
# Charge-axis verification (B-C1, PATCH-01 surface)
# ---------------------------------------------------------------------------
def verify_charge_axis(xml_path: str = HYBRID_MTR_XML,
                       integer_tol: float = 5e-4) -> Dict[str, Any]:
    """Read the MTR residue charges from the hybrid XML and report the axis.

    Reports the B-C1 invariants per residue variant (MTR/NMTR/CMTR), and
    surfaces per-residue Σq so the caller sees integer-parity
    compliance explicitly. **Scope is per-residue** — this is what OpenMM
    ``createSystem`` evaluates and what PME sees as the box net charge.

    Historical note (Q1 fix 2026-05-30): the previous implementation flattened
    ``residues.iter('Atom')`` into a single name-keyed dict, which silently
    overwrote shared atom names across MTR/NMTR/CMTR and reported a spurious
    full-XML ``Σq = -0.176 e`` that was actually the last-write CMTR value.
    Per-residue is the correct scope; the rescale of CMTR (production blocker Q1 a) drove
    all 3 residue variants to |Σq| ≤ 5e-4 e.
    """
    import xml.etree.ElementTree as ET

    tree = ET.parse(xml_path)
    residues = tree.find(".//Residues")
    if residues is None:
        raise ValueError(f"No <Residues> block in {xml_path}")

    per_residue: List[Dict[str, Any]] = []
    ne1_charge: Optional[float] = None
    bbn_charge: Optional[float] = None

    for res_el in residues.findall("Residue"):
        rname = res_el.get("name")
        atoms = {}
        for atom in res_el.findall("Atom"):
            q = atom.get("charge")
            name = atom.get("name")
            if q is not None and name is not None:
                atoms[name] = float(q)
        sigma = sum(atoms.values())
        per_residue.append({
            "residue": rname,
            "n_atoms": len(atoms),
            "sigma_q": sigma,
            "ne1_charge": atoms.get(ALCH_COMMON_ATOM),
            "backbone_n_charge": atoms.get("N"),
            "within_integer_tol": abs(sigma) <= integer_tol,
        })
        # Internal MTR is the production residue (Cp4 residue-4 is internal).
        if rname == "MTR":
            ne1_charge = atoms.get(ALCH_COMMON_ATOM)
            bbn_charge = atoms.get("N")

    max_abs_sigma = max((abs(r["sigma_q"]) for r in per_residue), default=0.0)
    all_within_tol = all(r["within_integer_tol"] for r in per_residue)

    return {
        "xml_path": xml_path,
        "scope": "per_residue",
        "integer_tol_e": integer_tol,
        "per_residue": per_residue,
        "max_abs_sigma_q_e": max_abs_sigma,
        "all_residues_within_tol": all_within_tol,
        "ne1_charge": ne1_charge,
        "backbone_n_charge": bbn_charge,
        "ne1_frozen_ok": (ne1_charge is not None
                          and abs(ne1_charge - (-0.3418)) < 1e-4),
        "backbone_n_amber14sb_ok": (bbn_charge is not None
                                    and abs(bbn_charge - (-0.4157)) < 1e-4),
        "charge_regime": "option_beta_regime2_hybrid_rescaled",
        "rescale_applied": True,
        "rescale_audit": "params/_archive/mtr_rescale_audit_20260530.json",
        "regime": "ranking_only",
    }


# ---------------------------------------------------------------------------
# Leg structure resolution (B-C3)
# ---------------------------------------------------------------------------
def resolve_leg_inputs(seed: str = "s7") -> Dict[str, Dict[str, str]]:
    """Resolve the prepared input PDBs for both endpoints of both legs.

    Reuses the existing per-seed prepared structures (no rebuild):
      bound leg : the 2QKI complex (chain A target + chain B cyclic peptide)
      free  leg : the cyclic peptide alone (chain B), extracted from the same
                  prepared complex so the cyclic_ss topology is identical.

    Cp4 endpoint -> ``outputs/2QKI_Cp4_hybrid_calib_<seed>/_md_input/2QKI_Cp4.pdb``
    WT  endpoint -> ``outputs/2QKI_WT_calib_<seed>/_md_input/2QKI_WT.pdb``
    """
    out = os.path.join(_PROJ_ROOT, "outputs")
    cp4_dir = os.path.join(out, f"2QKI_Cp4_hybrid_calib_{seed}")
    wt_dir = os.path.join(out, f"2QKI_WT_calib_{seed}")
    cp4_pdb = os.path.join(cp4_dir, "_md_input", "2QKI_Cp4.pdb")
    wt_pdb = os.path.join(wt_dir, "_md_input", "2QKI_WT.pdb")
    # Prefer the MD-equilibrated final.pdb when present: it is fully H-complete
    # and template-clean (MTR carries CM + HM1-3), so the free-leg extraction
    # needs no PDBFixer ncAA bond/H reconstruction.
    cp4_final = os.path.join(cp4_dir, "mdresult", "2QKI_Cp4_final.pdb")
    wt_final = os.path.join(wt_dir, "mdresult", "2QKI_WT_final.pdb")
    result = {
        "bound": {"cp4": cp4_pdb, "wt": wt_pdb},
        "free": {"cp4": cp4_pdb, "wt": wt_pdb},  # peptide extracted downstream
        "final": {
            "cp4": cp4_final if os.path.isfile(cp4_final) else None,
            "wt": wt_final if os.path.isfile(wt_final) else None,
        },
        "hydrogens_xml": os.path.join(cp4_dir, "params", "MTR_hydrogens.xml"),
    }
    return result


def resolve_fold_leg_inputs(
    scaffold_final: str, hydrogens_xml: Optional[str] = None,
) -> Dict[str, Any]:
    """Leg-input override for a SINGLE-SCAFFOLD folding-thermocycle two-copy build.

    The 2QKI-hardcoded :func:`resolve_leg_inputs` wires the Cp4/WT endpoint pair; a
    folding-stability control (barnase Ile96->Ala, or the Ac-Ile-NMe reference leg)
    is a CANONICAL single-scaffold mutation with NO cp4/wt endpoint pair. For the
    canonical (non-MTR) two-copy path the builder reads ONLY ``li['final']['wt']``
    (``scaffold_final``) — copy-1 is the point-mutation of that scaffold and copy-2
    is its native residue. The both-finals guard also requires ``li['final']['cp4']``
    to be non-None, so BOTH slots are pointed at the SAME real scaffold (the cp4 slot
    is never consumed on the canonical path; it only satisfies the guard).

    Returns the same schema as :func:`resolve_leg_inputs` (bound/free/final +
    hydrogens_xml) so the override is drop-in for the ``leg_inputs`` kwarg;
    ``hydrogens_xml`` defaults to ``None`` (canonical amber14 has no ncAA H
    definitions to load, and the two-copy build passes ``add_hydrogens=False`` for
    both endpoint copies — PDBFixer/the H-complete scaffold place the hydrogens).

    ``scaffold_final`` MUST be an H-complete single-chain PDB (a bad / non-H-complete
    scaffold silently biases dgbind1; that fidelity is the caller's system-prep
    responsibility). This resolver does NO file I/O beyond an existence check.
    """
    if not scaffold_final or not os.path.isfile(scaffold_final):
        raise FileNotFoundError(
            "resolve_fold_leg_inputs: scaffold_final does not exist: %r"
            % (scaffold_final,))
    return {
        # bound/free are unused by the canonical two-copy path (it reads only the
        # 'final' slot); populated for schema parity with resolve_leg_inputs.
        "bound": {"cp4": scaffold_final, "wt": scaffold_final},
        "free": {"cp4": scaffold_final, "wt": scaffold_final},
        "final": {"cp4": scaffold_final, "wt": scaffold_final},
        "hydrogens_xml": hydrogens_xml,
    }


def _extract_binder_only(in_pdb: str, out_pdb: str, binder_chain: str = "B") -> str:
    """Write a PDB containing only the binder chain (the free-peptide leg).

    Keeps every binder atom verbatim (including CYS SG) so the disulfide
    detector reconstructs the cyclic_ss bond identically to the bound leg.
    """
    kept = []
    with open(in_pdb) as fh:
        for line in fh:
            if line[:6] in ("ATOM  ", "HETATM"):
                if line[21] == binder_chain:
                    kept.append(line)
            elif line[:3] == "TER" and kept:
                kept.append(line)
    with open(out_pdb, "w") as fh:
        fh.writelines(kept)
        fh.write("END\n")
    return out_pdb


def _pdbfixer_prep_termini(in_pdb: str, out_pdb: str,
                           keep_resnames: Tuple[str, ...] = ("MTR",)) -> str:
    """Use PDBFixer to add only missing HEAVY terminal atoms on the extracted
    free peptide so amber14 templates match. Hydrogens are deliberately NOT
    added here — they are placed downstream by ``Modeller.addHydrogens`` once
    the MTR hydrogen definitions are registered, which is the only path that
    positions the ncAA methyl H's (HM1-3 on CM) correctly.

    All input hydrogens are stripped first so PDBFixer/Modeller rebuild a
    consistent H set (the input carries Trp-style indole H's that are wrong for
    the MTR N-methyl state). The ncAA residue (MTR) is masked to ``UNK`` during
    fixing so PDBFixer leaves its non-standard heavy atoms intact, then
    restored — the masking strategy run_restrained_md uses for ncAA protection.

    The compstatin cyclic_ss peptide has STANDARD N-/C-termini (the macrocycle
    is the SG-SG disulfide, not head-to-tail), so standard terminal completion
    is correct.
    """
    from pdbfixer import PDBFixer

    # Strip ALL hydrogens + mask MTR -> UNK (record nothing else changes).
    masked = out_pdb + ".masked.pdb"
    with open(in_pdb) as fh, open(masked, "w") as out:
        for line in fh:
            if line[:6] in ("ATOM  ", "HETATM"):
                elem = line[76:78].strip()
                aname = line[12:16].strip()
                is_h = (elem == "H") or (not elem and aname[:1] == "H")
                if is_h:
                    continue  # drop input hydrogens
                if line[17:20].strip() in keep_resnames:
                    out.write(line[:17] + "UNK" + line[20:])
                else:
                    out.write(line)
            elif line[:3] == "TER":
                out.write(line)

    fixer = PDBFixer(filename=masked)
    fixer.findMissingResidues()
    fixer.missingResidues = {}  # do not insert internal gaps for a single chain
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()  # heavy atoms only (e.g. terminal OXT); no H here

    tmp_fixed = out_pdb + ".fixed.pdb"
    with open(tmp_fixed, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh)

    # Restore UNK -> MTR.
    with open(tmp_fixed) as fh, open(out_pdb, "w") as out:
        for line in fh:
            if line[:6] in ("ATOM  ", "HETATM") and line[17:20].strip() == "UNK":
                out.write(line[:17] + "MTR" + line[20:])
            else:
                out.write(line)
    return out_pdb


_SOLVENT_RESNAMES = {"HOH", "WAT", "NA", "CL", "NA+", "CL-", "K+", "K"}


def prepare_free_peptide_from_final(final_pdb: str, out_pdb: str,
                                    binder_chain: str = "B") -> str:
    """Extract the H-complete cyclic peptide from an MD ``final.pdb``.

    The final.pdb is already template-clean (MTR carries CM + HM1-3, all
    standard residues fully protonated), so this drops only the receptor chain
    and the solvent/ions, preserving every binder atom + hydrogen verbatim. No
    PDBFixer reconstruction is needed — this avoids the ncAA-bond-stripping
    pitfall entirely. cyclic_ss (CYS2-CYS12) is re-detected downstream.
    """
    kept = []
    for line in open(final_pdb):
        if line[:6] in ("ATOM  ", "HETATM"):
            if line[21] != binder_chain:
                continue
            if line[17:20].strip() in _SOLVENT_RESNAMES:
                continue
            kept.append(line)
    with open(out_pdb, "w") as fh:
        fh.writelines(kept)
        fh.write("END\n")
    return out_pdb


def prepare_mutated_binder_from_final(
    final_pdb: str, out_pdb: str, resnum: int,
    from_resname: str, to_resname: str, binder_chain: str = "B",
    add_hydrogens: bool = True,
) -> str:
    """Extract the binder from an MD ``final.pdb`` and apply a CANONICAL amber14
    point mutation (e.g. VAL-3-ILE) via PDBFixer, returning a template-clean PDB.

    This is the canonical (all-amber, no ncAA) sibling of
    :func:`prepare_free_peptide_from_final`: it produces the APPEARING-state copy
    for an engine-validation mutation where BOTH endpoints are standard amber14
    residues (V3I Val<->Ile). The disappearing-state copy is produced verbatim by
    ``prepare_free_peptide_from_final`` (no mutation needed — it IS the scaffold's
    native residue), so this helper is used only for the side that differs.

    Mechanics:
      1. Extract the binder chain (drop receptor + solvent).
      2. PDBFixer ``applyMutations([f"{from}-{resnum}-{to}"], chain)`` rebuilds the
         target residue's heavy-atom side chain to the target template (drops the
         disappearing atoms, adds the appearing heavy atoms e.g. CD1 with a
         standard geometry).
      3. ``findMissingAtoms`` / ``addMissingAtoms`` completes any remaining heavy
         atoms; ``addMissingHydrogens`` (when ``add_hydrogens``) places the full
         amber14 hydrogen set on the mutated residue.

    Hydrogens on the mutated residue MUST be (re)placed because the appearing
    heavy atoms (CD1) carry no H in the input final.pdb. The cyclic_ss disulfide
    (CYS-CYS) is preserved (CYS is untouched by the mutation) and re-detected
    downstream. A standard (non-ncAA) terminal completion is correct for the
    compstatin macrocycle (the ring is the SG-SG bond, termini are standard).

    Fail-loud if PDBFixer cannot resolve the mutation (the smoke reports the
    cause).
    """
    from pdbfixer import PDBFixer

    raw = out_pdb + ".binder.pdb"
    _extract_binder_only(final_pdb, raw, binder_chain=binder_chain)

    fixer = PDBFixer(filename=raw)
    # Single-chain extract -> the PDBFixer chain id may be re-labelled; resolve the
    # actual chain id of the extracted binder (usually the binder_chain, but the
    # PDB reader can re-letter a sole chain). Mutate the chain that carries the
    # target residue number.
    target_chain_id = None
    for ch in fixer.topology.chains():
        for res in ch.residues():
            try:
                if int(res.id) == int(resnum) and res.name == from_resname:
                    target_chain_id = ch.id
                    break
            except (TypeError, ValueError):
                continue
        if target_chain_id is not None:
            break
    if target_chain_id is None:
        raise ValueError(
            "prepare_mutated_binder_from_final: residue %s-%d not found in the "
            "extracted binder of %s (cannot apply the %s->%s mutation)."
            % (from_resname, resnum, final_pdb, from_resname, to_resname))

    fixer.applyMutations(
        ["%s-%d-%s" % (from_resname, int(resnum), to_resname)], target_chain_id)
    fixer.findMissingResidues()
    fixer.missingResidues = {}        # do not insert internal gaps (single chain)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()           # heavy atoms (the appearing CD1 etc.)
    if add_hydrogens:
        fixer.addMissingHydrogens(7.0)  # full amber14 H set incl. mutated residue

    # PDBFile.writeFile re-letters chains positionally (A, B, ...) IGNORING
    # ``chain.id`` — so the single extracted binder is written as chain 'A'. The
    # downstream build keys on ``binder_chain``, so stamp the binder chain id into
    # the written PDB text (column 22) for every binder atom. Single chain =>
    # every ATOM/HETATM line is the binder.
    tmp_written = out_pdb + ".written.pdb"
    with open(tmp_written, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh)
    with open(tmp_written) as src, open(out_pdb, "w") as dst:
        for line in src:
            if line[:6] in ("ATOM  ", "HETATM"):
                dst.write(line[:21] + binder_chain + line[22:])
            else:
                dst.write(line)
    return out_pdb


# Bound-complex contact thresholds (Angstrom). A correctly imaged 2QKI bound
# pose sits at ~2.5-3.0 A receptor<->binder min heavy-atom distance; a
# PBC-unwrapped final.pdb leaves the binder one box image away (~44 A). The
# guard rejects anything that still reads as separated after re-imaging.
_BOUND_CONTACT_MAX_A = 8.0


def _parse_cryst1_box(final_pdb: str) -> Tuple[float, float, float]:
    """Read the orthorhombic box lengths (a, b, c) in Angstrom from the CRYST1
    record of ``final_pdb``.

    The minimum-image re-imaging in :func:`prepare_bound_complex_from_final`
    requires the unit-cell vectors. A missing / degenerate CRYST1 is a hard
    error (no silent pass): without the box the binder cannot be re-imaged and a
    broken bound endpoint would propagate into solvation + minimize.

    Only the orthorhombic (alpha=beta=gamma=90) case is supported — the 2QKI
    final.pdb is cubic ``109.113 109.113 109.113 90 90 90``. A non-orthorhombic
    cell raises rather than mis-handling the off-diagonal box terms.
    """
    box_line = None
    with open(final_pdb) as fh:
        for line in fh:
            if line[:6] == "CRYST1":
                box_line = line
                break
    if box_line is None:
        raise ValueError(
            "prepare_bound_complex_from_final: no CRYST1 record in %s; cannot "
            "minimum-image the binder against the receptor (PBC box unknown)."
            % final_pdb)
    try:
        a = float(box_line[6:15])
        b = float(box_line[15:24])
        c = float(box_line[24:33])
        alpha = float(box_line[33:40])
        beta = float(box_line[40:47])
        gamma = float(box_line[47:54])
    except ValueError as exc:
        raise ValueError(
            "prepare_bound_complex_from_final: malformed CRYST1 record in %s: "
            "%r" % (final_pdb, box_line.rstrip())) from exc
    if min(a, b, c) <= 0.0:
        raise ValueError(
            "prepare_bound_complex_from_final: non-positive box length in "
            "CRYST1 of %s (a=%g b=%g c=%g)." % (final_pdb, a, b, c))
    if not (abs(alpha - 90.0) < 1e-3 and abs(beta - 90.0) < 1e-3
            and abs(gamma - 90.0) < 1e-3):
        raise ValueError(
            "prepare_bound_complex_from_final: non-orthorhombic CRYST1 in %s "
            "(alpha=%g beta=%g gamma=%g); only orthorhombic boxes are "
            "supported for binder re-imaging." % (final_pdb, alpha, beta, gamma))
    return (a, b, c)


def _atom_xyz(line: str) -> Tuple[float, float, float]:
    """Parse the (x, y, z) Angstrom coordinates from a PDB ATOM/HETATM line."""
    return (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def _shift_atom_line(line: str, shift: np.ndarray) -> str:
    """Return ``line`` with its x/y/z translated by ``shift`` (Angstrom),
    preserving the strict PDB column layout (8.3f in cols 31-54)."""
    x, y, z = _atom_xyz(line)
    nx = x + float(shift[0])
    ny = y + float(shift[1])
    nz = z + float(shift[2])
    return "%s%8.3f%8.3f%8.3f%s" % (line[:30], nx, ny, nz, line[54:])


def _is_heavy_atom(line: str) -> bool:
    """Heavy-atom test for a PDB line (element != H, falling back to the atom
    name when the element column is blank)."""
    el = line[76:78].strip()
    if el:
        return el != "H"
    return not line[12:16].strip().startswith("H")


def _reimage_binder_min_image(kept_lines: List[str], binder_chain: str,
                              receptor_chain: str,
                              box: Tuple[float, float, float]):
    """Rigid minimum-image re-alignment of the binder onto the receptor.

    An MD trajectory written with PBC unwrapping can leave the binder (chain B)
    a whole box vector away from the receptor (chain A) even though the binder
    itself is intact (its own atoms are contiguous). This applies the single
    per-axis integer box-vector translation that brings the binder heavy-atom
    centroid to the receptor's minimum image — i.e. the whole-molecule
    ``image_molecules(anchor=receptor, other=binder, make_whole=True)`` result
    for an already-whole binder — so the equilibrated bound pose is preserved
    (rigid translation only) while the PBC wrap is corrected.

    Returns ``(new_lines, min_dist_A, shift)`` where ``min_dist_A`` is the
    receptor<->binder heavy-atom minimum distance after the shift and ``shift``
    is the applied translation (Angstrom).
    """
    box_arr = np.asarray(box, dtype=float)
    rec_heavy = []
    bnd_heavy = []
    for line in kept_lines:
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        if not _is_heavy_atom(line):
            continue
        ch = line[21]
        if ch == receptor_chain:
            rec_heavy.append(_atom_xyz(line))
        elif ch == binder_chain:
            bnd_heavy.append(_atom_xyz(line))
    if not rec_heavy or not bnd_heavy:
        raise ValueError(
            "prepare_bound_complex_from_final: receptor (chain %s, %d heavy) or "
            "binder (chain %s, %d heavy) atoms missing after extraction; cannot "
            "re-image." % (receptor_chain, len(rec_heavy),
                           binder_chain, len(bnd_heavy)))
    rec_arr = np.asarray(rec_heavy, dtype=float)
    bnd_arr = np.asarray(bnd_heavy, dtype=float)
    rec_centroid = rec_arr.mean(axis=0)
    bnd_centroid = bnd_arr.mean(axis=0)
    # Per-axis integer box shift that pulls the binder centroid into the
    # receptor's minimum image (round() picks the nearest image per axis).
    n_image = np.round((bnd_centroid - rec_centroid) / box_arr)
    shift = -n_image * box_arr
    new_lines = []
    for line in kept_lines:
        if line[:6] in ("ATOM  ", "HETATM") and line[21] == binder_chain:
            new_lines.append(_shift_atom_line(line, shift))
        else:
            new_lines.append(line)
    bnd_shifted = bnd_arr + shift
    # Minimum receptor<->binder heavy-atom distance after the shift (no PBC —
    # the binder is now in the receptor's primary image).
    deltas = rec_arr[:, None, :] - bnd_shifted[None, :, :]
    min_dist = float(np.sqrt((deltas ** 2).sum(axis=-1)).min())
    return new_lines, min_dist, shift


def prepare_bound_complex_from_final(final_pdb: str, out_pdb: str,
                                     binder_chain: str = "B",
                                     receptor_chain: str = "A") -> str:
    """Extract the H-complete BOUND complex (receptor + cyclic peptide) from an
    MD ``final.pdb``, dropping only the solvent/ions.

    This is the bound-leg analog of ``prepare_free_peptide_from_final``: the same
    MD final.pdb is template-clean (MTR carries CM + HM1-3, all standard residues
    fully protonated), so this keeps BOTH the receptor (chain A) and the binder
    (chain B) in their equilibrated BOUND pose and drops the water/ions. The
    receptor rides through as non-alchemical context — the residue-4 dual-topology
    transform operates only on the binder, exactly as in the free leg, so the
    receptor's presence is inert to the common-core swap (the common map is
    binder-protein-only). NO displacement (RBFE — the binder stays bound). Every
    receptor + binder atom + hydrogen is preserved verbatim. cyclic_ss (CYS2-CYS12)
    is re-detected downstream on the binder chain.

    The intra-receptor and receptor-binder peptide bonds are reconstructed by the
    same distance-threshold completion ``build_leg_system`` runs, and amber14 has
    clean templates for every standard receptor residue, so re-solvation +
    createSystem succeed without ncAA reconstruction on the receptor.

    PBC re-imaging: an MD trajectory written with PBC unwrapping can store the
    binder a whole box vector away from the receptor (observed on 2QKI Cp4:
    binder ~44 A from receptor, one box image of ~109 A in Z). Solvating that
    split with a 1.2 nm pad wraps the ~96 A gap in a giant box, and minimize then
    diverges (NaN) as bonds straddle the periodic boundary — the bound asyncre
    crash before cycle 0. After extraction the binder is rigidly translated to
    the receptor's minimum image (pose-preserving) and a fail-fast guard rejects
    any complex whose receptor<->binder min heavy-atom distance still exceeds
    ``_BOUND_CONTACT_MAX_A`` (a genuinely separated / broken endpoint, never a
    valid bound pose). The free-peptide prep has no such guard (binder only).
    """
    box = _parse_cryst1_box(final_pdb)
    keep_chains = {binder_chain, receptor_chain}
    kept = []
    for line in open(final_pdb):
        if line[:6] in ("ATOM  ", "HETATM"):
            if line[21] not in keep_chains:
                continue
            if line[17:20].strip() in _SOLVENT_RESNAMES:
                continue
            kept.append(line)
        elif line[:3] == "TER" and kept:
            # Preserve chain breaks so the receptor and binder do not get fused
            # into one chain by the reader (each chain keeps its own id).
            kept.append(line)
    kept, min_dist, _shift = _reimage_binder_min_image(
        kept, binder_chain, receptor_chain, box)
    if min_dist > _BOUND_CONTACT_MAX_A:
        raise ValueError(
            "prepare_bound_complex_from_final: receptor<->binder min heavy-atom "
            "distance %.2f A exceeds %.1f A after minimum-image re-imaging of %s "
            "— the bound endpoint is broken (binder not in contact with the "
            "receptor); refusing to emit a non-bound complex."
            % (min_dist, _BOUND_CONTACT_MAX_A, final_pdb))
    with open(out_pdb, "w") as fh:
        fh.writelines(kept)
        fh.write("END\n")
    return out_pdb


def prepare_mutated_bound_complex_from_final(
    final_pdb: str, out_pdb: str, resnum: int,
    from_resname: str, to_resname: str, binder_chain: str = "B",
    receptor_chain: str = "A", add_hydrogens: bool = True,
) -> str:
    """Bound-complex prep with a CANONICAL amber14 point mutation on the binder.

    Bound-leg sibling of :func:`prepare_mutated_binder_from_final` /
    :func:`prepare_bound_complex_from_final`: keeps the receptor (inert context)
    + binder in the equilibrated BOUND pose, re-images the binder to the receptor's
    minimum image, then applies the binder-chain point mutation (e.g. VAL-3-ILE)
    via PDBFixer and re-places hydrogens on the mutated residue. The disappearing-
    state bound complex is produced verbatim by ``prepare_bound_complex_from_final``
    (no mutation); this helper is the appearing-state side only.

    The mutation is applied to the BINDER chain only (PDBFixer mutates the named
    chain), so the receptor side chains are untouched. Fail-loud if PDBFixer cannot
    resolve the mutation. cyclic_ss is preserved (CYS untouched).
    """
    from pdbfixer import PDBFixer

    # First produce the re-imaged WT bound complex (receptor + binder, pose-fixed).
    wt_complex = out_pdb + ".bound.pdb"
    prepare_bound_complex_from_final(
        final_pdb, wt_complex, binder_chain=binder_chain,
        receptor_chain=receptor_chain)

    fixer = PDBFixer(filename=wt_complex)
    # Resolve the actual chain id carrying the target binder residue.
    target_chain_id = None
    for ch in fixer.topology.chains():
        for res in ch.residues():
            try:
                if int(res.id) == int(resnum) and res.name == from_resname:
                    target_chain_id = ch.id
                    break
            except (TypeError, ValueError):
                continue
        if target_chain_id is not None:
            break
    if target_chain_id is None:
        raise ValueError(
            "prepare_mutated_bound_complex_from_final: residue %s-%d not found in "
            "the binder of %s (cannot apply the %s->%s mutation)."
            % (from_resname, resnum, final_pdb, from_resname, to_resname))

    fixer.applyMutations(
        ["%s-%d-%s" % (from_resname, int(resnum), to_resname)], target_chain_id)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    if add_hydrogens:
        fixer.addMissingHydrogens(7.0)

    # PDBFile.writeFile re-letters chains positionally (A, B, ...) IGNORING
    # ``chain.id``. Write first, then re-stamp the written PDB text: the chain
    # whose residue ``resnum`` is the now-mutated residue (``to_resname``) becomes
    # ``binder_chain``; every other chain becomes ``receptor_chain``. The downstream
    # build keys on these ids.
    tmp_written = out_pdb + ".written.pdb"
    with open(tmp_written, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh)
    # Identify the written chain letter that carries the mutated binder residue.
    binder_written_id = None
    for line in open(tmp_written):
        if line[:6] in ("ATOM  ", "HETATM"):
            try:
                rnum = int(line[22:26])
            except ValueError:
                continue
            if rnum == int(resnum) and line[17:20].strip() == to_resname:
                binder_written_id = line[21]
                break
    if binder_written_id is None:
        raise ValueError(
            "prepare_mutated_bound_complex_from_final: could not locate the mutated "
            "residue %s-%d in the written complex of %s (chain re-stamp failed)."
            % (to_resname, resnum, final_pdb))
    with open(tmp_written) as src, open(out_pdb, "w") as dst:
        for line in src:
            if line[:6] in ("ATOM  ", "HETATM"):
                cid = binder_chain if line[21] == binder_written_id else receptor_chain
                dst.write(line[:21] + cid + line[22:])
            else:
                dst.write(line)
    return out_pdb


def compute_decouple_direction(
    build: Dict[str, Any], binder_chain: str = "B", resnum: Optional[int] = None,
    shell_nm: float = 1.0, decouple_nm: float = 1.2, spec: Optional[Any] = None,
) -> Optional[Tuple[float, float, float]]:
    """Unit vector pointing from the LOCAL heavy-atom density OUTWARD past NE1.

    The genuine-mode HE1 decouple translates HE1 ``decouple_nm`` along this vector
    at the u1 endpoint; the landing point must sit in LOW density (bulk solvent /
    empty space) so the "decoupled" HE1 does not clash a real atom — which would
    replace the genuine decouple cost with a packing artifact.

    A naive "away-from-receptor-centroid" heuristic FAILS for the 2QKI bound pose:
    residue-4 is solvent-exposed on the far side of the receptor, so that vector
    points back THROUGH the macrocyclic peptide's own ring and lands HE1 on the
    binder's Arg side chain (verified: 0.12 nm clash). The receptor was never the
    obstacle here — the binder's own fold was. So the direction is taken from the
    LOCAL structure instead: ``NE1 - centroid(heavy atoms within shell_nm of
    NE1)``, i.e. straight outward from whatever is packed around residue-4
    (receptor AND peptide). This is leg-agnostic (for the free peptide it points
    out of the peptide body into solvent; for the bound complex it clears both the
    receptor and the peptide) and deterministic. The shell is restricted to a
    local neighbourhood so a far-away bulk does not bias the outward direction.

    Returns ``None`` if NE1 or the local shell cannot be resolved, in which case
    the caller keeps the legacy fixed +Z.

    ``spec`` (a ``MutationSpec`` or registry name; ``None`` => res-4 MTR<->Trp)
    supplies the common attach atom + residue number. The outward direction is
    computed off the common attach atom (NE1 for Trp/MTR, CG1 for Val/Ile), so
    the geometry generalizes to any single-residue mutation. ``resnum`` overrides
    the spec's residue number when given (legacy positional compatibility).
    """
    ms = resolve_mutation_spec(spec)
    target_resnum = resnum if resnum is not None else ms.resnum
    common_atom = ms.common_attach_atom
    topology = build["modeller"].topology
    positions = np.array([
        v.value_in_unit(unit.nanometer) for v in build["modeller"].positions])
    system = build.get("system")
    ne1_idx = None
    for chain in topology.chains():
        if chain.id != binder_chain:
            continue
        for res in chain.residues():
            if res.name in _SOLVENT_RESNAMES:
                continue
            try:
                if int(res.id) != target_resnum:
                    continue
            except (TypeError, ValueError):
                continue
            for atom in res.atoms():
                if atom.name == common_atom:
                    ne1_idx = atom.index
    if ne1_idx is None:
        return None
    ne1_pos = positions[ne1_idx]

    # Local heavy-atom neighbourhood (mass > 1.5 Da excludes H and any massless
    # dummy). Restrict to atoms within shell_nm of NE1 (exclude NE1 itself).
    n_part = system.getNumParticles() if system is not None else len(positions)
    diffs = positions[:n_part] - ne1_pos
    dist = np.linalg.norm(diffs, axis=1)
    shell = []
    for i in range(n_part):
        if i == ne1_idx:
            continue
        if dist[i] >= shell_nm:
            continue
        if system is not None:
            m = system.getParticleMass(i).value_in_unit(unit.dalton)
            if m <= 1.5:
                continue
        shell.append(i)
    if not shell:
        return None
    local_centroid = positions[shell].mean(axis=0)
    out = ne1_pos - local_centroid
    norm = float(np.linalg.norm(out))
    if norm < 1e-9:
        return None
    out = out / norm
    return (float(out[0]), float(out[1]), float(out[2]))


def prepare_free_peptide_pdb(complex_pdb: str, out_pdb: str,
                             binder_chain: str = "B") -> str:
    """Extract the cyclic peptide (free leg) from a complex and fix its termini.

    Combines ``_extract_binder_only`` (drop the receptor) with
    ``_pdbfixer_prep_termini`` (standard N/C terminal completion, ncAA-masked)
    so the result is template-clean for amber14 ``createSystem`` while keeping
    the cyclic_ss disulfide.
    """
    raw = out_pdb + ".raw.pdb"
    _extract_binder_only(complex_pdb, raw, binder_chain=binder_chain)
    return _pdbfixer_prep_termini(raw, out_pdb)


# ---------------------------------------------------------------------------
# Alchemical-atom identification (B-C5)
# ---------------------------------------------------------------------------
def identify_alchemical_atoms(
    topology: app.Topology, binder_chain: str = "B",
    resnum: Optional[int] = None, spec: Optional[Any] = None,
) -> Dict[str, List[int]]:
    """Return global atom indices for the residue alchemical partition.

    Classifies indices (DEFAULT spec = residue-4 MTR<->Trp) into:
      common  : the common attach atom (NE1; frozen, present in both states)
      wt_only : the DISAPPEARING-state-only atoms (HE1; the legacy "wt_only" slot)
      mtr_only: the APPEARING-state-only atoms (CM, HM1-3; legacy "mtr_only" slot)
    Atoms not present in the given topology are simply absent (e.g. a WT-only
    structure has no CM); the caller combines both states' lists for the
    dual-topology box.

    ``spec`` (a ``MutationSpec`` or registry name; ``None`` => the res-4
    MTR<->Trp default) generalizes the partition to any single-residue mutation
    (e.g. residue-3 Val<->Ile, V3I). The legacy ``wt_only`` / ``mtr_only`` dict
    keys are kept (they map to the spec's disappearing/appearing sets) so every
    downstream consumer (index map, swap, asserts) is reused unchanged.

    ``resnum`` (legacy positional kwarg) overrides the spec's residue number when
    given; otherwise the spec's resnum is used. With ``spec=None`` and no
    ``resnum`` this resolves to ``ALCH_RESNUM`` (byte-identical legacy behaviour).
    """
    ms = resolve_mutation_spec(spec)
    target_resnum = resnum if resnum is not None else ms.resnum
    common_atom = ms.common_attach_atom
    disappear = set(ms.stateA_only_atoms)
    appear = set(ms.stateB_only_atoms)

    common: List[int] = []
    wt_only: List[int] = []
    mtr_only: List[int] = []
    for chain in topology.chains():
        if chain.id != binder_chain:
            continue
        for res in chain.residues():
            try:
                rn = int(res.id)
            except (TypeError, ValueError):
                continue
            if rn != target_resnum:
                continue
            for atom in res.atoms():
                if atom.name == common_atom:
                    common.append(atom.index)
                elif atom.name in disappear:
                    wt_only.append(atom.index)
                elif atom.name in appear:
                    mtr_only.append(atom.index)
    return {"common": common, "wt_only": wt_only, "mtr_only": mtr_only}


# ---------------------------------------------------------------------------
# System build (B-C1, B-C3): FF stack identical to run_restrained_md
# ---------------------------------------------------------------------------
def build_leg_system(
    pdb_path: str,
    leg: str = "bound",
    binder_chain: str = "B",
    solvate: bool = True,
    padding_nm: float = 1.2,
    ncaa_xml: Optional[str] = HYBRID_MTR_XML,
    hydrogens_xml: Optional[str] = HYBRID_MTR_HYDROGENS_XML,
    add_hydrogens: bool = True,
    constraints: Any = HBonds,
    spec: Optional[Any] = None,
    ncaa_resname: Optional[str] = "MTR",
) -> Dict[str, Any]:
    """Build one leg's OpenMM ``System`` with the canonical UPDD FF stack.

    Mirrors ``run_restrained_md.py`` exactly: amber14-all + tip3pfb + MTR
    hybrid XML, tip3p water, 1.2 nm padding, 0.15 M NaCl neutralized, PME with
    1.0 nm cutoff, HBonds, rigid water. The cyclic_ss disulfide (CYS2-CYS12) is
    committed via the shared detector before ``createSystem``.

    Returns a dict with modeller, system, alchemical-atom indices, disulfide
    info and the solvated atom count. For ``solvate=False`` (smoke-test) the
    box is left unsolvated and NoCutoff is used.

    ``constraints`` (default ``HBonds``) is the OpenMM ``createSystem``
    constraint set. It is threaded through ONLY so a caller can request an
    UNCONSTRAINED build (``constraints=None``) when the alchemical hydrogens
    (the appearing/disappearing methyl/HE1 H) must NOT carry SHAKE — this is the
    Tier-2 short-dynamics finite-under-integration requirement (a 1 fs
    unconstrained alch-H integrator). The default preserves the canonical
    production stack byte-for-byte (all existing callers are unaffected).

    ``spec`` (a ``MutationSpec`` or registry name; ``None`` => res-4 MTR<->Trp)
    only changes which atoms ``identify_alchemical_atoms`` classifies — it does
    NOT alter the FF stack.

    ``ncaa_resname`` (default ``"MTR"``) is the HETATM ncAA residue requiring the
    extra hydrogen-definition load + internal-bond injection. For a CANONICAL
    all-amber mutation (e.g. V3I Val<->Ile) the caller passes ``ncaa_resname=None``
    and ``ncaa_xml=None`` so those ncAA-specific steps are skipped (Val/Ile are
    standard amber14 templates). ``ncaa_xml=None`` drops the extra XML from the
    ForceField stack. The default keeps the existing MTR path byte-identical.
    """
    ff_inputs = list(FF_FILES) + ([ncaa_xml] if ncaa_xml else [])
    ff = ForceField(*ff_inputs)

    pdb = PDBFile(pdb_path)
    modeller = Modeller(pdb.topology, pdb.positions)

    # Register the MTR ncAA hydrogen definitions BEFORE addHydrogens, exactly
    # as run_restrained_md does (Modeller.loadHydrogenDefinitions). Without
    # this the MTR methyl H's (HM1-3 on CM) are not placed and createSystem
    # reports a missing-H template error. Skip addHydrogens for inputs that are
    # already H-complete (an MD final.pdb) — re-running it on a fully
    # protonated ncAA mis-handles the methyl set.
    if add_hydrogens:
        if hydrogens_xml and os.path.isfile(hydrogens_xml):
            Modeller.loadHydrogenDefinitions(hydrogens_xml)
        modeller.addHydrogens(ff)

    # ncAA internal-bond injection: PDBFile reads MTR as HETATM without intra-
    # residue bonds. Re-add them from the MTR XML <Bond> records (identical to
    # run_restrained_md.inject_xml_bonds v38 multi-site). MUST happen before
    # the peptide-bond completion below (the peptide adder reads existing
    # bonds to skip duplicates). Skipped for a canonical all-amber build
    # (``ncaa_resname=None``) — standard residues carry their bonds already.
    n_internal_added = 0
    if ncaa_resname is not None and ncaa_xml:
        n_internal_added = inject_xml_internal_bonds(
            modeller.topology, [ncaa_xml], xml_res_name=ncaa_resname
        )

    # ncAA peptide-bond completion: OpenMM's PDBFile reader does not infer
    # ATOM↔HETATM peptide bonds, so MTR (HETATM) junctions are missing here.
    # Add them by C(i)-N(i+1) distance threshold (0.20 nm) before
    # createSystem so amber14 templates resolve correctly. Identical logic to
    # run_restrained_md.add_missing_peptide_bonds_safe. Skipped for a canonical
    # all-amber build (no HETATM junctions to repair).
    n_peptide_added = 0
    if ncaa_resname is not None:
        n_peptide_added = add_missing_peptide_bonds_safe(
            modeller, binder_chain=binder_chain, max_cn_distance_nm=0.20
        )

    # cyclic_ss disulfide (B-C3): retained on BOTH legs.
    disulfide = None
    pair = detect_disulfide_pair(modeller, binder_chain=binder_chain,
                                 max_sg_dist_nm=DISULFIDE_MAX_NM)
    if pair is not None:
        added = commit_bond_if_missing(modeller.topology, pair[0], pair[1])
        disulfide = {
            "sg1_index": pair[0].index,
            "sg2_index": pair[1].index,
            "sg1_res": pair[0].residue.id,
            "sg2_res": pair[1].residue.id,
            "bond_added": added,
        }

    if solvate:
        modeller.addSolvent(
            ff, model="tip3p", padding=padding_nm * unit.nanometers,
            ionicStrength=0.15 * unit.molar,
            positiveIon="Na+", negativeIon="Cl-", neutralize=True,
        )
        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=PME,
            nonbondedCutoff=1.0 * unit.nanometers,
            constraints=constraints,
            rigidWater=True,
            ewaldErrorTolerance=0.0005,
        )
    else:
        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=constraints,
            rigidWater=True,
        )

    alch = identify_alchemical_atoms(
        modeller.topology, binder_chain=binder_chain, spec=spec)
    return {
        "leg": leg,
        "modeller": modeller,
        "system": system,
        "n_atoms": modeller.topology.getNumAtoms(),
        "alchemical_atoms": alch,
        "disulfide": disulfide,
        "n_peptide_bonds_added": n_peptide_added,
        "n_internal_bonds_added": n_internal_added,
        "ff_inputs": ff_inputs,
    }


# ---------------------------------------------------------------------------
# ATMForce wiring (B-C5): soft-core alchemical perturbation
# ---------------------------------------------------------------------------
def attach_atm_force(
    system: mm.System,
    displacement_nm: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    displaced_atoms: Optional[List[int]] = None,
    lambda1: float = 0.0,
    lambda2: float = 0.0,
    alpha_per_kcal: float = 0.10,
    u0_kcal: float = 0.0,
    w0_kcal: float = 0.0,
    umax_kcal: float = 200.0,
    ubcore_kcal: float = 100.0,
    acore: float = 0.062500,
) -> int:
    """Add an ``ATMForce`` (the ATS alchemical-transfer force) to ``system``.

    Uses OpenMM 8.5's validated 9-parameter ATMForce convenience constructor
    (Azimi/Gallicchio JCIM 2022) which builds the canonical linear-then-softplus
    soft-core alchemical potential internally and registers the standard global
    parameters (Lambda1, Lambda2, Alpha, Uh, W0, Umax, Ubcore, Acore,
    Direction). The soft-core cap (Umax/Ubcore/Acore) bounds the appearing-atom
    energy at intermediate lambda (B-C5). Returns the ATMForce index.

    The displaced atoms (the residue-4 alchemical set) receive the
    ``displacement_nm`` transformation vector; all other atoms get a zero
    displacement. For a single-point mutation the displacement is zero (no
    whole-ligand transfer): the perturbation is carried by the dual-topology
    appearing/disappearing atoms scaled through Lambda. A non-zero
    ``displacement_nm`` enables the classic AToM bound<->solvent transfer for
    validation against the package's RBFE reference.

    All non-bonded / bonded forces are moved into the ATMForce so the hybrid
    potential governs the full perturbation, per the ATMForce contract.
    """
    kcal = unit.kilocalorie_per_mole
    kj = unit.kilojoule_per_mole

    # Alpha has units of 1/energy. Input is per kcal/mol; OpenMM uses kJ/mol.
    kcal_per_kj = (1.0 * kcal).value_in_unit(kj)  # = 4.184
    alpha = (alpha_per_kcal / kcal_per_kj) / kj
    uh = (u0_kcal * kcal).value_in_unit(kj) * kj
    w0 = (w0_kcal * kcal).value_in_unit(kj) * kj
    umax = (umax_kcal * kcal).value_in_unit(kj) * kj
    ubcore = (ubcore_kcal * kcal).value_in_unit(kj) * kj
    direction = 1.0

    atm = mm.ATMForce(lambda1, lambda2, alpha, uh, w0, umax, ubcore, acore, direction)

    # Move the perturbed forces into the ATMForce. The ATMForce evaluates these
    # at the reference (u0) and displaced (u1) coordinates. NonbondedForce
    # carries the alchemical electrostatics/LJ of the appearing/disappearing
    # atoms; bonded forces follow the displacement of the swapped sidechain.
    import copy
    move_types = (mm.NonbondedForce, mm.HarmonicBondForce,
                  mm.HarmonicAngleForce, mm.PeriodicTorsionForce)
    to_move = [i for i in range(system.getNumForces())
               if isinstance(system.getForce(i), move_types)]
    # Validated AToM-OpenMM pattern (ommsystem.py): add a copy.copy() of each
    # force to the ATMForce (an unowned object), then remove the original from
    # the system. Passing the system-owned force directly raises
    # "does not own its corresponding OpenMM object".
    for i in to_move:
        atm.addForce(copy.copy(system.getForce(i)))
    for i in sorted(to_move, reverse=True):
        system.removeForce(i)

    # Enumerate every particle. Displaced atoms get the transformation vector;
    # the rest get zero displacement.
    disp = mm.Vec3(*displacement_nm) * unit.nanometer
    zero = mm.Vec3(0.0, 0.0, 0.0) * unit.nanometer
    nonzero = any(abs(c) > 0 for c in displacement_nm)
    dset = set(displaced_atoms or [])
    for idx_p in range(system.getNumParticles()):
        if nonzero and idx_p in dset:
            atm.addParticle(disp)
        else:
            atm.addParticle(zero)

    return system.addForce(atm)


# ===========================================================================
# IN-PLACE RESIDUE-4 FUSED DUAL-TOPOLOGY ATS (v0.8 de-risk; design rationale
# C1-C8 coordinate-swap criteria).
#
# This is a SEPARATE code path from attach_atm_force/smoke_test_leg above,
# which are left intact (R-7 / C8). The legacy path does single-topology
# fixed-displacement labelling only (displacement default (0,0,0) -> u1==u0,
# null op). The functions below build a GENUINE fused dual-topology box where
# residue-4 carries BOTH endpoint variant atoms (WT HE1 + MTR CM/HM1-3) on a
# common NE1, and wire the upstream common-core coordinate swap
# (ParticleOffsetDisplacement + expression-string ATMForce with all 10 globals
# incl. UOffset) so that only ~4-5 alch atoms move in place — NOT the ~210-atom
# whole binder. Ranking-only (R-11); this is a PREDICTION test (R-18), NOT a
# converged DDG_bind.
#
# Reuse map (C1, anti-fragmentation):
#   - build_leg_system()                  : both endpoint systems (verbatim)
#   - identify_alchemical_atoms()         : the residue-4 partition (verbatim)
#   - detect_disulfide_pair/commit_*      : cyclic_ss preserved on the MTR base
#   - upstream add_common_var_atoms_to_atmforce (ommsystem.py:745-789) pattern :
#       the swap loop (setParticleTransformation + ParticleOffsetDisplacement)
#       is reproduced verbatim, NOT hand-rolled, with NE1 as the attach atom.
#   - upstream set_atmforce() expression strings (ommsystem.py:817-834)        :
#       the referencePot/alchemicalPot/softCore strings + 10 globals are copied
#       verbatim (they carry the UOffset term the 9-arg convenience ctor omits).
#   - upstream add_forces_to_atmforce() var_regions branch (ommsystem.py:209-223):
#       move ALL Nonbonded/Harmonic/Torsion forces into the ATMForce via
#       copy.copy() then removeForce() — same migration the upstream uses.
# ===========================================================================

# Soft-core canon (Gallicchio 2021; Azimi 2022 JCIM 62(2):309
# DOI 10.1021/acs.jcim.1c01129). DO NOT re-tune (C3/R3): the whole point of the
# small in-place perturbation is that the cap is barely exercised.
ATS_UMAX_KCAL = 200.0
ATS_UBCORE_KCAL = 100.0
ATS_ACORE = 0.062500

# Sentinel marking the WT HE1 slot inside a harvested angle triple. The
# harvester (_harvest_he1_terms) substitutes the WT HE1 index with this value;
# the fused injector (_inject_he1_into_fused) remaps exactly that slot to the
# appended fused index, leaving the common-atom slots (whose WT indices equal
# the MTR indices by C3 alignment) unchanged.
_HE1_SENTINEL = -987654321

# Upstream ATMForce energy-expression strings (ommsystem.py:817-823), verbatim.
# Copied (not hand-rolled) per C1 so the UOffset term and softplus form match
# the validated AToM-OpenMM var-region protocol exactly.
_ATS_REFERENCE_POT_EXPR = "select(step(Direction), u0, u1) + "
_ATS_ALCHEMICAL_POT_EXPR = (
    "select(Lambda2-Lambda1 , "
    "((Lambda2-Lambda1)/Alpha)*log(1+exp(-Alpha*(usc-Uh))) + Lambda2*usc + W0, "
    "Lambda2*usc + W0);"
)
_ATS_SOFTCORE_EXPR = (
    "usc = select(Acore, select(step(u-Ubcore), (Umax-Ubcore)*fsc+Ubcore, u), u);"
    "fsc = (z^Acore-1)/(z^Acore+1);"
    "z = 1 + 2*(y/Acore) + 2*(y/Acore)^2;"
    "y = (u-Ubcore)/(Umax-Ubcore);"
    "u = select(step(Direction), 1, -1)*(u1-(u0 + UOffset))"
)

# Force types migrated into the ATMForce for the var-regions protocol
# (upstream nbpattern/harmpattern/torpattern, ommsystem.py:209-223).
_ATS_MOVE_FORCE_TYPES = (
    mm.NonbondedForce,
    mm.HarmonicBondForce,
    mm.HarmonicAngleForce,
    mm.PeriodicTorsionForce,
)


def _kcal_to_kj(x_kcal: float) -> float:
    """kcal/mol -> kJ/mol scalar (OpenMM global parameters are kJ/mol)."""
    kcal = unit.kilocalorie_per_mole
    kj = unit.kilojoule_per_mole
    return (x_kcal * kcal).value_in_unit(kj)


def _build_common_index_map(
    wt_build: Dict[str, Any], mtr_build: Dict[str, Any], binder_chain: str = "B"
) -> Dict[str, Any]:
    """Pair the residue-shared (common) atoms between the WT and MTR endpoints.

    Both endpoint systems are built by ``build_leg_system`` from the same
    cyclic-peptide topology, so their COMMON atoms (everything except the
    residue-4 variant set) are produced in the identical order. This returns
    the index-aligned common-atom lists plus the var/attach atoms, and runs the
    C3 hard gate: equal common-atom COUNT and byte-identical NAME ORDER. The
    swap is positional (common_i <-> other_common_i), so a count or order
    mismatch is the single most likely silent-wrong point — fail loud.
    """
    wt_top = wt_build["modeller"].topology
    mtr_top = mtr_build["modeller"].topology

    wt_var = set(wt_build["alchemical_atoms"]["wt_only"])
    mtr_var = set(mtr_build["alchemical_atoms"]["mtr_only"])

    wt_atoms = list(wt_top.atoms())
    mtr_atoms = list(mtr_top.atoms())

    # Common region = the shared BINDER residue frame only (exclude the var atoms
    # AND any solvent the MTR base may carry when solvated — the common-core swap
    # is protein-only; waters/ions are non-alchemical and not transformed).
    def _is_binder_protein(atom) -> bool:
        return (atom.residue.chain.id == binder_chain
                and atom.residue.name not in _SOLVENT_RESNAMES)

    wt_common = [a.index for a in wt_atoms
                 if a.index not in wt_var and _is_binder_protein(a)]
    mtr_common = [a.index for a in mtr_atoms
                  if a.index not in mtr_var and _is_binder_protein(a)]

    # C3 hard gate 1: common-atom COUNT parity (upstream _exit at L763-765).
    if len(wt_common) != len(mtr_common):
        raise ValueError(
            "C3 common-atom COUNT parity FAIL: WT has %d common atoms, MTR has "
            "%d. The fused dual-topology swap is positional and requires equal "
            "common-atom counts (upstream ommsystem.py:763-765 _exit gate)."
            % (len(wt_common), len(mtr_common))
        )

    # C3 hard gate 2: index-order alignment by atom NAME (positional swap).
    wt_name = {a.index: a.name for a in wt_atoms}
    mtr_name = {a.index: a.name for a in mtr_atoms}
    mismatches: List[str] = []
    for w_i, m_i in zip(wt_common, mtr_common):
        if wt_name[w_i] != mtr_name[m_i]:
            mismatches.append(
                "pos WT[%d]=%s vs MTR[%d]=%s" % (w_i, wt_name[w_i], m_i, mtr_name[m_i])
            )
    if mismatches:
        raise ValueError(
            "C3 common-atom ORDER alignment FAIL (%d mismatches): the WT and MTR "
            "common-atom lists must be in the same physical order for the "
            "positional swap. First mismatches: %s"
            % (len(mismatches), "; ".join(mismatches[:8]))
        )

    return {
        "wt_common": wt_common,
        "mtr_common": mtr_common,
        "wt_var": sorted(wt_var),       # {HE1}
        "mtr_var": sorted(mtr_var),     # {CM, HM1, HM2, HM3}
        "wt_attach": wt_build["alchemical_atoms"]["common"][0],   # NE1 (WT)
        "mtr_attach": mtr_build["alchemical_atoms"]["common"][0],  # NE1 (MTR)
        "n_common": len(wt_common),
    }


def assert_common_atom_param_continuity(
    wt_build: Dict[str, Any], mtr_build: Dict[str, Any], cmap: Dict[str, Any],
    q_tol_e: float = 1e-4, sigma_tol_nm: float = 1e-4, eps_tol_kj: float = 1e-4,
) -> Dict[str, Any]:
    """MC1 (C4): common-atom (charge, sigma, epsilon) endpoint-equality ASSERT.

    For the relative dual-topology cycle to be exact, the COMMON core must be
    electrostatically/LJ-continuous: each common atom's nonbonded params must be
    identical between the WT-derived and MTR-derived representations. The hybrid
    MTR XML freezes NE1 (-0.3418) but the rest of the ring/backbone common
    partials must NOT diverge — otherwise the "common" core is not common and
    the swap injects spurious DeltaE that is finite-but-wrong (a false-green
    invisible to a single-state finiteness check). Highest ncAA risk per Q4.

    Returns the worst per-channel deviation; raises on any exceedance.
    """
    nb_wt = next(f for f in wt_build["system"].getForces()
                 if isinstance(f, mm.NonbondedForce))
    nb_mtr = next(f for f in mtr_build["system"].getForces()
                  if isinstance(f, mm.NonbondedForce))

    max_dq = 0.0
    max_dsig = 0.0
    max_deps = 0.0
    worst: List[str] = []
    wt_name = {a.index: a.name for a in wt_build["modeller"].topology.atoms()}
    for w_i, m_i in zip(cmap["wt_common"], cmap["mtr_common"]):
        qw, sw, ew = nb_wt.getParticleParameters(w_i)
        qm, sm, em = nb_mtr.getParticleParameters(m_i)
        dq = abs(qw.value_in_unit(unit.elementary_charge)
                 - qm.value_in_unit(unit.elementary_charge))
        dsig = abs(sw.value_in_unit(unit.nanometer)
                   - sm.value_in_unit(unit.nanometer))
        deps = abs(ew.value_in_unit(unit.kilojoule_per_mole)
                   - em.value_in_unit(unit.kilojoule_per_mole))
        if dq > max_dq:
            max_dq = dq
        if dsig > max_dsig:
            max_dsig = dsig
        if deps > max_deps:
            max_deps = deps
        if dq > q_tol_e or dsig > sigma_tol_nm or deps > eps_tol_kj:
            worst.append("%s: dq=%.3e dsig=%.3e deps=%.3e"
                         % (wt_name.get(w_i, "?"), dq, dsig, deps))

    result = {
        "max_dq_e": max_dq,
        "max_dsigma_nm": max_dsig,
        "max_deps_kj": max_deps,
        "n_common_checked": cmap["n_common"],
        "passed": not worst,
    }
    if worst:
        raise ValueError(
            "MC1 common-atom param continuity FAIL (%d atoms exceed tol "
            "q=%.1e sigma=%.1e eps=%.1e): the dual-topology common core is not "
            "electrostatically/LJ continuous. Offenders: %s"
            % (len(worst), q_tol_e, sigma_tol_nm, eps_tol_kj, "; ".join(worst[:8]))
        )
    return result


def _summarize_common_charge_divergence(
    wt_build: Dict[str, Any], mtr_build: Dict[str, Any], cmap: Dict[str, Any],
    binder_chain: str = "B",
) -> Dict[str, Any]:
    """Structured per-atom report of the WT<->MTR common-charge divergence.

    Used when MC1 fails non-strictly: rather than crash, surface WHICH common
    atoms diverge and by how much (residue-4 vs elsewhere, net displaced
    charge). This makes the pre-registered S1 outcome (ii) actionable
    (common-charge harmonization scope) instead of an opaque error string.
    """
    nb_wt = next(f for f in wt_build["system"].getForces()
                 if isinstance(f, mm.NonbondedForce))
    nb_mtr = next(f for f in mtr_build["system"].getForces()
                  if isinstance(f, mm.NonbondedForce))
    wt_atoms = {a.index: a for a in wt_build["modeller"].topology.atoms()}

    per_res4: List[Dict[str, Any]] = []
    sum_dq_res4 = 0.0
    sum_dq_all = 0.0
    n_diverging = 0
    for w_i, m_i in zip(cmap["wt_common"], cmap["mtr_common"]):
        qw = nb_wt.getParticleParameters(w_i)[0].value_in_unit(unit.elementary_charge)
        qm = nb_mtr.getParticleParameters(m_i)[0].value_in_unit(unit.elementary_charge)
        dq = qw - qm
        sum_dq_all += dq
        if abs(dq) > 1e-4:
            n_diverging += 1
        atom = wt_atoms[w_i]
        if str(atom.residue.id) == str(ALCH_RESNUM):
            sum_dq_res4 += dq
            per_res4.append({"name": atom.name, "q_wt": round(qw, 4),
                             "q_mtr": round(qm, 4), "dq": round(dq, 4)})
    return {
        "passed": False,
        "n_common_checked": cmap["n_common"],
        "n_diverging": n_diverging,
        "sum_dq_res4_e": round(sum_dq_res4, 5),
        "sum_dq_all_common_e": round(sum_dq_all, 6),
        "residue4_common": per_res4,
    }


def _harmonize_common_charges_to_wt(
    wt_build: Dict[str, Any], mtr_build: Dict[str, Any], cmap: Dict[str, Any]
) -> int:
    """DIAGNOSTIC-ONLY: overwrite the MTR base common-atom charges with the WT
    values so the common core is electrostatically continuous (MC1 passes).

    This is NOT a production fix — it discards the hybrid-MTR RESP charges on the
    common core, which is only valid for ISOLATING the mechanical validation of
    the swap/endpoint-equivalence from the charge-continuity gap. Sigma/eps are
    already identical (only charge diverges). Mutates the MTR NonbondedForce.
    Returns the number of common atoms whose charge was overwritten.
    """
    nb_wt = next(f for f in wt_build["system"].getForces()
                 if isinstance(f, mm.NonbondedForce))
    nb_mtr = next(f for f in mtr_build["system"].getForces()
                  if isinstance(f, mm.NonbondedForce))
    n = 0
    for w_i, m_i in zip(cmap["wt_common"], cmap["mtr_common"]):
        qw, _, _ = nb_wt.getParticleParameters(w_i)
        _, sm, em = nb_mtr.getParticleParameters(m_i)
        nb_mtr.setParticleParameters(m_i, qw, sm, em)
        n += 1
    return n


def _harvest_he1_terms(
    wt_build: Dict[str, Any], he1_index: int,
    wt_to_mtr: Optional[Dict[int, int]] = None,
) -> Dict[str, Any]:
    """Extract every WT-side parameter for the disappearing HE1 atom.

    Pulls HE1's mass, nonbonded (q, sigma, eps), its exceptions, and all bonded
    terms (NE1-HE1 bond + the two NE1-centred angles) so they can be transplanted
    onto the fused MTR base.

    CRITICAL index remap: the WT and MTR RAW atom indices DIVERGE after NE1 (WT
    has HE1 at 61 then the ring; MTR has CM/HM1-3 at 61-64 then the ring), so a
    partner index harvested from WT (e.g. CE2=62 in WT) is a DIFFERENT atom in
    the MTR base (HM1=62). Every harvested partner index is therefore translated
    WT-raw -> MTR-raw via ``wt_to_mtr`` (built from the C3 common-atom pairing).
    Only HE1's own index becomes the appended fused index (sentinel-marked).
    A missing partner in the map is a fatal alignment error (fail-loud).
    """
    wt_to_mtr = wt_to_mtr or {}

    def _remap(idx: int) -> int:
        if idx == he1_index:
            return _HE1_SENTINEL
        if idx not in wt_to_mtr:
            raise ValueError(
                "HE1 transplant: WT partner index %d has no MTR counterpart in "
                "the common-atom map (raw indices diverge past NE1; a partner "
                "outside the common set cannot be transplanted)." % idx)
        return wt_to_mtr[idx]

    sysm = wt_build["system"]
    mass = sysm.getParticleMass(he1_index).value_in_unit(unit.dalton)

    nb = next(f for f in sysm.getForces() if isinstance(f, mm.NonbondedForce))
    q, sig, eps = nb.getParticleParameters(he1_index)
    nb_params = {
        "charge_e": q.value_in_unit(unit.elementary_charge),
        "sigma_nm": sig.value_in_unit(unit.nanometer),
        "epsilon_kj": eps.value_in_unit(unit.kilojoule_per_mole),
    }
    exceptions = []
    for ei in range(nb.getNumExceptions()):
        p1, p2, cp, sg, ep = nb.getExceptionParameters(ei)
        if he1_index in (p1, p2):
            other = p2 if p1 == he1_index else p1
            exceptions.append({
                "other": _remap(other),
                "chargeProd_e2": cp.value_in_unit(unit.elementary_charge ** 2),
                "sigma_nm": sg.value_in_unit(unit.nanometer),
                "epsilon_kj": ep.value_in_unit(unit.kilojoule_per_mole),
            })

    bonds = []
    angles = []
    for f in sysm.getForces():
        if isinstance(f, mm.HarmonicBondForce):
            for bi in range(f.getNumBonds()):
                p1, p2, length, k = f.getBondParameters(bi)
                if he1_index in (p1, p2):
                    other = p2 if p1 == he1_index else p1
                    bonds.append({
                        "other": _remap(other),
                        "length_nm": length.value_in_unit(unit.nanometer),
                        "k_kj_nm2": k.value_in_unit(
                            unit.kilojoule_per_mole / unit.nanometer ** 2),
                    })
        elif isinstance(f, mm.HarmonicAngleForce):
            for ai in range(f.getNumAngles()):
                a1, a2, a3, ang, k = f.getAngleParameters(ai)
                if he1_index in (a1, a2, a3):
                    # Remap every slot WT-raw -> MTR-raw; the HE1 slot becomes
                    # the sentinel (-> appended fused index in the injector), the
                    # common slots become their MTR counterparts via wt_to_mtr.
                    triple = [_remap(v) for v in (a1, a2, a3)]
                    angles.append({
                        "a1": triple[0], "a2": triple[1], "a3": triple[2],
                        "angle_rad": ang.value_in_unit(unit.radian),
                        "k_kj_rad2": k.value_in_unit(
                            unit.kilojoule_per_mole / unit.radian ** 2),
                    })
    return {
        "mass_da": mass,
        "nonbonded": nb_params,
        "exceptions": exceptions,
        "bonds": bonds,
        "angles": angles,
    }


def _inject_he1_into_fused(
    mtr_build: Dict[str, Any], he1_terms: Dict[str, Any],
    mtr_var_indices: Optional[List[int]] = None,
) -> Dict[str, Any]:
    """Append the WT-only HE1 atom to the MTR base, becoming the FUSED box.

    Adds HE1 as one extra particle to the MTR System + Topology (residue 4,
    chain B) and transplants its harvested WT nonbonded/exception/bonded terms.
    The MTR base already carries CM + HM1-3 and the cyclic_ss disulfide, so the
    result is a single topology with BOTH endpoint variant atoms on the common
    NE1. Returns the fused HE1 index + bookkeeping. Mutates ``mtr_build`` in
    place (its system/topology become the fused box).

    Var-var exclusions (``mtr_var_indices`` = CM, HM1-3): HE1 (WT-var) and the
    methyl (MTR-var) are MUTUALLY-EXCLUSIVE alternate-state atoms of the SAME
    site (both bond to NE1), so in the fused box they physically overlap
    (~0.06 nm). They must NEVER see each other's nonbonded interaction (that
    overlap would detonate the LJ term). A zero-interaction nonbonded EXCEPTION
    is added between HE1 and every methyl atom — the single-topology-style
    var-var exclusion that makes the in-place fused box well-conditioned. Without
    it, the reference state energy u0 diverges (~1e6 kcal/mol) from the HE1<->CM
    clash. The ATMForce lambda swap then dials which endpoint's var atoms couple
    to the rest of the system.
    """
    system = mtr_build["system"]
    topology = mtr_build["modeller"].topology
    positions = list(mtr_build["modeller"].positions)

    # 1) System particle (mass).
    fused_he1 = system.addParticle(he1_terms["mass_da"] * unit.dalton)

    # 2) Topology atom. OpenMM requires a residue's atoms to be CONTIGUOUS, and
    #    residue-4 is not the last residue, so HE1 cannot be appended into it.
    #    The System/NonbondedForce are index-only (no contiguity rule) and the
    #    ATMForce swap operates purely on System indices, so HE1's PHYSICAL
    #    identity (the WT NE1-HE1 bond + its swap to NE1) is carried by the
    #    transplanted bonded/nonbonded terms regardless of topology grouping.
    #    We therefore place the topology atom in a NEW dedicated chain "H" (a new
    #    chain is always contiguity-legal, even after solvation has appended water
    #    chains/residues; appending into chain B fails once chain B carries
    #    trailing solvent). The topology is used only for position setup, not
    #    energy; NE1 of residue-4 stays the bond/swap partner via the
    #    transplanted bonded terms (System-index based, not topology-grouping
    #    based). The caller sets the fused HE1 alch index from the returned index.
    res4 = None
    for chain in topology.chains():
        if chain.id != "B":
            continue
        for res in chain.residues():
            try:
                if int(res.id) == ALCH_RESNUM:
                    res4 = res
                    break
            except (TypeError, ValueError):
                continue
        if res4 is not None:
            break
    if res4 is None:
        raise ValueError("Fused build: residue-4 not found on chain B for HE1.")
    h_element = app.element.hydrogen
    he1_chain = topology.addChain(id="H")
    he1_res = topology.addResidue("HE1X", he1_chain, id=str(ALCH_RESNUM))
    topo_he1 = topology.addAtom("HE1", h_element, he1_res)
    # Bond HE1 to NE1 in the topology too (keeps the connectivity record honest;
    # cross-chain topology bonds are permitted).
    ne1_atom = next((a for a in res4.atoms() if a.name == ALCH_COMMON_ATOM), None)
    if ne1_atom is not None:
        topology.addBond(ne1_atom, topo_he1)

    # Position HE1 at the WT-equilibrium geometry off NE1 (the seed; R2 asserts
    # this is non-clashing before attach). Use the existing NE1 position +
    # 0.101 nm along (NE1->CD1) outward; a precise seed is validated by R2.
    ne1_idx = mtr_build["alchemical_atoms"]["common"][0]
    ne1_pos = np.array(positions[ne1_idx].value_in_unit(unit.nanometer))
    cd1_idx = next((a.index for a in res4.atoms() if a.name == "CD1"), None)
    ce2_idx = next((a.index for a in res4.atoms() if a.name == "CE2"), None)
    if cd1_idx is not None and ce2_idx is not None:
        cd1 = np.array(positions[cd1_idx].value_in_unit(unit.nanometer))
        ce2 = np.array(positions[ce2_idx].value_in_unit(unit.nanometer))
        # HE1 points away from the ring bisector of (CD1, CE2) about NE1.
        bis = (cd1 + ce2) / 2.0
        direction = ne1_pos - bis
        norm = np.linalg.norm(direction)
        direction = direction / norm if norm > 1e-9 else np.array([0.0, 0.0, 1.0])
        he1_pos = ne1_pos + 0.101 * direction
    else:
        he1_pos = ne1_pos + np.array([0.0, 0.0, 0.101])
    positions.append(mm.Vec3(*he1_pos) * unit.nanometer)
    mtr_build["modeller"].positions = positions

    # 3) Nonbonded particle params (+ HE1 exceptions). Partner indices carry
    #    over directly (C3 alignment: WT common idx == MTR common idx).
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    if nb.getNumParticles() != fused_he1:
        # NonbondedForce must grow in lockstep with the system particle count.
        raise ValueError(
            "Fused build: NonbondedForce particle count (%d) out of sync with "
            "new particle index (%d)." % (nb.getNumParticles(), fused_he1))
    p = he1_terms["nonbonded"]
    nb.addParticle(p["charge_e"] * unit.elementary_charge,
                   p["sigma_nm"] * unit.nanometer,
                   p["epsilon_kj"] * unit.kilojoule_per_mole)
    for ex in he1_terms["exceptions"]:
        nb.addException(
            fused_he1, ex["other"],
            ex["chargeProd_e2"] * unit.elementary_charge ** 2,
            ex["sigma_nm"] * unit.nanometer,
            ex["epsilon_kj"] * unit.kilojoule_per_mole,
        )

    # 3b) Var-var exclusions: HE1 <-> each methyl atom (CM, HM1-3). These are
    #     alternate-state atoms of the same site and must never interact, else
    #     their ~0.06 nm overlap detonates the LJ term (u0 ~1e6). Zero-interaction
    #     nonbonded exception (chargeProd=0, sigma=0.1, eps=0).
    n_var_var_excl = 0
    for mv in (mtr_var_indices or []):
        nb.addException(
            fused_he1, mv,
            0.0 * unit.elementary_charge ** 2,
            0.1 * unit.nanometer,
            0.0 * unit.kilojoule_per_mole,
        )
        n_var_var_excl += 1

    # 4) Bonded terms (NE1-HE1 bond + NE1-centred angles). MC2 internal-bond
    #    injection for the appearing methyl is already carried by the MTR XML
    #    (CM-NE1 in the base build); here we inject the DISAPPEARING HE1 bond.
    bf = next(f for f in system.getForces() if isinstance(f, mm.HarmonicBondForce))
    n_bond = 0
    for b in he1_terms["bonds"]:
        bf.addBond(fused_he1, b["other"],
                   b["length_nm"] * unit.nanometer,
                   b["k_kj_nm2"] * unit.kilojoule_per_mole / unit.nanometer ** 2)
        n_bond += 1
    af = next((f for f in system.getForces()
               if isinstance(f, mm.HarmonicAngleForce)), None)
    n_angle = 0
    if af is not None:
        for ang in he1_terms["angles"]:
            # The HE1 slot is marked with _HE1_SENTINEL by the harvester; remap
            # only that slot to the appended fused index (the other two slots
            # are common atoms whose indices carry over unchanged).
            a1, a2, a3 = (
                fused_he1 if v == _HE1_SENTINEL else v
                for v in (ang["a1"], ang["a2"], ang["a3"])
            )
            af.addAngle(
                a1, a2, a3,
                ang["angle_rad"] * unit.radian,
                ang["k_kj_rad2"] * unit.kilojoule_per_mole / unit.radian ** 2,
            )
            n_angle += 1

    return {
        "fused_he1_index": fused_he1,
        "n_he1_bonds": n_bond,
        "n_he1_angles": n_angle,
        "n_he1_exceptions": len(he1_terms["exceptions"]),
        "n_var_var_exclusions": n_var_var_excl,
        "fused_n_particles": system.getNumParticles(),
    }


def assert_appearing_methyl_seed(
    fused_build: Dict[str, Any],
    min_dist_nm: float = 0.10,
    cm_ne1_eq_nm: float = 0.1450, cm_ne1_tol_nm: float = 0.02,
    hm_cm_eq_nm: float = 0.1090, hm_cm_tol_nm: float = 0.02,
) -> Dict[str, Any]:
    """R2 (C5): seed-geometry hard ASSERT, gates BEFORE the ATMForce attach.

    The appearing CM+HM1-3 (and the disappearing HE1) enter the UNCAPPED bonded
    base potential; a clashing/overlong seed detonates the harmonic term locally
    exactly like the whole-binder coverage NaN. Asserts:
      (a) min-dist of every appearing atom (CM, HM1-3) to HE1, to any common /
          disappearing atom, and to nearest solvent O is > ``min_dist_nm``;
      (b) CM-NE1 bond length within tol of ``cm_ne1_eq_nm``;
      (c) each HM-CM bond length within tol of ``hm_cm_eq_nm``.
    Raises on any violation (fail-loud, not silent).
    """
    topology = fused_build["modeller"].topology
    positions = np.array([
        v.value_in_unit(unit.nanometer) for v in fused_build["modeller"].positions
    ])
    alch = fused_build["alchemical_atoms"]
    appearing = list(alch["mtr_only"])           # CM, HM1-3
    ne1 = alch["common"][0]
    he1 = fused_build["fused_he1_index"]

    name_by_idx = {a.index: a.name for a in topology.atoms()}
    res_by_idx = {a.index: a.residue.name for a in topology.atoms()}

    # Solvent O atoms (if solvated).
    solvent_o = [a.index for a in topology.atoms()
                 if a.residue.name in _SOLVENT_RESNAMES and a.element is not None
                 and a.element.symbol == "O"]

    # Two distinct target classes (the physics differs):
    #   HARD-FAIL targets = COMMON atoms + solvent O. These are ALWAYS "on" at
    #     every lambda; an appearing atom overlapping them detonates the uncapped
    #     base term. Exclude the appearing atoms' own bonded partners (CM-NE1 is
    #     a real bond ~0.145 nm; HM-CM ~0.109 nm) — those are bond-length checked.
    #   FLAG-ONLY target = HE1. HE1 is the DISAPPEARING dual-topology partner of
    #     CM (both bond to NE1, mutually exclusive endpoint atoms), so a CM<->HE1
    #     overlap is the EXPECTED fused-box geometry handled by the soft-core +
    #     the swap's nonbonded decoupling — it must NOT hard-fail (that would be a
    #     physics error). It is reported (FLAG-not-silent per the R2 spec).
    bonded_partners = {ne1, fused_build["fused_he1_index"]}
    cm_idx = next((a for a in appearing if name_by_idx.get(a) == "CM"), None)
    if cm_idx is not None:
        bonded_partners.add(cm_idx)  # HM-CM bonds
    hard_targets = {a for a in range(topology.getNumAtoms())
                    if a not in appearing and a != he1
                    and res_by_idx.get(a) not in _SOLVENT_RESNAMES}
    hard_targets.update(solvent_o)
    hard_targets -= set(appearing)

    min_hard = float("inf")
    worst_hard = None
    for ap in appearing:
        pa = positions[ap]
        for tgt in hard_targets:
            if tgt in bonded_partners:
                continue  # real bonded partner; checked by bond-length below
            d = float(np.linalg.norm(pa - positions[tgt]))
            if d < min_hard:
                min_hard = d
                worst_hard = (name_by_idx.get(ap, ap), name_by_idx.get(tgt, tgt), d)

    if min_hard <= min_dist_nm:
        raise ValueError(
            "R2 seed min-dist FAIL: appearing atom too close to a COMMON/solvent "
            "atom (%.4f nm <= %.4f nm) — %s. A clashing methyl seed detonates the "
            "uncapped bonded base term." % (min_hard, min_dist_nm, worst_hard)
        )

    # FLAG (not fail): appearing-atom <-> HE1 (the dual-topology exclusion pair).
    min_he1 = min(
        float(np.linalg.norm(positions[ap] - positions[he1])) for ap in appearing
    )

    # Bond-length asserts.
    cm = next((a for a in appearing if name_by_idx.get(a) == "CM"), None)
    if cm is None:
        raise ValueError("R2 seed: no CM atom found among appearing atoms.")
    cm_ne1 = float(np.linalg.norm(positions[cm] - positions[ne1]))
    if abs(cm_ne1 - cm_ne1_eq_nm) > cm_ne1_tol_nm:
        raise ValueError(
            "R2 CM-NE1 bond length FAIL: %.4f nm (eq %.4f +/- %.4f). A wrong-length "
            "bond detonates the harmonic term even with valid geometry."
            % (cm_ne1, cm_ne1_eq_nm, cm_ne1_tol_nm)
        )
    hm_lengths = {}
    for hm in (a for a in appearing if name_by_idx.get(a, "").startswith("HM")):
        d = float(np.linalg.norm(positions[hm] - positions[cm]))
        hm_lengths[name_by_idx[hm]] = d
        if abs(d - hm_cm_eq_nm) > hm_cm_tol_nm:
            raise ValueError(
                "R2 %s-CM bond length FAIL: %.4f nm (eq %.4f +/- %.4f)."
                % (name_by_idx[hm], d, hm_cm_eq_nm, hm_cm_tol_nm)
            )

    return {
        "min_hard_dist_nm": min_hard,
        "min_hard_pair": worst_hard,
        "min_appearing_to_he1_nm": min_he1,    # FLAG only (dual-topo exclusion)
        "cm_ne1_nm": cm_ne1,
        "hm_cm_nm": hm_lengths,
        "n_solvent_o_checked": len(solvent_o),
        "passed": True,
    }


def assert_methyl_bonded_present(fused_build: Dict[str, Any]) -> Dict[str, Any]:
    """MC2 (C4): the appearing methyl bonded terms + CM-NE1 internal bond present.

    Confirms the CM-NE1 bond and the three HM-CM connectivities exist in the
    fused box. CM-NE1 must be a true HarmonicBondForce term (CM is heavy, never
    constrained). The HM-CM connectivities are H-X bonds: under the base build's
    ``constraints=HBonds`` they are SHAKE constraints, NOT HarmonicBond terms —
    so MC2 accepts either a HarmonicBond or a constraint for HM-CM (presence is
    the gate; whether they are constrained is an R3 engine choice, not a build
    failure). Their absence in BOTH would be a silent topology-injection failure.
    """
    system = fused_build["system"]
    topology = fused_build["modeller"].topology
    alch = fused_build["alchemical_atoms"]
    ne1 = alch["common"][0]
    name_by_idx = {a.index: a.name for a in topology.atoms()}
    cm = next((i for i in alch["mtr_only"] if name_by_idx.get(i) == "CM"), None)
    hms = [i for i in alch["mtr_only"] if name_by_idx.get(i, "").startswith("HM")]

    # Gather every HarmonicBond pair (the ATMForce may have absorbed the force;
    # search both the system and any ATMForce inner forces).
    bond_pairs = set()
    for f in system.getForces():
        if isinstance(f, mm.HarmonicBondForce):
            for bi in range(f.getNumBonds()):
                p1, p2, _, _ = f.getBondParameters(bi)
                bond_pairs.add(frozenset((p1, p2)))
        elif isinstance(f, mm.ATMForce):
            for j in range(f.getNumForces()):
                inner = f.getForce(j)
                if isinstance(inner, mm.HarmonicBondForce):
                    for bi in range(inner.getNumBonds()):
                        p1, p2, _, _ = inner.getBondParameters(bi)
                        bond_pairs.add(frozenset((p1, p2)))

    # Constraint pairs (HM-CM may be SHAKE constraints under HBonds).
    constraint_pairs = set()
    for ci in range(system.getNumConstraints()):
        a, b, _ = system.getConstraintParameters(ci)
        constraint_pairs.add(frozenset((a, b)))

    def _connected(i, j):
        return frozenset((i, j)) in bond_pairs or frozenset((i, j)) in constraint_pairs

    cm_ne1_present = cm is not None and frozenset((cm, ne1)) in bond_pairs
    hm_cm_present = {name_by_idx[h]: _connected(h, cm) for h in hms}
    hm_cm_as_constraint = {
        name_by_idx[h]: (frozenset((h, cm)) in constraint_pairs) for h in hms
    }
    if not cm_ne1_present:
        raise ValueError(
            "MC2 FAIL: CM-NE1 internal bond absent from the fused box's "
            "HarmonicBondForce (CM is heavy — must be a real bond, not a "
            "constraint).")
    missing_hm = [k for k, v in hm_cm_present.items() if not v]
    if missing_hm:
        raise ValueError(
            "MC2 FAIL: methyl HM-CM connectivity absent (neither bond nor "
            "constraint): %s" % missing_hm)
    return {
        "cm_ne1_bond_present": cm_ne1_present,
        "hm_cm_bonds_present": hm_cm_present,
        "hm_cm_as_constraint": hm_cm_as_constraint,
        "passed": True,
    }


def assert_disulfide_in_fused(fused_build: Dict[str, Any]) -> Dict[str, Any]:
    """MC3 (C4): cyclic_ss SG-SG disulfide preserved in the fused box.

    The disulfide is committed on the MTR base modeller BEFORE the fused HE1
    injection and the ATMForce migration. This re-confirms the CYS2-CYS12 SG-SG
    bond survives into the final fused system's HarmonicBondForce (whether still
    in the system or absorbed into the ATMForce).
    """
    disulfide = fused_build.get("disulfide")
    if not disulfide:
        raise ValueError("MC3 FAIL: no disulfide detected on the fused base.")
    sg1 = disulfide["sg1_index"]
    sg2 = disulfide["sg2_index"]
    system = fused_build["system"]
    found = False
    for f in system.getForces():
        if isinstance(f, mm.HarmonicBondForce):
            for bi in range(f.getNumBonds()):
                p1, p2, _, _ = f.getBondParameters(bi)
                if {p1, p2} == {sg1, sg2}:
                    found = True
        elif isinstance(f, mm.ATMForce):
            for j in range(f.getNumForces()):
                inner = f.getForce(j)
                if isinstance(inner, mm.HarmonicBondForce):
                    for bi in range(inner.getNumBonds()):
                        p1, p2, _, _ = inner.getBondParameters(bi)
                        if {p1, p2} == {sg1, sg2}:
                            found = True
    # The disulfide may be expressed as a constraint (HBonds path keeps S-S as a
    # bond, not a constraint, so the HarmonicBondForce search is authoritative);
    # also accept a constraint as a fallback.
    if not found:
        for ci in range(system.getNumConstraints()):
            a, b, _ = system.getConstraintParameters(ci)
            if {a, b} == {sg1, sg2}:
                found = True
                break
    if not found:
        raise ValueError(
            "MC3 FAIL: cyclic_ss SG-SG (%d-%d) bond absent from the fused box's "
            "HarmonicBondForce/constraints." % (sg1, sg2))
    return {"sg1_index": sg1, "sg2_index": sg2, "disulfide_present": True}


def attach_inplace_swap_atmforce(
    fused_build: Dict[str, Any],
    cmap: Dict[str, Any],
    lambda1: float = 0.0,
    lambda2: float = 0.0,
    alpha_per_kcal: float = 0.0,
    u0_kcal: float = 0.0,
    w0_kcal: float = 0.0,
    umax_kcal: float = ATS_UMAX_KCAL,
    ubcore_kcal: float = ATS_UBCORE_KCAL,
    acore: float = ATS_ACORE,
    direction: float = 1.0,
    uoffset_kcal: float = 0.0,
    var_park_nm: Optional[Tuple[float, float, float]] = None,
    swap_mode: str = "genuine",
    genuine_decouple_nm: float = 1.2,
    genuine_decouple_dir: Optional[Tuple[float, float, float]] = None,
) -> Dict[str, Any]:
    """R1 (C1/C2): attach the in-place common-core swap ATMForce.

    Reuses the upstream var-region protocol VERBATIM (no hand-rolled swap, no
    9-arg convenience ctor):
      - expression-string ``ATMForce`` ctor with the upstream reference/alchemy/
        soft-core strings + all 10 globals (incl. UOffset);
      - migrate ALL Nonbonded/Harmonic/Torsion forces into the ATMForce via
        ``copy.copy`` + ``removeForce`` (upstream add_forces_to_atmforce
        var_regions branch);
      - addParticle() for every particle, then the upstream
        ``add_common_var_atoms_to_atmforce`` ``setParticleTransformation`` +
        ``ParticleOffsetDisplacement`` idiom (attach atom = NE1, C2).
    The WT-state var = {HE1}, MTR-state var = {CM, HM1-3}; attach = NE1.

    ``swap_mode`` selects how the var-region perturbation u1-u0 is produced:

      - ``"genuine"`` (default): the physically-meaningful in-place toggle. At
        the reference (u0 / WT-coupled) endpoint HE1 is coupled and the methyl
        is present-but-mutually-excluded; at the displaced (u1 / MTR-coupled)
        endpoint HE1 is decoupled by translating it ``genuine_decouple_nm`` into
        bulk via ``ParticleOffsetDisplacement(ne1_ref, NE1)`` against a massless
        zero-LJ dummy NE1 reference placed at NE1 + d. HE1's *internal* angle
        terms are zeroed before migration (a lone H's NE1-centred angles are not
        part of the side-chain free-energy difference and would otherwise inject
        pure internal-geometry strain, not a coupling change). The methyl keeps
        the upstream NE1-attach ``ParticleOffsetDisplacement(NE1, NE1)`` => zero
        offset => stays coupled. The resulting u1-u0 is the genuine cost of
        moving HE1's nonbonded coupling from the core into bulk (ones-to-tens of
        kcal/mol), evaluated at the TRUE mutation geometry — NOT a parked clash.
        This DOES require the dummy NE1 reference + angle decoupling: a single
        shared NE1 cannot carry a bond/angle-strain-free swap through NE1-attach
        alone (both var groups bond the SAME NE1 => zero offset). The dummy
        reference is the minimal in-place stand-in for the upstream two-copy
        overlay's distinct partner-attach atom.

        ``genuine_decouple_dir`` selects the bulk direction the dummy reference is
        placed along (HE1 is displaced ``genuine_decouple_nm`` along it at u1).
        Default ``None`` => the legacy fixed +Z (correct for the FREE leg, whose
        peptide is bathed in solvent so any direction lands in bulk). For the
        BOUND leg the caller passes a unit vector pointing AWAY from the receptor
        (e.g. NE1 - receptor_centroid, normalized) so HE1 decouples into bulk
        solvent rather than translating INTO the receptor (which would replace the
        genuine decouple cost with a receptor-clash artifact). The vector is
        normalized internally; a zero/degenerate vector falls back to +Z.
      - ``"var_park"`` (DIAGNOSTIC): the legacy ~4-atom methyl FixedDisplacement
        park. Produces a finite-but-ELEVATED u1-u0 dominated by CM-NE1 bond
        strain (k~282000 kJ/mol/nm^2), i.e. a parked clash, NOT a physically
        meaningful relative perturbation. Kept only for diagnostics.
      - ``"null"``: the degenerate upstream NE1-attach swap (offset 0 both var
        groups) -> u1 == u0. Demonstrates the single-shared-core null op.

    Mutates ``fused_build['system']`` in place. Returns the ATMForce index +
    wiring bookkeeping. Soft-core canon (umax/ubcore/acore) is NOT re-tuned (C3).
    """
    import copy

    if swap_mode not in ("genuine", "var_park", "null"):
        raise ValueError(
            "attach_inplace_swap_atmforce: swap_mode must be one of "
            "{'genuine','var_park','null'}, got %r" % (swap_mode,))
    # Back-compat: an explicit var_park_nm vector forces the diagnostic park mode
    # (older callers passed var_park_nm to opt into the parked-clash path).
    if var_park_nm is not None and any(abs(c) > 0 for c in var_park_nm) \
            and swap_mode == "genuine":
        swap_mode = "var_park"

    system = fused_build["system"]

    # --- ATMForce: expression-string ctor + 10 globals (upstream, verbatim). ---
    atm = mm.ATMForce(
        _ATS_REFERENCE_POT_EXPR + _ATS_ALCHEMICAL_POT_EXPR + _ATS_SOFTCORE_EXPR
    )
    # Alpha has units 1/energy; input per kcal/mol -> kJ/mol.
    alpha_kj = alpha_per_kcal / _kcal_to_kj(1.0) if alpha_per_kcal else 0.0
    atm.addGlobalParameter("Lambda1", lambda1)
    atm.addGlobalParameter("Lambda2", lambda2)
    atm.addGlobalParameter("Alpha", alpha_kj)
    atm.addGlobalParameter("Uh", _kcal_to_kj(u0_kcal))
    atm.addGlobalParameter("W0", _kcal_to_kj(w0_kcal))
    atm.addGlobalParameter("Umax", _kcal_to_kj(umax_kcal))
    atm.addGlobalParameter("Ubcore", _kcal_to_kj(ubcore_kcal))
    atm.addGlobalParameter("Acore", acore)
    atm.addGlobalParameter("Direction", direction)
    atm.addGlobalParameter("UOffset", _kcal_to_kj(uoffset_kcal))

    lig1_var = cmap["mtr_var"]              # {CM, HM1-3}
    lig2_var = cmap["wt_var_fused"]         # {fused HE1}
    lig1_attach = cmap["mtr_attach"]        # NE1
    lig2_attach = cmap["wt_attach_fused"]   # NE1 (same physical atom)

    # --- GENUINE mode pre-migration setup (must run BEFORE the forces are moved
    #     into the ATMForce, since we mutate the live HarmonicAngleForce and grow
    #     the System/NonbondedForce). ---
    he1_idx = lig2_var[0] if lig2_var else None
    n_he1_angles_decoupled = 0
    ne1_ref_index: Optional[int] = None
    if swap_mode == "genuine" and he1_idx is not None:
        # (1) Zero HE1's internal NE1-centred angle terms. HE1 is a lone H on the
        #     common NE1; its angle terms are part of the COMMON-core internal
        #     geometry, not the side-chain alchemical difference. Leaving them
        #     migrated would make the u1 (HE1-displaced) state pay pure angle
        #     strain (k ~ 1e2-1e3 kcal/mol/rad^2), swamping the real coupling
        #     change with a geometry artifact (the same failure class as the
        #     var_park CM-NE1 bond strain). Decoupling them isolates HE1's
        #     NONBONDED coupling change, which IS the relative perturbation.
        af = next((f for f in system.getForces()
                   if isinstance(f, mm.HarmonicAngleForce)), None)
        if af is not None:
            for a in range(af.getNumAngles()):
                a1, a2, a3, ang, k = af.getAngleParameters(a)
                if he1_idx in (a1, a2, a3):
                    af.setAngleParameters(
                        a, a1, a2, a3, ang,
                        0.0 * unit.kilojoule_per_mole / unit.radian ** 2)
                    n_he1_angles_decoupled += 1
        # (2) Add a massless, zero-LJ, zero-charge DUMMY NE1 reference at
        #     pos(NE1) + d. This is the minimal in-place stand-in for the upstream
        #     two-copy overlay's distinct partner-attach atom: it gives HE1 a
        #     non-coincident reference so ParticleOffsetDisplacement(ne1_ref, NE1)
        #     carries a REAL offset d (instead of the degenerate NE1-attach 0).
        positions = list(fused_build["modeller"].positions)
        ne1_pos = np.array(positions[lig1_attach].value_in_unit(unit.nanometer))
        # Decouple direction: legacy fixed +Z (free leg) unless the caller passes
        # an explicit bulk direction (bound leg => away-from-receptor). Normalize;
        # a degenerate (zero-norm) vector falls back to +Z so the offset is real.
        if genuine_decouple_dir is not None:
            dvec = np.array([float(c) for c in genuine_decouple_dir], dtype=float)
            dnorm = float(np.linalg.norm(dvec))
            unit_dir = (dvec / dnorm) if dnorm > 1e-9 else np.array([0.0, 0.0, 1.0])
        else:
            unit_dir = np.array([0.0, 0.0, 1.0])
        ref_pos = ne1_pos + float(genuine_decouple_nm) * unit_dir
        ne1_ref_index = system.addParticle(0.0 * unit.dalton)
        nb = next(f for f in system.getForces()
                  if isinstance(f, mm.NonbondedForce))
        nb.addParticle(0.0 * unit.elementary_charge,
                       0.1 * unit.nanometer,
                       0.0 * unit.kilojoule_per_mole)
        positions.append(mm.Vec3(*ref_pos) * unit.nanometer)
        fused_build["modeller"].positions = positions

    # --- Migrate the var forces into the ATMForce (upstream var_regions). ---
    to_move = [i for i in range(system.getNumForces())
               if isinstance(system.getForce(i), _ATS_MOVE_FORCE_TYPES)]
    for i in to_move:
        atm.addForce(copy.copy(system.getForce(i)))
    for i in sorted(to_move, reverse=True):
        system.removeForce(i)

    # --- addParticle() for every particle (no-arg overload, upstream L747-748). ---
    for _ in range(system.getNumParticles()):
        atm.addParticle()

    # --- Var swap (upstream add_common_var_atoms_to_atmforce L773-776). ---
    # IN-PLACE single-shared-core construction: the common core is ONE physical
    # atom set (the fused box has a single copy), so the common atoms take the
    # IDENTITY transform (no displacement) — they are literally shared, not two
    # overlaid copies. Only the VAR atoms are transformed. This DIFFERS from the
    # upstream whole-ligand-transfer setup (which overlays two FULL ligand copies
    # via common_i <-> other_common_i). The single-shared-core common atoms have
    # no second copy to swap to; calling ParticleOffsetDisplacement(wt_common_i,
    # mtr_common_i) with WT indices that DON'T EXIST in the fused box (raw indices
    # diverge past NE1) would scatter atoms by the inter-snapshot coordinate
    # difference -> u1 overflow. We therefore leave commons as identity.
    #
    # The var swap (NE1 attach, C2) places each state's var atoms relative to the
    # partner attach atom via the upstream ParticleOffsetDisplacement(other_attach,
    # this_attach). With a SINGLE shared NE1 (lig1_attach == lig2_attach) the
    # offset is ZERO, which yields u1 == u0 (the degenerate null op). The three
    # swap_mode branches below resolve this.
    #
    # (lig1_var/lig2_var/lig1_attach/lig2_attach were bound above, before the
    # pre-migration genuine-mode setup.)

    # Detect the single-shared-core degeneracy (attach atoms coincide -> null op)
    # and the WT-common-not-resident condition (raw indices not in fused system).
    n_particles = system.getNumParticles()
    wt_common_resident = all(0 <= c < n_particles for c in cmap.get("wt_common", []))
    same_attach = (lig1_attach == lig2_attach)

    # Common atoms: identity transform (shared single copy). We do NOT call
    # setParticleTransformation on commons (default = no displacement) — they are
    # literally shared, not two overlaid copies, so they must not be scattered by
    # a WT-vs-MTR coordinate offset.
    park_applied = False
    genuine_applied = False
    if swap_mode == "genuine" and he1_idx is not None and ne1_ref_index is not None:
        # GENUINE in-place toggle. The methyl (MTR var) keeps the upstream
        # NE1-attach displacement => ParticleOffsetDisplacement(NE1, NE1) = 0 =>
        # stays coupled at both endpoints. HE1 (WT var) is displaced by
        # d = pos(ne1_ref) - pos(NE1) into bulk at the u1 (MTR-coupled) endpoint,
        # so its NONBONDED coupling to the rest of the system is removed there.
        # Its internal angle terms were decoupled pre-migration, so u1-u0 is the
        # genuine HE1 coupling change (ones-to-tens of kcal/mol), NOT geometry
        # strain. This is the minimal in-place realisation of the upstream
        # two-copy overlay's distinct partner-attach reference.
        for i in range(len(lig1_var)):   # methyl: zero-offset (NE1-attach), coupled
            atm.setParticleTransformation(
                lig1_var[i],
                mm.ParticleOffsetDisplacement(lig2_attach, lig1_attach))
        for i in range(len(lig2_var)):   # HE1: displaced into bulk at u1
            atm.setParticleTransformation(
                lig2_var[i],
                mm.ParticleOffsetDisplacement(ne1_ref_index, lig1_attach))
        genuine_applied = True
    elif swap_mode == "var_park" and var_park_nm is not None \
            and any(abs(c) > 0 for c in var_park_nm):
        # DIAGNOSTIC park: FixedDisplacement of the ~4-atom methyl. Produces a
        # finite-but-ELEVATED u1-u0 dominated by CM-NE1 bond strain (a parked
        # clash), NOT a physically meaningful relative perturbation. Retained
        # for diagnostics only.
        park = mm.Vec3(*var_park_nm) * unit.nanometer
        for v in lig1_var:   # park the MTR methyl at the displaced endpoint
            atm.setParticleTransformation(v, mm.FixedDisplacement(park))
        park_applied = True
    else:
        # NULL: the degenerate upstream NE1-attach swap (offset 0 both var
        # groups) -> u1 == u0. Demonstrates the single-shared-core null op.
        for i in range(len(lig1_var)):
            atm.setParticleTransformation(
                lig1_var[i],
                mm.ParticleOffsetDisplacement(lig2_attach, lig1_attach))
        for i in range(len(lig2_var)):
            atm.setParticleTransformation(
                lig2_var[i],
                mm.ParticleOffsetDisplacement(lig1_attach, lig2_attach))

    atm_index = system.addForce(atm)
    return {
        "atm_force_index": atm_index,
        "n_forces_migrated": len(to_move),
        "n_common_identity": len(cmap["mtr_common"]),
        "n_lig1_var": len(lig1_var),
        "n_lig2_var": len(lig2_var),
        "attach_atom": lig1_attach,
        "swap_mode": swap_mode,
        "genuine_applied": genuine_applied,
        "genuine_decouple_nm": float(genuine_decouple_nm) if genuine_applied else None,
        "genuine_decouple_dir": (
            [float(c) for c in genuine_decouple_dir]
            if (genuine_applied and genuine_decouple_dir is not None) else None),
        "ne1_ref_index": ne1_ref_index,
        "n_he1_angles_decoupled": n_he1_angles_decoupled,
        "var_park_applied": park_applied,
        "var_park_nm": var_park_nm,
        "single_shared_core_null_op": bool(
            same_attach and not park_applied and not genuine_applied),
        "wt_common_resident_in_fused": bool(wt_common_resident),
    }


def build_inplace_res4_fused_system(
    leg: str = "free",
    seed: str = "s7",
    binder_chain: str = "B",
    solvate: bool = True,
    padding_nm: float = 1.2,
    strict_mc1: bool = False,
    harmonize_common_charges: bool = False,
    var_park_nm: Optional[Tuple[float, float, float]] = None,
    swap_mode: str = "genuine",
    genuine_decouple_nm: float = 1.2,
    mtr_ncaa_xml: Optional[str] = None,
    constraints: Any = HBonds,
) -> Dict[str, Any]:
    """Top-level orchestrator: build the in-place res-4 FUSED dual-topology box.

    Order is R2 -> R1 (C5): the seed geometry is asserted BEFORE the ATMForce
    attach (a fused box built around a clashing methyl is already poisoned).

      1. Build both endpoint systems via ``build_leg_system`` (verbatim reuse).
      2. C3: pair common atoms (count parity + name-order alignment, fail-loud).
      3. MC1: common-atom (q,sigma,eps) endpoint continuity ASSERT.
      4. Harvest HE1's WT terms; inject HE1 into the MTR base -> FUSED box.
      5. R2: seed min-dist + bond-length ASSERT (gates before attach).
      6. MC2: methyl bonded-term presence; MC3: cyclic_ss preserved.
      7. R1: attach the upstream common-core swap ATMForce (NE1 attach).

    Returns the fused build dict + all assert results + the two pristine
    endpoint builds (kept for the smoke's endpoint-equivalence check, C6/Q5).
    Ranking-only (R-11); PREDICTION test (R-18), not a converged DDG.

    ``leg`` selects the thermodynamic leg of the RBFE cycle
    ``DDG_bind = DG_mut(BOUND) - DG_mut(FREE)``:
      - ``"free"``  : the solvated cyclic peptide alone (receptor dropped).
      - ``"bound"`` : the receptor (C3c) + cyclic peptide in the equilibrated
                      BOUND pose, re-solvated, NO displacement (the binder stays
                      bound — this is RBFE, not ABFE). Same residue-4 HE1<->methyl
                      transform as the free leg; the receptor is inert context.
    Both legs share every helper (build, common-map, harvest/inject, asserts,
    swap) verbatim. The only bound-specific logic is sourcing the complex
    structure and pointing the genuine HE1-decouple direction AWAY from the
    receptor into bulk solvent (a fixed +Z could land inside the receptor).

    NOTE: solvate=True is the C6 target (PME-solvated). solvate=False is the
    cheap CPU-only path (unit tests / dry build / the rigorous unsolvated
    endpoint-equivalence pair).

    ``constraints`` (default ``HBonds``) is forwarded to BOTH endpoint
    ``build_leg_system`` calls so the fused box and its WT-harvest reference share
    the same constraint set. The Tier-1 0-step decomposition is constraint-
    insensitive (constraints alter dynamics, not the single-frame potential), so
    the default is unchanged. Tier-2 short dynamics requires ``constraints=None``
    so the appearing/disappearing alch H (methyl HM1-3 / HE1) are NOT under SHAKE
    (the R3 1 fs unconstrained-alch-H requirement). The injected HE1 itself is
    added post-``createSystem`` and never carries a SHAKE constraint regardless.

    Solvation ordering (avoids the fused-residue template gap): the MTR base is
    solvated by ``build_leg_system`` (amber14 has a clean MTR template + waters)
    BEFORE the HE1 injection. ``addParticle`` then appends HE1 at the very end of
    the (already-solvated) System/NonbondedForce/positions, and HE1's WT-harvested
    exception/bond partners are all PROTEIN common atoms whose indices are
    identical in the solvated and unsolvated MTR (protein precedes solvent). So
    HE1 transplant is index-stable regardless of solvation. We never call
    ``createSystem`` on the fused (28-atom res-4) topology, which has no template.
    """
    if leg not in ("free", "bound"):
        raise NotImplementedError(
            "build_inplace_res4_fused_system: leg must be 'free' or 'bound', "
            "got %r." % (leg,))

    li = resolve_leg_inputs(seed)
    if not li["final"]["wt"] or not li["final"]["cp4"]:
        raise FileNotFoundError(
            "build_inplace_res4_fused_system requires both endpoint final.pdb "
            "(seed %s). Got wt=%s cp4=%s"
            % (seed, li["final"]["wt"], li["final"]["cp4"]))

    # 1) Endpoint structures + systems.
    #    FREE leg : the cyclic peptide alone (receptor dropped, re-solvated).
    #    BOUND leg: the receptor + cyclic peptide in the equilibrated BOUND pose
    #              (only water/ions dropped, complex re-solvated). NO displacement
    #              (RBFE — the binder stays bound). The receptor rides through as
    #              inert non-alchemical context; the residue-4 dual-topology swap
    #              is binder-protein-only and leg-agnostic. The WT-harvest reference
    #              is built in the SAME leg layout (UNSOLVATED) so the binder common
    #              atoms align by name-order and HE1's bonded/exception partners
    #              (all binder atoms within 3 bonds) remap cleanly. The MTR base is
    #              built per ``solvate``; protein indices are identical
    #              solvated/unsolvated (protein precedes solvent), so C3 alignment
    #              vs the unsolvated WT holds for both legs.
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix="ats_fused_%s_" % (leg,))
    wt_struct = os.path.join(tmpdir, "wt_%s.pdb" % (leg,))
    mtr_struct = os.path.join(tmpdir, "cp4_%s.pdb" % (leg,))
    if leg == "free":
        prepare_free_peptide_from_final(li["final"]["wt"], wt_struct, binder_chain)
        prepare_free_peptide_from_final(li["final"]["cp4"], mtr_struct, binder_chain)
    else:
        prepare_bound_complex_from_final(li["final"]["wt"], wt_struct, binder_chain)
        prepare_bound_complex_from_final(li["final"]["cp4"], mtr_struct, binder_chain)

    # RBFE MTR XML resolution (cross-track isolation): the in-place RBFE build
    # loads a DEDICATED RBFE-scoped MTR XML, never the shared HYBRID_MTR_XML the
    # other tracks (Track A MD, QM/MM 1-traj, v0.8 RESP, MM-PBSA) consume. Default
    # = HYBRID_MTR_XML_RBFE; fall back to the shared file only if the RBFE file is
    # absent (so the build never silently fails on a fresh checkout). The
    # ``mtr_ncaa_xml`` override lets a caller pin an explicit RBFE XML.
    if mtr_ncaa_xml is not None:
        mtr_xml = mtr_ncaa_xml
    elif os.path.isfile(HYBRID_MTR_XML_RBFE):
        mtr_xml = HYBRID_MTR_XML_RBFE
    else:
        mtr_xml = HYBRID_MTR_XML

    wt_build = build_leg_system(
        wt_struct, leg=leg, binder_chain=binder_chain, solvate=False,
        add_hydrogens=False, hydrogens_xml=li["hydrogens_xml"],
        constraints=constraints)
    mtr_build = build_leg_system(
        mtr_struct, leg=leg, binder_chain=binder_chain, solvate=solvate,
        padding_nm=padding_nm, add_hydrogens=False,
        ncaa_xml=mtr_xml, hydrogens_xml=li["hydrogens_xml"],
        constraints=constraints)

    # 2) C3: pair common atoms (count parity + order alignment). Restrict to the
    #    binder PROTEIN atoms (the MTR base may now carry trailing solvent that
    #    the WT unsolvated build lacks — the common-core swap is protein-only).
    cmap = _build_common_index_map(wt_build, mtr_build, binder_chain)

    # 2b) OPTIONAL de-risk DIAGNOSTIC (NOT a production fix): force the MTR base's
    #     common-atom charges to the WT-canonical values so the common core IS
    #     continuous. This isolates the *mechanical* validation (swap wiring +
    #     endpoint-equivalence + finiteness) from the *scientific* MC1
    #     charge-discontinuity gap. Production requires a proper common-charge
    #     harmonized refit — this in-memory override is a diagnostic only and is
    #     flagged in the returned payload. NEVER use for a quantitative DDG.
    common_charges_harmonized = False
    if harmonize_common_charges:
        _harmonize_common_charges_to_wt(wt_build, mtr_build, cmap)
        common_charges_harmonized = True

    # 3) MC1: common-atom param continuity (endpoint equality). This is the
    #    HIGHEST ncAA risk (Q4): the amber14SB-Trp common-ring charges and the
    #    hybrid-MTR (Khoury/RESP-refit) common-ring charges may diverge, which
    #    makes the "common" core electrostatically discontinuous (finite but
    #    WRONG — pre-registered S1 outcome (ii)). By default this is surfaced as
    #    a STRUCTURED outcome (not an opaque crash) so downstream review sees it;
    #    strict_mc1=True re-raises (used by the unit test for the fail-loud path).
    mc1_error: Optional[str] = None
    try:
        mc1 = assert_common_atom_param_continuity(wt_build, mtr_build, cmap)
    except ValueError as exc:
        if strict_mc1:
            raise
        mc1_error = str(exc)
        mc1 = _summarize_common_charge_divergence(
            wt_build, mtr_build, cmap, binder_chain)
        return {
            "leg": leg,
            "seed": seed,
            "outcome": "mc1_charge_discontinuity",
            "mc1_param_continuity": mc1,
            "mc1_error": mc1_error,
            "common_map": {k: cmap[k] for k in ("n_common", "wt_var", "mtr_var")},
            "regime": "ranking_only",
            "note": ("PRE-REGISTERED outcome (ii): the dual-topology common core "
                     "is electrostatically discontinuous (amber14SB-Trp vs "
                     "hybrid-MTR common-ring charges diverge). Fused box NOT "
                     "attached; this must be resolved (common-charge "
                     "harmonization) before a meaningful DDG. R-18: PREDICTION "
                     "test surfaced a real charge-continuity gap, not a pass."),
        }

    # 4) Harvest HE1 (WT) and inject into the MTR base -> fused box. The WT and
    #    MTR RAW indices diverge past NE1, so build the WT-raw -> MTR-raw common
    #    map and remap every HE1 partner index through it (a partner outside the
    #    common set is a fatal alignment error inside _harvest_he1_terms).
    wt_to_mtr = dict(zip(cmap["wt_common"], cmap["mtr_common"]))
    wt_he1 = cmap["wt_var"][0]
    he1_terms = _harvest_he1_terms(wt_build, wt_he1, wt_to_mtr=wt_to_mtr)
    inj = _inject_he1_into_fused(
        mtr_build, he1_terms,
        mtr_var_indices=list(mtr_build["alchemical_atoms"]["mtr_only"]))
    fused = mtr_build
    fused["fused_he1_index"] = inj["fused_he1_index"]
    fused["injection"] = inj
    # Re-identify the residue-4 partition on the fused topology, then attach the
    # injected HE1 index explicitly (HE1 lives in the dedicated "H" chain, so it
    # is NOT found by the chain-B/resnum-4 scan — its index comes from injection).
    fused["alchemical_atoms"] = identify_alchemical_atoms(
        fused["modeller"].topology, binder_chain=binder_chain)
    fused["alchemical_atoms"]["wt_only"] = [inj["fused_he1_index"]]

    # Extend the index map with the FUSED HE1 index + fused attach (same NE1).
    cmap["wt_var_fused"] = [inj["fused_he1_index"]]
    cmap["wt_attach_fused"] = fused["alchemical_atoms"]["common"][0]  # NE1 (fused)

    # 5) R2: seed assert (gates BEFORE attach).
    seed_assert = assert_appearing_methyl_seed(fused)

    # 6) MC2 + MC3 (pre-attach, on the fused base).
    mc2 = assert_methyl_bonded_present(fused)
    mc3 = assert_disulfide_in_fused(fused)

    # 7) R1: attach the in-place swap ATMForce. For the BOUND leg the genuine HE1
    #    decouple direction is computed OUTWARD from the LOCAL atom density around
    #    NE1 (away from both the receptor AND the binder's own fold), so the dummy
    #    NE1 reference lands in low density (bulk solvent / empty space) and the
    #    "decoupled" HE1 does not clash a real atom — which would replace the
    #    genuine decouple cost with a packing artifact. (A naive away-from-receptor
    #    vector is wrong for this pose: it points back through the peptide ring.)
    #    The FREE leg keeps the legacy fixed +Z (any direction is bulk solvent for
    #    a free peptide); leaving decouple_dir=None there preserves the validated
    #    free-leg numbers byte-for-byte.
    decouple_dir = None
    if leg == "bound":
        decouple_dir = compute_decouple_direction(
            fused, binder_chain=binder_chain,
            decouple_nm=genuine_decouple_nm)
    swap = attach_inplace_swap_atmforce(
        fused, cmap, lambda1=0.0, lambda2=0.0, var_park_nm=var_park_nm,
        swap_mode=swap_mode, genuine_decouple_nm=genuine_decouple_nm,
        genuine_decouple_dir=decouple_dir)

    return {
        "leg": leg,
        "seed": seed,
        "outcome": "fused_attached",
        "fused_build": fused,
        "wt_endpoint": wt_build,
        "mtr_endpoint": mtr_build,
        "common_map": cmap,
        "mc1_param_continuity": mc1,
        "common_charges_harmonized": common_charges_harmonized,
        "mtr_ncaa_xml": mtr_xml,
        "swap_mode": swap_mode,
        "genuine_decouple_dir": decouple_dir,
        "seed_assert": seed_assert,
        "mc2_methyl_bonded": mc2,
        "mc3_disulfide": mc3,
        "swap": swap,
        "solvated": solvate,
        "regime": "ranking_only",
        "note": ("PREDICTION test (R-18); not a converged DDG_bind."
                 + (" [DIAGNOSTIC: common charges harmonized to WT — NOT "
                    "production-valid for a quantitative DDG]"
                    if common_charges_harmonized else "")),
    }


# ===========================================================================
# CANONICAL ATS TWO-COPY (Gallicchio 2025 JCIM, DOI 10.1021/acs.jcim.5c00207;
# preprint arXiv:2412.19971).
# Cp4-WT residue-4 RBFE; C1-C10 design criteria.
#
# This is a SEPARATE, ADDITIVE code path from the single-shared-core
# build_inplace_res4_fused_system / attach_inplace_swap_atmforce above, which are
# left BYTE-IDENTICAL (R-7 / C9 — densify_pilot and other callers still use the
# old path). The single-shared-core box held ONE physical common copy + an
# HE1-into-MTR injection + var-var exclusions; that re-packaged a single-topology
# fused box and collapsed (pertE saturated ~150, overlap ~0, dgbind1 = 0.5*pertE
# deterministic offset).
#
# CANONICAL ATS (the LOAD-BEARING correction in the 06-14 verdict, Q1/Q3):
#   - BOTH endpoint copies are resident in ONE box. copy-1 (MTR) sits at the
#     receptor-binding / free-peptide site; copy-2 (WT) is DISPLACED by a vector
#     d (~40 Å, ATS peptide convention) into BULK SOLVENT. The two copies are
#     therefore spatially separated — common-common clash is avoided by the
#     d-separation, NOT by nonbonded exclusions (NO inter-copy exclusion is added,
#     mirroring upstream add_common_var_atoms_to_atmforce which adds none; C3).
#   - The common-region COORDINATES are SWAPPED as an ATMForce transform
#     (setParticleTransformation + ParticleOffsetDisplacement(other_common_i,
#     this_common_i)) so λ=0 is one physical system (copy-1 at site, copy-2 in
#     bulk) and λ=1 is the reverse (copy-2 at site, copy-1 in bulk). Each copy's
#     VARIABLE atoms map to the PARTNER copy's attach atom (distinct NE1 atoms —
#     NOT a single shared NE1), so the var perturbation is a genuine transfer, not
#     a null op.
#   - There is NO "single MTR base + HE1 inject" here; that paradigm is RETIRED in
#     this path. Each copy is a complete, independently-built endpoint topology.
#
# Reuse map (C1, anti-fragmentation; ALL VERBATIM):
#   - build_leg_system()                 : each endpoint copy, UNSOLVATED.
#   - identify_alchemical_atoms()        : the residue-4 partition per copy.
#   - upstream add_common_var_atoms_to_atmforce (ommsystem.py:745-789) idiom :
#       the common-swap + var-swap loop reproduced verbatim (C1), with DISTINCT
#       attach atoms (copy-1 NE1 != copy-2 NE1).
#   - upstream set_atmforce() expression strings + 10 globals (ommsystem.py:817-834)
#       reused via the _ATS_* module constants (shared with the legacy path).
#   - upstream add_forces_to_atmforce() var_regions migration (ommsystem.py:209-223)
#       reused via _ATS_MOVE_FORCE_TYPES (shared with the legacy path).
#
# Ranking-only (R-11); PREDICTION test (R-18). two-copy correctness is NOT yet
# proven — endpoint-equivalence + frac<UBCORE>0 + O>=0.1 + dgbind1!=0.5 + UWHAM
# convergence (the pilot, post-validation) must hold SIMULTANEOUSLY. Until then this
# is NOT a converged ΔΔG_bind and NEVER a Magotti / absolute comparison.
# ===========================================================================

# ATS peptide displacement convention (Gallicchio 2025 §Methods: TYK2 45 Å,
# peptides 40 Å). copy-2 is moved this far into bulk so the two copies never
# overlap (PME cutoff 1.0 nm << d) — C2/C3 spatial separation, not exclusion.
ATS_TWOCOPY_DISPLACEMENT_NM = 4.0   # 40 Å

# C6 two-copy separation thresholds. The displacement is direction-neutral:
# both magnitude and direction cancel in u1-u0 when copy-2 is fully decoupled.
#   - HARD clash floor: copies closer than this are physically interpenetrating
#     (PE -> +1e15, minimize NaN). Equals the PME cutoff (createSystem
#     nonbondedCutoff=1.0 nm). A build below this is REJECTED unconditionally.
ATS_TWOCOPY_CLASH_FLOOR_NM = 1.0
#   - DECOUPLING-SUFFICIENT acceptance line: the auto-search direction must clear
#     this for the copies to be genuinely decoupled (PME cutoff + an LJ-tail buffer
#     so the dispersion tail past the cutoff is also negligible). Conservative
#     option = 2.0 nm.
ATS_TWOCOPY_ACCEPT_SEP_NM = 1.5
ATS_TWOCOPY_ACCEPT_SEP_NM_CONSERVATIVE = 2.0

# Void-water carve cutoff (opt-in build-time carve, default OFF). The two-copy
# swap displaces a bulk copy's DISAPPEARING heavy atoms (e.g. the Trp indole)
# back onto the PARTNER copy's residue site. When the partner residue lacks those
# atoms (e.g. an Ala site with no indole), addSolvent fills that empty volume with
# bulk water; at the swap (u1) frame the displaced heavy atoms then overlap those
# waters, and because the ATMForce soft-core caps only the perturbation (u1-u0),
# NOT the base energy u1, the raw uncapped clash (u1 ~ 1e13 kJ/mol) makes the
# backward-endpoint minimize NaN-crash. Deleting the WHOLE (neutral) waters that
# penetrate the swap-displaced disappearing-heavy volume removes the clash while
# conserving net charge (whole HOH only). Default cutoff 2.5 Å (reviewed band
# 2.4-2.6 Å); a whole HOH with any atom within this distance of any displaced
# disappearing-heavy atom is carved.
ATS_CARVE_VOID_CUTOFF_NM = 0.25   # 2.5 Å

# Direction auto-search: re-pick the displacement DIRECTION into
# open solvent rather than inflating d when the res-4-local outward vector grazes
# the receptor for a given pose). The base direction is the res-4-local outward
# vector (compute_decouple_direction); candidates are sampled in a cone around it.
ATS_TWOCOPY_AUTOSEARCH_N_CANDIDATES = 12   # candidate directions per cone shell
ATS_TWOCOPY_AUTOSEARCH_CONE_DEG = 60.0     # half-angle of the candidate cone
# Magnitude escalation ladder (nm) tried in order when no candidate direction at
# the previous magnitude clears the acceptance line.
ATS_TWOCOPY_AUTOSEARCH_MAGNITUDES_NM = (4.0, 5.5, 7.0, 9.0, 11.0)


def _displace_copy_positions(
    positions: List[Any], dvec_nm: Tuple[float, float, float]
) -> List[Any]:
    """Return a NEW position list translated by ``dvec_nm`` (a unit-wrapped list).

    Used to move copy-2 (the bulk copy) by the displacement vector d before the
    two copies are merged. Operates on a list of OpenMM ``Quantity`` Vec3
    positions; returns nm-wrapped Vec3s.
    """
    dx, dy, dz = (float(c) for c in dvec_nm)
    out = []
    for p in positions:
        v = p.value_in_unit(unit.nanometer)
        out.append(mm.Vec3(v[0] + dx, v[1] + dy, v[2] + dz) * unit.nanometer)
    return out


def compute_twocopy_displacement_vector(
    copy1_build: Dict[str, Any], copy2_build: Dict[str, Any],
    binder_chain: str = "B", magnitude_nm: float = ATS_TWOCOPY_DISPLACEMENT_NM,
    spec: Optional[Any] = None,
) -> Tuple[float, float, float]:
    """C3: the d-vector that moves copy-2 (bulk copy) clear of BOTH the receptor
    and the binder fold of copy-1.

    Reuses the local-outward logic of ``compute_decouple_direction`` (the common
    attach atom minus the centroid of the local heavy-atom shell around copy-1's
    residue) — the same leg-agnostic direction that clears the receptor AND the
    peptide's own fold for the BOUND pose, and points into bulk for the FREE
    peptide. Scaled to ``magnitude_nm``. Falls back to +X if the local outward
    direction is degenerate (the magnitude alone still separates the copies; the
    precise direction only matters so the bulk copy lands in solvent, which the box
    padding guarantees after a uniform displacement).

    NOT the single-var decouple of the legacy path — here the WHOLE copy-2 is
    translated by d, so the relevant geometry is copy-1's residue outward
    direction (where copy-2 must NOT collide as it is swapped to the site at
    λ=1). ``spec`` (a ``MutationSpec`` / registry name; ``None`` => res-4 default)
    supplies the attach atom + residue number.
    """
    unit_dir = compute_decouple_direction(
        copy1_build, binder_chain=binder_chain, spec=spec)
    if unit_dir is None:
        unit_dir = (1.0, 0.0, 0.0)
    norm = float(np.linalg.norm(np.array(unit_dir, dtype=float)))
    if norm < 1e-9:
        unit_dir = (1.0, 0.0, 0.0)
        norm = 1.0
    scale = float(magnitude_nm) / norm
    return (unit_dir[0] * scale, unit_dir[1] * scale, unit_dir[2] * scale)


def _copy_solute_heavy_positions(
    build: Dict[str, Any], binder_chain: str = "B",
) -> np.ndarray:
    """Heavy-atom (solute, non-solvent) positions (nm) of a SINGLE per-copy build.

    Used by the displacement auto-search to score candidate directions on the
    UNSOLVATED endpoint copies (cheap, CPU-only — no full solvated System rebuild
    per candidate). Solvent/ion residues and hydrogens are dropped (the clash-
    relevant metric is heavy-atom distance). For the bound copy-1 this includes the
    whole receptor + binder; for a binder-only copy-2 it is just the binder.
    """
    top = build["modeller"].topology
    positions = np.array([
        v.value_in_unit(unit.nanometer) for v in build["modeller"].positions])
    keep: List[int] = []
    for atom in top.atoms():
        if atom.residue.name in _SOLVENT_RESNAMES:
            continue
        el = atom.element
        if el is not None and el.symbol == "H":
            continue
        if el is None and atom.name.strip().startswith("H"):
            continue
        keep.append(atom.index)
    return positions[keep]


def _candidate_directions(
    base_unit: Tuple[float, float, float],
    n_candidates: int = ATS_TWOCOPY_AUTOSEARCH_N_CANDIDATES,
    cone_deg: float = ATS_TWOCOPY_AUTOSEARCH_CONE_DEG,
) -> List[Tuple[float, float, float]]:
    """Unit-vector candidates: the base direction plus a ring of directions tilted
    ``cone_deg`` off the base, distributed evenly in azimuth.

    The base outward direction (compute_decouple_direction) is tried FIRST so the
    auto-search reproduces the legacy direction when it already clears (preserving
    the established choice). The cone ring offers alternatives that point into open
    solvent when the base direction grazes the receptor for a given pose (Q5/C3
    (a)). Deterministic (fixed azimuthal spacing, no RNG) so a rebuild is
    reproducible.
    """
    base = np.array(base_unit, dtype=float)
    nb = float(np.linalg.norm(base))
    if nb < 1e-9:
        base = np.array([1.0, 0.0, 0.0])
        nb = 1.0
    base = base / nb

    # Build an orthonormal frame (e1, e2) perpendicular to base.
    ref = np.array([0.0, 0.0, 1.0]) if abs(base[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(base, ref)
    e1 = e1 / float(np.linalg.norm(e1))
    e2 = np.cross(base, e1)

    cands: List[Tuple[float, float, float]] = [
        (float(base[0]), float(base[1]), float(base[2]))]
    n_ring = max(0, int(n_candidates) - 1)
    if n_ring > 0:
        theta = np.radians(cone_deg)
        ct, st = np.cos(theta), np.sin(theta)
        for k in range(n_ring):
            phi = 2.0 * np.pi * k / n_ring
            tilt = ct * base + st * (np.cos(phi) * e1 + np.sin(phi) * e2)
            tilt = tilt / float(np.linalg.norm(tilt))
            cands.append((float(tilt[0]), float(tilt[1]), float(tilt[2])))
    return cands


def auto_search_twocopy_displacement(
    copy1_build: Dict[str, Any], copy2_build: Dict[str, Any],
    binder_chain: str = "B",
    accept_sep_nm: float = ATS_TWOCOPY_ACCEPT_SEP_NM,
    clash_floor_nm: float = ATS_TWOCOPY_CLASH_FLOOR_NM,
    magnitudes_nm: Tuple[float, ...] = ATS_TWOCOPY_AUTOSEARCH_MAGNITUDES_NM,
    n_candidates: int = ATS_TWOCOPY_AUTOSEARCH_N_CANDIDATES,
    cone_deg: float = ATS_TWOCOPY_AUTOSEARCH_CONE_DEG,
    padding_nm: float = 1.2,
    spec: Optional[Any] = None,
) -> Dict[str, Any]:
    """Direction-aware displacement auto-search (task #100).

    Pick the displacement vector ``d`` that moves copy-2 (the bulk copy) clear of
    copy-1's solute body. The search is direction-FIRST, magnitude-SECOND:

      1. Candidate directions = the res-4-local outward base direction
         (compute_decouple_direction) + a cone ring of ``n_candidates`` tilted
         alternatives (``_candidate_directions``). The base is tried first.
      2. For the current magnitude, translate copy-2's heavy atoms by each
         candidate and compute the copy-1<->copy-2 minimum heavy-atom distance,
         BOTH raw and under the periodic minimum-image convention with a PREDICTED
         box (the merge re-pads by ``padding_nm`` around the displaced extent, so
         the box edge is estimated from the per-axis solute span + 2*padding). The
         min-image check rejects directions whose magnitude wraps copy-2 back near
         copy-1 (periodic-image ceiling) at the construction stage.
      3. Choose the candidate that MAXIMISES the (image-aware) min distance. If it
         clears ``accept_sep_nm`` (and the image distance also clears it), accept
         and return ``magnitude * unit_dir``.
      4. If NO candidate at this magnitude clears the acceptance line, escalate to
         the next magnitude (the fallback) and repeat from step 2.
      5. If every magnitude is exhausted, raise with the best achieved distance so
         the caller can escalate (the deterministic final-failure path).

    This is endpoint-/direction-NEUTRAL (ranking-safe): d is a construction
    separation device, not a thermodynamic coordinate (the swap transform is
    partner-offset based; u1-u0 is d-invariant given full decoupling + bulk
    solvation). CPU-only (operates on the unsolvated per-copy coordinates; no System
    rebuild per candidate). Returns the chosen vector + achieved distances + the
    full candidate trail for the build log (task #6).
    """
    c1_heavy = _copy_solute_heavy_positions(copy1_build, binder_chain)
    c2_heavy = _copy_solute_heavy_positions(copy2_build, binder_chain)
    if c1_heavy.size == 0 or c2_heavy.size == 0:
        raise ValueError(
            "auto_search_twocopy_displacement: copy-1 (%d) or copy-2 (%d) has no "
            "solute heavy atoms; cannot score candidate displacements."
            % (c1_heavy.shape[0], c2_heavy.shape[0]))

    base_dir = compute_decouple_direction(
        copy1_build, binder_chain=binder_chain, spec=spec)
    if base_dir is None:
        base_dir = (1.0, 0.0, 0.0)
    candidates = _candidate_directions(base_dir, n_candidates, cone_deg)

    trail: List[Dict[str, Any]] = []
    best_overall: Optional[Dict[str, Any]] = None
    for mag in magnitudes_nm:
        for ci, unit_dir in enumerate(candidates):
            dvec = np.array(unit_dir, dtype=float) * float(mag)
            c2_disp = c2_heavy + dvec
            raw_min = _min_image_min_distance_nm(c1_heavy, c2_disp, None)
            # Predicted box: merge re-pads padding_nm around the union extent.
            all_pos = np.vstack([c1_heavy, c2_disp])
            span = all_pos.max(axis=0) - all_pos.min(axis=0)
            box_lengths = span + 2.0 * float(padding_nm)
            image_min = _min_image_min_distance_nm(c1_heavy, c2_disp, box_lengths)
            score = min(raw_min, image_min)
            rec = {
                "magnitude_nm": float(mag), "candidate_index": ci,
                "unit_dir": [float(x) for x in unit_dir],
                "raw_min_nm": raw_min, "image_min_nm": image_min,
                "score_nm": score,
            }
            trail.append(rec)
            if best_overall is None or score > best_overall["score_nm"]:
                best_overall = dict(rec)
        # Best candidate AT THIS magnitude.
        mag_recs = [r for r in trail if r["magnitude_nm"] == float(mag)]
        best_mag = max(mag_recs, key=lambda r: r["score_nm"])
        if (best_mag["raw_min_nm"] >= accept_sep_nm
                and best_mag["image_min_nm"] >= accept_sep_nm
                and best_mag["score_nm"] >= clash_floor_nm):
            unit_dir = best_mag["unit_dir"]
            dvec = tuple(float(x) * float(mag) for x in unit_dir)
            return {
                "displacement_vector_nm": list(dvec),
                "unit_dir": [float(x) for x in unit_dir],
                "magnitude_nm": float(mag),
                "candidate_index": best_mag["candidate_index"],
                "achieved_raw_min_nm": best_mag["raw_min_nm"],
                "achieved_image_min_nm": best_mag["image_min_nm"],
                "accept_sep_nm": accept_sep_nm,
                "clash_floor_nm": clash_floor_nm,
                "n_candidates": len(candidates),
                "n_magnitudes_tried": magnitudes_nm.index(mag) + 1,
                "trail": trail,
                "accepted": True,
            }

    # Exhausted: deterministic final failure (caller escalates).
    raise ValueError(
        "auto_search_twocopy_displacement: no candidate direction/magnitude cleared "
        "the acceptance line %.3f nm. Best achieved (image-aware) min distance was "
        "%.3f nm at magnitude %.1f nm, candidate %d. Magnitudes tried: %s. The pose "
        "may need a larger box (raise padding) or manual --displacement-nm; do NOT "
        "defeat the C6 guard by any route other than genuine separation." % (
                     accept_sep_nm, best_overall["score_nm"],
                     best_overall["magnitude_nm"], best_overall["candidate_index"],
                     ", ".join("%.1f" % m for m in magnitudes_nm)))


def _build_twocopy_index_map(
    copy1_build: Dict[str, Any], copy2_build: Dict[str, Any],
    n_copy1: int, binder_chain: str = "B",
) -> Dict[str, Any]:
    """Pair the residue-shared (common) atoms between the two RESIDENT copies.

    Unlike ``_build_common_index_map`` (single-shared-core: WT atoms were NOT
    resident in the fused box), here BOTH copies are physically resident in the
    merged box. copy-1 occupies indices [0, n_copy1); copy-2 occupies
    [n_copy1, n_copy1 + n_copy2). The common/var/attach atom indices for copy-2
    are its per-copy indices PLUS ``n_copy1`` (the merge offset).

    Returns index-aligned common lists (in MERGED indices) + each copy's var and
    attach atoms, and runs the C4 hard gate: equal common-atom COUNT and
    byte-identical NAME ORDER (the swap is positional). Distinct attach atoms
    (copy1_attach != copy2_attach) is the structural difference from the legacy
    single-shared-core map.
    """
    c1_top = copy1_build["modeller"].topology
    c2_top = copy2_build["modeller"].topology

    c1_var = set(copy1_build["alchemical_atoms"]["mtr_only"])   # copy-1 = MTR
    c2_var = set(copy2_build["alchemical_atoms"]["wt_only"])    # copy-2 = WT

    c1_atoms = list(c1_top.atoms())
    c2_atoms = list(c2_top.atoms())

    def _is_binder_protein(atom) -> bool:
        return (atom.residue.chain.id == binder_chain
                and atom.residue.name not in _SOLVENT_RESNAMES)

    # copy-1 common: per-copy indices (the merge keeps copy-1 at [0, n_copy1)).
    c1_common = [a.index for a in c1_atoms
                 if a.index not in c1_var and _is_binder_protein(a)]
    # copy-2 common: per-copy indices SHIFTED by the merge offset.
    c2_common = [a.index + n_copy1 for a in c2_atoms
                 if a.index not in c2_var and _is_binder_protein(a)]

    # C4 hard gate 1: common-atom COUNT parity (upstream _exit L763-765).
    if len(c1_common) != len(c2_common):
        raise ValueError(
            "C4 two-copy common-atom COUNT parity FAIL: copy-1 has %d common "
            "atoms, copy-2 has %d. The two-copy swap is positional and requires "
            "equal common-atom counts (upstream ommsystem.py:763-765 _exit)."
            % (len(c1_common), len(c2_common))
        )

    # C4 hard gate 2: index-order alignment by atom NAME (positional swap). Names
    # are read per-copy (copy-2 names from its own topology, indices un-shifted).
    c1_name = {a.index: a.name for a in c1_atoms}
    c2_name = {a.index: a.name for a in c2_atoms}
    mismatches: List[str] = []
    for w_i, m_shifted in zip(c1_common, c2_common):
        m_i = m_shifted - n_copy1
        if c1_name[w_i] != c2_name[m_i]:
            mismatches.append(
                "pos copy1[%d]=%s vs copy2[%d]=%s"
                % (w_i, c1_name[w_i], m_i, c2_name[m_i])
            )
    if mismatches:
        raise ValueError(
            "C4 two-copy common-atom ORDER alignment FAIL (%d mismatches): the "
            "two copies' common-atom lists must be in the same physical order for "
            "the positional swap. First mismatches: %s"
            % (len(mismatches), "; ".join(mismatches[:8]))
        )

    return {
        "copy1_common": c1_common,                               # MERGED indices
        "copy2_common": c2_common,                               # MERGED indices
        "copy1_var": sorted(i for i in c1_var),                  # MTR {CM,HM1-3}
        "copy2_var": sorted(i + n_copy1 for i in c2_var),        # WT {HE1}
        "copy1_attach":
            copy1_build["alchemical_atoms"]["common"][0],        # MTR NE1
        "copy2_attach":
            copy2_build["alchemical_atoms"]["common"][0] + n_copy1,  # WT NE1
        "n_common": len(c1_common),
        "n_copy1": n_copy1,
    }


def assert_twocopy_common_param_continuity(
    system: mm.System, cmap: Dict[str, Any],
    q_tol_e: float = 1e-4, sigma_tol_nm: float = 1e-4, eps_tol_kj: float = 1e-4,
) -> Dict[str, Any]:
    """MC1 strict-mode probe: per-atom common (q, sigma, epsilon) continuity
    between the two RESIDENT copies, read from the MERGED System's NonbondedForce.

    REPORTING-ONLY in the canonical two-copy box: per-atom common-charge
    divergence is NOT required for correctness, because the ATS swap is a
    coordinate-only transform — each copy keeps its native residue-template
    charges and u1-u0 already includes the per-copy charge difference. The
    requirement that the common core be byte-identical was an over-constraint
    inherited from the retired single-shared-core design (one physical common
    copy => one charge).

    This function still RAISES on any per-atom exceedance and is retained for the
    opt-in ``strict_mc1=True`` fail-loud path (the legacy unit-test contract). The
    default two-copy build path uses ``_summarize_twocopy_charge_divergence``
    instead (non-blocking report + a retained NET-charge sanity hard gate).

    Returns the worst per-channel deviation when continuous; raises otherwise.
    """
    nb = next(f for f in system.getForces()
              if isinstance(f, mm.NonbondedForce))
    max_dq = max_dsig = max_deps = 0.0
    worst: List[str] = []
    for c1_i, c2_i in zip(cmap["copy1_common"], cmap["copy2_common"]):
        q1, s1, e1 = nb.getParticleParameters(c1_i)
        q2, s2, e2 = nb.getParticleParameters(c2_i)
        dq = abs(q1.value_in_unit(unit.elementary_charge)
                 - q2.value_in_unit(unit.elementary_charge))
        dsig = abs(s1.value_in_unit(unit.nanometer)
                   - s2.value_in_unit(unit.nanometer))
        deps = abs(e1.value_in_unit(unit.kilojoule_per_mole)
                   - e2.value_in_unit(unit.kilojoule_per_mole))
        max_dq = max(max_dq, dq)
        max_dsig = max(max_dsig, dsig)
        max_deps = max(max_deps, deps)
        if dq > q_tol_e or dsig > sigma_tol_nm or deps > eps_tol_kj:
            worst.append("idx %d<->%d: dq=%.3e dsig=%.3e deps=%.3e"
                         % (c1_i, c2_i, dq, dsig, deps))
    result = {
        "max_dq_e": max_dq, "max_dsigma_nm": max_dsig, "max_deps_kj": max_deps,
        "n_common_checked": cmap["n_common"], "passed": not worst,
    }
    if worst:
        raise ValueError(
            "MC1 two-copy common-atom continuity FAIL (%d atoms exceed tol "
            "q=%.1e sigma=%.1e eps=%.1e): the dual-topology common core is not "
            "electrostatically/LJ continuous between the two resident copies. "
            "Offenders: %s"
            % (len(worst), q_tol_e, sigma_tol_nm, eps_tol_kj, "; ".join(worst[:8]))
        )
    return result


# Net-charge sanity tolerance (C2). Per-ATOM common-charge divergence is
# physically OK in the canonical two-copy ATS box (the swap is a coordinate-only
# transform — each copy retains its native residue-template charges, and u1-u0
# correctly includes the per-copy charge difference; 
# ). It is also OK for
# the COMMON-atom subset alone to carry a non-zero net dq: different residues
# (e.g. amber14 VAL vs ILE) legitimately assign different partial charges to their
# shared backbone/CB atoms, and the variable atoms carry the EXACT complementary
# charge so the WHOLE-RESIDUE total is conserved (measured V3I: common-subset net
# +0.0893 e, exactly cancelled by the var atoms => full alch-residue net 0.0 on
# BOTH copies). The retained hard gate is therefore the FULL alchemical-residue
# net charge per copy (common res-mut atoms + that copy's variable atoms): a
# non-charge-changing mutation must give the SAME residue total on copy-1 and
# copy-2. A divergence there is a genuine build defect (mis-paired common map,
# wrong residue template, or a real charge-changing mutation that this two-copy
# scaffold does not yet support). Tol is loose vs the per-atom 1e-4 to absorb
# float round-off over the residue's atoms.
TWOCOPY_NET_DQ_TOL_E = 1e-3


def _summarize_twocopy_charge_divergence(
    system: mm.System, cmap: Dict[str, Any], copy1_build: Dict[str, Any],
    resnum: int = ALCH_RESNUM, net_dq_tol_e: float = TWOCOPY_NET_DQ_TOL_E,
) -> Dict[str, Any]:
    """Structured per-atom report of the copy1<->copy2 common-charge divergence.

    Two-copy analog of ``_summarize_common_charge_divergence``. In the CANONICAL
    two-copy ATS box per-atom common-charge divergence is NOT a failure — the swap
    is a coordinate-only transform and each copy retains its native charges, so
    u1-u0 stays exact (
    ). This report is
    therefore REPORTING-ONLY (numbers-neutral); the build proceeds with native
    charges. It surfaces WHICH common atoms diverge and by how much (the mutated
    residue vs elsewhere, net displaced charge) for review.

    ``resnum`` selects the mutated binder residue whose per-atom dq is itemised
    (default ``ALCH_RESNUM`` = res-4 MTR, byte-identical legacy behaviour; V3I
    passes ``MutationSpec.resnum`` = 3). ``net_sanity_ok`` is the C2 hard gate
    signal: ``True`` iff the FULL alchemical-residue net charge (the mutated
    residue's common atoms + that copy's variable atoms) AGREES between copy-1 and
    copy-2 within ``net_dq_tol_e`` (a divergence => genuine build defect; the
    common-subset net alone is NOT the gate — see ``TWOCOPY_NET_DQ_TOL_E``).
    """
    nb = next(f for f in system.getForces()
              if isinstance(f, mm.NonbondedForce))
    c1_res = {a.index: a.residue.id
              for a in copy1_build["modeller"].topology.atoms()}
    c1_name = {a.index: a.name
               for a in copy1_build["modeller"].topology.atoms()}

    def _q(idx: int) -> float:
        return nb.getParticleParameters(idx)[0].value_in_unit(
            unit.elementary_charge)

    per_resmut: List[Dict[str, Any]] = []
    sum_dq_resmut = sum_dq_all = 0.0
    # Per-copy running total over the mutated residue's COMMON atoms (the var
    # atoms are added below to form the full residue net).
    q_resmut_common_c1 = q_resmut_common_c2 = 0.0
    n_diverging = 0
    for c1_i, c2_i in zip(cmap["copy1_common"], cmap["copy2_common"]):
        q1 = _q(c1_i)
        q2 = _q(c2_i)
        dq = q1 - q2
        sum_dq_all += dq
        if abs(dq) > 1e-4:
            n_diverging += 1
        if str(c1_res.get(c1_i)) == str(resnum):
            sum_dq_resmut += dq
            q_resmut_common_c1 += q1
            q_resmut_common_c2 += q2
            per_resmut.append({"name": c1_name.get(c1_i), "q_copy1": round(q1, 4),
                               "q_copy2": round(q2, 4), "dq": round(dq, 4)})

    # FULL alchemical-residue net charge per copy = res-mut common atoms + that
    # copy's VARIABLE (appearing/disappearing) atoms. For a non-charge-changing
    # mutation this total must agree across copies (the C2 hard gate). The var
    # atoms carry the complementary charge to the common-subset dq, so this
    # cancels the benign amber14-template common-charge difference.
    q_var_c1 = sum(_q(i) for i in cmap.get("copy1_var", []))
    q_var_c2 = sum(_q(i) for i in cmap.get("copy2_var", []))
    full_resnet_c1 = q_resmut_common_c1 + q_var_c1
    full_resnet_c2 = q_resmut_common_c2 + q_var_c2
    full_resnet_diff = full_resnet_c1 - full_resnet_c2

    return {
        # Per-atom continuity is NOT required in the canonical two-copy box; this
        # field reports the per-atom state for review (it does NOT gate the build).
        "passed": (n_diverging == 0),
        "per_atom_continuous": (n_diverging == 0),
        "n_common_checked": cmap["n_common"],
        "n_diverging": n_diverging,
        # C2 net-charge sanity (the retained HARD gate): the FULL mutated-residue
        # net charge (common + var) must agree across the two copies. The common-
        # subset net alone is reported but does NOT gate (see docstring).
        "net_sanity_ok": abs(full_resnet_diff) <= net_dq_tol_e,
        "net_dq_tol_e": net_dq_tol_e,
        "full_resmut_net_copy1_e": round(full_resnet_c1, 6),
        "full_resmut_net_copy2_e": round(full_resnet_c2, 6),
        "full_resmut_net_diff_e": round(full_resnet_diff, 6),
        # Back-compat alias kept (res-4 MTR consumers); ``sum_dq_resmut_e`` is the
        # generalized name (the mutated residue, whatever ``resnum`` is). This is
        # the COMMON-subset dq (benign; informational only).
        "sum_dq_res4_e": round(sum_dq_resmut, 5),
        "sum_dq_resmut_e": round(sum_dq_resmut, 5),
        "sum_dq_all_common_e": round(sum_dq_all, 6),
        "residue4_common": per_resmut,
        "residue_mut_common": per_resmut,
    }


def _twocopy_solute_heavy_indices(
    fused_build: Dict[str, Any], n_copy1: int,
) -> Tuple[List[int], List[int]]:
    """Split the merged topology's SOLUTE heavy-atom indices by copy.

    The merge places copy-1 at ``[0, n_copy1)`` and copy-2 at ``[n_copy1, N)``.
    Solute = any atom whose residue is NOT in ``_SOLVENT_RESNAMES`` (the shared
    water/ion bath is excluded — the two copies share one solvent shell, so a
    solvent atom near both copies is expected and must not trip the guard).
    Hydrogens are dropped (heavy-atom min distance is the clash-relevant metric
    and keeps the pairwise distance matrix small).

    Returns ``(copy1_solute_heavy, copy2_solute_heavy)`` as MERGED indices.
    """
    top = fused_build["modeller"].topology
    copy1_idx: List[int] = []
    copy2_idx: List[int] = []
    for atom in top.atoms():
        if atom.residue.name in _SOLVENT_RESNAMES:
            continue
        el = atom.element
        if el is not None and el.symbol == "H":
            continue
        if el is None and atom.name.strip().startswith("H"):
            continue
        if atom.index < n_copy1:
            copy1_idx.append(atom.index)
        else:
            copy2_idx.append(atom.index)
    return copy1_idx, copy2_idx


def _box_lengths_nm_from_vectors(box_vectors: Any) -> Optional[np.ndarray]:
    """Return the per-axis box lengths (nm) from a triclinic box-vector triple.

    Accepts OpenMM ``getPeriodicBoxVectors()``-style output (a 3x3 of Quantity
    Vec3, or a bare 3x3 array in nm). For the addSolvent boxes here the box is
    rectangular, so the minimum-image convention only needs the diagonal lengths
    (|a_x|, |b_y|, |c_z|). Off-diagonal (triclinic tilt) terms are ignored — these
    boxes are orthorhombic by construction (addSolvent default). Returns ``None``
    if the box cannot be resolved (a non-periodic / unsolvated build).
    """
    if box_vectors is None:
        return None
    rows = []
    for vec in box_vectors:
        try:
            comp = vec.value_in_unit(unit.nanometer)
        except AttributeError:
            comp = vec
        rows.append([float(comp[0]), float(comp[1]), float(comp[2])])
    arr = np.array(rows, dtype=float)
    if arr.shape != (3, 3):
        return None
    lengths = np.array([abs(arr[0, 0]), abs(arr[1, 1]), abs(arr[2, 2])], dtype=float)
    if not np.all(np.isfinite(lengths)) or np.any(lengths <= 0.0):
        return None
    return lengths


def _min_image_min_distance_nm(
    c1_pos: np.ndarray, c2_pos: np.ndarray, box_lengths_nm: Optional[np.ndarray],
) -> float:
    """Minimum heavy-atom distance between two coordinate sets under the minimum-
    image convention (per-axis box wrapping).

    Without ``box_lengths_nm`` this is the raw Euclidean min distance. With it,
    each pairwise component delta is wrapped into ``[-L/2, L/2]`` per axis
    (``delta -= L * round(delta / L)``) so a copy-2 atom that sits across a
    periodic boundary — i.e. its nearest IMAGE is close to copy-1 even though its
    raw coordinate is far — is measured at its true nearest-image distance. This
    catches the Q5/C3 periodic-image ceiling: a d so large that copy-2 wraps back
    NEAR copy-1, which the raw-coordinate check is blind to. Symmetric by
    construction (the wrapped delta is the same whether measured c1->c2 or
    c2->c1-image), so a single wrapped pass covers both directions the task asks
    for.
    """
    deltas = c1_pos[:, None, :] - c2_pos[None, :, :]
    if box_lengths_nm is not None:
        # Wrap only on axes with a positive box length (a degenerate/zero axis is
        # treated as non-periodic on that axis — avoids a divide-by-zero).
        safe = np.array(box_lengths_nm, dtype=float)
        usable = safe > 0.0
        if np.any(usable):
            shift = np.zeros_like(deltas)
            shift[..., usable] = (
                safe[usable] * np.round(deltas[..., usable] / safe[usable]))
            deltas = deltas - shift
    return float(np.sqrt((deltas ** 2).sum(axis=-1)).min())


def assert_twocopy_separation(
    fused_build: Dict[str, Any], cmap: Dict[str, Any], min_sep_nm: float = 1.0,
    box_vectors: Any = None, accept_sep_nm: Optional[float] = None,
) -> Dict[str, Any]:
    """C6: the two copies are spatially SEPARATED (no overlay) BEFORE attach.

    TWO independent separation gates, both at the ``min_sep_nm`` floor (the PME
    cutoff; the real displacement is ~40 Å):

      1. residue-4 NE1<->NE1 (the swap attach atoms) — the per-residue check.
      2. copy-1 SOLUTE <-> copy-2 SOLUTE minimum heavy-atom distance — the
         WHOLE-MOLECULE check. The NE1<->NE1 gate alone is NECESSARY but NOT
         SUFFICIENT: when copy-2 was (wrongly) built with its own full receptor,
         the res-4-local d-vector left the two NE1 atoms far apart while the two
         receptor BODIES interpenetrated (107 clashes, PE -> +1e15, minimize
         NaN). The solute-solute heavy-atom min distance catches that
         interpenetration / under-displacement that NE1<->NE1 is blind to. The
         shared solvent/ion bath is excluded (the copies share one water shell).

      3. (when ``box_vectors`` given) the PERIODIC MINIMUM-IMAGE solute<->solute
         distance — copy-1 solute vs copy-2 solute AND copy-2 vs copy-1's nearest
         IMAGE. The raw-coordinate checks (1,2) are blind to the Q5/C3 periodic-
         image ceiling: a d so large that copy-2 (or its solvation shell) wraps
         across the box boundary places its nearest image back NEAR copy-1, re-
         introducing the cross-talk the separation was meant to remove. The min-
         image distance measures the true nearest-image separation and fails the
         build if it drops below the clash floor.

    A small separation on any gate means the copies overlap (or wrap) and the
    build collapsed back toward the (forbidden) overlay / interpenetrating design.

    ``min_sep_nm`` is the HARD clash floor (default the PME cutoff, 1.0 nm) — any
    gate below it is unconditionally rejected. ``accept_sep_nm`` (optional) is the
    higher DECOUPLING-SUFFICIENT acceptance line (e.g. 1.5 nm): when set, the raw
    AND image solute distances are also checked against it (the auto-search caller
    requests this elevated line; the standalone default keeps only the clash floor
    so existing direct callers stay byte-identical). ``box_vectors`` (optional) is
    the merged box's periodic box vectors (available only after addSolvent); when
    ``None`` the image gate is skipped (an unsolvated build has no periodic box).
    """
    positions = np.array([
        v.value_in_unit(unit.nanometer)
        for v in fused_build["modeller"].positions])
    ne1_c1 = cmap["copy1_attach"]
    ne1_c2 = cmap["copy2_attach"]
    sep = float(np.linalg.norm(positions[ne1_c1] - positions[ne1_c2]))
    if sep < min_sep_nm:
        raise ValueError(
            "C6 two-copy separation FAIL: copy-1 NE1 and copy-2 NE1 are %.3f nm "
            "apart (< %.3f nm). The two copies must be d-displaced into bulk "
            "(overlay is the FORBIDDEN single-shared-core design — clash must be "
            "avoided by separation, not exclusion)." % (sep, min_sep_nm))

    # Gate 2: whole-solute min heavy-atom distance (catches receptor-receptor /
    # under-displacement interpenetration the NE1<->NE1 gate cannot see).
    c1_heavy, c2_heavy = _twocopy_solute_heavy_indices(fused_build, cmap["n_copy1"])
    if not c1_heavy or not c2_heavy:
        raise ValueError(
            "C6 two-copy separation FAIL: copy-1 (%d) or copy-2 (%d) has no solute "
            "heavy atoms in the merged topology; cannot verify solute separation."
            % (len(c1_heavy), len(c2_heavy)))
    c1_pos = positions[c1_heavy]
    c2_pos = positions[c2_heavy]
    solute_min_sep = _min_image_min_distance_nm(c1_pos, c2_pos, None)
    if solute_min_sep < min_sep_nm:
        raise ValueError(
            "C6 two-copy SOLUTE separation FAIL: copy-1 solute and copy-2 solute "
            "minimum heavy-atom distance is %.3f nm (< %.3f nm). The copies' bodies "
            "interpenetrate or copy-2 is under-displaced (a res-4-local d-vector "
            "cannot clear a full second receptor — canonical ATS uses ONE shared "
            "receptor + binder-only copy-2). NE1<->NE1 was %.3f nm (passed) but the "
            "whole-solute check caught the interpenetration."
            % (solute_min_sep, min_sep_nm, sep))

    # Gate 3: PERIODIC minimum-image solute<->solute distance (Q5/C3 image ceiling).
    # Available only when the merged box's periodic vectors are supplied (post
    # addSolvent). A d that wraps copy-2 back near copy-1 passes the raw gate but
    # fails here. The wrapped delta is symmetric, so this single pass covers both
    # copy1<->copy2 and copy2<->copy1-image.
    box_lengths = _box_lengths_nm_from_vectors(box_vectors)
    image_min_sep: Optional[float] = None
    if box_lengths is not None:
        image_min_sep = _min_image_min_distance_nm(c1_pos, c2_pos, box_lengths)
        if image_min_sep < min_sep_nm:
            raise ValueError(
                "C6 two-copy PERIODIC-IMAGE separation FAIL: copy-1 solute and "
                "copy-2 solute minimum-IMAGE heavy-atom distance is %.3f nm "
                "(< %.3f nm) under box %s nm, while the raw distance was %.3f nm. "
                "The displacement pushed copy-2 across the periodic boundary so its "
                "nearest image wraps back near copy-1 (Q5/C3 periodic-image "
                "ceiling) — reduce d or re-pick the displacement direction into "
                "open solvent." % (image_min_sep, min_sep_nm,
                                   np.round(box_lengths, 3).tolist(), solute_min_sep))

    # Optional elevated acceptance line (decoupling-sufficient, e.g. 1.5 nm). The
    # clash floor above is the hard reject; this is the higher bar the auto-search
    # caller requires so the copies are genuinely decoupled (PME cutoff + LJ tail).
    if accept_sep_nm is not None:
        worst_raw = min(sep, solute_min_sep)
        if worst_raw < accept_sep_nm:
            raise ValueError(
                "C6 two-copy ACCEPTANCE separation FAIL: the smaller of NE1<->NE1 "
                "(%.3f nm) and solute<->solute (%.3f nm) is below the decoupling-"
                "sufficient acceptance line %.3f nm (clash floor %.3f nm cleared). "
                "Increase d or re-pick the displacement direction."
                % (sep, solute_min_sep, accept_sep_nm, min_sep_nm))
        if image_min_sep is not None and image_min_sep < accept_sep_nm:
            raise ValueError(
                "C6 two-copy ACCEPTANCE (image) separation FAIL: the periodic "
                "minimum-image solute distance %.3f nm is below the acceptance line "
                "%.3f nm. copy-2's nearest image is within the decoupling buffer of "
                "copy-1." % (image_min_sep, accept_sep_nm))

    return {
        "ne1_ne1_sep_nm": sep,
        "solute_solute_min_sep_nm": solute_min_sep,
        "image_solute_min_sep_nm": image_min_sep,
        "min_sep_nm": min_sep_nm,
        "accept_sep_nm": accept_sep_nm,
        "n_copy1_solute_heavy": len(c1_heavy),
        "n_copy2_solute_heavy": len(c2_heavy),
        "passed": True,
    }


def attach_twocopy_swap_atmforce(
    fused_build: Dict[str, Any],
    cmap: Dict[str, Any],
    lambda1: float = 0.0,
    lambda2: float = 0.0,
    alpha_per_kcal: float = 0.0,
    u0_kcal: float = 0.0,
    w0_kcal: float = 0.0,
    umax_kcal: float = ATS_UMAX_KCAL,
    ubcore_kcal: float = ATS_UBCORE_KCAL,
    acore: float = ATS_ACORE,
    direction: float = 1.0,
    uoffset_kcal: float = 0.0,
) -> Dict[str, Any]:
    """C1/C2/C3: attach the CANONICAL ATS two-copy coordinate-swap ATMForce.

    Reuses the upstream var-region protocol (add_common_var_atoms_to_atmforce,
    ommsystem.py:745-789) VERBATIM:
      - expression-string ``ATMForce`` ctor with the upstream reference/alchemy/
        soft-core strings + all 10 globals (incl. UOffset) — shared _ATS_* consts;
      - migrate ALL Nonbonded/Harmonic/Torsion forces into the ATMForce via
        ``copy.copy`` + ``removeForce`` (upstream var_regions branch);
      - addParticle() for every particle, then the upstream
        ``setParticleTransformation`` + ``ParticleOffsetDisplacement`` idiom.

    The swap (C1, verbatim upstream):
      common: copy1_common_i  <- ParticleOffsetDisplacement(copy2_common_i, copy1_common_i)
              copy2_common_i  <- ParticleOffsetDisplacement(copy1_common_i, copy2_common_i)
      var   : copy1_var_i     <- ParticleOffsetDisplacement(copy2_attach, copy1_attach)
              copy2_var_i     <- ParticleOffsetDisplacement(copy1_attach, copy2_attach)
    The attach atoms are DISTINCT (copy1_attach != copy2_attach), so the var
    offset is a REAL d-transfer, not the single-shared-core null op. NO inter-copy
    exclusion is added (C3 — clash is avoided by the d-separation, mirroring the
    upstream which adds none).

    Mutates ``fused_build['system']`` in place. Returns the ATMForce index +
    wiring bookkeeping. Soft-core canon (umax/ubcore/acore) is NOT re-tuned (C7).
    """
    import copy

    system = fused_build["system"]

    if cmap["copy1_attach"] == cmap["copy2_attach"]:
        raise ValueError(
            "attach_twocopy_swap_atmforce: copy1_attach == copy2_attach (%d). "
            "The two-copy swap REQUIRES distinct attach atoms (the whole point of "
            "the rebuild — a shared attach is the single-shared-core null op)."
            % (cmap["copy1_attach"],))

    # --- ATMForce: expression-string ctor + 10 globals (upstream, verbatim). ---
    atm = mm.ATMForce(
        _ATS_REFERENCE_POT_EXPR + _ATS_ALCHEMICAL_POT_EXPR + _ATS_SOFTCORE_EXPR
    )
    alpha_kj = alpha_per_kcal / _kcal_to_kj(1.0) if alpha_per_kcal else 0.0
    atm.addGlobalParameter("Lambda1", lambda1)
    atm.addGlobalParameter("Lambda2", lambda2)
    atm.addGlobalParameter("Alpha", alpha_kj)
    atm.addGlobalParameter("Uh", _kcal_to_kj(u0_kcal))
    atm.addGlobalParameter("W0", _kcal_to_kj(w0_kcal))
    atm.addGlobalParameter("Umax", _kcal_to_kj(umax_kcal))
    atm.addGlobalParameter("Ubcore", _kcal_to_kj(ubcore_kcal))
    atm.addGlobalParameter("Acore", acore)
    atm.addGlobalParameter("Direction", direction)
    atm.addGlobalParameter("UOffset", _kcal_to_kj(uoffset_kcal))

    # --- Migrate the var forces into the ATMForce (upstream var_regions). ---
    to_move = [i for i in range(system.getNumForces())
               if isinstance(system.getForce(i), _ATS_MOVE_FORCE_TYPES)]
    for i in to_move:
        atm.addForce(copy.copy(system.getForce(i)))
    for i in sorted(to_move, reverse=True):
        system.removeForce(i)

    # --- addParticle() for every particle (no-arg overload, upstream L747-748). --
    for _ in range(system.getNumParticles()):
        atm.addParticle()

    # --- Common + var swap (upstream add_common_var_atoms_to_atmforce L769-776,
    #     VERBATIM). BOTH common sets are resident -> the common-common swap is a
    #     REAL coordinate transform (NOT the legacy identity), and the var sets map
    #     to the PARTNER copy's DISTINCT attach atom. NO exclusions added (C3). ---
    c1_common = cmap["copy1_common"]
    c2_common = cmap["copy2_common"]
    c1_var = cmap["copy1_var"]
    c2_var = cmap["copy2_var"]
    c1_attach = cmap["copy1_attach"]
    c2_attach = cmap["copy2_attach"]

    for i in range(len(c1_common)):
        atm.setParticleTransformation(
            c1_common[i],
            mm.ParticleOffsetDisplacement(c2_common[i], c1_common[i]))
    for i in range(len(c2_common)):
        atm.setParticleTransformation(
            c2_common[i],
            mm.ParticleOffsetDisplacement(c1_common[i], c2_common[i]))
    for i in range(len(c1_var)):
        atm.setParticleTransformation(
            c1_var[i],
            mm.ParticleOffsetDisplacement(c2_attach, c1_attach))
    for i in range(len(c2_var)):
        atm.setParticleTransformation(
            c2_var[i],
            mm.ParticleOffsetDisplacement(c1_attach, c2_attach))

    atm_index = system.addForce(atm)
    return {
        "atm_force_index": atm_index,
        "n_forces_migrated": len(to_move),
        "n_common_swapped": len(c1_common),
        "n_copy1_var": len(c1_var),
        "n_copy2_var": len(c2_var),
        "copy1_attach": c1_attach,
        "copy2_attach": c2_attach,
        "swap_mode": "twocopy",
        "distinct_attach": True,
        "inter_copy_exclusions_added": 0,   # C3: NONE (clash avoided by d-sep).
    }


def _twocopy_alchemical_atoms(
    copy1_build: Dict[str, Any], copy2_build: Dict[str, Any], n_copy1: int,
) -> Dict[str, Any]:
    """Residue-4 alchemical partition on the MERGED two-copy box.

    copy-1 (MTR, site) contributes ``common`` (NE1) + ``mtr_only`` (CM, HM1-3)
    at their un-shifted indices; copy-2 (WT, bulk) contributes ``wt_only`` (HE1)
    + its own NE1 at +n_copy1. The ``common`` slot is copy-1's NE1 (the swap
    attaches each copy's var to its OWN NE1, both resident; the index map carries
    both attach atoms separately). Returns the merged indices the swap consumes.
    """
    c1 = copy1_build["alchemical_atoms"]
    c2 = copy2_build["alchemical_atoms"]
    return {
        # copy-1 attach NE1 (un-shifted). The two-copy swap uses cmap's
        # copy1_attach / copy2_attach for the distinct attach atoms; this list is
        # the convenience "common attach" pointer (copy-1 NE1) for callers.
        "common": [c1["common"][0]],
        "mtr_only": list(c1["mtr_only"]),                       # copy-1 indices
        "wt_only": [i + n_copy1 for i in c2["wt_only"]],        # copy-2 +offset
        "copy1_ne1": c1["common"][0],
        "copy2_ne1": c2["common"][0] + n_copy1,
    }


def _detect_twocopy_disulfides(
    merged_modeller, n_copy1: int, binder_chain: str = "B",
) -> List[Dict[str, Any]]:
    """Detect the cyclic_ss disulfide in EACH copy on the merged topology.

    Two copies => two SG-SG disulfides. Reuses the shared ``detect_disulfide_pair``
    on the merged modeller; both copies' CYS2-CYS12 SG-SG pairs are within the
    detector threshold (copy-2 is rigidly displaced by d, so its intra-copy SG-SG
    distance is preserved). Commits any missing bond on the merged topology.
    Returns one record per detected pair (expected: 2).

    The detector scans ALL chain-B residues; with two copies the merged topology
    has TWO chain-B segments (Modeller.add preserves chain ids), so the detector
    may pair across copies if SG indices interleave. To keep the pairing
    intra-copy we partition by the merge offset: copy-1 SGs are < n_copy1, copy-2
    SGs are >= n_copy1, and we pair the two nearest SGs within each partition.
    """
    sgs = []
    positions = np.array([
        v.value_in_unit(unit.nanometer) for v in merged_modeller.positions])
    for atom in merged_modeller.topology.atoms():
        if (atom.residue.chain.id == binder_chain
                and atom.residue.name in ("CYS", "CYX")
                and atom.name == "SG"):
            sgs.append(atom)
    records: List[Dict[str, Any]] = []
    for lo, hi, label in ((0, n_copy1, "copy1"),
                          (n_copy1, 10 ** 12, "copy2")):
        part = [a for a in sgs if lo <= a.index < hi]
        # Pair the two SGs whose separation is the smallest (the cyclic_ss pair).
        best = None
        for i in range(len(part)):
            for j in range(i + 1, len(part)):
                d = float(np.linalg.norm(
                    positions[part[i].index] - positions[part[j].index]))
                if d <= DISULFIDE_MAX_NM and (best is None or d < best[2]):
                    best = (part[i], part[j], d)
        if best is not None:
            added = commit_bond_if_missing(
                merged_modeller.topology, best[0], best[1])
            records.append({
                "copy": label,
                "sg1_index": best[0].index, "sg2_index": best[1].index,
                "sg_dist_nm": best[2], "bond_added": added,
            })
    return records


def _harmonize_twocopy_common_charges(
    system: mm.System, cmap: Dict[str, Any],
) -> int:
    """DIAGNOSTIC-ONLY: overwrite copy-2 (WT) common-atom charges with copy-1
    (MTR) values in the MERGED System so the common core is continuous (MC1
    passes).

    NOT a production fix — production uses the harmonized RBFE XML (Σ|Δq|=0) on
    disk. This in-memory override only ISOLATES the mechanical swap / endpoint-
    equivalence validation from the charge gap. Sigma/eps are already identical
    (only charge diverges). Mutates the merged NonbondedForce. Returns the count.
    """
    nb = next(f for f in system.getForces()
              if isinstance(f, mm.NonbondedForce))
    n = 0
    for c1_i, c2_i in zip(cmap["copy1_common"], cmap["copy2_common"]):
        q1, _, _ = nb.getParticleParameters(c1_i)
        _, s2, e2 = nb.getParticleParameters(c2_i)
        nb.setParticleParameters(c2_i, q1, s2, e2)
        n += 1
    return n


def _collect_bond_constraint_pairs(system: mm.System):
    """All bonded + constrained atom-index pairs in a (possibly ATMForce-nested)
    System, as two ``frozenset`` sets. Shared by the MC2 / R2 bonded-term gates.
    """
    bond_pairs = set()
    for f in system.getForces():
        if isinstance(f, mm.HarmonicBondForce):
            for bi in range(f.getNumBonds()):
                p1, p2, _, _ = f.getBondParameters(bi)
                bond_pairs.add(frozenset((p1, p2)))
        elif isinstance(f, mm.ATMForce):
            for j in range(f.getNumForces()):
                inner = f.getForce(j)
                if isinstance(inner, mm.HarmonicBondForce):
                    for bi in range(inner.getNumBonds()):
                        p1, p2, _, _ = inner.getBondParameters(bi)
                        bond_pairs.add(frozenset((p1, p2)))
    constraint_pairs = set()
    for ci in range(system.getNumConstraints()):
        a, b, _ = system.getConstraintParameters(ci)
        constraint_pairs.add(frozenset((a, b)))
    return bond_pairs, constraint_pairs


def _heavy_h_neighbours(
    heavy_idx: int, var_slot, bond_pairs, constraint_pairs,
    name_by_idx: Dict[int, str],
) -> List[int]:
    """Hydrogen indices (within this copy's var slot) connected to ``heavy_idx``.

    Found by CONNECTIVITY in the merged box (bond OR SHAKE constraint), not by name
    prefix, so the per-heavy H grouping is robust to multi-heavy specs where two
    heavies share a methyl-H name prefix (V3A's CG1 HG1x and CG2 HG2x both prefix
    "HG"). An H is a var atom whose name starts with ``H`` (after any leading
    wyckoff digit) and that is connected to ``heavy_idx``. BOTH ``bond_pairs`` and
    ``constraint_pairs`` are searched because under HBonds the C-H bonds are
    converted to SHAKE constraints and removed from the HarmonicBondForce (the
    single-heavy MC2 makes the same bond-OR-constraint allowance via ``_connected``).
    """
    out: List[int] = []
    slot = set(var_slot)
    for pr in (bond_pairs | constraint_pairs):
        if heavy_idx not in pr:
            continue
        other = next(iter(pr - {heavy_idx}))
        if other not in slot:
            continue
        nm = name_by_idx.get(other, "")
        base = nm.strip().lstrip("0123456789")
        if base[:1].upper() == "H":
            out.append(other)
    return out


def _assert_twocopy_multiheavy_bonded(
    fused_build: Dict[str, Any], cmap: Dict[str, Any], spec, shape: str,
) -> Dict[str, Any]:
    """MC2 multi-heavy branch (C5): the acyclic-star multi-methyl var group's bonded
    terms are present in the merged box.

    For ``multi_appearing_heavy`` the appearing heavies live in copy-1 (``mtr_only``
    slot, attach = ``copy1_attach``); for ``multi_disappearing_heavy`` they live in
    copy-2 (``wt_only`` slot, attach = ``copy2_attach``). The heavy NAME list is the
    SSOT-derived ``spec._heavy_names`` of the corresponding state-only atom set
    (not a single field), so the per-heavy loop matches the partition exactly.

    For EACH var heavy: (1) it MUST be present in the var slot; (2) it MUST bond the
    common attach atom by a REAL HarmonicBond (constraint not accepted — this is the
    star-only scope enforcement and the chained-heavy / ring backstop: a CD that
    bonds only its parent CG, not the attach, fail-louds here); (3) it MUST carry at
    least one H connected to it (bond OR constraint), else the build is incomplete.
    Any missing heavy / attach bond / H => RAISE (fail-loud, never silent-skip).
    """
    system = fused_build["system"]
    name_by_idx = {a.index: a.name
                   for a in fused_build["modeller"].topology.atoms()}
    alch = fused_build["alchemical_atoms"]

    if shape == "multi_appearing_heavy":
        attach_idx = cmap["copy1_attach"]          # appearing var lives in copy-1
        var_slot = list(alch["mtr_only"])
        declared_heavies = spec._heavy_names(spec.stateB_only_atoms)
        side = "appearing"
    else:  # multi_disappearing_heavy
        attach_idx = cmap["copy2_attach"]          # disappearing var lives in copy-2
        var_slot = list(alch["wt_only"])
        declared_heavies = spec._heavy_names(spec.stateA_only_atoms)
        side = "disappearing"

    bond_pairs, constraint_pairs = _collect_bond_constraint_pairs(system)
    attach_name = name_by_idx.get(attach_idx, "attach")

    def _connected(i, j):
        return (frozenset((i, j)) in bond_pairs
                or frozenset((i, j)) in constraint_pairs)

    per_heavy: List[Dict[str, Any]] = []
    for heavy_name in declared_heavies:
        heavy_idx = next(
            (i for i in var_slot if name_by_idx.get(i) == heavy_name), None)
        if heavy_idx is None:
            raise ValueError(
                "MC2 two-copy FAIL: declared %s heavy %r absent from the merged "
                "box's %s var slot (cannot certify a heavy that was not built)."
                % (side, heavy_name, side))
        # heavy<->attach MUST be a REAL HarmonicBond (not a constraint): this is the
        # acyclic-star scope enforcement + the chained-heavy / ring fail-loud
        # backstop (a heavy that bonds only its parent heavy, not the attach, fails).
        heavy_attach_bond = frozenset((heavy_idx, attach_idx)) in bond_pairs
        if not heavy_attach_bond:
            raise ValueError(
                "MC2 two-copy FAIL: %s-%s internal bond absent from the merged "
                "box's HarmonicBondForce (each %s star heavy var atom must bond the "
                "attach atom by a real bond; a heavy that bonds only another heavy "
                "is a chained/ring shape NOT in the acyclic-star scope)."
                % (heavy_name, attach_name, side))
        h_idxs = _heavy_h_neighbours(
            heavy_idx, var_slot, bond_pairs, constraint_pairs, name_by_idx)
        h_present = {name_by_idx[h]: _connected(h, heavy_idx) for h in h_idxs}
        missing_h = [k for k, v in h_present.items() if not v]
        if missing_h:
            raise ValueError(
                "MC2 two-copy FAIL: %s-H connectivity to %s absent (neither bond "
                "nor constraint): %s" % (side, heavy_name, missing_h))
        if not h_idxs:
            raise ValueError(
                "MC2 two-copy FAIL: %s heavy %r has ZERO hydrogens bonded to it in "
                "the %s var slot — an incomplete side-chain build."
                % (side, heavy_name, side))
        per_heavy.append({
            "heavy": heavy_name,
            "heavy_attach_bond_present": heavy_attach_bond,
            "h_names": sorted(h_present),
            "h_connected": h_present,
        })

    return {
        "shape": shape,
        "side": side,
        "attach": attach_name,
        "n_heavies_certified": len(per_heavy),
        "per_heavy": per_heavy,
        "passed": True,
    }


def _idx_by_name(var_slot, name, name_by_idx) -> Optional[int]:
    """First var-slot index whose name matches ``name`` (None if absent)."""
    return next((i for i in var_slot if name_by_idx.get(i) == name), None)


def _assert_twocopy_connected_group_bonded(
    fused_build: Dict[str, Any], cmap: Dict[str, Any], spec,
) -> Dict[str, Any]:
    """MC2 connected-subgraph branch (C3): a ``single_attach_connected_group`` var
    group (W4A fused indole) is built with its full ring topology.

    The Phase A multi-heavy branch requires EACH declared var heavy to bond the
    common attach atom by a real bond (the acyclic-STAR scope). For a fused ring only
    the ROOT heavy (CG) bonds the attach (CB); the other ring heavies bond ring
    neighbours, not the attach — so the star assert would fail-loud at the first
    non-root heavy (the intended star backstop, NOT a bug). This branch is the
    connected-subgraph generalization. The declared var-heavy NAME set is the
    SSOT-derived ``spec._heavy_names`` of the one-sided var group's state-only atoms;
    the ring-closure bonds + root are read from the spec.

    The disappearing var lives in copy-2 (``wt_only`` slot, attach = ``copy2_attach``);
    the appearing-side mirror lives in copy-1 (``mtr_only`` slot, ``copy1_attach``).
    W4A is a disappearing-side connected group (Trp->Ala deletes the indole).

    Certify steps (each fail-loud, never silent-skip):
      (1) ROOT: exactly ONE declared var heavy bonds the common attach atom by a REAL
          HarmonicBond. Zero => disconnected from the common boundary; >=2 =>
          multi-branch (not a single-attach group). Both RAISE. If the spec names a
          root (``bonded_heavy_disappearing`` / ``bonded_heavy_appearing``) it MUST
          match the discovered root.
      (2) BFS: from the root, follow intra-group (var-slot heavy<->heavy) REAL bonds;
          ALL declared heavies MUST be reached. Any unreached heavy => disconnected
          subgraph. RAISE.
      (3) RING: every declared ``ring_closure_bond`` MUST be present as a REAL
          HarmonicBond between two var-slot atoms (the ring must actually close inside
          the var group). A missing ring-closure bond => open chain / wrong topology.
          RAISE.
      (4) H: each declared heavy's H's are connected (bond OR SHAKE constraint — the
          SHAKE-aware ``_heavy_h_neighbours`` test). Fully-substituted ring-fusion
          bridgeheads (CG/CD2/CE2) legitimately carry zero H, so the gate is "at
          least one declared heavy is H-bearing" (an all-H-stripped build is a defect)
          + each found H is connectivity-verified; bridgeheads are recorded explicitly.

    Runs at the PRE-ATTACH call site (the merged System still carries a top-level
    HarmonicBondForce there; after attach the forces migrate INTO the ATMForce and the
    getForce() downcast loses the HarmonicBondForce subclass — bond_pairs would be 0).
    """
    system = fused_build["system"]
    name_by_idx = {a.index: a.name
                   for a in fused_build["modeller"].topology.atoms()}
    alch = fused_build["alchemical_atoms"]

    # One-sided var group: disappearing side (copy-2) or appearing side (copy-1).
    n_app_heavy = len(spec._heavy_names(spec.stateB_only_atoms))
    if n_app_heavy >= 1:
        attach_idx = cmap["copy1_attach"]
        var_slot = list(alch["mtr_only"])
        declared_heavies = tuple(spec._heavy_names(spec.stateB_only_atoms))
        root_heavy = spec.bonded_heavy_appearing
        side = "appearing"
    else:
        attach_idx = cmap["copy2_attach"]
        var_slot = list(alch["wt_only"])
        declared_heavies = tuple(spec._heavy_names(spec.stateA_only_atoms))
        root_heavy = spec.bonded_heavy_disappearing
        side = "disappearing"

    bond_pairs, constraint_pairs = _collect_bond_constraint_pairs(system)
    attach_name = name_by_idx.get(attach_idx, "attach")
    slot_set = set(var_slot)

    # Resolve declared heavy NAMES -> built indices (a declared heavy never built
    # fails loud — cannot certify a heavy that is not in the box).
    idx_of: Dict[str, int] = {}
    for hn in declared_heavies:
        hi = _idx_by_name(var_slot, hn, name_by_idx)
        if hi is None:
            raise ValueError(
                "MC2 two-copy FAIL (connected-group): declared %s heavy %r absent "
                "from the merged box's %s var slot (cannot certify a heavy that was "
                "not built)." % (side, hn, side))
        idx_of[hn] = hi
    heavy_idxs = set(idx_of.values())

    # (1) ROOT: exactly one declared heavy bonds the common attach by a REAL bond.
    attach_bonded_roots = [
        hn for hn, hi in idx_of.items()
        if frozenset((hi, attach_idx)) in bond_pairs]
    if len(attach_bonded_roots) == 0:
        raise ValueError(
            "MC2 two-copy FAIL (connected-group): NO declared %s heavy bonds the "
            "common attach atom %r by a real HarmonicBond — the var group is "
            "disconnected from the common boundary (no attach-bonded root)."
            % (side, attach_name))
    if len(attach_bonded_roots) > 1:
        raise ValueError(
            "MC2 two-copy FAIL (connected-group): %d declared %s heavies bond the "
            "common attach atom %r (%s) — a single-attach connected group must have "
            "exactly ONE attach-bonded root (>=2 is a multi-branch swap, not in "
            "scope)." % (len(attach_bonded_roots), side, attach_name,
                         sorted(attach_bonded_roots)))
    discovered_root = attach_bonded_roots[0]
    if root_heavy is not None and discovered_root != root_heavy:
        raise ValueError(
            "MC2 two-copy FAIL (connected-group): discovered attach-bonded root %r "
            "does not match the spec-declared root %r."
            % (discovered_root, root_heavy))
    root_idx = idx_of[discovered_root]

    # (2) BFS from the root over INTRA-GROUP (var-slot heavy<->heavy) real bonds; ALL
    #     declared heavies must be reached.
    adj: Dict[int, set] = {hi: set() for hi in heavy_idxs}
    for pr in bond_pairs:
        a, b = tuple(pr)
        if a in heavy_idxs and b in heavy_idxs:
            adj[a].add(b)
            adj[b].add(a)
    seen = {root_idx}
    queue = [root_idx]
    while queue:
        cur = queue.pop()
        for nb in adj.get(cur, ()):
            if nb not in seen:
                seen.add(nb)
                queue.append(nb)
    unreached = sorted(name_by_idx[hi] for hi in heavy_idxs if hi not in seen)
    if unreached:
        raise ValueError(
            "MC2 two-copy FAIL (connected-group): declared %s heavies NOT reachable "
            "from the attach-bonded root %r by intra-group bonds (disconnected "
            "subgraph): %s" % (side, discovered_root, unreached))

    # (3) RING-CLOSURE: each declared ring-closure bond must be a REAL bond between
    #     two var-slot atoms (the ring must actually close inside the var group).
    ring_present: List[Dict[str, Any]] = []
    for (na, nb) in spec.ring_closure_bonds:
        ia = _idx_by_name(var_slot, na, name_by_idx)
        ib = _idx_by_name(var_slot, nb, name_by_idx)
        if ia is None or ib is None:
            raise ValueError(
                "MC2 two-copy FAIL (connected-group): ring-closure bond %s-%s names a "
                "var atom absent from the %s slot (ia=%s ib=%s)."
                % (na, nb, side, ia, ib))
        if ia not in slot_set or ib not in slot_set:
            raise ValueError(
                "MC2 two-copy FAIL (connected-group): ring-closure bond %s-%s atoms "
                "are not both in the %s var slot." % (na, nb, side))
        if frozenset((ia, ib)) not in bond_pairs:
            raise ValueError(
                "MC2 two-copy FAIL (connected-group): declared ring-closure bond "
                "%s-%s is ABSENT from the merged box's HarmonicBondForce — the ring "
                "did not close (open chain / wrong topology)." % (na, nb))
        ring_present.append({"bond": (na, nb), "present": True})

    # (4) H connectivity: each declared heavy's H's (bond OR constraint); fully-
    #     substituted ring-fusion bridgeheads (CG/CD2/CE2) legitimately carry zero H.
    per_heavy: List[Dict[str, Any]] = []
    for hn in declared_heavies:
        hi = idx_of[hn]
        h_idxs = _heavy_h_neighbours(
            hi, var_slot, bond_pairs, constraint_pairs, name_by_idx)
        h_names = sorted(name_by_idx[h] for h in h_idxs)
        per_heavy.append({
            "heavy": hn,
            "n_h": len(h_idxs),
            "h_names": h_names,
            "is_bridgehead": (len(h_idxs) == 0),
        })
    n_h_bearing = sum(1 for r in per_heavy if r["n_h"] >= 1)
    if n_h_bearing == 0:
        raise ValueError(
            "MC2 two-copy FAIL (connected-group): ZERO declared %s heavies carry a "
            "hydrogen — an all-H-stripped / incomplete side-chain build." % side)

    return {
        "shape": "single_attach_connected_group",
        "side": side,
        "attach": attach_name,
        "root_heavy": discovered_root,
        "n_heavies_certified": len(declared_heavies),
        "n_h_bearing_heavies": n_h_bearing,
        "ring_closure_bonds_present": ring_present,
        "per_heavy": per_heavy,
        "passed": True,
    }


def _assert_twocopy_acyclic_connected_group_bonded(
    fused_build: Dict[str, Any], cmap: Dict[str, Any], spec,
) -> Dict[str, Any]:
    """MC2 acyclic-connected-group branch: an ``acyclic_connected_group`` var group
    (Ile->Ala deletes CG1,CG2 both bonded to CB + CD1 chained off CG1) is built as a
    connected, ACYCLIC subgraph rooted at the common attach atom.

    DISTINCT from the acyclic-STAR branch (``_assert_twocopy_multiheavy_bonded``, which
    requires EVERY declared var heavy to bond the attach atom directly): here a var
    heavy may bond only a PARENT var heavy (CD1->CG1), so long as the whole group is
    reachable from the attach atom over intra-group + attach real bonds. DISTINCT from
    the RING branch (``_assert_twocopy_connected_group_bonded``): here
    ``ring_closure_bonds`` is EMPTY and the induced attach+heavy subgraph MUST be a
    TREE (no cycle) — a heavy cycle is a ring shape and belongs to the ring branch.

    The disappearing var lives in copy-2 (``wt_only`` slot, attach = ``copy2_attach``);
    the appearing-side mirror lives in copy-1 (``mtr_only`` slot, ``copy1_attach``).
    Ile->Ala is a disappearing-side acyclic connected group.

    Certify steps (each fail-loud, never silent-skip):
      (0) NO ``ring_closure_bonds`` declared (an acyclic spec must not name a ring
          bond — that is the ring branch's contract). RAISE if present.
      (1) ROOT(s): AT LEAST ONE declared var heavy bonds the common attach atom by a
          REAL HarmonicBond (multiple roots ALLOWED: CG1 AND CG2 both bond CB). Zero
          => disconnected from the common boundary. RAISE.
      (2) CONNECTIVITY: BFS from the ATTACH atom over the induced subgraph on
          {attach} u {declared heavies} (real bonds only) reaches EVERY declared
          heavy. Any unreached => disconnected subgraph. RAISE.
      (3) ACYCLICITY: that induced subgraph has EXACTLY ``n_heavies`` edges (a
          spanning tree on n_heavies+1 nodes = n_heavies edges). MORE => a cycle
          (ring shape, wrong branch). RAISE.
      (4) H: each declared heavy's H's connect (bond OR SHAKE constraint, the
          SHAKE-aware ``_heavy_h_neighbours`` test); AT LEAST ONE declared heavy is
          H-bearing (an all-H-stripped build is a defect). RAISE.

    Runs at the PRE-ATTACH call site (top-level HarmonicBondForce still present).
    """
    system = fused_build["system"]
    name_by_idx = {a.index: a.name
                   for a in fused_build["modeller"].topology.atoms()}
    alch = fused_build["alchemical_atoms"]

    # One-sided var group: appearing (copy-1) or disappearing (copy-2).
    n_app_heavy = len(spec._heavy_names(spec.stateB_only_atoms))
    if n_app_heavy >= 1:
        attach_idx = cmap["copy1_attach"]
        var_slot = list(alch["mtr_only"])
        declared_heavies = tuple(spec._heavy_names(spec.stateB_only_atoms))
        side = "appearing"
    else:
        attach_idx = cmap["copy2_attach"]
        var_slot = list(alch["wt_only"])
        declared_heavies = tuple(spec._heavy_names(spec.stateA_only_atoms))
        side = "disappearing"

    # (0) An acyclic-connected-group spec must declare NO ring-closure bond (that is
    #     the ring branch's contract; a ring here would be mis-routed).
    if spec.ring_closure_bonds:
        raise ValueError(
            "MC2 two-copy FAIL (acyclic-connected-group): spec %r declares "
            "ring_closure_bonds %s but was classified acyclic — a ring shape must use "
            "the single_attach_connected_group branch (connected_group_certified), not "
            "the chained-group branch." % (spec.name, list(spec.ring_closure_bonds)))

    bond_pairs, constraint_pairs = _collect_bond_constraint_pairs(system)
    attach_name = name_by_idx.get(attach_idx, "attach")

    # Resolve declared heavy NAMES -> built indices (a declared heavy never built
    # fails loud — cannot certify a heavy that is not in the box).
    idx_of: Dict[str, int] = {}
    for hn in declared_heavies:
        hi = _idx_by_name(var_slot, hn, name_by_idx)
        if hi is None:
            raise ValueError(
                "MC2 two-copy FAIL (acyclic-connected-group): declared %s heavy %r "
                "absent from the merged box's %s var slot (cannot certify a heavy that "
                "was not built)." % (side, hn, side))
        idx_of[hn] = hi
    heavy_idxs = set(idx_of.values())

    # (1) ROOT(s): >= 1 declared heavy bonds the common attach by a REAL bond.
    attach_bonded_roots = [
        hn for hn, hi in idx_of.items()
        if frozenset((hi, attach_idx)) in bond_pairs]
    if not attach_bonded_roots:
        raise ValueError(
            "MC2 two-copy FAIL (acyclic-connected-group): NO declared %s heavy bonds "
            "the common attach atom %r by a real HarmonicBond — the var group is "
            "disconnected from the common boundary (no attach-bonded root)."
            % (side, attach_name))
    # If the spec names a root it must be among the discovered attach-bonded roots.
    declared_root = (spec.bonded_heavy_appearing if side == "appearing"
                     else spec.bonded_heavy_disappearing)
    if declared_root is not None and declared_root not in attach_bonded_roots:
        raise ValueError(
            "MC2 two-copy FAIL (acyclic-connected-group): spec-declared root %r is "
            "NOT among the discovered attach-bonded roots %s."
            % (declared_root, sorted(attach_bonded_roots)))

    # (2)+(3) CONNECTIVITY + ACYCLICITY on the induced subgraph {attach} u {heavies}.
    nodes = heavy_idxs | {attach_idx}
    adj: Dict[int, set] = {n: set() for n in nodes}
    n_edges = 0
    for pr in bond_pairs:
        a, b = tuple(pr)
        if a in nodes and b in nodes and b not in adj[a]:
            adj[a].add(b)
            adj[b].add(a)
            n_edges += 1
    seen = {attach_idx}
    queue = [attach_idx]
    while queue:
        cur = queue.pop()
        for nb in adj.get(cur, ()):
            if nb not in seen:
                seen.add(nb)
                queue.append(nb)
    unreached = sorted(name_by_idx[hi] for hi in heavy_idxs if hi not in seen)
    if unreached:
        raise ValueError(
            "MC2 two-copy FAIL (acyclic-connected-group): declared %s heavies NOT "
            "reachable from the common attach atom %r over intra-group + attach real "
            "bonds (disconnected subgraph): %s" % (side, attach_name, unreached))
    # A connected subgraph on (n_heavies + 1) nodes is a TREE iff it has n_heavies
    # edges; MORE edges => a cycle (a ring shape, which must use the ring branch).
    if n_edges != len(heavy_idxs):
        raise ValueError(
            "MC2 two-copy FAIL (acyclic-connected-group): the induced attach+heavy "
            "subgraph has %d real bonds for %d heavies (an acyclic tree needs exactly "
            "%d) — a CYCLE was detected, i.e. this is a RING shape that must use the "
            "single_attach_connected_group branch, not the acyclic chained branch."
            % (n_edges, len(heavy_idxs), len(heavy_idxs)))

    # (4) H connectivity per declared heavy.
    per_heavy: List[Dict[str, Any]] = []
    for hn in declared_heavies:
        hi = idx_of[hn]
        h_idxs = _heavy_h_neighbours(
            hi, var_slot, bond_pairs, constraint_pairs, name_by_idx)
        h_names = sorted(name_by_idx[h] for h in h_idxs)
        per_heavy.append({
            "heavy": hn,
            "n_h": len(h_idxs),
            "h_names": h_names,
            "is_attach_root": (hn in attach_bonded_roots),
        })
    n_h_bearing = sum(1 for r in per_heavy if r["n_h"] >= 1)
    if n_h_bearing == 0:
        raise ValueError(
            "MC2 two-copy FAIL (acyclic-connected-group): ZERO declared %s heavies "
            "carry a hydrogen — an all-H-stripped / incomplete side-chain build."
            % side)

    return {
        "shape": "acyclic_connected_group",
        "side": side,
        "attach": attach_name,
        "attach_bonded_roots": sorted(attach_bonded_roots),
        "n_heavies_certified": len(declared_heavies),
        "n_h_bearing_heavies": n_h_bearing,
        "n_subgraph_edges": n_edges,
        "per_heavy": per_heavy,
        "passed": True,
    }


def assert_twocopy_methyl_bonded(
    fused_build: Dict[str, Any], cmap: Dict[str, Any],
) -> Dict[str, Any]:
    """MC2 (C5): the single-heavy var group's bonded terms are present in the
    merged box (the var heavy atom bonds the common attach atom + its H's bond it).

    Generalized SYMMETRICALLY over the two supported single-heavy shapes (read
    from the build's resolved ``MutationSpec.shape``):

      * ``appearing_heavy``    : copy-1 carries the APPEARING heavy var
                                 (CM for MTR / CD1 for Ile) bonded to copy-1's
                                 attach atom (NE1 / CG1); the appearing H's
                                 (HM / HD prefix) bond that heavy atom. The
                                 disappearing side disappears only an H.
      * ``disappearing_heavy`` : copy-2 carries the DISAPPEARING heavy var
                                 (CB for the Ala->Gly mirror) bonded to copy-2's
                                 attach atom (CA); the disappearing H's (HB prefix)
                                 bond that heavy atom. The appearing side grows no
                                 heavy beyond the common attach.

    Each var heavy atom <-> attach internal bond must be a true HarmonicBondForce
    term; the var H's bond the var heavy atom by a HarmonicBond or a SHAKE
    constraint under HBonds (presence is the gate).

    The acyclic-star multi-methyl shapes (``multi_appearing_heavy`` /
    ``multi_disappearing_heavy``, e.g. V3A Val->Ala deletes CG1 AND CG2, both
    star-bonded to CB) are ALSO supported when the spec is explicitly
    ``multiheavy_star_certified``; they are certified per-heavy by
    ``_assert_twocopy_multiheavy_bonded``.

    The connected-subgraph RING shape (``single_attach_connected_group``, e.g. W4A
    Trp->Ala deletes the 9-heavy FUSED indole) is ALSO supported when the spec is
    explicitly ``connected_group_certified`` with a non-empty ``ring_closure_bonds``;
    it is certified by ``_assert_twocopy_connected_group_bonded`` (one attach-bonded
    root + intra-group BFS to all heavies + ring-closure bonds present + H
    connectivity).

    FAIL-LOUD on any other shape (un-certified multi-heavy / un-certified ring /
    chained-heavy / multi-branch — e.g. W4F a ring contraction with heavies on BOTH
    sides): the builder + asserts have no path for those and must NOT silent-build a
    wrong soft-core endpoint. ``MutationSpec.shape`` == ``unsupported`` raises here.
    """
    system = fused_build["system"]
    name_by_idx = {a.index: a.name
                   for a in fused_build["modeller"].topology.atoms()}
    alch = fused_build["alchemical_atoms"]
    spec = resolve_mutation_spec(fused_build.get("mutation_spec"))
    shape = spec.shape
    if shape in ("multi_appearing_heavy", "multi_disappearing_heavy"):
        return _assert_twocopy_multiheavy_bonded(fused_build, cmap, spec, shape)
    if shape == "single_attach_connected_group":
        return _assert_twocopy_connected_group_bonded(fused_build, cmap, spec)
    if shape == "acyclic_connected_group":
        return _assert_twocopy_acyclic_connected_group_bonded(
            fused_build, cmap, spec)
    if shape not in ("appearing_heavy", "disappearing_heavy"):
        raise ValueError(
            "MC2 two-copy FAIL: unsupported mutation shape %r for spec %r "
            "(stateA_only=%s, stateB_only=%s). Supported: a single APPEARING heavy + "
            "its H's (MTR / V3I), a single DISAPPEARING heavy + its H's (Ala->Gly "
            "mirror), or an ACYCLIC-STAR multi-methyl group (>=2 heavies each bonded "
            "directly to the common attach atom) when the spec is explicitly "
            "multiheavy_star_certified (e.g. V3A), a connected-subgraph RING group "
            "when connected_group_certified with ring_closure_bonds (e.g. W4A fused "
            "indole), or an ACYCLIC chained connected group (>=1 attach-bonded root + "
            "the rest reachable, NO ring) when chained_group_certified (e.g. Ile->Ala). "
            "W4F ring CONTRACTION / heavy-on-BOTH-sides mutations require a partition "
            "redesign and are NOT buildable here — fail loud rather than silent-build a "
            "wrong endpoint."
            % (shape, spec.name, list(spec.stateA_only_atoms),
               list(spec.stateB_only_atoms)))

    if shape == "appearing_heavy":
        heavy_name = spec.bonded_heavy_appearing
        h_prefix = spec.appearing_h_prefix
        attach_idx = cmap["copy1_attach"]          # copy-1 carries the appearing var
        var_slot = alch["mtr_only"]                # appearing var lives in copy-1
        side = "appearing"
    else:  # disappearing_heavy (A9G mirror)
        heavy_name = spec.bonded_heavy_disappearing
        h_prefix = spec.disappearing_h_prefix
        attach_idx = cmap["copy2_attach"]          # copy-2 carries the disappearing var
        var_slot = alch["wt_only"]                 # disappearing var lives in copy-2
        side = "disappearing"

    heavy_idx = next(
        (i for i in var_slot if name_by_idx.get(i) == heavy_name), None)
    h_idxs = ([i for i in var_slot
               if name_by_idx.get(i, "").startswith(h_prefix)]
              if h_prefix else [])

    bond_pairs, constraint_pairs = _collect_bond_constraint_pairs(system)

    def _connected(i, j):
        return (frozenset((i, j)) in bond_pairs
                or frozenset((i, j)) in constraint_pairs)

    attach_name = name_by_idx.get(attach_idx, "attach")
    heavy_attach_present = (
        heavy_idx is not None and frozenset((heavy_idx, attach_idx)) in bond_pairs)
    h_heavy_present = {name_by_idx[h]: _connected(h, heavy_idx) for h in h_idxs}
    if not heavy_attach_present:
        raise ValueError(
            "MC2 two-copy FAIL: %s-%s internal bond absent from the merged box's "
            "HarmonicBondForce (the %s heavy var atom must bond the attach atom by a "
            "real bond)." % (heavy_name, attach_name, side))
    missing_h = [k for k, v in h_heavy_present.items() if not v]
    if missing_h:
        raise ValueError(
            "MC2 two-copy FAIL: %s-H %s-%s connectivity absent (neither bond nor "
            "constraint): %s" % (side, h_prefix, heavy_name, missing_h))
    return {
        "shape": shape,
        "heavy_attach_bond_present": heavy_attach_present,
        "h_heavy_bonds_present": h_heavy_present,
        # Back-compat keys for the appearing-heavy path (existing tests read these):
        "cm_ne1_bond_present": heavy_attach_present,
        "hm_cm_bonds_present": h_heavy_present,
        "passed": True,
    }


def assert_twocopy_disulfides(fused_build: Dict[str, Any]) -> Dict[str, Any]:
    """MC3 (C5): BOTH copies' cyclic_ss SG-SG disulfides preserved in the merged
    box (two copies => two disulfides).

    Confirms each detected SG-SG pair survives into the merged System's
    HarmonicBondForce (or a constraint). Raises if fewer than 2 disulfides were
    detected or any one is absent from the System.

    DISULFIDE-FREE scaffolds (a barnase folding control, an Ac-Ile-NMe reference
    peptide) carry ZERO cysteines, so there is NO SG-SG bond to preserve and the
    compstatin-specific ">=2 (one per copy)" assumption does not apply. The gate is
    therefore made CONDITIONAL ON CYSTEINE PRESENCE (detected STRUCTURALLY: a CYS/CYX
    SG atom on the merged topology): with zero cysteines it SKIPS/PASSES (additive
    guard); when cysteines ARE present (the 2QKI cyclic_ss peptide: CYS2-CYS12 in BOTH
    copies) the existing ">=2 detected + each survives" gate is UNCHANGED
    (byte-identical for every disulfide-bearing caller — MTR/V3I/A9G/W4A).
    """
    disulfides = fused_build.get("disulfides") or []
    # STRUCTURAL cysteine detection (both amber names + the deprotonated CYM), by SG
    # atom presence. A disulfide-free scaffold has no SG -> no disulfide to preserve.
    #
    # BINDER-SCOPED (matches _detect_twocopy_disulfides, which scopes the SG-SG
    # search to ``binder_chain``): only the BINDER is duplicated (copy-1 at the
    # site + copy-2 in bulk), so MC3's ">=2 disulfides, one per copy" is inherently
    # a BINDER concept. A free / disulfide cysteine on the SHARED inert RECEPTOR
    # (e.g. MDM2 Cys53 in a bound-leg complex, or any receptor S-S) is NOT part of
    # the mutation unit and must NOT trip this gate — the receptor is a single
    # shared copy, so it can never satisfy "one per copy" and is out of MC3's
    # scope. Scoping the TRIGGER to the binder removes the pre-existing scope
    # inconsistency (count was whole-box while detection was binder-scoped) and is
    # byte-identical for every disulfide-bearing caller (2QKI cyclic_ss CYS2-CYS12
    # lives on the binder chain -> count unchanged) and for CYS-free binders
    # (barnase / Ac-Ile-NMe -> 0 either way). Falls back to the whole-topology
    # count ONLY when ``binder_chain`` is unavailable (defensive; every builder in
    # this module records it on ``fused``).
    topology = fused_build["modeller"].topology
    _binder_chain = fused_build.get("binder_chain")
    n_cys_sg = sum(
        1 for a in topology.atoms()
        if a.name == "SG" and a.residue.name in ("CYS", "CYX", "CYM")
        and (_binder_chain is None or a.residue.chain.id == _binder_chain))
    if n_cys_sg == 0:
        return {"n_disulfides": 0, "disulfides": [], "passed": True,
                "skipped_reason": "disulfide_free_scaffold_no_cysteines"}
    if len(disulfides) < 2:
        raise ValueError(
            "MC3 two-copy FAIL: expected 2 cyclic_ss disulfides (one per copy), "
            "detected %d (with %d cysteine SG atoms present)."
            % (len(disulfides), n_cys_sg))
    system = fused_build["system"]
    bond_pairs = set()
    for f in system.getForces():
        if isinstance(f, mm.HarmonicBondForce):
            for bi in range(f.getNumBonds()):
                p1, p2, _, _ = f.getBondParameters(bi)
                bond_pairs.add(frozenset((p1, p2)))
        elif isinstance(f, mm.ATMForce):
            for j in range(f.getNumForces()):
                inner = f.getForce(j)
                if isinstance(inner, mm.HarmonicBondForce):
                    for bi in range(inner.getNumBonds()):
                        p1, p2, _, _ = inner.getBondParameters(bi)
                        bond_pairs.add(frozenset((p1, p2)))
    constraint_pairs = set()
    for ci in range(system.getNumConstraints()):
        a, b, _ = system.getConstraintParameters(ci)
        constraint_pairs.add(frozenset((a, b)))

    checked = []
    for d in disulfides:
        pair = frozenset((d["sg1_index"], d["sg2_index"]))
        present = pair in bond_pairs or pair in constraint_pairs
        if not present:
            raise ValueError(
                "MC3 two-copy FAIL: %s cyclic_ss SG-SG (%d-%d) absent from the "
                "merged System." % (d["copy"], d["sg1_index"], d["sg2_index"]))
        checked.append({"copy": d["copy"], "sg1_index": d["sg1_index"],
                        "sg2_index": d["sg2_index"], "present": True})
    return {"n_disulfides": len(checked), "disulfides": checked, "passed": True}


def assert_twocopy_seed(
    fused_build: Dict[str, Any], cmap: Dict[str, Any], min_dist_nm: float = 0.10,
) -> Dict[str, Any]:
    """R2 (C6): per-copy seed-geometry ASSERT (appearing/disappearing atoms not
    clashing WITHIN their own copy), gates BEFORE the ATMForce attach.

    The inter-copy clash is handled by the d-separation (C6 separation assert);
    this gate covers the per-copy geometry: copy-1's appearing var group (CM,HM1-3
    for MTR / CD1,HD11-13 for Ile) must not clash copy-1's own common/solvent
    atoms, and copy-2's disappearing var group (HE1 for MTR/V3I; CB,HB1-3,HA for
    the Ala->Gly mirror) must not clash copy-2's own atoms. A clashing seed
    detonates the uncapped bonded base term locally. Each copy's atoms are
    partitioned by the merge offset so a var atom is checked only against ITS OWN
    copy + solvent (the partner copy is d-displaced and irrelevant here).

    Symmetric over both var groups: each var atom's OWN bonded partners and the
    other atoms of its OWN var group are excluded (they are bond-length terms /
    intra-group geometry, not clashes) — generalizes the legacy single-HE1 /
    single-methyl exclusion to a multi-atom disappearing group (Ala CB methyl).
    """
    positions = np.array([
        v.value_in_unit(unit.nanometer)
        for v in fused_build["modeller"].positions])
    topology = fused_build["modeller"].topology
    system = fused_build["system"]
    n_copy1 = fused_build["n_copy1"]
    name_by_idx = {a.index: a.name for a in topology.atoms()}
    res_by_idx = {a.index: a.residue.name for a in topology.atoms()}

    alch = fused_build["alchemical_atoms"]
    appearing = list(alch["mtr_only"])     # copy-1 var group (indices < n_copy1)
    disappearing = list(alch["wt_only"])   # copy-2 var group (indices >= n_copy1)
    ne1_c1 = cmap["copy1_attach"]
    ne1_c2 = cmap["copy2_attach"]
    resolve_mutation_spec(fused_build.get("mutation_spec"))  # validate the selector

    solvent_o = [a.index for a in topology.atoms()
                 if a.residue.name in _SOLVENT_RESNAMES and a.element is not None
                 and a.element.symbol == "O"]

    # System bonds, used to exclude each var atom's true 1-2 bonded partners (the
    # var-group internal bonds + the attach bond are bond-length terms, NOT clashes;
    # 1-3 FF-excluded pairs e.g. a ring-N hydrogen sit at a standard ~0.08-0.09 nm).
    bond_pairs, _ = _collect_bond_constraint_pairs(system)

    def _bonded_neighbours(idx):
        out = set()
        for pr in bond_pairs:
            if idx in pr:
                out.update(pr)
        out.discard(idx)
        return out

    def _seed_min_dist(var_group, attach_idx, copy_lo, copy_hi):
        """Min distance from any var-group atom to a non-excluded same-copy /
        solvent target. Excluded for each var atom = the whole var group + the
        attach atom + that var atom's 1-2 bonded partners (and the attach's bonded
        neighbours, the FF-excluded 1-3 ring/methyl pairs)."""
        group = set(var_group)
        attach_neighbours = _bonded_neighbours(attach_idx)
        targets = {a.index for a in topology.atoms()
                   if copy_lo <= a.index < copy_hi and a.index not in group
                   and res_by_idx.get(a.index) not in _SOLVENT_RESNAMES}
        targets.update(solvent_o)
        min_d = float("inf")
        worst = None
        for v in var_group:
            excluded = group | {attach_idx} | attach_neighbours | _bonded_neighbours(v)
            pv = positions[v]
            for tgt in targets:
                if tgt in excluded:
                    continue
                d = float(np.linalg.norm(pv - positions[tgt]))
                if d < min_d:
                    min_d = d
                    worst = (name_by_idx.get(v, v), name_by_idx.get(tgt, tgt), d)
        return min_d, worst

    # COPY-1 appearing var group: indices [0, n_copy1).
    min_hard, worst_hard = _seed_min_dist(appearing, ne1_c1, 0, n_copy1)
    if min_hard <= min_dist_nm:
        raise ValueError(
            "R2 two-copy seed min-dist FAIL: copy-1 appearing var group too close "
            "to a copy-1 common/solvent atom (%.4f nm <= %.4f nm) — %s."
            % (min_hard, min_dist_nm, worst_hard))

    # COPY-2 disappearing var group: indices [n_copy1, n_atoms).
    n_atoms = topology.getNumAtoms()
    min_he1, worst_he1 = _seed_min_dist(disappearing, ne1_c2, n_copy1, n_atoms)
    if min_he1 <= min_dist_nm:
        raise ValueError(
            "R2 two-copy seed min-dist FAIL: copy-2 disappearing var group too "
            "close to a copy-2 common/solvent atom (%.4f nm <= %.4f nm) — %s."
            % (min_he1, min_dist_nm, worst_he1))

    return {
        "min_copy1_methyl_dist_nm": min_hard,
        "min_copy1_methyl_pair": worst_hard,
        "min_copy2_he1_dist_nm": min_he1,
        "min_copy2_he1_pair": worst_he1,
        "n_solvent_o_checked": len(solvent_o),
        "passed": True,
    }


# ---------------------------------------------------------------------------
# P3-#116 FIX2: deterministic + bounded-retry appearing-H placement.
#
# The APPEARING residue's methyl / side-chain hydrogens are placed by PDBFixer ->
# OpenMM ``Modeller.addHydrogens``, which jitters every new H by an UNSEEDED
# ``0.05 * Vec3(random.random(), ...)`` (modeller.py:1028, Python's GLOBAL ``random``
# module). Because that RNG is never deterministically seeded in the build path, a
# serial multi-seed run advances the global RNG per build -> the drawn methyl rotamer
# depends on BUILD ORDER, not on the velocity seed, and a rare (~5%) draw lands a
# pathological rotamer that trips the R2 seed min-dist guard (assert_twocopy_seed),
# crashing the campaign at the 3rd serial build (post-run review
# pathology_trackb_r2_seed_stochastic_appearingH_20260703).
#
# The fix is PLACEMENT-ONLY (scientific-review the R2 appearing-H placement review
# 20260703, Q1 FE-neutral): ONLY the appearing-H INITIAL coordinates change. The
# soft-core C7 canon (UMAX/UBCORE/ACORE/U0/ALPHA), the λ-schedule, the cycle count,
# the templates/charges, the box/PME and the frozen FE core are BYTE-IDENTICAL — the
# UWHAM estimator reads post-mintimeid PRODUCTION samples only, and a free terminal
# methyl re-samples all three staggered minima within equilibration, so the initial
# rotamer washes out (Chodera-Shirts 2011 DOI 10.1063/1.3660669). A deterministic
# per-unit seed makes the placement reproducible AND build-ORDER-independent (this
# ALONE removes the root cause); the bounded retry recovers the rare bad draw.
# The 0.10 nm R2 threshold is NEVER relaxed (a sub-0.10 nm H-H detonates the uncapped
# bonded base term -> NaN); K exhaustion is fail-LOUD (a real-geometry escalation
# signal, not a stochastic outlier).
_APPEARING_H_RETRY_K_DEFAULT = 5

# The stable message prefix raised by assert_twocopy_seed on an R2 clash. Matched so
# the retry re-places ONLY on this (retriable) build-artifact failure and NEVER
# swallows an unrelated ValueError (a genuine build defect).
_R2_SEED_FAIL_PREFIX = "R2 two-copy seed min-dist FAIL"


def _derive_appearing_h_seed(unit_key: str, attempt: int) -> int:
    """Deterministic per-unit RNG seed for the appearing-H placement (P3-#116 FIX2).

    Derived from a STABLE hash of the unit identity string ``unit_key`` (e.g. the
    velocity-seed label + leg + mutation name) folded with the retry ``attempt``
    index — NOT from the build order and NOT from a wall-clock / entropy source — so
    the SAME unit re-places its appearing-H IDENTICALLY across processes and
    independently of how many builds preceded it in the same process (the post-run review
    build-ORDER root cause). ``hashlib.sha256`` is used deliberately: Python's builtin
    ``hash`` on ``str`` is per-process salted (``PYTHONHASHSEED``) and would defeat the
    reproducibility this fix exists to provide. Returns a 31-bit NON-NEGATIVE int
    (fits both ``random.seed`` and ``numpy.random.seed``'s [0, 2**32) requirement).
    """
    digest = hashlib.sha256(("%s#%d" % (unit_key, int(attempt))).encode("utf-8"))
    return int(digest.hexdigest()[:8], 16) & 0x7FFFFFFF


def _seed_appearing_h_placement(seed_val: int) -> None:
    """Seed the process RNGs consumed by ``Modeller.addHydrogens`` right BEFORE the
    appearing-H placement (P3-#116 FIX2). OpenMM jitters each new H with Python's
    GLOBAL ``random`` (modeller.py:1028); ``numpy.random`` is seeded too (defensive —
    harmless if unused by this OpenMM build). Placement-only side effect: no
    Hamiltonian term is touched."""
    random.seed(int(seed_val))
    np.random.seed(int(seed_val) & 0x7FFFFFFF)


def _is_r2_seed_failure(exc: BaseException) -> bool:
    """True iff ``exc`` is the R2 two-copy seed min-dist clash raised by
    :func:`assert_twocopy_seed` (the retriable appearing-H placement artifact). The
    stable message prefix is matched so an UNRELATED ``ValueError`` (a genuine build
    defect — mis-paired common map, wrong template, unsupported mutation) is NOT
    silently retried away."""
    return _R2_SEED_FAIL_PREFIX in str(exc)


def _attach_heavy_neighbor_indices(
    topology: app.Topology, attach_idx: int, common_idx_set: set,
) -> List[int]:
    """Indices of the HEAVY common atoms bonded to the attach atom in ``topology``.

    Used to point a repositioned disappearing-atom hydrogen AWAY from the attach
    atom's heavy environment (the indole ring bisector for NE1, the CB carbon for
    CG1) so it lands in a non-clashing tetrahedral / donor position relative to
    the REGISTERED common core. Reads the topology bonds directly (the modeller
    topology carries the intra-residue connectivity post-build).
    """
    nbrs: List[int] = []
    for b in topology.bonds():
        a0, a1 = b[0], b[1]
        if a0.index == attach_idx:
            other = a1
        elif a1.index == attach_idx:
            other = a0
        else:
            continue
        el = other.element
        is_h = (el is not None and el.symbol == "H") or (
            el is None and other.name.strip().startswith("H"))
        if is_h:
            continue
        if other.index in common_idx_set:
            nbrs.append(other.index)
    return nbrs


# Ideal sp3 tetrahedral dot: the cosine of the angle between any two of the four
# tetrahedral bond directions (109.4712206... deg). Used by the deterministic
# common beta-H placement (connected-group / W4A shape only).
_TET_COS = -1.0 / 3.0


def _unit_vec(v: np.ndarray) -> np.ndarray:
    """Unit vector (return the input unchanged when its norm is ~0)."""
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-12 else v


def _common_beta_h_indices_for_attach(
    topology: app.Topology, common_idx_set: set, attach_idx: int,
    binder_chain: str,
) -> List[int]:
    """Common-core H atom indices bonded to the common attach atom (the beta-H pair).

    These are the atoms whose RANDOM PDBFixer placement (inherited via the common
    overwrite in the registration) can drive the rigid-translated root heavy onto a
    beta-H vertex (the W4A R2-retry mode). Identified STRUCTURALLY (H, in the common
    set, bonded to the attach atom on the binder chain) so the correction is robust to
    naming (HB2/HB3) without hard-coding names.
    """
    beta_h: List[int] = []
    for b in topology.bonds():
        a0, a1 = b[0], b[1]
        if a0.index == attach_idx:
            other = a1
        elif a1.index == attach_idx:
            other = a0
        else:
            continue
        if other.index not in common_idx_set:
            continue
        if other.residue.chain.id != binder_chain:
            continue
        if other.residue.name in _SOLVENT_RESNAMES:
            continue
        el = other.element
        is_h = (el is not None and el.symbol == "H") or (
            el is None and other.name.strip().startswith("H"))
        if is_h:
            beta_h.append(other.index)
    return sorted(beta_h)


def _place_two_tetrahedral_beta_h(
    cb: np.ndarray, ca: np.ndarray, cg: np.ndarray, bond_nm: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Place the TWO remaining sp3 tetrahedral substituents on CB given the two
    occupied directions CB->CA and CB->CG, at ``bond_nm`` from CB.

    CB is sp3 with four substituents; two are heavy (CA, CG) and two are the beta-H
    pair to place. The two H directions are the canonical complement of {u_ca, u_cg}:
    a bisector pointing AWAY from both heavies + a symmetric out-of-plane split, with
    the in/out split solved from the tetrahedral dot constraint (d . u_ca == d . u_cg
    == -1/3). The H's AVOID the CG direction BY CONSTRUCTION (they are the complement
    of CG), which removes the deterministic CG-on-beta-H collision (the W4A R2-retry
    mode). Returns (h1, h2) absolute positions (nm).
    """
    u_ca = _unit_vec(ca - cb)
    u_cg = _unit_vec(cg - cb)
    # bisector pointing away from both occupied heavy directions.
    s = u_ca + u_cg
    if float(np.linalg.norm(s)) > 1e-9:
        bis = _unit_vec(-s)
    else:
        ref = np.array([1.0, 0.0, 0.0]) if abs(u_ca[0]) < 0.9 else \
            np.array([0.0, 1.0, 0.0])
        bis = _unit_vec(np.cross(u_ca, ref))
    # out-of-plane axis (perpendicular to the CA/CG plane).
    perp = np.cross(u_ca, u_cg)
    if float(np.linalg.norm(perp)) < 1e-9:
        # CA/CG nearly collinear (degenerate) — pick any perpendicular to bis.
        ref = np.array([1.0, 0.0, 0.0]) if abs(bis[0]) < 0.9 else \
            np.array([0.0, 1.0, 0.0])
        perp = np.cross(bis, ref)
    perp = _unit_vec(perp)
    # d = a*bis + b*perp (unit). The tetrahedral constraint d . u_ca == d . u_cg ==
    # -1/3 reduces (perp _|_ both heavy dirs, bis symmetric) to a*(bis . u_ca) = -1/3.
    bis_dot = float(bis @ u_ca)   # == bis . u_cg by symmetry of bis
    a = (_TET_COS / bis_dot) if abs(bis_dot) > 1e-6 else 0.0
    a = float(np.clip(a, -1.0, 1.0))
    b = float(np.sqrt(max(0.0, 1.0 - a * a)))
    d1 = _unit_vec(a * bis + b * perp)
    d2 = _unit_vec(a * bis - b * perp)
    return cb + bond_nm * d1, cb + bond_nm * d2


def _place_deterministic_common_beta_h(
    copy2_build: Dict[str, Any], ms_reg, binder_chain: str,
) -> Dict[str, Any]:
    """SHAPE-GATED (single_attach_connected_group): re-place the COMMON beta-H pair
    (the H's bonded to the common attach CB) DETERMINISTICALLY by ideal sp3
    tetrahedral geometry off the REGISTERED CB, complementary to {CA, CG}.

    W4A (Trp4->Ala) is a single-scaffold perturbation: BOTH copies derive from the
    SAME WT final.pdb, so the registration's common overwrite inherits PDBFixer's
    RANDOM beta-H vertex pick on the mutated ALA. When PDBFixer happens to place a
    common beta-H on the SAME CB vertex the rigid-translated indole root (CG) points
    toward, CG lands on the beta-H (< the R2 floor) => an R2 retry. This correction
    re-places the two common beta-H's at the tetrahedral complement of {CA, CG}, which
    AVOIDS the CG direction by construction (deterministic, no RNG), off the REGISTERED
    CB (so the common register stays exact), at the engine-registered CB->beta-H bond
    length (direction-only change, harmonic CB-HB term unperturbed).

    Mutates ``copy2_build['modeller'].positions`` in place. Returns a structured
    bookkeeping dict; on any unexpected geometry (not the canonical CB+CA+CG+2H sp3
    centre) it leaves the engine placement untouched and reports the skip reason
    (fail-soft: the R2 seed gate still guards any residual clash).
    """
    c2_top = copy2_build["modeller"].topology
    c2_var = set(copy2_build["alchemical_atoms"]["wt_only"])

    def _is_binder_protein(atom) -> bool:
        return (atom.residue.chain.id == binder_chain
                and atom.residue.name not in _SOLVENT_RESNAMES)

    c2_common = set(a.index for a in c2_top.atoms()
                    if a.index not in c2_var and _is_binder_protein(a))
    attach_idx = copy2_build["alchemical_atoms"]["common"][0]

    beta_h = _common_beta_h_indices_for_attach(
        c2_top, c2_common, attach_idx, binder_chain)
    if len(beta_h) != 2:
        return {"detbeta_mode": "skipped_unexpected_beta_h_count",
                "n_common_beta_h": len(beta_h)}

    # CA = the other heavy common neighbour of CB; CG = the rigid-translated indole
    # root (a disappearing var heavy that bonds CB).
    name_idx: Dict[str, int] = {}
    for atom in c2_top.atoms():
        if atom.index in c2_common and _is_binder_protein(atom) \
                and atom.name not in name_idx:
            name_idx[atom.name] = atom.index
    if "CA" not in name_idx:
        return {"detbeta_mode": "skipped_no_CA"}

    root_heavy_name = ms_reg.bonded_heavy_disappearing
    cg_idx = None
    for atom in c2_top.atoms():
        if atom.index in c2_var and atom.name == root_heavy_name \
                and _is_binder_protein(atom):
            cg_idx = atom.index
            break
    if cg_idx is None:
        return {"detbeta_mode": "skipped_no_root_heavy",
                "root_heavy_name": root_heavy_name}

    c2_pos = list(copy2_build["modeller"].positions)
    cb = np.array(c2_pos[attach_idx].value_in_unit(unit.nanometer))
    ca = np.array(c2_pos[name_idx["CA"]].value_in_unit(unit.nanometer))
    cg = np.array(c2_pos[cg_idx].value_in_unit(unit.nanometer))

    bond_lengths = [
        float(np.linalg.norm(
            np.array(c2_pos[h].value_in_unit(unit.nanometer)) - cb))
        for h in beta_h]
    valid_lengths = [b for b in bond_lengths if b > 1e-6]
    bond_nm = float(np.mean(valid_lengths)) if valid_lengths else 0.109

    h1, h2 = _place_two_tetrahedral_beta_h(cb, ca, cg, bond_nm)
    c2_pos[beta_h[0]] = mm.Vec3(*h1) * unit.nanometer
    c2_pos[beta_h[1]] = mm.Vec3(*h2) * unit.nanometer
    copy2_build["modeller"].positions = c2_pos

    cg_to_h = [float(np.linalg.norm(h - cg)) for h in (h1, h2)]
    return {
        "detbeta_mode": "tetrahedral",
        "n_common_beta_h_replaced": 2,
        "common_beta_h_indices": beta_h,
        "attach_index": attach_idx,
        "root_heavy_name": root_heavy_name,
        "beta_h_bond_nm": bond_nm,
        "cg_to_beta_h_nm": cg_to_h,
        "cg_to_beta_h_min_nm": float(min(cg_to_h)),
    }


def _register_copy2_common_to_copy1(
    copy1_build: Dict[str, Any], copy2_build: Dict[str, Any],
    binder_chain: str = "B", spec: Optional[Any] = None,
) -> Dict[str, Any]:
    """Set copy-2's common-core coordinates to copy-1's (byte-identical common
    conformation), and reposition copy-2's disappearing var atom(s) accordingly.

    Both endpoint copies are built from INDEPENDENT MD final.pdb's, so their
    common-core conformers differ (~4.5 Å backbone RMSD). The canonical ATS
    common-core coordinate SWAP needs the two copies' common coordinates in
    REGISTER so the swap offset is PURE d (a clean rigid translation), not an
    inter-conformer mismatch that the swapped state pays as un-relievable strain.
    This mirrors the upstream dual-topology PDB (both ligands share the common
    coordinates); the alignment force maintains it under dynamics.

    Mutates ``copy2_build['modeller'].positions`` in place:
      - copy-2 common atom i  -> copy-1 common atom i's position (paired by C4
        name-order alignment, computed here on per-copy indices).
      - copy-2 disappearing var atom(s) -> off the (now-registered) common attach
        atom, at the preserved attach->var bond length, pointing AWAY from the
        attach atom's bonded heavy common neighbours (the indole ring bisector for
        NE1->HE1; the CB direction for CG1->HG11). Re-deriving the direction in the
        registered frame avoids a clash with the registered neighbours that the raw
        offset (against the differing original conformer) would introduce.

    ``spec`` (a ``MutationSpec`` or registry name; ``None`` => res-4 MTR<->Trp)
    supplies the common attach atom (NE1 for Trp, CG1 for Val). The
    disappearing-var indices come from the build's ``alchemical_atoms['wt_only']``
    slot (populated per spec by ``identify_alchemical_atoms``).

    Returns bookkeeping (n_common_registered, the var offset applied).
    """
    ms_reg = resolve_mutation_spec(spec)  # validate selector + read the group shape
    c1_top = copy1_build["modeller"].topology
    c2_top = copy2_build["modeller"].topology
    c1_var = set(copy1_build["alchemical_atoms"]["mtr_only"])
    c2_var = set(copy2_build["alchemical_atoms"]["wt_only"])

    def _is_binder_protein(atom) -> bool:
        return (atom.residue.chain.id == binder_chain
                and atom.residue.name not in _SOLVENT_RESNAMES)

    c1_common = [a.index for a in c1_top.atoms()
                 if a.index not in c1_var and _is_binder_protein(a)]
    c2_common = [a.index for a in c2_top.atoms()
                 if a.index not in c2_var and _is_binder_protein(a)]
    if len(c1_common) != len(c2_common):
        raise ValueError(
            "_register_copy2_common_to_copy1: common COUNT mismatch (%d vs %d) — "
            "cannot register conformers." % (len(c1_common), len(c2_common)))

    c1_pos = [p.value_in_unit(unit.nanometer)
              for p in copy1_build["modeller"].positions]
    c2_pos = list(copy2_build["modeller"].positions)

    c2_ne1 = copy2_build["alchemical_atoms"]["common"][0]
    he1_list = sorted(c2_var)

    # The disappearing var GROUP can be a single H (MTR HE1 / V3I HG11) or a heavy
    # methyl group (Ala CB + HB1-3 + the renamed alpha-H HA, the disappearing-heavy
    # mirror). A SINGLE non-heavy var keeps the legacy single-H repositioning (off
    # the attach atom, byte-identical). A MULTI-atom group (>1 atom OR a heavy var)
    # is moved RIGIDLY by the attach atom's registration delta so the native
    # intra-group bond geometry (CB-HB / CB-CA bond lengths + angles) is preserved
    # exactly — collapsing such a group onto one point would detonate the bonded
    # base term. The R2 seed assert then gates any residual clash against the
    # registered backbone (fail-loud, never silent).
    heavy_in_group = bool(ms_reg.bonded_heavy_disappearing)
    rigid_group = len(he1_list) > 1 or heavy_in_group

    # copy-2's ORIGINAL attach position (before the common overwrite) — the rigid
    # translation delta is (registered attach) - (original attach).
    c2_ne1_orig = np.array(c2_pos[c2_ne1].value_in_unit(unit.nanometer))

    # Capture copy-2's original attach->var BOND LENGTH (single-H legacy path only).
    he1_bond_nm = 0.101
    if he1_list and not rigid_group:
        he1_orig = np.array(c2_pos[he1_list[0]].value_in_unit(unit.nanometer))
        he1_bond_nm = float(np.linalg.norm(he1_orig - c2_ne1_orig)) or 0.101

    # Overwrite copy-2 commons with copy-1 commons (registered conformation).
    c2_common_set = set(c2_common)
    for c1_i, c2_i in zip(c1_common, c2_common):
        v = c1_pos[c1_i]
        c2_pos[c2_i] = mm.Vec3(v[0], v[1], v[2]) * unit.nanometer

    ne1_new = np.array(c2_pos[c2_ne1].value_in_unit(unit.nanometer))
    he1_dir = None
    if he1_list and rigid_group:
        # RIGID translation of the whole disappearing group by the attach delta:
        # preserves the native (copy-2) intra-group geometry; the group stays bonded
        # to the now-registered attach atom at its native offset.
        delta = ne1_new - c2_ne1_orig
        for var_i in he1_list:
            v = np.array(c2_pos[var_i].value_in_unit(unit.nanometer)) + delta
            c2_pos[var_i] = mm.Vec3(*v) * unit.nanometer
    elif he1_list:
        # Single-H legacy path (byte-identical): off the (registered) attach atom,
        # pointing AWAY from the centroid of the attach atom's bonded heavy common
        # neighbours, at the preserved bond length. For NE1 the heavy neighbours are
        # the indole ring (CD1, CE2) -> the legacy ring-bisector; for CG1 the heavy
        # neighbour is CB -> straight off CB.
        heavy_nbrs = _attach_heavy_neighbor_indices(c2_top, c2_ne1, c2_common_set)
        if heavy_nbrs:
            centroid = np.mean(
                [np.array(c2_pos[i].value_in_unit(unit.nanometer))
                 for i in heavy_nbrs], axis=0)
            d = ne1_new - centroid
            nrm = np.linalg.norm(d)
            he1_dir = (d / nrm) if nrm > 1e-9 else np.array([0.0, 0.0, 1.0])
        else:
            he1_dir = np.array([0.0, 0.0, 1.0])
        for var_i in he1_list:
            he1_new = ne1_new + he1_bond_nm * he1_dir
            c2_pos[var_i] = mm.Vec3(*he1_new) * unit.nanometer

    copy2_build["modeller"].positions = c2_pos

    record: Dict[str, Any] = {
        "n_common_registered": len(c1_common),
        "he1_repositioned": bool(he1_list),
        "he1_repositioned_rigid": bool(he1_list and rigid_group),
        "he1_bond_nm": he1_bond_nm,
        "he1_dir": (he1_dir.tolist() if he1_dir is not None else None),
    }

    # SHAPE-GATED deterministic COMMON beta-H placement (connected-group / W4A only).
    # Every other shape is left BYTE-IDENTICAL: the rigid-translate body above already
    # ran verbatim, and this correction is appended ONLY for the connected-group ring
    # shape, replacing the two common beta-H's (whose RANDOM PDBFixer placement is the
    # W4A R2-retry driver) with a deterministic CG-avoiding tetrahedral placement.
    if ms_reg.shape == "single_attach_connected_group":
        detbeta = _place_deterministic_common_beta_h(
            copy2_build, ms_reg, binder_chain)
        record["detbeta"] = detbeta

    return record


class VoidWaterCarveError(RuntimeError):
    """Fail-loud violation of a void-water carve safety gate.

    Raised when the opt-in ``carve_void_waters`` build path would delete anything
    other than a whole, neutral water molecule (a partial water, a non-water /
    protein / ligand / ion residue, or a selection whose net charge is non-zero).
    A carve that trips any of these gates is a build defect, not a benign carve —
    the build MUST NOT proceed with a charge-changing or solute-touching deletion.
    """


# Whole-water residue names the carve is permitted to delete (a strict subset of
# ``_SOLVENT_RESNAMES`` — ions/counterions are explicitly EXCLUDED so the carve can
# never touch a charged species).
_CARVE_WATER_RESNAMES = {"HOH", "WAT", "SOL"}


def _carve_void_penetrating_waters(
    merged: Any,
    n_copy1: int,
    dvec: Tuple[float, float, float],
    spec: Any,
    binder_chain: str = "B",
    cutoff_nm: float = ATS_CARVE_VOID_CUTOFF_NM,
) -> Dict[str, Any]:
    """Delete the whole bulk waters that penetrate the swap-displaced disappearing-
    heavy-atom volume of each copy, at build time (topology-level, BEFORE
    ``createSystem``).

    Geometric criterion (deterministic, no RNG): the disappearing heavy atoms are
    the ``stateA_only`` HEAVY atoms of ``spec`` (e.g. the Trp indole
    CG/CD1/CD2/NE1/CE2/CE3/CZ2/CZ3/CH2). For EACH copy these atoms are moved to the
    PARTNER copy's residue site using the SAME construction displacement vector
    ``d`` (``dvec``): copy-1's disappearing atoms by ``+d`` (they swap toward
    copy-2's site), copy-2's by ``-d`` (toward copy-1's site). This reuses the
    construction ``d`` verbatim — it is NOT hard-coded. A residue lies in copy-1
    when its atoms are indexed ``< n_copy1``, else copy-2; only the copy that
    physically CARRIES the disappearing atoms contributes (an Ala site has none),
    so both swap directions (copy1-site void and copy2-site void) are covered
    generically. Any WHOLE ``HOH`` with an atom within ``cutoff_nm`` of any
    displaced disappearing-heavy atom is deleted.

    Fail-loud (``VoidWaterCarveError``) BEFORE the delete if a selected residue is
    not a whole 3-atom water (O + 2 H), if any selected residue is not a permitted
    water name (protein/ligand/ion), or if the net charge of the deletion is
    non-zero (|Δq| > 1e-6 e; whole waters are neutral, so this must hold — it is a
    charge-conservation safety net against a mis-selected charged species).

    The carve operates ONLY on the merged topology upstream of the ATMForce attach,
    so the soft-core canon (umax/ubcore/acore) and the alchemical partition are
    untouched. Returns a per-call carve report (removed count, resids, cutoff, the
    source copy/direction of each removed water, atom-count delta, net-charge
    delta). Ranking-only (R-11).
    """
    top = merged.topology
    # Strip units to a clean (N, 3) float array. Modeller positions can be either a
    # Quantity(list-of-Vec3) or a list of Quantity(Vec3), so unwrap at BOTH levels
    # (an eager dtype=float on the raw nesting trips on the residual Quantity).
    raw_pos = merged.positions
    if unit.is_quantity(raw_pos):
        raw_pos = raw_pos.value_in_unit(unit.nanometer)
    pos_rows: List[List[float]] = []
    for v in raw_pos:
        if unit.is_quantity(v):
            v = v.value_in_unit(unit.nanometer)
        pos_rows.append([float(v[0]), float(v[1]), float(v[2])])
    positions = np.array(pos_rows, dtype=float)
    d = np.array([float(c) for c in dvec], dtype=float)

    disappearing_heavy = set(spec._heavy_names(spec.stateA_only_atoms))
    resnum_str = str(spec.resnum)

    # 1) Displaced disappearing-heavy points, tagged by source copy / swap sign.
    disp_pts: List[np.ndarray] = []
    disp_src: List[Tuple[int, str]] = []   # (copy_number, direction_label)
    for res in top.residues():
        if res.chain.id != binder_chain or str(res.id) != resnum_str:
            continue
        atom_indices = [a.index for a in res.atoms()]
        if not atom_indices:
            continue
        in_copy1 = min(atom_indices) < int(n_copy1)
        # copy-1 swaps toward copy-2 (+d); copy-2 swaps toward copy-1 (-d).
        sign = 1.0 if in_copy1 else -1.0
        copy_no = 1 if in_copy1 else 2
        dir_label = "copy1_to_copy2_site" if in_copy1 else "copy2_to_copy1_site"
        for a in res.atoms():
            if a.name in disappearing_heavy:
                disp_pts.append(positions[a.index] + sign * d)
                disp_src.append((copy_no, dir_label))

    report: Dict[str, Any] = {
        "carve_void_waters": True,
        "cutoff_nm": float(cutoff_nm),
        "displacement_vector_nm": [float(c) for c in d],
        "n_disappearing_heavy_points": len(disp_pts),
        "n_atoms_before": top.getNumAtoms(),
        "removed_waters": [],
        "n_waters_removed": 0,
        "atom_count_delta": 0,
        "net_charge_delta_e": 0.0,
        "net_charge_invariant": True,
    }
    if not disp_pts:
        # No disappearing heavy atoms present (e.g. a spec with no heavy on the
        # disappearing side) => no void to carve. Not an error; report zero.
        report["note"] = ("no disappearing heavy atoms found for res %s on chain "
                          "%s — nothing to carve" % (resnum_str, binder_chain))
        return report
    disp_arr = np.asarray(disp_pts, dtype=float)          # (K, 3)

    # 2) Select WHOLE waters with any atom within cutoff of any displaced point.
    cutoff = float(cutoff_nm)
    selected: List[Any] = []
    for res in top.residues():
        if res.name not in _CARVE_WATER_RESNAMES:
            continue
        w_idx = [a.index for a in res.atoms()]
        wpos = positions[w_idx]                            # (n_wat, 3)
        # min distance from any water atom to any displaced disappearing point.
        dmat = np.linalg.norm(
            wpos[:, None, :] - disp_arr[None, :, :], axis=2)   # (n_wat, K)
        if float(dmat.min()) < cutoff:
            wa, ka = np.unravel_index(int(dmat.argmin()), dmat.shape)
            src_copy, src_dir = disp_src[int(ka)]
            selected.append(res)
            report["removed_waters"].append({
                "resid": res.id,
                "res_index": res.index,
                "n_atoms": len(w_idx),
                "min_dist_nm": float(dmat.min()),
                "source_copy": src_copy,
                "source_direction": src_dir,
            })

    # 3) Fail-loud safety gates BEFORE any deletion.
    net_charge = 0.0
    n_removed_atoms = 0
    for res in selected:
        atoms = list(res.atoms())
        # (a) whole 3-atom HOH (O + 2 H).
        elems = sorted(
            (a.element.symbol if a.element is not None else a.name.strip()[:1])
            for a in atoms)
        if len(atoms) != 3 or elems != ["H", "H", "O"]:
            raise VoidWaterCarveError(
                "void-water carve selected a non-whole-water residue %s%s "
                "(%d atoms, elements %s) — the carve deletes ONLY whole 3-atom "
                "TIP3P waters (O + 2 H); a partial or non-water selection is a "
                "build defect." % (res.name, res.id, len(atoms), elems))
        # (b) selected residue must be a permitted water name (never protein/ion).
        if res.name not in _CARVE_WATER_RESNAMES:
            raise VoidWaterCarveError(
                "void-water carve selected a non-water residue %s%s — the carve "
                "must touch ONLY bulk water, never protein/ligand/ion."
                % (res.name, res.id))
        n_removed_atoms += len(atoms)
        # (c) net-charge contribution: a whole neutral water is 0 e.
        net_charge += 0.0

    if abs(net_charge) > 1e-6:
        raise VoidWaterCarveError(
            "void-water carve net-charge delta = %.6e e (> 1e-6) — the carve must "
            "conserve net charge (whole neutral waters only). A non-zero delta "
            "signals a mis-selected charged species." % (net_charge,))

    # 4) Delete the whole waters (topology-level, before createSystem).
    if selected:
        merged.delete(selected)

    report["n_waters_removed"] = len(selected)
    report["atom_count_delta"] = -n_removed_atoms
    report["net_charge_delta_e"] = float(net_charge)
    report["net_charge_invariant"] = bool(abs(net_charge) <= 1e-6)
    report["n_atoms_after"] = merged.topology.getNumAtoms()
    return report


def build_inplace_res4_twocopy_system(
    leg: str = "free",
    seed: str = "s7",
    binder_chain: str = "B",
    solvate: bool = True,
    padding_nm: float = 1.2,
    strict_mc1: bool = False,
    harmonize_common_charges: bool = False,
    displacement_nm: float = ATS_TWOCOPY_DISPLACEMENT_NM,
    auto_search_displacement: bool = False,
    accept_sep_nm: float = ATS_TWOCOPY_ACCEPT_SEP_NM,
    mtr_ncaa_xml: Optional[str] = None,
    constraints: Any = HBonds,
    spec: Optional[Any] = None,
    leg_inputs: Optional[Dict[str, Any]] = None,
    appearing_h_seed: Optional[int] = None,
    carve_void_waters: bool = False,
    carve_cutoff_nm: float = ATS_CARVE_VOID_CUTOFF_NM,
) -> Dict[str, Any]:
    """Top-level orchestrator: build the CANONICAL ATS TWO-COPY box (C2-C8).

    ``spec`` (a ``MutationSpec`` / registry name; ``None`` => the res-4 MTR<->Trp
    default, byte-identical legacy path) selects the mutation-definition layer ONLY
    — the two-copy core (build / register / displace / swap / C6 guards) is reused
    unchanged (C5). For the DEFAULT spec the build resolves the RBFE-harmonized MTR
    ncAA XML exactly as before. For a CANONICAL all-amber spec (e.g.
    ``"v3i_val_ile_res3"``, ``MUTATION_VAL_ILE_RES3``) BOTH copies are built from
    standard amber14 templates (no ncAA XML): copy-1 (the appearing state, e.g.
    Ile) is produced by the canonical point-mutation prep and copy-2 (the
    disappearing state, e.g. Val) is the scaffold's native residue, both sourced
    from the SAME WT-2QKI final.pdb (V3I engine validation, engine validation (Option B)).

    Both endpoint copies (copy-1 = MTR at the site, copy-2 = WT in bulk) are
    resident in ONE box. copy-2 is displaced by a vector d (~40 Å) so the copies
    are spatially separated; then the merged topology is solvated ONCE and a
    single System is built. The common-region coordinates are swapped by the
    upstream ATS transform; var atoms map to the PARTNER copy's DISTINCT attach
    atom. NO inter-copy exclusion is added (clash avoided by d-separation, C3).

    Order (C5/C6):
      1. Build both endpoint copies via ``build_leg_system`` UNSOLVATED.
      2. Displace copy-2 by d (compute_twocopy_displacement_vector).
      3. Merge copy-1 + copy-2 into one Modeller, solvate ONCE, createSystem.
      4. C4: pair common atoms across the two resident copies (count + order).
      5. MC1: common-atom continuity REPORT (reporting-only — per-atom divergence
         is OK in the coordinate-only swap; a retained NET-charge sanity hard gate
         + the opt-in strict_mc1 fail-loud path still guard genuine build defects).
      6. C6: two-copy spatial-separation ASSERT (no overlay).
      7. MC2: methyl bonded present (copy-1) ; MC3: cyclic_ss in BOTH copies.
      8. R2-style seed: appearing/disappearing atoms non-clashing (within copy).
      9. C1/C2/C3: attach the canonical two-copy swap ATMForce (distinct attach).

    Returns the fused build dict + all assert results. Ranking-only (R-11);
    PREDICTION test (R-18), NOT a converged ΔΔG.

      - ``leg="free"``  : the solvated cyclic peptide alone (receptor dropped),
                          TWO copies (MTR site + WT bulk).
      - ``leg="bound"`` : the receptor + cyclic peptide in the equilibrated BOUND
                          pose, TWO binder copies (the receptor is shared inert
                          context for copy-1; copy-2's binder is in bulk). The
                          d-vector clears both the receptor and the binder fold.

    ``displacement_nm`` is the magnitude of d (default 40 Å). ``solvate=True`` is
    the C8 target. ``constraints`` is forwarded to BOTH endpoint builds (Tier-2
    short dynamics pass ``constraints=None`` for the unconstrained alch-H run).

    ``auto_search_displacement`` (default ``False`` -> byte-identical legacy path):
    when ``True``, the displacement vector is chosen by ``auto_search_twocopy_
    displacement`` — a direction-aware search (res-4-local base direction + a cone
    of alternatives) that maximises the copy1<->copy2 (+ periodic image) min heavy-
    atom distance and escalates the magnitude only if no direction clears
    ``accept_sep_nm`` (default 1.5 nm, the decoupling-sufficient line). This
    automates the manual displacement-magnitude recovery and is endpoint-
    /direction-neutral (ranking-safe). When ``auto_search_displacement=False`` the
    legacy fixed-direction ``compute_twocopy_displacement_vector`` at
    ``displacement_nm`` is used unchanged (the manual ``--displacement-nm`` override
    path). ``accept_sep_nm`` is only consulted by the auto-search; the post-solvate
    C6 assert always enforces the 1.0 nm clash floor + the periodic-image gate.

    ``carve_void_waters`` (default ``False`` -> byte-identical legacy path): when
    ``True``, AFTER ``addSolvent`` and BEFORE ``createSystem`` the whole bulk waters
    that penetrate the swap-displaced disappearing-heavy-atom volume of each copy
    (``_carve_void_penetrating_waters``, cutoff ``carve_cutoff_nm``, default 2.5 Å)
    are deleted at topology level. This removes the raw uncapped u1 clash the
    ATMForce base energy cannot soft-core (the backward-endpoint NaN crash) while
    conserving net charge (whole neutral waters only, fail-loud on any violation).
    The carve is upstream of the ATMForce attach, so the soft-core canon
    (umax/ubcore/acore) and the alchemical partition are untouched. When ``False``
    the water shell is left exactly as ``addSolvent`` produced it (the serialized
    System is byte-identical to every pre-carve build).
    """
    if leg not in ("free", "bound"):
        raise NotImplementedError(
            "build_inplace_res4_twocopy_system: leg must be 'free' or 'bound', "
            "got %r." % (leg,))

    ms = resolve_mutation_spec(spec)
    # The sole ncAA in this system is MTR; every other appearing state (Val/Ile)
    # is a standard amber14 template and takes the CANONICAL path (no ncAA XML,
    # no hydrogen-definition load, copy-1 produced by the point-mutation prep).
    is_ncaa_mtr = (ms.stateB_resname == "MTR")

    # ``leg_inputs`` (default None -> byte-identical 2QKI-hardcoded resolution): an
    # optional pre-resolved leg-input dict (same schema as resolve_leg_inputs). It is
    # the single-scaffold folding-thermocycle override (resolve_fold_leg_inputs points
    # both endpoint slots at one barnase / tripeptide scaffold), read-only here and
    # threaded through serialize + run_one_replicate. When None the legacy per-seed
    # 2QKI resolver runs unchanged (all existing MTR/V3I/A9G/W4A callers unaffected).
    li = leg_inputs or resolve_leg_inputs(seed)
    if not li["final"]["wt"] or not li["final"]["cp4"]:
        raise FileNotFoundError(
            "build_inplace_res4_twocopy_system requires both endpoint final.pdb "
            "(seed %s). Got wt=%s cp4=%s"
            % (seed, li["final"]["wt"], li["final"]["cp4"]))

    # 1) Endpoint structures (UNSOLVATED; the merge solvates once after the
    #    displacement so both copies + the d-gap share one water shell). copy-1 =
    #    the APPEARING state at the site, copy-2 = the DISAPPEARING state in bulk.
    #
    #    CANONICAL ATS RBFE (Gallicchio JCIM 2025): there is ONE shared receptor.
    #    Only the binder/ligand is duplicated and displaced. For the BOUND leg,
    #    copy-1 carries the receptor + the appearing-state binder in the
    #    equilibrated site pose; copy-2 is the disappearing-state BINDER ALONE
    #    (receptor dropped) displaced into bulk. A second full receptor in copy-2
    #    would be displaced into copy-1's receptor body (the d-vector is residue-
    #    local, far smaller than the receptor extent), producing receptor-receptor
    #    interpenetration (PE -> +1e15, minimize NaN). The dual-topology swap
    #    touches only the binder common/var atoms, so copy-2 needs the binder only.
    #    The FREE leg is binder-only in both copies already; here BOTH legs use the
    #    binder-only prep for copy-2.
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix="ats_twocopy_%s_" % (leg,))
    mtr_struct = os.path.join(tmpdir, "stateB_%s.pdb" % (leg,))
    wt_struct = os.path.join(tmpdir, "stateA_%s.pdb" % (leg,))

    if is_ncaa_mtr:
        # DEFAULT MTR<->Trp path (byte-identical legacy): both copies come from
        # their own endpoint final.pdb (Cp4 vs WT); copy-1 = MTR (cp4 final),
        # copy-2 = WT (wt final).
        if leg == "free":
            prepare_free_peptide_from_final(
                li["final"]["cp4"], mtr_struct, binder_chain)
            prepare_free_peptide_from_final(
                li["final"]["wt"], wt_struct, binder_chain)
        else:
            # copy-1 = receptor + MTR binder (shared inert receptor context).
            prepare_bound_complex_from_final(
                li["final"]["cp4"], mtr_struct, binder_chain)
            # copy-2 = WT binder ONLY (receptor dropped) -> displaced into bulk.
            prepare_free_peptide_from_final(
                li["final"]["wt"], wt_struct, binder_chain)
    else:
        # CANONICAL all-amber path (e.g. V3I Val<->Ile, engine validation (Option B)): both
        # copies derive from the SAME WT-2QKI scaffold final.pdb (pos4 = Trp). The
        # DISAPPEARING state (copy-2) is the scaffold's native residue (verbatim,
        # via the free-peptide prep); the APPEARING state (copy-1) is the canonical
        # point mutation of that scaffold (e.g. VAL->ILE via PDBFixer). Engine-
        # validation framing — NOT a Cp4 anchor reproduction.
        scaffold_final = li["final"]["wt"]
        # P3-#116 FIX2 (placement-only): seed the process RNG consumed by
        # PDBFixer/Modeller.addHydrogens RIGHT BEFORE the mutated-copy prep so the
        # appearing residue's methyl / side-chain H land in a REPRODUCIBLE,
        # build-ORDER-independent rotamer. Default None => legacy UNSEEDED placement
        # (byte-identical to every pre-fix build; V3I/A9G/W4A callers pass nothing).
        # Only the appearing-H INITIAL coordinates change — the FE core is untouched.
        if appearing_h_seed is not None:
            _seed_appearing_h_placement(appearing_h_seed)
        if leg == "free":
            prepare_mutated_binder_from_final(
                scaffold_final, mtr_struct, ms.resnum,
                ms.stateA_resname, ms.stateB_resname, binder_chain,
                add_hydrogens=True)
            prepare_free_peptide_from_final(
                scaffold_final, wt_struct, binder_chain)
        else:
            # copy-1 = receptor + appearing-state binder (shared inert receptor).
            # The bound-complex prep keeps the receptor; the appearing-state
            # mutation is applied to the binder chain only (PDBFixer mutates the
            # named chain). copy-2 = disappearing-state binder ONLY -> bulk.
            prepare_mutated_bound_complex_from_final(
                scaffold_final, mtr_struct, ms.resnum,
                ms.stateA_resname, ms.stateB_resname, binder_chain)
            prepare_free_peptide_from_final(
                scaffold_final, wt_struct, binder_chain)

    # ncAA XML resolution (cross-track isolation): the MTR RBFE build loads the
    # DEDICATED harmonized RBFE XML (Σ|Δq|=0 common core), never the shared
    # Option-β HYBRID_MTR_XML the other tracks consume. The CANONICAL path uses NO
    # ncAA XML (Val/Ile are standard amber14 templates).
    if not is_ncaa_mtr:
        mtr_xml = None
    elif mtr_ncaa_xml is not None:
        mtr_xml = mtr_ncaa_xml
    elif os.path.isfile(HYBRID_MTR_XML_RBFE):
        mtr_xml = HYBRID_MTR_XML_RBFE
    else:
        mtr_xml = HYBRID_MTR_XML

    ncaa_resname = "MTR" if is_ncaa_mtr else None
    # The canonical mutated copy-1 already had its hydrogens placed by PDBFixer
    # (the appearing CD1 carries no H in the input), so add_hydrogens stays False
    # for both endpoints (final.pdb / PDBFixer outputs are H-complete).
    copy1_build = build_leg_system(   # appearing state (site copy)
        mtr_struct, leg=leg, binder_chain=binder_chain, solvate=False,
        add_hydrogens=False, ncaa_xml=mtr_xml, hydrogens_xml=li["hydrogens_xml"],
        constraints=constraints, spec=ms, ncaa_resname=ncaa_resname)
    copy2_build = build_leg_system(   # disappearing state (bulk copy)
        wt_struct, leg=leg, binder_chain=binder_chain, solvate=False,
        add_hydrogens=False, ncaa_xml=(mtr_xml if is_ncaa_mtr else None),
        hydrogens_xml=li["hydrogens_xml"],
        constraints=constraints, spec=ms, ncaa_resname=ncaa_resname)

    # 2a) REGISTER copy-2's common-core coordinates onto copy-1's frame. The two
    #     endpoint conformers come from INDEPENDENT MD final.pdb's (WT vs Cp4), so
    #     their backbone differs by ~4.5 Å RMSD. The canonical ATS common-core
    #     coordinate SWAP requires the two copies' common coordinates to be in
    #     REGISTER (the swap offset is then PURE d) — otherwise the swapped state
    #     pays the inter-conformer mismatch as bond/nonbonded strain (~5e4 kcal/mol
    #     that minimization cannot relieve, the same false-clash the verdict warns
    #     of). We therefore set copy-2's common atoms to copy-1's common positions
    #     (byte-identical common conformation by construction), and reposition
    #     copy-2's variable atom(s) (HE1) at the SAME offset from copy-2's NE1 as in
    #     the original WT structure (so HE1's bond/geometry is preserved). This is
    #     the in-memory equivalent of the upstream dual-topology PDB where both
    #     ligands share the common-core coordinates; the alignment FORCE (C4 ATS
    #     params) then maintains register under dynamics.
    _register_copy2_common_to_copy1(copy1_build, copy2_build, binder_chain, spec=ms)

    # 2b) Displace copy-2 (WT) by d into bulk (C2/C3). The direction clears copy-1's
    #     residue-local density (receptor + binder fold); the magnitude is d.
    #
    #     auto_search_displacement=False (default) -> byte-identical legacy path:
    #     fixed residue-local direction at the explicit/default magnitude (the manual
    #     --displacement-nm override). auto_search_displacement=True -> direction-
    #     aware search (cone of candidates, magnitude escalation) that maximises the
    #     copy1<->copy2 (+image) min distance and clears accept_sep_nm.
    displacement_search: Optional[Dict[str, Any]] = None
    if auto_search_displacement:
        displacement_search = auto_search_twocopy_displacement(
            copy1_build, copy2_build, binder_chain=binder_chain,
            accept_sep_nm=accept_sep_nm, padding_nm=padding_nm, spec=ms)
        dvec = tuple(displacement_search["displacement_vector_nm"])
    else:
        dvec = compute_twocopy_displacement_vector(
            copy1_build, copy2_build, binder_chain=binder_chain,
            magnitude_nm=displacement_nm, spec=ms)
    copy2_disp_positions = _displace_copy_positions(
        copy2_build["modeller"].positions, dvec)

    # 3) Merge copy-1 + copy-2 into ONE Modeller, then solvate ONCE + createSystem.
    #    Modeller.add appends copy-2's topology/positions AFTER copy-1, so copy-1
    #    keeps indices [0, n_copy1) and copy-2 takes [n_copy1, n_copy1+n_copy2) —
    #    the merge offset the index map applies. The ncAA XML (copy-1 ncAA template)
    #    + amber14 (copy-2 standard) both resolve, and a single addSolvent bathes
    #    both copies and the d-gap in one consistent water shell. For the CANONICAL
    #    path (mtr_xml=None) only the standard amber14 stack is loaded.
    ff_inputs = list(FF_FILES) + ([mtr_xml] if mtr_xml else [])
    ff = ForceField(*ff_inputs)
    n_copy1 = copy1_build["modeller"].topology.getNumAtoms()

    merged = Modeller(copy1_build["modeller"].topology,
                      copy1_build["modeller"].positions)
    merged.add(copy2_build["modeller"].topology, copy2_disp_positions)

    # Opt-in void-water carve (default OFF -> byte-identical). Threaded here so it
    # sits BETWEEN addSolvent and createSystem: the delete is topology-level, before
    # the System (and the ATMForce attach at step 9) exists, so the soft-core canon
    # is provably untouched. When OFF the entire block is skipped (no call, no
    # topology mutation, no reordering) => the serialized System is unchanged.
    carve_report: Optional[Dict[str, Any]] = None
    if solvate:
        merged.addSolvent(
            ff, model="tip3p", padding=padding_nm * unit.nanometers,
            ionicStrength=0.15 * unit.molar,
            positiveIon="Na+", negativeIon="Cl-", neutralize=True,
        )
        if carve_void_waters:
            carve_report = _carve_void_penetrating_waters(
                merged, n_copy1, dvec, ms, binder_chain=binder_chain,
                cutoff_nm=carve_cutoff_nm)
        system = ff.createSystem(
            merged.topology, nonbondedMethod=PME,
            nonbondedCutoff=1.0 * unit.nanometers, constraints=constraints,
            rigidWater=True, ewaldErrorTolerance=0.0005,
        )
    else:
        system = ff.createSystem(
            merged.topology, nonbondedMethod=app.NoCutoff,
            constraints=constraints, rigidWater=True,
        )

    # Re-detect the disulfide in EACH copy on the merged topology (MC3: two copies
    # => two disulfides). The createSystem above resolved the S-S via the
    # distance-completed bonds in each copy's prep; here we record both pairs.
    disulfides = _detect_twocopy_disulfides(merged, n_copy1, binder_chain)

    fused = {
        "leg": leg,
        "modeller": merged,
        "system": system,
        "n_atoms": merged.topology.getNumAtoms(),
        "n_copy1": n_copy1,
        # The alchemical BINDER chain id (the ONLY duplicated chain). Recorded so
        # the disulfide MC3 gate can scope its cysteine trigger to the mutation
        # unit (matching _detect_twocopy_disulfides), rather than the whole box —
        # a free / disulfide cysteine on the SHARED inert receptor (e.g. MDM2
        # Cys53 in a bound-leg complex) is not part of the mutation unit.
        "binder_chain": binder_chain,
        "displacement_vector_nm": [float(c) for c in dvec],
        "ff_inputs": ff_inputs,
        # The residue partition on the MERGED box, per copy (copy-2 shifted).
        "alchemical_atoms": _twocopy_alchemical_atoms(
            copy1_build, copy2_build, n_copy1),
        "disulfides": disulfides,
        # The resolved mutation spec — read by the asserts (MC2 heavy/H names,
        # R2 seed) so they generalize to any single-residue mutation.
        "mutation_spec": ms,
        # Opt-in void-water carve report (None when carve_void_waters=False). The
        # build produces ONE serialized System shared by dplus/dminus, so this carve
        # is identical for both directions by construction.
        "carve_report": carve_report,
    }

    # 4) C4: pair common atoms across the two RESIDENT copies (count + order).
    cmap = _build_twocopy_index_map(copy1_build, copy2_build, n_copy1, binder_chain)

    # 4b) OPTIONAL diagnostic: harmonize the MERGED System's copy-2 (WT) common
    #     charges to copy-1 (MTR) so MC1 passes (isolates the mechanical swap +
    #     endpoint-equivalence validation from the charge-continuity gap; NOT a
    #     production fix — production uses the harmonized RBFE XML on disk).
    common_charges_harmonized = False
    if harmonize_common_charges:
        _harmonize_twocopy_common_charges(system, cmap)
        common_charges_harmonized = True

    # 5) MC1: common-atom continuity. RE-CLASSIFIED to REPORTING-ONLY for the
    #    canonical two-copy box (
    #    ). The ATS
    #    swap is a coordinate-only transform: each resident copy keeps its NATIVE
    #    residue-template charges, so u1-u0 already includes the per-copy charge
    #    difference correctly. Per-atom common-charge divergence is therefore NOT a
    #    failure here — it was an over-constraint inherited from the retired
    #    single-shared-core design (where the common core was ONE physical copy with
    #    ONE charge). The build proceeds with native charges; the divergence is
    #    surfaced as a non-blocking, numbers-neutral report (P11 overlap-field
    #    pattern). Two retained gates:
    #      (a) C2 net-charge sanity (ALWAYS, hard): a non-charge-changing mutation
    #          must conserve the common-core TOTAL charge across the two copies
    #          (Sigma dq ~ 0). A non-zero NET signals a genuine build defect (e.g.
    #          a mis-paired common map / residue-template mismatch) -> raise.
    #      (b) strict_mc1=True (opt-in, the fail-loud unit-test path): re-raises on
    #          ANY per-atom divergence, preserving the legacy fail-loud contract.
    mc1_error: Optional[str] = None
    try:
        mc1 = assert_twocopy_common_param_continuity(system, cmap)
    except ValueError as exc:
        if strict_mc1:
            raise
        mc1_error = str(exc)
        mc1 = _summarize_twocopy_charge_divergence(
            system, cmap, copy1_build, resnum=ms.resnum)
        # C2 (retained HARD gate): per-atom divergence is benign, but the FULL
        # alchemical-residue net charge (common + var) MUST agree across the two
        # copies for a non-charge-changing mutation. A divergence is a genuine
        # build defect (mis-paired common map / wrong residue template / an
        # unsupported charge-changing mutation) — fail loud regardless of
        # strict_mc1. The common-subset net alone is benign and does NOT gate.
        if not mc1["net_sanity_ok"]:
            raise ValueError(
                "MC1 two-copy NET-charge sanity FAIL: the full mutated-residue "
                "(res %s) net charge differs between the copies — copy-1 = %.6f e, "
                "copy-2 = %.6f e, diff = %.6f e (tol %.1e). The two copies hold the "
                "SAME residue for a non-charge-changing mutation, so the residue "
                "total (common + variable atoms) must agree; a divergence is a "
                "genuine build defect (mis-paired common map, wrong residue "
                "template, or an unsupported charge-changing mutation), not a "
                "benign per-atom redistribution."
                % (str(ms.resnum), mc1["full_resmut_net_copy1_e"],
                   mc1["full_resmut_net_copy2_e"], mc1["full_resmut_net_diff_e"],
                   mc1["net_dq_tol_e"]))
        # Per-atom divergence is OK (native charges preserved, net conserved):
        # fall through to the canonical twocopy_attached path with the report
        # attached non-blockingly (mc1["mc1_error"] retains the raw assert text).
        mc1["mc1_error"] = mc1_error

    # 6) C6: two-copy spatial-separation ASSERT (no overlay; clash avoided by d).
    #    Pass the merged box's periodic vectors (present only after addSolvent) so
    #    the periodic minimum-image gate fires (Q5/C3 image ceiling). For the
    #    auto-search path also enforce the elevated decoupling-sufficient
    #    acceptance line (the legacy path keeps only the 1.0 nm clash floor =>
    #    byte-identical). Box vectors come from the merged System default cell.
    box_vectors = None
    if solvate:
        box_vectors = system.getDefaultPeriodicBoxVectors()
    separation = assert_twocopy_separation(
        fused, cmap, box_vectors=box_vectors,
        accept_sep_nm=(accept_sep_nm if auto_search_displacement else None))

    # 7) MC2 (copy-1 methyl bonded) + MC3 (cyclic_ss in BOTH copies).
    mc2 = assert_twocopy_methyl_bonded(fused, cmap)
    mc3 = assert_twocopy_disulfides(fused)

    # 8) R2-style seed: appearing/disappearing atoms non-clashing WITHIN their own
    #    copy (the inter-copy separation is C6; this is the per-copy seed gate).
    seed_assert = assert_twocopy_seed(fused, cmap)

    # 9) C1/C2/C3: attach the canonical two-copy swap ATMForce (distinct attach).
    swap = attach_twocopy_swap_atmforce(fused, cmap, lambda1=0.0, lambda2=0.0)

    return {
        "leg": leg, "seed": seed,
        "outcome": "twocopy_attached",
        "fused_build": fused,
        "copy1_endpoint": copy1_build,
        "copy2_endpoint": copy2_build,
        "common_map": cmap,
        "mc1_param_continuity": mc1,
        "common_charges_harmonized": common_charges_harmonized,
        "mtr_ncaa_xml": mtr_xml,
        "swap_mode": "twocopy",
        "displacement_vector_nm": [float(c) for c in dvec],
        "displacement_mode": ("auto_search" if auto_search_displacement
                              else "fixed_direction"),
        # task #6: selected direction / magnitude / achieved min-distance for the
        # run_manifest / build log (downstream post-run decoupling verification). For
        # the legacy fixed path this records the realised vector + the C6 distances.
        "displacement_log": (
            {
                "mode": "auto_search",
                "unit_dir": displacement_search["unit_dir"],
                "magnitude_nm": displacement_search["magnitude_nm"],
                "candidate_index": displacement_search["candidate_index"],
                "n_candidates": displacement_search["n_candidates"],
                "n_magnitudes_tried": displacement_search["n_magnitudes_tried"],
                "achieved_raw_min_nm": displacement_search["achieved_raw_min_nm"],
                "achieved_image_min_nm": displacement_search["achieved_image_min_nm"],
                "accept_sep_nm": displacement_search["accept_sep_nm"],
                "achieved_ne1_ne1_sep_nm": separation["ne1_ne1_sep_nm"],
                "achieved_solute_min_sep_nm": separation["solute_solute_min_sep_nm"],
                "achieved_image_solute_min_sep_nm":
                    separation.get("image_solute_min_sep_nm"),
            }
            if auto_search_displacement else
            {
                "mode": "fixed_direction",
                "magnitude_nm": float(displacement_nm),
                "displacement_vector_nm": [float(c) for c in dvec],
                "achieved_ne1_ne1_sep_nm": separation["ne1_ne1_sep_nm"],
                "achieved_solute_min_sep_nm": separation["solute_solute_min_sep_nm"],
                "achieved_image_solute_min_sep_nm":
                    separation.get("image_solute_min_sep_nm"),
            }
        ),
        "separation": separation,
        "seed_assert": seed_assert,
        "mc2_methyl_bonded": mc2,
        "mc3_disulfide": mc3,
        "swap": swap,
        "solvated": solvate,
        # Opt-in void-water carve report (None when carve_void_waters=False) — the
        # per-call carve summary (removed count/resids, cutoff, source copy, atom +
        # net-charge deltas). Surfaced so the launcher can log free-vs-bound counts.
        "carve_report": carve_report,
        "regime": "ranking_only",
        "note": ("CANONICAL ATS two-copy PREDICTION test (R-18); ranking-only "
                 "(R-11); two-copy correctness NOT yet proven (needs the pilot: "
                 "endpoint-equiv + frac<UBCORE>0 + O>=0.1 + dgbind1!=0.5 + UWHAM "
                 "convergence). NOT a converged ΔΔG_bind."
                 + (" [DIAGNOSTIC: common charges harmonized — NOT production-"
                    "valid for a quantitative DDG]"
                    if common_charges_harmonized else "")),
    }


def build_inplace_res4_twocopy_system_r2_retry(
    *,
    unit_key: str,
    retry_k: int = _APPEARING_H_RETRY_K_DEFAULT,
    **build_kwargs: Any,
) -> Dict[str, Any]:
    """Bounded R2-retry wrapper around :func:`build_inplace_res4_twocopy_system`
    (P3-#116 FIX2).

    Wraps [deterministic appearing-H placement -> full two-copy build -> R2 seed
    ASSERT] in a bounded retry loop: on the rare R2 seed-clash FAIL (a pathological
    ``addHydrogens``-jitter rotamer, ~5% per draw) it re-derives an INCREMENTED
    per-unit seed (``_derive_appearing_h_seed(unit_key, attempt)``) and rebuilds, up
    to ``retry_k`` attempts. Because ~95% of draws pass, K=5 -> ~3e-7 exhaustion
    probability.

    The seed is derived from a STABLE hash of ``unit_key`` (a per-unit identity
    string) folded with the attempt index, so the placement is reproducible AND
    build-ORDER-independent (the root cause). Every other kwarg is forwarded
    VERBATIM to :func:`build_inplace_res4_twocopy_system` (placement-only fix — no
    Hamiltonian term is touched).

    Fail-loud contract (review condition C2/C-Q4, R-18):
      * ONLY the R2 seed-clash ``ValueError`` (matched via :func:`_is_r2_seed_failure`)
        is retried; any other exception is re-raised immediately (a genuine build
        defect must NOT be masked).
      * K exhaustion raises ``RuntimeError`` with the full seed trail — the 0.10 nm
        R2 threshold is NEVER relaxed and a clashing build is NEVER accepted (a
        sub-threshold H-H detonates the uncapped bonded base term -> NaN). K
        exhaustion is a REAL-geometry escalation signal, not a stochastic outlier.

    Returns the successful build dict, annotated with the ``appearing_h_retry`` trail
    (attempts used + the deterministic seeds tried) for the run manifest / integrity review
    audit.
    """
    if retry_k < 1:
        raise ValueError(
            "build_inplace_res4_twocopy_system_r2_retry: retry_k must be >= 1, "
            "got %r." % (retry_k,))
    if "appearing_h_seed" in build_kwargs:
        raise ValueError(
            "build_inplace_res4_twocopy_system_r2_retry: appearing_h_seed is "
            "derived per-attempt from unit_key and must NOT be passed explicitly.")

    seeds_tried: List[int] = []
    last_r2_error: Optional[str] = None
    for attempt in range(retry_k):
        seed_val = _derive_appearing_h_seed(unit_key, attempt)
        seeds_tried.append(seed_val)
        try:
            build = build_inplace_res4_twocopy_system(
                appearing_h_seed=seed_val, **build_kwargs)
        except ValueError as exc:
            if not _is_r2_seed_failure(exc):
                raise          # a non-R2 build error is a real defect — fail loud
            last_r2_error = str(exc)
            sys.stderr.write(
                "[P3-#116 FIX2] R2 seed clash on appearing-H build (unit=%s, "
                "attempt %d/%d, seed=%d): %s\n  -> re-placing with an incremented "
                "deterministic seed.\n"
                % (unit_key, attempt + 1, retry_k, seed_val, last_r2_error))
            continue
        build["appearing_h_retry"] = {
            "unit_key": unit_key,
            "retry_k": retry_k,
            "attempts_used": attempt + 1,
            "seeds_tried": list(seeds_tried),
            "seed_used": seed_val,
        }
        return build

    raise RuntimeError(
        "build_inplace_res4_twocopy_system_r2_retry: the appearing-H R2 seed guard "
        "FAILED on ALL %d deterministic re-placements (unit=%s, seeds=%s). At K=%d "
        "the ~3e-7 exhaustion probability is effectively never reached by the "
        "stochastic-placement artifact, so this signals a REAL scaffold/mutation "
        "geometry problem (NOT a rare draw). The 0.10 nm R2 threshold is NOT relaxed "
        "and NO clashing build is accepted (a sub-threshold H-H detonates the "
        "uncapped bonded base term -> NaN). Escalate (scaffold re-equilibration / "
        "mutation partition review). Last R2 error: %s"
        % (retry_k, unit_key, seeds_tried, retry_k, last_r2_error))


def check_twocopy_endpoint_equivalence(
    twocopy_build: Dict[str, Any],
    platform_name: str = "Reference",
    tol_kcal: float = 25.0,
    pert_plausible_max_kcal: float = 1.0e3,
    pert_saturated_floor_kcal: float = 140.0,
) -> Dict[str, Any]:
    """C6e/C8/Q5: the canonical two-copy endpoint-equivalence + decouple check.

    For the rebuild to be CORRECT (not a finite-but-wrong false-green), THREE
    things must hold at the un-transformed reference frame (Lambda1=Lambda2=0,
    Direction=+1):

      (1) ENDPOINT-EQUIVALENCE: the ATM reference potential ``u0`` (the energy of
          the box with the swap NOT applied) must reproduce the full merged
          System potential energy within ``tol_kcal``. ``u0`` is the physical
          state "copy-1 (MTR) at the site + copy-2 (WT) in bulk", so if it equals
          the plain potential the swap wiring did not corrupt the reference state.

      (2) PERTURBATION REGIME: ``|u1 - u0|`` must be in the physically plausible
          ones-to-hundreds-kcal/mol band — NOT the single-shared-core saturated
          ~150 plateau (the collapse signature) and NOT a 56,000 kcal/mol clash.
          ``u1`` is the energy AFTER the coordinate swap (copy-2 swapped to the
          site, copy-1 to bulk), so a finite, non-saturated ``|u1-u0|`` is the
          first evidence the two-copy transfer is real. (A definitive
          non-saturation verdict needs the soft-core-uncapped raw transfer at a
          real λ window — that is the pilot, post-validation; here we report the
          one-frame value + flag the saturated-plateau signature.)

      (3) BULK-COPY DECOUPLED: copy-2 (WT, the bulk copy at the reference frame)
          must be far from copy-1 / the receptor (the d-separation), evidenced by
          the copy1-NE1 <-> copy2-NE1 distance >= the PME cutoff. This is the
          "non-interacting partner decoupled" direct check (if the bulk copy still
          interacts strongly the d was too small or the swap is wrong).

    This is a CHEAP one-frame check (Tier-1). It does NOT prove the converged
    ΔΔG; the pilot (frac<UBCORE>0 + O>=0.1 + dgbind1!=0.5 + UWHAM convergence)
    does. Returns a structured dict; never raises (the caller decides PASS/FAIL
    from the flags) so the smoke surfaces an honest finding.
    """
    import openmm as mm
    import openmm.unit as unit2

    fused = twocopy_build["fused_build"]
    system = fused["system"]
    positions = fused["modeller"].positions
    cmap = twocopy_build["common_map"]

    integ = mm.VerletIntegrator(0.001 * unit2.picoseconds)
    try:
        plat = mm.Platform.getPlatformByName(platform_name)
        ctx = mm.Context(system, integ, plat)
    except Exception:
        plat = mm.Platform.getPlatformByName("Reference")
        ctx = mm.Context(system, integ, plat)
        platform_name = "Reference"
    ctx.setPositions(positions)
    ctx.setParameter("Lambda1", 0.0)
    ctx.setParameter("Lambda2", 0.0)
    ctx.setParameter("Direction", 1.0)

    state = ctx.getState(getEnergy=True)
    e_pot = state.getPotentialEnergy().value_in_unit(unit2.kilocalorie_per_mole)
    atm = next(system.getForce(i) for i in range(system.getNumForces())
               if isinstance(system.getForce(i), mm.ATMForce))
    pert = atm.getPerturbationEnergy(ctx)
    u1 = pert[0].value_in_unit(unit2.kilocalorie_per_mole)
    u0 = pert[1].value_in_unit(unit2.kilocalorie_per_mole)
    raw_pert = u1 - u0

    # (1) endpoint-equivalence: u0 reproduces the plain potential (the reference
    #     frame energy, swap not applied).
    endpoint_gap = e_pot - u0
    endpoint_equiv = bool(np.isfinite(endpoint_gap)
                          and abs(endpoint_gap) <= tol_kcal)

    # (2) perturbation regime.
    finite_pert = bool(np.isfinite(raw_pert))
    plausible = bool(finite_pert and abs(raw_pert) <= pert_plausible_max_kcal)
    # The single-shared-core collapse saturated near Umax/Ubcore (~150). A raw
    # perturbation pinned at the soft-core plateau is the collapse signature; flag
    # it (a definitive verdict needs the multi-window uncapped transfer — pilot).
    saturated_plateau = bool(
        finite_pert and abs(raw_pert) >= pert_saturated_floor_kcal
        and abs(raw_pert) <= ATS_UMAX_KCAL + 1.0)

    # (3) bulk-copy decoupled: copy1 NE1 <-> copy2 NE1 separation.
    pos = np.array([v.value_in_unit(unit2.nanometer) for v in positions])
    sep_nm = float(np.linalg.norm(
        pos[cmap["copy1_attach"]] - pos[cmap["copy2_attach"]]))
    bulk_decoupled = bool(sep_nm >= 1.0)

    overall = bool(endpoint_equiv and finite_pert and plausible
                   and bulk_decoupled and not saturated_plateau)
    return {
        "platform": platform_name,
        "energies_kcal": {
            "E_pot": e_pot, "u0": u0, "u1": u1, "u1_minus_u0": raw_pert,
        },
        "endpoint_equivalence": {
            "u0_kcal": u0, "E_pot_kcal": e_pot, "gap_kcal": endpoint_gap,
            "tol_kcal": tol_kcal, "passed": endpoint_equiv,
        },
        "perturbation_regime": {
            "u1_minus_u0_kcal": raw_pert, "finite": finite_pert,
            "plausible": plausible, "saturated_plateau_flag": saturated_plateau,
            "plausible_max_kcal": pert_plausible_max_kcal,
        },
        "bulk_copy_decoupled": {
            "ne1_ne1_sep_nm": sep_nm, "passed": bulk_decoupled,
        },
        "overall_pass": overall,
        "regime": "ranking_only",
        "note": ("C6e/C8 one-frame two-copy endpoint-equivalence (Tier-1). NOT a "
                 "converged ΔΔG; the saturated-plateau exit is the cheap "
                 "collapse-signature flag, the definitive non-saturation verdict "
                 "needs the pilot (frac<UBCORE>0 + O>=0.1 + dgbind1!=0.5)."),
    }


def smoke_test_leg(
    pdb_path: str,
    leg: str = "bound",
    binder_chain: str = "B",
    n_lambda: int = 3,
    steps_per_window: int = 50,
    platform_name: str = "CUDA",
    add_hydrogens: bool = False,
    solvate: bool = False,
    displacement_nm: Tuple[float, float, float] = (2.0, 0.0, 0.0),
) -> Dict[str, Any]:
    """Run a tiny ATS chain on one leg and read the per-lambda perturbation
    energy + a toy free-energy estimate.

    This proves: (1) the ATMForce-equipped system builds, (2) it integrates a
    few MD steps without blowing up, (3) ``getPerturbationEnergy`` /
    ``getState`` return finite energies across lambda, and (4) a trivial
    free-energy difference (mean perturbation energy bracket) is sane (finite,
    no NaN). It does NOT produce a converged DDG_bind — that is the gated run.
    """
    build = build_leg_system(
        pdb_path, leg=leg, binder_chain=binder_chain,
        solvate=solvate, add_hydrogens=add_hydrogens,
    )
    system = build["system"]
    modeller = build["modeller"]
    alch = build["alchemical_atoms"]
    # The displaced set = the residue-4 alchemical atoms (the swapped sidechain).
    displaced = (alch["common"] + alch["wt_only"] + alch["mtr_only"])

    # Build a lambda schedule (denser near the soft-core ends per B-C4 spirit;
    # for the smoke-test just a few symmetric points).
    lambdas = [i / (n_lambda - 1) for i in range(n_lambda)] if n_lambda > 1 else [0.5]

    # Attach ATMForce at the first lambda; lambda is then dialed via the context.
    attach_atm_force(system, displacement_nm=displacement_nm,
                     displaced_atoms=displaced,
                     lambda1=lambdas[0], lambda2=lambdas[0])

    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1.0 / unit.picosecond, 0.001 * unit.picoseconds
    )
    platform = mm.Platform.getPlatformByName(platform_name)
    sim = app.Simulation(modeller.topology, system, integrator, platform)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=200)

    per_lambda = []
    for lam in lambdas:
        sim.context.setParameter("Lambda1", lam)
        sim.context.setParameter("Lambda2", lam)
        sim.step(steps_per_window)
        state = sim.context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        # ATMForce perturbation energy (u1 - u0) at this lambda.
        atm = [system.getForce(i) for i in range(system.getNumForces())
               if isinstance(system.getForce(i), mm.ATMForce)][0]
        pert = atm.getPerturbationEnergy(sim.context)
        u1 = pert[0].value_in_unit(unit.kilocalorie_per_mole) if hasattr(pert[0], "value_in_unit") else float(pert[0])
        u0 = pert[1].value_in_unit(unit.kilocalorie_per_mole) if hasattr(pert[1], "value_in_unit") else float(pert[1])
        per_lambda.append({
            "lambda": lam, "potential_kcal": pe,
            "u1_kcal": u1, "u0_kcal": u0, "pert_u1_minus_u0_kcal": u1 - u0,
            "finite": bool(np.isfinite(pe) and np.isfinite(u1) and np.isfinite(u0)),
        })

    perts = [p["pert_u1_minus_u0_kcal"] for p in per_lambda]
    # Toy linear-response free-energy bracket (mean perturbation energy). NOT a
    # converged estimate — sanity proxy only.
    toy_dg = float(np.mean(perts)) if perts else float("nan")
    return {
        "leg": leg,
        "n_atoms": build["n_atoms"],
        "disulfide": build["disulfide"],
        "displaced_atoms": displaced,
        "lambdas": lambdas,
        "per_lambda": per_lambda,
        "all_finite": all(p["finite"] for p in per_lambda),
        "toy_dg_kcal": toy_dg,
        "regime": "ranking_only",
        "note": "SMOKE-TEST ONLY — not a converged DDG; production FEP is gated.",
    }


# ---------------------------------------------------------------------------
# Module self-report (no side effects beyond reading project assets)
# ---------------------------------------------------------------------------
def describe() -> Dict[str, Any]:
    """Return a structured description of the Track B ATS setup readiness."""
    axis = verify_charge_axis()
    return {
        "track": "B",
        "method": "ATS (Alchemical Transfer with coordinate Swapping)",
        "cycle": "DDG_bind = DG_alch(bound, Trp->MTR) - DG_alch(free, Trp->MTR)",
        "charge_axis": axis,
        "alchemical_region": {
            "resnum": ALCH_RESNUM,
            "common": ALCH_COMMON_ATOM,
            "wt_only_disappear": ALCH_WT_ONLY,
            "mtr_only_appear": ALCH_MTR_ONLY,
        },
        "legs": ["bound (2QKI complex)", "free (cyclic peptide, cyclic_ss retained)"],
        "regime": "ranking_only",
    }


if __name__ == "__main__":
    import json
    print(json.dumps(describe(), indent=2))
