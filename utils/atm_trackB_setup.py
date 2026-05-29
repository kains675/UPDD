#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B — ATS (Alchemical Transfer with coordinate Swapping) system setup.

Builds the Cp4<->WT single-point relative-binding alchemical system for the
2QKI compstatin complex, per SciVal verdict
``verdict_trackAB_qm1traj_fep_prevalidation_20260529.md`` (conditions
B-C1 .. B-C7). INFRASTRUCTURE / SMOKE-TEST scope only — the production FEP
(>=11 lambda x >=5 ns x >=3 replicas) is gated on Track A's read-out and is
NOT launched here.

Thermodynamic cycle (B-C3):
    DDG_bind = DG_alch(bound, Trp->MTR) - DG_alch(free-peptide, Trp->MTR)
Both legs retain the compstatin ``cyclic_ss`` disulfide (CYS2-CYS12 SG-SG).

Charge axis (B-C1): Option-beta / regime-2 hybrid MTR
(``params/MTR_gaff2_hybrid.xml``): amber14SB-frozen backbone, NE1 frozen at
-0.3418, Khoury 2014 OMW sidechain. Formal charge 0 at both endpoints
(charge-conserving perturbation). See ``open_issues`` re: the -0.176 e partial
residual (Keeper PATCH-01 / SciVal gate before production).

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

import os
import sys
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
# Charge-axis verification (B-C1, Keeper PATCH-01 surface)
# ---------------------------------------------------------------------------
def verify_charge_axis(xml_path: str = HYBRID_MTR_XML) -> Dict[str, Any]:
    """Read the MTR residue charges from the hybrid XML and report the axis.

    Confirms the B-C1 invariants that are *checkable from the file*: NE1 frozen
    at -0.3418 and backbone N at amber14SB -0.4157. Returns the partial-charge
    sum so the caller (and Keeper) can see the -0.176 e residual explicitly
    rather than trusting the ``sigmaq0`` label.
    """
    import xml.etree.ElementTree as ET

    tree = ET.parse(xml_path)
    residues = tree.find(".//Residues")
    if residues is None:
        raise ValueError(f"No <Residues> block in {xml_path}")
    charges = {}
    for atom in residues.iter("Atom"):
        q = atom.get("charge")
        name = atom.get("name")
        if q is not None and name is not None:
            charges[name] = float(q)

    sigma_q = sum(charges.values())
    ne1 = charges.get(ALCH_COMMON_ATOM)
    bbn = charges.get("N")
    return {
        "xml_path": xml_path,
        "n_atoms": len(charges),
        "sigma_q": sigma_q,
        "ne1_charge": ne1,
        "backbone_n_charge": bbn,
        "ne1_frozen_ok": (ne1 is not None and abs(ne1 - (-0.3418)) < 1e-4),
        "backbone_n_amber14sb_ok": (bbn is not None and abs(bbn - (-0.4157)) < 1e-4),
        "charge_regime": "option_beta_regime2_hybrid",
        # The label is "sigmaq0" (formal charge 0, charge-conserving) but the
        # partial residual is NOT exactly 0 — surface it for the production gate.
        "partial_residual_nonzero": abs(sigma_q) > 1e-3,
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
    topology: app.Topology, binder_chain: str = "B", resnum: int = ALCH_RESNUM
) -> Dict[str, List[int]]:
    """Return global atom indices for the residue-4 alchemical partition.

    Classifies indices into:
      common  : NE1 (frozen, present in both states)
      wt_only : HE1 (TRP indole donor, disappears Trp->MTR)
      mtr_only: CM, HM1-3 (N-methyl, appears Trp->MTR)
    Atoms not present in the given topology are simply absent (e.g. a WT-only
    structure has no CM); the caller combines both states' lists for the
    dual-topology box.
    """
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
            if rn != resnum:
                continue
            for atom in res.atoms():
                if atom.name == ALCH_COMMON_ATOM:
                    common.append(atom.index)
                elif atom.name in ALCH_WT_ONLY:
                    wt_only.append(atom.index)
                elif atom.name in ALCH_MTR_ONLY:
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
    ncaa_xml: str = HYBRID_MTR_XML,
    hydrogens_xml: Optional[str] = HYBRID_MTR_HYDROGENS_XML,
    add_hydrogens: bool = True,
) -> Dict[str, Any]:
    """Build one leg's OpenMM ``System`` with the canonical UPDD FF stack.

    Mirrors ``run_restrained_md.py`` exactly: amber14-all + tip3pfb + MTR
    hybrid XML, tip3p water, 1.2 nm padding, 0.15 M NaCl neutralized, PME with
    1.0 nm cutoff, HBonds, rigid water. The cyclic_ss disulfide (CYS2-CYS12) is
    committed via the shared detector before ``createSystem``.

    Returns a dict with modeller, system, alchemical-atom indices, disulfide
    info and the solvated atom count. For ``solvate=False`` (smoke-test) the
    box is left unsolvated and NoCutoff is used.
    """
    ff_inputs = list(FF_FILES) + [ncaa_xml]
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
            constraints=HBonds,
            rigidWater=True,
            ewaldErrorTolerance=0.0005,
        )
    else:
        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=HBonds,
            rigidWater=True,
        )

    alch = identify_alchemical_atoms(modeller.topology, binder_chain=binder_chain)
    return {
        "leg": leg,
        "modeller": modeller,
        "system": system,
        "n_atoms": modeller.topology.getNumAtoms(),
        "alchemical_atoms": alch,
        "disulfide": disulfide,
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


# ---------------------------------------------------------------------------
# Smoke-test (B-C2): TINY ATS run proving the pipeline is end-to-end sane.
# NOT the production FEP (>=11 lambda x >=5 ns x >=3 replicas — gated run).
# ---------------------------------------------------------------------------
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
