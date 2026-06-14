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
    """
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
    with open(out_pdb, "w") as fh:
        fh.writelines(kept)
        fh.write("END\n")
    return out_pdb


def compute_decouple_direction(
    build: Dict[str, Any], binder_chain: str = "B", resnum: int = ALCH_RESNUM,
    shell_nm: float = 1.0, decouple_nm: float = 1.2,
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
    """
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
                if int(res.id) != resnum:
                    continue
            except (TypeError, ValueError):
                continue
            for atom in res.atoms():
                if atom.name == ALCH_COMMON_ATOM:
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
    constraints: Any = HBonds,
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

    # ncAA internal-bond injection: PDBFile reads MTR as HETATM without intra-
    # residue bonds. Re-add them from the MTR XML <Bond> records (identical to
    # run_restrained_md.inject_xml_bonds v38 multi-site). MUST happen before
    # the peptide-bond completion below (the peptide adder reads existing
    # bonds to skip duplicates).
    n_internal_added = inject_xml_internal_bonds(
        modeller.topology, [ncaa_xml], xml_res_name="MTR"
    )

    # ncAA peptide-bond completion: OpenMM's PDBFile reader does not infer
    # ATOM↔HETATM peptide bonds, so MTR (HETATM) junctions are missing here.
    # Add them by C(i)-N(i+1) distance threshold (0.20 nm) before
    # createSystem so amber14 templates resolve correctly. Identical logic to
    # run_restrained_md.add_missing_peptide_bonds_safe.
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

    alch = identify_alchemical_atoms(modeller.topology, binder_chain=binder_chain)
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
) -> Tuple[float, float, float]:
    """C3: the d-vector that moves copy-2 (bulk copy) clear of BOTH the receptor
    and the binder fold of copy-1.

    Reuses the local-outward logic of ``compute_decouple_direction`` (NE1 minus
    the centroid of the local heavy-atom shell around copy-1's residue-4 NE1) —
    the same leg-agnostic direction that clears the receptor AND the peptide's own
    fold for the BOUND pose, and points into bulk for the FREE peptide. Scaled to
    ``magnitude_nm``. Falls back to +X if the local outward direction is
    degenerate (the magnitude alone still separates the copies; the precise
    direction only matters so the bulk copy lands in solvent, which the box
    padding guarantees after a uniform displacement).

    NOT the single-HE1 decouple of the legacy path — here the WHOLE copy-2 is
    translated by d, so the relevant geometry is copy-1's residue-4 outward
    direction (where copy-2 must NOT collide as it is swapped to the site at
    λ=1).
    """
    unit_dir = compute_decouple_direction(copy1_build, binder_chain=binder_chain)
    if unit_dir is None:
        unit_dir = (1.0, 0.0, 0.0)
    norm = float(np.linalg.norm(np.array(unit_dir, dtype=float)))
    if norm < 1e-9:
        unit_dir = (1.0, 0.0, 0.0)
        norm = 1.0
    scale = float(magnitude_nm) / norm
    return (unit_dir[0] * scale, unit_dir[1] * scale, unit_dir[2] * scale)


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
    """MC1 (C5, HIGHEST RISK): common-atom (q, sigma, epsilon) continuity between
    the two RESIDENT copies, read from the MERGED System's NonbondedForce.

    The two-copy box EXPOSES the common-charge gap the single-shared-core box hid
    (it held one physical common copy). For the relative cycle to be exact each
    common atom's nonbonded params must be byte-identical between copy-1 (MTR) and
    copy-2 (WT). If they diverge the "common" core is not common and the swap
    injects spurious ΔE that is finite-but-wrong (a false-green). Production
    requires the harmonized RBFE XML (Σ|Δq| = 0, NMTR = amber14SB-Trp).

    Raises on any exceedance; returns the worst per-channel deviation otherwise.
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


def _summarize_twocopy_charge_divergence(
    system: mm.System, cmap: Dict[str, Any], copy1_build: Dict[str, Any],
) -> Dict[str, Any]:
    """Structured per-atom report of the copy1<->copy2 common-charge divergence.

    Two-copy analog of ``_summarize_common_charge_divergence``: when MC1 fails
    non-strictly, surface WHICH common atoms diverge and by how much (residue-4
    vs elsewhere, net displaced charge) so the pre-registered outcome (ii) is
    actionable instead of an opaque crash.
    """
    nb = next(f for f in system.getForces()
              if isinstance(f, mm.NonbondedForce))
    c1_res = {a.index: a.residue.id
              for a in copy1_build["modeller"].topology.atoms()}
    c1_name = {a.index: a.name
               for a in copy1_build["modeller"].topology.atoms()}
    per_res4: List[Dict[str, Any]] = []
    sum_dq_res4 = sum_dq_all = 0.0
    n_diverging = 0
    for c1_i, c2_i in zip(cmap["copy1_common"], cmap["copy2_common"]):
        q1 = nb.getParticleParameters(c1_i)[0].value_in_unit(unit.elementary_charge)
        q2 = nb.getParticleParameters(c2_i)[0].value_in_unit(unit.elementary_charge)
        dq = q1 - q2
        sum_dq_all += dq
        if abs(dq) > 1e-4:
            n_diverging += 1
        if str(c1_res.get(c1_i)) == str(ALCH_RESNUM):
            sum_dq_res4 += dq
            per_res4.append({"name": c1_name.get(c1_i), "q_copy1": round(q1, 4),
                             "q_copy2": round(q2, 4), "dq": round(dq, 4)})
    return {
        "passed": False,
        "n_common_checked": cmap["n_common"],
        "n_diverging": n_diverging,
        "sum_dq_res4_e": round(sum_dq_res4, 5),
        "sum_dq_all_common_e": round(sum_dq_all, 6),
        "residue4_common": per_res4,
    }


def assert_twocopy_separation(
    fused_build: Dict[str, Any], cmap: Dict[str, Any], min_sep_nm: float = 1.0,
) -> Dict[str, Any]:
    """C6: the two copies are spatially SEPARATED (no overlay) BEFORE attach.

    Asserts the minimum distance between copy-1's residue-4 NE1 and copy-2's
    residue-4 NE1 is >= ``min_sep_nm`` (the PME cutoff floor; the real
    displacement is ~40 Å). This is the structural proof that clash is avoided by
    d-separation, NOT by exclusion. A small separation means the copies overlap
    and the build collapsed back toward the (forbidden) overlay design.
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
    return {"ne1_ne1_sep_nm": sep, "min_sep_nm": min_sep_nm, "passed": True}


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


def assert_twocopy_methyl_bonded(
    fused_build: Dict[str, Any], cmap: Dict[str, Any],
) -> Dict[str, Any]:
    """MC2 (C5): copy-1 (MTR) appearing-methyl bonded terms present in the merged
    box (CM-NE1 internal bond + HM-CM connectivities).

    Two-copy analog of ``assert_methyl_bonded_present`` but reads the merged
    System: CM-NE1 must be a true HarmonicBondForce term; HM-CM may be a
    HarmonicBond or a SHAKE constraint under HBonds (presence is the gate). The
    methyl lives only in copy-1; copy-2 (WT) carries HE1 instead.
    """
    system = fused_build["system"]
    name_by_idx = {a.index: a.name
                   for a in fused_build["modeller"].topology.atoms()}
    alch = fused_build["alchemical_atoms"]
    ne1 = cmap["copy1_attach"]
    cm = next((i for i in alch["mtr_only"] if name_by_idx.get(i) == "CM"), None)
    hms = [i for i in alch["mtr_only"] if name_by_idx.get(i, "").startswith("HM")]

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

    def _connected(i, j):
        return (frozenset((i, j)) in bond_pairs
                or frozenset((i, j)) in constraint_pairs)

    cm_ne1_present = cm is not None and frozenset((cm, ne1)) in bond_pairs
    hm_cm_present = {name_by_idx[h]: _connected(h, cm) for h in hms}
    if not cm_ne1_present:
        raise ValueError(
            "MC2 two-copy FAIL: CM-NE1 internal bond absent from the merged box's "
            "HarmonicBondForce (CM is heavy — must be a real bond).")
    missing_hm = [k for k, v in hm_cm_present.items() if not v]
    if missing_hm:
        raise ValueError(
            "MC2 two-copy FAIL: methyl HM-CM connectivity absent (neither bond "
            "nor constraint): %s" % missing_hm)
    return {
        "cm_ne1_bond_present": cm_ne1_present,
        "hm_cm_bonds_present": hm_cm_present,
        "passed": True,
    }


def assert_twocopy_disulfides(fused_build: Dict[str, Any]) -> Dict[str, Any]:
    """MC3 (C5): BOTH copies' cyclic_ss SG-SG disulfides preserved in the merged
    box (two copies => two disulfides).

    Confirms each detected SG-SG pair survives into the merged System's
    HarmonicBondForce (or a constraint). Raises if fewer than 2 disulfides were
    detected or any one is absent from the System.
    """
    disulfides = fused_build.get("disulfides") or []
    if len(disulfides) < 2:
        raise ValueError(
            "MC3 two-copy FAIL: expected 2 cyclic_ss disulfides (one per copy), "
            "detected %d." % (len(disulfides),))
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
    this gate covers the per-copy geometry: copy-1's appearing methyl (CM, HM1-3)
    must not clash copy-1's own common/solvent atoms, and copy-2's disappearing
    HE1 must not clash copy-2's own atoms. A clashing seed detonates the uncapped
    bonded base term locally. Each copy's atoms are partitioned by the merge
    offset so an appearing atom is only checked against ITS OWN copy + solvent
    (the partner copy is d-displaced and irrelevant here).
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
    appearing = list(alch["mtr_only"])     # copy-1 methyl (indices < n_copy1)
    he1 = alch["wt_only"][0]               # copy-2 HE1 (index >= n_copy1)
    ne1_c1 = cmap["copy1_attach"]
    ne1_c2 = cmap["copy2_attach"]

    solvent_o = [a.index for a in topology.atoms()
                 if a.residue.name in _SOLVENT_RESNAMES and a.element is not None
                 and a.element.symbol == "O"]

    # HARD targets for the COPY-1 appearing methyl = copy-1's own non-appearing
    # atoms (< n_copy1) + solvent O. Exclude the methyl's own bonded partners
    # (CM-NE1, HM-CM) which are bond-length terms, not clashes.
    cm_idx = next((a for a in appearing if name_by_idx.get(a) == "CM"), None)
    bonded_partners = {ne1_c1}
    if cm_idx is not None:
        bonded_partners.add(cm_idx)
    copy1_targets = {a.index for a in topology.atoms()
                     if a.index < n_copy1 and a.index not in appearing
                     and res_by_idx.get(a.index) not in _SOLVENT_RESNAMES}
    copy1_targets.update(solvent_o)

    min_hard = float("inf")
    worst_hard = None
    for ap in appearing:
        pa = positions[ap]
        for tgt in copy1_targets:
            if tgt in bonded_partners:
                continue
            d = float(np.linalg.norm(pa - positions[tgt]))
            if d < min_hard:
                min_hard = d
                worst_hard = (name_by_idx.get(ap, ap), name_by_idx.get(tgt, tgt), d)
    if min_hard <= min_dist_nm:
        raise ValueError(
            "R2 two-copy seed min-dist FAIL: copy-1 appearing methyl too close to "
            "a copy-1 common/solvent atom (%.4f nm <= %.4f nm) — %s."
            % (min_hard, min_dist_nm, worst_hard))

    # COPY-2 HE1: check against copy-2's own atoms (>= n_copy1) + solvent O,
    # EXCLUDING HE1's own 1-2/1-3 neighbours (NE1 bond + the ring atoms CD1/CE2
    # bonded to NE1). These are nonbonded-EXCLUDED in the FF (1-3 pairs sit at a
    # standard ~0.08-0.09 nm from a ring-N hydrogen), so counting them as a clash
    # is a false positive — the legacy seed assert likewise excludes HE1's bonded
    # partners. We read the System bonds to find NE1's bonded ring neighbours.
    he1_excluded = {ne1_c2}
    for f in system.getForces():
        if isinstance(f, mm.HarmonicBondForce):
            for bi in range(f.getNumBonds()):
                p1, p2, _, _ = f.getBondParameters(bi)
                if ne1_c2 in (p1, p2):
                    he1_excluded.add(p2 if p1 == ne1_c2 else p1)
        elif isinstance(f, mm.ATMForce):
            for j in range(f.getNumForces()):
                inner = f.getForce(j)
                if isinstance(inner, mm.HarmonicBondForce):
                    for bi in range(inner.getNumBonds()):
                        p1, p2, _, _ = inner.getBondParameters(bi)
                        if ne1_c2 in (p1, p2):
                            he1_excluded.add(p2 if p1 == ne1_c2 else p1)
    copy2_targets = {a.index for a in topology.atoms()
                     if a.index >= n_copy1 and a.index != he1
                     and res_by_idx.get(a.index) not in _SOLVENT_RESNAMES}
    copy2_targets.update(solvent_o)
    min_he1 = float("inf")
    worst_he1 = None
    for tgt in copy2_targets:
        if tgt in he1_excluded:  # NE1 bond + 1-3 ring neighbours (FF-excluded).
            continue
        d = float(np.linalg.norm(positions[he1] - positions[tgt]))
        if d < min_he1:
            min_he1 = d
            worst_he1 = (name_by_idx.get(he1, he1), name_by_idx.get(tgt, tgt), d)
    if min_he1 <= min_dist_nm:
        raise ValueError(
            "R2 two-copy seed min-dist FAIL: copy-2 HE1 too close to a copy-2 "
            "common/solvent atom (%.4f nm <= %.4f nm) — %s."
            % (min_he1, min_dist_nm, worst_he1))

    return {
        "min_copy1_methyl_dist_nm": min_hard,
        "min_copy1_methyl_pair": worst_hard,
        "min_copy2_he1_dist_nm": min_he1,
        "min_copy2_he1_pair": worst_he1,
        "n_solvent_o_checked": len(solvent_o),
        "passed": True,
    }


def _register_copy2_common_to_copy1(
    copy1_build: Dict[str, Any], copy2_build: Dict[str, Any],
    binder_chain: str = "B",
) -> Dict[str, Any]:
    """Set copy-2's common-core coordinates to copy-1's (byte-identical common
    conformation), and reposition copy-2's variable atom(s) accordingly.

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
      - copy-2 var atom (HE1) -> copy-2 NE1's NEW position + (HE1 - NE1) original
        offset, so HE1 keeps its WT bond geometry off the (now-registered) NE1.

    Returns bookkeeping (n_common_registered, the HE1 offset applied).
    """
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

    # Capture copy-2's original NE1->HE1 BOND LENGTH (preserve the bond magnitude;
    # the DIRECTION is recomputed in the registered ring frame below so HE1 does
    # not clash the registered ring — using the raw WT offset against the MTR-frame
    # ring would mis-place HE1 since the two ring conformers differ).
    he1_bond_nm = 0.101
    if he1_list:
        ne1_orig = np.array(c2_pos[c2_ne1].value_in_unit(unit.nanometer))
        he1_orig = np.array(c2_pos[he1_list[0]].value_in_unit(unit.nanometer))
        he1_bond_nm = float(np.linalg.norm(he1_orig - ne1_orig)) or 0.101

    # Overwrite copy-2 commons with copy-1 commons (registered conformation).
    name_c2 = {a.index: a.name for a in c2_top.atoms()}
    c2_common_by_name = {name_c2[i]: i for i in c2_common}
    for c1_i, c2_i in zip(c1_common, c2_common):
        v = c1_pos[c1_i]
        c2_pos[c2_i] = mm.Vec3(v[0], v[1], v[2]) * unit.nanometer

    # Reposition HE1 in the REGISTERED ring frame: off the (now copy-1-framed) NE1,
    # pointing AWAY from the (CD1,CE2) ring bisector at the preserved bond length.
    # This is the same indole-donor geometry the legacy HE1 injector uses, so HE1
    # lands at the real Trp NE1-HE1 site relative to the registered ring (no clash
    # with the registered CD1/CE2).
    he1_dir = None
    if he1_list:
        ne1_new = np.array(c2_pos[c2_ne1].value_in_unit(unit.nanometer))
        cd1_i = c2_common_by_name.get("CD1")
        ce2_i = c2_common_by_name.get("CE2")
        if cd1_i is not None and ce2_i is not None:
            cd1 = np.array(c2_pos[cd1_i].value_in_unit(unit.nanometer))
            ce2 = np.array(c2_pos[ce2_i].value_in_unit(unit.nanometer))
            d = ne1_new - (cd1 + ce2) / 2.0
            nrm = np.linalg.norm(d)
            he1_dir = (d / nrm) if nrm > 1e-9 else np.array([0.0, 0.0, 1.0])
        else:
            he1_dir = np.array([0.0, 0.0, 1.0])
        he1_new = ne1_new + he1_bond_nm * he1_dir
        c2_pos[he1_list[0]] = mm.Vec3(*he1_new) * unit.nanometer

    copy2_build["modeller"].positions = c2_pos
    return {
        "n_common_registered": len(c1_common),
        "he1_repositioned": bool(he1_list),
        "he1_bond_nm": he1_bond_nm,
        "he1_dir": (he1_dir.tolist() if he1_dir is not None else None),
    }


def build_inplace_res4_twocopy_system(
    leg: str = "free",
    seed: str = "s7",
    binder_chain: str = "B",
    solvate: bool = True,
    padding_nm: float = 1.2,
    strict_mc1: bool = False,
    harmonize_common_charges: bool = False,
    displacement_nm: float = ATS_TWOCOPY_DISPLACEMENT_NM,
    mtr_ncaa_xml: Optional[str] = None,
    constraints: Any = HBonds,
) -> Dict[str, Any]:
    """Top-level orchestrator: build the CANONICAL ATS TWO-COPY box (C2-C8).

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
      5. MC1: common-atom (q,sigma,eps) continuity ASSERT (highest risk).
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
    """
    if leg not in ("free", "bound"):
        raise NotImplementedError(
            "build_inplace_res4_twocopy_system: leg must be 'free' or 'bound', "
            "got %r." % (leg,))

    li = resolve_leg_inputs(seed)
    if not li["final"]["wt"] or not li["final"]["cp4"]:
        raise FileNotFoundError(
            "build_inplace_res4_twocopy_system requires both endpoint final.pdb "
            "(seed %s). Got wt=%s cp4=%s"
            % (seed, li["final"]["wt"], li["final"]["cp4"]))

    # 1) Endpoint structures (UNSOLVATED; the merge solvates once after the
    #    displacement so both copies + the d-gap share one water shell). copy-1 =
    #    MTR (site), copy-2 = WT (bulk).
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix="ats_twocopy_%s_" % (leg,))
    mtr_struct = os.path.join(tmpdir, "cp4_%s.pdb" % (leg,))
    wt_struct = os.path.join(tmpdir, "wt_%s.pdb" % (leg,))
    if leg == "free":
        prepare_free_peptide_from_final(li["final"]["cp4"], mtr_struct, binder_chain)
        prepare_free_peptide_from_final(li["final"]["wt"], wt_struct, binder_chain)
    else:
        prepare_bound_complex_from_final(li["final"]["cp4"], mtr_struct, binder_chain)
        prepare_bound_complex_from_final(li["final"]["wt"], wt_struct, binder_chain)

    # RBFE MTR XML resolution (cross-track isolation): the RBFE build loads the
    # DEDICATED harmonized RBFE XML (Σ|Δq|=0 common core), never the shared
    # Option-β HYBRID_MTR_XML the other tracks consume.
    if mtr_ncaa_xml is not None:
        mtr_xml = mtr_ncaa_xml
    elif os.path.isfile(HYBRID_MTR_XML_RBFE):
        mtr_xml = HYBRID_MTR_XML_RBFE
    else:
        mtr_xml = HYBRID_MTR_XML

    copy1_build = build_leg_system(   # MTR, site copy
        mtr_struct, leg=leg, binder_chain=binder_chain, solvate=False,
        add_hydrogens=False, ncaa_xml=mtr_xml, hydrogens_xml=li["hydrogens_xml"],
        constraints=constraints)
    copy2_build = build_leg_system(   # WT, bulk copy
        wt_struct, leg=leg, binder_chain=binder_chain, solvate=False,
        add_hydrogens=False, hydrogens_xml=li["hydrogens_xml"],
        constraints=constraints)

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
    _register_copy2_common_to_copy1(copy1_build, copy2_build, binder_chain)

    # 2b) Displace copy-2 (WT) by d into bulk (C2/C3). The direction clears copy-1's
    #     residue-4 local density (receptor + binder fold); the magnitude is d.
    dvec = compute_twocopy_displacement_vector(
        copy1_build, copy2_build, binder_chain=binder_chain,
        magnitude_nm=displacement_nm)
    copy2_disp_positions = _displace_copy_positions(
        copy2_build["modeller"].positions, dvec)

    # 3) Merge copy-1 + copy-2 into ONE Modeller, then solvate ONCE + createSystem.
    #    Modeller.add appends copy-2's topology/positions AFTER copy-1, so copy-1
    #    keeps indices [0, n_copy1) and copy-2 takes [n_copy1, n_copy1+n_copy2) —
    #    the merge offset the index map applies. The ncAA XML (copy-1 MTR template)
    #    + amber14 (copy-2 WT standard) both resolve, and a single addSolvent
    #    bathes both copies and the d-gap in one consistent water shell.
    ff_inputs = list(FF_FILES) + [mtr_xml]
    ff = ForceField(*ff_inputs)
    n_copy1 = copy1_build["modeller"].topology.getNumAtoms()

    merged = Modeller(copy1_build["modeller"].topology,
                      copy1_build["modeller"].positions)
    merged.add(copy2_build["modeller"].topology, copy2_disp_positions)

    if solvate:
        merged.addSolvent(
            ff, model="tip3p", padding=padding_nm * unit.nanometers,
            ionicStrength=0.15 * unit.molar,
            positiveIon="Na+", negativeIon="Cl-", neutralize=True,
        )
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
        "displacement_vector_nm": [float(c) for c in dvec],
        "ff_inputs": ff_inputs,
        # The residue-4 partition on the MERGED box, per copy (copy-2 shifted).
        "alchemical_atoms": _twocopy_alchemical_atoms(
            copy1_build, copy2_build, n_copy1),
        "disulfides": disulfides,
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

    # 5) MC1: common-atom continuity (highest ncAA risk). Surface as a STRUCTURED
    #    outcome (not an opaque crash) by default so review sees the charge gap;
    #    strict_mc1=True re-raises (the fail-loud unit-test path).
    mc1_error: Optional[str] = None
    try:
        mc1 = assert_twocopy_common_param_continuity(system, cmap)
    except ValueError as exc:
        if strict_mc1:
            raise
        mc1_error = str(exc)
        mc1 = _summarize_twocopy_charge_divergence(system, cmap, copy1_build)
        return {
            "leg": leg, "seed": seed,
            "outcome": "mc1_charge_discontinuity",
            "mc1_param_continuity": mc1,
            "mc1_error": mc1_error,
            "common_map": {k: cmap[k] for k in ("n_common", "copy1_var", "copy2_var")},
            "displacement_vector_nm": [float(c) for c in dvec],
            "regime": "ranking_only",
            "note": ("PRE-REGISTERED outcome (ii): the two-copy common core is "
                     "electrostatically discontinuous (copy-1 MTR vs copy-2 WT "
                     "common charges diverge). The two-copy box EXPOSED the gap "
                     "the single-shared-core box hid. Resolve via the harmonized "
                     "RBFE XML (Σ|Δq|=0) before a meaningful ΔΔG. R-18: a real "
                     "charge-continuity finding, not a pass."),
        }

    # 6) C6: two-copy spatial-separation ASSERT (no overlay; clash avoided by d).
    separation = assert_twocopy_separation(fused, cmap)

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
        "separation": separation,
        "seed_assert": seed_assert,
        "mc2_methyl_bonded": mc2,
        "mc3_disulfide": mc3,
        "swap": swap,
        "solvated": solvate,
        "regime": "ranking_only",
        "note": ("CANONICAL ATS two-copy PREDICTION test (R-18); ranking-only "
                 "(R-11); two-copy correctness NOT yet proven (needs the pilot: "
                 "endpoint-equiv + frac<UBCORE>0 + O>=0.1 + dgbind1!=0.5 + UWHAM "
                 "convergence). NOT a converged ΔΔG_bind."
                 + (" [DIAGNOSTIC: common charges harmonized — NOT production-"
                    "valid for a quantitative DDG]"
                    if common_charges_harmonized else "")),
    }


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
