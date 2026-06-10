#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B Phase IV — OpenFE 1.11 RBFE setup + smoke for Cp4 vs WT compstatin.

Q5a PRIMARY: OpenFE 1.11 ``RelativeHybridTopologyProtocol`` — Cp4 Trp4 ↔
WT Trp4 (technically a Trp → 1-MeW perturbation at residue 4 of cyclic
compstatin). Cycle preserves CYS2-CYS12 disulfide + linear backbone (NOT
head-to-tail cyclic — verified SG-SG 2.04 Å, N1-OXT13 12.0 Å).

# Architecture decision (engineering)

OpenFE 1.11's ``RelativeHybridTopologyProtocol`` REQUIRES the
alchemical/transforming entity to be a ``SmallMoleculeComponent`` (see
``openfe/protocols/openmm_rfe/hybridtop_protocols.py`` L245-249). For a
13-res cyclic peptide ligand, we therefore convert the binder PDB → RDKit
Mol → SmallMoleculeComponent. The receptor (MASP-2, 631 res) stays as a
``ProteinComponent`` (static, both states).

Atom mapping: ``LomapAtomMapper`` with element + ring constraints. The
Cp4↔WT perturbation differs by exactly +4 atoms (CM + HM1-3 on indole N1
of Trp4), so the MCS yields a common-core of ~194 atoms with 4
non-mapped atoms on stateB (Cp4 side). MTR is essentially "Trp + methyl
on N1" — no atom deletion, only addition.

# Disulfide handling

The CYS2-CYS12 SG-SG bond is detected by RDKit during PDB→Mol conversion
(within 2.2 Å). It becomes part of the molecular graph and is preserved
in both stateA (WT) and stateB (Cp4) — i.e., it's in the common-core.

# Charge / partial charge strategy

OpenFE 1.11 default ``partial_charge_settings.partial_charge_method =
'am1bcc'`` would re-assign charges via openff-toolkit AmberTools sqm. For
a 198-atom peptide this is slow (~5-10 min) but tractable. The Cp4 (MTR)
side needs OFF-compatible parameters; we let openff-2.2.1 + AM1-BCC
assign — this is consistent across both endpoints and the relative ΔΔG
benefits from cancellation.

Alternative (not implemented here, deferred to Phase 4 if smoke PASS):
load custom MTR charges from ``params/MTR_gaff2_hybrid.xml`` via SDF
``<atom><charge>`` overrides. This would require RESP-equivalent
re-derivation for the WT side too, defeating the cancellation argument.

# Smoke test settings

Per condition C1:
- 1 seed × 11-λ × 250 ps equilibration on bound + free legs
- MBAR overlap off-diagonal > 0.05 ALL pairs → C1 PASS
- forward/reverse hysteresis < 1.0 kcal/mol → PASS
- FAIL → escalate to the Q5b/c fallback (GROMACS-PMX hybrid)

# Cross-references

* Prior pivot note: conda-forge openfe 0.15.0 lacks RBFE (now on 1.11)
* Cp4 source PDB: outputs/2QKI_Cp4_hybrid_calib_s7/_md_input/2QKI_Cp4.pdb
* WT source PDB:  outputs/2QKI_WT_calib_s7/_md_input/2QKI_WT.pdb
* AToM v2.1 free-leg in-progress: PID 1876426 (PRESERVED, archived post-completion)

# Output structure

    outputs/_trackb/openfe_rbfe/
    ├── {seed_tag}/
    │   ├── binder_cp4.sdf          # extracted + RDKit-parsed Cp4 binder
    │   ├── binder_wt.sdf           # extracted + RDKit-parsed WT binder
    │   ├── receptor.pdb            # MASP-2 receptor only
    │   ├── mapping.json            # LomapAtomMapper atom correspondences
    │   ├── mapping_score.txt       # LOMAP score (0-1, > 0.5 = good)
    │   ├── settings.json           # serialized protocol settings used
    │   ├── transform_complex.json  # bound-leg AlchemicalNetwork transformation
    │   ├── transform_solvent.json  # free-leg AlchemicalNetwork transformation
    │   ├── smoke_results.json      # MBAR overlap matrix + ΔΔG + hysteresis
    │   └── build_metadata.json     # seed + timestamp + version + paths

CLI:
    python scripts/trackb_openfe_rbfe_setup.py --build-only --seed-tag s7
    python scripts/trackb_openfe_rbfe_setup.py --smoke --seed-tag s7 \\
        --equilibration-ps 250 --production-ps 0 --cuda-device 0
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Environment hygiene: ensure conda env bin is in PATH for AmberTools (sqm,
# antechamber) discovery by openff-toolkit when launched outside `conda
# activate openfe11`. openff-toolkit's `shutil.which` PATH probe fails if the
# user invokes the script via the absolute python path without activation.
# ---------------------------------------------------------------------------
_ENV_BIN = os.path.dirname(sys.executable)
if _ENV_BIN not in os.environ.get("PATH", ""):
    os.environ["PATH"] = _ENV_BIN + os.pathsep + os.environ.get("PATH", "")

# ---------------------------------------------------------------------------
# Paths (project-scope absolute)
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path("/home/san/UPDD_proj")
CP4_COMPLEX_PDB = PROJECT_ROOT / "outputs" / "2QKI_Cp4_hybrid_calib_s7" / "_md_input" / "2QKI_Cp4.pdb"
WT_COMPLEX_PDB = PROJECT_ROOT / "outputs" / "2QKI_WT_calib_s7" / "_md_input" / "2QKI_WT.pdb"
OUTPUT_ROOT = PROJECT_ROOT / "outputs" / "_trackb" / "openfe_rbfe"

# MTR ncAA force field files (canonical; same paths used by AToM v2.1 launcher)
HYBRID_MTR_XML = PROJECT_ROOT / "params" / "MTR_gaff2_hybrid.xml"
MTR_HYDROGENS_XML = PROJECT_ROOT / "outputs" / "2QKI_Cp4_hybrid_calib_s7" / "params" / "MTR_hydrogens.xml"


# ---------------------------------------------------------------------------
# Step 1: Split complex PDB → receptor.pdb + binder.pdb
# ---------------------------------------------------------------------------

def split_complex_to_receptor_binder(
    complex_pdb: Path,
    out_receptor_pdb: Path,
    out_binder_pdb: Path,
    *,
    receptor_chain: str = "A",
    binder_chain: str = "B",
) -> Tuple[int, int]:
    """Split a 2-chain prepared complex PDB into receptor and binder files.

    Returns (n_receptor_atoms, n_binder_atoms).

    Reuses chain-split logic from ``phase4_trackB_v2_make_system.py`` (no
    mechanical re-implementation, but lighter — we don't need the
    multi-chain handling since prepared PDBs are guaranteed 2-chain).
    """
    if not complex_pdb.exists():
        raise FileNotFoundError(f"Complex PDB missing: {complex_pdb}")

    receptor_lines: List[str] = []
    binder_lines: List[str] = []
    other_lines: List[str] = []

    with complex_pdb.open() as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                chain_id = line[21:22]
                if chain_id == receptor_chain:
                    receptor_lines.append(line)
                elif chain_id == binder_chain:
                    binder_lines.append(line)
                else:
                    other_lines.append(line)
            elif line.startswith(("CRYST1", "HEADER", "TITLE", "REMARK")):
                # Preserve CRYST1 + header for both outputs
                receptor_lines.insert(0, line)
                binder_lines.insert(0, line)

    if not receptor_lines:
        raise ValueError(f"No chain {receptor_chain} atoms found in {complex_pdb}")
    if not binder_lines:
        raise ValueError(f"No chain {binder_chain} atoms found in {complex_pdb}")
    if other_lines:
        print(f"  [warn] {len(other_lines)} atoms with chain != A/B ignored")

    out_receptor_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_binder_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_receptor_pdb.write_text("".join(receptor_lines) + "END\n")
    out_binder_pdb.write_text("".join(binder_lines) + "END\n")

    n_rec = sum(1 for L in receptor_lines if L.startswith(("ATOM", "HETATM")))
    n_bnd = sum(1 for L in binder_lines if L.startswith(("ATOM", "HETATM")))
    return n_rec, n_bnd


# ---------------------------------------------------------------------------
# Step 1b: Protonate binder via OpenMM addHydrogens + MTR_hydrogens.xml
# ---------------------------------------------------------------------------

def protonate_binder_pdb(
    raw_pdb: Path,
    out_protonated_pdb: Path,
    *,
    is_cp4: bool,
) -> int:
    """Protonate binder PDB using OpenMM Modeller.addHydrogens.

    Cp4 source PDBs ship with MTR4 having only heavy atoms (15 atoms, no H).
    WT source PDBs have full Trp protonation but addHydrogens is idempotent
    for fully-protonated residues — re-running normalizes hydrogen names
    across both endpoints.

    For Cp4, load MTR_hydrogens.xml definitions BEFORE addHydrogens so the
    HM1-3 (N-methyl on CM) + HB2/3 + HA + HD1 are added at correct positions.

    Returns the post-protonation atom count.
    """
    from openmm.app import PDBFile, Modeller, ForceField

    pdb = PDBFile(str(raw_pdb))
    modeller = Modeller(pdb.topology, pdb.positions)

    # Load hydrogen definitions BEFORE addHydrogens (Cp4 only; WT amber14SB built-in is enough)
    if is_cp4:
        if not MTR_HYDROGENS_XML.exists():
            raise FileNotFoundError(f"MTR_hydrogens.xml not found: {MTR_HYDROGENS_XML}")
        Modeller.loadHydrogenDefinitions(str(MTR_HYDROGENS_XML))

    # Use ForceField for template-aware addHydrogens (knows MTR residue topology)
    ff_files = ["amber14-all.xml"]
    if is_cp4:
        ff_files.append(str(HYBRID_MTR_XML))
    try:
        forcefield = ForceField(*ff_files)
        modeller.addHydrogens(forcefield, pH=7.0)
    except Exception as exc:
        # Fallback: addHydrogens without force field (uses default template lookup)
        print(f"  [warn] addHydrogens(ff) failed ({exc}); retry without ff")
        modeller.addHydrogens(pH=7.0)

    out_protonated_pdb.parent.mkdir(parents=True, exist_ok=True)
    with out_protonated_pdb.open("w") as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f, keepIds=True)

    return modeller.topology.getNumAtoms()


# ---------------------------------------------------------------------------
# Step 2: Binder PDB → RDKit Mol → SDF (OpenFE-compatible)
# ---------------------------------------------------------------------------

def _assign_peptide_formal_charges(mol) -> int:
    """Correct PDB-loaded peptide charges to chemically-true assignments.

    RDKit's ``MolFromPDBFile(sanitize=True)`` auto-assigns charges by
    valence: it correctly identifies N-terminal NH3+ (ILE1 N → +1) and ARG
    NH2 → +1 (whichever N got the C=N double bond in Kekulé perception).
    BUT:
      1. It does NOT assign -1 to ASP OD2 or C-term OXT (because the
         heuristic uses bond-count not H-absence).
      2. It treats every HIS as HIP (both N protonated) because the
         imidazole 4-bond rule treats both N atoms as charged candidates.
         In our amber14SB-protonated PDBs, HIS10 is HID (only HD1, no
         HE2) → neutral, so ND1 +1 is wrong.

    This function applies ADDITIVE corrections: leaves RDKit's defaults in
    place except for the cases above. The valence model stays consistent
    (we don't reset charges on already-charged Ns and break Kekulé
    perception).

    Returns the count of corrections applied.
    """
    n_corrections = 0

    # Index residues by (resnum, resname) → list of atoms
    residue_atoms: Dict[Tuple[int, str], List[Any]] = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        key = (info.GetResidueNumber(), info.GetResidueName().strip())
        residue_atoms.setdefault(key, []).append(atom)

    for (resnum, resname), atoms in residue_atoms.items():
        atom_by_name = {
            a.GetPDBResidueInfo().GetName().strip(): a for a in atoms
        }

        # ASP/GLU sidechain COO-: assign -1 to OD2 (or OE2)
        if resname == "ASP" and "OD2" in atom_by_name:
            od2 = atom_by_name["OD2"]
            if od2.GetFormalCharge() == 0 and not any(
                n.GetSymbol() == "H" for n in od2.GetNeighbors()
            ):
                od2.SetFormalCharge(-1)
                n_corrections += 1
        if resname == "GLU" and "OE2" in atom_by_name:
            oe2 = atom_by_name["OE2"]
            if oe2.GetFormalCharge() == 0 and not any(
                n.GetSymbol() == "H" for n in oe2.GetNeighbors()
            ):
                oe2.SetFormalCharge(-1)
                n_corrections += 1

        # C-terminal OXT: assign -1
        if "OXT" in atom_by_name:
            oxt = atom_by_name["OXT"]
            if oxt.GetFormalCharge() == 0:
                oxt.SetFormalCharge(-1)
                n_corrections += 1

        # HIS HID/HIE correction: if missing one of (HD1, HE2), the
        # imidazole is neutral. RDKit defaults to HIP (+1 on ND1) →
        # we need to flip ND1 back to 0 AND explicitly set bond orders.
        if resname == "HIS":
            from rdkit import Chem as _Chem
            h_names = {
                a.GetPDBResidueInfo().GetName().strip()
                for a in atoms
                if a.GetSymbol() == "H"
            }
            has_hd1 = "HD1" in h_names
            has_he2 = "HE2" in h_names
            if not (has_hd1 and has_he2):  # NOT HIP → neutral imidazole
                # HID convention (HD1 only): CG=CD2, CE1=NE2 double bonds
                # HIE convention (HE2 only): CG-CD2 single, CD2=CE1 ? actually
                # imidazole neutral HIS:
                #   HID: bonds CG-ND1(s), ND1-CE1(s), CE1=NE2(d), NE2-CD2(s), CD2=CG(d)
                #   HIE: bonds CG-CD2(s), CD2-NE2(s), NE2=CE1(d), CE1-ND1(s), ND1=CG(d)
                if "ND1" in atom_by_name and atom_by_name["ND1"].GetFormalCharge() == +1:
                    atom_by_name["ND1"].SetFormalCharge(0)
                    n_corrections += 1
                if "NE2" in atom_by_name and atom_by_name["NE2"].GetFormalCharge() == +1:
                    atom_by_name["NE2"].SetFormalCharge(0)
                    n_corrections += 1
                # Explicitly Kekulize the imidazole ring per HID/HIE convention
                cg = atom_by_name.get("CG")
                nd1 = atom_by_name.get("ND1")
                cd2 = atom_by_name.get("CD2")
                ce1 = atom_by_name.get("CE1")
                ne2 = atom_by_name.get("NE2")
                if all(x is not None for x in [cg, nd1, cd2, ce1, ne2]):
                    # Find bond objects
                    def _bond(a, b):
                        return mol.GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())

                    bonds_to_set = []
                    if has_hd1 and not has_he2:  # HID
                        bonds_to_set = [
                            (cg, nd1, _Chem.BondType.SINGLE),
                            (nd1, ce1, _Chem.BondType.SINGLE),
                            (ce1, ne2, _Chem.BondType.DOUBLE),
                            (ne2, cd2, _Chem.BondType.SINGLE),
                            (cd2, cg, _Chem.BondType.DOUBLE),
                        ]
                    elif has_he2 and not has_hd1:  # HIE
                        bonds_to_set = [
                            (cg, nd1, _Chem.BondType.DOUBLE),
                            (nd1, ce1, _Chem.BondType.SINGLE),
                            (ce1, ne2, _Chem.BondType.DOUBLE),
                            (ne2, cd2, _Chem.BondType.SINGLE),
                            (cd2, cg, _Chem.BondType.SINGLE),
                        ]
                    for a, b, btype in bonds_to_set:
                        bd = _bond(a, b)
                        if bd is not None:
                            bd.SetBondType(btype)
                            bd.SetIsAromatic(False)
                    # Clear aromatic flags on ring atoms (we set explicit Kekulé)
                    for a in [cg, nd1, cd2, ce1, ne2]:
                        a.SetIsAromatic(False)

    return n_corrections


def binder_pdb_to_sdf(
    binder_pdb: Path,
    out_sdf: Path,
    *,
    name: str,
    is_canonical_only: bool = True,
) -> Dict[str, Any]:
    """Convert binder PDB → SDF using openff polymer-aware loader (canonical AAs).

    For WT (all canonical residues) we use ``openff.toolkit.Topology.from_pdb``
    which leverages the openff residue substructure library to assign correct
    bond orders, formal charges, and tautomer state in one shot. Total charge
    comes out at 0.0 for amber14SB-protonated compstatin (verified).

    For Cp4 (contains MTR ncAA) the polymer loader FAILS — MTR isn't in the
    library and constructing ``_custom_substructures`` for a mid-chain ncAA
    requires SMARTS-label-to-PDB-atom-name correspondence that openff
    currently doesn't validate gracefully. We fall back to RDKit
    MolFromPDBFile + additive charge correction here, which produces a valid
    SDF for LomapAtomMapper consumption but FAILS at the openff toolkit's
    `_assign_partial_charges` step at MD-execution time.

    The MTR substructure registration is THE current blocker for production
    OpenFE RBFE on cyclic compstatin.

    Returns metadata dict.
    """
    from rdkit import Chem

    n_charge_assigned = 0
    used_openff_polymer = False

    if is_canonical_only:
        # Canonical-only path: openff polymer loader → SDF via openff (reliable)
        try:
            from openff.toolkit import Topology
            top = Topology.from_pdb(str(binder_pdb))
            if top.n_molecules != 1:
                raise RuntimeError(f"Expected 1 molecule, got {top.n_molecules}")
            offmol = next(top.molecules)
            offmol.to_file(str(out_sdf), file_format="sdf")
            mol = Chem.MolFromMolFile(str(out_sdf), removeHs=False, sanitize=True)
            if mol is None:
                raise RuntimeError("Re-parse of openff SDF via RDKit failed")
            used_openff_polymer = True
            n_charge_assigned = int(abs(float(offmol.total_charge.m_as("elementary_charge"))) > 0.001)
        except Exception as exc:
            print(f"  [info] openff polymer loader failed ({exc}); falling back to RDKit path")
            used_openff_polymer = False

    if not used_openff_polymer:
        # RDKit fallback for ncAA-containing binder
        mol = Chem.MolFromPDBFile(str(binder_pdb), sanitize=True, removeHs=False, proximityBonding=True)
        if mol is None:
            mol = Chem.MolFromPDBFile(str(binder_pdb), sanitize=True, removeHs=False, proximityBonding=False)
        if mol is None:
            raise ValueError(f"RDKit failed to parse {binder_pdb}")

        # Apply additive charge corrections (preserves RDKit's good defaults)
        n_charge_assigned = _assign_peptide_formal_charges(mol)

        # Re-sanitize after charge fixups (esp. HIS where we flipped ND1 from +1 to 0)
        try:
            Chem.SanitizeMol(mol)
        except Exception as exc:
            print(f"  [info] Post-correction sanitize failed ({exc}); applying Kekulize-only retry")
            try:
                Chem.Kekulize(mol, clearAromaticFlags=False)
            except Exception as exc2:
                print(f"  [warn] Kekulize also failed: {exc2}; writing un-kekulized SDF (may fail openfe load)")

    # Detect disulfide
    has_disulfide = False
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        if a1.GetSymbol() == "S" and a2.GetSymbol() == "S":
            has_disulfide = True
            break

    # Compute formal charge sum
    formal_charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())

    # Set name for OpenFE
    mol.SetProp("_Name", name)

    if not used_openff_polymer:
        # Already wrote via openff path; only write if we built via RDKit fallback
        out_sdf.parent.mkdir(parents=True, exist_ok=True)
        writer = Chem.SDWriter(str(out_sdf))
        try:
            writer.write(mol)
        except Exception as exc:
            print(f"  [error] SDF write failed: {exc}")
            raise
        finally:
            writer.close()

    return {
        "n_atoms": mol.GetNumAtoms(),
        "n_bonds": mol.GetNumBonds(),
        "has_disulfide": has_disulfide,
        "formal_charge_sum": formal_charge,
        "n_formal_charges_assigned": n_charge_assigned,
        "used_openff_polymer_loader": used_openff_polymer,
        "sdf_path": str(out_sdf),
    }


# ---------------------------------------------------------------------------
# Step 3: OpenFE ChemicalSystem + LomapAtomMapping
# ---------------------------------------------------------------------------

def build_openfe_systems_and_mapping(
    sdf_cp4: Path,
    sdf_wt: Path,
    receptor_pdb: Path,
    *,
    out_mapping_json: Path,
) -> Dict[str, Any]:
    """Build OpenFE Cp4/WT SmallMoleculeComponents + LomapAtomMapping.

    Returns dict with mapping score + n_common + n_unique_A + n_unique_B.
    The full openfe Component + Mapping objects are NOT returned (they're
    not pickle-stable across forked workers); they're rebuilt fresh in
    the smoke launcher.
    """
    from openfe import (
        SmallMoleculeComponent,
        ProteinComponent,
        SolventComponent,
        ChemicalSystem,
    )
    from openfe.setup.atom_mapping import LomapAtomMapper

    # Load components
    smc_wt = SmallMoleculeComponent.from_sdf_file(str(sdf_wt))
    smc_cp4 = SmallMoleculeComponent.from_sdf_file(str(sdf_cp4))
    protein = ProteinComponent.from_pdb_file(str(receptor_pdb))
    solvent = SolventComponent()  # default TIP3P + 0.15 M NaCl

    print(f"  WT  binder: {smc_wt.to_rdkit().GetNumAtoms()} atoms")
    print(f"  Cp4 binder: {smc_cp4.to_rdkit().GetNumAtoms()} atoms")
    print(f"  Receptor:   {protein.to_openmm_topology().getNumAtoms()} atoms")

    # LomapAtomMapper for Cp4↔WT
    mapper = LomapAtomMapper(
        time=20,
        threed=True,
        max3d=1.0,
        element_change=False,  # strict: Trp→1-MeW keeps all elements (just +CH3)
    )
    mappings = list(mapper.suggest_mappings(smc_wt, smc_cp4))
    if not mappings:
        raise RuntimeError("LomapAtomMapper produced no mapping")
    mapping = mappings[0]

    # Compute mapping score (LOMAP) — openfe 1.11 namespace
    from openfe.setup import lomap_scorers
    score = lomap_scorers.default_lomap_score(mapping)
    n_common = len(mapping.componentA_to_componentB)
    n_unique_a = mapping.componentA.to_rdkit().GetNumAtoms() - n_common
    n_unique_b = mapping.componentB.to_rdkit().GetNumAtoms() - n_common

    # Serialize mapping for inspection
    out_mapping_json.parent.mkdir(parents=True, exist_ok=True)
    out_mapping_json.write_text(
        json.dumps(
            {
                "componentA_name": "WT",
                "componentB_name": "Cp4",
                "componentA_to_componentB": {
                    str(k): int(v) for k, v in mapping.componentA_to_componentB.items()
                },
                "lomap_score": float(score),
                "n_common_atoms": n_common,
                "n_unique_componentA": n_unique_a,
                "n_unique_componentB": n_unique_b,
            },
            indent=2,
        )
    )

    # Build complex + solvent ChemicalSystems for both endpoints
    sys_complex_wt = ChemicalSystem(
        components={"protein": protein, "ligand": smc_wt, "solvent": solvent},
        name="complex_wt",
    )
    sys_complex_cp4 = ChemicalSystem(
        components={"protein": protein, "ligand": smc_cp4, "solvent": solvent},
        name="complex_cp4",
    )
    sys_solvent_wt = ChemicalSystem(
        components={"ligand": smc_wt, "solvent": solvent},
        name="solvent_wt",
    )
    sys_solvent_cp4 = ChemicalSystem(
        components={"ligand": smc_cp4, "solvent": solvent},
        name="solvent_cp4",
    )

    return {
        "lomap_score": float(score),
        "n_common_atoms": n_common,
        "n_unique_wt": n_unique_a,
        "n_unique_cp4": n_unique_b,
        "smc_wt": smc_wt,
        "smc_cp4": smc_cp4,
        "protein": protein,
        "solvent": solvent,
        "mapping": mapping,
        "sys_complex_wt": sys_complex_wt,
        "sys_complex_cp4": sys_complex_cp4,
        "sys_solvent_wt": sys_solvent_wt,
        "sys_solvent_cp4": sys_solvent_cp4,
    }


# ---------------------------------------------------------------------------
# Step 4: Build smoke protocol settings
# ---------------------------------------------------------------------------

def build_smoke_protocol(
    *,
    equilibration_ps: float = 250.0,
    production_ps: float = 0.0,
    n_lambda: int = 11,
    n_repeats: int = 1,
    temperature_K: float = 298.15,
    seed: int = 7,
):
    """Build RBFE protocol with smoke-test settings.

    250 ps equilibration + 0 ps production = pure smoke (overlap matrix
    can be computed from equilibration MBAR samples). For production
    runs set ``production_ps=5000`` (5 ns) per λ-window.
    """
    from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
    from openff.units import unit

    settings = RelativeHybridTopologyProtocol.default_settings()

    # λ-windows
    settings.lambda_settings.lambda_windows = n_lambda

    # Timing
    settings.simulation_settings.equilibration_length = equilibration_ps * unit.picosecond
    settings.simulation_settings.production_length = production_ps * unit.picosecond
    # Set real-time analysis interval below production length to avoid validator error
    rt_interval_ps = max(min(50.0, equilibration_ps / 2.0), 10.0)
    settings.simulation_settings.real_time_analysis_interval = rt_interval_ps * unit.picosecond
    settings.simulation_settings.real_time_analysis_minimum_time = 0.0 * unit.picosecond
    settings.simulation_settings.minimization_steps = 1000  # smoke: lighter than 5000 default

    # Sampler
    settings.simulation_settings.sampler_method = "repex"
    settings.simulation_settings.n_replicas = n_lambda

    # Repeats
    settings.protocol_repeats = n_repeats

    # Temperature
    settings.thermo_settings.temperature = temperature_K * unit.kelvin

    # Engine (CUDA, already default 'cuda' lowercased)
    settings.engine_settings.compute_platform = "cuda"

    # Charge backend: NAGL (graph-NN AM1-BCC). Sub-second on 200+ atom peptide
    # vs ambertools sqm which times out (~10 min on Cp4-sized molecule,
    # validated 2026-05-31). NAGL openff-gnn-am1bcc-1.0.0.pt produces
    # total_charge = 0.000 on amber14SB-protonated WT compstatin.
    settings.partial_charge_settings.partial_charge_method = "nagl"
    settings.partial_charge_settings.nagl_model = "openff-gnn-am1bcc-1.0.0.pt"

    return settings


# ---------------------------------------------------------------------------
# Step 5: Run smoke (bound + free legs)
# ---------------------------------------------------------------------------

def run_smoke_leg(
    *,
    leg: str,
    sys_stateA,
    sys_stateB,
    mapping,
    settings,
    out_dir: Path,
    cuda_device: int = 0,
) -> Dict[str, Any]:
    """Run one RBFE leg (complex or solvent) and collect MBAR overlap.

    Uses openfe's ``Transformation`` + ``execute`` pattern. We run it
    in-process (no openfe_quickrun CLI dispatch) so we can capture
    intermediate dataframe + free energy convergence in real-time.
    """
    from openfe import Transformation
    from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol
    from gufe.protocols import execute_DAG
    from gufe import AlchemicalNetwork

    out_dir.mkdir(parents=True, exist_ok=True)

    protocol = RelativeHybridTopologyProtocol(settings=settings)
    transformation = Transformation(
        stateA=sys_stateA,
        stateB=sys_stateB,
        mapping=mapping,
        protocol=protocol,
        name=f"{leg}_wt_to_cp4",
    )

    # Persist transformation JSON for inspection
    (out_dir / "transformation.json").write_text(
        json.dumps({"leg": leg, "stateA": sys_stateA.name, "stateB": sys_stateB.name}, indent=2)
    )

    # Set CUDA device
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", str(cuda_device))

    # Execute
    t0 = time.time()
    dag = transformation.create()
    scratch_dir = out_dir / "scratch"
    scratch_dir.mkdir(parents=True, exist_ok=True)
    dag_result = execute_DAG(dag, shared_basedir=out_dir, scratch_basedir=scratch_dir, n_retries=0)
    wall_s = time.time() - t0

    # Gather results
    proto_result = protocol.gather([dag_result])
    estimate = proto_result.get_estimate()
    uncertainty = proto_result.get_uncertainty()

    # MBAR overlap matrix
    overlap_matrix = None
    try:
        overlap_matrix = proto_result.get_overlap_matrix()
    except Exception as exc:
        print(f"  [warn] overlap matrix extraction failed: {exc}")

    if overlap_matrix is not None:
        # off-diagonal min
        n = overlap_matrix.shape[0]
        offdiag = []
        for i in range(n):
            for j in range(n):
                if i != j:
                    offdiag.append(float(overlap_matrix[i, j]))
        offdiag_min = float(np.min(offdiag)) if offdiag else 0.0
        offdiag_mean = float(np.mean(offdiag)) if offdiag else 0.0
    else:
        offdiag_min = -1.0
        offdiag_mean = -1.0

    result = {
        "leg": leg,
        "wall_s": wall_s,
        "delta_g_kcal_mol": float(estimate.to("kilocalorie/mol").magnitude) if estimate is not None else None,
        "uncertainty_kcal_mol": float(uncertainty.to("kilocalorie/mol").magnitude) if uncertainty is not None else None,
        "overlap_min_offdiag": offdiag_min,
        "overlap_mean_offdiag": offdiag_mean,
        "n_lambda": settings.lambda_settings.lambda_windows,
        "n_replicas": settings.simulation_settings.n_replicas,
    }
    (out_dir / "smoke_result.json").write_text(json.dumps(result, indent=2))
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def cmd_build(args):
    """Phase 2: build SDFs + receptor + mapping; no MD."""
    seed_dir = OUTPUT_ROOT / args.seed_tag
    seed_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: split complex PDBs
    print(f"[1/3] Split Cp4 complex → receptor + binder")
    rec_cp4_pdb = seed_dir / "receptor_from_cp4.pdb"
    bnd_cp4_pdb = seed_dir / "binder_cp4_raw.pdb"
    n_rec_cp4, n_bnd_cp4 = split_complex_to_receptor_binder(CP4_COMPLEX_PDB, rec_cp4_pdb, bnd_cp4_pdb)
    print(f"  Cp4: receptor {n_rec_cp4} atoms, binder {n_bnd_cp4} atoms")

    print(f"[1/3] Split WT complex → receptor + binder")
    rec_wt_pdb = seed_dir / "receptor_from_wt.pdb"
    bnd_wt_pdb = seed_dir / "binder_wt_raw.pdb"
    n_rec_wt, n_bnd_wt = split_complex_to_receptor_binder(WT_COMPLEX_PDB, rec_wt_pdb, bnd_wt_pdb)
    print(f"  WT:  receptor {n_rec_wt} atoms, binder {n_bnd_wt} atoms")

    # Use Cp4-derived receptor as canonical (WT one should be identical chain A)
    receptor_pdb = seed_dir / "receptor.pdb"
    shutil.copy(rec_cp4_pdb, receptor_pdb)

    # Step 1b: protonate binders (Cp4 MTR needs MTR_hydrogens.xml; WT is consistency-only)
    print(f"[1b/3] Protonate binders (Cp4 MTR + WT amber14SB)")
    bnd_cp4_pdb_h = seed_dir / "binder_cp4_protonated.pdb"
    bnd_wt_pdb_h = seed_dir / "binder_wt_protonated.pdb"
    n_cp4_h = protonate_binder_pdb(bnd_cp4_pdb, bnd_cp4_pdb_h, is_cp4=True)
    n_wt_h = protonate_binder_pdb(bnd_wt_pdb, bnd_wt_pdb_h, is_cp4=False)
    print(f"  Cp4 protonated: {n_cp4_h} atoms (raw was {n_bnd_cp4})")
    print(f"  WT  protonated: {n_wt_h} atoms (raw was {n_bnd_wt})")

    # Step 2: binder PDB → SDF (WT uses openff polymer loader; Cp4 uses RDKit fallback for MTR)
    print(f"[2/3] Binder PDB → SDF (WT via openff polymer; Cp4 via RDKit + charge correction)")
    sdf_cp4 = seed_dir / "binder_cp4.sdf"
    sdf_wt = seed_dir / "binder_wt.sdf"
    meta_wt = binder_pdb_to_sdf(bnd_wt_pdb_h, sdf_wt, name="WT", is_canonical_only=True)
    meta_cp4 = binder_pdb_to_sdf(bnd_cp4_pdb_h, sdf_cp4, name="Cp4", is_canonical_only=False)
    print(f"  Cp4 SDF: {meta_cp4['n_atoms']} atoms, {meta_cp4['n_bonds']} bonds, "
          f"disulfide={meta_cp4['has_disulfide']}, formal_charge={meta_cp4['formal_charge_sum']}")
    print(f"  WT  SDF: {meta_wt['n_atoms']} atoms, {meta_wt['n_bonds']} bonds, "
          f"disulfide={meta_wt['has_disulfide']}, formal_charge={meta_wt['formal_charge_sum']}")

    # Step 3: OpenFE systems + LomapAtomMapping
    print(f"[3/3] OpenFE systems + Lomap atom mapping")
    mapping_json = seed_dir / "mapping.json"
    build_info = build_openfe_systems_and_mapping(
        sdf_cp4=sdf_cp4,
        sdf_wt=sdf_wt,
        receptor_pdb=receptor_pdb,
        out_mapping_json=mapping_json,
    )
    print(f"  LOMAP score:   {build_info['lomap_score']:.3f} (>= 0.5 = favorable)")
    print(f"  Common atoms:  {build_info['n_common_atoms']}")
    print(f"  Unique on WT:  {build_info['n_unique_wt']}")
    print(f"  Unique on Cp4: {build_info['n_unique_cp4']}")

    build_meta = {
        "seed_tag": args.seed_tag,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "cp4_complex_pdb": str(CP4_COMPLEX_PDB),
        "wt_complex_pdb": str(WT_COMPLEX_PDB),
        "receptor_pdb": str(receptor_pdb),
        "sdf_cp4": str(sdf_cp4),
        "sdf_wt": str(sdf_wt),
        "binder_meta_cp4": meta_cp4,
        "binder_meta_wt": meta_wt,
        "lomap_score": build_info["lomap_score"],
        "n_common_atoms": build_info["n_common_atoms"],
        "n_unique_wt": build_info["n_unique_wt"],
        "n_unique_cp4": build_info["n_unique_cp4"],
        "openfe_version": _get_openfe_version(),
        "openmm_version": _get_openmm_version(),
    }
    (seed_dir / "build_metadata.json").write_text(json.dumps(build_meta, indent=2))
    print(f"\nBuild PASS. Metadata → {seed_dir / 'build_metadata.json'}")
    return 0


def cmd_smoke(args):
    """Phase 3: run smoke MD on bound + free legs."""
    seed_dir = OUTPUT_ROOT / args.seed_tag
    if not (seed_dir / "build_metadata.json").exists():
        print(f"[fatal] Build artifacts missing in {seed_dir}. Run --build first.")
        return 1

    sdf_cp4 = seed_dir / "binder_cp4.sdf"
    sdf_wt = seed_dir / "binder_wt.sdf"
    receptor_pdb = seed_dir / "receptor.pdb"

    print(f"[1/3] Rebuild OpenFE systems + mapping (fresh, for fork-safe execution)")
    build_info = build_openfe_systems_and_mapping(
        sdf_cp4=sdf_cp4,
        sdf_wt=sdf_wt,
        receptor_pdb=receptor_pdb,
        out_mapping_json=seed_dir / "mapping.json",
    )

    print(f"[2/3] Build smoke protocol settings")
    settings = build_smoke_protocol(
        equilibration_ps=args.equilibration_ps,
        production_ps=args.production_ps,
        n_lambda=args.n_lambda,
        n_repeats=args.n_repeats,
        seed=int(args.seed_tag.lstrip("s")),
    )
    # Persist settings dict
    try:
        import dataclasses
        settings_dict = {
            "equilibration_ps": float(args.equilibration_ps),
            "production_ps": float(args.production_ps),
            "n_lambda": int(args.n_lambda),
            "n_repeats": int(args.n_repeats),
            "sampler": "repex",
            "temperature_K": 298.15,
            "platform": "CUDA",
        }
        (seed_dir / "settings.json").write_text(json.dumps(settings_dict, indent=2))
    except Exception as e:
        print(f"  [warn] settings persistence failed: {e}")

    results = {}
    legs_to_run = []
    if "complex" in args.legs:
        legs_to_run.append("complex")
    if "solvent" in args.legs:
        legs_to_run.append("solvent")

    for leg in legs_to_run:
        print(f"\n[3/3] Run leg: {leg}")
        leg_dir = seed_dir / f"smoke_{leg}"
        if leg == "complex":
            sys_a = build_info["sys_complex_wt"]
            sys_b = build_info["sys_complex_cp4"]
        else:
            sys_a = build_info["sys_solvent_wt"]
            sys_b = build_info["sys_solvent_cp4"]
        try:
            result = run_smoke_leg(
                leg=leg,
                sys_stateA=sys_a,
                sys_stateB=sys_b,
                mapping=build_info["mapping"],
                settings=settings,
                out_dir=leg_dir,
                cuda_device=args.cuda_device,
            )
            results[leg] = result
            print(f"  {leg}: ΔG={result['delta_g_kcal_mol']:.3f} ± {result['uncertainty_kcal_mol']:.3f} kcal/mol")
            print(f"  {leg}: overlap min off-diag = {result['overlap_min_offdiag']:.3f}, mean = {result['overlap_mean_offdiag']:.3f}")
            print(f"  {leg}: wall = {result['wall_s']:.1f} s")
        except Exception as exc:
            print(f"  [error] leg {leg} failed: {exc}")
            import traceback
            traceback.print_exc()
            results[leg] = {"error": str(exc), "leg": leg}

    # Aggregate smoke result + C1 gate evaluation
    gate_pass = _evaluate_c1_gate(results)
    summary = {
        "seed_tag": args.seed_tag,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "legs": results,
        "c1_gate_pass": gate_pass["pass"],
        "c1_gate_details": gate_pass,
    }
    (seed_dir / "smoke_results.json").write_text(json.dumps(summary, indent=2))
    print(f"\n========= SMOKE GATE C1 =========")
    print(f"  PASS: {gate_pass['pass']}")
    for k, v in gate_pass.items():
        if k != "pass":
            print(f"  {k}: {v}")
    return 0 if gate_pass["pass"] else 2


def _evaluate_c1_gate(results: Dict[str, Any]) -> Dict[str, Any]:
    """C1 gate: MBAR overlap > 0.05 ALL pairs + hysteresis < 1.0 kcal/mol.

    Hysteresis is computed across legs (complex - solvent ΔΔG_bind cycle),
    but for a 1-rep smoke we use the uncertainty as proxy. Full
    forward/reverse hysteresis requires 2 repeats with reversed
    initialization; deferred to production.
    """
    details = {}
    overlap_pass = True
    for leg, r in results.items():
        if not isinstance(r, dict) or "overlap_min_offdiag" not in r:
            details[f"{leg}_overlap"] = None
            overlap_pass = False
            continue
        if r["overlap_min_offdiag"] < 0.05:
            overlap_pass = False
        details[f"{leg}_overlap_min"] = r["overlap_min_offdiag"]

    # ΔΔG_bind = ΔG_complex - ΔG_solvent
    if "complex" in results and "solvent" in results:
        cx = results["complex"]
        sv = results["solvent"]
        if (isinstance(cx, dict) and isinstance(sv, dict)
            and "delta_g_kcal_mol" in cx and cx["delta_g_kcal_mol"] is not None
            and "delta_g_kcal_mol" in sv and sv["delta_g_kcal_mol"] is not None):
            ddg = cx["delta_g_kcal_mol"] - sv["delta_g_kcal_mol"]
            details["ddg_bind_kcal_mol"] = ddg
            # Uncertainty (propagated)
            unc = ((cx.get("uncertainty_kcal_mol") or 0) ** 2 + (sv.get("uncertainty_kcal_mol") or 0) ** 2) ** 0.5
            details["ddg_uncertainty_kcal_mol"] = unc
            hysteresis_proxy = unc * 2.0  # 2σ as forward/reverse proxy
            details["hysteresis_proxy_kcal_mol"] = hysteresis_proxy
            hysteresis_pass = hysteresis_proxy < 1.0
        else:
            hysteresis_pass = False
            details["hysteresis_proxy_kcal_mol"] = None
    else:
        hysteresis_pass = None  # only one leg run
        details["hysteresis_proxy_kcal_mol"] = "single-leg"

    overall_pass = overlap_pass and (hysteresis_pass is not False)
    details["pass"] = overall_pass
    details["overlap_gate_threshold"] = 0.05
    details["hysteresis_gate_threshold"] = 1.0
    return details


def _get_openfe_version() -> str:
    try:
        import openfe
        return openfe.__version__
    except Exception:
        return "unknown"


def _get_openmm_version() -> str:
    try:
        import openmm
        return openmm.__version__
    except Exception:
        return "unknown"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest="mode", required=True)

    p_build = subparsers.add_parser("build", help="Phase 2: build SDFs + mapping (no MD)")
    p_build.add_argument("--seed-tag", default="s7", help="Cohort seed tag (s7|s19|s23|s101|s127|s163|s199|s251)")

    p_smoke = subparsers.add_parser("smoke", help="Phase 3: run smoke MD on bound + free legs")
    p_smoke.add_argument("--seed-tag", default="s7")
    p_smoke.add_argument("--equilibration-ps", type=float, default=250.0)
    p_smoke.add_argument("--production-ps", type=float, default=0.0)
    p_smoke.add_argument("--n-lambda", type=int, default=11)
    p_smoke.add_argument("--n-repeats", type=int, default=1)
    p_smoke.add_argument("--cuda-device", type=int, default=0)
    p_smoke.add_argument("--legs", default="complex,solvent",
                         help="Comma-separated list: complex,solvent (default both)")

    p_both = subparsers.add_parser("both", help="Build + smoke in one shot")
    p_both.add_argument("--seed-tag", default="s7")
    p_both.add_argument("--equilibration-ps", type=float, default=250.0)
    p_both.add_argument("--production-ps", type=float, default=0.0)
    p_both.add_argument("--n-lambda", type=int, default=11)
    p_both.add_argument("--n-repeats", type=int, default=1)
    p_both.add_argument("--cuda-device", type=int, default=0)
    p_both.add_argument("--legs", default="complex,solvent")

    args = parser.parse_args()

    if args.mode == "build":
        return cmd_build(args)
    if args.mode == "smoke":
        return cmd_smoke(args)
    if args.mode == "both":
        rc = cmd_build(args)
        if rc != 0:
            return rc
        return cmd_smoke(args)
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
