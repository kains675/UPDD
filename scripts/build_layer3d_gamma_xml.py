#!/usr/bin/env python
"""
Build params/MTR_gaff2_layer3d_gamma.xml — v0.9 Layer 3D γ option XML.

Scientific design:
  - Bonded parameters: INHERIT all amber14SB Trp bonded parameters (bonds,
    angles, dihedrals, impropers) via amber14SB-compatible `protein-*` atom
    types. The amber14SB force-field loaded at runtime supplies the canonical
    Trp bonded definitions for 24 of the 27 atoms.
  - Methyl group junction: the new NE1-CM bond + 3 angles (HC-CT-NA,
    CT-NA-CW, CT-NA-CN) are added as explicit class-typed supplemental
    entries; values match the amber-class methyl-on-aromatic-N parameters
    already in `MTR_gaff2_resp.xml` (provenance: Khoury 2014 OMW XML
    bonded section, which itself maps onto ff14SB CT-NA family).
  - Methyl-junction ring-planarity improper: amber14SB Trp has the
    wildcard improper `protein-NA *-*-protein-H` (k=4.184, period=2,
    phase=π) which enforces NE1 sp2 planarity with the indole ring.
    With H → CM substitution this improper no longer matches; we add
    the equivalent `protein-NA *-*-protein-CT` improper with the same
    k value (scientifically defensible — the role is identical, only the
    substituent atom type changes).
  - Partial charges: Khoury 2014 OMW RESP-A2 (full RESP, not the v0.8 β
    hybrid sidechain-only override). All 27 atoms get Khoury values.
  - Σq = 0 enforced via least-significant-bit redistribution across all
    27 atoms.

This is the v0.9 PARALLEL FALLBACK XML.

Inputs:
  - `params/MTR_gaff2_resp.xml` (existing full-RESP XML; γ uses identical
    atom types + charges + 24 amber-class supplemental terms; γ adds the
    NA-*-*-CT improper for indole planarity).

Outputs:
  - `params/MTR_gaff2_layer3d_gamma.xml`
  - `outputs/analysis/d_gamma_xml_build_<TS>/MTR_gaff2_layer3d_gamma.xml`
  - `outputs/analysis/d_gamma_xml_build_<TS>/manifest.json`
  - `outputs/analysis/d_gamma_xml_build_<TS>/build_log.md`
  - `outputs/analysis/d_gamma_xml_build_<TS>/validation_report.md`

References:
  - Khoury 2014 ACS Synth Biol DOI 10.1021/sb400168u
  - Maier 2015 J. Chem. Theory Comput. 11, 3696 DOI 10.1021/acs.jctc.5b00255
    (amber14SB Trp canonical bonded parameter source)
  - amber14SB OpenMM XML: openmm/.../amber14/protein.ff14SB.xml (Trp residue
    L892-944)
  - Limitations roadmap §L1 (b'/b''/c split)
"""
from __future__ import annotations

import json
import sys
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

PROJ = Path(__file__).resolve().parent.parent
TEMPLATE_XML = PROJ / "params" / "MTR_gaff2_resp.xml"
OUT_XML = PROJ / "params" / "MTR_gaff2_layer3d_gamma.xml"
OUT_MANIFEST = PROJ / "params" / "MTR_gaff2_layer3d_gamma_manifest.json"

# Build timestamp marker dir for provenance copies
ANALYSIS_DIR_PARENT = PROJ / "outputs" / "analysis"


KHOURY_OMW_CHARGES: Dict[str, float] = {
    "N":   -0.280943,
    "H":    0.241023,
    "CA":  -0.124259,
    "HA":   0.120930,
    "C":    0.575115,
    "O":   -0.571038,
    "CB":  -0.157536,
    "HB2":  0.104569,
    "HB3":  0.104569,
    "CG":  -0.126422,
    "CD2":  0.083375,
    "CE2":  0.119586,
    "NE1":  0.077670,
    "CD1": -0.252655,
    "HD1":  0.204947,
    "CM":  -0.239062,
    "HM1":  0.112200,
    "HM2":  0.112200,
    "HM3":  0.112200,
    "CZ2": -0.241414,
    "HZ2":  0.146305,
    "CH2": -0.139177,
    "HH2":  0.136546,
    "CZ3": -0.206367,
    "HZ3":  0.143250,
    "CE3": -0.211658,
    "HE3":  0.156047,
}

# γ-option NE1 ring-planarity improper substitution
# amber14SB Trp ring-planarity improper is:
#     <Improper k1="4.184" periodicity1="2" phase1="3.141592653589793"
#               type1="protein-NA" type2="" type3="" type4="protein-H"/>
# With H → CM (protein-CT) substitution we need the equivalent:
#     type1="protein-NA" type2="" type3="" type4="protein-CT"
# k1 = 4.184 kcal/mol (= 17.5728 kJ/mol; identical to H variant)
# period = 2, phase = π (sp2 planarity)
NE1_CM_PLANARITY_IMPROPER: Dict[str, str] = {
    "type1": "protein-NA",
    "type2": "",
    "type3": "",
    "type4": "protein-CT",
    "periodicity1": "2",
    "phase1": "3.141592653589793",
    "k1": "4.184",
}


def find_residue(root: ET.Element, name: str) -> ET.Element:
    for res in root.findall(".//Residue"):
        if res.get("name") == name:
            return res
    raise RuntimeError(f"Residue {name!r} not found")


def apply_charges(residue: ET.Element) -> Dict[str, float]:
    """Override all atom charges with Khoury OMW values. Returns final dict."""
    overridden: Dict[str, float] = {}
    missing: List[str] = []
    for atom in residue.findall("Atom"):
        name = atom.get("name", "")
        if name in KHOURY_OMW_CHARGES:
            charge = KHOURY_OMW_CHARGES[name]
            atom.set("charge", f"{charge:.6f}")
            overridden[name] = charge
        else:
            missing.append(name)
    if missing:
        raise RuntimeError(f"Khoury OMW does not cover atoms: {missing}")
    return overridden


def neutralize_charges(residue: ET.Element) -> float:
    """Distribute residual Σq across all atoms to ensure exact 0."""
    atoms = list(residue.findall("Atom"))
    sigma = sum(float(a.get("charge", "0.0")) for a in atoms)
    if abs(sigma) > 1e-12:
        n_atoms = len(atoms)
        per_atom = sigma / n_atoms
        for atom in atoms:
            q = float(atom.get("charge", "0.0")) - per_atom
            atom.set("charge", f"{q:.6f}")
    return sum(float(a.get("charge", "0.0")) for a in atoms)


def add_planarity_improper(root: ET.Element) -> ET.Element:
    """Append the protein-NA *-*-protein-CT improper (γ-specific).

    Returns the inserted Improper element for manifest reporting.
    """
    torsion_force = root.find(".//PeriodicTorsionForce")
    if torsion_force is None:
        raise RuntimeError("PeriodicTorsionForce section not found in XML")

    # Avoid duplicate insertion if already present
    for impr in torsion_force.findall("Improper"):
        if (impr.get("type1") == "protein-NA"
                and impr.get("type4") == "protein-CT"
                and not impr.get("type2")
                and not impr.get("type3")):
            return impr

    improper = ET.SubElement(torsion_force, "Improper")
    for key, value in NE1_CM_PLANARITY_IMPROPER.items():
        improper.set(key, value)
    return improper


def collect_atom_types(residue: ET.Element) -> Dict[str, Dict[str, int]]:
    """Tally atom types in the MTR residue (for validation report)."""
    type_counts: Dict[str, int] = {}
    name_to_type: Dict[str, str] = {}
    for atom in residue.findall("Atom"):
        name = atom.get("name", "")
        atype = atom.get("type", "")
        name_to_type[name] = atype
        type_counts[atype] = type_counts.get(atype, 0) + 1
    return {"counts": type_counts, "name_to_type": name_to_type}


def parse_template_bonded(root: ET.Element) -> Dict[str, int]:
    """Count bonded terms in the template (pre-γ) for validation report."""
    counts: Dict[str, int] = {}
    counts["HarmonicBondForce_class"] = len(
        [b for b in root.findall(".//HarmonicBondForce/Bond")
         if b.get("class1") is not None])
    counts["HarmonicBondForce_type"] = len(
        [b for b in root.findall(".//HarmonicBondForce/Bond")
         if b.get("type1") is not None])
    counts["HarmonicAngleForce_class"] = len(
        [a for a in root.findall(".//HarmonicAngleForce/Angle")
         if a.get("class1") is not None])
    counts["HarmonicAngleForce_type"] = len(
        [a for a in root.findall(".//HarmonicAngleForce/Angle")
         if a.get("type1") is not None])
    counts["PeriodicTorsionForce_Proper_class"] = len(
        [p for p in root.findall(".//PeriodicTorsionForce/Proper")
         if p.get("class1") is not None])
    counts["PeriodicTorsionForce_Improper_class"] = len(
        [p for p in root.findall(".//PeriodicTorsionForce/Improper")
         if p.get("class1") is not None])
    counts["PeriodicTorsionForce_Improper_type"] = len(
        [p for p in root.findall(".//PeriodicTorsionForce/Improper")
         if p.get("type1") is not None])
    return counts


def write_xml(tree: ET.ElementTree, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    tree.write(target, encoding="utf-8", xml_declaration=False)


def build_manifest(
    final_charges: Dict[str, float],
    sigma_q: float,
    type_inventory: Dict[str, Any],
    bonded_counts_before: Dict[str, int],
    bonded_counts_after: Dict[str, int],
    template_path: Path,
) -> Dict[str, Any]:
    return {
        "schema": "mtr_layer3d_gamma_manifest/0.1",
        "generated_by": "scripts/build_layer3d_gamma_xml.py",
        "build_timestamp_utc": datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "rationale": (
            "v0.9 Layer 3D γ-option XML: amber14SB Trp bonded parameter "
            "inheritance + Khoury OMW full RESP-A2 partial charges + GAFF2 "
            "methyl-junction supplemental terms. PARALLEL FALLBACK option "
            "for v0.9 pilot. Ranking-regime production candidate when "
            "relaxed-scan β option (b' in limitations_roadmap §L1) misses "
            "acceptance criterion."
        ),
        "design_principles": {
            "bonded_inheritance": (
                "All 24 amber-equivalent atoms use protein-* types (protein-N, "
                "protein-H, protein-CX, protein-H1, protein-CT, protein-HC, "
                "protein-C*, protein-CB, protein-CW, protein-H4, protein-NA, "
                "protein-CN, protein-CA, protein-HA, protein-C, protein-O). "
                "amber14SB.xml loaded at runtime supplies all canonical Trp "
                "bond/angle/torsion parameters."
            ),
            "methyl_junction_bonded": (
                "CM uses protein-CT (CT–CT, CT–CW, CT–CN inherited from "
                "amber14SB). HM1/2/3 use protein-HC. The CT–NA bond + 3 "
                "junction angles (HC–CT–NA, CT–NA–CW, CT–NA–CN) are added as "
                "explicit class-typed supplemental entries to handle the "
                "amber14SB gap (standard Trp has H–NA, not CT–NA)."
            ),
            "ring_planarity_improper": (
                "amber14SB has wildcard improper protein-NA *–*–protein-H "
                "(k=4.184 kcal/mol, period=2, phase=π) enforcing NE1 sp2 "
                "planarity. With H→CM substitution this improper no longer "
                "matches. γ XML adds the equivalent protein-NA *–*–protein-CT "
                "improper with the same k value (the planarity role is "
                "identical; only the substituent atom type changes)."
            ),
            "partial_charges": (
                "Full Khoury OMW RESP-A2 fit applied to all 27 atoms. NOT "
                "the v0.8 β hybrid (which preserved N/H/NE1 baseline). "
                "γ vs β distinction: β = charges-only mitigation of NE1 sign "
                "flip; γ = full Khoury charges + bonded inheritance shift "
                "from Khoury's ff03-derived bonded to amber14SB ff14SB Trp."
            ),
        },
        "scientific_provenance": {
            "khoury_2014": {
                "doi": "10.1021/sb400168u",
                "block": "ffncaa.in OMW (1-methyl-tryptophan)",
                "block_line": "778-830",
                "convention": "internal-residue, ff03-compatible, HF/6-31G* ESP, RESP-A2 fit",
            },
            "amber14sb": {
                "reference": "Maier 2015 J. Chem. Theory Comput. 11, 3696",
                "doi": "10.1021/acs.jctc.5b00255",
                "source_file": "openmm/.../amber14/protein.ff14SB.xml",
                "trp_residue_lines": "892-944",
            },
            "limitations_roadmap": {
                "path": "plan/limitations_roadmap_20260519.md",
                "section": "§L1 path (c) γ parallel fallback",
            },
        },
        "atom_name_mapping_khoury_to_updd": {
            "CZ1": "CM",
            "HZ11": "HM1",
            "HZ12": "HM2",
            "HZ13": "HM3",
        },
        "template_source": {
            "path": str(template_path.relative_to(PROJ)),
            "convention": "amber14SB-typed (protein-*), Khoury OMW RESP-A2 charges",
            "note": (
                "Existing MTR_gaff2_resp.xml already implements the protein-* "
                "type pattern. γ XML differs from it ONLY by the added "
                "NE1 ring-planarity improper (protein-NA *–*–protein-CT)."
            ),
        },
        "n_atoms": len(final_charges),
        "n_atoms_amber_inherited": 23,  # 27 total - 4 methyl-specific (CM HM1 HM2 HM3)
        "n_atoms_methyl_junction": 4,
        "final_sigma_q": round(sigma_q, 9),
        "atom_type_inventory": type_inventory,
        "bonded_term_counts": {
            "before_gamma_additions": bonded_counts_before,
            "after_gamma_additions": bonded_counts_after,
            "gamma_additions": {
                "NE1_CM_planarity_improper": NE1_CM_PLANARITY_IMPROPER,
            },
        },
        "final_charges_per_atom": {k: round(v, 6) for k, v in final_charges.items()},
        "validation_gates": {
            "sigma_q_neutral": abs(sigma_q) < 1e-5,
            "all_27_atoms_present": len(final_charges) == 27,
            "all_atoms_have_khoury_charge": (
                sorted(final_charges) == sorted(KHOURY_OMW_CHARGES)
            ),
        },
    }


def main() -> int:
    if not TEMPLATE_XML.exists():
        print(f"[ERROR] Template not found: {TEMPLATE_XML}", file=sys.stderr)
        return 1

    print(f"[info] Template XML: {TEMPLATE_XML}")
    tree = ET.parse(TEMPLATE_XML)
    root = tree.getroot()

    residue = find_residue(root, "MTR")

    # Snapshot pre-modification state
    bonded_before = parse_template_bonded(root)
    print(f"[info] Template bonded term counts: {bonded_before}")

    # 1. Apply Khoury charges (idempotent for resp XML, but explicit)
    overridden = apply_charges(residue)
    print(f"[info] Khoury OMW charges applied to {len(overridden)} atoms")

    # 2. Normalize Σq (resp XML may already be neutral; this re-checks)
    sigma_q_final = neutralize_charges(residue)
    print(f"[info] Final Σq = {sigma_q_final:+.3e}")
    if abs(sigma_q_final) > 1e-5:
        print(f"[ERROR] Σq normalization failed: {sigma_q_final}", file=sys.stderr)
        return 1

    # 3. Add γ-specific NE1 ring-planarity improper for CM substitution
    impr = add_planarity_improper(root)
    print(f"[info] Added improper: NA-*-*-CT k={impr.get('k1')} (γ-specific)")

    # Snapshot post-modification state
    bonded_after = parse_template_bonded(root)
    print(f"[info] Final bonded term counts: {bonded_after}")

    # 4. Atom type inventory
    type_inventory = collect_atom_types(residue)
    print(f"[info] Atom type counts: {type_inventory['counts']}")

    # 5. Write XML
    write_xml(tree, OUT_XML)
    print(f"[OK] Wrote: {OUT_XML.relative_to(PROJ)}")

    # 6. Build + write manifest
    final_charges = {a.get("name", ""): float(a.get("charge", "0.0"))
                     for a in residue.findall("Atom")}
    manifest = build_manifest(
        final_charges=final_charges,
        sigma_q=sigma_q_final,
        type_inventory=type_inventory,
        bonded_counts_before=bonded_before,
        bonded_counts_after=bonded_after,
        template_path=TEMPLATE_XML,
    )
    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with OUT_MANIFEST.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"[OK] Wrote: {OUT_MANIFEST.relative_to(PROJ)}")

    # 7. Mirror copies to analysis output dir if a TS path was provided
    ts_dirs = sorted(ANALYSIS_DIR_PARENT.glob("d_gamma_xml_build_*"))
    if ts_dirs:
        latest = ts_dirs[-1]
        for src in (OUT_XML, OUT_MANIFEST):
            dst = latest / src.name.replace("_layer3d_gamma_manifest", "manifest")
            if src.name.endswith("manifest.json"):
                dst = latest / "manifest.json"
            else:
                dst = latest / src.name
            dst.write_bytes(src.read_bytes())
            print(f"[OK] Mirror: {dst.relative_to(PROJ)}")

    print()
    print("=== Validation summary ===")
    print(f"Atoms covered: {len(final_charges)}/27")
    print(f"Σq          = {sigma_q_final:+.2e}  (gate: |Σq| < 1e-5)")
    print(f"Atom types  = {sorted(set(type_inventory['name_to_type'].values()))}")
    print(f"γ additions = NE1 ring planarity (protein-NA *-*-protein-CT)")
    print()
    print("Next: smoke MD with this XML (1×5ns, 4 seeds).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
