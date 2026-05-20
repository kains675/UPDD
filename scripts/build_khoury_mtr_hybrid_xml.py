#!/usr/bin/env python
"""
Build params/MTR_gaff2_hybrid.xml — sidechain-only Khoury OMW RESP-A2 override
with backbone N, backbone amide H, and indole NE1 retained from the
amber14SB-patched baseline. Mitigates the MD-propagation-level junction
limitation identified in the v0.8 full-RESP pilot (NE1 sign flip from
-0.342 → +0.078 caused HD1 +/+ Coulomb repulsion → ring planarity
destabilization → NaN cascade in 2/8 seeds).

Preserved atoms (baseline amber14SB-patched values):
  - N    = -0.4157   (backbone N, ff03/ff14SB Trp anchor)
  - H    = +0.2719   (backbone amide H, paired with N)
  - NE1  = -0.3418   (indole-N, paired with HD1; sign-flip preventer)

Overridden atoms (24 atoms with Khoury OMW RESP-A2):
  - sidechain carbons + hydrogens (CA, HA, CB, HB2/3, CG, CD1, CD2, CE2/3,
    NE1 [excluded], CZ2/3, CH2, HD1 [paired with NE1, but retained Khoury per
    default scope], HE3, HZ2/3, HH2)
  - methyl group (CM, HM1/2/3) — the MTR functional substitution
  - C, O backbone carbonyl

Residual Σq distributed across the 24 overridden atoms only (baseline
atoms left untouched at exact amber14SB values).

UPDD ↔ Khoury atom-name mapping (same as full-RESP):
  - Khoury CZ1 → UPDD CM
  - Khoury HZ11/HZ12/HZ13 → UPDD HM1/HM2/HM3

Inputs:
  - Template (baseline): outputs/2QKI_Cp4_calib_s7/params/MTR_gaff2.xml

Outputs:
  - params/MTR_gaff2_hybrid.xml
  - params/MTR_gaff2_hybrid_manifest.json (provenance + Δq per atom)
"""
from __future__ import annotations

import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
TEMPLATE_XML = PROJ / "outputs/2QKI_Cp4_calib_s7/params/MTR_gaff2.xml"
OUT_XML = PROJ / "params/MTR_gaff2_hybrid.xml"
OUT_MANIFEST = PROJ / "params/MTR_gaff2_hybrid_manifest.json"

KHOURY_OMW_CHARGES = {
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

# Atoms preserved from baseline (NOT overridden by Khoury)
PRESERVE_BASELINE = {"N", "H", "NE1"}


def main() -> int:
    if not TEMPLATE_XML.exists():
        print(f"[ERROR] Template XML not found: {TEMPLATE_XML}", file=sys.stderr)
        return 1

    tree = ET.parse(TEMPLATE_XML)
    root = tree.getroot()

    residue = None
    for res in root.findall(".//Residue"):
        if res.get("name") == "MTR":
            residue = res
            break
    if residue is None:
        print("[ERROR] MTR Residue not found in template XML", file=sys.stderr)
        return 1

    # Snapshot baseline
    template_charges = {}
    for atom in residue.findall("Atom"):
        template_charges[atom.get("name")] = float(atom.get("charge", "0.0"))

    template_sigma_q = sum(template_charges.values())
    print(f"[info] Template baseline Σq = {template_sigma_q:+.6f}")
    if abs(template_sigma_q) > 1e-5:
        print(f"[WARN] Template not neutral: {template_sigma_q}")

    # Apply hybrid: preserve baseline for N/H/NE1, override others with Khoury
    overrides = 0
    preserved = 0
    missing = []
    overridden_atom_names = []
    for atom in residue.findall("Atom"):
        name = atom.get("name")
        if name in PRESERVE_BASELINE:
            preserved += 1
        elif name in KHOURY_OMW_CHARGES:
            atom.set("charge", f"{KHOURY_OMW_CHARGES[name]:.6f}")
            overrides += 1
            overridden_atom_names.append(name)
        else:
            missing.append(name)

    new_sigma_q = sum(float(a.get("charge")) for a in residue.findall("Atom"))
    print(f"[info] Preserved (baseline): {preserved} atoms — {sorted(PRESERVE_BASELINE)}")
    print(f"[info] Overridden (Khoury):  {overrides} atoms")
    print(f"[info] Σq pre-residual: {new_sigma_q:+.6f}")

    if abs(new_sigma_q) > 1e-6:
        # Distribute residual across overridden atoms only (preserve exact baseline)
        residual = -new_sigma_q
        per_atom_adj = residual / overrides
        for atom in residue.findall("Atom"):
            if atom.get("name") not in PRESERVE_BASELINE:
                q = float(atom.get("charge")) + per_atom_adj
                atom.set("charge", f"{q:.6f}")
        new_sigma_q = sum(float(a.get("charge")) for a in residue.findall("Atom"))
        print(f"[info] Residual {residual:+.6f} distributed across {overrides} overridden atoms ({per_atom_adj:+.6f} each)")

    print(f"[info] Final Σq: {new_sigma_q:+.2e}")

    if missing:
        print(f"[WARN] Atoms in template not covered by Khoury OMW: {missing}")

    OUT_XML.parent.mkdir(parents=True, exist_ok=True)
    tree.write(OUT_XML, encoding="utf-8", xml_declaration=False)
    print(f"[OK] Wrote: {OUT_XML}")

    final_charges = {a.get("name"): float(a.get("charge")) for a in residue.findall("Atom")}
    delta_vs_template = {n: final_charges[n] - template_charges[n] for n in final_charges}
    delta_vs_full_resp = {
        n: final_charges[n] - KHOURY_OMW_CHARGES.get(n, template_charges[n])
        for n in final_charges
    }
    rmsd_vs_template = (sum(d * d for d in delta_vs_template.values()) / len(delta_vs_template)) ** 0.5
    max_abs_vs_template = max(abs(d) for d in delta_vs_template.values())

    manifest = {
        "schema": "mtr_hybrid_manifest/0.1",
        "generated_by": "scripts/build_khoury_mtr_hybrid_xml.py",
        "rationale": (
            "Sidechain-only Khoury RESP-A2 override; backbone N/H + indole NE1 "
            "retained from amber14SB-patched baseline to prevent the MD-propagation-"
            "level junction limitation (NE1 sign flip → HD1 +/+ Coulomb repulsion → "
            "ring planarity destabilization → NaN cascade) observed in 2/8 seeds of "
            "the v0.8 full-RESP pilot."
        ),
        "preserved_baseline_atoms": sorted(PRESERVE_BASELINE),
        "overridden_atoms": sorted(overridden_atom_names),
        "missing_from_khoury": missing,
        "source_khoury": {
            "doi": "10.1021/sb400168u",
            "block": "ffncaa.in OMW (1-methyl-tryptophan)",
            "block_line": "778-830",
            "convention": "internal-residue, ff03-compatible, HF/6-31G* ESP, RESP-A2 fit",
        },
        "atom_name_mapping_khoury_to_updd": {
            "CZ1": "CM",
            "HZ11": "HM1",
            "HZ12": "HM2",
            "HZ13": "HM3",
        },
        "template_baseline": {
            "path": str(TEMPLATE_XML),
            "convention": "amber14SB-patched (UPDD_NCAA_AMBER14_PATCH=1 baked into XML)",
            "sigma_q": round(template_sigma_q, 6),
        },
        "n_atoms": len(final_charges),
        "n_preserved": preserved,
        "n_overridden": overrides,
        "final_sigma_q": round(new_sigma_q, 6),
        "vs_amber14sb_patched_baseline": {
            "rmsd_e": round(rmsd_vs_template, 4),
            "max_abs_delta_e": round(max_abs_vs_template, 4),
            "delta_per_atom": {k: round(v, 4) for k, v in delta_vs_template.items()},
        },
        "vs_full_khoury_resp": {
            "delta_per_atom": {k: round(v, 4) for k, v in delta_vs_full_resp.items()},
            "comment": "Non-zero entries identify the preserved baseline atoms (N, H, NE1) + residual redistribution drift.",
        },
        "final_charges_per_atom": {k: round(v, 6) for k, v in final_charges.items()},
    }
    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with OUT_MANIFEST.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"[OK] Wrote: {OUT_MANIFEST}")
    print()
    print(f"=== Δq vs amber14SB-patched baseline ===")
    print(f"RMSD     = {rmsd_vs_template:.4f} e")
    print(f"max |Δq| = {max_abs_vs_template:.4f} e")
    print(f"Σq (final) = {new_sigma_q:+.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
