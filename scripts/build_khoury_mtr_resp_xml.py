#!/usr/bin/env python
"""
Build params/MTR_gaff2_resp.xml from the existing MTR_gaff2.xml template by
overriding partial charges with Khoury 2014 OMW (1-methyl-tryptophan) RESP-A2
values. Khoury OMW block (ffncaa.in L778-830) provides 27 atoms, Σq = 0
under the internal-residue convention (ff03-compatible, HF/6-31G* ESP fit).

UPDD ↔ Khoury atom-name mapping:
  - Khoury CZ1  → UPDD CM
  - Khoury HZ11 → UPDD HM1
  - Khoury HZ12 → UPDD HM2
  - Khoury HZ13 → UPDD HM3
  - All other 23 atom names match.

Inputs:
  - Template: outputs/2QKI_Cp4_calib_s7/params/MTR_gaff2.xml

Outputs:
  - params/MTR_gaff2_resp.xml
  - params/MTR_gaff2_resp_manifest.json (provenance + Δq vs baseline)
"""
from __future__ import annotations

import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
TEMPLATE_XML = PROJ / "outputs/2QKI_Cp4_calib_s7/params/MTR_gaff2.xml"
OUT_XML = PROJ / "params/MTR_gaff2_resp.xml"
OUT_MANIFEST = PROJ / "params/MTR_gaff2_resp_manifest.json"

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

KHOURY_SIGMA_Q = sum(KHOURY_OMW_CHARGES.values())


def main() -> int:
    if not TEMPLATE_XML.exists():
        print(f"[ERROR] Template XML not found: {TEMPLATE_XML}", file=sys.stderr)
        return 1

    print(f"[info] Σq (Khoury OMW) = {KHOURY_SIGMA_Q:+.6f}")
    if abs(KHOURY_SIGMA_Q) > 1e-5:
        print(f"[ERROR] Khoury OMW Σq not neutral: {KHOURY_SIGMA_Q}", file=sys.stderr)
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

    overrides = 0
    missing = []
    template_charges = {}
    for atom in residue.findall("Atom"):
        name = atom.get("name")
        template_charges[name] = float(atom.get("charge", "0.0"))
        if name in KHOURY_OMW_CHARGES:
            atom.set("charge", f"{KHOURY_OMW_CHARGES[name]:.6f}")
            overrides += 1
        else:
            missing.append(name)

    new_sigma_q = sum(float(a.get("charge")) for a in residue.findall("Atom"))

    if abs(new_sigma_q) > 1e-6:
        residual = -new_sigma_q
        atoms_list = list(residue.findall("Atom"))
        per_atom_adj = residual / len(atoms_list)
        for atom in atoms_list:
            q = float(atom.get("charge")) + per_atom_adj
            atom.set("charge", f"{q:.6f}")
        new_sigma_q = sum(float(a.get("charge")) for a in residue.findall("Atom"))
        print(f"[info] Residual {residual:+.2e} distributed across {len(atoms_list)} atoms")

    print(f"[info] Overrides applied: {overrides}")
    print(f"[info] New Σq (post-adjust): {new_sigma_q:+.2e}")
    if missing:
        print(f"[WARN] Atoms in template not covered by Khoury OMW: {missing}")

    OUT_XML.parent.mkdir(parents=True, exist_ok=True)
    tree.write(OUT_XML, encoding="utf-8", xml_declaration=False)
    print(f"[OK] Wrote: {OUT_XML}")

    final_charges = {a.get("name"): float(a.get("charge")) for a in residue.findall("Atom")}
    delta = {n: final_charges[n] - template_charges[n] for n in final_charges}
    rmsd = (sum(d * d for d in delta.values()) / len(delta)) ** 0.5
    max_abs = max(abs(d) for d in delta.values())

    manifest = {
        "schema": "mtr_resp_manifest/0.1",
        "generated_by": "scripts/build_khoury_mtr_resp_xml.py",
        "source": {
            "khoury_2014_doi": "10.1021/sb400168u",
            "block": "ffncaa.in OMW (= 1-methyl-tryptophan)",
            "block_line": "778-830",
            "convention": "internal-residue, ff03-compatible, HF/6-31G* ESP, RESP-A2 fit",
        },
        "atom_name_mapping": {
            "CZ1 (Khoury)": "CM (UPDD)",
            "HZ11 (Khoury)": "HM1 (UPDD)",
            "HZ12 (Khoury)": "HM2 (UPDD)",
            "HZ13 (Khoury)": "HM3 (UPDD)",
        },
        "khoury_sigma_q": round(KHOURY_SIGMA_Q, 6),
        "n_atoms": len(final_charges),
        "overrides_applied": overrides,
        "missing_from_khoury": missing,
        "final_sigma_q": round(new_sigma_q, 6),
        "vs_template_baseline": {
            "template_path": str(TEMPLATE_XML),
            "rmsd_delta_q_e": round(rmsd, 4),
            "max_abs_delta_q_e": round(max_abs, 4),
            "delta_q_per_atom": {k: round(v, 4) for k, v in delta.items()},
        },
        "final_charges_per_atom": {k: round(v, 6) for k, v in final_charges.items()},
    }
    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with OUT_MANIFEST.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"[OK] Wrote: {OUT_MANIFEST}")
    print()
    print(f"=== Δq vs baseline ===")
    print(f"RMSD     = {rmsd:.4f} e")
    print(f"max |Δq| = {max_abs:.4f} e")
    print(f"Σq (final) = {new_sigma_q:+.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
