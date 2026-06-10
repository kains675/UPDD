#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Track B production blocker Q1 — Σq integer rescale on MTR/NMTR/CMTR.

Production blocker Q1 (a):
- Rescale scope: MTR sidechain heavy + H atoms ONLY
- Frozen: backbone N/H/CA/HA/C/O (amber14SB) + NE1 (B-C1 anchor, −0.3418)
  + OXT/HXT (C-term oxygen/hydroxyl, amber14SB Trp C-term values)
- Distribution: uniform per-atom Δq = −Σq_residue / N_rescale_atoms
- Acceptance: post-rescale per-residue Σq ≤ |5e-4| e (integer-parity)

Three residue variants are processed independently (residue-scope summation):
  MTR  (internal,  27 atoms) — pre-existing Σq ≈ −5e-6 e (already integer)
  NMTR (N-term,    27 atoms) — pre-existing Σq ≈  0     e (already integer)
  CMTR (C-term,    29 atoms) — pre-existing Σq ≈ −0.176 e (needs rescale)

Per-atom Δ for CMTR ≈ +0.0088 e, well inside Khoury (2014) ±0.02 e per-atom
parameterization tolerance.

NOTE on the −0.176 e provenance: the previous launcher/verifier used a flat
``Residues.iter('Atom')`` walk with a name-keyed dict that overwrote (CMTR is
last in iteration → its values won, OXT/HXT extra atoms contributed). The
true per-residue scope is what production FF code (OpenMM ``createSystem``)
sees, so per-residue is the correct scope to enforce integer-parity.
"""

from __future__ import annotations

import argparse
import os
import sys
import xml.etree.ElementTree as ET
from typing import Dict, List, Tuple

# Per-residue acceptance: |Σq| ≤ 5e-4 e (integer-parity)
ACCEPT_TOL = 5e-4

# Frozen atoms by name (not rescale-eligible):
#   backbone amber14SB:  N, H, CA, HA, C, O
#   NE1 (B-C1 anchor):   NE1
#   C-term cap atoms:    OXT, HXT (amber14SB Trp C-term values)
FROZEN_ATOM_NAMES = {"N", "H", "CA", "HA", "C", "O", "NE1", "OXT", "HXT"}


def _read_residues(xml_path: str) -> ET.ElementTree:
    parser = ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))
    tree = ET.parse(xml_path, parser=parser)
    return tree


def _residue_charge_table(res_el: ET.Element) -> List[Tuple[str, str, float]]:
    """Return (name, type, charge) tuples for every <Atom> in this residue."""
    out: List[Tuple[str, str, float]] = []
    for a in res_el.findall("Atom"):
        out.append((a.get("name"), a.get("type"), float(a.get("charge"))))
    return out


def _compute_residue_sigma_q(res_el: ET.Element) -> Dict[str, float]:
    """Per-residue Σq summary with frozen/rescale partition."""
    atoms = _residue_charge_table(res_el)
    sigma_all = sum(q for _, _, q in atoms)
    frozen = [(n, t, q) for n, t, q in atoms if n in FROZEN_ATOM_NAMES]
    rescale = [(n, t, q) for n, t, q in atoms if n not in FROZEN_ATOM_NAMES]
    return {
        "n_atoms": len(atoms),
        "sigma_q_all": sigma_all,
        "n_frozen": len(frozen),
        "sigma_q_frozen": sum(q for _, _, q in frozen),
        "n_rescale": len(rescale),
        "sigma_q_rescale_current": sum(q for _, _, q in rescale),
    }


def rescale_residue(res_el: ET.Element, target_sigma: float = 0.0) -> Dict[str, object]:
    """Apply uniform per-atom Δq on rescale-eligible atoms to drive Σq → 0.

    Returns an audit log with per-atom Δq deltas, pre/post Σq, and rescale scope.
    """
    rname = res_el.get("name")
    pre = _compute_residue_sigma_q(res_el)
    delta_total = target_sigma - pre["sigma_q_all"]

    rescale_atoms = [a for a in res_el.findall("Atom")
                     if a.get("name") not in FROZEN_ATOM_NAMES]
    n_rescale = len(rescale_atoms)
    if n_rescale == 0:
        return {
            "residue": rname,
            "pre": pre,
            "delta_total_e": delta_total,
            "per_atom_delta_e": 0.0,
            "post_sigma_q": pre["sigma_q_all"],
            "atom_changes": [],
            "skipped": True,
            "reason": "no_rescale_eligible_atoms",
        }

    per_atom_delta = delta_total / n_rescale

    atom_changes: List[Dict[str, object]] = []
    for a in rescale_atoms:
        name = a.get("name")
        q_old = float(a.get("charge"))
        q_new = q_old + per_atom_delta
        a.set("charge", repr(q_new))
        atom_changes.append({
            "name": name,
            "type": a.get("type"),
            "q_old": q_old,
            "q_new": q_new,
            "delta": per_atom_delta,
        })

    post = _compute_residue_sigma_q(res_el)

    return {
        "residue": rname,
        "pre": pre,
        "target_sigma_q": target_sigma,
        "delta_total_e": delta_total,
        "per_atom_delta_e": per_atom_delta,
        "n_rescale_atoms": n_rescale,
        "post_sigma_q": post["sigma_q_all"],
        "post_within_tol": abs(post["sigma_q_all"]) <= ACCEPT_TOL,
        "atom_changes": atom_changes,
        "skipped": False,
    }


def main():
    p = argparse.ArgumentParser(
        description="Rescale MTR/NMTR/CMTR per-residue Σq to integer (production blocker Q1 a)"
    )
    p.add_argument("--xml", default="params/MTR_gaff2_hybrid.xml",
                   help="MTR hybrid XML to patch in-place")
    p.add_argument("--dry-run", action="store_true",
                   help="Report rescale plan + audit but do not write XML")
    p.add_argument("--audit-out",
                   default="params/_archive/mtr_rescale_audit_20260530.json",
                   help="Write JSON audit log here")
    args = p.parse_args()

    xml_path = os.path.abspath(args.xml)
    if not os.path.isfile(xml_path):
        print(f"ERROR: {xml_path} not found", file=sys.stderr)
        sys.exit(2)

    tree = _read_residues(xml_path)
    root = tree.getroot()
    residues_el = root.find(".//Residues")
    if residues_el is None:
        print(f"ERROR: no <Residues> block in {xml_path}", file=sys.stderr)
        sys.exit(2)

    print(f"Track B Q1 rescale — input: {xml_path}")
    print(f"Acceptance tolerance: |Σq| ≤ {ACCEPT_TOL:.0e} e per residue")
    print(f"Frozen atom names: {sorted(FROZEN_ATOM_NAMES)}")
    print()

    audit: List[Dict[str, object]] = []
    for res_el in residues_el.findall("Residue"):
        rname = res_el.get("name")
        pre = _compute_residue_sigma_q(res_el)
        print(f"--- {rname} ---")
        print(f"  N_atoms = {pre['n_atoms']}  N_frozen = {pre['n_frozen']}  "
              f"N_rescale = {pre['n_rescale']}")
        print(f"  Σq_all     = {pre['sigma_q_all']:+.6f} e")
        print(f"  Σq_frozen  = {pre['sigma_q_frozen']:+.6f} e")
        print(f"  Σq_rescale = {pre['sigma_q_rescale_current']:+.6f} e")
        if abs(pre["sigma_q_all"]) <= ACCEPT_TOL:
            print(f"  STATUS: already within tol ({pre['sigma_q_all']:+.2e}). "
                  f"Skip (delta would be 0).")
            audit.append({
                "residue": rname, "pre": pre,
                "skipped": True, "reason": "within_tolerance",
                "post_sigma_q": pre["sigma_q_all"],
                "post_within_tol": True,
            })
            print()
            continue

        log = rescale_residue(res_el, target_sigma=0.0)
        print(f"  RESCALE: Δ_total = {log['delta_total_e']:+.6f} e  "
              f"per-atom Δ = {log['per_atom_delta_e']:+.8f} e "
              f"(|Δ| = {abs(log['per_atom_delta_e']):.2e}, "
              f"Khoury tol ±2e-2)")
        print(f"  POST   Σq = {log['post_sigma_q']:+.2e} e   "
              f"(within tol: {log['post_within_tol']})")
        if log["per_atom_delta_e"] and abs(log["per_atom_delta_e"]) > 2e-2:
            print(f"  WARNING: per-atom |Δ| > Khoury 2e-2 e tolerance")
        audit.append(log)
        print()

    # Verify whole-XML integer parity after patch (per-residue sum of pre-patch
    # diff vs post-patch — should be exactly the planned Δ_total).
    print("=== Post-patch per-residue Σq (final check) ===")
    all_pass = True
    for res_el in residues_el.findall("Residue"):
        post = _compute_residue_sigma_q(res_el)
        rname = res_el.get("name")
        ok = abs(post["sigma_q_all"]) <= ACCEPT_TOL
        all_pass &= ok
        mark = "OK" if ok else "FAIL"
        print(f"  {rname:5s}: Σq = {post['sigma_q_all']:+.6e} e  [{mark}]")
    print()
    print(f"OVERALL: {'PASS' if all_pass else 'FAIL'} "
          f"(all residue |Σq| ≤ {ACCEPT_TOL:.0e})")

    # Persist audit
    import json
    os.makedirs(os.path.dirname(os.path.abspath(args.audit_out)), exist_ok=True)
    with open(args.audit_out, "w") as fh:
        json.dump({
            "xml_path": xml_path,
            "tolerance_e": ACCEPT_TOL,
            "frozen_atom_names": sorted(FROZEN_ATOM_NAMES),
            "scope": "per_residue",
            "method": "uniform_per_atom_delta_on_rescale_eligible",
            "audit": audit,
            "all_pass": all_pass,
            "regime": "ranking_only",
        }, fh, indent=2)
    print(f"Audit written: {args.audit_out}")

    if args.dry_run:
        print("\nDRY RUN — XML NOT modified.")
        return 0 if all_pass else 1

    # Persist patched XML
    tree.write(xml_path, encoding="utf-8", xml_declaration=False)
    print(f"\nXML patched in-place: {xml_path}")
    return 0 if all_pass else 1


if __name__ == "__main__":
    sys.exit(main())
