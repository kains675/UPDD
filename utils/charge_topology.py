#!/usr/bin/env python
"""utils/charge_topology.py — Stage 1 (R-15/R-16) chemistry-true charge helpers.

This module contains the pure-stdlib implementation of the chemistry-true
charge computation used by:

    * ``run_qmmm.py`` R-15 / R-16 guards (build_qm_mol path)
    * ``audit_charge_consistency.py`` historical JSON audit (Stage 1)
    * ``tests/test_charge_consistency.py`` unit tests

Keeping the implementation in its own module avoids the heavy pyscf/openmm
imports of ``run_qmmm.py`` when we only need to parse a PDB and classify
residue protonation states (SciVal verdict v2 §5).

References:
    - SciVal verdict_target_card_charge_rationale_20260419_v2.md §1.1, §5
    - Senn & Thiel 2009 Angew 48, 1198 (H-link charge = 0 premise)
    - Maier 2015 JCTC 11, 3696 (ff14SB HIS tautomer parameters)

Constitution: R-15 (magnitude) + R-16 (binder SSOT) — see CLAUDE.md Rev. 2.
"""

from collections import OrderedDict
from typing import Any, Dict, Iterable, List, Optional, Tuple


# ==========================================
# Residue-name encoded protonation variants (AMBER ff14SB)
# ==========================================
_ASP_PROTONATED_NAMES = ("ASH",)
_GLU_PROTONATED_NAMES = ("GLH",)
_LYS_DEPROTONATED_NAMES = ("LYN",)
_HIS_POSITIVE_NAMES = ("HIP",)
_HIS_NEUTRAL_NAMES = ("HID", "HIE")


def _residue_group_atom_names(residue_atoms: Iterable[Dict[str, Any]]):
    """Return the set of atom-name strings (stripped) for a residue record."""
    return {str(a.get("name", "")).strip() for a in residue_atoms}


def classify_residue_charge(resname: str, atom_names: Iterable[str]) -> int:
    """Canonical net charge at pH 7.4 for a single residue record.

    Pure function: 3-letter residue name plus the set of atom names present.
    No bond/connectivity heuristics beyond AMBER/OpenMM naming conventions.

    Args:
        resname: 3-letter residue name (case-insensitive).
        atom_names: iterable of atom-name strings (whitespace stripped).

    Returns:
        int: canonical residue net charge (−1, 0, or +1) at pH 7.4.

    Notes:
        - TYR pKa ~10 → 0 at pH 7.4 (PROPKA integration TODO).
        - CYS pKa ~8.3 → 0 at pH 7.4 (disulfide CYX / cysteinate CYM not
          yet detected).
        - HIS tautomer inference: HD1+HE2 both present → HIP (+1); HD1 only
          → HID (0); HE2 only → HIE (0).
    """
    rn = (resname or "").upper().strip()
    names = set(atom_names)

    if rn in _ASP_PROTONATED_NAMES:
        return 0
    if rn in _GLU_PROTONATED_NAMES:
        return 0
    if rn in _LYS_DEPROTONATED_NAMES:
        return 0
    if rn in _HIS_POSITIVE_NAMES:
        return +1
    if rn in _HIS_NEUTRAL_NAMES:
        return 0

    if rn == "ARG":
        hh_set = {"HH11", "HH12", "HH21", "HH22"}
        return +1 if hh_set.issubset(names) else 0

    if rn == "LYS":
        hz_set = {"HZ1", "HZ2", "HZ3"}
        return +1 if hz_set.issubset(names) else 0

    if rn == "ASP":
        # Residue name "ASP" (not "ASH") → deprotonated carboxylate at pH 7.4.
        return -1

    if rn == "GLU":
        return -1

    if rn == "HIS":
        has_hd1 = "HD1" in names
        has_he2 = "HE2" in names
        if has_hd1 and has_he2:
            return +1
        return 0

    if rn == "TYR":
        return 0

    if rn == "CYS":
        return 0

    return 0


def residues_by_chain_resnum(
    atoms: Iterable[Dict[str, Any]],
    chain_filter: Optional[str] = None,
) -> "OrderedDict[Tuple[str, int], Tuple[str, List[Dict[str, Any]]]]":
    """Group atoms by (chain, resnum) preserving insertion order."""
    out: "OrderedDict[Tuple[str, int], Tuple[str, List[Dict[str, Any]]]]" = OrderedDict()
    for a in atoms:
        if chain_filter is not None and a.get("chain") != chain_filter:
            continue
        key = (a.get("chain"), int(a.get("resnum", 0)))
        if key not in out:
            out[key] = (a.get("resname", ""), [])
        out[key][1].append(a)
    return out


def detect_cyclic_peptide(binder_group) -> Tuple[bool, Dict[str, Any]]:
    """Heuristic cyclic-peptide detection.

    Linear peptide:
        - First residue carries 2-3 hydrogens on backbone N (NH3+ terminus).
        - Last residue carries an OXT atom (free C-terminal carboxylate).
    Head-to-tail cyclic peptide:
        - First residue has exactly 1 backbone N-H (peptide-bonded).
        - Last residue has no OXT (C=O bonded to N of residue 1).

    Args:
        binder_group: OrderedDict from :func:`residues_by_chain_resnum`.

    Returns:
        (is_cyclic, diag) with diag containing the measurements used.
    """
    if not binder_group:
        return False, {"reason": "empty_binder_group"}

    resnums = sorted({k[1] for k in binder_group.keys()})
    if not resnums:
        return False, {"reason": "no_resnums"}

    first_resnum = resnums[0]
    last_resnum = resnums[-1]
    first_key = next(k for k in binder_group.keys() if k[1] == first_resnum)
    last_key = next(k for k in binder_group.keys() if k[1] == last_resnum)

    first_rname, first_atoms = binder_group[first_key]
    last_rname, last_atoms = binder_group[last_key]

    first_names = _residue_group_atom_names(first_atoms)
    last_names = _residue_group_atom_names(last_atoms)

    nh_candidates = {"H", "H1", "H2", "H3"}
    n_nh_first = len(first_names & nh_candidates)
    has_oxt_last = "OXT" in last_names

    is_cyclic = (n_nh_first <= 1) and (not has_oxt_last)
    return is_cyclic, {
        "first_residue": {
            "chain": first_key[0], "resnum": first_resnum,
            "resname": first_rname, "n_backbone_NH": n_nh_first,
        },
        "last_residue": {
            "chain": last_key[0], "resnum": last_resnum,
            "resname": last_rname, "has_OXT": bool(has_oxt_last),
        },
    }


def compute_binder_chem_charge(
    atoms: Iterable[Dict[str, Any]],
    binder_chain: str,
    ncaa_charge_map: Optional[Dict[str, int]] = None,
) -> Tuple[int, Dict[str, Any]]:
    """Compute chemistry-true binder net charge from a parsed atom list.

    R-16: binder net charge is derived at runtime from the snapshot PDB;
    ``target_card.binder_net_charge`` is only a hint.

    Args:
        atoms: iterable of atom dicts (chain, resnum, resname, name, element).
        binder_chain: chain letter of the binder (e.g. ``"B"``).
        ncaa_charge_map: optional {residue_name: net_charge}. Defaults
            include NME / NMA / ACE = 0 (neutral caps / N-methyl alanine).

    Returns:
        (binder_chem_charge, diagnostic).
    """
    ncaa_charge_map = dict(ncaa_charge_map or {})
    # R-16: seed defaults from the ncAA registry (single source of truth) so
    # charged ncAAs (ORN/DAB/HAR/NMK/NMR/DAR/TPO/PTR/…) are recognized
    # without requiring callers to repeat the table.
    try:
        from ncaa_registry import get_ncaa_charge_map  # type: ignore
        for _rn, _q in get_ncaa_charge_map().items():
            ncaa_charge_map.setdefault(_rn, _q)
    except ImportError:
        pass
    ncaa_charge_map.setdefault("NME", 0)
    ncaa_charge_map.setdefault("NMA", 0)
    ncaa_charge_map.setdefault("ACE", 0)

    binder_group = residues_by_chain_resnum(atoms, chain_filter=binder_chain)
    is_cyclic, cyclic_diag = detect_cyclic_peptide(binder_group)

    per_residue = []
    total = 0
    for (chain, resnum), (resname, res_atoms) in binder_group.items():
        names = _residue_group_atom_names(res_atoms)
        rn = (resname or "").upper().strip()
        if rn in ncaa_charge_map:
            q = int(ncaa_charge_map[rn])
        else:
            q = classify_residue_charge(rn, names)
        total += q
        per_residue.append({
            "chain": chain, "resnum": resnum, "resname": rn, "charge": q,
        })

    termini_delta = 0
    if not is_cyclic and binder_group:
        # Neutral linear peptide: +1 NH3+ and −1 COO− cancel.
        termini_delta = 0
    total += termini_delta

    return int(total), {
        "cyclic": bool(is_cyclic),
        "cyclic_detection": cyclic_diag,
        "per_residue": per_residue,
        "termini_delta": termini_delta,
        "ncaa_charge_map": ncaa_charge_map,
    }


def compute_qm_net_charge_topology(
    atoms: Iterable[Dict[str, Any]],
    binder_chain: str,
    target_chain: str,
    target_contact_residues: Iterable[int],
    whole_residue_exceptions: Optional[Iterable[int]] = None,
    pH: float = 7.4,
    ncaa_charge_map: Optional[Dict[str, int]] = None,
) -> Tuple[int, int, int, Dict[str, Any]]:
    """Compute chemistry-true QM net charge for the topology-mode QM region.

    Returns a 4-tuple ``(binder_charge, target_charge, total_charge, diag)``.

    The function is side-effect free; all logging is left to the caller.
    """
    contact_set = {int(r) for r in (target_contact_residues or [])}
    whole_set = {int(r) for r in (whole_residue_exceptions or [])}

    binder_charge, binder_diag = compute_binder_chem_charge(
        atoms, binder_chain, ncaa_charge_map=ncaa_charge_map,
    )

    target_group = residues_by_chain_resnum(atoms, chain_filter=target_chain)
    target_charge = 0
    target_per_residue = []
    for (chain, resnum), (resname, res_atoms) in target_group.items():
        if int(resnum) not in contact_set:
            continue
        names = _residue_group_atom_names(res_atoms)
        rn = (resname or "").upper().strip()
        q = classify_residue_charge(rn, names)
        target_charge += q
        target_per_residue.append({
            "chain": chain,
            "resnum": int(resnum),
            "resname": rn,
            "charge": q,
            "whole_residue": int(resnum) in whole_set,
        })

    total = int(binder_charge) + int(target_charge)
    diag = {
        "pH": float(pH),
        "binder_chain": binder_chain,
        "target_chain": target_chain,
        "target_contact_residues": sorted(contact_set),
        "whole_residue_exceptions": sorted(whole_set),
        "binder_charge": int(binder_charge),
        "target_charge": int(target_charge),
        "total_charge": int(total),
        "binder_diag": binder_diag,
        "target_per_residue": target_per_residue,
    }
    return int(binder_charge), int(target_charge), int(total), diag


# ==========================================
# Minimal PDB parser — no numpy, no pyscf
# ==========================================
# Extracted so the audit module and unit tests can load snapshots without
# importing the heavy ``run_qmmm`` stack.
_WATER_RESNAMES = {"HOH", "WAT", "TIP3", "TIP", "SOL"}


def parse_pdb_atoms_lite(pdb_path: str) -> List[Dict[str, Any]]:
    """Stream-parse a PDB file into a list of atom dicts (stdlib only).

    Fields returned per atom:
        record, serial, name, resname, chain, resnum, x, y, z, element.

    Water residues are filtered out. If the element column (76-78) is
    blank, we fall back to the first alphabetic character of the atom name.
    """
    atoms: List[Dict[str, Any]] = []
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            try:
                resname = line[17:20].strip()
                if resname in _WATER_RESNAMES:
                    continue
                name = line[12:16].strip()
                element = line[76:78].strip() if len(line) > 76 else ""
                if not element:
                    stripped = name.lstrip("0123456789")
                    element = stripped[:1].upper() if stripped else "C"
                atoms.append({
                    "record": rec,
                    "serial": int(line[6:11]),
                    "name": name,
                    "resname": resname,
                    "chain": line[21].strip(),
                    "resnum": int(line[22:26]),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "element": element,
                })
            except (ValueError, IndexError):
                continue
    return atoms
