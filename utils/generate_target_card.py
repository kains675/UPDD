#!/usr/bin/env python3
"""
generate_target_card.py — UPDD target card generator (v0.6.1)

MD snapshot PDB에서 binder-target contact residues를 자동 추출하여
target_card.json을 생성합니다. 어떤 타겟에도 범용으로 사용 가능합니다.

사용법:
    python utils/generate_target_card.py --pdb <snap.pdb> --target-chain A \
        --binder-chain B --cutoff 4.5 --target-id 6WGN --output target_cards/6WGN.json
"""

import argparse
import json
import os
import sys
from collections import defaultdict
from datetime import date

import numpy as np


_BACKBONE_NAMES = frozenset({'N', 'CA', 'C', 'O', 'H', 'HA',
                              'OXT', 'H1', 'H2', 'H3', 'HA2', 'HA3'})
_SOLVENT_RESNAMES = frozenset({'HOH', 'WAT', 'TIP', 'TIP3', 'TIP4',
                               'SOL', 'NA', 'CL', 'K', 'MG', 'ZN'})
_STANDARD_AA = frozenset({
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
    'TYR', 'VAL', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'CYX', 'CYM',
    'ASH', 'GLH', 'LYN',  # protonation variants
    'NME', 'ACE',          # capping groups
})


# ==========================================
# [v0.6.5 Phase 4] target_iso 전하 자동 계산 (pH 7.4 기준)
# ==========================================
# Reference protonation state at physiological pH 7.4:
#   - ASP/GLU: deprotonated carboxylate (-1)
#   - ARG/LYS: protonated amine/guanidinium (+1)
#   - HIS: pKa ~6.0 → mostly neutral; HIP variant explicitly +1
#   - CYS sidechain (sulfhydryl): neutral at pH 7.4 (pKa ~8.3); CYM (-1) explicit
#   - All others: 0
# Note: protonation variants (ASH, GLH, LYN) are explicitly listed.
# This is an *approximation* for target_iso_net_charge — for high-precision work,
# perform PROPKA / H++ analysis on the actual MD topology.
_RESIDUE_CHARGE_AT_PH7 = {
    "ASP": -1, "GLU": -1, "ASH": 0, "GLH": 0,        # Asp/Glu (protonated → 0)
    "ARG": +1, "LYS": +1, "LYN": 0,                  # Arg/Lys (deprotonated Lys → 0)
    "HIS": 0, "HIE": 0, "HID": 0, "HIP": +1,         # Histidine variants
    "HSD": 0, "HSE": 0, "HSP": +1,                   # CHARMM histidine variants
    "CYS": 0, "CYX": 0, "CYM": -1,                   # Cys (CYM = thiolate)
    "TYR": 0,                                          # Tyr neutral at pH 7.4 (pKa ~10)
}


def compute_target_iso_charge(contact_residues_with_names):
    """Compute net charge of target sidechain subset at pH 7.4.

    Used by Phase 4 (qm_int_kcal_frozen) to seed E_target_iso SCF with a
    physiologically reasonable charge. The resulting value is *approximate* —
    titratable residues with pKa shifts in the binding pocket may deviate.

    Args:
        contact_residues_with_names: list of dicts with keys 'resnum', 'resname'.

    Returns:
        Tuple[int, str]: (net_charge, human-readable rationale).
    """
    total = 0
    parts = []
    for r in contact_residues_with_names:
        rn = (r.get("resname") or "").upper()
        q = _RESIDUE_CHARGE_AT_PH7.get(rn, 0)
        if q != 0:
            parts.append("{}{}({:+d})".format(rn, r.get("resnum", "?"), q))
            total += q
    rationale = " + ".join(parts) if parts else "all neutral"
    return total, "pH7.4: {} = {:+d}".format(rationale, total)


# ==========================================
# Stage 2 (2026-04-19) — schema 0.6.4 auto-rationale generator
# ==========================================
def generate_target_iso_rationale(
    contact_residues_with_names,
    md_to_crystal_offset=0,
    schema_version="0.6.4",
    today_iso=None,
):
    """Generate an R-15 compliant target_iso_charge_rationale string.

    The string includes both MD and crystal numbering for every contributing
    residue (per Stage 2 numbering_convention spec) and a breakdown of the
    ionization contribution (ASP/GLU/ARG/LYS/HIP counts).

    Args:
        contact_residues_with_names: list of {resnum, resname} dicts in MD
            numbering (as produced by ``extract_contacts``).
        md_to_crystal_offset: integer offset (crystal_resnum = md_resnum +
            offset). 0 means MD and crystal share a coordinate system.
        schema_version: emitted verbatim into the rationale for traceability.
        today_iso: ISO date string (YYYY-MM-DD). Defaults to today.

    Returns:
        Tuple[int, str]: (net_charge, rationale_string).
    """
    if today_iso is None:
        today_iso = str(date.today())

    total = 0
    ionized_details = []   # list of (resname, md, crystal, charge)
    neutral_details = []   # list of (md, crystal, resname)
    tautomer_notes = []    # e.g. HIS92/95 = HID (neutral)

    for r in contact_residues_with_names:
        md_num = r.get("resnum")
        rn = (r.get("resname") or "").upper()
        crystal_num = (md_num + md_to_crystal_offset
                       if md_num is not None else None)
        q = _RESIDUE_CHARGE_AT_PH7.get(rn, 0)
        if q != 0:
            ionized_details.append((rn, md_num, crystal_num, q))
            total += q
        else:
            neutral_details.append((md_num, crystal_num, rn))
        # Mark explicit HIS tautomer variants for human auditability.
        if rn in ("HID", "HIE"):
            tautomer_notes.append(
                "HIS{}/{} = {} (neutral)".format(md_num, crystal_num, rn)
            )
        elif rn == "HIP":
            tautomer_notes.append(
                "HIS{}/{} = HIP (+1 explicit)".format(md_num, crystal_num)
            )
        elif rn == "HIS":
            # Ambiguous HIS in card — flagged but treated as neutral default.
            tautomer_notes.append(
                "HIS{}/{} = HIS (tautomer resolved at runtime from atom names)"
                .format(md_num, crystal_num)
            )

    # Breakdown counts by residue family.
    by_family = {"GLU": 0, "ASP": 0, "ARG": 0, "LYS": 0, "HIP": 0}
    for rn, _, _, q in ionized_details:
        if rn in by_family:
            by_family[rn] += 1
    breakdown_parts = []
    for fam, sign, label in (
        ("GLU", -1, "GLU"), ("ASP", -1, "ASP"),
        ("ARG", +1, "ARG"), ("LYS", +1, "LYS"),
        ("HIP", +1, "HIP"),
    ):
        if by_family[fam]:
            breakdown_parts.append(
                "{}x{}({:+d})".format(by_family[fam], label, sign)
            )
    if not breakdown_parts:
        breakdown_parts.append("all neutral")
    breakdown_str = " + ".join(breakdown_parts)

    # Ionized residue list (MD/crystal pairs).
    ion_list = ", ".join(
        "{}{}/{}".format(rn, md, crys)
        for rn, md, crys, _ in ionized_details
    ) or "none"

    # Non-ionizable contact list (MD/crystal pairs, resname in parens).
    neutral_list_parts = []
    for md, crys, rn in neutral_details:
        # Skip HIS variants: they already appear in tautomer_notes.
        if rn in ("HIS", "HID", "HIE", "HIP"):
            continue
        neutral_list_parts.append("{}/{}".format(md, crys))
    neutral_list = (", ".join(neutral_list_parts)
                    if neutral_list_parts else "none")

    tautomer_str = "; ".join(tautomer_notes) if tautomer_notes else "none"

    rationale = (
        "Auto-generated {today} (schema {schema}): "
        "target QM region net charge = {total:+d} = {breakdown}. "
        "Contributing residues (MD/crystal): {ionized}. "
        "Tautomers: {tautomers}. "
        "Non-ionizable contacts (MD/crystal): {neutrals}."
    ).format(
        today=today_iso,
        schema=schema_version,
        total=total,
        breakdown=breakdown_str,
        ionized=ion_list,
        tautomers=tautomer_str,
        neutrals=neutral_list,
    )
    return total, rationale


def validate_numbering_convention(atoms, contact_residues_with_names,
                                  md_to_crystal_offset, target_chain):
    """Cross-check that `md_to_crystal_offset` is consistent with PDB identities.

    We can't verify crystal identity without a reference PDB, but we do
    verify per-contact residue identity exists at MD numbering (since the
    target card is MD-based). This is a minimal sanity check.

    Args:
        atoms: parsed PDB atom list (from :func:`parse_pdb`).
        contact_residues_with_names: list of {resnum, resname}.
        md_to_crystal_offset: int (unused at this check, reserved for future
            PDB cross-validation).
        target_chain: chain ID of the target protein.

    Returns:
        list[str]: warnings (empty if clean).
    """
    warnings_out = []
    target_by_rn = {}
    for a in atoms:
        if a["chain"] != target_chain:
            continue
        target_by_rn.setdefault(a["resnum"], a["resname"])
    for r in contact_residues_with_names:
        md_num = r.get("resnum")
        declared = (r.get("resname") or "").upper()
        found = (target_by_rn.get(md_num) or "").upper()
        if not found:
            warnings_out.append(
                "MD {} not present in PDB chain {}".format(md_num, target_chain)
            )
        elif found != declared:
            warnings_out.append(
                "MD {} identity mismatch: card={} pdb={}".format(
                    md_num, declared, found
                )
            )
    return warnings_out


def parse_pdb(pdb_path):
    """Parse ATOM/HETATM records. Returns list of atom dicts."""
    atoms = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            try:
                record = line[:6].strip()
                name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21]
                resnum = int(line[22:26].strip())
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                atoms.append({
                    'record': record,
                    'name': name,
                    'resname': resname,
                    'chain': chain,
                    'resnum': resnum,
                    'x': x, 'y': y, 'z': z,
                })
            except (ValueError, IndexError):
                continue
    return atoms


def extract_contacts(atoms, target_chain, binder_chain, cutoff_A):
    """
    Find target residues within cutoff_A of any binder atom.
    Returns sorted list of (resnum, resname) tuples.
    """
    target_atoms = [a for a in atoms
                    if a['chain'] == target_chain
                    and a['resname'] not in _SOLVENT_RESNAMES]
    binder_atoms = [a for a in atoms
                    if a['chain'] == binder_chain
                    and a['resname'] not in _SOLVENT_RESNAMES]

    if not target_atoms or not binder_atoms:
        raise ValueError(
            f"No atoms found for chain {target_chain!r} or {binder_chain!r}. "
            f"Check --target-chain and --binder-chain."
        )

    binder_xyz = np.array([[a['x'], a['y'], a['z']] for a in binder_atoms])
    contact_residues = {}  # resnum -> resname

    for ta in target_atoms:
        ta_xyz = np.array([ta['x'], ta['y'], ta['z']])
        dists = np.linalg.norm(binder_xyz - ta_xyz, axis=1)
        if dists.min() <= cutoff_A:
            resnum = ta['resnum']
            resname = ta['resname']
            contact_residues[resnum] = resname  # last resname wins (consistent)

    return sorted(contact_residues.items())  # [(resnum, resname), ...]


def extract_contacts_multi(pdb_paths, target_chain, binder_chain, cutoff_A,
                            occupancy_threshold=0.6):
    """
    여러 snapshot 에서 occupancy 기반 contact residues 선택.

    Returns
    -------
    core : list of (resnum, resname)
        occupancy >= threshold 인 잔기 (main contact set).
    extended : list of (resnum, resname)
        0 < occupancy < threshold 인 경계 잔기 (edge / transient contacts).
    occupancy_map : dict[int, float]
        resnum -> fraction of snapshots where contact was observed.
    """
    contact_counts = defaultdict(int)
    contact_resnames = {}
    n_snaps = len(pdb_paths)

    for pdb_path in pdb_paths:
        atoms = parse_pdb(pdb_path)
        contacts = extract_contacts(atoms, target_chain, binder_chain, cutoff_A)
        for resnum, resname in contacts:
            contact_counts[resnum] += 1
            contact_resnames[resnum] = resname  # last wins (same residue across snaps)

    threshold_count = max(1, round(occupancy_threshold * n_snaps))

    core = sorted(
        [(r, contact_resnames[r]) for r, cnt in contact_counts.items()
         if cnt >= threshold_count],
        key=lambda x: x[0]
    )
    extended = sorted(
        [(r, contact_resnames[r]) for r, cnt in contact_counts.items()
         if 0 < cnt < threshold_count],
        key=lambda x: x[0]
    )
    occupancy_map = {r: contact_counts[r] / n_snaps for r in contact_counts}

    return core, extended, occupancy_map


def detect_cofactors(atoms, target_chain):
    """
    Detect non-standard HETATM residues in target chain (excluding solvent).
    Returns list of {resname, resnum, chain} dicts.
    """
    seen = {}
    for a in atoms:
        if (a['record'] == 'HETATM'
                and a['chain'] == target_chain
                and a['resname'] not in _SOLVENT_RESNAMES
                and a['resname'] not in _STANDARD_AA):
            key = (a['chain'], a['resnum'], a['resname'])
            if key not in seen:
                seen[key] = {
                    'resname': a['resname'],
                    'chain': a['chain'],
                    'resnum': a['resnum'],
                    'treatment': 'mm_point_charge',
                }
    return list(seen.values())


def qm_partition_rule(resname):
    """
    Determine QM partitioning rule for a target residue.
    GLY and PRO: whole_residue (backbone coupling or ring fusion)
    Others: sidechain_from_cb
    """
    if resname in ('GLY',):
        return 'whole_residue'
    if resname in ('PRO',):
        return 'whole_residue'
    return 'sidechain_from_cb'


def estimate_qm_atoms(contact_residues, binder_atoms):
    """Rough estimate of QM atom count for budget check."""
    sidechain_avg = 5  # heavy atoms per non-GLY/PRO sidechain
    whole_avg = 7      # heavy atoms per GLY/PRO whole residue
    h_factor = 1.8     # total atoms including H

    target_count = sum(
        whole_avg if resname in ('GLY', 'PRO') else sidechain_avg
        for _, resname in contact_residues
    )
    binder_count = len([a for a in binder_atoms if a['name'] != 'H'])
    return int((target_count + binder_count) * h_factor)


def build_target_card(
    pdb_path, target_chain, binder_chain, cutoff_A,
    target_id, binder_topology, max_n_qm,
    all_pdb_paths=None, occupancy_threshold=0.6, cofactor_mandatory=False,
    md_to_crystal_offset=0,
    numbering_source="MD",
    target_residues_mode="whole",
):
    """Build the full target_card dict.

    Parameters
    ----------
    all_pdb_paths : list[str] or None
        Full list of snapshot PDBs used for occupancy-based contact filtering.
        If None or length <= 1, falls back to single-snapshot mode (legacy).
    occupancy_threshold : float
        Fraction [0,1] of snapshots in which a residue must contact the binder
        to be included in `target_contact_residues`. Residues with
        0 < occupancy < threshold are recorded as `extended_contacts`.
    cofactor_mandatory : bool
        If True, raise when no cofactor HETATM is detected (for cofactor-
        dependent targets such as GTP-bound kinases).
    """
    # Always parse main PDB for cofactor detection and binder atom counting
    atoms_main = parse_pdb(pdb_path)
    cofactors = detect_cofactors(atoms_main, target_chain)
    binder_atoms = [a for a in atoms_main if a['chain'] == binder_chain]

    # 단일 vs 다중 snapshot 처리
    if all_pdb_paths and len(all_pdb_paths) > 1:
        contacts, extended, occ_map = extract_contacts_multi(
            all_pdb_paths, target_chain, binder_chain, cutoff_A, occupancy_threshold
        )
    else:
        contacts_raw = extract_contacts(atoms_main, target_chain, binder_chain, cutoff_A)
        contacts = contacts_raw
        extended = []
        occ_map = {r: 1.0 for r, _ in contacts_raw}
        print(
            "  [WARN] Single snapshot mode — occupancy filter not applied. "
            "Pass --snapshots for multi-snapshot occupancy filtering.",
            file=sys.stderr,
        )

    contact_resnums = [r for r, _ in contacts]
    extended_resnums = [r for r, _ in extended]
    # Schema 0.6.4 (Stage 2): target_residues_mode=whole → entire residue is QM.
    # whole_residue_exceptions list no longer has a functional meaning
    # (it only applied in 0.6.3 sidechain_from_cb mode). Kept as empty
    # for backward-compat field presence.
    whole_exceptions = []

    n_qm_est = estimate_qm_atoms(contacts, binder_atoms)

    # Cofactor handling with mandatory branch
    if not cofactors and cofactor_mandatory:
        raise ValueError(
            f"[cofactor_mandatory=True] No cofactor detected in MD snapshot {pdb_path}. "
            f"For cofactor-dependent targets (e.g. GTP-bound kinases), "
            f"rebuild MD with cofactor included before generating target_card."
        )
    elif cofactors:
        print(f"  [INFO] Cofactors detected: {[c['resname'] for c in cofactors]}", file=sys.stderr)
    else:
        print(
            "  [WARN] No cofactor detected. If target requires cofactor, "
            "use --cofactor-mandatory to enforce.",
            file=sys.stderr,
        )

    n_snapshots = len(all_pdb_paths) if all_pdb_paths else 1

    # Stage 2 (2026-04-19) — schema 0.6.4 auto-rationale with MD/crystal pairs.
    # contacts 는 [(resnum, resname), ...] (extract_contacts 출력 규약).
    contact_res_for_charge = [
        {"resnum": r, "resname": rname} for r, rname in contacts
    ]
    target_iso_charge, _legacy_rationale = compute_target_iso_charge(
        contact_res_for_charge
    )
    _, auto_rationale = generate_target_iso_rationale(
        contact_res_for_charge,
        md_to_crystal_offset=int(md_to_crystal_offset),
        schema_version="0.6.5",
    )

    # Validate numbering convention against PDB identities at each contact.
    nc_warnings = validate_numbering_convention(
        atoms_main, contact_res_for_charge,
        md_to_crystal_offset=int(md_to_crystal_offset),
        target_chain=target_chain,
    )
    for w in nc_warnings:
        print(
            "  [WARN] numbering_convention validator: " + w,
            file=sys.stderr,
        )

    numbering_convention = {
        "source": numbering_source,
        "md_to_crystal_offset": int(md_to_crystal_offset),
        "note": (
            "numbering_source={}; crystal_resnum = md_resnum + {}. "
            "Generated by generate_target_card.py at {}."
        ).format(numbering_source, int(md_to_crystal_offset),
                 str(date.today())),
    }

    # [R-17 v0.6.5] Enrich cofactor entries with R-17 fields so generated
    # cards are canonical v0.6.5. Well-known resnames get a parameter library
    # lookup (GNP → gaff2_resp, MG → li_merz_12_6_4). Unknown resnames emit
    # empty ff_parameters + pdb_literal source; operator must fill in.
    enriched_cofactors = _enrich_cofactor_entries_r17(cofactors) if cofactors else []

    card = {
        "schema_version": "0.6.5",
        "target_id": target_id,
        "target_chain": target_chain,
        "binder_chain": binder_chain,
        "binder_topology": binder_topology,
        "numbering_convention": numbering_convention,
        "qm_partitioning": (
            "whole_residue" if target_residues_mode == "whole"
            else "sidechain_only"
        ),
        "target_contact_residues": contact_resnums,
        "extended_contacts": extended_resnums,
        # R-16 hint — runtime PDB is the ground truth; R-16 guard fails
        # fast when this hint does not match. Default 0 is conservative.
        "binder_net_charge": 0,
        "binder_charge_note": (
            "R-16 hint — runtime PDB chemistry is the source of truth. "
            "Override with actual binder net charge if known (e.g. -1 for "
            "cyclic ARG+ASP+GLU)."
        ),
        "target_iso_net_charge": target_iso_charge,
        "target_iso_charge_rationale": auto_rationale,
        "cofactor_residues": enriched_cofactors,
        "partition_rules": {
            "binder_residues_mode": "whole",
            "target_residues_mode": target_residues_mode,
            "glycine_handling": "whole_residue",
            "proline_handling": "whole_residue",
            "link_atom_model": "senn_thiel_2009_hcap",
            "cb_link_bond_length_A": 1.09,
        },
        "qm_atom_budget": {
            "max_n_qm": max_n_qm,
            "expected_n_qm_estimate": n_qm_est,
            "enforcement": "hard_fail",
        },
        "provenance": {
            "generated_by": "generate_target_card.py",
            "schema_version": "0.6.5",
            "source_pdb": os.path.abspath(pdb_path),
            "contact_cutoff_A": cutoff_A,
            "generated_date": str(date.today()),
            "contact_residues_count": len(contact_resnums),
            "n_snapshots_used": n_snapshots,
            "occupancy_threshold": occupancy_threshold,
            "occupancy_map": {str(r): round(v, 2) for r, v in occ_map.items()},
            "numbering_validation_warnings": nc_warnings,
        },
        "backbone_keep_residues": [],
        "notes": (
            f"Auto-generated from MD snapshot{'s' if n_snapshots > 1 else ''}. "
            f"{len(contact_resnums)} core contacts "
            f"(occupancy >= {occupancy_threshold}, n_snap={n_snapshots}) "
            f"at {cutoff_A}Å cutoff. "
            f"{len(extended_resnums)} extended (transient) contacts. "
            f"n_qm estimate: {n_qm_est}. target_residues_mode={target_residues_mode}. "
            + ("Cofactor absent in MD snapshot — verify MD prep." if not cofactors else "")
        ),
    }
    return card


def main():
    parser = argparse.ArgumentParser(
        description="Generate UPDD target_card.json from MD snapshot PDB."
    )
    parser.add_argument('--pdb', required=True,
                        help='MD snapshot PDB file (ATOM/HETATM)')
    parser.add_argument('--target-chain', required=True,
                        help='Chain ID of the target protein (e.g. A)')
    parser.add_argument('--binder-chain', required=True,
                        help='Chain ID of the peptide binder (e.g. B)')
    parser.add_argument('--cutoff', type=float, default=4.5,
                        help='Heavy-atom contact cutoff in Angstrom (default: 4.5)')
    parser.add_argument('--target-id', required=True,
                        help='Target identifier (e.g. 6WGN). Used for output filename.')
    parser.add_argument('--binder-topology', default='linear',
                        choices=['linear', 'cyclic_htc', 'cyclic_ss', 'cyclic_cc'],
                        help='Binder ring topology (default: linear)')
    parser.add_argument('--max-n-qm', type=int, default=600,
                        help='Hard QM atom budget limit (default: 600)')
    parser.add_argument('--output', default=None,
                        help='Output path (default: target_cards/<target-id>.json)')
    parser.add_argument('--snapshots', nargs='+', default=None,
                        help='추가 snapshot PDB 파일들 (occupancy 계산용). '
                             '없으면 --pdb 단독 사용.')
    parser.add_argument('--occupancy-threshold', type=float, default=0.6,
                        help='Contact occupancy 임계값 (기본 0.6 = 5 snap 중 3개 이상).')
    parser.add_argument('--cofactor-mandatory', action='store_true', default=False,
                        help='cofactor 가 MD 에 없으면 오류 발생 (GTP-bound 타겟 등).')
    parser.add_argument('--md-to-crystal-offset', type=int, default=0,
                        help=('crystal_resnum = md_resnum + offset '
                              '(schema 0.6.4 numbering_convention). '
                              '예: 6WGN MD→crystal = +3.'))
    parser.add_argument('--numbering-source', default='MD',
                        choices=['MD', 'crystal', 'Uniprot'],
                        help='target_contact_residues 의 좌표계 (기본 MD).')
    parser.add_argument('--target-residues-mode', default='whole',
                        choices=['whole', 'sidechain_from_cb'],
                        help=('target contact residue QM 포함 방식 '
                              '(Stage 2 default = whole; legacy = sidechain_from_cb).'))
    args = parser.parse_args()

    if not os.path.exists(args.pdb):
        sys.exit(f"[ERROR] PDB not found: {args.pdb}")

    # --snapshots 와 --pdb 합치기 (중복 제거, --pdb 는 반드시 포함)
    all_pdbs = [args.pdb]
    if args.snapshots:
        for p in args.snapshots:
            if p != args.pdb and p not in all_pdbs:
                if not os.path.exists(p):
                    sys.exit(f"[ERROR] Snapshot PDB not found: {p}")
                all_pdbs.append(p)

    out_path = args.output or os.path.join(
        os.path.dirname(__file__), '..', 'target_cards', f"{args.target_id}.json"
    )
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)

    print(f"[generate_target_card] {args.target_id}", file=sys.stderr)
    print(f"  PDB        : {args.pdb}", file=sys.stderr)
    print(f"  chains     : target={args.target_chain}, binder={args.binder_chain}", file=sys.stderr)
    print(f"  cutoff     : {args.cutoff} Å", file=sys.stderr)
    print(f"  topology   : {args.binder_topology}", file=sys.stderr)
    print(f"  snapshots  : {len(all_pdbs)} (occupancy threshold={args.occupancy_threshold})",
          file=sys.stderr)
    if args.cofactor_mandatory:
        print(f"  cofactor   : MANDATORY (fail-fast if absent)", file=sys.stderr)

    card = build_target_card(
        pdb_path=args.pdb,
        target_chain=args.target_chain,
        binder_chain=args.binder_chain,
        cutoff_A=args.cutoff,
        target_id=args.target_id,
        binder_topology=args.binder_topology,
        max_n_qm=args.max_n_qm,
        all_pdb_paths=all_pdbs,
        occupancy_threshold=args.occupancy_threshold,
        cofactor_mandatory=args.cofactor_mandatory,
        md_to_crystal_offset=args.md_to_crystal_offset,
        numbering_source=args.numbering_source,
        target_residues_mode=args.target_residues_mode,
    )

    with open(out_path, 'w') as f:
        json.dump(card, f, indent=2)

    n = len(card['target_contact_residues'])
    est = card['qm_atom_budget']['expected_n_qm_estimate']
    mode = card['partition_rules']['target_residues_mode']
    print(f"  contacts   : {n} residues", file=sys.stderr)
    print(f"  mode       : target_residues_mode={mode}", file=sys.stderr)
    print(f"  n_qm est   : {est}", file=sys.stderr)
    print(f"  output     : {out_path}", file=sys.stderr)
    print(f"[generate_target_card] DONE", file=sys.stderr)

    # Also print the card as JSON to stdout (for piping/inspection)
    print(json.dumps(card, indent=2))


# ==========================================
# [R-17 v0.6.5] Cofactor schema enrichment + completeness validation
# ==========================================
# Known-cofactor parameter library. Lookup by resname → default
# ``ff_parameters`` template. Cards for new targets simply reference these
# paths and the pipeline (G3 parameterize_cofactor) resolves them at runtime.
#
# Rationale sources cite the primary parameterization literature so cards
# produced here are reviewable without external context.
_COFACTOR_FF_LIBRARY = {
    "GNP": {
        "method": "gaff2_resp",
        "frcmod": "params/gnp/gnp.frcmod",
        "mol2":   "params/gnp/gnp.mol2",
        "reference": "Bayly 1993 JPC 97,10269 + AMBER GAFF2 (Wang 2004 JCC 25,1157)",
    },
    "GTP": {
        "method": "gaff2_resp",
        "frcmod": "params/gtp/gtp.frcmod",
        "mol2":   "params/gtp/gtp.mol2",
        "reference": "Meagher 2003 JCC 24,1016 (phosphate parameters)",
    },
    "GDP": {
        "method": "gaff2_resp",
        "frcmod": "params/gdp/gdp.frcmod",
        "mol2":   "params/gdp/gdp.mol2",
        "reference": "Meagher 2003 JCC 24,1016",
    },
    "ATP": {
        "method": "gaff2_resp",
        "frcmod": "params/atp/atp.frcmod",
        "mol2":   "params/atp/atp.mol2",
        "reference": "Meagher 2003 JCC 24,1016",
    },
    "ADP": {
        "method": "gaff2_resp",
        "frcmod": "params/adp/adp.frcmod",
        "mol2":   "params/adp/adp.mol2",
        "reference": "Meagher 2003 JCC 24,1016",
    },
    "MG":  {"method": "li_merz_12_6_4", "reference": "Li & Merz 2014 JCTC 10,1518"},
    "ZN":  {"method": "li_merz_12_6_4", "reference": "Li & Merz 2014 JCTC 10,1518"},
    "CA":  {"method": "li_merz_12_6_4", "reference": "Li & Merz 2014 JCTC 10,1518"},
    "MN":  {"method": "li_merz_12_6_4", "reference": "Li & Merz 2014 JCTC 10,1518"},
    "FE":  {"method": "li_merz_12_6_4", "reference": "Li & Merz 2014 JCTC 10,1518"},
}

# Treatment defaults when the target_card generator detects a cofactor class.
# These are DEFAULTS — operators should review per-target (e.g. catalytic Zn
# typically wants treatment="qm_region" in v0.7+).
_DEFAULT_TREATMENT_BY_CLASS = {
    "nucleotide_triphosphate": "mm_point_charge",      # GNP, GTP, ATP
    "nucleotide_diphosphate":  "mm_point_charge",      # GDP, ADP
    "divalent_ion":            "mm_point_charge",      # MG, CA, MN, FE (non-cat)
    "transition_metal":        "mm_point_charge",      # ZN (cat → qm_region v0.7+)
}

_COFACTOR_CLASS_BY_RESNAME = {
    "GNP": "nucleotide_triphosphate",
    "GTP": "nucleotide_triphosphate",
    "ATP": "nucleotide_triphosphate",
    "GDP": "nucleotide_diphosphate",
    "ADP": "nucleotide_diphosphate",
    "MG":  "divalent_ion",
    "CA":  "divalent_ion",
    "MN":  "divalent_ion",
    "FE":  "divalent_ion",
    "ZN":  "transition_metal",
}


def _enrich_cofactor_entries_r17(cofactor_seeds):
    """[R-17 v0.6.5] Produce full cofactor_residues entries from seeds.

    Args:
        cofactor_seeds: list of dicts carrying (at minimum) ``resname``,
            ``chain``, ``resnum``, and optionally ``charge``. These come
            from the MD PDB scan earlier in the pipeline.

    Returns:
        list of dicts with R-17 v0.6.5 fields filled in. Well-known resnames
        resolve against the internal library; unknown resnames get
        ``ff_parameters={}`` + ``source='pdb_literal'`` so the operator must
        supply parameters before G3 will pass.
    """
    enriched = []
    for c in (cofactor_seeds or []):
        resname = str(c.get("resname", "")).strip().upper()
        entry = {
            "resname": resname,
            "chain":   c.get("chain"),
            "resnum":  c.get("resnum"),
            "required": True,
            "source":  "pdb_literal",
            "treatment": "mm_point_charge",
            "alternative_treatments": ["mm_point_charge"],
            "charge":  c.get("charge", 0),
            "ff_parameters": {},
        }
        # Library lookup for known cofactors.
        if resname in _COFACTOR_FF_LIBRARY:
            entry["ff_parameters"] = dict(_COFACTOR_FF_LIBRARY[resname])
        # Class-driven treatment override (divalent ions can also go qm_region
        # in v0.7+, flagged via alternative_treatments).
        klass = _COFACTOR_CLASS_BY_RESNAME.get(resname)
        if klass:
            entry["treatment"] = _DEFAULT_TREATMENT_BY_CLASS.get(
                klass, entry["treatment"]
            )
            if klass in ("divalent_ion", "transition_metal"):
                entry["alternative_treatments"] = [
                    "mm_point_charge", "qm_region",
                ]
        enriched.append(entry)
    return enriched


def validate_cofactor_completeness(card, repo_root=None):
    """[R-17 v0.6.5] Verify that every required cofactor's ff_parameters
    resolve to existing files on disk.

    Args:
        card: target_card dict (v0.6.5 shape).
        repo_root: optional absolute path to the project root. Defaults to
            ``<utils>/../``.

    Returns:
        dict with keys:
            ``ok``     (bool) — True iff every required entry is satisfied.
            ``issues`` (list[str]) — per-entry remediation hints.
            ``entries`` (list[dict]) — per-entry diagnostic.
    """
    if repo_root is None:
        repo_root = os.path.abspath(
            os.path.join(os.path.dirname(__file__), os.pardir)
        )
    issues = []
    report = []
    for entry in (card.get("cofactor_residues") or []):
        if not isinstance(entry, dict):
            continue
        if not entry.get("required"):
            continue
        resname = str(entry.get("resname", "")).strip().upper()
        ff = entry.get("ff_parameters") or {}
        method = str(ff.get("method", "")).strip().lower()
        entry_report = {
            "resname": resname, "method": method,
            "cache_hit": False, "files_checked": [],
        }
        # Ion-builtin methods need no per-target frcmod.
        if method in ("li_merz_12_6_4", "li_merz", "amber_ion", "ion_builtin"):
            entry_report["cache_hit"] = True
            entry_report["status"] = "ion_builtin"
            report.append(entry_report)
            continue
        # Any other method must point at real files (frcmod / mol2 / xml).
        found_any = False
        for key in ("frcmod", "mol2", "xml"):
            rel_path = ff.get(key)
            if not rel_path:
                continue
            abs_path = rel_path if os.path.isabs(rel_path) else os.path.join(
                repo_root, rel_path
            )
            exists = os.path.exists(abs_path)
            entry_report["files_checked"].append({"key": key, "path": abs_path,
                                                  "exists": exists})
            if exists:
                found_any = True
        entry_report["cache_hit"] = found_any
        entry_report["status"] = "cache_hit" if found_any else "missing"
        if not found_any:
            issues.append(
                "cofactor '{}' (method={}): no existing ff_parameters "
                "files. Regenerate via parameterize_cofactor.ensure_cofactor_params "
                "or ship the frcmod/mol2 under params/{}/.".format(
                    resname, method or "unspecified", resname.lower()
                )
            )
        report.append(entry_report)
    return {"ok": not issues, "issues": issues, "entries": report}


if __name__ == '__main__':
    main()
