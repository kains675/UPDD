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
    all_pdb_paths=None, occupancy_threshold=0.6, cofactor_mandatory=False
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
    whole_exceptions = [r for r, rname in contacts
                        if qm_partition_rule(rname) == 'whole_residue']

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

    # [v0.6.5 Phase 4] target_iso_net_charge — Phase 4 (qm_int_kcal_frozen) 활성화용.
    # contacts 는 [(resnum, resname), ...] 튜플 리스트 (extract_contacts 출력 규약).
    contact_res_for_charge = [
        {"resnum": r, "resname": rname} for r, rname in contacts
    ]
    target_iso_charge, target_iso_rationale = compute_target_iso_charge(
        contact_res_for_charge
    )

    card = {
        "schema_version": "0.6.3",
        "target_id": target_id,
        "target_chain": target_chain,
        "binder_chain": binder_chain,
        "binder_topology": binder_topology,
        "numbering_convention": "md_snapshot",
        "qm_partitioning": "sidechain_only",
        "target_contact_residues": contact_resnums,
        "whole_residue_exceptions": whole_exceptions,
        "extended_contacts": extended_resnums,
        # [v0.6.5 Phase 4] 초분자 분해용 전하 정보 (pH 7.4 approx).
        "binder_net_charge": 0,
        "binder_charge_note": (
            "computed at runtime from binder QM atom proton parity "
            "(default 0 if even electrons)"
        ),
        "target_iso_net_charge": target_iso_charge,
        "target_iso_charge_rationale": target_iso_rationale,
        "cofactor_residues": [
            dict(c, treatment="mm_point_charge") for c in cofactors
        ] if cofactors else [],
        "partition_rules": {
            "binder_residues_mode": "whole",
            "target_residues_mode": "sidechain_from_cb",
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
            "source_pdb": os.path.abspath(pdb_path),
            "contact_cutoff_A": cutoff_A,
            "generated_date": str(date.today()),
            "contact_residues_count": len(contact_resnums),
            "n_snapshots_used": n_snapshots,
            "occupancy_threshold": occupancy_threshold,
            "occupancy_map": {str(r): round(v, 2) for r, v in occ_map.items()},
            "whole_residue_exceptions_resnames": {
                str(r): rname for r, rname in contacts
                if r in whole_exceptions
            },
        },
        "backbone_keep_residues": [],
        "notes": (
            f"Auto-generated from MD snapshot{'s' if n_snapshots > 1 else ''}. "
            f"{len(contact_resnums)} core contacts "
            f"(occupancy >= {occupancy_threshold}, n_snap={n_snapshots}) "
            f"at {cutoff_A}Å cutoff. "
            f"{len(extended_resnums)} extended (transient) contacts. "
            f"n_qm estimate: {n_qm_est}. "
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
    )

    with open(out_path, 'w') as f:
        json.dump(card, f, indent=2)

    n = len(card['target_contact_residues'])
    exc = card['whole_residue_exceptions']
    est = card['qm_atom_budget']['expected_n_qm_estimate']
    print(f"  contacts   : {n} residues", file=sys.stderr)
    print(f"  whole-res  : {exc} (GLY/PRO)", file=sys.stderr)
    print(f"  n_qm est   : {est}", file=sys.stderr)
    print(f"  output     : {out_path}", file=sys.stderr)
    print(f"[generate_target_card] DONE", file=sys.stderr)

    # Also print the card as JSON to stdout (for piping/inspection)
    print(json.dumps(card, indent=2))


if __name__ == '__main__':
    main()
