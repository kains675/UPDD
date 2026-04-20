#!/usr/bin/env python
"""
ncaa_mutate.py
--------------
PDB에 ncAA 치환 적용 (AUTO_SASA 지원).
UPDD의 ncaa_registry와 연동하여 하드코딩을 제거하고 Manifest.json을 출력합니다.

[Architecture Zenith Update]
- CLI 인자를 --ncaa_code 로 표준화하여 UPDD.py와 완벽 동기화
- [CRITICAL FIX] AUTO_SASA 탐색 시 펩타이드가 깊게 결합되어 BSA 필터를 통과하는 후보가 0개일 때, 필터를 완전히 무시하고 가장 표면에 노출된 잔기를 강제 선택하는 3단계 폴백(Ultimate Fallback) 완벽 복원
"""

import os
import sys
import glob
import copy
import shutil  # [v35 AUDIT] Principle 3: required by the fallback pdb-copy path in main() (line ~297) — prior absence produced a latent NameError on every non-substituted design.
import argparse
import json
import logging
import numpy as np

# utils 디렉토리 내부에서 직접 registry 참조
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    from ncaa_registry import resolve_ncaa_definition
except ImportError:
    sys.exit("[!] ncaa_registry.py 파일을 찾을 수 없습니다. utils 폴더에 위치해야 합니다.")

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("ncaa_mutate")

MAX_SASA = {
    "ALA": 121.0, "ARG": 265.0, "ASN": 187.0, "ASP": 187.0,
    "CYS": 148.0, "GLU": 214.0, "GLN": 214.0, "GLY":  97.0,
    "HIS": 216.0, "ILE": 195.0, "LEU": 191.0, "LYS": 230.0,
    "MET": 203.0, "PHE": 228.0, "PRO": 154.0, "SER": 143.0,
    "THR": 163.0, "TRP": 264.0, "TYR": 255.0, "VAL": 165.0,
}

def write_manifest(output_dir, args, resolved_targets):
    manifest = {
        "input_positions": args.residues,
        "resolved_targets": [{"chain": c, "resid": r} for c, r in resolved_targets],
        "ncaa_code": args.ncaa_code,
        "status": "SUCCESS"
    }
    # [v35 AUDIT] Principle 7/v27: encoding="utf-8" — output_dir may contain 한글 path segments.
    with open(os.path.join(output_dir, "mutation_manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

def _make_cm_atom_line(serial: int, chain_id: str, res_num: int, res_name_padded: str, n_coord: tuple, ca_coord: tuple) -> str:
    """N-메틸 ncAA(NMA) 의 추가 CM 탄소 원자에 대한 PDB HETATM 레코드를 생성한다.

    [v4 CRITICAL FIX] 이전 구현은 f-string 으로 `{res_name_padded}{chain_id}` 를
    직접 이어 붙여 PDB col 21 의 blank 를 생략했다. 이로 인해 chain_id 가 col 21
    (index 20) 에 배치되고, 파서가 col 22 (index 21) 로 chain 을 읽을 때 항상 blank
    가 반환되어 CM 원자가 "chain=' '" 의 유령 잔기로 분리되었다. 이는 OpenMM 의
    addHydrogens() 가 NMA 잔기를 ALA-like (5 heavy atoms) 로 오인하고 "No template
    found for residue NMA" 오류로 즉시 실패하는 근본 원인이었다.

    본 수정은 80-char PDB ATOM/HETATM 레코드를 col-exact 위치로 조립하여 포맷
    오류의 재발을 차단한다. 필드 레이아웃 (1-based col / 0-based Python index):

        - HETATM       cols 1-6   / [0:6]
        - serial       cols 7-11  / [6:11]     right-justified
        - blank        col 12     / [11]
        - atom name    cols 13-16 / [12:16]    " CM " (C 원소, 1-letter padding)
        - altLoc       col 17     / [16]       blank
        - resName      cols 18-20 / [17:20]    "NMA" left-justified
        - blank        col 21     / [20]       ← 이전 구현이 생략했던 위치
        - chainID      col 22     / [21]       chain_id 의 첫 글자
        - resSeq       cols 23-26 / [22:26]    right-justified
        - iCode        col 27     / [26]       blank
        - blank        cols 28-30 / [27:30]
        - x            cols 31-38 / [30:38]
        - y            cols 39-46 / [38:46]
        - z            cols 47-54 / [46:54]
        - occupancy    cols 55-60 / [54:60]
        - tempFactor   cols 61-66 / [60:66]
        - element      cols 77-78 / [76:78]    right-justified " C"
    """
    n  = np.array(n_coord,  dtype=float)
    ca = np.array(ca_coord, dtype=float)
    v  = n - ca
    norm = np.linalg.norm(v)
    if norm < 1e-6:
        v = np.array([1.0, 0.0, 0.0])
        norm = 1.0
    v /= norm
    cm = n + v * 1.46

    # PDB 레코드를 col-exact 로 조립 (80-char 고정폭)
    line_buf = list(" " * 80)
    line_buf[0:6]   = list("HETATM")
    line_buf[6:11]  = list(f"{serial:5d}"[-5:])
    line_buf[12:16] = list(" CM ")                         # 탄소(1-letter element) 명명 규약
    line_buf[17:20] = list(f"{res_name_padded[:3]:<3}")
    line_buf[21]    = chain_id[0] if chain_id else " "     # ← 올바른 col 22 위치
    line_buf[22:26] = list(f"{res_num:4d}"[-4:])
    line_buf[30:38] = list(f"{cm[0]:8.3f}"[-8:])
    line_buf[38:46] = list(f"{cm[1]:8.3f}"[-8:])
    line_buf[46:54] = list(f"{cm[2]:8.3f}"[-8:])
    line_buf[54:60] = list("  1.00")
    line_buf[60:66] = list("  0.00")
    line_buf[76:78] = list(" C")                           # element 우측 정렬 (col 77-78)
    return "".join(line_buf).rstrip() + "\n"

def _compute_complex_and_unbound_sasa(pdb_path, binder_chain, target_chain):
    """[v3 9-1] 복합체 및 단독 binder 의 잔기별 SASA 를 BioPython ShrakeRupley 로 계산한다.

    Returns:
        Tuple[BioPython Structure, dict, dict]:
            (struct_complex, complex_sasa, unbound_sasa)
            - complex_sasa: {res_num: sasa} — binder_chain 잔기만 (복합체 환경)
            - unbound_sasa: {res_num: sasa} — target_chain 제거 후 binder 단독 SASA.
              target chain 이 존재하지 않으면 빈 dict.
    """
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.SASA import ShrakeRupley
    except ImportError:
        raise RuntimeError("BioPython이 설치되어 있지 않습니다. (conda install -c conda-forge biopython)")

    parser = PDBParser(QUIET=True)
    struct_complex = parser.get_structure("complex", pdb_path)
    sr = ShrakeRupley()

    sr.compute(struct_complex, level="R")
    complex_sasa: dict = {}
    for chain in struct_complex[0].get_chains():
        if chain.id == binder_chain:
            for res in chain.get_residues():
                if res.get_id()[0] == " ":
                    complex_sasa[res.get_id()[1]] = res.sasa

    struct_binder_only = copy.deepcopy(struct_complex)
    try:
        struct_binder_only[0].detach_child(target_chain)
    except KeyError:
        struct_binder_only = None

    unbound_sasa: dict = {}
    if struct_binder_only is not None:
        sr.compute(struct_binder_only, level="R")
        for chain in struct_binder_only[0].get_chains():
            for res in chain.get_residues():
                if res.get_id()[0] == " ":
                    unbound_sasa[res.get_id()[1]] = res.sasa

    return struct_complex, complex_sasa, unbound_sasa


def _residue_metric(res, complex_sasa, unbound_sasa):
    """단일 잔기 → (res_num, rsasa, plddt, bsa) 튜플 (또는 None)."""
    if res.get_id()[0] != " " or "CA" not in res:
        return None
    res_num = res.get_id()[1]
    res_name = res.get_resname()
    plddt = res["CA"].get_bfactor()
    sasa_abs = complex_sasa.get(res_num, 0.0)
    max_ref = MAX_SASA.get(res_name, 200.0)
    rsasa = sasa_abs / max_ref if max_ref > 0 else 0.0
    bsa = (unbound_sasa.get(res_num, 0.0) - sasa_abs) if unbound_sasa else 0.0
    return res_num, rsasa, plddt, bsa


def _identify_terminal_resids(binder_chain_obj):
    """[v5 해법 A] chain 의 첫/마지막 표준 잔기 ID(set) 를 반환한다.

    근거 (CLAUDE.md v5 §해법 A):
        ncAA 가 펩타이드 chain 의 N-말단 또는 C-말단에 배치되면 OpenMM
        addHydrogens() 가 자유산/자유아민 형태의 terminal 템플릿(N{RES}/C{RES})
        을 요구한다. 현재 parameterize_ncaa.py 가 출력하는 XML 은 internal 잔기
        템플릿 한 종류이므로 (해법 B 적용 전), 말단 배치 시 즉시 template
        매칭이 실패한다. 본 헬퍼는 AUTO_SASA 후보 산출 단계에서 chain 의
        첫/마지막 표준 잔기를 사전에 배제하기 위한 키 집합을 반환한다.

        말단 잔기는 구조적으로 유연(고 RMSF)하여 결합 자유에너지 기여도가
        본질적으로 낮으므로, 후보에서 배제하더라도 데이터 신뢰성이 저하되지
        않는다 (Principle 1: 신뢰성 향상 유지).

    Args:
        binder_chain_obj: BioPython Chain 객체.

    Returns:
        set[int]: 첫/마지막 표준 잔기의 resid. 표준 잔기가 2개 미만이면 빈 set.
    """
    standard_residues = [r for r in binder_chain_obj.get_residues() if r.get_id()[0] == " "]
    if len(standard_residues) < 2:
        return set()
    return {standard_residues[0].get_id()[1], standard_residues[-1].get_id()[1]}


def _apply_strict_filter(binder_chain_obj, complex_sasa, unbound_sasa,
                         plddt_cutoff=70.0, bsa_exclusion=20.0):
    """[v3 9-1] 1차 엄격 필터: pLDDT ≥ plddt_cutoff, BSA ≤ bsa_exclusion.

    [v5 해법 A] chain 첫/마지막 표준 잔기는 후보에서 제외한다.
    """
    terminal_resids = _identify_terminal_resids(binder_chain_obj)
    candidates = []
    for res in binder_chain_obj.get_residues():
        metric = _residue_metric(res, complex_sasa, unbound_sasa)
        if metric is None:
            continue
        res_num, rsasa, plddt, bsa = metric
        if res_num in terminal_resids:
            continue
        if plddt < plddt_cutoff:
            continue
        if unbound_sasa and bsa > bsa_exclusion:
            continue
        candidates.append((res_num, rsasa, plddt))
    return candidates


def _apply_relaxed_filter(binder_chain_obj, complex_sasa, unbound_sasa, bsa_exclusion=20.0):
    """[v3 9-1] 2차 완화 필터: 두 단계 (60/bsa_excl, 50/2.5*bsa_excl) 순차 적용.

    [v5 해법 A] chain 첫/마지막 표준 잔기는 후보에서 제외한다.
    """
    terminal_resids = _identify_terminal_resids(binder_chain_obj)
    for fb_plddt, fb_bsa in [(60.0, bsa_exclusion), (50.0, bsa_exclusion * 2.5)]:
        fb_candidates = []
        for res in binder_chain_obj.get_residues():
            metric = _residue_metric(res, complex_sasa, unbound_sasa)
            if metric is None:
                continue
            res_num, rsasa, plddt, bsa = metric
            if res_num in terminal_resids:
                continue
            if plddt < fb_plddt:
                continue
            if unbound_sasa and bsa > fb_bsa:
                continue
            fb_candidates.append((res_num, rsasa, plddt))
        if fb_candidates:
            return fb_candidates
    return []


def _apply_sasa_fallback(binder_chain_obj, complex_sasa, unbound_sasa):
    """[v3 9-1] 3차 폴백: 모든 필터를 무시하고 rSASA 만으로 후보 추출.

    [v5 해법 A] chain 첫/마지막 표준 잔기는 폴백에서도 제외한다 — 말단 배치는
    OpenMM 템플릿 매칭 자체가 실패하므로 데이터 신뢰성이 아닌 파이프라인
    동작 가능성의 문제이며, 폴백 단계에서도 동일하게 우회되어야 한다.
    """
    log.warning("[AUTO_SASA] 엄격한 필터(pLDDT, BSA)를 만족하는 타겟이 없어, 필터를 무시하고 표면 노출도가 가장 높은 잔기를 강제로 선택합니다.")
    terminal_resids = _identify_terminal_resids(binder_chain_obj)
    fallback = []
    for res in binder_chain_obj.get_residues():
        metric = _residue_metric(res, complex_sasa, unbound_sasa)
        if metric is None:
            continue
        res_num, rsasa, plddt, _bsa = metric
        if res_num in terminal_resids:
            continue
        fallback.append((res_num, rsasa, plddt))
    return fallback


def _select_by_bsa_and_plddt(pdb_path, n, binder_chain, target_chain,
                             plddt_cutoff=70.0, bsa_exclusion=20.0):
    """[v3 9-1] 4단계 필터(엄격 → 완화 → 폴백) 를 순차 적용하는 얇은 오케스트레이터.

    외부 인터페이스 (인자, 반환값) 는 기존과 완벽히 동일하다.
    """
    _struct, complex_sasa, unbound_sasa = _compute_complex_and_unbound_sasa(
        pdb_path, binder_chain, target_chain
    )
    binder_chain_obj = _struct[0][binder_chain]

    candidates = _apply_strict_filter(
        binder_chain_obj, complex_sasa, unbound_sasa, plddt_cutoff, bsa_exclusion
    )
    if not candidates:
        candidates = _apply_relaxed_filter(
            binder_chain_obj, complex_sasa, unbound_sasa, bsa_exclusion
        )
    if not candidates:
        candidates = _apply_sasa_fallback(binder_chain_obj, complex_sasa, unbound_sasa)

    if not candidates:
        # 이 에러가 난다면 진짜로 PDB 체인 안에 아미노산이 1개도 없는 비정상 파일입니다.
        raise RuntimeError(f"바인더 체인({binder_chain}) 내에 평가할 수 있는 표준 아미노산 잔기가 전혀 없습니다.")

    # rSASA(용매 노출 비율)가 가장 높은 순서대로 정렬하여 추출
    candidates.sort(key=lambda x: x[1], reverse=True)
    return [(binder_chain, item[0]) for item in candidates[:n]]

def parse_residue_targets(residues_str, input_dir=None, binder_chain="B", target_chain="A", plddt_cutoff=70.0, bsa_exclusion=20.0):
    residues_str = residues_str.strip()
    if residues_str.upper().startswith("AUTO_SASA:"):
        if not input_dir or not os.path.exists(input_dir):
            raise FileNotFoundError(f"[!] SASA 계산을 위한 레퍼런스 폴더를 찾을 수 없습니다: {input_dir}")
            
        pdb_files = sorted(glob.glob(os.path.join(input_dir, "*.pdb")))
        if not pdb_files:
            raise FileNotFoundError(f"[!] {input_dir} 폴더에 PDB 파일이 없습니다. 이전 단계(AlphaFold2)가 정상적으로 완료되었는지 확인하세요.")
            
        rank001 = [p for p in pdb_files if "rank_001" in os.path.basename(p)]
        ref_pdb = rank001[0] if rank001 else pdb_files[0]
        n = int(residues_str.split(":")[1])
        return _select_by_bsa_and_plddt(ref_pdb, n, binder_chain, target_chain, plddt_cutoff, bsa_exclusion)
    
    targets = []
    for token in residues_str.split(","):
        token = token.strip()
        if not token: continue
        if token[0].isalpha() and token[1:].lstrip("-").isdigit():
            targets.append((token[0].upper(), int(token[1:])))
        elif token.lstrip("-").isdigit():
            targets.append((binder_chain, int(token)))
    return targets

def _generate_cb_from_backbone(n_coord, ca_coord, c_coord):
    """N, CA, C 좌표로부터 CB 위치를 이상적 정사면체 기하학으로 계산한다.

    Engh & Huber (1991): CA-CB = 1.521 Å, N-CA-CB = 110.5°.
    CB 는 N 과 C 의 반대쪽 (L-amino acid 입체배치) 에 위치한다.
    GLY → NMA 치환 시 원본에 CB 가 없으므로 이 함수로 자동 생성한다.
    """
    import numpy as _np
    n  = _np.array(n_coord,  dtype=float)
    ca = _np.array(ca_coord, dtype=float)
    c  = _np.array(c_coord,  dtype=float)
    ca_n = n - ca
    ca_c = c - ca
    n_can = ca_n / max(_np.linalg.norm(ca_n), 1e-8)
    n_cac = ca_c / max(_np.linalg.norm(ca_c), 1e-8)
    cb_dir = -(n_can + n_cac)
    cb_dir = cb_dir / max(_np.linalg.norm(cb_dir), 1e-8)
    cb = ca + cb_dir * 1.521
    return (float(cb[0]), float(cb[1]), float(cb[2]))


# [v0.6.7] Generic parent-topology helpers moved to utils.parent_topology.
# The legacy hardcoded _PARENT_NEIGHBOR_HINT / _PARENT_RESIDUE_BONDS dicts are
# replaced by functions that query amber14 ff14SB.xml templates for any parent.
# Import lazy to avoid hard dependency at module load if parent_topology itself
# fails to import (e.g., networkx not available).

def _get_neighbor_hint(parent_residue: str, attach_name: str):
    try:
        from parent_topology import parent_neighbor_hint
        return parent_neighbor_hint(parent_residue, attach_name)
    except ImportError:
        return ()


def _get_parent_bonds(parent_residue: str):
    try:
        from parent_topology import parent_heavy_bonds
        return parent_heavy_bonds(parent_residue)
    except ImportError:
        return []


def _place_extension_atom_sp2(attach_coord, neighbor1_coord, neighbor2_coord, bond_length):
    """sp2 attach atom (예: TRP NE1) 에서 두 ring neighbor 의 bisector 반대 방향
    으로 extension atom 을 배치한다. TRP NE1 의 경우 CM 이 indole 평면 위,
    ring 바깥쪽으로 1.466 Å (amber N-alkyl bond).

    Ref: Capece 2012 (1-Me-Trp), 정상값 확인용 SMILES 3D 기하 비교.
    """
    import numpy as _np
    a  = _np.array(attach_coord,   dtype=float)
    n1 = _np.array(neighbor1_coord, dtype=float)
    n2 = _np.array(neighbor2_coord, dtype=float)
    v1 = a - n1
    v2 = a - n2
    v1 = v1 / max(_np.linalg.norm(v1), 1e-8)
    v2 = v2 / max(_np.linalg.norm(v2), 1e-8)
    bisector = v1 + v2
    bisector = bisector / max(_np.linalg.norm(bisector), 1e-8)
    ext = a + bisector * bond_length
    return (float(ext[0]), float(ext[1]), float(ext[2]))


def _make_extension_atom_line(serial, chain_id, res_num, res_name_padded, atom_name, element, coord):
    """[v0.6.6] extension_atoms 전용 PDB ATOM line writer."""
    buf = list(" " * 80)
    buf[0:6]   = list("HETATM")
    buf[6:11]  = list(f"{serial:5d}")
    buf[12:16] = list(f" {atom_name:<3s}" if len(atom_name) <= 3 else atom_name[:4])
    buf[17:20] = list(res_name_padded[:3])
    buf[21]    = chain_id
    buf[22:26] = list(f"{res_num:4d}")
    buf[30:38] = list(f"{coord[0]:8.3f}")
    buf[38:46] = list(f"{coord[1]:8.3f}")
    buf[46:54] = list(f"{coord[2]:8.3f}")
    buf[54:60] = list("  1.00")
    buf[60:66] = list("  0.00")
    sym = element.rjust(2)
    buf[76:78] = list(sym)
    return "".join(buf).rstrip() + "\n"


def _substitute_residue(atom_lines, chain_id, res_num, ncaa_def):
    new_res_name_padded = ncaa_def.pdb_resname.ljust(3)
    backbone_atoms = {"N", "CA", "C", "O", "OXT", "H", "HA"}
    keep_atoms = set(ncaa_def.keep_atoms)
    extension_atoms = getattr(ncaa_def, "extension_atoms", ()) or ()
    parent_residue = getattr(ncaa_def, "parent_residue", None)

    result, substituted, in_target_res, cm_inserted = [], False, False, False
    n_coord, ca_coord, last_serial = None, None, 0
    # [v0.6.6] Strategy A: source residue 에서 extension anchor + neighbor 좌표 수집
    source_heavy_coords = {}   # atom_name → (x,y,z)
    source_resname = None

    for line in atom_lines:
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM"):
            try: last_serial = max(last_serial, int(line[6:11]))
            except ValueError: pass

        if rec not in ("ATOM", "HETATM"):
            result.append(line)
            continue

        try:
            line_chain = line[21]
            line_resnum = int(line[22:26].strip())
        except (ValueError, IndexError):
            result.append(line)
            continue

        atom_name = line[12:16].strip()
        is_target = (line_chain == chain_id and line_resnum == res_num)

        if in_target_res and not is_target and not cm_inserted:
            if ncaa_def.mutation_type == "N-M" and n_coord and ca_coord:
                last_serial += 1
                result.append(_make_cm_atom_line(last_serial, chain_id, res_num, new_res_name_padded, n_coord, ca_coord))
            cm_inserted = True
            in_target_res = False

        if is_target:
            in_target_res, substituted = True, True
            if source_resname is None:
                source_resname = line[17:20].strip()
            try:
                coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                if atom_name == "N": n_coord = coord
                elif atom_name == "CA": ca_coord = coord
                # [v0.6.6] extension_atoms 처리용 source heavy coords 캐시
                if extension_atoms and line[76:78].strip() != "H":
                    source_heavy_coords[atom_name] = coord
            except ValueError: pass

            if atom_name not in keep_atoms: continue

            new_line = line[:17] + new_res_name_padded + line[20:]
            if atom_name not in backbone_atoms and new_line.startswith("ATOM  "):
                new_line = "HETATM" + new_line[6:]
            result.append(new_line)
        else:
            result.append(line)

    if in_target_res and not cm_inserted:
        if ncaa_def.mutation_type == "N-M" and n_coord and ca_coord:
            last_serial += 1
            result.append(_make_cm_atom_line(last_serial, chain_id, res_num, new_res_name_padded, n_coord, ca_coord))

    # [v45 FIX — Error C] CB 가 keep_atoms 에 포함되어 있으나 결과에 없으면 자동 생성.
    # GLY 위치에 NMA 를 치환하면 원본에 CB 가 없으므로 backbone 기하학에서 생성한다.
    if substituted and "CB" in keep_atoms and n_coord and ca_coord:
        has_cb = any(
            l[12:16].strip() == "CB"
            for l in result
            if (l.startswith("ATOM") or l.startswith("HETATM"))
            and len(l) > 26 and l[21] == chain_id
            and l[22:26].strip() == str(res_num)
        )
        if not has_cb:
            c_coord_local = None
            for l in result:
                if (l.startswith("ATOM") or l.startswith("HETATM")) \
                   and len(l) > 54 and l[21] == chain_id \
                   and l[22:26].strip() == str(res_num) \
                   and l[12:16].strip() == "C":
                    c_coord_local = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
                    break
            if c_coord_local:
                cb_xyz = _generate_cb_from_backbone(n_coord, ca_coord, c_coord_local)
                last_serial += 1
                buf = list(" " * 80)
                buf[0:6]   = list("HETATM")
                buf[6:11]  = list(f"{last_serial:5d}")
                buf[12:16] = list(" CB ")
                buf[17:20] = list(new_res_name_padded[:3])
                buf[21]    = chain_id
                buf[22:26] = list(f"{res_num:4d}")
                buf[30:38] = list(f"{cb_xyz[0]:8.3f}")
                buf[38:46] = list(f"{cb_xyz[1]:8.3f}")
                buf[46:54] = list(f"{cb_xyz[2]:8.3f}")
                buf[54:60] = list("  1.00")
                buf[60:66] = list("  0.00")
                buf[76:78] = list(" C")
                cb_line = "".join(buf).rstrip() + "\n"
                # CA 원자 바로 다음에 삽입
                for idx_l in range(len(result)):
                    rl = result[idx_l]
                    if (rl.startswith("ATOM") or rl.startswith("HETATM")) \
                       and len(rl) > 26 and rl[21] == chain_id \
                       and rl[22:26].strip() == str(res_num) \
                       and rl[12:16].strip() == "CA":
                        result.insert(idx_l + 1, cb_line)
                        break

    # [v0.6.6] Strategy A: extension_atoms 를 기하학적 배치로 추가.
    # Capece 2012 (doi:10.1021/jp2082825): 1-Me-Trp 의 CM 은 NE1 에 1.466 Å,
    # indole 평면의 ring-바깥쪽 bisector 방향. SciVal 2026-04-21 approved.
    if substituted and extension_atoms and parent_residue:
        # Parent residue 가 실제로 source 와 일치하는지 검증 (fail-fast)
        if source_resname and source_resname != parent_residue:
            raise RuntimeError(
                f"[ncAA mutate] parent_residue 불일치: ncaa_def.parent_residue="
                f"'{parent_residue}' 인데 source 는 '{source_resname}' (chain "
                f"{chain_id}, resnum {res_num}). MTR 같은 parent-bootstrap ncAA "
                f"는 source 잔기가 parent 와 동일해야 합니다."
            )
        # [v0.6.7] Support chained extensions: attach_name may refer to a
        # previously-placed extension atom (e.g., phosphate O1P attaches to
        # extension P, not a parent atom). placed_ext_coords tracks them.
        placed_ext_coords: dict = {}
        declared_ext_names = {ed[0] for ed in extension_atoms}
        for (ext_name, ext_elem, attach_name, bond_len) in extension_atoms:
            if attach_name in source_heavy_coords:
                attach_coord = source_heavy_coords[attach_name]
            elif attach_name in placed_ext_coords:
                attach_coord = placed_ext_coords[attach_name]
            else:
                raise RuntimeError(
                    f"[ncAA mutate] extension_atom '{ext_name}' 의 attach "
                    f"'{attach_name}' 가 source 잔기 또는 이전에 배치된 "
                    f"extension 에 없음. Extension 순서 (root 먼저) 확인 + "
                    f"parent '{parent_residue}' 일치 확인."
                )
            # v0.6.7: amber14 template 로부터 attach 의 heavy neighbor hint 조회
            # (chained extension 의 경우 hint 는 empty — bond extension 로만 fallback)
            if attach_name in declared_ext_names:
                hint = ()  # chained: no amber template hint
            else:
                hint = _get_neighbor_hint(parent_residue, attach_name)
            if len(hint) >= 2 and hint[0] in source_heavy_coords and hint[1] in source_heavy_coords:
                ext_coord = _place_extension_atom_sp2(
                    attach_coord,
                    source_heavy_coords[hint[0]],
                    source_heavy_coords[hint[1]],
                    bond_len,
                )
            elif len(hint) == 1 and hint[0] in source_heavy_coords:
                # Single neighbor (e.g., backbone N in internal form has only CA).
                # For N-alkyl extension (CM on backbone N), compute along the
                # trans-peptide direction: away from CA, in the backbone plane.
                # Fallback: place opposite the single neighbor.
                import numpy as _np
                a = _np.array(attach_coord, dtype=float)
                n1 = _np.array(source_heavy_coords[hint[0]], dtype=float)
                v = a - n1
                v = v / max(_np.linalg.norm(v), 1e-8)
                ext_xyz = a + v * bond_len
                ext_coord = (float(ext_xyz[0]), float(ext_xyz[1]), float(ext_xyz[2]))
            elif len(hint) == 0:
                # Chained extension (e.g., phosphate O-P-O tetrahedral).
                # Place radially outward from attach atom along a canonical
                # tetrahedral direction, staggered among siblings.
                import numpy as _np
                a = _np.array(attach_coord, dtype=float)
                # Count prior chained siblings already placed
                siblings_placed = sum(
                    1 for nm in placed_ext_coords
                    if any(ed2[0] == nm and ed2[2] == attach_name for ed2 in extension_atoms)
                )
                # Tetrahedral vertices around a central atom (unit vectors)
                tet_dirs = [
                    (1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1),
                ]
                dx, dy, dz = tet_dirs[siblings_placed % 4]
                v = _np.array([dx, dy, dz], dtype=float)
                v = v / _np.linalg.norm(v)
                ext_xyz = a + v * bond_len
                ext_coord = (float(ext_xyz[0]), float(ext_xyz[1]), float(ext_xyz[2]))
            else:
                raise RuntimeError(
                    f"[ncAA mutate] parent='{parent_residue}' attach="
                    f"'{attach_name}' 에 대한 neighbor hint 획득 실패. "
                    f"amber14 ff14SB.xml 에 parent 또는 attach atom 이 "
                    f"정의되어 있는지 확인."
                )
            placed_ext_coords[ext_name] = ext_coord
            last_serial += 1
            ext_line = _make_extension_atom_line(
                last_serial, chain_id, res_num, new_res_name_padded,
                ext_name, ext_elem, ext_coord,
            )
            # attach atom 바로 다음에 삽입 (readability)
            inserted = False
            for idx_l in range(len(result)):
                rl = result[idx_l]
                if (rl.startswith("ATOM") or rl.startswith("HETATM")) \
                   and len(rl) > 26 and rl[21] == chain_id \
                   and rl[22:26].strip() == str(res_num) \
                   and rl[12:16].strip() == attach_name:
                    result.insert(idx_l + 1, ext_line)
                    inserted = True
                    break
            if not inserted:
                result.append(ext_line)

    return result, substituted

def _emit_conect_records_for_mutation(atom_lines, chain_id, res_num, ncaa_def):
    """[v0.6.6] Parent-residue 기반 ncAA 의 heavy-atom bonds 를 CONECT records 로 변환.

    run_restrained_md.py 의 graph_policy=strict 는 거리 기반 edge 추정을 금지하므로
    non-standard 잔기 (MTR, SEP, TPO 등) 는 PDB 에 explicit CONECT 를 가져야 한다.

    Args:
        atom_lines: 치환 완료된 ATOM/HETATM 라인 리스트
        chain_id, res_num: 대상 잔기
        ncaa_def: NCAADef (parent_residue 와 extension_atoms 참조)

    Returns:
        list[str]: CONECT 라인 리스트. 빈 리스트면 parent_residue 미지정.
    """
    parent = getattr(ncaa_def, "parent_residue", None)
    if not parent:
        return []
    bonds = _get_parent_bonds(parent)
    if not bonds:
        return []

    # atom name → serial 매핑 (target 잔기에 한정)
    name_to_serial = {}
    for line in atom_lines:
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
            continue
        if line[21] != chain_id or line[22:26].strip() != str(res_num):
            continue
        try:
            serial = int(line[6:11])
        except ValueError:
            continue
        name_to_serial[line[12:16].strip()] = serial

    # [v0.6.7] External peptide bonds: mutated residue is mid-chain, so the
    # backbone C must bond to next_residue.N and backbone N must bond to
    # prev_residue.C. Without these CONECT records, OpenMM templates for the
    # adjacent standard residues (e.g., VAL3 flanking MTR4) have mismatched
    # bond patterns and fail template matching.
    prev_c_serial = None
    next_n_serial = None
    for line in atom_lines:
        if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
            continue
        if line[21] != chain_id:
            continue
        try:
            line_resnum = int(line[22:26].strip())
            line_serial = int(line[6:11])
        except ValueError:
            continue
        atom_nm = line[12:16].strip()
        if line_resnum == res_num - 1 and atom_nm == "C":
            prev_c_serial = line_serial
        elif line_resnum == res_num + 1 and atom_nm == "N":
            next_n_serial = line_serial

    bonds = list(bonds)
    # Extension atoms 의 attach bond 추가 (예: NE1-CM for MTR)
    for (ext_name, _ext_elem, attach_name, _bond_len) in (getattr(ncaa_def, "extension_atoms", ()) or ()):
        bonds.append((attach_name, ext_name))

    conect_lines = []
    emitted_pairs = set()
    for a1, a2 in bonds:
        s1 = name_to_serial.get(a1)
        s2 = name_to_serial.get(a2)
        if s1 is None or s2 is None:
            continue
        key = frozenset((s1, s2))
        if key in emitted_pairs:
            continue
        emitted_pairs.add(key)
        # CONECT 표준: CONECT <serial1> <serial2>  (columns 7-11, 12-16)
        conect_lines.append(f"CONECT{s1:5d}{s2:5d}\n")
    # v0.6.7: emit external peptide-bond CONECT
    mutated_n = name_to_serial.get("N")
    mutated_c = name_to_serial.get("C")
    if prev_c_serial is not None and mutated_n is not None:
        key = frozenset((prev_c_serial, mutated_n))
        if key not in emitted_pairs:
            emitted_pairs.add(key)
            conect_lines.append(f"CONECT{prev_c_serial:5d}{mutated_n:5d}\n")
    if mutated_c is not None and next_n_serial is not None:
        key = frozenset((mutated_c, next_n_serial))
        if key not in emitted_pairs:
            emitted_pairs.add(key)
            conect_lines.append(f"CONECT{mutated_c:5d}{next_n_serial:5d}\n")
    return conect_lines


def mutate_pdb(pdb_path, targets, ncaa_def, output_path):
    atom_lines = []
    other_lines = []
    # [Refactor] PDB I/O는 한글 path 및 비-ASCII REMARK 안전성을 위해 utf-8 강제.
    with open(pdb_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_lines.append(line)
            else:
                other_lines.append(line)

    total_substituted = 0
    conect_lines_all = []
    for chain_id, res_num in targets:
        atom_lines, ok = _substitute_residue(atom_lines, chain_id, res_num, ncaa_def)
        if ok:
            total_substituted += 1
            # [v0.6.6] Strategy A: parent_residue 기반 non-standard 잔기는 CONECT 명시
            conect_lines_all.extend(
                _emit_conect_records_for_mutation(atom_lines, chain_id, res_num, ncaa_def)
            )

    with open(output_path, "w", encoding="utf-8") as f:
        for line in atom_lines:
            f.write(line)
        for line in conect_lines_all:
            f.write(line)
        if not any(l.startswith("END") for l in other_lines):
            f.write("END\n")
        else:
            for line in other_lines:
                if line.startswith("END") or line.startswith("REMARK"):
                    f.write(line)

    return total_substituted

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputdir", required=True)
    parser.add_argument("--outputdir", required=True)
    parser.add_argument("--ncaa_code", required=True, help="Registry의 Canonical Code")
    parser.add_argument("--residues", required=True)
    parser.add_argument("--binder_chain", default="B")
    parser.add_argument("--target_chain", default="A")
    parser.add_argument("--plddt_cutoff", type=float, default=70.0)
    parser.add_argument("--bsa_exclusion", type=float, default=20.0)
    args = parser.parse_args()

    os.makedirs(args.outputdir, exist_ok=True)
    
    ncaa_def = resolve_ncaa_definition(args.ncaa_code)
    if ncaa_def is None:
        log.error(f"[!] 알 수 없는 ncAA Code: {args.ncaa_code}")
        sys.exit(1)

    try:
        targets = parse_residue_targets(args.residues, input_dir=args.inputdir, binder_chain=args.binder_chain, target_chain=args.target_chain, plddt_cutoff=args.plddt_cutoff, bsa_exclusion=args.bsa_exclusion)
    except Exception as e:
        log.error(str(e))
        sys.exit(1)
        
    if not targets:
        log.error("[!] 치환할 대상이 없습니다.")
        sys.exit(1)

    log.info(f"최종 치환 위치: {[(c+str(r)) for c, r in targets]}")

    pdb_files = sorted(glob.glob(os.path.join(args.inputdir, "*.pdb")))
    success, failed = 0, 0

    for pdb_path in pdb_files:
        stem = os.path.splitext(os.path.basename(pdb_path))[0]
        output_path = os.path.join(args.outputdir, f"{stem}_ncaa.pdb")
        try:
            n_sub = mutate_pdb(pdb_path, targets, ncaa_def, output_path)
            if n_sub > 0: success += 1
            else:
                shutil.copy(pdb_path, output_path)
                failed += 1
        except Exception as e:
            log.error(f"[실패] {stem}: {e}")
            failed += 1

    write_manifest(args.outputdir, args, targets)
    # [Reporting] success/failed 카운터를 명시적으로 출력하여 배치 단위 추적성 확보.
    log.info(f"[Mutation Summary] success={success}, failed={failed}, total={success + failed}")

if __name__ == "__main__":
    main()
