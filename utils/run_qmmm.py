#!/usr/bin/env python
"""
run_qmmm.py
-----------
PySCF + OpenMM을 이용해 QM/MM 에너지 계산을 수행합니다.
- QM 영역: ncAA 잔기 + 인접 5Å 이내 잔기
- MM 영역: 나머지 단백질 + 용매
- 수준: ωB97X-D3 / 6-31G* (기본값)
사용 환경: qmmm (conda)
"""

import os
import sys
import glob
import argparse
import json
import re
import numpy as np

from pyscf import gto, dft, lib
from pyscf.qmmm import itrf as qmmm_itrf

# [FIX5] PySCF OpenMP 스레드 수 명시 설정
lib.num_threads(os.cpu_count() or 8)

from openmm import app, unit
import openmm as mm
from openmm.app import ForceField, PDBFile

# [v3 3-3, 3-4] 도메인 상수와 PDB 파서 헬퍼를 utils_common(SSoT) 에서 임포트.
# CLI 모드(`python utils/run_qmmm.py`)에서도 동작하도록 utils 경로를 sys.path 에 등록한다.
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from utils_common import KNOWN_COFACTORS, parse_pdb_atom_line  # noqa: E402


# ==========================================
# 상수 및 원소 정보
# ==========================================
ELEMENT_CHARGE = {
    "H": 1,  "C": 6,  "N": 7,  "O": 8,  "S": 16,
    "F": 9,  "P": 15, "Cl": 17,"Br": 35,"Si": 14,
}

STD_AA = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS","MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL","HID","HIE","HIP","ACE","NME"
}


# ==========================================
# PDB 파싱 — utils_common.parse_pdb_atom_line 위임
# ==========================================
def parse_pdb_atoms(pdb_path):
    """PDB 파일을 utils_common.parse_pdb_atom_line 으로 파싱하여 원자 dict 리스트를 반환한다.

    QM/MM 단계 고유의 후처리:
        - 물 (HOH/WAT/TIP3) 잔기는 즉시 제외한다.
        - element 컬럼이 비어 있거나 atom name 이 숫자로 시작 (예: 1HB) 하는
          경우, 정규식으로 숫자 접두를 제거하고 첫 알파벳을 대문자화한 결과를
          element 로 폴백한다 (탄소 'C' 가 최종 폴백).
        - is_ncaa 플래그를 부착: HETATM 이면서 KNOWN_COFACTORS 에 포함되지 않은
          잔기는 ncAA 후보로 분류된다.
    """
    atoms = []
    with open(pdb_path, encoding="utf-8") as f:
        for line in f:
            base = parse_pdb_atom_line(line)
            if base is None:
                continue
            resname = base["resname"]
            if resname in ("HOH", "WAT", "TIP3"):
                continue

            # element 폴백: 비어있거나 atom name 첫 글자가 숫자인 경우
            element = base["element"]
            if not element:
                element = re.sub(r"^[0-9]+", "", line[12:14].strip())[:1].upper() or "C"

            atoms.append({
                "serial":  base["serial"],
                "name":    base["name"],
                "resname": resname,
                "chain":   base["chain"],
                "resnum":  base["resnum"],
                "x":       base["x"],
                "y":       base["y"],
                "z":       base["z"],
                "element": element,
                "is_ncaa": (
                    base["record"] == "HETATM"
                    and resname not in ("HOH", "WAT", "TIP3")
                    and resname not in KNOWN_COFACTORS
                ),
            })
    return atoms


# ==========================================
# QM/MM 영역 분할
# ==========================================
def _partition_by_interface_distance(atoms, cutoff, binder_chain):
    """
    [Refactor] fast/plaid 모드에서 중복되던 인터페이스 잔기 추출 로직을
    공통 헬퍼로 추출한 함수.

    1. binder_chain과 비-binder 사이의 거리(원자 간 최소 거리)가 ``cutoff``
       이하인 binder 잔기를 1차 추출한다.
    2. 1차로 추출된 binder 잔기 좌표 군집과 ``cutoff`` 이하 거리에 있는
       비-binder 잔기를 2차로 추출한다.
    3. 두 단계에서 누적된 (chain, resnum) 집합을 QM 영역으로 분류한다.

    fast(4.0 Å)와 plaid(3.0 Å) 모드의 임계값 차이는 호출부에서 cutoff 인자로
    전달하므로, 본 헬퍼는 정책 결정에 관여하지 않는다 (Principle 1: 기존
    임계값 보존).
    """
    coords_b = np.array([[a["x"], a["y"], a["z"]] for a in atoms if a["chain"] == binder_chain])
    coords_a = np.array([[a["x"], a["y"], a["z"]] for a in atoms if a["chain"] != binder_chain])

    if len(coords_b) == 0 or len(coords_a) == 0:
        qm_atoms = [a for a in atoms if a["chain"] == binder_chain]
        mm_atoms = [a for a in atoms if a["chain"] != binder_chain]
        return qm_atoms, mm_atoms

    qm_resnums = set()
    for a in atoms:
        if a["chain"] == binder_chain:
            dist = np.min(np.linalg.norm(coords_a - np.array([a["x"], a["y"], a["z"]]), axis=1))
            if dist <= cutoff:
                qm_resnums.add((a["chain"], a["resnum"]))

    b_core_coords = np.array([[a["x"], a["y"], a["z"]] for a in atoms if (a["chain"], a["resnum"]) in qm_resnums])
    if len(b_core_coords) > 0:
        for a in atoms:
            if a["chain"] != binder_chain:
                dist = np.min(np.linalg.norm(b_core_coords - np.array([a["x"], a["y"], a["z"]]), axis=1))
                if dist <= cutoff:
                    qm_resnums.add((a["chain"], a["resnum"]))

    qm_atoms = [a for a in atoms if (a["chain"], a["resnum"]) in qm_resnums]
    mm_atoms = [a for a in atoms if (a["chain"], a["resnum"]) not in qm_resnums]
    return qm_atoms, mm_atoms


# ==========================================
# [v4 S-2] QM/MM Link Atom (H-cap) 경계 보정
# ==========================================
def _add_link_atoms(qm_atoms, mm_atoms, all_atoms, cn_max_distance=1.7):
    """[v4 S-2] QM/MM 경계에서 절단된 펩타이드 결합 (C-N) 에 H-cap link atom 을 삽입한다.

    배경:
        잔기 단위로 QM 영역을 분할하면 backbone 의 C(i) - N(i+1) 펩타이드 결합이
        QM-MM 경계를 가로지르는 경우가 발생한다. 이때 QM 측 원자 (예: C) 는
        본래 결합 partner (N) 가 사라져 *공중 부유* 상태가 되어, RKS/DFT 가
        partial valence 구조로 발산하거나 spurious dipole moment 가 생긴다.

    해결 (Senn & Thiel, Angew. Chem. Int. Ed. 2009):
        절단된 결합의 QM 측 원자에서 1.09 Å (C-H 평균 결합 길이) 떨어진 위치에
        link H 원자를 삽입하여 QM 영역의 valence 를 satisfy 한다. link H 의
        좌표는 원래 결합 방향을 따라 배치되어 기하 왜곡을 최소화한다.

    Args:
        qm_atoms:  partition_qmmm 가 산출한 QM 원자 dict 리스트
        mm_atoms:  partition_qmmm 가 산출한 MM 원자 dict 리스트
        all_atoms: 전체 원자 (parse_pdb_atoms 결과)
        cn_max_distance: 펩타이드 결합으로 인정할 C-N 거리 임계값 (Å). 기본 1.7 Å.

    Returns:
        list[dict]: 추가될 link H 원자 dict 리스트. build_qm_mol 가 본 리스트를
        qm_atoms 에 합쳐 PySCF 분자 정의에 포함시킨다.
    """
    if not qm_atoms or not mm_atoms:
        return []

    # QM/MM 멤버십을 (chain, resnum, name) 키로 빠른 조회
    qm_keys = {(a["chain"], a["resnum"], a["name"]) for a in qm_atoms}
    mm_keys = {(a["chain"], a["resnum"], a["name"]) for a in mm_atoms}

    link_atoms = []
    BOND_LEN_CH = 1.09  # C-H 평균 결합 길이 (Å)

    # 펩타이드 결합 후보: backbone C 와 그 C 가 속한 잔기 i 의 다음 잔기 i+1 의 N
    # PDB 의 동일 chain 에서 resnum 이 1 차이나는 모든 (C, N) 쌍을 검토
    by_chain_res = {}
    for a in all_atoms:
        by_chain_res.setdefault((a["chain"], a["resnum"]), []).append(a)

    for (chain, resnum), atoms_in_res in by_chain_res.items():
        c_atom = next((a for a in atoms_in_res if a["name"] == "C"), None)
        if c_atom is None:
            continue
        next_res_atoms = by_chain_res.get((chain, resnum + 1))
        if not next_res_atoms:
            continue
        n_next = next((a for a in next_res_atoms if a["name"] == "N"), None)
        if n_next is None:
            continue

        # C-N 거리 검증 (펩타이드 결합으로 인정 가능한 거리?)
        d = np.linalg.norm(np.array([c_atom["x"], c_atom["y"], c_atom["z"]])
                           - np.array([n_next["x"], n_next["y"], n_next["z"]]))
        if d > cn_max_distance:
            continue

        c_key = (c_atom["chain"], c_atom["resnum"], "C")
        n_key = (n_next["chain"], n_next["resnum"], "N")

        # 경계 케이스: C ∈ QM, N ∈ MM → QM 측 (C) 에 link H 추가, 방향은 C→N
        if c_key in qm_keys and n_key in mm_keys:
            link = _make_link_h(c_atom, n_next, BOND_LEN_CH)
            link_atoms.append(link)

        # 반대 경계: N ∈ QM, C ∈ MM → QM 측 (N) 에 link H 추가, 방향은 N→C
        if n_key in qm_keys and c_key in mm_keys:
            link = _make_link_h(n_next, c_atom, BOND_LEN_CH)
            link_atoms.append(link)

    if link_atoms:
        print(f"  [Link Atom] QM/MM 경계 절단 결합 {len(link_atoms)}개 → H-cap 삽입 (Senn & Thiel 2009)")
    return link_atoms


def _make_link_h(qm_atom, mm_atom, bond_len):
    """qm_atom 위치에서 mm_atom 방향으로 bond_len Å 거리에 H link atom dict 를 생성한다."""
    p_qm = np.array([qm_atom["x"], qm_atom["y"], qm_atom["z"]])
    p_mm = np.array([mm_atom["x"], mm_atom["y"], mm_atom["z"]])
    v = p_mm - p_qm
    n = np.linalg.norm(v)
    if n < 1e-6:
        v = np.array([1.0, 0.0, 0.0])
        n = 1.0
    v /= n
    h_pos = p_qm + v * bond_len
    return {
        "serial":  -1,                      # link atom 은 PDB 일련번호 미부여
        "name":    "HL",                    # Hydrogen Link
        "resname": qm_atom["resname"],
        "chain":   qm_atom["chain"],
        "resnum":  qm_atom["resnum"],
        "x": float(h_pos[0]),
        "y": float(h_pos[1]),
        "z": float(h_pos[2]),
        "element": "H",
        "is_ncaa": False,
        "is_link": True,                    # build_qm_mol 진단용 플래그
    }


def partition_qmmm(atoms, qm_cutoff=5.0, mode="full", binder_chain="B"):
    """
    mode="full": 바인더 전체를 QM으로.
    mode="fast": 타겟과 바인더가 맞닿은 4.0 Å 이내 핵심 인터페이스 잔기만 추출하여 QM으로.
    mode="plaid" : 타겟과 바인더가 맞닿은 3.0 Å 이내 핵심 인터페이스 잔기만 추출하여 QM으로.
    """
    ncaa_atoms = [a for a in atoms if a["is_ncaa"]]

    # [Plaid 초고속 다이어트 모드] — 3.0 Å 임계값
    if mode == "plaid":
        # [v0.3 #4] cutoff 값을 단일 변수로 관리하여 print/연산 일치 보장
        _cutoff = 3.0
        print(f"  [Calculate Mode: PLAID] {_cutoff} Å 이내 인터페이스 핵심 잔기만 QM 영역으로 추출합니다.")
        qm_atoms, mm_atoms = _partition_by_interface_distance(atoms, cutoff=_cutoff, binder_chain=binder_chain)
        print(f"  [PLAID] 추출된 핵심 QM 원자 수: {len(qm_atoms)}개")
        if len(qm_atoms) == 0:
            print("  [!] 경고: 펩타이드가 단백질에서 90Å 이상 멀어졌습니다! (결합 붕괴)")
        return qm_atoms, mm_atoms

    # [Fast 다이어트 모드]
    if mode == "fast":
        # [v0.3 #4] cutoff 값을 단일 변수로 관리하여 print/연산 일치 보장
        _cutoff = 4.0
        print(f"  [Calculate Mode: FAST] {_cutoff} Å 이내 인터페이스 핵심 잔기만 QM 영역으로 추출합니다.")
        qm_atoms, mm_atoms = _partition_by_interface_distance(atoms, cutoff=_cutoff, binder_chain=binder_chain)
        print(f"  [FAST] 추출된 핵심 QM 원자 수: {len(qm_atoms)}개")
        return qm_atoms, mm_atoms

    # [Full 모드]
    if len(ncaa_atoms) == 0:
        qm_atoms = [a for a in atoms if a["chain"] == binder_chain]
        mm_atoms = [a for a in atoms if a["chain"] != binder_chain]
        print(f"  [Calculate Mode: FULL] Chain {binder_chain} 전체를 QM 영역으로 ({len(qm_atoms)}개 원자)")
        return qm_atoms, mm_atoms

    qm_resnums = set()
    ncaa_coords = np.array([[a["x"], a["y"], a["z"]] for a in ncaa_atoms])
    for a in atoms:
        if a["is_ncaa"]:
            qm_resnums.add((a["chain"], a["resnum"]))
            continue
        dist = np.min(np.linalg.norm(ncaa_coords - np.array([a["x"], a["y"], a["z"]]), axis=1))
        if dist <= qm_cutoff:
            qm_resnums.add((a["chain"], a["resnum"]))

    qm_atoms = [a for a in atoms if (a["chain"], a["resnum"]) in qm_resnums]
    mm_atoms = [a for a in atoms if (a["chain"], a["resnum"]) not in qm_resnums]
    print(f"  [Mode: FULL] QM 영역: {len(qm_atoms)}개 원자 ({len(qm_resnums)}개 잔기)")
    return qm_atoms, mm_atoms


# ==========================================
# PySCF QM 분자 빌드
# ==========================================
def build_qm_mol(qm_atoms, charge=0, spin=0, basis="6-31G*"):
    """QM 원자 목록 → PySCF gto.M 객체.

    홀수 전자 보정 정책:
        총 양성자 수가 홀수인 경우 charge 를 -1 로 조정하여 singlet (spin=0)
        조건을 유지한다.

        음이온(-1) 선택 근거:
            펩타이드 결합 환경에서 C-N 절단으로 단편화된 cluster 가 RKS
            폐각 가정과 충돌할 때, 전자 수용(음이온 형성)이 양성자 방출보다
            국소 정전기적으로 유리한 경우가 일반적이다 (백본 카르보닐 산소
            및 측쇄 카르복실 그룹의 음전하 친화도). 따라서 charge=-1 가
            대다수의 인터페이스 단편에 대해 더 안정적인 수렴 경로를 제공한다.

        주의:
            - 양이온 라디칼 환경, 또는 금속 이온 (e.g. Zn2+, Mg2+) 인접
              cluster 에서는 charge=+1 이 더 물리적으로 정확할 수 있다.
            - 본 보정은 fragment cluster 의 SCF 수렴을 위한 *수치적 보정* 이며,
              ESP/에너지 절대값이 아닌 *상대적 비교* (동일 cutoff 모드 내 design
              간 순위 매김) 용도로 결과를 해석해야 한다.
            - 결합 절단 단편의 chemical accuracy 가 필요한 경우, 사용자는 QM
              cutoff 를 늘리거나 link-atom (수소 캡) 방식으로 전환을 고려해야 한다.
    """
    atom_str = ""
    total_protons = 0  # 총 양성자 수 계산용 변수
    
    for a in qm_atoms:
        elem = a["element"] if a["element"] else "C"
        # [FIX1] H 포함: PySCF는 H를 자동 추가하지 않음
        atom_str += f"{elem} {a['x']:.6f} {a['y']:.6f} {a['z']:.6f};"
        # 원소별 양성자 수 누적
        total_protons += ELEMENT_CHARGE.get(elem.capitalize(), 0)

    if not atom_str:
        raise ValueError("QM 영역에 원자가 없습니다.")

    # 홀수 전자 보정 (정책은 함수 docstring 참조)
    if total_protons % 2 != 0:
        charge = -1
        print(f"    [!] 홀수 전자 감지 (Protons: {total_protons}) -> Charge를 -1로 자동 조정하여 Spin=0 붕괴 방어")

    # Si 등 비표준 원소는 def2-SVP로 fallback
    # [FIX1] 모든 원소 포함 (H 포함), H는 6-31G** 사용
    all_elements = set(a["element"] for a in qm_atoms if a["element"])
    basis_dict = {}
    for elem in all_elements:
        if elem in ("Si", "Br", "I"):
            basis_dict[elem] = "def2-SVP"
        elif elem == "H":
            basis_dict[elem] = "6-31G**"
        else:
            basis_dict[elem] = basis

    mol = gto.M(
        atom    = atom_str.rstrip(";"),
        basis   = basis_dict if basis_dict else basis,
        charge  = charge,
        spin    = spin,
        unit    = "Angstrom",
        verbose = 3,
    )
    return mol

# ==========================================
# MM 점전하 배열 생성
# ==========================================
# [FIX2] AMBER ff14SB 전 잔기 원자별 부분전하 테이블 (Maier et al. 2015)
# [Refactor] AMBER ff14SB 표준 잔기·원자별 부분전하 테이블과 백본 폴백을
# utils/amber_charges.py 로 분리하여 본 모듈의 가독성과 재사용성을 확보한다.
# 데이터 의미와 변경 절차는 amber_charges.py 의 모듈 docstring을 참조한다.
# [v3 4-2] _BACKBONE_FALLBACK 은 BACKBONE_FALLBACK_CHARGES 로 정식 이름 변경되었다.
from amber_charges import AMBER_FF14SB, BACKBONE_FALLBACK_CHARGES  # noqa: E402


def get_mm_charges(mm_atoms):
    """
    [FIX2] AMBER ff14SB 잔기별·원자별 전하 테이블로 MM 점전하 할당
    Returns: np.array shape (n_mm, 4) — [x, y, z, q] (Angstrom)

    [#13 v0.2] PySCF add_mm_charges 는 mol.unit(Å) 기준으로 자동 Bohr 변환을
    수행하므로 (https://pyscf.org/pyscf_api_docs/pyscf.qmmm.html) 여기에서
    Å→Bohr 수동 변환을 수행하면 이중 변환이 발생한다. 좌표를 Å 그대로 반환.
    """
    mm_charge_array = []
    n_unknown = 0
    for a in mm_atoms:
        res_q = AMBER_FF14SB.get(a["resname"], {})
        q = res_q.get(a["name"], BACKBONE_FALLBACK_CHARGES.get(a["name"], 0.0))
        if q == 0.0 and a["name"] not in BACKBONE_FALLBACK_CHARGES and not res_q:
            n_unknown += 1
        # [#13] Å 좌표 그대로 전달 — PySCF 가 mol.unit 기준 자동 변환.
        mm_charge_array.append([a["x"], a["y"], a["z"], q])
    if n_unknown:
        print(f"  [MM charges] 테이블 미등록 원자 {n_unknown}개 → 전하=0.0")
    return np.array(mm_charge_array)


# ==========================================
# QM/MM 계산 실행
# ==========================================
def run_qmmm_calc(pdb_path, output_dir, qm_basis, qm_xc, ncaa_elem, qm_cutoff=5.0, mode="full", binder_chain="B"):
    """단일 스냅샷 PDB에 대한 QM/MM 계산"""

    # [#40 v0.2] max_memory 를 환경 변수 UPDD_QM_MAX_MEMORY(MB) 로 주입 (기본 16000 = 16GB).
    _max_mem = int(os.environ.get("UPDD_QM_MAX_MEMORY", "16000"))

    basename = os.path.basename(pdb_path).replace(".pdb", "")
    out_json = os.path.join(output_dir, f"{basename}_qmmm_{mode}.json")

    if os.path.exists(out_json):
        # [v0.3 #1+#5] resume 무결성 강화:
        #   - converged=True AND energy_total_hartree 유효 → 스킵 (정상 완료)
        #   - converged=True 인데 interaction_kcal=None → QM-only 만 실패한 케이스
        #     → JSON 삭제 후 전체 재계산 (SciVal 승인: 독립 SCF 재계산이 신뢰성 보장)
        #   - 그 외 (converged=False / energy_total=None / 손상) → 삭제 후 재계산
        try:
            with open(out_json, "r", encoding="utf-8") as f:
                existing = json.load(f)
            if existing.get("converged") and existing.get("energy_total_hartree") is not None:
                # QM-only 만 실패한 경우 → 전체 재계산
                if existing.get("interaction_kcal") is None:
                    print(f"\n  [Retry] QM-only 미완료 → 전체 재계산: {basename}")
                    os.remove(out_json)
                else:
                    print(f"\n  [Skip] 이미 완료: {basename}")
                    return existing
            else:
                reason = existing.get("reason", "null_energy")
                print(f"\n  [Cleanup] 미수렴 결과 삭제 ({reason}): {basename}")
                os.remove(out_json)
        except (json.JSONDecodeError, IOError):
            print(f"\n  [Cleanup] 손상된 JSON 삭제: {basename}")
            try:
                os.remove(out_json)
            except OSError:
                pass

    print(f"\n  계산: {basename}")

    atoms = parse_pdb_atoms(pdb_path)
    if not atoms:
        print(f"  [!] 원자 없음: {pdb_path}")
        return None

    coords_b = np.array([[a["x"], a["y"], a["z"]] for a in atoms if a["chain"] == binder_chain])
    coords_a = np.array([[a["x"], a["y"], a["z"]] for a in atoms if a["chain"] != binder_chain])
    
    if len(coords_b) > 0 and len(coords_a) > 0:
        min_dist = np.min([np.min(np.linalg.norm(coords_a - b_pos, axis=1)) for b_pos in coords_b])
        if min_dist > 8.0:
            # [v39 Fix 5] PBC 래핑 가능성을 진단하기 위해 CRYST1 박스 크기를 함께 출력
            pbc_hint = ""
            try:
                with open(pdb_path, "r", encoding="utf-8", errors="ignore") as _cpf:
                    for _cline in _cpf:
                        if _cline.startswith("CRYST1"):
                            _box = float(_cline[6:15])
                            if _box > 0 and min_dist > _box * 0.8:
                                pbc_hint = (f" (PBC 박스={_box:.1f}Å — 래핑 의심, "
                                            f"extract_snapshots 의 image_molecules 확인)")
                            break
            except Exception:
                pass
            print(f"  [!] 도망자 발견: 펩타이드가 {min_dist:.1f}Å 우주로 떠나갔습니다.{pbc_hint}")
            # [#36 v0.2] 도망자는 converged=False + interaction_kcal=None 으로
            # 랭킹에서 제외. reason/min_dist_A 로 사후 분석 가능.
            return {
                "snapshot": basename, "energy_total_hartree": None, "energy_qmmm_kcal": None,
                "interaction_kcal": None, "converged": False, "is_run_mmgbsa": False,
                "reason": "binder_escaped", "min_dist_A": round(float(min_dist), 1),
            }

    qm_atoms, mm_atoms = partition_qmmm(atoms, qm_cutoff, mode=mode, binder_chain=binder_chain)

    # [#38 v0.2] QM 원자 수 상한 가드 — QM-only 계산 발산 및 GPU OOM 예방.
    MAX_QM_ATOMS = 500
    if len(qm_atoms) == 0:
        # [#36] QM 원자 0개는 converged=False 로 표시하여 랭킹 제외.
        return {
            "snapshot": basename, "energy_total_hartree": None, "energy_qmmm_kcal": None,
            "interaction_kcal": None, "converged": False, "is_run_mmgbsa": False,
            "reason": "no_qm_atoms",
        }
    if len(qm_atoms) > MAX_QM_ATOMS:
        print(f"  [!] QM 원자 {len(qm_atoms)}개 > 상한 {MAX_QM_ATOMS}")
        return {
            "snapshot": basename, "energy_total_hartree": None, "energy_qmmm_kcal": None,
            "interaction_kcal": None, "converged": False, "is_run_mmgbsa": False,
            "reason": "qm_atoms_exceeded", "n_qm_atoms": len(qm_atoms),
        }

    # [v4 S-2] QM/MM 경계 H-cap link atom 추가 (Senn & Thiel 2009)
    link_atoms = _add_link_atoms(qm_atoms, mm_atoms, atoms)
    if link_atoms:
        qm_atoms = qm_atoms + link_atoms

    try:
        mol = build_qm_mol(qm_atoms, basis=qm_basis)
    except Exception as e:
        print(f"  [!] QM 분자 빌드 실패: {e}")
        return None

    print(f"  QM 기저함수: {qm_basis} | 범함수: {qm_xc}")
    print(f"  QM 원자 수 (중원자): {mol.natm}")

    mm_coords_charges = get_mm_charges(mm_atoms)

    # ── 4. DFT 계산 (QM/MM 임베딩) ──────────────────────────
    mf = dft.RKS(mol)
    mf.xc = qm_xc
    # [v0.3.3] verbose=4: v0.1 롤백
    mf.verbose = 4
    mf.max_memory = _max_mem
    mf.grids.level = 3
    mf.max_cycle   = 300
    mf.conv_tol    = 1e-8

    gpu_success = False
    try:
        # 1. 먼저 평범한 계산기를 GPU 로켓으로 변신!
        import gpu4pyscf
        from gpu4pyscf.qmmm import itrf as gpu_qmmm_itrf
        mf = mf.to_gpu()
        
        # 2. 변신이 성공하면, GPU 전용 짐칸에 점전하 탑재!
        if len(mm_coords_charges) > 0:
            mf = gpu_qmmm_itrf.add_mm_charges(mf, mm_coords_charges[:, :3], mm_coords_charges[:, 3])
            
        gpu_success = True
        print(f"  [+] MM 점전하: {len(mm_coords_charges)}개 포함 (GPU 상태)")
        print("  🚀 [GPU 가속 - DFT계산]")
        
    except Exception as e:
        print(f"  [i] GPU 가속 실패 (사유: {e}). CPU 멀티코어로 우회합니다.")
        # GPU 변신이 실패하면 CPU 장갑차(QMMMRKS)로 우회 탑재
        mf = dft.RKS(mol)
        mf.xc = qm_xc
        # [v0.3.3] verbose=4: v0.1 롤백
        mf.verbose = 4
        mf.max_memory = _max_mem
        mf.grids.level = 3
        mf.max_cycle   = 300
        mf.conv_tol    = 1e-8
        if len(mm_coords_charges) > 0:
            mf = qmmm_itrf.mm_charge(mf, mm_coords_charges[:, :3], mm_coords_charges[:, 3])
            print(f"  [+] MM 점전하: {len(mm_coords_charges)}개 포함 (CPU 상태)")

    try:
        energy_total = mf.kernel()
        print(f"  QM/MM 에너지: {energy_total:.8f} Hartree")
    except Exception as e:
        print(f"  [!] DFT 수렴 실패: {e}")
        energy_total = None

    # ── 5. QM 단독 에너지 계산 (상호작용 에너지 분해) ──────────────────────────────
    e_qm_only = None
    if energy_total is not None and len(mm_coords_charges) > 0:
        mf_qm = dft.RKS(mol)
        mf_qm.xc = qm_xc
        # [v0.3.3] verbose=4: v0.1 롤백
        mf_qm.verbose = 4
        mf_qm.max_cycle = 300
        mf_qm.conv_tol = 1e-8  
        mf_qm.max_memory = _max_mem
                
        # 2라운드는 점전하 짐칸이 없으므로 그냥 변신만 시키면 됩니다!
        if gpu_success:
            try:
                mf_qm = mf_qm.to_gpu()
                print("  🚀 [GPU 가속 - QM 단독계산]")
            except Exception as e:
                print(f"  [!] 2라운드 GPU 전환 실패 (사유: {e})")

        mf_qm.grids.level = 3
        try:
            e_qm_only  = mf_qm.kernel()
            e_interact = (energy_total - e_qm_only) * 627.509
            print(f"  QM 단독 에너지: {e_qm_only:.8f} Hartree")
            print(f"  QM-MM 상호작용: {e_interact:.4f} kcal/mol")
        except Exception as e_qm_err:
            # [v0.3 #5] QM-only 실패 진단 — QM/MM E_total 은 유효하므로 converged=True 유지.
            # interaction_kcal=None 으로 기록되며, 다음 resume 시 QM-only 재시도 분기로
            # 진입한다 (run_qmmm_calc 의 resume 로직 #1+#5 참조).
            print(f"  [!] QM-only 계산 실패: {e_qm_err}")
            e_interact = None
    else:
        e_interact = None

    # ── 6. 결과 저장 ─────────────────────────────────────────
    result = {
        "snapshot":       basename,
        "pdb_path":       pdb_path,
        "qm_method":      f"{qm_xc}/{qm_basis}",
        "n_qm_atoms":     mol.natm,
        "n_mm_atoms":     len(mm_atoms),
        "ncaa_element":   ncaa_elem,
        "energy_total_hartree":   energy_total,
        "energy_qm_hartree":      e_qm_only,
        "energy_qmmm_kcal":       float(energy_total * 627.509) if energy_total else None,
        "interaction_kcal":       e_interact,
        "converged":              energy_total is not None,
    }

    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    print(f"  결과 저장: {out_json}")
    return result


# ==========================================
# CLI 메인
# ==========================================
def main():
    parser = argparse.ArgumentParser(
        description="QM/MM — PySCF ωB97X-D3 + OpenMM 점전하 임베딩"
    )
    parser.add_argument("--snapdir",   required=True,         help="스냅샷 PDB 디렉토리")
    parser.add_argument("--outputdir", required=True,         help="QM/MM 결과 출력 디렉토리")
    parser.add_argument("--qm_basis",  default="6-31G*",      help="QM 기저함수 (기본: 6-31G*)")
    parser.add_argument("--qm_xc",     default="wb97x-d3",    help="DFT 범함수 (기본: wb97x-d3)")
    parser.add_argument("--ncaa_elem", default="none",        help="ncAA 핵심 원소 (예: Si)")
    parser.add_argument("--cutoff",    type=float, default=5.0, help="QM 영역 cutoff Å (기본 5.0)")
    parser.add_argument("--mode", default="fast", choices=["full", "fast", "plaid"], help="QM/MM 계산 모드")
    parser.add_argument("--filter", default=None, help="특정 디자인 ID만 골라서 계산 (예: design_w2_2_s1)")
    parser.add_argument("--binder_chain", default="B", help="Binder chain ID (기본 B)")
    args = parser.parse_args()

    os.makedirs(args.outputdir, exist_ok=True)

    all_pdbs = sorted(glob.glob(os.path.join(args.snapdir, "*.pdb")))
    if args.filter:
        pdb_files = [f for f in all_pdbs if os.path.basename(f).startswith(args.filter)]
        print(f"  [Filter] '{args.filter}'로 시작하는 스냅샷 {len(pdb_files)}개를 찾았습니다.")
    else:
        pdb_files = all_pdbs
    if not pdb_files:
        print(f"[!] PDB 파일 없음: {args.snapdir}")
        return

    print(f"\n[QM/MM Calculation]")
    print(f"  스냅샷 수 : {len(pdb_files)}개")
    print(f"  수준      : {args.qm_xc}/{args.qm_basis}")
    print(f"  QM cutoff : {args.cutoff} Å")
    print(f"  ncAA 원소 : {args.ncaa_elem}")

    all_results = []
    failed      = []

    for pdb_path in pdb_files:
        result = run_qmmm_calc(
            pdb_path   = pdb_path,
            output_dir = args.outputdir,
            qm_basis   = args.qm_basis,
            qm_xc      = args.qm_xc,
            ncaa_elem  = args.ncaa_elem,
            qm_cutoff  = args.cutoff,
            mode       = args.mode, 
            binder_chain = args.binder_chain
        )
        if result:
            all_results.append(result)
        else:
            failed.append(os.path.basename(pdb_path))

    # [v55 FIX] --filter 모드(병렬 워커)에서는 전체 summary를 작성하지 않는다.
    # 여러 워커가 동시에 summary를 덮어쓰면 마지막 워커의 결과만 남는다.
    # 통합 summary는 orchestrator(UPDD.py execute_qmmm)에서 생성한다.
    if args.filter:
        summary_path = "(개별 JSON만 저장 — 통합 summary는 orchestrator에서 생성)"
    else:
        summary_path = os.path.join(args.outputdir, "qmmm_summary.json")
        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump({
                "method":    f"{args.qm_xc}/{args.qm_basis}",
                "n_calc":    len(all_results),
                "n_failed":  len(failed),
                "failed":    failed,
                "results":   all_results,
            }, f, indent=2)

    # 수렴 결과 에너지 순 정렬 출력
    converged = [r for r in all_results if r["converged"]]
    converged.sort(key=lambda r: r["energy_total_hartree"])

    print(f"\n  완료: {len(converged)}/{len(pdb_files)} 수렴")
    if converged:
        print(f"\n  [에너지 순위 (낮을수록 안정)]")
        for i, r in enumerate(converged, 1):
            e_kcal = r["energy_qmmm_kcal"]
            ia     = r["interaction_kcal"]
            ia_str = f"{ia:+.2f} kcal/mol" if ia is not None else "N/A"
            print(f"  {i:2}. {r['snapshot'][:40]:<40}  "
                  f"E={e_kcal:>14.4f} kcal/mol  QM-MM={ia_str}")

    print(f"\n  요약 저장: {summary_path}")


if __name__ == "__main__":
    main()
