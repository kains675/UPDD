#!/usr/bin/env python
"""
utils/utils_common.py
---------------------
UPDD 파이프라인의 utils 패키지 전반에서 중복 구현되던 보조 함수와 도메인
상수를 통합한 공통 모듈. 신규 항목을 추가할 때는 반드시 다음 두 원칙을 준수한다.

1. 본 모듈의 함수는 사이드 이펙트(파일 I/O, 로깅, OS 호출)를 최소화하고
   수학적/위상수학적 변환에만 집중한다.
2. 기존 호출 모듈의 로컬 함수는 삭제하지 않고 본 모듈의 함수를 호출하는
   wrapper 형태로 유지한다 (CLAUDE.md 원칙 1: 로직 보존).

[v1] resolve_chainid_by_letter 통합
- graft_ligand_generic._resolve_chainid_by_letter
- extract_snapshots._resolve_binder_chainid
두 곳에서 동일한 알고리즘으로 mdtraj Topology의 chain letter → chain index
를 해석하던 로직을 본 모듈의 resolve_chainid_by_letter()로 통합하였다.

[v3] 도메인 상수 및 PDB 파싱 헬퍼 통합
- run_qmmm.KNOWN_COFACTORS 와 preprocess_target.COMMON_IONS / COMMON_COFACTORS
  를 단일 진실 공급원(SSoT)으로 통합한다 (CLAUDE.md 3-3).
- run_qmmm / preprocess_target / ncaa_mutate 의 인라인 PDB 컬럼 파서를
  parse_pdb_atom_line() 으로 통합한다 (CLAUDE.md 3-4).
- preprocess_target.dist 함수를 atom_distance() 로 통합한다 (CLAUDE.md 12).
"""

from typing import Optional, Any, Dict, Tuple, List
import math


# ==========================================
# 도메인 상수 — 단일 진실 공급원 (Single Source of Truth)
# ==========================================
# [3-3] PDB 보조인자/이온 명단. KNOWN_COFACTORS는 QM/MM, MM-GBSA, preprocess
# 단계 모두에서 "리간드가 아닌 보조분자"로 분류해야 하는 잔기 식별에 사용된다.
# 본 상수를 수정할 때는 반드시 모든 호출부의 정합성을 함께 확인한다.
KNOWN_COFACTORS = {
    "GDP", "GTP", "ATP", "ADP", "AMP", "GMP", "CTP", "UTP",
    "NAD", "NADH", "NAP", "NADP", "FAD", "FMN", "HEM", "HEC",
    "MG", "ZN", "CA", "FE", "MN", "CU", "CO", "NI", "FE2", "FE3",
    "SO4", "PO4", "GOL", "EDO", "PEG", "MPD", "ACT", "IMD",
    "CL", "NA", "K", "BR", "IOD",
}

COMMON_IONS = {
    "NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE", "CU", "CO", "NI",
}

# preprocess_target 의 COMMON_COFACTORS 와 KNOWN_COFACTORS 의 합집합과 동치인
# 호환용 상수. preprocess_target 호출 규약 보존을 위해 별도로 노출한다.
COMMON_COFACTORS = {
    "GDP", "GTP", "ADP", "ATP", "SAM", "SAH", "FAD", "FMN", "NAD", "NAP",
    "HEM", "HEME", "PLP", "PQN", "COA",
}

WATER_RESNAMES = {"HOH", "WAT", "TIP", "TIP3", "SOL"}


# ==========================================
# 기하 헬퍼
# ==========================================
def atom_distance(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    """3차원 공간 두 점(angstrom 단위)의 유클리드 거리를 반환한다.

    preprocess_target.dist 와 정확히 동일한 수치 계약을 유지한다 (numpy 의존성
    없이 stdlib math 만 사용).
    """
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


# ==========================================
# PDB 파싱 — 컬럼 기반 표준 파서
# ==========================================
def parse_pdb_atom_line(line: str) -> Optional[Dict[str, Any]]:
    """PDB ATOM/HETATM 레코드 한 줄을 파싱하여 표준 dict 를 반환한다.

    PDB v3 컬럼 규약에 따라 고정 폭으로 파싱한다. 파싱이 불가능한(필드 누락,
    숫자 변환 실패 등) 경우 None을 반환하며, 호출부는 None 반환을 "이 줄은
    원자 레코드가 아니다" 라는 정상 신호로 처리해야 한다.

    Returns:
        Optional[Dict[str, Any]]: 다음 키를 가진 dict 또는 None
            - record  (str): "ATOM" 또는 "HETATM"
            - serial  (int): 원자 일련번호
            - name    (str): 원자 이름 (앞뒤 공백 제거)
            - resname (str): 잔기 이름 (앞뒤 공백 제거)
            - chain   (str): chain ID (공백이면 빈 문자열)
            - resnum  (int): 잔기 번호
            - x, y, z (float): 좌표 (Å)
            - element (str): 원소 기호 (col 76-78). 미지정 시 빈 문자열
    """
    rec = line[:6].strip()
    if rec not in ("ATOM", "HETATM"):
        return None
    try:
        return {
            "record":  rec,
            "serial":  int(line[6:11]),
            "name":    line[12:16].strip(),
            "resname": line[17:20].strip(),
            "chain":   line[21].strip(),
            "resnum":  int(line[22:26]),
            "x":       float(line[30:38]),
            "y":       float(line[38:46]),
            "z":       float(line[46:54]),
            "element": line[76:78].strip() if len(line) > 76 and line[76:78].strip() else "",
        }
    except (ValueError, IndexError):
        return None


def resolve_chainid_by_letter(topology: Any, chain_letter: Optional[str]) -> Optional[int]:
    """mdtraj Topology의 chain letter(A/B/C...)를 chain.index 정수로 해석한다.

    Args:
        topology: mdtraj.Topology 객체. ``chains`` 이터러블을 노출하며 각
            chain은 ``chain_id`` 속성과 ``index`` 속성을 가져야 한다.
        chain_letter: 해석 대상 체인 문자(대소문자 무관). None 또는 빈 문자열을
            전달하면 해석을 시도하지 않고 즉시 None을 반환한다.

    Returns:
        Optional[int]: 매칭되는 chain의 index. 매칭 실패 시 None을 반환하며,
        호출부가 폴백 경로(legacy chainid heuristic 등)로 진행하도록 위임한다.

    Notes:
        - 본 함수는 사이드 이펙트가 없으며, 빈 chain_id 컬럼을 가진 PDB 파일
          (' '으로 출력하는 일부 writer)을 안전하게 처리한다.
        - 호출부는 본 함수의 None 반환을 정상적인 폴백 신호로 처리해야 한다.
    """
    if not chain_letter:
        return None
    target = str(chain_letter).strip().upper()
    for chain in topology.chains:
        chain_letter_found = getattr(chain, "chain_id", None) or ""
        if chain_letter_found.upper() == target:
            return chain.index
    return None
