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
# [DIAG] 구조화 로그 헬퍼 — Path 에이전트 자동 진단용
# ==========================================
# Path 가 `grep "\[DIAG\]"` 한 번으로 모든 결정점을 수집할 수 있도록 JSON-line
# 형식으로 출력한다. 계산 결과는 변경하지 않으며 순수 관찰성(observability)
# 개선만 수행한다. payload 는 JSON-serializable 타입만 허용
# (None/str/int/float/bool/list/dict).
def _diag(event, **kwargs):
    """Emit a single-line structured diagnostic log.

    Format: ``[DIAG] {"event": "<event>", ...kwargs}``

    Args:
        event: 결정점 식별자 (snake_case 권장).
        **kwargs: payload 필드. 모두 JSON-serializable 이어야 한다.
    """
    payload = {"event": event}
    payload.update(kwargs)
    print("[DIAG] " + json.dumps(payload, ensure_ascii=False), flush=True)


# ==========================================
# 예외 클래스 — Policy A fail-fast
# ==========================================
class OddElectronError(ValueError):
    """Raised when the QM region has an odd number of electrons.

    Policy A (v0.6): fail-fast — no silent charge=-1 correction.
    The caller (run_qmmm_calc) should catch this, record the failure in the
    snapshot JSON, and skip this snapshot without aborting the whole design.

    Inherits from ``ValueError`` for backward compatibility with existing
    ``except ValueError`` handlers upstream.
    """
    pass


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
# Sidechain atom selector — topology mode helper
# ==========================================
_BACKBONE_NAMES = frozenset({
    'N', 'CA', 'C', 'O',
    'H', 'HA', 'OXT',
    'H1', 'H2', 'H3',
    'HA2', 'HA3',
})


# ==========================================
# Atomic proton lookup — supermolecular subsystem parity check
# ==========================================
# Used by Phase 4 (qm_int_kcal_frozen) to determine binder subsystem charge
# from proton parity (Policy A: fail-fast on odd-electron subsystems).
# Includes biologically/structurally relevant elements; ELEMENT_CHARGE above is
# the legacy QM-region table (kept verbatim for backward compatibility).
_ELEMENT_PROTONS = {
    "H": 1, "C": 6, "N": 7, "O": 8, "S": 16, "P": 15,
    "F": 9, "Cl": 17, "Br": 35, "I": 53, "Si": 14, "Se": 34,
    "Fe": 26, "Zn": 30, "Mg": 12, "Ca": 20, "Na": 11, "K": 19,
}


def _select_sidechain_atoms(residue_atoms):
    """Return the sidechain-only subset of a residue's atom list.

    Used by ``partition_qmmm(mode="topology")`` to place target contact residue
    sidechains into the QM region while leaving their backbone in MM. This keeps
    the QM cluster chemically well-defined (Cα–Cβ bond is the clean break) and
    avoids dragging the whole backbone chain through the QM/MM boundary.

    Special cases:
        - GLY has no chemically meaningful sidechain → empty list.
        - PRO's pyrrolidine ring is fused to the backbone; splitting it creates
          a non-physical fragment → keep the whole residue.

    Args:
        residue_atoms: List of atom dicts belonging to a single residue
            (order/membership as produced by :func:`parse_pdb_atoms`).

    Returns:
        list[dict]: Subset of ``residue_atoms`` containing only sidechain atoms.
    """
    resname = residue_atoms[0].get('resname', '') if residue_atoms else ''
    if resname == 'GLY':
        return []
    if resname == 'PRO':
        return list(residue_atoms)
    return [a for a in residue_atoms if a.get('name', '').strip() not in _BACKBONE_NAMES]


# ==========================================
# [v4 S-2] QM/MM Link Atom (H-cap) 경계 보정
# ==========================================
def _add_link_atoms(qm_atoms, mm_atoms, all_atoms, cn_max_distance=1.7,
                    include_cb_breaks=False, binder_topology="linear", binder_chain="B"):
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
        include_cb_breaks: True 이면 Cα-Cβ 경계에도 H-cap 을 삽입한다.
            topology mode (sidechain-only QM) 에서 target contact residue 의
            Cβ (+sidechain) 가 QM, Cα (+backbone) 가 MM 인 경우에 사용된다.
            방향은 Cβ→Cα 이며 거리는 C-H 평균 1.09 Å 이다.
        binder_topology: Binder 의 위상 구조.
            - ``"linear"`` (기본): 선형 peptide, 추가 처리 없음.
            - ``"cyclic_htc"`` / ``"cyclic_cc"``: Head-to-Tail / C-to-C cyclization.
              binder chain 내 마지막 잔기의 C 와 첫 잔기의 N 사이 peptide bond 를
              검사하여 QM/MM 경계에 걸치면 H-cap link atom 을 삽입한다.
            - ``"cyclic_ss"``: Disulfide bond. binder chain 내 Cys Sγ-Sγ
              (≤ 2.1 Å) 쌍이 QM/MM 경계를 가로지르면 H-cap 을 삽입한다.
            topology mode (binder whole-QM) 에서는 cyclic bond 가 QM 내부에 있어
            추가 link atom 이 생성되지 않는다. legacy FAST/PLAID 모드에서 binder
            가 분할될 때에만 유효하다.
        binder_chain: Binder chain ID. cyclic topology 처리 시
            대상 chain 을 특정하기 위해 사용된다 (기본 "B").

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

    # Cyclic topology: last residue C → first residue N bond
    #   Head-to-Tail (cyclic_htc) / C-to-C (cyclic_cc) cyclization 에서 선형 루프가
    #   놓치는 마지막→첫 peptide bond 를 보정한다. binder 가 QM/MM 경계를 가로지르는
    #   legacy FAST/PLAID 모드에서만 유효 (topology mode 는 binder whole-QM).
    if binder_topology in ("cyclic_htc", "cyclic_cc"):
        binder_resnums = sorted(
            set(a["resnum"] for a in all_atoms if a["chain"] == binder_chain)
        )
        if len(binder_resnums) >= 2:
            first_res = binder_resnums[0]
            last_res = binder_resnums[-1]
            c_last = next(
                (a for a in all_atoms
                 if a["chain"] == binder_chain and a["resnum"] == last_res
                 and a["name"] == "C"),
                None,
            )
            n_first = next(
                (a for a in all_atoms
                 if a["chain"] == binder_chain and a["resnum"] == first_res
                 and a["name"] == "N"),
                None,
            )
            if c_last and n_first:
                d_cyc = np.linalg.norm(
                    np.array([c_last["x"], c_last["y"], c_last["z"]])
                    - np.array([n_first["x"], n_first["y"], n_first["z"]])
                )
                if d_cyc <= cn_max_distance:
                    c_key = (c_last["chain"], c_last["resnum"], "C")
                    n_key = (n_first["chain"], n_first["resnum"], "N")
                    # C ∈ QM, N ∈ MM 경계
                    if c_key in qm_keys and n_key in mm_keys:
                        lk = _make_link_h(c_last, n_first, BOND_LEN_CH)
                        lk["link_type"] = "cyclic_htc_c_to_n"
                        link_atoms.append(lk)
                    # N ∈ QM, C ∈ MM 경계
                    if n_key in qm_keys and c_key in mm_keys:
                        lk = _make_link_h(n_first, c_last, BOND_LEN_CH)
                        lk["link_type"] = "cyclic_htc_n_to_c"
                        link_atoms.append(lk)

    # Disulfide cyclization: Cys Sγ-Sγ 경계 처리
    #   binder 내부 SS bond 이 QM 내부에 완전히 포함되면 추가 처리 불필요.
    #   QM/MM 경계를 가로지르는 SS bond 에만 H-cap 삽입 (C-S 근사 거리 1.34 Å).
    if binder_topology == "cyclic_ss":
        qm_sg = [a for a in qm_atoms
                 if a["chain"] == binder_chain and a["name"] == "SG"]
        mm_sg = [a for a in mm_atoms
                 if a["chain"] == binder_chain and a["name"] == "SG"]
        for qm_s in qm_sg:
            for mm_s in mm_sg:
                d_ss = np.linalg.norm(
                    np.array([qm_s["x"], qm_s["y"], qm_s["z"]])
                    - np.array([mm_s["x"], mm_s["y"], mm_s["z"]])
                )
                if d_ss <= 2.1:  # SS bond 상한 2.1 Å
                    lk = _make_link_h(qm_s, mm_s, 1.34)  # C-S 근사 결합길이
                    lk["link_type"] = "cyclic_ss"
                    link_atoms.append(lk)

    # Cα-Cβ 경계 H-cap (topology mode 에서 sidechain-only QM 시 사용)
    #   - Target contact residue 의 Cβ+sidechain 이 QM, Cα+backbone 이 MM 인 경우
    #   - 방향: Cβ→Cα, 거리: 1.09 Å (C-H 평균)
    if include_cb_breaks:
        n_cb_breaks = 0
        for atom in qm_atoms:
            if atom.get("name", "").strip() != "CB":
                continue
            cb_chain = atom["chain"]
            cb_resnum = atom["resnum"]
            # 같은 잔기의 CA 를 MM 측에서 탐색
            ca = next(
                (a for a in mm_atoms
                 if a["chain"] == cb_chain and a["resnum"] == cb_resnum
                 and a.get("name", "").strip() == "CA"),
                None,
            )
            if ca is None:
                continue
            # _make_link_h 와 동일한 규약: qm_atom(=CB) 에서 mm_atom(=CA) 방향으로 bond_len
            link = _make_link_h(atom, ca, BOND_LEN_CH)
            link["link_type"] = "ca_cb"  # 진단용 flag (기존 peptide C-N 절단과 구분)
            link_atoms.append(link)
            n_cb_breaks += 1
        if n_cb_breaks:
            print(f"  [Link Atom] Cα-Cβ 경계 {n_cb_breaks}개 → H-cap 삽입 (topology mode)")

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
        "is_link": True,                    # build_qm_mol 진단용 플래그 (legacy)
        "is_link_atom": True,               # subsystem 분리용 표준 플래그
    }


def partition_qmmm(atoms, qm_cutoff=5.0, mode="full", binder_chain="B", target_card=None):
    """
    mode="full": 바인더 전체를 QM으로.
    mode="fast": 타겟과 바인더가 맞닿은 4.0 Å 이내 핵심 인터페이스 잔기만 추출하여 QM으로.
    mode="plaid" : 타겟과 바인더가 맞닿은 3.0 Å 이내 핵심 인터페이스 잔기만 추출하여 QM으로.
    mode="topology" : ``target_card`` 에 사전 정의된 contact residue 만 QM 으로.
        - binder chain 의 모든 잔기는 whole-residue QM.
        - target chain 중 ``target_contact_residues`` 에 나열된 잔기만 sidechain-only QM
          (GLY 는 제외, PRO 는 whole residue).
        - ``cofactor_residues`` 에 선언된 잔기는 MM point charge 처리.
        - ``qm_atom_budget.max_n_qm`` 을 초과하면 RuntimeError 로 즉시 실패한다.
        - Snapshot 간 QM region 이 topology 에 고정되어 n_qm variance = 0 을 보장한다.
    """
    # Topology mode — 사전 정의 contact residue 기반 QM 영역 (snapshot-invariant)
    if mode == "topology":
        if target_card is None:
            raise ValueError("topology mode requires target_card argument")

        rules = target_card.get("partition_rules", {})
        contact_resnums = set(target_card.get("target_contact_residues", []))
        # whole_residue_exceptions: target contact residues that
        # must be treated as whole-residue QM despite the default sidechain_only
        # policy (e.g. GLY_60 DxxGQ motif backbone N-H contributes to KD2 binding).
        whole_res_exceptions = set(target_card.get("whole_residue_exceptions", []))
        cofactor_residues = target_card.get("cofactor_residues", [])
        cofactor_resnames = {c.get("resname") for c in cofactor_residues if c.get("resname")}
        tc_binder_chain = target_card.get("binder_chain", binder_chain)
        target_chain = target_card.get("target_chain", "A")
        max_n_qm = target_card.get("qm_atom_budget", {}).get("max_n_qm", 600)
        if "qm_atom_budget" not in target_card:
            import warnings
            warnings.warn(
                f"target_card '{target_card.get('target_id')}' missing qm_atom_budget — "
                "using fallback max_n_qm=600. Regenerate with generate_target_card.py.",
                stacklevel=2,
            )

        # (chain, resnum) 로 원자를 잔기 단위로 그룹화
        from collections import defaultdict
        by_res = defaultdict(list)
        for a in atoms:
            by_res[(a["chain"], a["resnum"])].append(a)

        qm_atoms = []
        mm_atoms = []

        for (chain, resnum), res_atoms in by_res.items():
            resname = res_atoms[0].get("resname", "")

            if resname in cofactor_resnames:
                # Cofactor → MM point charge (GNP, MG 등)
                mm_atoms.extend(res_atoms)
            elif chain == tc_binder_chain:
                # Binder 는 whole residue QM (ncAA 포함 전체 잔기)
                qm_atoms.extend(res_atoms)
            elif chain == target_chain and resnum in contact_resnums:
                if resnum in whole_res_exceptions:
                    # Whole-residue QM (e.g., GLY backbone-mediated DxxGQ contact).
                    # Entire residue (backbone + sidechain) enters QM region.
                    qm_atoms.extend(res_atoms)
                else:
                    # Target contact 잔기 → sidechain 만 QM, backbone 은 MM
                    sc = _select_sidechain_atoms(res_atoms)
                    sc_ids = {id(a) for a in sc}
                    backbone = [a for a in res_atoms if id(a) not in sc_ids]
                    qm_atoms.extend(sc)
                    mm_atoms.extend(backbone)
            else:
                # 그 외 target 잔기 / 다른 chain → 전부 MM
                mm_atoms.extend(res_atoms)

        # QM atom budget 검증 — 초과 시 명시적 실패 (v0.6 Fail-fast 정책)
        if len(qm_atoms) > max_n_qm:
            raise RuntimeError(
                "QM atom count {} exceeds max_n_qm={} (target_card: {}). "
                "Reduce target_contact_residues or adjust qm_atom_budget.".format(
                    len(qm_atoms), max_n_qm, target_card.get("target_id")
                )
            )

        print(
            "  [Calculate Mode: TOPOLOGY] target_card={} | binder_chain={} | "
            "contacts={} | n_qm={} (budget {})".format(
                target_card.get("target_id"), tc_binder_chain,
                len(contact_resnums), len(qm_atoms), max_n_qm,
            )
        )
        return qm_atoms, mm_atoms

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

    # Policy A (fail-fast): 홀수 전자 → 즉시 예외 발생.
    # v0.5 까지의 silent charge=-1 보정은 서로 다른 quantum system 의 에너지
    # 혼용 원인이었으므로 제거한다. caller 는 OddElectronError 를 catch 하여
    # snapshot FAILED 로 기록하고 skip 해야 한다 (design abort 금지).
    if total_protons % 2 != 0:
        raise OddElectronError(
            f"Odd-electron QM region detected (total_protons={total_protons}). "
            f"n_qm_atoms={len(qm_atoms)}. "
            f"QM/MM boundary produces a chemically ambiguous open-shell fragment. "
            f"Recommended actions: (1) adjust target_contact_residues in target_card.json, "
            f"(2) check binder link-atom placement (binder_topology setting), "
            f"(3) verify titratable residue protonation states. "
            f"[Policy A: fail-fast — no silent charge correction applied]"
        )

    # Si 등 비표준 원소는 def2-SVP로 fallback
    # [FIX1] 모든 원소 포함 (H 포함), H는 6-31G** 사용
    # DF mode (Phase 5a): def2-* 계열 basis 가 전달되면 모든 원소에 일괄 적용.
    # 이유: def2-universal-jfit auxbasis 는 Karlsruhe family 와 매칭되어야 하며,
    # H 만 6-31G** 로 섞으면 Pople/Karlsruhe 혼합으로 DF 보조함수 정합성이 깨진다.
    # 비-DF (Pople) 경로는 기존 H=6-31G**, Si/Br/I=def2-SVP 특수처리 유지 (backward-compat).
    all_elements = set(a["element"] for a in qm_atoms if a["element"])
    if basis.startswith("def2"):
        basis_dict = {elem: basis for elem in all_elements}
    else:
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
# DFT 객체 설정 헬퍼 — 중복 제거
# ==========================================
def _make_dft_mf(mol, qm_xc, max_memory, use_df, df_auxbasis):
    """Configure RKS DFT object with optional density fitting.

    Phase 4 (qm_int_kcal_frozen) 의 binder/target 단독 SCF 와 main QM/MM
    SCF 모두에 동일한 DFT 설정을 적용하기 위해 추출한 헬퍼.

    GPU 변환(``to_gpu()``) 및 MM 점전하 부착(``mm_charge``)은 호출부 책임이며,
    본 함수는 순수 RKS 객체만 반환한다. 이는 ``density_fit()`` decorator chain
    이 GPU 변환·점전하 부착 *이전* 에 적용되어야 한다는 Phase 5a 제약 (UPDATE.md
    해당 버전) 을 호출부에서 명시 제어할 수 있도록 하기 위함이다.

    Args:
        mol: PySCF ``gto.M`` 분자 객체
        qm_xc: 교환-상관 범함수 이름 (예: ``"wb97x-d3"``)
        max_memory: PySCF ``max_memory`` (MB)
        use_df: density fitting 활성화 여부
        df_auxbasis: 보조 기저함수 (use_df=True 시 사용)

    Returns:
        pyscf.dft.RKS 객체 (density_fit 적용 후, GPU 변환·MM 부착 미적용)
    """
    mf = dft.RKS(mol)
    mf.xc = qm_xc
    mf.verbose = 1
    mf.max_memory = max_memory
    mf.grids.level = 3
    mf.max_cycle = 300
    mf.conv_tol = 1e-8
    if use_df:
        mf = mf.density_fit(auxbasis=df_auxbasis)
    return mf


# ==========================================
# QM/MM 계산 실행
# ==========================================
def run_qmmm_calc(pdb_path, output_dir, qm_basis, qm_xc, ncaa_elem, qm_cutoff=5.0, mode="full", binder_chain="B", target_id=None, use_df=True, df_auxbasis="def2-universal-jfit"):
    """단일 스냅샷 PDB에 대한 QM/MM 계산.

    ``target_id`` 가 지정되면 target_card JSON 을 로드하여
    topology mode 로 QM region 을 구성한다. 이 경우 ``mode`` 인자는 무시되고
    항상 topology 경로를 사용한다 (Phase 1: snapshot-invariant QM region).
    ``target_id=None`` 이면 기존 distance-based mode (full/fast/plaid) 유지.
    """

    # [#40 v0.2] max_memory 를 환경 변수 UPDD_QM_MAX_MEMORY(MB) 로 주입 (기본 16000 = 16GB).
    _max_mem = int(os.environ.get("UPDD_QM_MAX_MEMORY", "16000"))

    # Target card 로드 (topology mode 활성화 조건)
    target_card = None
    if target_id is not None:
        from utils_common import load_target_card
        target_card = load_target_card(target_id)

    # topology mode 여부 결정 — target_card 있으면 topology, 없으면 legacy mode 유지
    use_mode = "topology" if target_card is not None else mode

    basename = os.path.basename(pdb_path).replace(".pdb", "")
    out_json = os.path.join(output_dir, f"{basename}_qmmm_{use_mode}.json")

    if os.path.exists(out_json):
        # [v0.3 #1+#5] resume 무결성 강화:
        #   - converged=True AND energy_total_hartree 유효 → 스킵 (정상 완료)
        #   - converged=True 인데 interaction_kcal=None → QM-only 만 실패한 케이스
        #     → JSON 삭제 후 전체 재계산 (독립 SCF 재계산으로 신뢰성 보장)
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

    # topology mode: target_card 의 binder_chain 을 우선 사용
    effective_binder_chain = target_card["binder_chain"] if target_card else binder_chain
    qm_atoms, mm_atoms = partition_qmmm(
        atoms,
        qm_cutoff=qm_cutoff,
        mode=use_mode,
        binder_chain=effective_binder_chain,
        target_card=target_card,
    )

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
    # topology mode 에서는 Cα-Cβ 경계에도 H-cap 삽입 (sidechain-only QM)
    # binder_topology 전달 — cyclic peptide 의 last→first bond 처리
    include_cb = (use_mode == "topology")
    _binder_topo = (target_card.get("binder_topology", "linear")
                    if target_card else "linear")
    _binder_chain_for_link = (target_card.get("binder_chain", effective_binder_chain)
                              if target_card else effective_binder_chain)
    link_atoms = _add_link_atoms(
        qm_atoms, mm_atoms, atoms,
        include_cb_breaks=include_cb,
        binder_topology=_binder_topo,
        binder_chain=_binder_chain_for_link,
    )
    if link_atoms:
        qm_atoms = qm_atoms + link_atoms

    # [Phase 6 DIAG] QM region 통계 — Path 진단용 (n_qm/link/mode 추적)
    _n_link_diag = sum(1 for a in qm_atoms if a.get("is_link_atom"))
    _diag(
        "qm_region",
        snapshot=basename,
        n_qm=len(qm_atoms),
        n_qm_real=len(qm_atoms) - _n_link_diag,
        n_link=_n_link_diag,
        n_mm=len(mm_atoms),
        partitioning_mode=use_mode,
        target_id=(target_card.get("target_id") if target_card else None),
        binder_topology=_binder_topo,
    )

    try:
        mol = build_qm_mol(qm_atoms, basis=qm_basis)
    except OddElectronError as e:
        # [Phase 6 DIAG] OddElectronError 신호 — Path 가 grep 으로 패턴 추적
        _diag(
            "odd_electron_error",
            snapshot=basename,
            detail=str(e)[:200],
            policy="fail_fast_policy_a",
        )
        # Policy A: snapshot 실패 기록 후 skip (design abort 금지).
        # 결과 JSON 에 status/reason 을 명시 기록하여 downstream ranking 에서 격리.
        print(f"  [POLICY-A] OddElectronError → snapshot skipped: {basename}")
        print(f"  [POLICY-A] detail: {e}")
        failure_result = {
            "snapshot":            basename,
            "pdb_path":            pdb_path,
            "qm_method":           f"{qm_xc}/{qm_basis}",
            "n_qm_atoms":          len(qm_atoms),
            "n_mm_atoms":          len(mm_atoms),
            "ncaa_element":        ncaa_elem,
            "energy_total_hartree": None,
            "energy_qm_hartree":    None,
            "energy_qmmm_kcal":     None,
            "interaction_kcal":     None,
            "converged":            False,
            "is_run_mmgbsa":        False,
            "status":               "FAILED",
            "reason":               "odd_electron",
            "odd_electron_detail":  str(e),
        }
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(failure_result, f, indent=2)
        return failure_result
    except Exception as e:
        print(f"  [!] QM 분자 빌드 실패: {e}")
        return None

    print(f"  QM 기저함수: {qm_basis} | 범함수: {qm_xc}")
    print(f"  QM 원자 수 (중원자): {mol.natm}")

    # [Phase 6 DIAG] PySCF 분자 빌드 결과 — 전하/스핀/기저 추적
    _diag(
        "qm_mol_built",
        snapshot=basename,
        n_atoms_mol=int(mol.natm),
        charge=int(mol.charge),
        spin=int(mol.spin),
        basis=qm_basis,
    )

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

    # Phase 5a: DF-J/K density fitting (16GB VRAM 실행 가능성 복구).
    # 560-atom QM region direct-SCF → ~266GB VRAM 필요 → 16GB RTX 5070 Ti 불가.
    # density_fit() 는 반드시 to_gpu() / mm_charge() 이전에 호출해야 한다 (decorator chain 순서).
    if use_df:
        mf = mf.density_fit(auxbasis=df_auxbasis)
        print(f"  [DF-J/K] density_fit 활성화 (auxbasis={df_auxbasis})")

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

        # [Phase 6 DIAG] 백엔드 확정 신호 (GPU 경로)
        _diag(
            "compute_backend",
            snapshot=basename,
            backend="GPU",
            df_mode=bool(use_df),
            df_auxbasis=(df_auxbasis if use_df else None),
        )

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
        # Phase 5a: CPU fallback 에도 동일하게 DF-J/K 적용 (mm_charge 이전).
        # GPU/CPU 경로 간 method 동일성 보장 → interaction_kcal 일관성.
        if use_df:
            mf = mf.density_fit(auxbasis=df_auxbasis)
        if len(mm_coords_charges) > 0:
            mf = qmmm_itrf.mm_charge(mf, mm_coords_charges[:, :3], mm_coords_charges[:, 3])
            print(f"  [+] MM 점전하: {len(mm_coords_charges)}개 포함 (CPU 상태)")

        # [Phase 6 DIAG] 백엔드 확정 신호 (CPU fallback 경로)
        _diag(
            "compute_backend",
            snapshot=basename,
            backend="CPU",
            df_mode=bool(use_df),
            df_auxbasis=(df_auxbasis if use_df else None),
        )

    try:
        energy_total = mf.kernel()
        print(f"  QM/MM 에너지: {energy_total:.8f} Hartree")
    except Exception as e:
        print(f"  [!] DFT 수렴 실패: {e}")
        energy_total = None

    # [Phase 6 DIAG] SCF 수렴 결과 — converged/에너지 추적
    _diag(
        "scf_result",
        snapshot=basename,
        converged=energy_total is not None,
        energy_hartree=(round(float(energy_total), 8) if energy_total is not None else None),
        energy_kcal=(round(float(energy_total) * 627.509, 3) if energy_total is not None else None),
    )

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
        # Phase 5a: QM-only 도 동일 DF-J/K 적용 (mf 와 method 일치 보장).
        # interaction_kcal = (E_QMMM - E_QM) 의 일관된 차분 → DF 오차 cancel.
        if use_df:
            mf_qm = mf_qm.density_fit(auxbasis=df_auxbasis)

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

    # ── 6. Phase 4 qm_int_kcal_frozen — 초분자 Morokuma-style 분해 ──
    # ΔE_int = E(AB) − E(A) − E(B) (frozen geometry, no relaxation).
    # Ref: Boys & Bernardi 1970 Mol Phys 19,553; Kitaura & Morokuma 1976 IJQC 10,325.
    # topology mode 한정 — distance-based mode 에서는 binder/target 분리가 모호.
    qm_int_kcal_frozen = None
    e_binder_iso = None
    e_target_iso = None
    phase4_diagnostic = {
        "bsse_corrected": False,
        "bsse_uncertainty_kcal": "+3 to +8 (def2-SVP estimate)",
        "geometry": "frozen_md_snapshot",
        "relaxation_neglected": True,
        "link_atom_consistency": "shared_between_complex_and_target_iso",
        "interaction_definition": "supermolecular_morokuma_frozen",
    }

    if e_qm_only is not None and target_card is not None:
        tc_binder_chain_p4 = target_card.get("binder_chain", binder_chain)
        tc_target_chain_p4 = target_card.get("target_chain", "A")
        target_iso_net_charge = target_card.get("target_iso_net_charge", None)

        # qm_atoms 는 link_atoms 합산본 (L792). subsystem 분리 시 link 원자는
        # E_complex 와 E_target_iso 에만 포함, E_binder_iso 에는 포함 안 함.
        binder_qm = [a for a in qm_atoms
                     if a.get("chain") == tc_binder_chain_p4 and not a.get("is_link_atom")]
        target_qm_raw = [a for a in qm_atoms
                         if a.get("chain") == tc_target_chain_p4 and not a.get("is_link_atom")]
        link_atoms_in_qm = [a for a in qm_atoms if a.get("is_link_atom")]

        if target_iso_net_charge is None:
            import warnings as _warnings
            _warnings.warn(
                "target_card missing 'target_iso_net_charge' — qm_int_kcal_frozen skipped. "
                "Add this field to target_cards/{}.json".format(target_card.get("target_id")),
                stacklevel=2,
            )
        elif len(binder_qm) > 0 and len(target_qm_raw) > 0:
            try:
                # Binder subsystem charge — proton parity (Policy A: fail-fast)
                binder_protons = sum(
                    _ELEMENT_PROTONS.get(a["element"].capitalize(), 0)
                    for a in binder_qm if a.get("element")
                )
                if binder_protons % 2 != 0:
                    raise OddElectronError(
                        f"Binder subsystem odd-electron (protons={binder_protons}). "
                        "Cannot compute qm_int_kcal_frozen."
                    )
                binder_charge = 0  # even protons → standard whole-peptide charge 0

                # E_binder_iso: binder atoms only, no link atoms (whole-residue,
                # no Cα-Cβ cuts → no H-cap needed for binder side).
                mol_binder = build_qm_mol(binder_qm, charge=binder_charge,
                                           spin=0, basis=qm_basis)
                mf_binder = _make_dft_mf(mol_binder, qm_xc, _max_mem, use_df, df_auxbasis)
                if gpu_success:
                    try:
                        mf_binder = mf_binder.to_gpu()
                    except Exception:
                        pass
                print(f"  [Phase4] E_binder_iso SCF "
                      f"({mol_binder.natm} atoms, charge={binder_charge})")
                e_binder_iso = mf_binder.kernel()

                # E_target_iso: target sidechains + SAME H-cap link atoms used in
                # E_complex. Link atom set must be identical between E_complex and
                # E_target_iso to cancel boundary
                # artefacts in the difference.
                target_with_hcaps = target_qm_raw + link_atoms_in_qm
                mol_target = build_qm_mol(target_with_hcaps,
                                           charge=target_iso_net_charge, spin=0,
                                           basis=qm_basis)
                mf_target = _make_dft_mf(mol_target, qm_xc, _max_mem, use_df, df_auxbasis)
                if gpu_success:
                    try:
                        mf_target = mf_target.to_gpu()
                    except Exception:
                        pass
                print(f"  [Phase4] E_target_iso SCF "
                      f"({mol_target.natm} atoms, charge={target_iso_net_charge})")
                e_target_iso = mf_target.kernel()

                if e_binder_iso is not None and e_target_iso is not None:
                    qm_int_kcal_frozen = (e_qm_only - e_binder_iso - e_target_iso) * 627.509
                    print(f"  [Phase4] qm_int_kcal_frozen = "
                          f"{qm_int_kcal_frozen:.4f} kcal/mol")

            except OddElectronError as _oe:
                print(f"  [Phase4] OddElectronError in subsystem → "
                      f"qm_int_kcal_frozen=None: {_oe}")
                _diag(
                    "odd_electron_error",
                    snapshot=basename,
                    detail=str(_oe)[:200],
                    policy="fail_fast_policy_a",
                    subsystem="binder_iso",
                )
            except Exception as _e4:
                print(f"  [Phase4] qm_int_kcal_frozen failed: {_e4}")

    # [Phase 6 DIAG] 상호작용 에너지 — frozen / legacy / BSSE 진단
    _diag(
        "int_energy",
        snapshot=basename,
        qm_int_kcal_frozen=(
            round(float(qm_int_kcal_frozen), 4) if qm_int_kcal_frozen is not None else None
        ),
        interaction_kcal=(round(float(e_interact), 4) if e_interact is not None else None),
        e_binder_iso_hartree=(
            round(float(e_binder_iso), 8) if e_binder_iso is not None else None
        ),
        e_target_iso_hartree=(
            round(float(e_target_iso), 8) if e_target_iso is not None else None
        ),
        bsse_uncertainty="3_to_8_kcal_def2svp",
    )

    # ── 7. 결과 저장 ─────────────────────────────────────────
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
        # Phase 5a: DF-J/K 진단 필드 — Path 진단 투명성 (어떤 DF mode/auxbasis 인지 추적).
        "df_mode":                use_df,
        "df_auxbasis":            df_auxbasis if use_df else None,
        # Phase 4: 초분자 QM interaction energy (frozen geometry).
        # topology mode 에서만 계산됨 — legacy mode 는 None.
        "qm_int_kcal_frozen":     qm_int_kcal_frozen,
        "e_binder_iso_hartree":   e_binder_iso,
        "e_target_iso_hartree":   e_target_iso,
        **phase4_diagnostic,
    }

    # Topology mode diagnostic fields — Path 진단 투명성 확보.
    # target_card 로드 경로를 명시 기록하여 snapshot 별 QM region 구성 추적 가능.
    if target_card is not None:
        tc_target_chain = target_card.get("target_chain", "A")
        tc_binder_chain = target_card["binder_chain"]
        result.update({
            "partitioning_mode":      "topology",
            "target_card_id":         target_card["target_id"],
            "target_card_version":    target_card["schema_version"],
            "n_qm_binder":            sum(
                1 for a in qm_atoms
                if a.get("chain") == tc_binder_chain and not a.get("is_link")
            ),
            "n_qm_target_sidechain":  sum(
                1 for a in qm_atoms
                if a.get("chain") == tc_target_chain and not a.get("is_link")
            ),
            "n_link_atoms":           sum(1 for a in qm_atoms if a.get("is_link")),
            "cofactor_treatment":     "mm_point_charge",
        })

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
    # Phase 5a: --qm_basis 기본값 6-31G* → def2-SVP (DF-J/K 매칭).
    # 기존 6-31G* 재현이 필요하면 명시적으로 `--no-df --qm_basis 6-31G*` 사용.
    parser.add_argument("--qm_basis",  default="def2-SVP",    help="QM 기저함수 (기본: def2-SVP, DF mode)")
    parser.add_argument("--qm_xc",     default="wb97x-d3",    help="DFT 범함수 (기본: wb97x-d3)")
    parser.add_argument("--ncaa_elem", default="none",        help="ncAA 핵심 원소 (예: Si)")
    parser.add_argument("--cutoff",    type=float, default=5.0, help="QM 영역 cutoff Å (기본 5.0)")
    parser.add_argument("--mode", default="fast", choices=["full", "fast", "plaid"], help="QM/MM 계산 모드")
    parser.add_argument("--filter", default=None, help="특정 디자인 ID만 골라서 계산 (예: design_w2_2_s1)")
    parser.add_argument("--binder_chain", default="B", help="Binder chain ID (기본 B)")
    # Phase 1: target_card 기반 topology mode 활성화. 지정 시 --mode 는 무시되고
    # target_cards/{target_id}.json 의 QM partitioning 정책을 따른다 (snapshot-invariant QM region).
    parser.add_argument("--target-id", dest="target_id", type=str, default=None,
                        help="Target card ID (e.g., 6WGN). Activates topology mode.")
    # Phase 5a: DF-J/K density fitting 제어. 기본 활성 (16GB VRAM 필수).
    # `--no-df` 는 소규모 QM region (예: <100 atoms) 에서 direct-SCF 비교 검증용.
    parser.add_argument("--no-df", dest="use_df", action="store_false",
                        help="DF-J/K density fitting 비활성화 (소규모 QM region에서만 사용)")
    parser.set_defaults(use_df=True)
    parser.add_argument("--df-auxbasis", dest="df_auxbasis", default="def2-universal-jfit",
                        help="DF 보조 기저함수 (기본: def2-universal-jfit)")
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
            binder_chain = args.binder_chain,
            target_id  = args.target_id,
            # Phase 5a: DF-J/K 인자 전달.
            use_df       = args.use_df,
            df_auxbasis  = args.df_auxbasis,
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
