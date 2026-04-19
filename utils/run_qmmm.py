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

# ==========================================
# [B8] GPU4PySCF 16GB VRAM boundary config — SciVal-approved (2026-04-18)
# ==========================================
# RTX 5070 Ti 16GB 등 16GB VRAM 경계에서 gpu4pyscf 의 기본 ``min_ao_blksize`` (보통
# 128-256) 가 J/K 중간 버퍼 blksize 를 너무 크게 잡아 OOM 을 유발한다. 32GB+ GPU 에서만
# 256 이 안전하므로 16GB 에서는 64 로 하향.
# 반드시 ``import gpu4pyscf`` / ``cupy`` **이전** 에 설정되어야 한다
# (config 가 import 시점에 read-once 됨).
import pyscf as _early_pyscf  # noqa: E402 — config pre-set, other pyscf imports follow
_early_pyscf.__config__.min_ao_blksize = int(os.environ.get("UPDD_MIN_AO_BLKSIZE", "64"))

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
# DF auto-mode — VRAM 기반 자동 결정 (v0.6.x)
# ==========================================
def _decide_df_auto():
    """VRAM 총량 기반 DF 자동 결정.

    Threshold: ``UPDD_DF_AUTO_VRAM_THRESHOLD_GB`` 환경변수 (default 20 GiB).

    - VRAM 총량 >= threshold: DF on (density fitting 이 VRAM 수용)
    - VRAM 총량 <  threshold: DF off (direct-SCF 로 VRAM peak 회피)
    - VRAM 감지 실패         : DF off (보수적 선택, direct-SCF 는 항상 작동)

    Note:
        ``import cupy`` 은 함수 내부에서만 수행한다. cupy 가 없는 환경
        (예: CPU-only smoke test, `--help` 출력) 에서도 argparse 가
        import 오류 없이 정상 동작해야 하기 때문.

    Returns:
        Tuple of (use_df: bool, reason: str).
    """
    threshold_gb = float(os.environ.get("UPDD_DF_AUTO_VRAM_THRESHOLD_GB", "20"))
    try:
        import cupy as _cp
        props = _cp.cuda.runtime.getDeviceProperties(0)
        total_bytes = props["totalGlobalMem"]
        total_gib = total_bytes / (1024 ** 3)
        if total_gib >= threshold_gb:
            return True, f"auto: VRAM {total_gib:.1f} GiB >= threshold {threshold_gb} GiB → DF on"
        else:
            return False, f"auto: VRAM {total_gib:.1f} GiB < threshold {threshold_gb} GiB → DF off"
    except Exception as e:
        return False, f"auto: VRAM 감지 실패 ({type(e).__name__}: {str(e)[:60]}) → DF off (conservative)"


# DF 모드에서 J/K 중간 버퍼 제한 (gpu4pyscf 가 blksize 결정에 참고)
# gpu4pyscf.df.df.cholesky_eri_gpu() L140 에서
# ``cupy.zeros([blksize, nao, nao])`` peak buffer 를 생성하는데,
# blksize 는 with_df.max_memory 와 nao 로 결정된다.
# n_qm=466, nao~5000 기준 4000 MB → blksize ~40 → VRAM peak <6 GB 기대.
# 또한 ``with_df.use_gpu_memory = False`` 를 병행 강제하여 CDERI 를 CPU
# pinned RAM 에 저장한다 (auto heuristic 재확인; 소형 시스템 보호).
UPDD_DF_JK_MAX_MEMORY_MB = int(os.environ.get("UPDD_DF_JK_MAX_MEMORY_MB", "4000"))

# ==========================================
# SciVal 승인 축소 번들 (2026-04-18) — ranking preservation 기준
# ==========================================
# gpu4pyscf DF 가 n_qm=466 에서 공식 지원 범위 (Wu 2024 arXiv:2404.09452 벤치 <=168
# atoms) 밖으로 direct-SCF 경로 wall-time 감소 (164 s/cycle → 130 s/cycle 내외) 목표.
# 절대 binding E 는 설계 간 동일하게 최대 0.5 kcal/mol shift 하지만 ranking 은 유지
# (SciVal verified). 각 항목은 환경변수 override 로 rollback/debug 가능.
#
# - B2 ``grids.level=2``   : NWCHEM prune / ORCA DefGrid1 수준 (default 3 → 2)
# - B3 ``conv_tol=1e-7``   : binding E ranking 에 충분 (PySCF default 1e-9 → 1e-7)
# - B5 ``direct_scf_tol=1e-12`` : ERI screening 완화 (default 1e-13 → 1e-12)
# - B6 ``init_guess='huckel'``  : Protein+peptide convergence 가속 (default 'minao')
# - B9 ``UPDD_VRAM_POOL_FRACTION=0.65`` : CuPy pool 상한 (VRAM 대비, default 0.65)
#
# B1 (rks_lowmem) 과 B7 (mo_coeff snapshot 재사용) 은 SciVal 거부로 미구현:
#   - B1: QM/MM 구조 비호환 (gpu4pyscf/qmmm/itrf.py L85 square Fock vs hf_lowmem
#         tril 반환 → shape mismatch).
#   - B7: snapshot 간 statistical independence 위반.
UPDD_DFT_GRIDS_LEVEL = int(os.environ.get("UPDD_DFT_GRIDS_LEVEL", "2"))
UPDD_SCF_CONV_TOL = float(os.environ.get("UPDD_SCF_CONV_TOL", "1e-7"))
UPDD_SCF_DIRECT_TOL = float(os.environ.get("UPDD_SCF_DIRECT_TOL", "1e-12"))
UPDD_SCF_INIT_GUESS = os.environ.get("UPDD_SCF_INIT_GUESS", "minao")
UPDD_VRAM_POOL_FRACTION = float(os.environ.get("UPDD_VRAM_POOL_FRACTION", "0.65"))


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


class ChargeDeclarationMismatch(ValueError):
    """Raised when target_card declared QM net charge does NOT match the
    chemistry-computed (topology-based) value.

    R-15 (CLAUDE.md, 2026-04-19): declared charge must equal the net charge
    derived from residue identity + atom-name pattern at the pipeline's
    declared pH. A magnitude mismatch implies at least a 1-electron
    difference that parity check alone cannot detect.

    The caller should fail fast (do NOT silent-correct) and record the
    residue-level breakdown in the diagnostic JSON for later audit.
    """
    pass


class BinderChargeMismatch(ValueError):
    """Raised when ``target_card.binder_net_charge`` (hint) does NOT match
    the value computed from the runtime snapshot PDB.

    R-16 (CLAUDE.md, 2026-04-19): binder net charge is derived at runtime
    from snapshot PDB (residue identity + cyclic/linear detection + atom
    patterns). ``target_card.binder_net_charge`` is an "expected" hint.
    Runtime PDB is the ground truth; mismatch is a SSOT violation.
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


# ==========================================
# R-15 / R-16 — Chemistry-true charge computation (delegated to charge_topology)
# ==========================================
# Implementation lives in ``utils/charge_topology.py`` so the audit module
# and unit tests can reuse it without importing pyscf. See SciVal
# verdict_target_card_charge_rationale_20260419_v2.md §1.1, §5.
from charge_topology import (  # noqa: E402
    compute_binder_chem_charge,
    compute_qm_net_charge_topology,
)


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
        # Stage 2 (Option B3) — partition_rules.target_residues_mode 지원:
        #   "whole"             : 모든 contact residue 를 whole-residue QM 로
        #                         (chemistry-consistent; cyclic peptide binder 의
        #                          unequal ARG/ASP/GLU 수를 허용).
        #   "sidechain_from_cb" : legacy (Stage 1 이전). sidechain 만 QM, backbone MM.
        # whole_residue_exceptions 는 schema 0.6.3 (sidechain mode) 의 per-residue
        # override 였다. 0.6.4+ whole mode 에서는 의미 없음 → 무시/경고.
        target_residues_mode = rules.get("target_residues_mode", "sidechain_from_cb")
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
        if target_residues_mode == "whole" and whole_res_exceptions:
            import warnings as _w2
            _w2.warn(
                "target_residues_mode=whole makes whole_residue_exceptions "
                "redundant; field is being ignored. Regenerate target_card "
                "with schema 0.6.4 (generate_target_card.py).",
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
                # Stage 2 (2026-04-19) — target_residues_mode dispatch.
                if target_residues_mode == "whole":
                    # Option B3 — 모든 contact residue 전체가 QM.
                    # Sidechain_from_cb 의 Cα-Cβ artificial cut 없음 →
                    # 전하 consistency + link atom 수 대폭 감소.
                    qm_atoms.extend(res_atoms)
                elif resnum in whole_res_exceptions:
                    # Legacy per-residue override (schema 0.6.3 path).
                    qm_atoms.extend(res_atoms)
                else:
                    # Legacy sidechain-only path: sidechain QM, backbone MM.
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
            "contacts={} | target_residues_mode={} | n_qm={} (budget {})".format(
                target_card.get("target_id"), tc_binder_chain,
                len(contact_resnums), target_residues_mode,
                len(qm_atoms), max_n_qm,
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

    charge 는 caller 가 binder_net_charge + target_iso_net_charge 합산하여 전달.
    홀수 전자 감지 시 OddElectronError (Policy A fail-fast): n_electrons = total_protons - charge.
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
    # 전자 수 = 양성자 수 - QM 전하. charge 가 올바르게 전달되어야
    # 이온화된 잔기(ASP/GLU 등)를 포함한 QM region 에서 closed-shell 판정이 정확하다.
    n_electrons = total_protons - charge
    if n_electrons % 2 != 0:
        raise OddElectronError(
            f"Odd-electron QM region detected "
            f"(total_protons={total_protons}, charge={charge}, n_electrons={n_electrons}). "
            f"n_qm_atoms={len(qm_atoms)}. "
            f"QM/MM boundary produces a chemically ambiguous open-shell fragment. "
            f"Recommended actions: (1) verify binder_net_charge + target_iso_net_charge in target_card.json, "
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
# DF OOM 회피 가드 — gpu4pyscf 버전 호환 try/except
# ==========================================
def _apply_df_oom_guard(mf):
    """Configure mf.with_df to avoid gpu4pyscf GPU OOM at large n_qm.

    Two knobs (both guarded for gpu4pyscf version compatibility):

    1. ``with_df.use_gpu_memory = False``
       Force CDERI (3-center Cholesky integrals) to CPU pinned RAM instead of
       VRAM. gpu4pyscf auto heuristic picks this for large systems but we pin
       it explicitly to protect smaller systems from regressions when the
       auto threshold shifts between releases.

    2. ``with_df.max_memory = UPDD_DF_JK_MAX_MEMORY_MB``
       Bound the J/K-build intermediate buffer ``cupy.zeros([blksize, nao,
       nao])`` by capping ``max_memory`` so gpu4pyscf picks a smaller
       ``blksize``. At n_qm=466, nao~5000, 4000 MB → blksize ~40 → VRAM peak
       <6 GB (instead of 10-20 GB at default blksize).

    Both writes are ``try/except`` guarded because ``mf.with_df`` may not
    expose these attributes on every gpu4pyscf 1.x patch release. The
    numerical result of the SCF is invariant under these settings — only
    memory layout and block size change (Cholesky decomposition itself is
    deterministic).

    Args:
        mf: A PySCF / gpu4pyscf SCF object *after* ``density_fit()`` has been
            applied (so ``mf.with_df`` exists).
    """
    try:
        mf.with_df.use_gpu_memory = False
    except Exception:
        pass  # 속성 없는 버전은 무시
    try:
        mf.with_df.max_memory = UPDD_DF_JK_MAX_MEMORY_MB
    except Exception:
        pass


# ==========================================
# SciVal 승인 축소 번들 적용 — B2/B3/B5/B6 중앙화
# ==========================================
def _apply_scf_bundle(mf):
    """Apply SciVal-approved SCF/DFT tuning bundle (B2/B3/B5/B6) to ``mf``.

    Ranking preservation 기준으로 검증됨 (2026-04-18). 절대 binding E 는 모든
    design 에 동일 shift 하므로 상대 순위는 유지.

    - B2: ``mf.grids.level`` = ``UPDD_DFT_GRIDS_LEVEL`` (default 2)
    - B3: ``mf.conv_tol`` = ``UPDD_SCF_CONV_TOL`` (default 1e-7)
    - B5: ``mf.direct_scf_tol`` = ``UPDD_SCF_DIRECT_TOL`` (default 1e-12)
    - B6: ``mf.init_guess`` = ``UPDD_SCF_INIT_GUESS`` (default 'minao' —
          'huckel' 은 HOMO-LUMO near-degeneracy 시 발산 유발, opt-in 으로 분리)

    각 속성 쓰기는 ``try/except`` 가드로 감싸, PySCF/gpu4pyscf 의 특정 method
    객체 (예: 일부 사용자 정의 wrapper) 가 속성을 노출하지 않는 경우에도
    import-time 오류 없이 통과. 이 함수는 ``idempotent`` — 동일 mf 에 반복
    호출해도 부작용 없음.

    Args:
        mf: PySCF/gpu4pyscf mean-field 객체. ``density_fit()`` / ``to_gpu()``
            등 decorator 체인 뒤에 호출해도 무방 (idempotent).

    Returns:
        The same ``mf`` object for convenient chaining.
    """
    try:
        mf.grids.level = UPDD_DFT_GRIDS_LEVEL
    except Exception:
        pass
    try:
        mf.conv_tol = UPDD_SCF_CONV_TOL
    except Exception:
        pass
    try:
        mf.direct_scf_tol = UPDD_SCF_DIRECT_TOL
    except Exception:
        pass
    try:
        mf.init_guess = UPDD_SCF_INIT_GUESS
    except Exception:
        pass
    return mf


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
    mf.verbose = 4
    mf.max_memory = max_memory
    mf.max_cycle = 300
    # SciVal 승인 축소 번들 (B2/B3/B5/B6) — grids.level / conv_tol /
    # direct_scf_tol / init_guess 를 환경변수 기반으로 일괄 적용.
    # _make_dft_mf 는 Phase 4 binder/target 단독 SCF 의 진입점이기도 하므로
    # 이 한 곳에서 번들을 걸면 모든 하위 SCF 가 일관된 설정으로 수렴.
    _apply_scf_bundle(mf)
    if use_df:
        mf = mf.density_fit(auxbasis=df_auxbasis)
        # gpu4pyscf DF OOM 회피:
        # (1) CDERI 를 pinned RAM 에 명시 저장 (auto heuristic 재확인, 소형 시스템 보호)
        # (2) max_memory 축소로 J/K 중간 버퍼 [blksize, nao, nao] blksize 강제 축소
        _apply_df_oom_guard(mf)
        # density_fit decorator 가 grids / conv_tol / init_guess 를 감싸 wrapper
        # 객체로 교체할 수 있어 번들 재적용 (idempotent).
        _apply_scf_bundle(mf)
    return mf


# ==========================================
# QM/MM 계산 실행
# ==========================================
def run_qmmm_calc(pdb_path, output_dir, qm_basis, qm_xc, ncaa_elem, qm_cutoff=5.0, mode="full", binder_chain="B", target_id=None, use_df=True, df_auxbasis="weigend"):
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

        # Stage 2 (2026-04-19) — numbering_convention 소비 + 진단.
        # schema 0.6.4 부터 필수. 0.6.3 이하 card 는 offset=0 default + 경고.
        nc = target_card.get("numbering_convention")
        if isinstance(nc, dict):
            _nc_source = str(nc.get("source", "MD"))
            _nc_offset = int(nc.get("md_to_crystal_offset", 0))
            _nc_note = str(nc.get("note", ""))
            _nc_missing = False
        else:
            import warnings as _wnc
            _wnc.warn(
                "target_card '%s' is missing the v0.6.4 numbering_convention "
                "object — defaulting to source=MD, md_to_crystal_offset=0. "
                "Regenerate with generate_target_card.py to remove this warning."
                % target_card.get("target_id"),
                stacklevel=2,
            )
            _nc_source = "MD"
            _nc_offset = 0
            _nc_note = "missing in pre-0.6.4 card — assumed zero offset, verify!"
            _nc_missing = True
        _diag(
            "numbering_convention",
            target_id=target_card.get("target_id"),
            source=_nc_source,
            md_to_crystal_offset=_nc_offset,
            missing_from_card=bool(_nc_missing),
            schema_version=target_card.get("schema_version", "unknown"),
            note=_nc_note[:200],
        )

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
    try:
        qm_atoms, mm_atoms = partition_qmmm(
            atoms,
            qm_cutoff=qm_cutoff,
            mode=use_mode,
            binder_chain=effective_binder_chain,
            target_card=target_card,
        )
    except RuntimeError as e:
        # max_n_qm 초과 등 partition 단계 실패 — 프로세스 crash 방지
        print(f"  [!] QM region 구성 실패: {e}")
        _diag("qm_region", snapshot=basename, n_qm=-1, n_qm_real=-1, n_link=-1,
              n_mm=-1, partitioning_mode=use_mode,
              target_id=(target_card.get("target_id") if target_card else None),
              binder_topology="unknown", error=str(e)[:120])
        failure_result = {
            "snapshot": basename, "pdb_path": pdb_path,
            "qm_method": f"{qm_xc}/{qm_basis}",
            "energy_total_hartree": None, "energy_qmmm_kcal": None,
            "interaction_kcal": None, "converged": False, "is_run_mmgbsa": False,
            "status": "FAILED", "reason": "partition_error", "partition_error": str(e),
            "regime": "ranking_only",
        }
        if out_json:
            with open(out_json, "w", encoding="utf-8") as f:
                json.dump(failure_result, f, indent=2)
        return failure_result

    if len(qm_atoms) == 0:
        return {
            "snapshot": basename, "energy_total_hartree": None, "energy_qmmm_kcal": None,
            "interaction_kcal": None, "converged": False, "is_run_mmgbsa": False,
            "reason": "no_qm_atoms",
        }

    # [v4 S-2] QM/MM 경계 H-cap link atom 추가 (Senn & Thiel 2009)
    # topology mode 에서는 target_residues_mode 에 따라 경계가 달라진다:
    #   - "sidechain_from_cb" (legacy): Cα-Cβ 절단 → include_cb_breaks=True
    #   - "whole"               (B3) : peptide-bond 경계만 → include_cb_breaks=False
    # binder_topology 전달 — cyclic peptide 의 last→first bond 처리.
    _t_rules = (target_card.get("partition_rules", {}) if target_card else {})
    _target_res_mode_for_link = _t_rules.get("target_residues_mode",
                                             "sidechain_from_cb")
    include_cb = (use_mode == "topology"
                  and _target_res_mode_for_link != "whole")
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
        target_residues_mode=_target_res_mode_for_link,
        target_id=(target_card.get("target_id") if target_card else None),
        binder_topology=_binder_topo,
    )

    # QM region 순 전하: binder + target sidechain 이온화 상태 (pH 7.4)
    # target_iso_net_charge 는 target_card 에서 사전 결정 (generate_target_card.py).
    # 부재 시 0 으로 fallback 하되 경고 발행.
    _binder_charge_declared = int(target_card.get("binder_net_charge", 0)) if target_card else 0
    _target_iso_charge = target_card.get("target_iso_net_charge") if target_card else None
    if _target_iso_charge is None:
        import warnings as _w
        _w.warn(
            "target_iso_net_charge missing from target_card — using 0. "
            "Regenerate target_card with generate_target_card.py for correct QM charge.",
            RuntimeWarning, stacklevel=2,
        )
        _target_iso_charge = 0

    # ─── R-15 / R-16 guards (chemistry-true charge topology) ─────────────
    # SciVal v2/v3 verdict (2026-04-19): declared charges from target_card
    # must match the runtime chemistry-true computation from the snapshot
    # PDB. R-16 covers the binder (BinderChargeMismatch); R-15 covers the
    # *magnitude* of the total QM net charge (ChargeDeclarationMismatch).
    # The existing OddElectronError parity check in build_qm_mol() remains
    # an independent safety net (SciVal: "parity 와 magnitude 는 독립 check").
    _binder_charge_computed = _binder_charge_declared
    _target_charge_computed = int(_target_iso_charge)
    _total_charge_computed = _binder_charge_computed + _target_charge_computed
    _charge_topology_diag = None
    _declared_total = _binder_charge_declared + int(_target_iso_charge)
    try:
        if target_card is not None:
            _tc_binder_chain = target_card.get("binder_chain", binder_chain)
            _tc_target_chain = target_card.get("target_chain", "A")
            _tc_contacts = target_card.get("target_contact_residues", [])
            _tc_whole = target_card.get("whole_residue_exceptions", [])
            try:
                (
                    _binder_charge_computed,
                    _target_charge_computed,
                    _total_charge_computed,
                    _charge_topology_diag,
                ) = compute_qm_net_charge_topology(
                    atoms,
                    binder_chain=_tc_binder_chain,
                    target_chain=_tc_target_chain,
                    target_contact_residues=_tc_contacts,
                    whole_residue_exceptions=_tc_whole,
                    pH=7.4,
                )
            except Exception as _chem_err:  # noqa: BLE001 — always emit diagnostic
                _diag(
                    "chemistry_charge_compute_error",
                    snapshot=basename,
                    detail=str(_chem_err)[:200],
                )
                _charge_topology_diag = {"error": str(_chem_err)[:200]}

            _diag(
                "chemistry_charge_computed",
                snapshot=basename,
                pH=7.4,
                binder_chain=_tc_binder_chain,
                target_chain=_tc_target_chain,
                binder_chem=int(_binder_charge_computed),
                target_chem=int(_target_charge_computed),
                total_chem=int(_total_charge_computed),
                cyclic=bool(
                    (_charge_topology_diag or {}).get("binder_diag", {}).get("cyclic", False)
                ),
            )

            # ── R-16: binder hint vs runtime chemistry ─────────────────
            if _binder_charge_computed != _binder_charge_declared:
                _diag(
                    "binder_charge_mismatch",
                    snapshot=basename,
                    binder_card=int(_binder_charge_declared),
                    binder_chem=int(_binder_charge_computed),
                    policy="fail_fast_r14",
                )
                _msg = (
                    f"R-16 BinderChargeMismatch: target_card.binder_net_charge="
                    f"{_binder_charge_declared} but runtime PDB computes "
                    f"{_binder_charge_computed} for chain "
                    f"{_tc_binder_chain}. Per-residue breakdown: "
                    f"{(_charge_topology_diag or {}).get('binder_diag', {}).get('per_residue')}"
                )
                raise BinderChargeMismatch(_msg)

            # ── R-15: declared magnitude vs runtime chemistry ──────────
            _declared_total = _binder_charge_declared + int(_target_iso_charge)
            if _declared_total != _total_charge_computed:
                _diag(
                    "charge_declaration_mismatch",
                    snapshot=basename,
                    declared=int(_declared_total),
                    computed=int(_total_charge_computed),
                    binder_declared=int(_binder_charge_declared),
                    target_declared=int(_target_iso_charge),
                    binder_chem=int(_binder_charge_computed),
                    target_chem=int(_target_charge_computed),
                    policy="fail_fast_r13",
                )
                _msg = (
                    f"R-15 ChargeDeclarationMismatch: declared qm_net_charge="
                    f"{_declared_total} (binder={_binder_charge_declared} + target="
                    f"{_target_iso_charge}) but chemistry-true topology computes "
                    f"{_total_charge_computed} (binder={_binder_charge_computed} + "
                    f"target={_target_charge_computed}). "
                    "Parity-invariant silent pass prevented. "
                    "See SciVal verdict_target_card_charge_rationale_20260419_v2.md."
                )
                raise ChargeDeclarationMismatch(_msg)

            # Equal → verification passed.
            _diag(
                "charge_magnitude_verified",
                snapshot=basename,
                declared=int(_declared_total),
                computed=int(_total_charge_computed),
            )
    except (BinderChargeMismatch, ChargeDeclarationMismatch) as _chem_guard_err:
        _reason = (
            "binder_charge_mismatch_r14"
            if isinstance(_chem_guard_err, BinderChargeMismatch)
            else "charge_declaration_mismatch_r13"
        )
        print(f"  [R-15/R-16] charge guard raised → skipped: {basename}")
        print(f"  [R-15/R-16] detail: {_chem_guard_err}")
        failure_result = {
            "snapshot":               basename,
            "pdb_path":               pdb_path,
            "qm_method":              f"{qm_xc}/{qm_basis}",
            "n_qm_atoms":             len(qm_atoms),
            "n_mm_atoms":             len(mm_atoms),
            "ncaa_element":           ncaa_elem,
            "energy_total_hartree":   None,
            "energy_qm_hartree":      None,
            "energy_qmmm_kcal":       None,
            "interaction_kcal":       None,
            "converged":              False,
            "is_run_mmgbsa":          False,
            "status":                 "FAILED",
            "reason":                 _reason,
            "charge_guard_detail":    str(_chem_guard_err),
            "qm_net_charge_declared": int(_declared_total),
            "qm_net_charge_computed": int(_total_charge_computed),
            "binder_charge_declared": int(_binder_charge_declared),
            "binder_charge_computed": int(_binder_charge_computed),
            "charge_consistency_audit": "failed",
            "charge_topology_diag":   _charge_topology_diag,
            "regime":                 "ranking_only",
        }
        with open(out_json, "w", encoding="utf-8") as f:
            json.dump(failure_result, f, indent=2)
        return failure_result

    # v0.6.x legacy path: ``_binder_charge`` name retained for downstream
    # diagnostic messages below. Use computed value as the source of truth
    # (equal to declared when guards pass; guards raise on mismatch).
    _binder_charge = int(_binder_charge_computed)
    qm_net_charge = _binder_charge + int(_target_iso_charge)

    try:
        mol = build_qm_mol(qm_atoms, charge=qm_net_charge, basis=qm_basis)
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
            "qm_net_charge_declared": int(_declared_total),
            "qm_net_charge_computed": int(_total_charge_computed),
            "binder_charge_declared": int(_binder_charge_declared),
            "binder_charge_computed": int(_binder_charge_computed),
            "charge_consistency_audit": "failed",
            "regime":               "ranking_only",
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
    mf.max_cycle   = 300
    # SciVal 승인 축소 번들 (B2/B3/B5/B6) — _apply_scf_bundle 로 중앙화된 설정
    # 적용. grids.level / conv_tol / direct_scf_tol / init_guess 를 환경변수로 관리.
    _apply_scf_bundle(mf)
    print(f"  [SCF-Bundle] grids={UPDD_DFT_GRIDS_LEVEL} "
          f"conv_tol={UPDD_SCF_CONV_TOL} "
          f"direct_scf_tol={UPDD_SCF_DIRECT_TOL} "
          f"init_guess={UPDD_SCF_INIT_GUESS}")
    _diag(
        "scf_bundle",
        snapshot=basename,
        grids_level=UPDD_DFT_GRIDS_LEVEL,
        conv_tol=UPDD_SCF_CONV_TOL,
        direct_scf_tol=UPDD_SCF_DIRECT_TOL,
        init_guess=UPDD_SCF_INIT_GUESS,
        min_ao_blksize=int(os.environ.get("UPDD_MIN_AO_BLKSIZE", "64")),
        vram_pool_fraction=UPDD_VRAM_POOL_FRACTION,
    )

    # Phase 5a: DF-J/K density fitting (16GB VRAM 실행 가능성 복구).
    # 560-atom QM region direct-SCF → ~266GB VRAM 필요 → 16GB RTX 5070 Ti 불가.
    # density_fit() 는 반드시 to_gpu() / mm_charge() 이전에 호출해야 한다 (decorator chain 순서).
    if use_df:
        mf = mf.density_fit(auxbasis=df_auxbasis)
        print(f"  [DF-J/K] density_fit 활성화 (auxbasis={df_auxbasis})")
        print(f"  [DF-J/K] use_gpu_memory=False (CPU pinned RAM) | "
              f"max_memory={UPDD_DF_JK_MAX_MEMORY_MB} MB (J/K blksize 제한)")
        _apply_df_oom_guard(mf)
        # density_fit wrapper 가 번들 속성을 감싸는 경우에 대비해 재적용 (idempotent).
        _apply_scf_bundle(mf)

    gpu_success = False
    try:
        # 1. 먼저 평범한 계산기를 GPU 로켓으로 변신!
        import gpu4pyscf
        from gpu4pyscf.qmmm import itrf as gpu_qmmm_itrf
        # CuPy 기본 동작: VRAM 전체를 pool로 선점 → 디스플레이·pinned memory 공간 없어짐.
        # to_gpu() 전에 pool 상한을 10GB로 고정해 나머지 ~6GB를 확보한다.
        try:
            import cupy as _cp_pre
            _cp_pre.get_default_memory_pool().free_all_blocks()
            _cp_pre.get_default_pinned_memory_pool().free_all_blocks()
            _vram_free, _vram_total = _cp_pre.cuda.Device(0).mem_info
            # B9: VRAM pool 상한 (default 0.65 = 65%). 환경변수
            # UPDD_VRAM_POOL_FRACTION 으로 override — 0.875 등 상향 시 디스플레이
            # / pinned memory 공간 축소 주의.
            _pool_limit = int(_vram_total * UPDD_VRAM_POOL_FRACTION)
            _cp_pre.get_default_memory_pool().set_limit(size=_pool_limit)
            print(f"  [VRAM] GPU 시도 전: used={(_vram_total-_vram_free)//1024**2}MB "
                  f"free={_vram_free//1024**2}MB total={_vram_total//1024**2}MB "
                  f"pool_limit={_pool_limit//1024**2}MB")
        except Exception:
            pass
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
        mf.max_cycle   = 300
        # SciVal 승인 축소 번들 (B2/B3/B5/B6) — GPU/CPU 경로 method 동일성 유지.
        _apply_scf_bundle(mf)
        # Phase 5a: CPU fallback 에도 동일하게 DF-J/K 적용 (mm_charge 이전).
        # GPU/CPU 경로 간 method 동일성 보장 → interaction_kcal 일관성.
        if use_df:
            mf = mf.density_fit(auxbasis=df_auxbasis)
            # CPU 경로에서도 use_gpu_memory=False / max_memory 를 동일 적용
            # (GPU/CPU 양경로의 DF 객체 상태를 대칭적으로 구성하기 위함).
            _apply_df_oom_guard(mf)
            _apply_scf_bundle(mf)  # density_fit wrapper 대비 재적용 (idempotent)
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
        _oom = any(kw in str(e).lower() for kw in ("cudaerror", "out of memory", "cudamemory"))
        print(f"  [!] GPU kernel 예외: {type(e).__name__}: {str(e)[:200]}")
        if gpu_success and _oom:
            # ── Tier 1 fallback: GPU + CUDA Unified Memory (VRAM 부족분 → RAM 자동 overflow) ──
            # VRAM 한계를 초과할 때 PCIe를 통해 시스템 RAM을 VRAM 연장으로 사용.
            # 계산 결과는 순수 VRAM과 동일; 대역폭 차이(GDDR7 ~600 GB/s vs PCIe ~32 GB/s)로
            # 느려질 수 있으나 디스크 I/O 없음.
            print(f"  [!] GPU VRAM OOM (n_qm={len(qm_atoms)}) — CUDA Unified Memory(VRAM+RAM)로 재시도.")
            energy_total = None
            # GPU mf 객체가 GPU 텐서를 참조하고 있으면 free_all_blocks()로도 해제 불가.
            # UM 시도 전에 mf 명시 삭제 + GC로 VRAM 완전 반납.
            try:
                import cupy as _cp_gc, gc as _gc
                del mf
                _gc.collect()
                _cp_gc.get_default_memory_pool().free_all_blocks()
                _cp_gc.get_default_pinned_memory_pool().free_all_blocks()
                _free_after, _total_after = _cp_gc.cuda.Device(0).mem_info
                print(f"  [VRAM] mf 해제 후: used={(_total_after-_free_after)//1024**2}MB "
                      f"free={_free_after//1024**2}MB")
            except Exception:
                pass
            try:
                import cupy as _cp
                from gpu4pyscf.qmmm import itrf as _gpu_qmmm_retry

                # 실패한 1차 GPU 시도로 단편화된 VRAM 블록을 먼저 해제.
                # free_all_blocks()는 CuPy 풀이 캐시한 미사용 블록을 OS에 반환하여
                # Unified Memory 풀이 연속된 큰 블록을 확보할 수 있게 한다.
                _cp.get_default_memory_pool().free_all_blocks()
                _cp.get_default_pinned_memory_pool().free_all_blocks()

                # VRAM-first managed allocator:
                # malloc_managed 로 할당된 각 블록에 cudaMemAdviseSetPreferredLocation(GPU)
                # 힌트를 부여해 CUDA가 처음부터 VRAM에 배치하도록 지시.
                # VRAM이 부족해질 때만 LRU 페이지가 RAM으로 축출됨 (on-demand eviction).
                # 힌트 없는 순수 malloc_managed 는 첫 접근 시점에 마이그레이션하므로
                # 초기 page-fault overhead 가 크다.
                _GPU_DEVICE = 0
                class _VRAMFirstAllocator:
                    def __init__(self):
                        self._pool = _cp.cuda.MemoryPool(_cp.cuda.malloc_managed)
                    def malloc(self, size):
                        mem = self._pool.malloc(size)
                        if size > 0:
                            try:
                                # SetPreferredLocation: CUDA야, 이 블록은 GPU에 두는 걸 선호해
                                _cp.cuda.runtime.memAdvise(
                                    mem.ptr, size,
                                    _cp.cuda.runtime.cudaMemAdviseSetPreferredLocation,
                                    _GPU_DEVICE,
                                )
                                # SetAccessedBy: GPU가 자주 접근함 → 프리페치 최적화
                                _cp.cuda.runtime.memAdvise(
                                    mem.ptr, size,
                                    _cp.cuda.runtime.cudaMemAdviseSetAccessedBy,
                                    _GPU_DEVICE,
                                )
                            except Exception:
                                pass  # 힌트 실패는 성능 저하일 뿐, 정확도에 무관
                        return mem

                # RAM-first managed allocator:
                # 각 블록의 PreferredLocation 을 cudaCpuDeviceId(=-1, RAM)로 지정하여
                # 페이지가 기본적으로 시스템 RAM 에 상주하도록 한다. GPU 커널이 해당
                # 페이지에 접근하면 CUDA 가 fault-driven 으로 VRAM 에 migrate 해 준다.
                # AccessedBy=_GPU_DEVICE 힌트는 GPU 쪽 page table 을 사전 매핑해
                # page fault 시 VRAM→RAM round-trip 없이 직접 접근 가능하게 한다.
                # RAM-first: 페이지 기본 위치 RAM, GPU 접근 시 on-demand migrate,
                # PCIe 대역폭 bound (RTX 5070 Ti + AM5 PCIe 5.0 x16 환경에서 유의미).
                class _RAMFirstAllocator:
                    """RAM-first Unified Memory allocator.

                    페이지 기본 위치 RAM, GPU 접근 시 on-demand migrate,
                    PCIe 대역폭 bound. VRAM-first 가 OOM 일 때 Working-set
                    전체가 RAM 에 안착할 수 있는 fallback 경로로 사용.
                    """
                    def __init__(self):
                        self._pool = _cp.cuda.MemoryPool(_cp.cuda.malloc_managed)
                    def malloc(self, size):
                        mem = self._pool.malloc(size)
                        if size > 0:
                            try:
                                # SetPreferredLocation(CpuDeviceId): 이 블록은 RAM 에 두어라
                                _cp.cuda.runtime.memAdvise(
                                    mem.ptr, size,
                                    _cp.cuda.runtime.cudaMemAdviseSetPreferredLocation,
                                    _cp.cuda.runtime.cudaCpuDeviceId,
                                )
                                # SetAccessedBy(GPU): GPU 쪽 page table 사전 매핑 →
                                # fault 시 VRAM round-trip 없이 PCIe 직접 접근
                                _cp.cuda.runtime.memAdvise(
                                    mem.ptr, size,
                                    _cp.cuda.runtime.cudaMemAdviseSetAccessedBy,
                                    _GPU_DEVICE,
                                )
                            except Exception:
                                pass  # 힌트 실패는 성능 저하일 뿐, 정확도에 무관
                        return mem

                _vram_first = _VRAMFirstAllocator()
                _ram_first = None  # RAM-first 진입 시에만 인스턴스화 (finally 해제용)
                _cp.cuda.set_allocator(_vram_first.malloc)
                print(f"  [+] CUDA Unified Memory 활성화 (VRAM 우선, 부족분만 RAM으로 overflow)")
                _diag("compute_backend", snapshot=basename, backend="GPU_unified_memory",
                      df_mode=bool(use_df), df_auxbasis=(df_auxbasis if use_df else None))
                mf_um = _make_dft_mf(mol, qm_xc, _max_mem, use_df, df_auxbasis)
                mf_um.diis_space = 4  # 8→4: DIIS 히스토리 절반, VRAM ~1.2 GB 절약
                mf_um = mf_um.to_gpu()
                if len(mm_coords_charges) > 0:
                    mf_um = _gpu_qmmm_retry.add_mm_charges(
                        mf_um, mm_coords_charges[:, :3], mm_coords_charges[:, 3]
                    )
                    print(f"  [+] MM 점전하: {len(mm_coords_charges)}개 포함 (Unified Memory GPU)")
                try:
                    energy_total = mf_um.kernel()
                    mf = mf_um
                    print(f"  QM/MM 에너지 (GPU+RAM): {energy_total:.8f} Hartree")
                except Exception as e_vf:
                    _oom_vf = any(kw in str(e_vf).lower() for kw in ("cudaerror", "out of memory", "cudamemory"))
                    if not _oom_vf:
                        # 비-OOM 예외는 기존 경로(CPU direct fallback)로 위임.
                        raise
                    # ── Tier 1.5 fallback: RAM-first Unified Memory ──────────────
                    # VRAM-first 가 OOM → 작업 세트가 VRAM+가끔 RAM 으로 불충분.
                    # RAM-first 는 대량 텐서를 RAM 에 상주시키고 GPU 는 PCIe 5.0 을
                    # 통해 필요한 페이지만 pull 한다. 수치 결과는 VRAM-first 와 동일
                    # (Managed memory semantics — placement 힌트만 다름).
                    print(f"  [!] VRAM-first 도 OOM — RAM 우선 Unified Memory 재시도 (PCIe 5.0 page migration).")
                    # VRAM-first pool 해제 → RAM-first pool 이 큰 연속 블록 확보 가능
                    try:
                        del mf_um
                    except Exception:
                        pass
                    try:
                        _vram_first._pool.free_all_blocks()
                    except Exception:
                        pass
                    _cp.get_default_memory_pool().free_all_blocks()
                    _cp.get_default_pinned_memory_pool().free_all_blocks()

                    _ram_first = _RAMFirstAllocator()
                    _cp.cuda.set_allocator(_ram_first.malloc)
                    _diag("compute_backend", snapshot=basename,
                          backend="GPU_unified_memory_ram_first",
                          df_mode=bool(use_df),
                          df_auxbasis=(df_auxbasis if use_df else None))
                    mf_rf = _make_dft_mf(mol, qm_xc, _max_mem, use_df, df_auxbasis)
                    mf_rf.diis_space = 4  # VRAM-first 와 동일 설정 유지
                    mf_rf = mf_rf.to_gpu()
                    if len(mm_coords_charges) > 0:
                        mf_rf = _gpu_qmmm_retry.add_mm_charges(
                            mf_rf, mm_coords_charges[:, :3], mm_coords_charges[:, 3]
                        )
                        print(f"  [+] MM 점전하: {len(mm_coords_charges)}개 포함 (RAM-first Unified Memory GPU)")
                    try:
                        energy_total = mf_rf.kernel()
                        mf = mf_rf
                        print(f"  QM/MM 에너지 (GPU, RAM-first): {energy_total:.8f} Hartree")
                    except Exception as e_rf:
                        _oom_rf = any(kw in str(e_rf).lower() for kw in ("cudaerror", "out of memory", "cudamemory"))
                        if _oom_rf:
                            print(f"  [!] RAM-first 도 OOM — CPU direct-SCF 로 최종 재시도.")
                        else:
                            print(f"  [!] RAM-first Unified Memory 실패: {e_rf} — CPU direct-SCF 로 최종 재시도.")
                        # RAM-first mf 및 pool 정리 후 CPU direct 경로로 fall-through
                        try:
                            del mf_rf
                        except Exception:
                            pass
                        try:
                            _ram_first._pool.free_all_blocks()
                        except Exception:
                            pass
                        # Re-raise 하여 outer CPU direct fallback 로 넘긴다.
                        raise
            except Exception as e_um:
                _oom2 = any(kw in str(e_um).lower() for kw in ("cudaerror", "out of memory", "cudamemory"))
                if _oom2:
                    print(f"  [!] Unified Memory도 OOM — CPU direct-SCF로 최종 재시도.")
                else:
                    print(f"  [!] Unified Memory GPU 실패: {e_um} — CPU direct-SCF로 최종 재시도.")
                # ── Tier 2 fallback: CPU direct-SCF (DF 없음, 디스크 I/O 없음) ──
                _diag("compute_backend", snapshot=basename, backend="CPU_direct_final_fallback",
                      df_mode=False, df_auxbasis=None)
                mf_cpu = dft.RKS(mol)
                mf_cpu.xc         = qm_xc
                mf_cpu.verbose    = 4
                mf_cpu.max_memory = _max_mem
                mf_cpu.max_cycle  = 300
                mf_cpu.direct_scf = True
                # SciVal 승인 축소 번들 (B2/B3/B5/B6) — 최종 CPU direct-SCF 에도 동일 적용.
                # GPU 경로와 method 비교 일관성 (ranking preservation).
                _apply_scf_bundle(mf_cpu)
                if len(mm_coords_charges) > 0:
                    mf_cpu = qmmm_itrf.mm_charge(
                        mf_cpu, mm_coords_charges[:, :3], mm_coords_charges[:, 3]
                    )
                    print(f"  [+] MM 점전하: {len(mm_coords_charges)}개 포함 (CPU direct)")
                # gpu4pyscf 임포트 시 PySCF 내부까지 cupy allocator로 패치됨.
                # CPU kernel 실행 전 device/pinned allocator 완전 해제하지 않으면
                # grid partition 등 CPU 코드에서도 pinned_memory OOM 발생.
                try:
                    import cupy as _cp_detach
                    _cp_detach.get_default_memory_pool().free_all_blocks()
                    _cp_detach.get_default_pinned_memory_pool().free_all_blocks()
                    _cp_detach.cuda.set_allocator(None)
                    _cp_detach.cuda.set_pinned_memory_allocator(None)
                except Exception:
                    pass
                try:
                    energy_total = mf_cpu.kernel()
                    mf = mf_cpu
                    print(f"  QM/MM 에너지 (CPU-direct): {energy_total:.8f} Hartree")
                except Exception as e_cpu:
                    print(f"  [!] CPU direct-SCF 실패: {e_cpu}")
                    energy_total = None
            finally:
                # _VRAMFirstAllocator / _RAMFirstAllocator pool 이 managed memory 를
                # VRAM/RAM 에 잡아두는 것을 방지. 스냅샷 간 누적 방지를 위해 pool 명시
                # 해제 후 allocator 복원. 두 pool 모두 존재할 수 있으므로 각각 독립 해제.
                try:
                    _vram_first._pool.free_all_blocks()
                    del _vram_first
                except Exception:
                    pass
                try:
                    if _ram_first is not None:
                        _ram_first._pool.free_all_blocks()
                        del _ram_first
                except Exception:
                    pass
                try:
                    import cupy as _cp2, gc
                    _cp2.get_default_memory_pool().free_all_blocks()
                    _cp2.get_default_pinned_memory_pool().free_all_blocks()
                    _cp2.cuda.set_allocator(_cp2.get_default_memory_pool().malloc)
                    _cp2.cuda.set_pinned_memory_allocator(
                        _cp2.get_default_pinned_memory_pool().malloc)
                    gc.collect()
                except Exception:
                    pass
        else:
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
        mf_qm.max_memory = _max_mem
        # SciVal 승인 축소 번들 (B2/B3/B5/B6) — main SCF 와 동일 tuning 으로
        # interaction_kcal = (E_QMMM - E_QM) 차분의 method 일관성 유지.
        _apply_scf_bundle(mf_qm)
        # Phase 5a: QM-only 도 동일 DF-J/K 적용 (mf 와 method 일치 보장).
        # interaction_kcal = (E_QMMM - E_QM) 의 일관된 차분 → DF 오차 cancel.
        if use_df:
            mf_qm = mf_qm.density_fit(auxbasis=df_auxbasis)
            # QM-only SCF 에도 동일 OOM 가드 (main SCF 와 DF 메모리 정책 일치).
            _apply_df_oom_guard(mf_qm)
            _apply_scf_bundle(mf_qm)  # density_fit wrapper 대비 재적용 (idempotent)

        # 2라운드는 점전하 짐칸이 없으므로 그냥 변신만 시키면 됩니다!
        if gpu_success:
            try:
                mf_qm = mf_qm.to_gpu()
                print("  🚀 [GPU 가속 - QM 단독계산]")
            except Exception as e:
                print(f"  [!] 2라운드 GPU 전환 실패 (사유: {e})")

        # 이전에는 여기서 mf_qm.grids.level=3 하드 오버라이드가 있었으나
        # SciVal 축소 번들(B2)로 중앙화 (UPDD_DFT_GRIDS_LEVEL 환경변수 참조).
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
                # Sub-issue #3 (Option A) + R-16: use chemistry-true binder net
                # charge from snapshot PDB rather than hardcoding 0. The previous
                # comment ("even protons → standard whole-peptide charge 0") was
                # misleading — proton parity is independent of net charge (charge
                # ±2 is also even parity). Binder net charge is derived from
                # residue identity + atom-name pattern (ARG HH11-HH22, ASP OD
                # no-HD, GLU OE no-HE, LYS HZ*, HIS HD1/HE2 tautomer, termini
                # detection for cyclic vs linear). See SciVal v2 §1.1.
                binder_charge, _binder_subsys_diag = compute_binder_chem_charge(
                    binder_qm, binder_chain=tc_binder_chain_p4,
                )
                _diag(
                    "binder_chemistry_verified",
                    snapshot=basename,
                    binder_card=int(_binder_charge_declared),
                    binder_chem=int(binder_charge),
                    cyclic=bool(_binder_subsys_diag.get("cyclic", False)),
                    subsystem="binder_iso",
                )

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
    # R-11 (CLAUDE.md Rev. 2): "regime" 필드를 명시하여 이 결과가 ranking-only
    # 해석용인지 absolute-accuracy 해석용인지 downstream 소비자가 구분하게 한다.
    # v0.7 calibration 진입 전까지 "ranking_only" 고정. calibration 이후
    # 결과는 "absolute_calibrated" 로 교체될 예정 (calibrate_results.py 에서 set).
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
        # R-15 / R-16 audit fields (2026-04-19 Stage 1) — always present when
        # topology mode resolves charge via compute_qm_net_charge_topology().
        # "passed" means declared == computed at runtime (magnitude + binder
        # hint). Guards raise before reaching this path for mismatches.
        "qm_net_charge_declared": int(_declared_total),
        "qm_net_charge_computed": int(_total_charge_computed),
        "binder_charge_declared": int(_binder_charge_declared),
        "binder_charge_computed": int(_binder_charge_computed),
        "charge_consistency_audit": "passed",
        # R-11: Regime tag — v0.7 calibration 진입 전까지 "ranking_only" 고정.
        "regime":                 "ranking_only",
        "regime_note":            "Absolute binding energy 의 정량 해석 금지; design 간 상대 순위만 validated.",
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
    # v0.6.x — Per-snapshot subprocess isolation 지원.
    # UPDD.py 가 snapshot 단위로 subprocess 를 spawn 할 때 사용한다.
    # 단일 stem (예: "snap02") 을 지정하면 해당 stem 이 파일명에 포함된 PDB 만 처리.
    parser.add_argument("--snapshot-filter", dest="snapshot_filter", default=None,
                        help="특정 snapshot stem(예: snap02)만 처리 — per-snapshot subprocess isolation용")
    parser.add_argument("--binder_chain", default="B", help="Binder chain ID (기본 B)")
    # Phase 1: target_card 기반 topology mode 활성화. 지정 시 --mode 는 무시되고
    # target_cards/{target_id}.json 의 QM partitioning 정책을 따른다 (snapshot-invariant QM region).
    parser.add_argument("--target-id", dest="target_id", type=str, default=None,
                        help="Target card ID (e.g., 6WGN). Activates topology mode.")
    # Phase 5a: DF-J/K density fitting 제어. 기본 활성 (16GB VRAM 필수).
    # `--no-df` 는 소규모 QM region (예: <100 atoms) 에서 direct-SCF 비교 검증용.
    # v0.6.x: `--df-mode` 도입. `--no-df` 는 `--df-mode off` 의 alias 로 보존.
    parser.add_argument("--no-df", dest="use_df", action="store_false",
                        help="DF 비활성 (direct-SCF). --df-mode off 와 동일.")
    parser.set_defaults(use_df=True)
    parser.add_argument("--df-mode", dest="df_mode",
                        choices=["auto", "on", "off"], default="auto",
                        help="DF 모드: auto (VRAM 감지, default), on (강제), off (강제). "
                             "환경변수 UPDD_DF_AUTO_VRAM_THRESHOLD_GB 로 auto 임계값 조정 "
                             "(default 20 GiB). --no-df 는 off 의 alias. "
                             "환경변수 UPDD_DF_JK_MAX_MEMORY_MB (default 4000) 으로 "
                             "DF J/K 버퍼 조절 가능. VRAM OOM 시 하향, 여유 시 상향. "
                             "SCF 수렴/성능 튜닝 (SciVal 승인 축소 번들): "
                             "UPDD_SCF_CONV_TOL (default 1e-7), "
                             "UPDD_DFT_GRIDS_LEVEL (default 2), "
                             "UPDD_SCF_INIT_GUESS (default minao; "
                             "'huckel' 은 opt-in — near-degenerate HOMO-LUMO 에서 발산 위험), "
                             "UPDD_SCF_DIRECT_TOL (default 1e-12), "
                             "UPDD_MIN_AO_BLKSIZE (default 64, gpu4pyscf 16GB 경계 버그 우회), "
                             "UPDD_VRAM_POOL_FRACTION (default 0.65).")
    parser.add_argument("--df-auxbasis", dest="df_auxbasis", default="def2-universal-jfit",
                        help="DF 보조 기저함수 (기본: def2-universal-jfit, def2-SVP orbital basis 와 일치). "
                             "def2 family 와 매칭 안 되는 구형 aux (weigend 등) 는 명시 지정 시에만 사용.")
    args = parser.parse_args()

    # ==========================================
    # DF 결정 로직 (v0.6.x — auto/on/off)
    # ==========================================
    # Precedence: `--no-df` (explicit) > `--df-mode`.
    # argparse 기본값(`use_df=True`) 과 `--no-df`(store_false) 만으로는 사용자가
    # 명시적으로 `--no-df` 를 썼는지 구분 못함 → `sys.argv` 직접 스캔.
    _no_df_explicit = ("--no-df" in sys.argv)
    if _no_df_explicit:
        final_use_df = False
        _df_decision_reason = "explicit --no-df"
    elif args.df_mode == "on":
        final_use_df = True
        _df_decision_reason = "--df-mode on (forced)"
    elif args.df_mode == "off":
        final_use_df = False
        _df_decision_reason = "--df-mode off (forced)"
    else:  # auto
        final_use_df, _df_decision_reason = _decide_df_auto()

    args.use_df = final_use_df
    _threshold_gb = float(os.environ.get("UPDD_DF_AUTO_VRAM_THRESHOLD_GB", "20"))
    print(f"[DF-Mode] use_df={final_use_df} | reason: {_df_decision_reason}")
    # Main 진입부이므로 `_diag` 대신 `[DIAG]` JSON-line 을 직접 print.
    # snapshot 루프 이전 단계라 snapshot 필드는 없다.
    print("[DIAG] " + json.dumps({
        "event": "df_decision",
        "use_df": bool(final_use_df),
        "reason": _df_decision_reason,
        "df_mode_cli": args.df_mode,
        "no_df_explicit": _no_df_explicit,
        "threshold_gb": _threshold_gb,
    }, ensure_ascii=False), flush=True)

    os.makedirs(args.outputdir, exist_ok=True)

    all_pdbs = sorted(glob.glob(os.path.join(args.snapdir, "*.pdb")))
    if args.filter:
        pdb_files = [f for f in all_pdbs if os.path.basename(f).startswith(args.filter)]
        print(f"  [Filter] '{args.filter}'로 시작하는 스냅샷 {len(pdb_files)}개를 찾았습니다.")
    else:
        pdb_files = all_pdbs

    # v0.6.x — Per-snapshot subprocess isolation: stem 기반 post-filter.
    # SciVal Q6: legacy mode(--target-id 미지정) 에서도 동작하지만 topology mode 와의
    # 조합을 권장. n_qm 이 snapshot 간 동일해야 snapshot subprocess 의 의미가 살아난다.
    if args.snapshot_filter:
        before = len(pdb_files)
        _stem = args.snapshot_filter
        pdb_files = [
            f for f in pdb_files
            if f"_{_stem}_" in os.path.basename(f)
            or os.path.basename(f).endswith(f"_{_stem}.pdb")
        ]
        print(f"  [SnapFilter] stem='{_stem}' 매칭: {before} → {len(pdb_files)}개")
        if not args.target_id:
            print("  [WARN] --snapshot-filter는 topology mode(--target-id) 와 함께 사용 권장. "
                  "legacy mode 에서는 snapshot 간 n_qm 이 변동하여 isolation 효과가 제한적.")

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
