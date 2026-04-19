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
import os
import json
import time
import glob
import tempfile


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


# ==========================================
# [v0.5] Per-design MD Status SSOT
# ==========================================

_MD_STATUS_FILENAME = "_md_status.json"

EXPLOSION_MARKER_PATTERNS = [
    "_EXPLODED_dt2fs.log", "_EXPLODED_dt1fs.log",
    "_EXPLODED_MD.log", "_EXPLODED_NVT.log", "_EXPLODED_COLD_NVT.log",
    "_EXPLODED_MIN.log", "_EXPLODED_CYCLIC_MIN.log", "_EXPLODED_CYCLIC_RELAX.log",
]

PARTIAL_PLACEHOLDERS = ["_partial_md.pdb", "_partial_nvt.pdb"]


def load_md_status(project_dir):
    # type: (str) -> Optional[dict]
    """Read outputs/{project}/_md_status.json and return as dict.

    Returns None if the file does not exist or is corrupted.
    """
    path = os.path.join(project_dir, _MD_STATUS_FILENAME)
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError, OSError):
        return None


def write_md_status(project_dir, data):
    # type: (str, dict) -> bool
    """Atomic write of _md_status.json using tempfile + os.replace.

    Returns True on success, False on failure.
    """
    path = os.path.join(project_dir, _MD_STATUS_FILENAME)
    data["last_updated"] = time.time()
    try:
        fd, tmp_path = tempfile.mkstemp(
            dir=project_dir, suffix=".tmp", prefix="_md_status_"
        )
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        os.replace(tmp_path, path)
        return True
    except (IOError, OSError):
        if "tmp_path" in locals() and os.path.exists(tmp_path):
            try:
                os.remove(tmp_path)
            except OSError:
                pass
        return False


def get_design_tier(project_dir, design_stem):
    # type: (str, str) -> str
    """Return the tier string for a design, or 'UNKNOWN' if unavailable."""
    md_status = load_md_status(project_dir)
    if md_status is None:
        return "UNKNOWN"
    entry = md_status.get("designs", {}).get(design_stem, {})
    return entry.get("tier", "UNKNOWN")


def classify_md_outcome(md_dir, design_stem):
    # type: (str, str) -> dict
    """Analyze mdresult/ file patterns for a single design and return tier classification.

    Used by both the in-pipeline producer (Phase 2) and backfill tool (Phase 7).

    Args:
        md_dir: Path to mdresult/ directory.
        design_stem: Design base name (e.g. 'design_w4_0_s2').

    Returns:
        dict with keys: tier, ranking_include, md_stability_tier, passes,
        stale_markers, extract_snapshots_action.
    """
    final_pdb = os.path.join(md_dir, design_stem + "_final.pdb")
    has_final = os.path.exists(final_pdb)

    found_markers = []  # type: List[str]
    for pat in EXPLOSION_MARKER_PATTERNS:
        marker_path = os.path.join(md_dir, design_stem + pat)
        if os.path.exists(marker_path):
            found_markers.append(pat)

    found_partials = []  # type: List[str]
    for pat in PARTIAL_PLACEHOLDERS:
        pp = os.path.join(md_dir, design_stem + pat)
        if os.path.exists(pp):
            found_partials.append(pat)

    has_dt2fs_marker = "_EXPLODED_dt2fs.log" in found_markers
    has_dt1fs_marker = "_EXPLODED_dt1fs.log" in found_markers
    has_any_marker = len(found_markers) > 0

    passes = []  # type: List[dict]
    tier = "FAIL"
    extract_action = "skip"

    if not has_any_marker and has_final:
        tier = "SUCCESS"
        extract_action = "normal_trajectory"
        passes.append({
            "pass": 1,
            "dt_fs": 2.0,
            "status": "SUCCESS",
            "completion_ratio": 1.0,
        })

    elif has_any_marker and has_final:
        tier = "MARGINAL"
        extract_action = "normal_trajectory"

        pass1_status = "FAIL_NAN_EARLY"
        if has_dt2fs_marker:
            pass1_status = "FAIL_NAN_EARLY"
        elif "_EXPLODED_NVT.log" in found_markers:
            pass1_status = "FAIL_NVT"
        elif "_EXPLODED_MIN.log" in found_markers or \
             "_EXPLODED_CYCLIC_MIN.log" in found_markers:
            pass1_status = "FAIL_MIN"
        else:
            pass1_status = "FAIL_NAN_LATE"

        passes.append({
            "pass": 1,
            "dt_fs": 2.0,
            "status": pass1_status,
            "exploded_marker": found_markers[0] if found_markers else None,
            "partial_placeholder": found_partials[0] if found_partials else None,
        })
        passes.append({
            "pass": 2,
            "dt_fs": 1.0,
            "status": "SUCCESS",
            "completion_ratio": 1.0,
            "final_pdb": os.path.basename(final_pdb),
        })

    elif has_any_marker and not has_final:
        if found_partials:
            tier = "PARTIAL"
            extract_action = "trim_exploded"
            passes.append({
                "pass": 1,
                "dt_fs": 2.0,
                "status": "PARTIAL_SUCCESS",
                "exploded_marker": found_markers[0] if found_markers else None,
                "partial_placeholder": found_partials[0] if found_partials else None,
            })
        else:
            tier = "FAIL"
            extract_action = "skip"
            passes.append({
                "pass": 1,
                "dt_fs": 2.0,
                "status": "FAIL_NAN_EARLY",
                "exploded_marker": found_markers[0] if found_markers else None,
            })
            if has_dt1fs_marker:
                passes.append({
                    "pass": 2,
                    "dt_fs": 1.0,
                    "status": "FAIL_NAN_EARLY",
                    "exploded_marker": "_EXPLODED_dt1fs.log",
                })

    else:
        if has_final:
            tier = "SUCCESS"
            extract_action = "normal_trajectory"
            passes.append({
                "pass": 1,
                "dt_fs": 2.0,
                "status": "SUCCESS",
                "completion_ratio": 1.0,
            })
        else:
            tier = "FAIL"
            extract_action = "skip"

    ranking_include = (tier == "SUCCESS")

    return {
        "tier": tier,
        "ranking_include": ranking_include,
        "md_stability_tier": tier,
        "passes": passes,
        "stale_markers": found_markers,
        "extract_snapshots_action": extract_action,
    }


# ==========================================
# [v0.6.1] Target Card loader — Topology-based QM region selection
# ==========================================
def load_target_card(target_id, cards_dir=None):
    # type: (str, Optional[str]) -> dict
    """Load and validate a target_card JSON for topology-based QM/MM partitioning.

    Target card schema encodes the pre-curated QM region selection policy for a
    given target (chain assignments, contact residues, cofactor treatment,
    fragment-boundary rules). Pipeline callers pass the returned dict directly to
    ``run_qmmm.partition_qmmm(..., mode="topology", target_card=<dict>)``.

    Supported schema versions:
        - 0.6.0 — initial topology mode (Phase 1)
        - 0.6.1 — whole_residue_exceptions added
        - 0.6.2 — multi-snapshot occupancy filtering
        - 0.6.3 — [v0.6.5 Phase 4] target_iso_net_charge / binder_net_charge added
                   for supermolecular qm_int_kcal_frozen decomposition.
        - 0.6.4 — Stage 2 (2026-04-19): numbering_convention object, auto-generated
                   target_iso_charge_rationale, target_residues_mode=whole as default
                   (Option B3), max_n_qm 500→600. whole_residue_exceptions removed
                   (redundant under whole mode; kept optional for 0.6.3 backward compat).

    Args:
        target_id: Target identifier (e.g. "6WGN"). Must match the JSON filename
            stem (case-sensitive).
        cards_dir: Optional override for the target_cards directory. Defaults to
            ``<utils>/../target_cards`` when not provided.

    Returns:
        dict: Parsed target card with ``schema_version`` validated.

    Raises:
        FileNotFoundError: Card JSON not located at the resolved path.
        ValueError: ``schema_version`` missing or not in the supported set.
    """
    if cards_dir is None:
        cards_dir = os.path.join(os.path.dirname(__file__), "..", "target_cards")
    path = os.path.join(cards_dir, "{}.json".format(target_id))
    if not os.path.exists(path):
        raise FileNotFoundError("target_card not found: {}".format(path))
    with open(path, "r", encoding="utf-8") as f:
        card = json.load(f)
    supported = {"0.6.0", "0.6.1", "0.6.2", "0.6.3", "0.6.4"}
    if card.get("schema_version") not in supported:
        raise ValueError(
            "Unsupported target_card schema_version: {} (supported: {})".format(
                card.get("schema_version"), sorted(supported)
            )
        )

    # [v0.6.5 Phase 4] Phase 4 (qm_int_kcal_frozen) 활성화 조건: target_iso_net_charge.
    # 누락 시 qm_int_kcal_frozen 은 계산되지 않고 None 으로 기록된다 (backward-compat).
    if "target_iso_net_charge" not in card:
        import warnings
        warnings.warn(
            "target_card '{}' missing 'target_iso_net_charge' — qm_int_kcal_frozen "
            "will be skipped. Regenerate with generate_target_card.py (v0.6.5+) "
            "or add the field manually.".format(target_id),
            stacklevel=2,
        )

    return card
