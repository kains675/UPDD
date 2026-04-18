#!/usr/bin/env python
"""
rank_results_qmmm.py
--------------------
AF2 pLDDT + QM/MM 에너지 + MM-GBSA ΔG를 통합하여
최종 후보 펩타이드를 랭킹합니다.

점수 구성 (기본 가중치):
  - pLDDT         : 30%  (높을수록 구조 안정성↑)
  - QM/MM 상호작용: 40%  (낮을수록 전자 수준 결합↑)
  - MM-GBSA ΔG    : 30%  (낮을수록 결합 에너지↑)
사용 환경: qmmm (conda)
"""

import os
import sys
import glob
import json
import argparse
import csv
import numpy as np
from typing import Optional, Tuple

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    from utils_common import load_md_status  # noqa: E402
except ImportError:
    def load_md_status(project_dir):  # type: ignore
        return None


# ==========================================
# [DIAG] 구조화 로그 헬퍼 — Path 에이전트 자동 진단용
# ==========================================
# Path 가 `grep "\[DIAG\]"` 한 번으로 ranking 단계의 모든 결정점을 수집할 수 있도록
# JSON-line 형식으로 출력한다. 점수/순위 계산은 변경하지 않으며 순수 관찰성만 추가한다.
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
# 데이터 로드 함수들
# ==========================================

def load_af2_scores(af2_dir):
    """
    AF2 결과 디렉토리에서 pLDDT 점수 로드
    ColabFold 출력 JSON (*_scores*.json) 파싱
    Returns: dict {design_name: pLDDT}
    """
    scores = {}

    # ColabFold JSON 파일 탐색
    json_files = glob.glob(os.path.join(af2_dir, "*_scores*.json"))
    if not json_files:
        json_files = glob.glob(os.path.join(af2_dir, "**", "*scores*.json"), recursive=True)

    for jf in json_files:
        try:
            with open(jf, encoding="utf-8") as f:
                data = json.load(f)

            # pLDDT 추출 (ColabFold 포맷)
            plddt_list = data.get("plddt", [])
            if not plddt_list:
                plddt_list = data.get("mean_plddt", [])

            if isinstance(plddt_list, list) and plddt_list:
                plddt_val = float(np.mean(plddt_list))
            elif isinstance(plddt_list, (int, float)):
                plddt_val = float(plddt_list)
            else:
                continue

            # 바인더 체인(B)만의 pLDDT 추출 시도
            ptm  = data.get("ptm",  None)
            iptm = data.get("iptm", None)

            design_name = os.path.basename(jf).split("_scores")[0]
            scores[design_name] = {
                "plddt": plddt_val,
                "ptm":   float(ptm)  if ptm  is not None else None,
                "iptm":  float(iptm) if iptm is not None else None,
            }
        except Exception as e:
            print(f"  [!] AF2 JSON 파싱 실패 ({os.path.basename(jf)}): {e}")

    # af2_ranking.csv fallback
    if not scores:
        csv_path = os.path.join(af2_dir, "af2_ranking.csv")
        if os.path.exists(csv_path):
            with open(csv_path, encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    name = row.get("design", row.get("name", ""))
                    plddt = float(row.get("plddt", row.get("mean_plddt", 0)))
                    scores[name] = {"plddt": plddt, "ptm": None, "iptm": None}

    print(f"  AF2 pLDDT 로드: {len(scores)}개 설계")
    return scores


def load_qmmm_scores(qmmm_dir):
    """
    QM/MM 결과에서 에너지 로드
    Returns: dict {snapshot_name: {energy_total, interaction_kcal, qm_int_kcal_frozen}}

    qm_int_kcal_frozen 필드 추가:
        Phase 4 (supermolecular Morokuma frozen-geometry decomposition) 결과를
        함께 로드한다. 구 JSON 에는 필드가 없으므로 None 으로 처리.
    """
    scores = {}
    summary = os.path.join(qmmm_dir, "qmmm_summary.json")

    if os.path.exists(summary):
        with open(summary, encoding="utf-8") as f:
            data = json.load(f)
        for r in data.get("results", []):
            name = r.get("snapshot", "")
            scores[name] = {
                "qmmm_energy_kcal":   r.get("energy_qmmm_kcal"),
                "qmmm_interact_kcal": r.get("interaction_kcal"),
                "qm_int_kcal_frozen": r.get("qm_int_kcal_frozen"),
                "converged":          r.get("converged", False),
            }
    else:
        # 개별 JSON 폴백 — run_qmmm.py 출력 패턴: {snap}_qmmm_{mode}.json
        for jf in glob.glob(os.path.join(qmmm_dir, "*_qmmm_*.json")):
            with open(jf, encoding="utf-8") as f:
                r = json.load(f)
            import re as _re
            name = r.get("snapshot") or _re.sub(r"_qmmm_\w+\.json$", "", os.path.basename(jf))
            scores[name] = {
                "qmmm_energy_kcal":   r.get("energy_qmmm_kcal"),
                "qmmm_interact_kcal": r.get("interaction_kcal"),
                "qm_int_kcal_frozen": r.get("qm_int_kcal_frozen"),
                "converged":          r.get("converged", False),
            }

    print(f"  QM/MM 결과 로드: {len(scores)}개 스냅샷")

    # [Phase 6 DIAG] QM/MM 로드 통계 — frozen 가용성 / 실패 수 추적
    n_with_frozen = sum(
        1 for v in scores.values() if v.get("qm_int_kcal_frozen") is not None
    )
    n_failed = sum(1 for v in scores.values() if not v.get("converged", False))
    _diag(
        "rank_load",
        source="qmmm",
        n_total=len(scores),
        n_with_frozen=n_with_frozen,
        n_failed=n_failed,
    )

    return scores


def load_mmgbsa_scores(mmgbsa_dir):
    """
    MM-GBSA 결과에서 ΔG 로드. [v4 S-1] interaction_entropy_per_design 가 summary 에
    있으면 design 별 IE 보정 ΔG 를 함께 부착하여, compute_ranking 이 보정 ΔG 를
    랭킹에 사용할 수 있도록 한다.

    Returns: dict {snapshot_name: {delta_g_kcal, delta_g_corrected_kcal, std_delta_g_kcal, minus_TdS_kcal, favorable}}
    """
    scores = {}
    summary = os.path.join(mmgbsa_dir, "mmgbsa_summary.json")

    if os.path.exists(summary):
        with open(summary, encoding="utf-8") as f:
            data = json.load(f)
        ie_per_design = data.get("interaction_entropy_per_design", {})

        for r in data.get("results", []):
            name = r.get("snapshot", "")
            # snapshot 이름의 "_snap" 이전 prefix = design key
            design_key = name.split("_snap")[0] if "_snap" in name else name
            ie = ie_per_design.get(design_key, {})
            scores[name] = {
                "delta_g_kcal":           r.get("delta_g_kcal"),
                "favorable":              r.get("favorable", False),
                "delta_g_corrected_kcal": ie.get("delta_g_corrected_kcal"),
                "std_delta_g_kcal":       ie.get("std_delta_g_kcal"),
                "minus_TdS_kcal":         ie.get("minus_TdS_kcal"),
            }
    else:
        for jf in glob.glob(os.path.join(mmgbsa_dir, "*_mmgbsa.json")):
            with open(jf, encoding="utf-8") as f:
                r = json.load(f)
            name = r.get("snapshot", os.path.basename(jf).replace("_mmgbsa.json",""))
            scores[name] = {
                "delta_g_kcal": r.get("delta_g_kcal"),
                "favorable":    r.get("favorable", False),
                "delta_g_corrected_kcal": None,
                "std_delta_g_kcal":       None,
                "minus_TdS_kcal":         None,
            }

    print(f"  MM-GBSA 결과 로드: {len(scores)}개 스냅샷")
    return scores


# ==========================================
# 이름 매칭 유틸
# ==========================================
def fuzzy_match(name, candidates):
    """
    설계 이름을 스냅샷 이름과 매칭
    design_0 ↔ design_0_snap01_f1234 혹은 design_0_NMA_snap01 등과 연결
    """
    if name in candidates:
        return name
    for cand in candidates:
        if cand.startswith(name + "_snap") or cand == name:
            return cand
    
    # 공통 prefix (design_0 등)
    for cand in candidates:
        if "_snap" in cand:
            cand_prefix = cand.split("_snap")[0]
            # cand_prefix usually looks like "design_0_NMA"
            if cand_prefix.startswith(name) or name.startswith(cand_prefix):
                return cand
    return None


def fuzzy_match_all(name, candidates):
    """[#29 v0.2] 설계 이름에 매칭되는 **모든** 스냅샷을 리스트로 반환.

    fuzzy_match 는 첫 매치 1개만 반환하므로 디자인당 스냅샷 E_int 평균 랭킹을
    위해 다중 매칭 버전을 추가했다. 매칭 규칙은 fuzzy_match 와 동일:
        - 정확히 이름이 일치하거나
        - "{name}_snap" 으로 시작하거나
        - "_snap" 이전 prefix 가 서로 startswith 관계
    """
    matches = []
    # 정확 일치
    if name in candidates:
        matches.append(name)
    # _snap 접두 일치
    for cand in candidates:
        if cand == name:
            continue
        if cand.startswith(name + "_snap"):
            if cand not in matches:
                matches.append(cand)
    # prefix 양방향 매칭
    for cand in candidates:
        if cand in matches or cand == name:
            continue
        if "_snap" in cand:
            cand_prefix = cand.split("_snap")[0]
            if cand_prefix.startswith(name) or name.startswith(cand_prefix):
                matches.append(cand)
    return matches


# ==========================================
# 정규화
# ==========================================
def robust_percentile_normalize(values, invert=False, lo_pct=5.0, hi_pct=95.0):
    """[v4 S-4] 백분위(percentile) 기반 robust 정규화 — 이상치 방어용.

    Min-Max 정규화는 단일 outlier (극값 1개) 가 전체 범위를 잡아먹어 나머지 데이터의
    동적 범위를 압축한다. 본 함수는 [lo_pct, hi_pct] 구간을 [0, 1] 로 매핑하고
    구간 밖의 값은 clamp 한다 (이상치는 0 또는 1 로 saturating).

    Args:
        values:  정규화 대상 (NaN 포함 가능)
        invert:  True 이면 낮을수록 우수한 지표를 반전
        lo_pct:  하단 백분위 (기본 5)
        hi_pct:  상단 백분위 (기본 95)

    Returns:
        np.ndarray: 정규화된 값 [0, 1]. NaN 입력은 NaN 출력으로 보존.
    """
    arr = np.array(values, dtype=float)
    if np.all(np.isnan(arr)):
        return arr
    lo = float(np.nanpercentile(arr, lo_pct))
    hi = float(np.nanpercentile(arr, hi_pct))
    if hi == lo:
        return np.where(np.isnan(arr), np.nan, 0.5)
    norm = (arr - lo) / (hi - lo)
    norm = np.clip(norm, 0.0, 1.0)
    norm = np.where(np.isnan(arr), np.nan, norm)
    return (1.0 - norm) if invert else norm


def minmax_normalize(values, invert=False):
    """Min-Max 정규화 → [0, 1] 범위로 매핑한다.

    Args:
        values: 정규화 대상 수치 리스트. NaN 항목을 포함할 수 있다.
        invert: True 이면 낮을수록 우수한 지표(에너지, ΔG)를 반전하여
            "높을수록 우수" 척도로 변환한다.

    Returns:
        np.ndarray: shape == len(values). 모든 값이 동일하면 0.5 로 채워진
        균일 배열을 반환한다.

    [v3 7-2] NaN 처리 정책:
        - 정규화 단계: np.nanmin / np.nanmax 가 사용되므로 NaN 항목은
          min/max 계산에서 자연스럽게 제외된다.
        - 출력 단계: NaN 입력은 NaN 출력으로 그대로 전파된다 (수치 보존).
        - 호출부 (compute_ranking) 는 NaN 출력을 만나면 해당 지표를 가중
          평균에서 제외하고, 남은 지표의 가중치를 비례 재분배한다. 즉,
          본 함수는 NaN 을 0 으로 대체하지 않으며, 이는 결측 데이터가 점수에
          긍정적/부정적 편향을 일으키지 않도록 하는 의도된 설계이다.
    """
    arr = np.array(values, dtype=float)
    mn, mx = np.nanmin(arr), np.nanmax(arr)
    if mx == mn:
        return np.ones_like(arr) * 0.5
    norm = (arr - mn) / (mx - mn)
    return (1.0 - norm) if invert else norm


# ==========================================
# 통합 랭킹 계산
# ==========================================
def compute_ranking(af2_scores, qmmm_scores, mmgbsa_scores,
                    w_plddt=0.30, w_qmmm=0.40, w_dg=0.30,
                    topology="linear", normalizer="robust",
                    md_status=None,
                    include_tiers=("SUCCESS",)):
    """세 지표(pLDDT, QM/MM 상호작용, MM-GBSA ΔG) 를 정규화·가중합산하여 최종 복합 점수를 계산한다.

    Args:
        af2_scores:    {design_name: {"plddt", "ptm", "iptm"}} dict
        qmmm_scores:   {snapshot_name: {"qmmm_interact_kcal", "converged"}} dict
        mmgbsa_scores: {snapshot_name: {"delta_g_kcal", "favorable"}} dict
        w_plddt, w_qmmm, w_dg: 정규화된 가중치 (합 = 1.0). 합이 1.0 ± 0.01 을
            벗어나면 자동으로 비례 재정규화된다.
        topology: "linear" / "cyclic_*" / ... — cyclic 토폴로지에 +0.05 보정 부여.

    Returns:
        List[dict]: composite_score 내림차순으로 정렬된 entry 리스트. 각 entry 는
        rank, design, plddt, qmmm_interact_kcal, delta_g_kcal, composite_score 등을 포함.

    [v3 7-2] NaN / 결측 처리 정책:
        - QM/MM 또는 MM-GBSA 매칭이 실패한 entry 는 해당 필드를 None 으로 둔다.
        - minmax_normalize 단계에서 None → np.nan 변환되어 min/max 계산에서 제외.
        - 가중 평균 단계에서 결측 지표는 가중치 합 재분배로 자동 보정된다
          (예: pLDDT 만 있는 entry 는 pLDDT 단독 평균으로 계산).
        - 모든 지표가 결측이면 composite_score = 0.0 (사실상 최하위로 정렬됨).
    """
    # [Validation] 가중치 합이 1.0 ± 0.01에서 벗어나면 라이브러리 호출 경로
    # (CLI 없이 compute_ranking을 직접 호출하는 외부 코드 포함)에서도 자동
    # 정규화가 적용되도록 한다. 이는 main()의 CLI 검증과 동일한 안전망을
    # 함수 진입점에 중첩시켜 데이터 신뢰도를 보장한다.
    total_w = w_plddt + w_qmmm + w_dg
    if abs(total_w - 1.0) > 0.01:
        print(f"  [!] compute_ranking 가중치 합({total_w:.3f}) ≠ 1.0 — 자동 정규화")
        w_plddt /= total_w
        w_qmmm  /= total_w
        w_dg    /= total_w

    entries = []
    excluded_entries = []  # type: list

    # 모든 설계 이름 수집
    all_names = set(af2_scores.keys())

    for name in all_names:
        # [v0.5 Bug 4 fix] Tier filter — UNKNOWN은 허용 (legacy run 보호)
        tier = "UNKNOWN"
        if md_status:
            entry_md = md_status.get("designs", {}).get(name, {})
            tier = entry_md.get("tier", "UNKNOWN")

        if tier != "UNKNOWN" and tier not in include_tiers:
            excluded_entries.append({
                "design": name, "tier": tier, "reason": "tier_not_included",
            })
            continue

        entry = {"design": name, "md_tier": tier}

        # AF2 pLDDT
        af2 = af2_scores.get(name, {})
        entry["plddt"]  = af2.get("plddt")
        entry["ptm"]    = af2.get("ptm")
        entry["iptm"]   = af2.get("iptm")

        # [#29 + #36 v0.2] QM/MM: 다중 스냅샷 평균 랭킹
        # - converged=True AND (qm_int_kcal_frozen 또는 interaction_kcal) 유효 스냅샷만 필터
        # - N >= 2: 단순 산술평균 (등간격 스냅샷이므로 가중평균 금지)
        # - N >= 3: σ = np.std(ddof=1)
        # - N < 2: 랭킹 제외 (converged=False, qmmm_interact_kcal=None)
        #
        # 랭킹 score 우선순위:
        #   1. qm_int_kcal_frozen (supermolecular frozen-geometry, 물리적 의미)
        #   2. interaction_kcal (QM-MM embedding coupling, legacy fallback)
        # 두 필드 모두 entry 에 기록되어 CSV 에서 비교 가능.
        qm_matches = fuzzy_match_all(name, list(qmmm_scores.keys()))
        qm_legacy_vals = []   # interaction_kcal (legacy)
        qm_frozen_vals = []   # qm_int_kcal_frozen (Phase 4)
        for m in qm_matches:
            qm = qmmm_scores.get(m, {})
            if not qm.get("converged", False):
                continue
            ie = qm.get("qmmm_interact_kcal")
            if ie is not None:
                qm_legacy_vals.append(float(ie))
            iqf = qm.get("qm_int_kcal_frozen")
            if iqf is not None:
                qm_frozen_vals.append(float(iqf))

        # Phase 4 우선: 모든 valid snap 에서 qm_int_kcal_frozen 가 있으면 그것을
        # 랭킹 score 로 사용. 일부라도 누락이면 legacy interaction_kcal 사용.
        n_legacy = len(qm_legacy_vals)
        n_frozen = len(qm_frozen_vals)
        entry["qmmm_n_snapshots"] = n_legacy

        if n_legacy >= 2:
            entry["qmmm_interact_kcal"] = float(np.mean(qm_legacy_vals))
            entry["qmmm_converged"]     = True
            entry["qmmm_energy_std"]    = (
                float(np.std(qm_legacy_vals, ddof=1)) if n_legacy >= 3 else None
            )
        else:
            entry["qmmm_interact_kcal"] = None
            entry["qmmm_converged"]     = False
            entry["qmmm_energy_std"]    = None

        # Phase 4 mean — 모든 snap 이 frozen 값을 가질 때만 의미 있음.
        if n_frozen >= 2 and n_frozen == n_legacy:
            entry["qm_int_kcal_frozen_mean"] = float(np.mean(qm_frozen_vals))
            entry["qm_int_kcal_frozen_std"]  = (
                float(np.std(qm_frozen_vals, ddof=1)) if n_frozen >= 3 else None
            )
        else:
            entry["qm_int_kcal_frozen_mean"] = None
            entry["qm_int_kcal_frozen_std"]  = None

        # [Phase 6 DIAG] design 별 ranking 집계 결과 — frozen vs legacy 사용 추적
        _diag(
            "rank_design",
            design=entry["design"],
            n_snaps_valid=n_legacy,
            qm_int_kcal_frozen_mean=entry.get("qm_int_kcal_frozen_mean"),
            qmmm_interact_kcal=entry.get("qmmm_interact_kcal"),
            using_frozen=entry.get("qm_int_kcal_frozen_mean") is not None,
        )

        # [Phase 6 DIAG] frozen 누락 → legacy fallback 시 별도 신호
        if (entry.get("qm_int_kcal_frozen_mean") is None
                and entry.get("qmmm_interact_kcal") is not None):
            _diag(
                "rank_fallback",
                design=entry["design"],
                reason="qm_int_kcal_frozen_unavailable",
                fallback_to="interaction_kcal",
            )

        # [#29 v0.2] MM-GBSA: 다중 스냅샷 평균 (QM/MM 과 동일 패턴)
        gb_matches = fuzzy_match_all(name, list(mmgbsa_scores.keys()))
        gb_valid_vals = []
        gb_valid_corr = []
        gb_favorable_any = False
        gb_minus_TdS = None
        for m in gb_matches:
            gb = mmgbsa_scores.get(m, {})
            dg = gb.get("delta_g_kcal")
            if dg is not None:
                gb_valid_vals.append(float(dg))
                if gb.get("favorable", False):
                    gb_favorable_any = True
                dgc = gb.get("delta_g_corrected_kcal")
                if dgc is not None:
                    gb_valid_corr.append(float(dgc))
                if gb_minus_TdS is None:
                    gb_minus_TdS = gb.get("minus_TdS_kcal")
        n_gb_snaps = len(gb_valid_vals)
        entry["mmgbsa_n_snapshots"] = n_gb_snaps
        if n_gb_snaps >= 2:
            entry["delta_g_kcal"] = float(np.mean(gb_valid_vals))
            entry["favorable"]    = gb_favorable_any
            entry["delta_g_corrected_kcal"] = (
                float(np.mean(gb_valid_corr)) if gb_valid_corr else None
            )
            entry["delta_g_std"] = (
                float(np.std(gb_valid_vals, ddof=1)) if n_gb_snaps >= 3 else None
            )
            entry["minus_TdS_kcal"] = gb_minus_TdS
        else:
            entry["delta_g_kcal"] = None
            entry["favorable"]    = False
            entry["delta_g_corrected_kcal"] = None
            entry["delta_g_std"] = None
            entry["minus_TdS_kcal"] = None

        entries.append(entry)

    if not entries:
        return []

    # [v3 7-1] fuzzy_match 매칭 통계 및 미연결 설계 진단
    n_qmmm_matched = sum(1 for e in entries if e["qmmm_interact_kcal"] is not None)
    n_dg_matched   = sum(1 for e in entries if e["delta_g_kcal"]      is not None)
    print(f"  [매칭] AF2 {len(entries)}개 중 QM/MM {n_qmmm_matched}개, MM-GBSA {n_dg_matched}개 연결 성공")

    if n_qmmm_matched < len(entries):
        unmatched_qmmm = [e["design"] for e in entries if e["qmmm_interact_kcal"] is None]
        head = unmatched_qmmm[:5]
        suffix = "..." if len(unmatched_qmmm) > 5 else ""
        print(f"  [WARNING] QM/MM 미연결 설계 ({len(unmatched_qmmm)}개): {head}{suffix}")
    if n_dg_matched < len(entries):
        unmatched_dg = [e["design"] for e in entries if e["delta_g_kcal"] is None]
        head = unmatched_dg[:5]
        suffix = "..." if len(unmatched_dg) > 5 else ""
        print(f"  [WARNING] MM-GBSA 미연결 설계 ({len(unmatched_dg)}개): {head}{suffix}")

    # ── 정규화 ────────────────────────────────────────────────
    # QM/MM 정규화 score 우선순위:
    #   1. qm_int_kcal_frozen_mean (supermolecular frozen-geometry — 진정한
    #      QM 상호작용 에너지 proxy, calibration 적합)
    #   2. qmmm_interact_kcal (legacy QM-MM embedding coupling, fallback)
    # entry 별로 frozen 값이 있으면 frozen 사용, 없으면 legacy 사용.
    plddt_vals  = [e["plddt"]             if e["plddt"]             is not None else np.nan for e in entries]
    qmmm_vals   = [
        (e.get("qm_int_kcal_frozen_mean") if e.get("qm_int_kcal_frozen_mean") is not None
         else (e["qmmm_interact_kcal"] if e["qmmm_interact_kcal"] is not None else np.nan))
        for e in entries
    ]
    # [v4 S-1] ΔG 정규화는 IE 보정 ΔG (delta_g_corrected_kcal) 가 있으면 우선 사용
    dg_vals     = [
        (e.get("delta_g_corrected_kcal") if e.get("delta_g_corrected_kcal") is not None
         else (e["delta_g_kcal"] if e["delta_g_kcal"] is not None else np.nan))
        for e in entries
    ]

    # [v4 S-4] 정규화 방식 선택: "robust" (백분위 기반, 이상치 방어) 또는 "minmax" (legacy)
    if normalizer == "robust":
        plddt_norm = robust_percentile_normalize(plddt_vals, invert=False)
        qmmm_norm  = robust_percentile_normalize(qmmm_vals,  invert=True)
        dg_norm    = robust_percentile_normalize(dg_vals,    invert=True)
    else:
        plddt_norm = minmax_normalize(plddt_vals, invert=False)
        qmmm_norm  = minmax_normalize(qmmm_vals,  invert=True)
        dg_norm    = minmax_normalize(dg_vals,    invert=True)

    # Cyclic 페널티 보정 (선택적)
    # NOTE [v4 S-1]: cyclic_bonus 0.05 임의값은 IE 가 활성화된 경우 (delta_g_corrected_kcal
    # 가 entries 에 채워진 경우) 자연스럽게 −TΔS 항으로 흡수되므로 IE-aware mode 에서는
    # 비활성화한다. 이는 IE 가 cyclic 펩타이드의 conformational 자유도 차이를 이미 반영하기 때문이다.
    has_ie_corrected = any(e.get("delta_g_corrected_kcal") is not None for e in entries)
    if has_ie_corrected:
        cyclic_bonus = 0.0
    else:
        cyclic_bonus = 0.05 if topology not in ("linear", None) else 0.0

    for i, entry in enumerate(entries):
        s_plddt = plddt_norm[i] if not np.isnan(plddt_norm[i]) else 0.0
        s_qmmm  = qmmm_norm[i]  if not np.isnan(qmmm_norm[i])  else 0.0
        s_dg    = dg_norm[i]    if not np.isnan(dg_norm[i])    else 0.0

        # 데이터 누락 시 해당 가중치 재분배
        available = []
        weights   = []
        if entry["plddt"] is not None:
            available.append(s_plddt)
            weights.append(w_plddt)
        # qm score 가용 조건: frozen 또는 legacy 어느 쪽이든 1 개 이상.
        if (entry.get("qm_int_kcal_frozen_mean") is not None
                or entry["qmmm_interact_kcal"] is not None):
            available.append(s_qmmm)
            weights.append(w_qmmm)
        if entry["delta_g_kcal"] is not None:
            available.append(s_dg)
            weights.append(w_dg)

        if not available:
            entry["composite_score"] = 0.0
        else:
            total_w = sum(weights)
            entry["composite_score"] = sum(
                a * w / total_w for a, w in zip(available, weights)
            ) + cyclic_bonus

        entry["score_plddt_norm"]  = round(float(s_plddt), 4)
        entry["score_qmmm_norm"]   = round(float(s_qmmm),  4)
        entry["score_dg_norm"]     = round(float(s_dg),    4)
        entry["composite_score"]   = round(entry["composite_score"], 4)

    # 복합 점수 내림차순 정렬
    entries.sort(key=lambda e: e["composite_score"], reverse=True)

    # 순위 부여
    for rank, entry in enumerate(entries, 1):
        entry["rank"] = rank

    # [v0.5] 제외된 entries를 함수 속성으로 전달 (시그니처 변경 최소화)
    compute_ranking._excluded_entries = excluded_entries  # type: ignore[attr-defined]

    # [Phase 6 DIAG] 최종 ranking 완료 — top1/매칭 수 추적
    _diag(
        "rank_complete",
        n_designs=len(entries),
        n_qmmm_matched=n_qmmm_matched,
        n_dg_matched=n_dg_matched,
        top1_design=(entries[0]["design"] if entries else None),
        top1_score=(entries[0].get("composite_score") if entries else None),
    )

    return entries


# ==========================================
# CSV 저장
# ==========================================
def save_bucket_reports(output_csv, excluded_entries):
    # type: (str, list) -> None
    """[v0.5] Tier별 bucket report CSV 생성.

    ranking_marginal.csv, ranking_partial.csv, ranking_excluded.csv
    """
    if not excluded_entries:
        return

    base_dir = os.path.dirname(output_csv) or "."
    base_name = os.path.basename(output_csv).replace(".csv", "")

    buckets = {}  # type: Dict[str, list]
    for ex in excluded_entries:
        tier = ex.get("tier", "UNKNOWN").lower()
        buckets.setdefault(tier, []).append(ex)

    fieldnames = ["design", "tier", "reason"]

    for tier_key, items in buckets.items():
        bucket_csv = os.path.join(base_dir, f"{base_name}_{tier_key}.csv")
        with open(bucket_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(items)
        print(f"  [v0.5] Bucket report: {bucket_csv} ({len(items)} entries)")

    all_excluded_csv = os.path.join(base_dir, f"{base_name}_excluded.csv")
    with open(all_excluded_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(excluded_entries)
    print(f"  [v0.5] All excluded: {all_excluded_csv} ({len(excluded_entries)} entries)")


def save_ranking_csv(entries, output_csv):
    if not entries:
        print("  [!] 저장할 데이터 없음")
        return

    fieldnames = [
        "rank", "design", "md_tier",
        "composite_score",
        "plddt", "ptm", "iptm",
        # [#29 v0.2] 다중 스냅샷 개수 필드 추가
        "qmmm_interact_kcal", "qmmm_converged", "qmmm_energy_std", "qmmm_n_snapshots",
        # supermolecular QM interaction energy (frozen geometry)
        "qm_int_kcal_frozen_mean", "qm_int_kcal_frozen_std",
        "delta_g_kcal", "delta_g_corrected_kcal", "delta_g_std", "minus_TdS_kcal", "favorable",
        "mmgbsa_n_snapshots",
        "score_plddt_norm", "score_qmmm_norm", "score_dg_norm",
    ]

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(entries)

    print(f"  CSV 저장: {output_csv}")


# ==========================================
# 최종 리포트 출력
# ==========================================
def print_report(entries, topology, w_plddt, w_qmmm, w_dg):
    if not entries:
        print(" [!] 출력할 랭킹 데이터가 없습니다.")
        return
    print(f"\n{'='*72}")
    print(f"  UPDD 최종 랭킹 리포트")
    print(f"  토폴로지: {topology}")
    print(f"  가중치  : pLDDT {w_plddt*100:.0f}%  |  QM/MM {w_qmmm*100:.0f}%  |  MM-GBSA {w_dg*100:.0f}%")
    print(f"{'='*72}")
    print(f"  {'순위':<4} {'설계명':<36} {'복합점수':>8}  {'pLDDT':>7}  {'QM-MM(kcal)':>12}  {'ΔG(kcal)':>10}")
    print(f"  {'-'*70}")

    for e in entries[:10]:  # 상위 10개 출력
        plddt  = f"{e['plddt']:.1f}"    if e['plddt']             is not None else "N/A"
        qmmm   = f"{e['qmmm_interact_kcal']:+.2f}" if e['qmmm_interact_kcal'] is not None else "N/A"
        dg     = f"{e['delta_g_kcal']:+.2f}"       if e['delta_g_kcal']       is not None else "N/A"
        flag   = "✅" if e.get("favorable") else "  "
        print(f"  {e['rank']:<4} {e['design'][:35]:<36} {e['composite_score']:>8.4f}  "
              f"{plddt:>7}  {qmmm:>12}  {dg:>10} {flag}")

    print(f"{'='*72}")
    print(f"  🏆 최우선 후보: {entries[0]['design']}")
    print(f"     복합 점수  : {entries[0]['composite_score']:.4f}")
    if entries[0]['plddt']:
        print(f"     pLDDT      : {entries[0]['plddt']:.1f}")
    if entries[0]['delta_g_kcal']:
        print(f"     ΔG_bind    : {entries[0]['delta_g_kcal']:+.4f} kcal/mol")
    print(f"{'='*72}\n")


# ==========================================
# CLI 메인
# ==========================================
def main():
    parser = argparse.ArgumentParser(
        description="Final Ranking — pLDDT + QM/MM + MM-GBSA 통합 랭킹"
    )
    parser.add_argument("--af2_dir",    required=True,  help="AF2 결과 디렉토리")
    parser.add_argument("--qmmm_dir",   required=True,  help="QM/MM 결과 디렉토리")
    parser.add_argument("--mmgbsa_dir", required=True,  help="MM-GBSA 결과 디렉토리")
    parser.add_argument("--outputcsv",  required=True,  help="출력 CSV 경로")
    parser.add_argument("--topology",   default="linear", help="펩타이드 토폴로지")
    parser.add_argument("--w_plddt",    type=float, default=0.30, help="pLDDT 가중치")
    parser.add_argument("--w_qmmm",     type=float, default=0.40, help="QM/MM 가중치")
    parser.add_argument("--w_dg",       type=float, default=0.30, help="MM-GBSA 가중치")
    parser.add_argument("--normalizer", choices=["robust", "minmax"], default="robust",
                        help="[v4 S-4] 점수 정규화 방식 (robust=백분위 5/95 clamp, minmax=legacy)")
    parser.add_argument("--md-status", type=str, default=None,
                        help="[v0.5] _md_status.json 경로 (없으면 자동 탐색)")
    parser.add_argument("--include-tiers", nargs="+",
                        default=["SUCCESS"],
                        choices=["SUCCESS", "MARGINAL", "PARTIAL", "FAIL"],
                        help="[v0.5] Ranking에 포함할 MD tier (기본: SUCCESS only)")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.outputcsv) or ".", exist_ok=True)

    # 가중치 합 검증
    total_w = args.w_plddt + args.w_qmmm + args.w_dg
    if abs(total_w - 1.0) > 0.01:
        print(f"  [!] 가중치 합({total_w:.2f}) ≠ 1.0 — 자동 정규화")
        args.w_plddt /= total_w
        args.w_qmmm  /= total_w
        args.w_dg    /= total_w

    print(f"\n[Final Ranking]")
    print(f"  토폴로지: {args.topology}")
    print(f"  가중치  : pLDDT {args.w_plddt:.0%} | QM/MM {args.w_qmmm:.0%} | MM-GBSA {args.w_dg:.0%}")

    # 데이터 로드
    af2_scores    = load_af2_scores(args.af2_dir)
    qmmm_scores   = load_qmmm_scores(args.qmmm_dir)
    mmgbsa_scores = load_mmgbsa_scores(args.mmgbsa_dir)

    if not af2_scores:
        print("  [!] AF2 점수 없음 — QM/MM + MM-GBSA만으로 랭킹")

    # [v0.5] MD status 로드
    md_status = None
    if args.md_status:
        try:
            with open(args.md_status, "r", encoding="utf-8") as f:
                md_status = json.load(f)
            print(f"  [v0.5] MD status: {args.md_status}")
        except (IOError, json.JSONDecodeError) as e:
            print(f"  [!] MD status 로드 실패: {e} — tier 필터 비활성화")
    else:
        # 자동 탐색: af2_dir 상위 디렉토리에서 _md_status.json 탐색
        _parent = os.path.dirname(args.af2_dir)
        _auto = os.path.join(_parent, "_md_status.json")
        if os.path.exists(_auto):
            md_status = load_md_status(_parent)
            if md_status:
                print(f"  [v0.5] MD status auto-detected: {_auto}")

    include_tiers = tuple(args.include_tiers)
    if md_status:
        print(f"  [v0.5] Include tiers: {include_tiers}")

    # 랭킹 계산
    entries = compute_ranking(
        af2_scores, qmmm_scores, mmgbsa_scores,
        w_plddt    = args.w_plddt,
        w_qmmm     = args.w_qmmm,
        w_dg       = args.w_dg,
        topology   = args.topology,
        normalizer = args.normalizer,
        md_status  = md_status,
        include_tiers = include_tiers,
    )

    if not entries:
        print("  [!] 랭킹할 데이터 없음")
        return

    # CSV 저장 + 리포트 출력
    save_ranking_csv(entries, args.outputcsv)
    print_report(entries, args.topology, args.w_plddt, args.w_qmmm, args.w_dg)

    # [v0.5] Bucket reports for excluded entries
    excluded = getattr(compute_ranking, "_excluded_entries", [])
    if excluded:
        save_bucket_reports(args.outputcsv, excluded)

    # 요약 JSON (상위 5개)
    summary_json = args.outputcsv.replace(".csv", "_summary.json")
    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump({
            "topology":  args.topology,
            "weights":   {"plddt": args.w_plddt, "qmmm": args.w_qmmm, "dg": args.w_dg},
            "n_total":   len(entries),
            "top5":      entries[:5],
        }, f, indent=2, ensure_ascii=False)

    print(f"  요약 JSON: {summary_json}")
    print(f"\n  파이프라인 완료! 🎉  최종 결과: {args.outputcsv}")


if __name__ == "__main__":
    main()
