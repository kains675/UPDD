#!/usr/bin/env python
"""
UPDD_benchmark.py — UPDD 파이프라인 정량적 벤치마킹 프레임워크

목적:
    UPDD 파이프라인이 산출하는 결합 친화도 순위가 실험적으로 측정된 친화도
    순위와 얼마나 일치하는지를 정량적으로 평가한다. 핵심 지표는 절대값이
    아닌 **순위 상관 (rank correlation)** 이다.

설계 원칙 (CLAUDE.md):
    - 기존 *.py 파일을 import 하지 않는다. subprocess 로만 유틸리티를 호출한다.
    - 벤치마크에서는 RFdiffusion/ProteinMPNN/AlphaFold2 단계를 건너뛴다.
      실험 복합체 구조 (PDB) 에서 직접 MD → QM/MM → MM-GBSA 를 실행하여
      스코어링 단계의 순위 예측 능력만 순수하게 평가한다.
    - 과학 라이브러리 (scipy, numpy, sklearn, matplotlib) 는 evaluate 시점에
      지연 import 한다.

사용법:
    python UPDD_benchmark.py                    # standard 모드 (기본)
    python UPDD_benchmark.py --mode quick       # 빠른 검증 (10분)
    python UPDD_benchmark.py --mode full        # 정밀 벤치마크 (수 일)
    python UPDD_benchmark.py --skip-run         # 기존 결과로 평가만 수행
    python UPDD_benchmark.py --topology linear  # linear 만 평가

디렉토리 구조:
    UPDD_benchmark/
    ├── UPDD_benchmark_dataset.csv    ← 입력 (사전 배치 필요)
    ├── pdbs/                         ← 자동 다운로드
    ├── results/                      ← 파이프라인 결과
    ├── report/                       ← 평가 리포트 + 시각화
    ├── benchmark_progress.json       ← 재시작 지원
    └── benchmark.log                 ← 실행 로그
"""
from __future__ import annotations

import os
import sys
import csv
import json
import glob
import math
import time
import shutil
import logging
import argparse
import subprocess
import urllib.request
import urllib.error
from collections import defaultdict


# ==========================================
# 1. 환경 및 경로 설정 (UPDD.py 와 동일 규약)
# ==========================================
HOME_DIR  = os.environ.get("UPDD_HOME", os.path.expanduser("~"))
UPDD_DIR  = os.environ.get("UPDD_DIR", os.path.join(HOME_DIR, "UPDD_proj"))
UTILS_DIR = os.path.join(UPDD_DIR, "utils")

ENV_NAMES = {
    "md":   "md_simulation",
    "qmmm": "qmmm",
}

SCRIPT_PATHS = {
    "preprocess_target": os.path.join(UTILS_DIR, "preprocess_target.py"),
    "run_md":            os.path.join(UTILS_DIR, "run_restrained_md.py"),
    "extract_snap":      os.path.join(UTILS_DIR, "extract_snapshots.py"),
    "run_qmmm":          os.path.join(UTILS_DIR, "run_qmmm.py"),
    "run_mmgbsa":        os.path.join(UTILS_DIR, "run_mmgbsa.py"),
    "run_cyclic_fold":   os.path.join(UTILS_DIR, "run_cyclic_fold.py"),
}

RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

# 모든 벤치마크 산출물은 UPDD_benchmark/ 하위에 자동 배치된다.
BENCHMARK_DIR   = os.path.join(UPDD_DIR, "UPDD_benchmark")
DATASET_CSV     = os.path.join(BENCHMARK_DIR, "UPDD_benchmark_dataset.csv")
PDB_DIR         = os.path.join(BENCHMARK_DIR, "pdbs")
PREPROCESS_DIR  = os.path.join(BENCHMARK_DIR, "preprocessed")
AF2_INPUT_DIR   = os.path.join(BENCHMARK_DIR, "af2_inputs")
AF2_RESULT_DIR  = os.path.join(BENCHMARK_DIR, "af2_results")
RESULTS_DIR     = os.path.join(BENCHMARK_DIR, "results")
REPORT_DIR      = os.path.join(BENCHMARK_DIR, "report")
FIGURES_DIR     = os.path.join(REPORT_DIR, "figures")
PROGRESS_FILE   = os.path.join(BENCHMARK_DIR, "benchmark_progress.json")
LOG_FILE        = os.path.join(BENCHMARK_DIR, "benchmark.log")

# UPDD.py 와 동일 규약으로 localcolabfold 경로 탐색을 수행한다.
PROJECT_DIR = os.environ.get("UPDD_PROJECT_DIR", os.path.join(HOME_DIR, "ai_projects"))

# 표준 아미노산 1-letter → 3-letter
AA1_TO_AA3 = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
}
AA3_TO_AA1 = {v: k for k, v in AA1_TO_AA3.items()}
STD_AA1 = set(AA1_TO_AA3.keys())

# 히스티딘 변이형을 정규 HIS 로 흡수
HIS_VARIANTS = {"HID", "HIE", "HIP", "HSD", "HSE", "HSP"}
for _v in HIS_VARIANTS:
    AA3_TO_AA1[_v] = "H"

# 모드별 MD/QMMM 파라미터 (CLAUDE.md 실행 모드 절 참조)
MODE_PRESETS = {
    "quick": {
        "md_steps":   10000,
        "run_qmmm":   False,
        "run_mmgbsa": True,
        "qmmm_mode":  "fast",
    },
    "standard": {
        "md_steps":   50000,
        "run_qmmm":   True,
        "run_mmgbsa": True,
        "qmmm_mode":  "fast",
    },
    "full": {
        "md_steps":   5000000,
        "run_qmmm":   True,
        "run_mmgbsa": True,
        "qmmm_mode":  "full",
    },
}


# ==========================================
# 2. 로거
# ==========================================
def setup_benchmark_logger(base_dir: str) -> logging.Logger:
    """benchmark.log 단일 파일을 base_dir 하위에 작성한다 (append 모드)."""
    os.makedirs(base_dir, exist_ok=True)
    log_file = os.path.join(base_dir, "benchmark.log")
    logger = logging.getLogger("UPDD_benchmark")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    for h in logger.handlers[:]:
        logger.removeHandler(h)
    fh = logging.FileHandler(log_file, mode="a", encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(logging.INFO)
    sh.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(fh)
    logger.addHandler(sh)
    logger.info(f"[LOG] benchmark 로그: {log_file}")
    return logger


# ------------------------------------------
# 2b. PBC unwrap 검증 헬퍼 (Fix 4: 방어적)
# ------------------------------------------
def _verify_snapshot_proximity(snap_dir: str,
                               target_chain: str = "A",
                               binder_chain: str = "B",
                               max_dist: float = 15.0,
                               logger: logging.Logger | None = None) -> bool:
    """스냅샷 PDB 에서 타겟-바인더 간 최소 원자 거리를 검증한다.

    extract_snapshots.py 의 PBC unwrap (v39 Fix 1) 이 정상 적용되었는지
    확인하는 방어적 단계이다. max_dist (기본 15 Å) 초과 시 경고를 출력하고
    False 를 반환한다 (하류 QM/MM 이 "도망자 발견" 으로 빠질 가능성이 높다).
    """
    snap_pdbs = glob.glob(os.path.join(snap_dir, "*.pdb"))
    if not snap_pdbs:
        return True  # 스냅샷이 없으면 검증 skip
    for sp in snap_pdbs:
        a_xyz, b_xyz = [], []
        with open(sp, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                chain = line[21:22]
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                except ValueError:
                    continue
                if chain == target_chain:
                    a_xyz.append((x, y, z))
                elif chain == binder_chain:
                    b_xyz.append((x, y, z))
        if not a_xyz or not b_xyz:
            continue
        import math as _m
        min_d = min(
            _m.sqrt((ax - bx)**2 + (ay - by)**2 + (az - bz)**2)
            for (ax, ay, az) in a_xyz[::10]  # 성능: 1/10 서브샘플
            for (bx, by, bz) in b_xyz
        )
        if min_d > max_dist:
            if logger:
                logger.warning(
                    f"  [!] PBC unwrap 미완: {os.path.basename(sp)}, "
                    f"min(A-B)={min_d:.1f} Å > {max_dist} Å (QM/MM '도망자' 예상)"
                )
            return False
    return True


# ------------------------------------------
# 2c. 진행 상황 표시 헬퍼 (실시간 단계 로그)
# ------------------------------------------
_STAGE_ICONS = {"start": "⏳", "ok": "✅", "fail": "❌", "skip": "⏭️"}
_TOTAL_STAGES = 6  # preprocess, AF2, MD, snapshot, QM/MM, MM-GBSA


def _log_stage(logger: logging.Logger,
               stage_num: int,
               total_stages: int,
               stage_name: str,
               status: str,
               elapsed: float | None = None,
               detail: str = "") -> None:
    """단계별 진행 상황을 일관된 형식으로 출력한다.

    Args:
        stage_num/total_stages: 예) 1/6, 2/6 ...
        stage_name: 한국어 단계명 (예: "타겟 전처리")
        status: start / ok / fail / skip
        elapsed: 초 단위 (ok/fail/skip 시에만 의미)
        detail: 보조 메시지 (예: "chain A, 85 aa")
    """
    icon = _STAGE_ICONS.get(status, "?")
    time_str = f"({elapsed:.1f}s)" if elapsed is not None else ""
    detail_str = f" {detail}" if detail else ""
    # stage_name 을 고정 폭 (20) 으로 맞추기 위해 패딩. 한글은 2칸 차지하므로
    # 한글 문자에 대해 한 칸을 추가해서 시각적 정렬을 맞춘다.
    han_chars = sum(1 for c in stage_name if ord(c) > 0x7F)
    pad = max(0, 20 - (len(stage_name) + han_chars))
    padded = stage_name + " " * pad
    logger.info(f"  [{stage_num}/{total_stages}] {padded} {icon} {time_str}{detail_str}")


def _fmt_eta(seconds: float) -> str:
    if seconds >= 3600:
        return f"{seconds/3600:.1f}h"
    if seconds >= 60:
        return f"{seconds/60:.0f}m"
    return f"{seconds:.0f}s"


def _print_progress_summary(logger: logging.Logger,
                            idx: int, total: int,
                            n_success: int, n_fail: int,
                            elapsed_total: float) -> None:
    """전체 진행률, 성공/실패 집계, ETA 를 한 줄로 출력한다."""
    pct = (idx / total * 100.0) if total > 0 else 0.0
    avg = (elapsed_total / idx) if idx > 0 else 0.0
    remaining = (total - idx) * avg
    eta_str = _fmt_eta(remaining)
    logger.info(
        f"\n{'─' * 4} 전체 진행률: {idx}/{total} ({pct:.1f}%) "
        f"| 성공: {n_success} | 실패: {n_fail} "
        f"| ETA: ~{eta_str} {'─' * 4}\n"
    )


def _print_final_summary(logger: logging.Logger,
                         progress: dict,
                         total_elapsed: float) -> None:
    """모든 레코드 완료 후 상태 집계 + 실패 단계 분포 + ΔG 통계를 출력한다."""
    statuses: dict[str, int] = defaultdict(int)
    failed_stages: dict[str, int] = defaultdict(int)
    for _, info in progress.items():
        statuses[info.get("status", "UNKNOWN")] += 1
        if info.get("failed_stage"):
            failed_stages[info["failed_stage"]] += 1

    hours = total_elapsed / 3600.0
    bar = "=" * 55
    logger.info(f"\n{bar}")
    logger.info("  UPDD Benchmark 최종 요약")
    logger.info(bar)
    logger.info(f"  총 소요 시간: {hours:.1f}h ({total_elapsed:.0f}s)")
    logger.info(f"  총 레코드:    {len(progress)} 개")
    for st, cnt in sorted(statuses.items()):
        logger.info(f"    {st:<12}: {cnt} 개")

    if failed_stages:
        logger.info("\n  실패 단계 분포:")
        for stage, cnt in sorted(failed_stages.items(), key=lambda x: -x[1]):
            logger.info(f"    {stage:<16}: {cnt} 개")

    dg_values = [v["delta_g_kcal"] for v in progress.values()
                 if v.get("delta_g_kcal") is not None]
    if dg_values:
        import statistics
        logger.info(f"\n  ΔG 통계 (성공 {len(dg_values)} 개):")
        logger.info(f"    평균: {statistics.mean(dg_values):+.2f} kcal/mol")
        logger.info(f"    범위: {min(dg_values):+.2f} ~ {max(dg_values):+.2f} kcal/mol")
    logger.info(f"{bar}\n")


# ==========================================
# 3. 데이터셋 로더
# ==========================================
def clean_peptide_sequence(raw: str) -> str:
    """CSV 서열 필드를 표준 1-letter 서열로 정규화한다.

    처리:
    - cyclo(...) 래퍼 제거
    - 소문자 D-아미노산 마커 (f, v 등) 는 동일 알파벳 대문자로 취급
    - pSer/pTyr 등 인산화 접두사 (p) 제거
    - 표준 20 AA 외 문자 (괄호, 하이픈, ncAA 이름) 는 제외
    - 결과가 비면 원본 대문자 필터로 폴백
    """
    if not raw:
        return ""
    s = raw.strip()
    # cyclo() 래퍼
    if s.lower().startswith("cyclo(") and s.endswith(")"):
        s = s[6:-1]
    # 인산화 접두사 (pSer, pTyr, pY)
    s = s.replace("pS", "S").replace("pT", "T").replace("pY", "Y")
    # 하이픈으로 연결된 ncAA 이름 (예: YFWLKPC-Dap-Gly-Ser) 는 하이픈 앞까지만 취한다
    if "-" in s:
        s = s.split("-")[0]
    # 표준 AA 만 추출 (대소문자 모두)
    result = "".join(c for c in s.upper() if c in STD_AA1)
    return result


def load_dataset(csv_path: str) -> list[dict]:
    """벤치마크 CSV 를 읽어 레코드 리스트를 반환한다.

    - `#` 으로 시작하는 줄은 섹션 주석이므로 건너뛴다.
    - peptide_id 가 비어있는 줄은 건너뛴다.
    - kd_nm 이 비어있으면 ic50_nm 을 Cheng-Prusoff 근사 (Kd ≈ IC50/2) 로 대체한다.
    - 실험값이 전혀 없으면 레코드에 None 으로 남겨둔다 (음성 대조군 처리는 run 시점에서 결정).
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"데이터셋 CSV 없음: {csv_path}")

    records = []
    with open(csv_path, "r", encoding="utf-8") as f:
        # 섹션 주석 제거하면서 DictReader 에 전달
        cleaned = (ln for ln in f if ln.strip() and not ln.lstrip().startswith("#"))
        reader = csv.DictReader(cleaned)
        for row in reader:
            pid = (row.get("peptide_id") or "").strip()
            if not pid:
                continue
            kd_nm = _safe_float(row.get("kd_nm"))
            ic50_nm = _safe_float(row.get("ic50_nm"))
            if kd_nm is None and ic50_nm is not None:
                kd_nm = ic50_nm / 2.0
            rec = {
                "peptide_id":     pid,
                "sequence_raw":   (row.get("sequence") or "").strip(),
                "sequence":       clean_peptide_sequence(row.get("sequence") or ""),
                "target_name":    (row.get("target_name") or "").strip(),
                "target_uniprot": (row.get("target_uniprot") or "").strip(),
                "pdb_id":         (row.get("pdb_id") or "").strip().upper(),
                "topology":       (row.get("topology") or "linear").strip().lower(),
                "kd_nm":          kd_nm,
                "ic50_nm":        ic50_nm,
                "assay_type":     (row.get("assay_type") or "").strip(),
                "reference":      (row.get("reference") or "").strip(),
                "quality_score":  _safe_int(row.get("quality_score"), default=0),
                "notes":          (row.get("notes") or "").strip(),
                "keep_hetatms":   (row.get("keep_hetatms") or "").strip(),
            }
            records.append(rec)
    return records


def _safe_float(val):
    if val is None:
        return None
    s = str(val).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _safe_int(val, default=0):
    if val is None:
        return default
    s = str(val).strip()
    if not s:
        return default
    try:
        return int(s)
    except ValueError:
        return default


def filter_dataset(records: list[dict],
                   topology: str | None = None,
                   min_quality: int = 0,
                   max_samples: int | None = None,
                   require_pdb: bool = False) -> list[dict]:
    """레코드 리스트를 CLI 인자에 맞게 필터링한다."""
    out = []
    for r in records:
        if topology and r["topology"] != topology.lower():
            continue
        if r["quality_score"] < min_quality:
            continue
        if require_pdb and not r["pdb_id"]:
            continue
        out.append(r)
        if max_samples is not None and len(out) >= max_samples:
            break
    return out


def assign_labels(records: list[dict]) -> None:
    """각 레코드에 binary label (1=binder, 0=non-binder) 을 할당한다.

    - peptide_id 가 neg_* / decoy_* 로 시작하면 0
    - Kd > 10 μM (10000 nM) 이면 0
    - 그 외 유효 Kd 있는 경우 1
    - Kd 정보 없으면 None (평가에서 제외)
    """
    for r in records:
        pid = r["peptide_id"].lower()
        if pid.startswith("neg_") or pid.startswith("decoy_"):
            r["label"] = 0
        elif r["kd_nm"] is None:
            r["label"] = None
        elif r["kd_nm"] > 10000.0:
            r["label"] = 0
        else:
            r["label"] = 1


# ==========================================
# 4. PDB 다운로드
# ==========================================
def download_all_pdbs(records: list[dict],
                      pdb_dir: str,
                      logger: logging.Logger) -> tuple[int, int, int]:
    """데이터셋의 모든 고유 PDB 를 RCSB 에서 다운로드한다.

    이미 존재하면 건너뛴다. 빈 pdb_id 레코드는 조용히 skip 한다.

    Returns:
        (n_downloaded, n_skipped, n_failed)
    """
    os.makedirs(pdb_dir, exist_ok=True)
    unique_ids = sorted({r["pdb_id"] for r in records if r.get("pdb_id")})
    logger.info(f"[download] 고유 PDB ID: {len(unique_ids)}")

    n_downloaded, n_skipped, n_failed = 0, 0, 0
    failed_ids: list[tuple[str, str]] = []

    for pdb_id in unique_ids:
        out_path = os.path.join(pdb_dir, f"{pdb_id}.pdb")
        if os.path.exists(out_path) and os.path.getsize(out_path) > 0:
            n_skipped += 1
            continue
        url = RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id)
        try:
            urllib.request.urlretrieve(url, out_path)
            if os.path.getsize(out_path) == 0:
                raise urllib.error.URLError("빈 응답")
            logger.info(f"  [DL] {pdb_id}.pdb 다운로드 완료")
            n_downloaded += 1
        except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
            logger.warning(f"  [DL] {pdb_id}.pdb 다운로드 실패: {e}")
            n_failed += 1
            failed_ids.append((pdb_id, str(e)))
            if os.path.exists(out_path) and os.path.getsize(out_path) == 0:
                os.remove(out_path)

    # 비어있는 pdb_id 를 가진 데이터 포인트는 보고만 한다
    missing_pdb = [r["peptide_id"] for r in records if not r.get("pdb_id")]
    if missing_pdb:
        logger.info(f"[download] PDB 미지정 레코드 ({len(missing_pdb)}): "
                    f"{', '.join(missing_pdb[:10])}{' …' if len(missing_pdb) > 10 else ''}")

    if failed_ids:
        fail_csv = os.path.join(pdb_dir, "_download_failures.csv")
        with open(fail_csv, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["pdb_id", "error"])
            w.writerows(failed_ids)
        logger.info(f"[download] 실패 목록: {fail_csv}")

    return n_downloaded, n_skipped, n_failed


# ==========================================
# 5. 체인 식별 + PDB 전처리 헬퍼
# ==========================================
def parse_pdb_chains(pdb_path: str) -> dict[str, list[tuple[int, str]]]:
    """PDB 파일에서 체인별 (resnum, resname3) 순서 리스트를 추출한다.

    CA 원자만 수집하여 잔기 순서를 보존한다. ALTLOC 은 첫 번째만 사용한다.
    """
    chains: dict[str, list[tuple[int, str]]] = defaultdict(list)
    seen: set[tuple[str, int, str]] = set()
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            altloc = line[16:17]
            if altloc not in (" ", "A"):
                continue
            resname = line[17:20].strip()
            chain_id = line[21:22].strip() or " "
            try:
                resnum = int(line[22:26])
            except ValueError:
                continue
            key = (chain_id, resnum, line[26:27])
            if key in seen:
                continue
            seen.add(key)
            chains[chain_id].append((resnum, resname))
    return chains


def chain_to_sequence(chain_residues: list[tuple[int, str]]) -> str:
    """체인의 (resnum, resname3) 리스트를 1-letter 서열로 변환한다."""
    out = []
    for _, resname in chain_residues:
        out.append(AA3_TO_AA1.get(resname.upper(), "X"))
    return "".join(out)


def sequence_similarity(a: str, b: str) -> float:
    """difflib.SequenceMatcher 기반 부분 일치율.

    벤치마크 목적상 두 서열 중 짧은 쪽 (통상 펩타이드) 이 긴 쪽의 부분
    서열로 포함되어 있는지 확인한다. cyclic 펩타이드가 선형 절단으로
    기록된 경우의 회전 대응을 위해 rotations 도 체크한다.
    """
    if not a or not b:
        return 0.0
    import difflib
    shorter, longer = (a, b) if len(a) <= len(b) else (b, a)
    best = difflib.SequenceMatcher(None, shorter, longer).ratio()
    # cyclic rotation: shorter 서열을 한 바퀴 돌려가며 최고값을 취한다
    if len(shorter) >= 4:
        for i in range(1, min(len(shorter), 10)):
            rot = shorter[i:] + shorter[:i]
            r = difflib.SequenceMatcher(None, rot, longer).ratio()
            if r > best:
                best = r
    return best


def identify_target_chain(pdb_path: str) -> tuple[str | None, int]:
    """raw PDB 에서 가장 긴 단백질 체인을 타겟으로 식별한다.

    벤치마크 AF2 파이프라인용 간소화:
        바인더 체인 식별은 AF2 가 담당하므로, 여기서는 타겟 체인 (가장 긴
        폴리펩타이드) 만 식별하면 된다. raw PDB 에 흔한 병리적 케이스를 방어한다.

    처리:
        - NMR 멀티모델: MODEL 1 만 사용 (MODEL 2 이상 등장 시 break)
        - ATOM/HETATM 모두 스캔하되 CA 원자만 카운트 (표준 AA 판별)
        - HIS 변이형 (HID/HIE/HIP/HSD 등) 을 표준 AA 로 인정

    Returns:
        (chain_id, n_residues) 또는 (None, 0)
    """
    chains: dict[str, set] = defaultdict(set)
    in_first_model = True
    seen_model = False

    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            rec = line[0:6]
            if rec == "MODEL ":
                tokens = line.split()
                if len(tokens) >= 2:
                    if seen_model and tokens[1] != "1":
                        in_first_model = False
                    elif tokens[1] == "1":
                        seen_model = True
                        in_first_model = True
                continue
            if rec == "ENDMDL":
                if not in_first_model:
                    break
                # 첫 MODEL 블록의 ENDMDL 이후에는 이후 블록을 무시
                in_first_model = False
                continue
            if not in_first_model and seen_model:
                continue
            if rec != "ATOM  " and rec != "HETATM":
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            altloc = line[16:17]
            if altloc not in (" ", "A"):
                continue
            resname = line[17:20].strip().upper()
            if resname not in AA3_TO_AA1:
                continue
            chain_id = line[21:22]
            resnum = line[22:26].strip()
            icode  = line[26:27]
            chains[chain_id].add((resnum, icode))

    if not chains:
        return None, 0
    target_chain = max(chains, key=lambda c: len(chains[c]))
    return target_chain, len(chains[target_chain])


def identify_chains(pdb_path: str, peptide_sequence: str) -> dict:
    """복합체 PDB 에서 펩타이드 체인과 타겟 체인을 자동 식별한다.

    Returns:
        dict: {
            "peptide_chain": str,
            "target_chain": str,
            "peptide_len": int,
            "target_len": int,
            "similarity": float,  # 0..1
        }
    """
    chains = parse_pdb_chains(pdb_path)
    if not chains:
        raise ValueError(f"PDB 에서 체인 추출 실패: {pdb_path}")

    # 각 체인에 대해 펩타이드 서열과의 유사도 계산
    scored = []
    for cid, residues in chains.items():
        seq = chain_to_sequence(residues)
        sim = sequence_similarity(peptide_sequence, seq)
        scored.append((cid, len(residues), seq, sim))

    # 펩타이드 후보: 유사도 최고 + 길이가 펩타이드와 근접한 것
    # 펩타이드 사전 길이를 ±50% 범위에서 우선순위 매김
    plen = len(peptide_sequence)
    def _peptide_score(item):
        cid, length, seq, sim = item
        len_penalty = abs(length - plen) / max(plen, 1)
        return sim - 0.1 * len_penalty
    peptide_candidate = max(scored, key=_peptide_score)

    # 타겟 후보: 펩타이드 아닌 체인 중 가장 긴 것
    target_candidates = [s for s in scored if s[0] != peptide_candidate[0]]
    if target_candidates:
        target_candidate = max(target_candidates, key=lambda s: s[1])
    else:
        # 단일 체인만 존재 (비정상) → 자기 자신을 타겟으로 처리하고 경고
        target_candidate = peptide_candidate

    return {
        "peptide_chain": peptide_candidate[0],
        "target_chain":  target_candidate[0],
        "peptide_len":   peptide_candidate[1],
        "target_len":    target_candidate[1],
        "peptide_seq":   peptide_candidate[2],
        "similarity":    peptide_candidate[3],
    }


def extract_complex_pdb(src_pdb: str,
                        dst_pdb: str,
                        target_chain: str,
                        peptide_chain: str) -> None:
    """복합체 PDB 에서 target/peptide 두 체인만 추출하고 체인 ID 를 A/B 로 표준화한다.

    - ATOM 레코드만 유지 (물, 이온, 리간드, HETATM 제거)
    - target_chain → 'A', peptide_chain → 'B'
    - ALTLOC 은 첫 번째만 유지
    - TER 레코드 재구성
    """
    target_lines: list[str] = []
    peptide_lines: list[str] = []
    with open(src_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            rec = line[0:6]
            if rec not in ("ATOM  ", "HETATM"):
                continue
            # 히스티딘 변이형은 ATOM 으로 취급되지만 일부 파일은 HETATM 로 기록. 표준 AA 만 유지.
            resname = line[17:20].strip().upper()
            if resname not in AA3_TO_AA1:
                continue
            altloc = line[16:17]
            if altloc not in (" ", "A"):
                continue
            chain_id = line[21:22]
            new_line = None
            if chain_id == target_chain:
                # ATOM 으로 강제 + chain A
                new_line = "ATOM  " + line[6:21] + "A" + line[22:]
                target_lines.append(new_line)
            elif chain_id == peptide_chain:
                new_line = "ATOM  " + line[6:21] + "B" + line[22:]
                peptide_lines.append(new_line)

    if not target_lines:
        raise ValueError(f"타겟 체인 '{target_chain}' 원자가 {src_pdb} 에 없음")
    if not peptide_lines:
        raise ValueError(f"펩타이드 체인 '{peptide_chain}' 원자가 {src_pdb} 에 없음")

    os.makedirs(os.path.dirname(dst_pdb) or ".", exist_ok=True)
    with open(dst_pdb, "w", encoding="utf-8") as f:
        f.write("REMARK   1 UPDD_benchmark extracted complex\n")
        f.writelines(target_lines)
        f.write("TER\n")
        f.writelines(peptide_lines)
        f.write("TER\nEND\n")


def prepare_benchmark_input(record: dict,
                            src_pdb: str,
                            work_dir: str,
                            logger: logging.Logger) -> tuple[str, dict]:
    """단일 벤치마크 레코드의 MD 입력 디렉토리를 준비한다.

    Returns:
        (inputdir, chain_info)
    """
    os.makedirs(work_dir, exist_ok=True)
    chain_info = identify_chains(src_pdb, record["sequence"])
    logger.info(f"  [prep] {record['peptide_id']}: target={chain_info['target_chain']} "
                f"({chain_info['target_len']} aa), peptide={chain_info['peptide_chain']} "
                f"({chain_info['peptide_len']} aa), sim={chain_info['similarity']:.2f}")

    complex_pdb = os.path.join(work_dir, f"{record['peptide_id']}_complex.pdb")
    extract_complex_pdb(
        src_pdb=src_pdb,
        dst_pdb=complex_pdb,
        target_chain=chain_info["target_chain"],
        peptide_chain=chain_info["peptide_chain"],
    )
    return work_dir, chain_info


# ==========================================
# 6. run 서브커맨드 (배치 파이프라인)
# ==========================================
def _find_colabfold_python() -> str | None:
    """localcolabfold pixi 환경의 Python 바이너리를 탐색한다.

    ColabDesign (AfCycDesign) 은 JAX 가 설치된 pixi 환경에서만 실행 가능하다.
    conda 환경이 아니므로 별도 탐색이 필요하다.
    """
    candidates = glob.glob(os.path.join(PROJECT_DIR, "localcolabfold",
                                        ".pixi", "envs", "*", "bin", "python3"))
    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return None


def _pixi_run(script: str, extra_args: list[str], log_file: str) -> int:
    """localcolabfold pixi 환경의 Python 으로 스크립트를 실행한다.

    ColabDesign (AfCycDesign) 전용. conda 환경이 아닌 pixi 환경을 사용한다.
    """
    pybin = _find_colabfold_python()
    if pybin is None:
        with open(log_file, "a", encoding="utf-8") as lf:
            lf.write("[ERROR] localcolabfold pixi Python 을 찾을 수 없음\n")
        return 126
    cmd = [pybin, script] + extra_args
    with open(log_file, "a", encoding="utf-8") as lf:
        lf.write("\n" + "=" * 60 + "\n")
        lf.write(f"[CMD] {' '.join(cmd)}\n")
        lf.write("=" * 60 + "\n")
        lf.flush()
        try:
            proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
            return proc.returncode
        except FileNotFoundError:
            lf.write(f"[ERROR] Python 바이너리 실행 불가: {pybin}\n")
            return 127
        except KeyboardInterrupt:
            lf.write("[INTERRUPT] 사용자 중단\n")
            raise


def _conda_run(env_key: str, script: str, extra_args: list[str],
               log_file: str) -> int:
    """conda run -n <env> python <script> <args> 를 실행하고 로그를 기록한다."""
    env_name = ENV_NAMES.get(env_key, env_key)
    cmd = ["conda", "run", "-n", env_name, "--no-capture-output",
           "python", script] + extra_args
    with open(log_file, "a", encoding="utf-8") as lf:
        lf.write("\n" + "=" * 60 + "\n")
        lf.write(f"[CMD] {' '.join(cmd)}\n")
        lf.write("=" * 60 + "\n")
        lf.flush()
        try:
            proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
            return proc.returncode
        except FileNotFoundError:
            lf.write("[ERROR] conda 명령을 찾을 수 없음\n")
            return 127
        except KeyboardInterrupt:
            lf.write("[INTERRUPT] 사용자 중단\n")
            raise


def _run_md_stage(record: dict, inputdir: str, outputdir: str,
                  md_steps: int, log_file: str) -> int:
    os.makedirs(outputdir, exist_ok=True)
    topology = record["topology"]
    valid_topos = {"linear", "cyclic_htc", "cyclic_ss", "cyclic_nm", "bicyclic"}
    # [v41] AfCycDesign 성공 시 cyclic_htc/cyclic_nm 의 N-C < 3 Å 이 보장되므로
    # 원래 토폴로지를 전달. ColabFold 폴백 시 N-C > 30 Å 이므로 linear 로 전달.
    topo_arg = topology if topology in valid_topos else "linear"
    args = [
        "--inputdir",     inputdir,
        "--outputdir",    outputdir,
        "--steps",        str(md_steps),
        "--topology",     topo_arg,
        "--ncaa_label",   "none",
        "--binder_chain", "B",
        "--seed",         "42",
        "--platform",     "CUDA",
    ]
    return _conda_run("md", SCRIPT_PATHS["run_md"], args, log_file)


def _run_snapshot_stage(md_dir: str, snap_dir: str,
                        n_snapshots: int, log_file: str) -> int:
    os.makedirs(snap_dir, exist_ok=True)
    args = [
        "--md_dir",       md_dir,
        "--outputdir",    snap_dir,
        "--n_snapshots",  str(n_snapshots),
        "--binder_chain", "B",
    ]
    return _conda_run("md", SCRIPT_PATHS["extract_snap"], args, log_file)


def _run_qmmm_stage(snap_dir: str, qmmm_dir: str,
                    mode: str, log_file: str) -> int:
    os.makedirs(qmmm_dir, exist_ok=True)
    args = [
        "--snapdir",      snap_dir,
        "--outputdir",    qmmm_dir,
        "--mode",         mode,
        "--cutoff",       "5.0",
        "--qm_basis",     "6-31G*",
        "--qm_xc",        "wb97x-d3",
        "--ncaa_elem",    "none",
        "--binder_chain", "B",
    ]
    return _conda_run("qmmm", SCRIPT_PATHS["run_qmmm"], args, log_file)


def _run_mmgbsa_stage(md_dir: str, mmgbsa_dir: str, log_file: str) -> int:
    os.makedirs(mmgbsa_dir, exist_ok=True)
    args = [
        "--md_dir",        md_dir,
        "--outputdir",     mmgbsa_dir,
        "--ncaa_elem",     "none",
        "--receptor_chain", "A",
        "--binder_chain",   "B",
    ]
    return _conda_run("md", SCRIPT_PATHS["run_mmgbsa"], args, log_file)


# ==========================================
# 6b. 전처리 + AF2 헬퍼 (신규: CLAUDE.md AF2 통합 파이프라인)
# ==========================================
def find_colabfold_bin() -> str | None:
    """UPDD.py 의 find_colabfold_bin 과 동일 규약으로 colabfold_batch 를 탐색한다."""
    candidates: list[str] = []
    w = shutil.which("colabfold_batch")
    if w:
        candidates.append(w)
    candidates += glob.glob(os.path.join(PROJECT_DIR, "localcolabfold",
                                         ".pixi", "envs", "*", "bin", "colabfold_batch"))
    candidates += glob.glob(os.path.join(PROJECT_DIR, ".venv*", "bin", "colabfold_batch"))
    candidates += glob.glob(os.path.join(PROJECT_DIR, "venv*", "bin", "colabfold_batch"))
    candidates += glob.glob(os.path.expanduser("~/.pixi/envs/*/bin/colabfold_batch"))
    for base in (UPDD_DIR, HOME_DIR):
        candidates += glob.glob(os.path.join(base, ".pixi", "envs", "*", "bin", "colabfold_batch"))
    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return None


def _run_preprocess_target(pdb_path: str,
                           output_path: str,
                           target_chain: str,
                           log_file: str,
                           keep_hetatms: str = "") -> int:
    """preprocess_target.py 를 subprocess 로 호출하여 타겟 PDB 를 전처리한다.

    keep_hetatms 는 CSV 포맷 (세미콜론 구분, 예: "GDP;MG") 으로 받아서
    preprocess_target.py 의 CLI 규약 (쉼표 구분, 예: "GDP,MG") 으로 변환한다.
    세미콜론 구분을 사용하는 이유: CSV 내부 구분자 (,) 와 충돌 방지.

    Returns: returncode (0=성공)
    """
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    args = [
        "--input",       pdb_path,
        "--output",      output_path,
        "--chains",      target_chain,
        "--hetatm_mode", "auto",
    ]
    if keep_hetatms:
        args += ["--keep_hetatms", keep_hetatms.replace(";", ",")]
    return _conda_run("md", SCRIPT_PATHS["preprocess_target"], args, log_file)


def extract_target_sequence(preprocessed_pdb: str, target_chain: str = "A") -> str:
    """전처리된 PDB 에서 타겟 체인의 아미노산 서열을 1-letter code 로 추출한다.

    CA 원자의 잔기명을 3-letter → 1-letter 로 변환한다. 동일 (chain, resnum, icode)
    튜플은 한 번만 방문하여 중복 (ALTLOC) 을 방지한다. HIS 변이형 (HID/HIE/HIP/HSD/HSE/HSP)
    은 H 로 흡수된다.
    """
    residues: list[str] = []
    seen: set[tuple[str, str, str]] = set()
    with open(preprocessed_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            altloc = line[16:17]
            if altloc not in (" ", "A"):
                continue
            chain = line[21:22]
            if chain != target_chain:
                continue
            resnum  = line[22:26].strip()
            icode   = line[26:27]
            resname = line[17:20].strip().upper()
            key = (chain, resnum, icode)
            if key in seen:
                continue
            seen.add(key)
            aa = AA3_TO_AA1.get(resname, "X")
            residues.append(aa)
    return "".join(residues)


def prepare_af2_input(peptide_id: str,
                      target_sequence: str,
                      binder_sequence_raw: str,
                      af2_input_dir: str) -> str:
    """ColabFold 입력 CSV 를 작성한다.

    ColabFold 형식: `id,sequence` 헤더 + 복합체는 `target_seq:binder_seq` (콜론 구분).
    바인더 서열은 clean_peptide_sequence 로 cyclic 래퍼/인산화 접두사를 제거하여
    linear 1-letter 로 정규화한다.
    """
    os.makedirs(af2_input_dir, exist_ok=True)
    csv_path = os.path.join(af2_input_dir, f"{peptide_id}_input.csv")
    clean_binder = clean_peptide_sequence(binder_sequence_raw)
    if not clean_binder:
        raise ValueError(f"바인더 서열 정규화 실패: raw={binder_sequence_raw!r}")
    complex_seq = f"{target_sequence}:{clean_binder}"
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["id", "sequence"])
        w.writerow([peptide_id, complex_seq])
    return csv_path


def _run_colabfold(input_csv: str,
                   output_dir: str,
                   af2_bin: str,
                   log_file: str,
                   logger: logging.Logger | None = None) -> int:
    """ColabFold 로 복합체 구조를 예측한다.

    CLAUDE.md 에서 지정한 플래그 (`--num-recycle 3 --pair-mode unpaired_paired
    --msa-mode mmseqs2_uniref_env`) 를 사용하여 UPDD 본 파이프라인과 동일한
    설정을 재현한다. af2_bin 이 없으면 즉시 returncode=1 을 반환한다.
    """
    if not af2_bin or not os.path.isfile(af2_bin):
        if logger:
            logger.warning(f"[AF2] colabfold_batch 바이너리를 찾을 수 없음: {af2_bin}")
        return 1

    os.makedirs(output_dir, exist_ok=True)
    cmd = [
        af2_bin, input_csv, output_dir,
        "--num-recycle", "3",
        "--pair-mode",   "unpaired_paired",
        "--msa-mode",    "mmseqs2_uniref_env",
    ]
    with open(log_file, "a", encoding="utf-8") as lf:
        lf.write("\n" + "=" * 60 + "\n")
        lf.write(f"[CMD] {' '.join(cmd)}\n")
        lf.write("=" * 60 + "\n")
        lf.flush()
        try:
            proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
            return proc.returncode
        except FileNotFoundError:
            lf.write("[ERROR] colabfold_batch 실행 불가\n")
            return 127


def select_best_af2_structure(af2_output_dir: str,
                              peptide_id: str) -> tuple[str | None, float | None]:
    """AF2 출력에서 평균 pLDDT 가 가장 높은 PDB 를 선택한다.

    ColabFold 기본 파일명 규약:
        {id}_unrelaxed_rank_001_alphafold2_multimer_v3_model_{k}_seed_{s}.pdb
        {id}_scores_rank_001_alphafold2_multimer_v3_model_{k}_seed_{s}.json

    각 scores JSON 의 `plddt` 필드는 잔기별 리스트 또는 scalar. 평균값을 비교
    기준으로 삼는다. PDB 파일 경로는 scores 파일명에서 "_scores_" → "_unrelaxed_"
    치환으로 추론한다. 파일이 실제로 존재할 때만 후보로 취급한다.
    """
    score_files = glob.glob(os.path.join(af2_output_dir, "*_scores_*.json"))
    if not score_files:
        return None, None

    best_pdb: str | None = None
    best_plddt: float = -1.0

    for sf in score_files:
        try:
            with open(sf, "r", encoding="utf-8") as f:
                data = json.load(f)
        except (OSError, json.JSONDecodeError):
            continue
        plddt_raw = data.get("plddt")
        if isinstance(plddt_raw, list) and plddt_raw:
            plddt = sum(plddt_raw) / len(plddt_raw)
        elif isinstance(plddt_raw, (int, float)):
            plddt = float(plddt_raw)
        else:
            continue
        if plddt > best_plddt:
            pdb_name = os.path.basename(sf).replace("_scores_", "_unrelaxed_").replace(".json", ".pdb")
            pdb_path = os.path.join(af2_output_dir, pdb_name)
            if os.path.exists(pdb_path):
                best_plddt = plddt
                best_pdb = pdb_path

    return best_pdb, (best_plddt if best_pdb else None)


def run_single_benchmark(record: dict,
                         src_pdb: str,
                         output_root: str,
                         mode_preset: dict,
                         logger: logging.Logger,
                         preprocess_dir: str | None = None,
                         af2_input_dir: str | None = None,
                         af2_result_dir: str | None = None,
                         af2_bin: str | None = None,
                         use_experimental: bool = False) -> dict:
    """단일 펩타이드-타겟 쌍에 대해 preprocess → AF2 → MD → snapshot → QM/MM → MM-GBSA 실행.

    기본 (AF2 통합) 플로우:
        1. identify_target_chain: raw PDB 에서 가장 긴 단백질 체인 식별
        2. preprocess_target.py: 타겟 체인만 추출 + 물/이온 제거 + PDBFixer 복구
        3. extract_target_sequence: 전처리된 PDB 에서 타겟 1-letter 서열 추출
        4. prepare_af2_input: ColabFold 입력 CSV 작성 (target:binder)
        5. _run_colabfold: ColabFold 실행
        6. select_best_af2_structure: 최고 pLDDT PDB 선택
        7. MD → Snapshot → QM/MM → MM-GBSA (AF2 출력을 입력으로)

    --use-experimental 플로우 (하위 호환):
        기존 raw PDB 복합체 추출 경로를 그대로 사용한다. 스코어링 단독 평가용.

    에러 처리 원칙:
        개별 단계 실패 시 그 지점에서 중단하고 status/failed_stage/error_msg 를
        기록한다. 다음 레코드로 파이프라인을 계속 진행할 수 있도록 예외를
        여기서 잡고 결과 dict 만 반환한다.
    """
    pid = record["peptide_id"]
    pdb_id = record.get("pdb_id", "")
    work_root    = os.path.join(output_root, pid)
    md_dir       = os.path.join(work_root, "md")
    snap_dir     = os.path.join(work_root, "snapshots")
    qmmm_dir     = os.path.join(work_root, "qmmm")
    mmgbsa_dir   = os.path.join(work_root, "mmgbsa")
    md_input_dir = os.path.join(work_root, "md_input")
    log_file     = os.path.join(work_root, "run.log")
    os.makedirs(work_root, exist_ok=True)

    result: dict = {
        "peptide_id":          pid,
        "status":              "FAILED",
        "delta_g_kcal":        None,
        "delta_g_corrected":   None,
        "qmmm_interact_kcal":  None,
        "rmsd_mean":           None,
        "plddt":               None,
        "failed_stage":        None,
        "error_msg":           None,
        "chain_info":          None,
        "pipeline":            "experimental" if use_experimental else "af2",
    }

    # ── 헤더 ──
    kd = record.get("kd_nm")
    kd_str = f", Kd={kd:g} nM" if kd else ""
    bar = "━" * 55
    logger.info(bar)
    logger.info(f"  {pid} (pdb={pdb_id or '?'}, topo={record.get('topology', '?')}{kd_str})")
    logger.info(bar)
    total = _TOTAL_STAGES  # 6 고정 (experimental 에서 stage 2 는 skip)

    # ── [1/6] 타겟 전처리 ──
    _log_stage(logger, 1, total, "타겟 전처리", "start")
    t0 = time.time()
    if use_experimental:
        try:
            _, chain_info = prepare_benchmark_input(record, src_pdb, md_input_dir, logger)
            result["chain_info"] = chain_info
            _log_stage(logger, 1, total, "타겟 전처리", "ok",
                       elapsed=time.time() - t0,
                       detail=(f"chain {chain_info['target_chain']}→A, "
                               f"{chain_info.get('target_len','?')} aa + "
                               f"{chain_info['peptide_chain']}→B, "
                               f"{chain_info.get('peptide_len','?')} aa"))
        except Exception as e:
            _log_stage(logger, 1, total, "타겟 전처리", "fail",
                       elapsed=time.time() - t0, detail=str(e)[:60])
            result["failed_stage"] = "prepare"
            result["error_msg"] = f"{type(e).__name__}: {e}"
            return result
    else:
        if preprocess_dir is None or af2_input_dir is None or af2_result_dir is None:
            _log_stage(logger, 1, total, "타겟 전처리", "fail",
                       elapsed=time.time() - t0, detail="AF2 디렉토리 인자 누락")
            result["failed_stage"] = "config"
            result["error_msg"] = "AF2 디렉토리 인자 누락 (preprocess_dir/af2_input_dir/af2_result_dir)"
            return result

        target_chain, n_res = identify_target_chain(src_pdb)
        if target_chain is None:
            _log_stage(logger, 1, total, "타겟 전처리", "fail",
                       elapsed=time.time() - t0, detail="타겟 체인 식별 실패")
            result["failed_stage"] = "chain_id"
            result["error_msg"] = "타겟 체인 식별 실패"
            return result
        result["chain_info"] = {
            "target_chain": target_chain, "target_len": n_res,
            "peptide_chain": None, "peptide_len": None, "similarity": None,
        }

        preprocessed_pdb = os.path.join(preprocess_dir, f"{pdb_id}_clean.pdb")
        preprocess_reused = os.path.exists(preprocessed_pdb)
        if not preprocess_reused:
            os.makedirs(preprocess_dir, exist_ok=True)
            keep_het = (record.get("keep_hetatms") or "").strip()
            rc = _run_preprocess_target(
                pdb_path=src_pdb,
                output_path=preprocessed_pdb,
                target_chain=target_chain,
                log_file=log_file,
                keep_hetatms=keep_het,
            )
            if rc != 0 or not os.path.exists(preprocessed_pdb):
                _log_stage(logger, 1, total, "타겟 전처리", "fail",
                           elapsed=time.time() - t0, detail=f"preprocess rc={rc}")
                result["failed_stage"] = "preprocess"
                result["error_msg"] = f"preprocess_target.py 실패 (rc={rc})"
                return result

        target_seq = extract_target_sequence(preprocessed_pdb, target_chain)
        if not target_seq or len(target_seq) < 10:
            _log_stage(logger, 1, total, "타겟 전처리", "fail",
                       elapsed=time.time() - t0,
                       detail=f"서열 추출 실패 (len={len(target_seq)})")
            result["failed_stage"] = "seq_extract"
            result["error_msg"] = f"타겟 서열 추출 실패 (len={len(target_seq)})"
            return result

        reuse_tag = " [reuse]" if preprocess_reused else ""
        _log_stage(logger, 1, total, "타겟 전처리", "ok",
                   elapsed=time.time() - t0,
                   detail=f"chain {target_chain}, {n_res} aa, seq_len={len(target_seq)}{reuse_tag}")

    # ── [2/6] AF2 구조 예측 ──
    if use_experimental:
        _log_stage(logger, 2, total, "AF2 구조 예측", "skip",
                   detail="--use-experimental (실험 복합체 직접 사용)")
    else:
        _log_stage(logger, 2, total, "AF2 구조 예측", "start")
        t0 = time.time()
        af2_out = os.path.join(af2_result_dir, pid)
        has_prior = os.path.isdir(af2_out) and bool(
            glob.glob(os.path.join(af2_out, "*_scores_*.json"))
        )
        if has_prior:
            best_pdb, plddt = select_best_af2_structure(af2_out, pid)
            if best_pdb is None:
                _log_stage(logger, 2, total, "AF2 구조 예측", "fail",
                           elapsed=time.time() - t0, detail="기존 출력에서 유효 구조 없음")
                result["failed_stage"] = "af2_select"
                result["error_msg"] = "AF2 기존 출력에서 유효한 구조를 찾지 못함"
                return result
            result["plddt"] = plddt
            _log_stage(logger, 2, total, "AF2 구조 예측", "skip",
                       detail=f"기존 결과 재사용 (pLDDT={plddt:.1f})")
        else:
            try:
                input_csv = prepare_af2_input(pid, target_seq, record.get("sequence_raw", ""),
                                              af2_input_dir)
            except Exception as e:
                _log_stage(logger, 2, total, "AF2 구조 예측", "fail",
                           elapsed=time.time() - t0, detail=str(e)[:60])
                result["failed_stage"] = "af2_input"
                result["error_msg"] = f"{type(e).__name__}: {e}"
                return result
            os.makedirs(af2_out, exist_ok=True)
            topo = record.get("topology", "linear")
            # [v44] 모든 토폴로지에서 ColabFold 사용. cyclic_htc/cyclic_nm 은
            # ColabFold + geometric gap closure 로 N-C 말단을 폐합한다.
            rc = _run_colabfold(input_csv, af2_out, af2_bin, log_file, logger)
            if rc != 0:
                _log_stage(logger, 2, total, "AF2 구조 예측", "fail",
                           elapsed=time.time() - t0, detail=f"ColabFold rc={rc}")
                result["failed_stage"] = "af2"
                result["error_msg"] = f"ColabFold 실패 (rc={rc})"
                return result

            # cyclic_htc/cyclic_nm: N-C gap closure 적용
            if topo in ("cyclic_htc", "cyclic_nm"):
                sys.path.insert(0, UTILS_DIR)
                from run_cyclic_fold import close_nc_gap as _close_nc_gap  # noqa: E402
                for _gap_pdb in glob.glob(os.path.join(af2_out, "*_unrelaxed_*.pdb")):
                    _close_nc_gap(_gap_pdb, _gap_pdb, binder_chain="B")
                sys.path.pop(0)

            best_pdb, plddt = select_best_af2_structure(af2_out, pid)
            if best_pdb is None:
                _log_stage(logger, 2, total, "AF2 구조 예측", "fail",
                           elapsed=time.time() - t0, detail="출력에서 유효 구조 없음")
                result["failed_stage"] = "af2_select"
                result["error_msg"] = "AF2 출력에서 유효한 구조를 찾지 못함"
                return result
            result["plddt"] = plddt
            result["pipeline"] = "af2+gapclosure" if topo in ("cyclic_htc", "cyclic_nm") else "af2"
            _log_stage(logger, 2, total, "AF2 구조 예측", "ok",
                       elapsed=time.time() - t0, detail=f"pLDDT={plddt:.1f}")

        # AF2 출력을 MD 입력 디렉토리로 복사
        os.makedirs(md_input_dir, exist_ok=True)
        md_input_pdb = os.path.join(md_input_dir, f"{pid}_af2.pdb")
        try:
            shutil.copyfile(best_pdb, md_input_pdb)
        except OSError as e:
            _log_stage(logger, 3, total, "MD 시뮬레이션", "fail", detail=f"AF2 PDB 복사 실패: {e}")
            result["failed_stage"] = "md_stage"
            result["error_msg"] = f"AF2 PDB 복사 실패: {e}"
            return result

    # [v44] gap closure 로 N-C < 4 Å 이 보장되므로 원래 토폴로지를 그대로 전달.
    md_record = dict(record)

    # ── [3/6] MD 시뮬레이션 ──
    _log_stage(logger, 3, total, "MD 시뮬레이션", "start",
               detail=f"steps={mode_preset['md_steps']}, topo={md_record['topology']}")
    t0 = time.time()
    try:
        rc = _run_md_stage(md_record, md_input_dir, md_dir, mode_preset["md_steps"], log_file)
        if rc != 0:
            raise RuntimeError(f"run_restrained_md returncode={rc}")
        final_pdbs = glob.glob(os.path.join(md_dir, "*_final.pdb"))
        if not final_pdbs:
            raise RuntimeError("MD final PDB 없음")
    except Exception as e:
        _log_stage(logger, 3, total, "MD 시뮬레이션", "fail",
                   elapsed=time.time() - t0, detail=str(e)[:60])
        result["failed_stage"] = "md"
        result["error_msg"] = f"{type(e).__name__}: {e}"
        return result
    _log_stage(logger, 3, total, "MD 시뮬레이션", "ok",
               elapsed=time.time() - t0,
               detail=f"{mode_preset['md_steps']*2/1000:.0f} ps")

    # ── [4/6] 스냅샷 추출 ──
    _log_stage(logger, 4, total, "스냅샷 추출", "start")
    t0 = time.time()
    snap_ok = True
    try:
        rc = _run_snapshot_stage(md_dir, snap_dir, n_snapshots=5, log_file=log_file)
        if rc != 0:
            raise RuntimeError(f"extract_snapshots returncode={rc}")
        n_snaps = len(glob.glob(os.path.join(snap_dir, "*.pdb")))
        _log_stage(logger, 4, total, "스냅샷 추출", "ok",
                   elapsed=time.time() - t0, detail=f"{n_snaps} 개 추출")
    except Exception as e:
        snap_ok = False
        _log_stage(logger, 4, total, "스냅샷 추출", "fail",
                   elapsed=time.time() - t0, detail=str(e)[:60])
        result["failed_stage"] = "snapshot"
        result["error_msg"] = f"{type(e).__name__}: {e}"
        # 스냅샷 실패해도 MM-GBSA 는 final PDB 로 계속 시도한다

    # [v39 Fix 4] PBC unwrap 검증 — QM/MM 전에 스냅샷 근접도를 방어적으로 확인
    if snap_ok:
        _verify_snapshot_proximity(snap_dir, logger=logger)

    # ── [5/6] QM/MM 계산 ──
    if not mode_preset["run_qmmm"]:
        _log_stage(logger, 5, total, "QM/MM 계산", "skip", detail=f"mode={mode_preset}")
    elif not snap_ok or not glob.glob(os.path.join(snap_dir, "*.pdb")):
        _log_stage(logger, 5, total, "QM/MM 계산", "skip", detail="스냅샷 부재")
    else:
        _log_stage(logger, 5, total, "QM/MM 계산", "start",
                   detail=f"mode={mode_preset['qmmm_mode']}")
        t0 = time.time()
        try:
            rc = _run_qmmm_stage(snap_dir, qmmm_dir, mode_preset["qmmm_mode"], log_file)
            if rc != 0:
                raise RuntimeError(f"run_qmmm returncode={rc}")
            qmmm_summary = os.path.join(qmmm_dir, "qmmm_summary.json")
            ia_mean = None
            if os.path.exists(qmmm_summary):
                with open(qmmm_summary, "r", encoding="utf-8") as f:
                    data = json.load(f)
                ia = [r["interaction_kcal"] for r in data.get("results", [])
                      if r.get("interaction_kcal") is not None]
                if ia:
                    ia_mean = float(sum(ia) / len(ia))
                    result["qmmm_interact_kcal"] = ia_mean
            detail_q = f"E_int={ia_mean:+.2f} kcal/mol" if ia_mean is not None else ""
            _log_stage(logger, 5, total, "QM/MM 계산", "ok",
                       elapsed=time.time() - t0, detail=detail_q)
        except Exception as e:
            _log_stage(logger, 5, total, "QM/MM 계산", "fail",
                       elapsed=time.time() - t0, detail=str(e)[:60])
            result["failed_stage"] = "qmmm"
            result["error_msg"] = f"{type(e).__name__}: {e}"

    # ── [6/6] MM-GBSA ΔG ──
    if not mode_preset["run_mmgbsa"]:
        _log_stage(logger, 6, total, "MM-GBSA ΔG", "skip")
    else:
        _log_stage(logger, 6, total, "MM-GBSA ΔG", "start")
        t0 = time.time()
        try:
            rc = _run_mmgbsa_stage(md_dir, mmgbsa_dir, log_file)
            if rc != 0:
                raise RuntimeError(f"run_mmgbsa returncode={rc}")
            mmgbsa_summary = os.path.join(mmgbsa_dir, "mmgbsa_summary.json")
            if os.path.exists(mmgbsa_summary):
                with open(mmgbsa_summary, "r", encoding="utf-8") as f:
                    data = json.load(f)
                result["delta_g_kcal"] = data.get("mean_dg")
                ie_per = data.get("interaction_entropy_per_design", {}) or {}
                corrected = [v.get("delta_g_corrected_kcal") for v in ie_per.values()
                             if v and v.get("reliable") and v.get("delta_g_corrected_kcal") is not None]
                if corrected:
                    result["delta_g_corrected"] = float(sum(corrected) / len(corrected))
            dg = result["delta_g_kcal"]
            detail_g = f"ΔG={dg:+.2f} kcal/mol" if dg is not None else ""
            _log_stage(logger, 6, total, "MM-GBSA ΔG", "ok",
                       elapsed=time.time() - t0, detail=detail_g)
        except Exception as e:
            _log_stage(logger, 6, total, "MM-GBSA ΔG", "fail",
                       elapsed=time.time() - t0, detail=str(e)[:60])
            result["failed_stage"] = "mmgbsa"
            result["error_msg"] = f"{type(e).__name__}: {e}"
            return result

    # ── 최종 상태 ──
    if result["delta_g_kcal"] is not None:
        result["status"] = "SUCCESS" if result["failed_stage"] is None else "PARTIAL"
    else:
        result["status"] = "FAILED"
    return result


def load_progress(progress_file: str) -> dict:
    if not os.path.exists(progress_file):
        return {}
    try:
        with open(progress_file, "r", encoding="utf-8") as f:
            return json.load(f)
    except (OSError, json.JSONDecodeError):
        return {}


def save_progress(progress_file: str, progress: dict) -> None:
    tmp = progress_file + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(progress, f, indent=2)
    os.replace(tmp, progress_file)


def run_all_benchmarks(records: list[dict],
                       pdb_dir: str,
                       results_dir: str,
                       mode: str,
                       progress_file: str,
                       logger: logging.Logger,
                       preprocess_dir: str | None = None,
                       af2_input_dir: str | None = None,
                       af2_result_dir: str | None = None,
                       af2_bin: str | None = None,
                       use_experimental: bool = False) -> int:
    """모든 데이터 포인트에 대해 파이프라인을 실행한다 (resume 지원).

    PDB 파일이 pdb_dir 에 없는 레코드는 'pdb_missing' 으로 기록하고 건너뛴다.
    이미 progress.json 에서 SUCCESS/PARTIAL 상태인 레코드는 재활용한다.

    기본적으로 AF2 통합 플로우를 사용한다. --use-experimental 이 True 이면
    raw PDB 복합체 추출 경로를 사용한다 (스코어링 단독 평가용).

    Returns:
        성공 (SUCCESS/PARTIAL) 레코드 수
    """
    os.makedirs(results_dir, exist_ok=True)
    mode_preset = MODE_PRESETS[mode]
    pipeline_label = "experimental" if use_experimental else "af2"
    logger.info(f"[run] 모드: {mode} → {mode_preset}")
    logger.info(f"[run] 파이프라인: {pipeline_label}")
    logger.info(f"[run] 대상 레코드: {len(records)}")
    if not use_experimental and not af2_bin:
        logger.warning("[run] colabfold_batch 바이너리를 찾지 못함 — AF2 단계가 실패합니다")

    errors_csv = os.path.join(results_dir, "benchmark_errors.csv")
    progress   = load_progress(progress_file)

    error_rows: list[list[str]] = []
    benchmark_start = time.time()
    n_total_records = len(records)
    n_success = 0
    n_fail    = 0

    # 기존 progress 에서 이미 완료된 레코드 집계 (resume 시 진행률 표시용)
    for pid_done in list(progress.keys()):
        st = progress[pid_done].get("status")
        if st in ("SUCCESS", "PARTIAL"):
            # 다만 이 카운트는 이번 세션에서 실제로 처리한 게 아니라 재사용이므로
            # progress summary 에서는 "성공" 카운트에는 포함하지 않고 skip 으로 취급
            pass

    for idx, rec in enumerate(records, start=1):
        pid = rec["peptide_id"]
        if progress.get(pid, {}).get("status") in ("SUCCESS", "PARTIAL"):
            logger.info(f"[{idx:>3}/{n_total_records}] {pid}: progress 에서 재사용 → skip")
            _print_progress_summary(logger, idx, n_total_records, n_success, n_fail,
                                    time.time() - benchmark_start)
            continue

        if not rec.get("pdb_id"):
            msg = "데이터셋에 pdb_id 비어있음"
            logger.warning(f"[{idx:>3}/{n_total_records}] {pid}: {msg}")
            progress[pid] = {"status": "FAILED", "failed_stage": "pdb_missing", "error_msg": msg}
            error_rows.append([pid, "pdb_missing", msg])
            n_fail += 1
            save_progress(progress_file, progress)
            _print_progress_summary(logger, idx, n_total_records, n_success, n_fail,
                                    time.time() - benchmark_start)
            continue

        src_pdb = os.path.join(pdb_dir, f"{rec['pdb_id']}.pdb")
        if not os.path.exists(src_pdb):
            msg = f"PDB 파일 없음: {src_pdb}"
            logger.warning(f"[{idx:>3}/{n_total_records}] {pid}: {msg}")
            progress[pid] = {"status": "FAILED", "failed_stage": "pdb_missing", "error_msg": msg}
            error_rows.append([pid, "pdb_missing", msg])
            n_fail += 1
            save_progress(progress_file, progress)
            _print_progress_summary(logger, idx, n_total_records, n_success, n_fail,
                                    time.time() - benchmark_start)
            continue

        logger.info(f"\n[{idx:>3}/{n_total_records}]")
        t0 = time.time()
        try:
            res = run_single_benchmark(
                rec, src_pdb, results_dir, mode_preset, logger,
                preprocess_dir=preprocess_dir,
                af2_input_dir=af2_input_dir,
                af2_result_dir=af2_result_dir,
                af2_bin=af2_bin,
                use_experimental=use_experimental,
            )
        except KeyboardInterrupt:
            logger.warning("[run] 사용자 중단 — progress 저장 후 종료")
            save_progress(progress_file, progress)
            raise
        except Exception as e:  # 방어적 최종 캐치
            res = {
                "peptide_id": pid, "status": "FAILED",
                "failed_stage": "unknown", "error_msg": f"{type(e).__name__}: {e}",
            }

        elapsed = time.time() - t0
        dg = res.get("delta_g_kcal")
        dg_str = f"{dg:+.2f} kcal/mol" if dg is not None else "N/A"
        status_icon = "✅" if res["status"] in ("SUCCESS", "PARTIAL") else "❌"
        logger.info(f"\n  → {status_icon} {res['status']}  ΔG={dg_str}  elapsed={elapsed:.1f}s")

        progress[pid] = {
            "status":             res["status"],
            "delta_g_kcal":       res.get("delta_g_kcal"),
            "delta_g_corrected":  res.get("delta_g_corrected"),
            "qmmm_interact_kcal": res.get("qmmm_interact_kcal"),
            "plddt":              res.get("plddt"),
            "pipeline":           res.get("pipeline"),
            "failed_stage":       res.get("failed_stage"),
            "error_msg":          res.get("error_msg"),
            "elapsed_s":          elapsed,
        }
        if res["status"] in ("SUCCESS", "PARTIAL"):
            n_success += 1
        else:
            n_fail += 1
            error_rows.append([pid, res.get("failed_stage") or "", res.get("error_msg") or ""])
        save_progress(progress_file, progress)
        _print_progress_summary(logger, idx, n_total_records, n_success, n_fail,
                                time.time() - benchmark_start)

    if error_rows:
        with open(errors_csv, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["peptide_id", "failed_stage", "error_msg"])
            w.writerows(error_rows)
        logger.info(f"[run] 실패 목록: {errors_csv}")

    _print_final_summary(logger, progress, time.time() - benchmark_start)

    succ = sum(1 for v in progress.values() if v.get("status") in ("SUCCESS", "PARTIAL"))
    return succ


# ==========================================
# 7. evaluate 서브커맨드
# ==========================================
def _import_sci():
    """scipy/numpy/sklearn/matplotlib 을 지연 import 한다."""
    import numpy as _np
    from scipy import stats as _stats
    try:
        from sklearn.metrics import roc_auc_score, roc_curve
        has_sklearn = True
    except ImportError:
        roc_auc_score, roc_curve = None, None
        has_sklearn = False
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as _plt
        has_mpl = True
    except ImportError:
        _plt = None
        has_mpl = False
    return {
        "np": _np, "stats": _stats,
        "roc_auc_score": roc_auc_score, "roc_curve": roc_curve,
        "plt": _plt, "has_sklearn": has_sklearn, "has_mpl": has_mpl,
    }


def _kd_to_pkd(kd_nm: float) -> float:
    """nM 단위 Kd → pKd ( -log10(Kd_M) ). 입력 보호."""
    if kd_nm is None or kd_nm <= 0:
        return float("nan")
    return -math.log10(kd_nm * 1e-9)


def collect_results(records: list[dict], results_dir: str) -> list[dict]:
    """각 레코드에 대해 progress/mmgbsa/qmmm summary 를 매칭하여 평가용 레코드로 만든다."""
    progress_file = os.path.join(results_dir, "benchmark_progress.json")
    progress = load_progress(progress_file)

    merged: list[dict] = []
    for r in records:
        pid = r["peptide_id"]
        entry = dict(r)
        p = progress.get(pid, {})
        entry["status"]             = p.get("status", "MISSING")
        entry["delta_g_kcal"]       = p.get("delta_g_kcal")
        entry["delta_g_corrected"]  = p.get("delta_g_corrected")
        entry["qmmm_interact_kcal"] = p.get("qmmm_interact_kcal")
        entry["failed_stage"]       = p.get("failed_stage")

        # summary JSON 에서 보강 (progress 가 축약일 수 있음)
        mmgbsa_summary = os.path.join(results_dir, pid, "mmgbsa", "mmgbsa_summary.json")
        if os.path.exists(mmgbsa_summary) and entry["delta_g_kcal"] is None:
            try:
                with open(mmgbsa_summary, "r", encoding="utf-8") as f:
                    data = json.load(f)
                entry["delta_g_kcal"] = data.get("mean_dg")
            except (OSError, json.JSONDecodeError):
                pass

        entry["pkd"] = _kd_to_pkd(r["kd_nm"]) if r["kd_nm"] else float("nan")
        merged.append(entry)
    return merged


def compute_metrics(sci, experimental_pkd, computed_scores, labels_binary):
    """핵심 평가 메트릭 계산.

    Args:
        experimental_pkd: 실험 pKd 배열
        computed_scores:  UPDD ΔG (낮을수록 강한 결합)
        labels_binary:    1=binder, 0=non-binder

    Note:
        - Spearman/Kendall 는 pKd (높을수록 강함) vs ΔG (낮을수록 강함) 이므로
          **반대 방향** 관계가 정답이다. 음수 부호는 붙이지 않고 상관계수의
          부호 자체로 해석한다: 이상적으로 ρ<0 (그리고 |ρ| 가 크다).
          본 함수는 그 원 부호를 유지하고, 리포트 단계에서 "rank_alignment"
          지표로 -ρ 를 함께 출력한다.
        - ROC-AUC 는 predict 방향을 맞추기 위해 -score 를 사용한다.
    """
    np = sci["np"]
    stats = sci["stats"]
    results: dict = {}

    exp = np.asarray(experimental_pkd, dtype=float)
    sco = np.asarray(computed_scores,  dtype=float)
    mask = np.isfinite(exp) & np.isfinite(sco)
    exp_m = exp[mask]
    sco_m = sco[mask]
    results["n_valid"] = int(mask.sum())

    if results["n_valid"] < 3:
        results["error"] = "유효 샘플이 3 개 미만"
        return results

    rho, p_spearman = stats.spearmanr(exp_m, sco_m)
    results["spearman_rho"]     = float(rho) if rho is not None else float("nan")
    results["spearman_p"]       = float(p_spearman) if p_spearman is not None else float("nan")
    results["rank_alignment"]   = -results["spearman_rho"]  # 높을수록 좋음

    tau, p_kendall = stats.kendalltau(exp_m, sco_m)
    results["kendall_tau"] = float(tau) if tau is not None else float("nan")
    results["kendall_p"]   = float(p_kendall) if p_kendall is not None else float("nan")

    r_val, p_pearson = stats.pearsonr(exp_m, sco_m)
    results["pearson_r"] = float(r_val)
    results["pearson_p"] = float(p_pearson)

    # 분류 메트릭
    lab = np.asarray(labels_binary, dtype=float)
    lab_mask = mask & np.isfinite(lab)
    lab_m = lab[lab_mask].astype(int)
    sco_l = sco[lab_mask]
    if sci["has_sklearn"] and lab_mask.sum() >= 4 and len(set(lab_m.tolist())) > 1:
        auc = sci["roc_auc_score"](lab_m, -sco_l)
        results["roc_auc"] = float(auc)
        fpr, tpr, _ = sci["roc_curve"](lab_m, -sco_l)
        results["roc_curve"] = {"fpr": fpr.tolist(), "tpr": tpr.tolist()}
    else:
        results["roc_auc"] = None

    # Enrichment Factor — 점수 (ΔG) 가 낮을수록 좋음
    for frac in (0.01, 0.05, 0.10):
        key = f"ef_{int(frac*100)}pct"
        results[key] = _enrichment_factor(lab_m.tolist() if len(lab_m) else [],
                                          sco_l.tolist() if len(sco_l) else [],
                                          frac)
    return results


def _enrichment_factor(labels: list[int], scores: list[float], top_fraction: float) -> float:
    n = len(labels)
    if n == 0:
        return 0.0
    n_top = max(1, int(n * top_fraction))
    n_active_total = sum(labels)
    if n_active_total == 0:
        return 0.0
    order = sorted(range(n), key=lambda i: scores[i])
    n_active_in_top = sum(labels[i] for i in order[:n_top])
    expected = n_active_total / n * n_top
    if expected == 0:
        return 0.0
    return n_active_in_top / expected


def run_kfold_cv(sci, records_valid: list[dict], n_folds: int = 5, seed: int = 42) -> dict:
    """K-Fold CV — 각 fold 에서 Spearman ρ 와 AUC 의 평균/표준편차를 산출."""
    import random
    np = sci["np"]
    n = len(records_valid)
    if n < n_folds * 2:
        return {"error": f"샘플 수 부족 (n={n}, folds={n_folds})"}

    rng = random.Random(seed)
    idx = list(range(n))
    rng.shuffle(idx)
    fold_size = n // n_folds

    rhos, aucs = [], []
    for f in range(n_folds):
        test_idx = set(idx[f * fold_size:(f + 1) * fold_size]) if f < n_folds - 1 \
                   else set(idx[f * fold_size:])
        fold_records = [records_valid[i] for i in test_idx]
        pkd = [r["pkd"]          for r in fold_records]
        sco = [r["delta_g_kcal"] for r in fold_records]
        lab = [r.get("label") if r.get("label") is not None else -1 for r in fold_records]
        m = compute_metrics(sci, pkd, sco, lab)
        if "spearman_rho" in m and math.isfinite(m["spearman_rho"]):
            rhos.append(m["spearman_rho"])
        if m.get("roc_auc") is not None:
            aucs.append(m["roc_auc"])

    out = {
        "n_folds": n_folds,
        "spearman_rho_mean": float(np.mean(rhos)) if rhos else float("nan"),
        "spearman_rho_std":  float(np.std(rhos))  if rhos else float("nan"),
        "roc_auc_mean":      float(np.mean(aucs)) if aucs else None,
        "roc_auc_std":       float(np.std(aucs))  if aucs else None,
    }
    return out


def run_loto_cv(sci, records_valid: list[dict]) -> dict:
    """Leave-One-Target-Out — 새로운 타겟에 대한 일반화 능력.

    타겟별로 한 그룹씩 제외하고 나머지 N-1 그룹의 Spearman ρ 를 계산.
    (기본적으로 제외 그룹 자체의 ρ 를 리포트한다.)
    """
    np = sci["np"]
    by_target: dict[str, list[dict]] = defaultdict(list)
    for r in records_valid:
        by_target[r["target_name"] or "UNKNOWN"].append(r)

    per_target: dict[str, dict] = {}
    for target, items in by_target.items():
        if len(items) < 3:
            continue
        pkd = [r["pkd"]          for r in items]
        sco = [r["delta_g_kcal"] for r in items]
        lab = [r.get("label") if r.get("label") is not None else -1 for r in items]
        m = compute_metrics(sci, pkd, sco, lab)
        per_target[target] = {
            "n":            len(items),
            "spearman_rho": m.get("spearman_rho"),
            "roc_auc":      m.get("roc_auc"),
        }

    valid_rho = [v["spearman_rho"] for v in per_target.values()
                 if v["spearman_rho"] is not None and math.isfinite(v["spearman_rho"])]
    return {
        "per_target": per_target,
        "mean_spearman_rho": float(np.mean(valid_rho)) if valid_rho else float("nan"),
        "n_targets": len(per_target),
    }


def run_subgroup_analysis(sci, records_valid: list[dict]) -> dict:
    """토폴로지별 하위그룹 메트릭."""
    out: dict[str, dict] = {}
    by_topo: dict[str, list[dict]] = defaultdict(list)
    for r in records_valid:
        by_topo[r["topology"]].append(r)
    for topo, items in by_topo.items():
        if len(items) < 3:
            out[topo] = {"n": len(items), "note": "샘플 부족"}
            continue
        pkd = [r["pkd"]          for r in items]
        sco = [r["delta_g_kcal"] for r in items]
        lab = [r.get("label") if r.get("label") is not None else -1 for r in items]
        out[topo] = compute_metrics(sci, pkd, sco, lab)
        out[topo]["n"] = len(items)
    return out


def generate_figures(sci, records_valid, metrics, subgroup, figures_dir):
    if not sci["has_mpl"]:
        return []
    os.makedirs(figures_dir, exist_ok=True)
    plt = sci["plt"]
    np  = sci["np"]
    produced = []

    # 1. 산점도
    pkd = np.array([r["pkd"]          for r in records_valid], dtype=float)
    sco = np.array([r["delta_g_kcal"] for r in records_valid], dtype=float)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(pkd, sco, alpha=0.7)
    ax.set_xlabel("Experimental pKd")
    ax.set_ylabel("Computed ΔG (kcal/mol)")
    ax.set_title(f"UPDD benchmark (ρ={metrics.get('spearman_rho', float('nan')):.2f})")
    ax.grid(True, alpha=0.3)
    p1 = os.path.join(figures_dir, "scatter_dg.png")
    fig.tight_layout()
    fig.savefig(p1, dpi=150)
    plt.close(fig)
    produced.append(p1)

    # 2. ROC curve
    if "roc_curve" in metrics and metrics["roc_curve"]:
        fig, ax = plt.subplots(figsize=(5, 5))
        fpr = metrics["roc_curve"]["fpr"]
        tpr = metrics["roc_curve"]["tpr"]
        ax.plot(fpr, tpr, label=f"AUC = {metrics.get('roc_auc', 0):.2f}")
        ax.plot([0, 1], [0, 1], "k--", alpha=0.5)
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.set_title("ROC Curve")
        ax.legend()
        ax.grid(True, alpha=0.3)
        p2 = os.path.join(figures_dir, "roc_curve.png")
        fig.tight_layout()
        fig.savefig(p2, dpi=150)
        plt.close(fig)
        produced.append(p2)

    # 3. Topology comparison
    topos = [t for t in subgroup if subgroup[t].get("spearman_rho") is not None]
    if topos:
        rhos = [subgroup[t]["spearman_rho"] for t in topos]
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.bar(topos, rhos)
        ax.axhline(0, color="k", linewidth=0.5)
        ax.set_ylabel("Spearman ρ")
        ax.set_title("Topology comparison")
        p3 = os.path.join(figures_dir, "topology_comparison.png")
        fig.tight_layout()
        fig.savefig(p3, dpi=150)
        plt.close(fig)
        produced.append(p3)

    return produced


def _tier_label(rho: float) -> str:
    if rho is None or not math.isfinite(rho):
        return "N/A"
    absrho = abs(rho)
    if absrho > 0.75:
        return "Tier 3 (우수)"
    if absrho > 0.60:
        return "Tier 2 (양호)"
    if absrho > 0.40:
        return "Tier 1 (유의)"
    return "Tier 0 (미달)"


def write_report(report_dir: str, records_all: list[dict], records_valid: list[dict],
                 metrics: dict, cv_results: dict, loto_results: dict,
                 subgroup: dict, figures: list[str],
                 pipeline_label: str = "af2") -> str:
    os.makedirs(report_dir, exist_ok=True)
    report_md = os.path.join(report_dir, "benchmark_report.md")

    n_total   = len(records_all)
    n_valid   = len(records_valid)
    n_failed  = sum(1 for r in records_all if r.get("status") == "FAILED")
    n_missing = sum(1 for r in records_all if r.get("status") == "MISSING")

    rho = metrics.get("spearman_rho", float("nan"))
    tier = _tier_label(rho)

    lines = []
    lines.append("# UPDD Benchmark Report")
    lines.append("")
    lines.append("## 벤치마크 방법론")
    lines.append("")
    if pipeline_label == "af2":
        lines.append("본 벤치마크는 UPDD 파이프라인의 **통합 성능** (구조 예측 + 스코어링) 을 평가한다.")
        lines.append("")
        lines.append("파이프라인: 타겟 PDB 전처리 → ColabFold 복합체 예측 → MD → QM/MM → MM-GBSA")
        lines.append("")
        lines.append("실험 구조 (PDB) 는 **타겟 체인 전처리에만** 사용되며, 복합체 구조는 ColabFold 가 "
                     "데이터셋의 바인더 서열로부터 예측한다. 따라서 본 평가에는 AF2 구조 예측 오차가 "
                     "포함되어 있으며, 스코어링 단독 성능과는 구분된다. 스코어링 단독 성능을 측정하려면 "
                     "`--use-experimental` 플래그를 사용한다.")
    else:
        lines.append("본 벤치마크는 **스코어링 단독 성능** 을 평가한다 (`--use-experimental`).")
        lines.append("")
        lines.append("파이프라인: 실험 복합체 구조 추출 → MD → QM/MM → MM-GBSA")
        lines.append("")
        lines.append("복합체 구조는 RCSB PDB 에서 직접 사용되며 AF2 단계를 건너뛴다. 따라서 본 평가는 "
                     "실제 운영 환경 (RFdiffusion → MPNN → AF2 → MD) 에서의 성능과는 차이가 있을 수 있다.")
    lines.append("")
    # [v40 FIX] cyclic_htc 방법론 주의사항
    lines.append("### cyclic_htc 토폴로지 주의사항")
    lines.append("")
    lines.append("데이터셋의 cyclic_htc 펩타이드는 ColabFold 로 구조 예측 후 **linear 토폴로지로 "
                 "MD 를 수행** 하였다. ColabFold 는 cyclic constraint 를 강제하지 않으므로, 예측 "
                 "구조는 linear 이며, 이 구조에서 head-to-tail cyclic bond 를 post-hoc 으로 형성하는 "
                 "것은 비물리적이다 (N-C 거리 30+ Å). cyclic_htc 하위그룹의 스코어링 결과는 "
                 '"linear 구조 기반 스코어링"으로 해석해야 하며, cyclic 구조 기반 스코어링은 '
                 "AfCycDesign 통합 후 가능하다. cyclic_ss (disulfide) 는 AF2 구조에서도 자연 "
                 "형성되므로 원래 토폴로지로 MD 를 수행한다.")
    lines.append("")
    lines.append(f"- Dataset 레코드: {n_total}")
    lines.append(f"- 평가 가능 (valid): {n_valid}")
    lines.append(f"- 실패: {n_failed}")
    lines.append(f"- 누락 (진행 불가): {n_missing}")
    lines.append("")
    lines.append("## Overall Performance")
    lines.append("")
    lines.append(f"| 지표 | 값 |")
    lines.append(f"|------|----|")
    lines.append(f"| Spearman ρ | {rho:+.3f} (p={metrics.get('spearman_p', float('nan')):.3g}) |")
    lines.append(f"| Kendall τ  | {metrics.get('kendall_tau', float('nan')):+.3f} |")
    lines.append(f"| Pearson r  | {metrics.get('pearson_r', float('nan')):+.3f} |")
    auc = metrics.get("roc_auc")
    lines.append(f"| ROC-AUC    | {auc if auc is None else f'{auc:.3f}'} |")
    lines.append(f"| EF 1%      | {metrics.get('ef_1pct', 0):.2f} |")
    lines.append(f"| EF 5%      | {metrics.get('ef_5pct', 0):.2f} |")
    lines.append(f"| EF 10%     | {metrics.get('ef_10pct', 0):.2f} |")
    lines.append(f"| **Tier**   | **{tier}** |")
    lines.append("")
    lines.append("해석: 핵심 지표는 Spearman ρ 이며, ΔG 가 낮을수록 pKd 가 높아야 하므로 "
                 "이상적으로 ρ<0 이다. |ρ|>0.4 = Tier1, |ρ|>0.6 = Tier2, |ρ|>0.75 = Tier3.")
    lines.append("")

    lines.append("## Cross-Validation")
    lines.append("")
    if "error" not in cv_results:
        lines.append(f"- {cv_results['n_folds']}-Fold Spearman ρ: "
                     f"{cv_results['spearman_rho_mean']:+.3f} ± {cv_results['spearman_rho_std']:.3f}")
        if cv_results.get("roc_auc_mean") is not None:
            lines.append(f"- {cv_results['n_folds']}-Fold ROC-AUC: "
                         f"{cv_results['roc_auc_mean']:.3f} ± {cv_results['roc_auc_std']:.3f}")
    else:
        lines.append(f"- (스킵: {cv_results['error']})")
    lines.append("")

    lines.append("## Leave-One-Target-Out")
    lines.append("")
    lines.append(f"- 타겟 수: {loto_results['n_targets']}")
    lines.append(f"- Spearman ρ 평균: {loto_results['mean_spearman_rho']:+.3f}")
    lines.append("")
    if loto_results["per_target"]:
        lines.append("| Target | n | Spearman ρ | ROC-AUC |")
        lines.append("|--------|---|------------|---------|")
        for t, v in sorted(loto_results["per_target"].items()):
            rho_v = v["spearman_rho"]
            auc_v = v["roc_auc"]
            lines.append(f"| {t} | {v['n']} | "
                         f"{rho_v:+.3f} | "
                         f"{auc_v if auc_v is None else f'{auc_v:.3f}'} |")
        lines.append("")

    lines.append("## Topology Subgroup")
    lines.append("")
    if subgroup:
        lines.append("| Topology | n | Spearman ρ | ROC-AUC | EF 10% |")
        lines.append("|----------|---|------------|---------|--------|")
        for t, m in sorted(subgroup.items()):
            if "spearman_rho" not in m:
                lines.append(f"| {t} | {m.get('n', 0)} | — | — | — |")
                continue
            auc_v = m.get("roc_auc")
            lines.append(f"| {t} | {m['n']} | {m.get('spearman_rho', 0):+.3f} | "
                         f"{auc_v if auc_v is None else f'{auc_v:.3f}'} | "
                         f"{m.get('ef_10pct', 0):.2f} |")
        lines.append("")

    if figures:
        lines.append("## Figures")
        lines.append("")
        for p in figures:
            rel = os.path.relpath(p, report_dir)
            lines.append(f"![{os.path.basename(p)}]({rel})")
        lines.append("")

    lines.append("## Failed / Missing Cases")
    lines.append("")
    bad = [r for r in records_all if r.get("status") not in ("SUCCESS", "PARTIAL")]
    if bad:
        lines.append(f"{len(bad)} 개의 데이터 포인트가 유효 점수를 산출하지 못했습니다.")
        lines.append("")
        lines.append("| peptide_id | status | stage |")
        lines.append("|------------|--------|-------|")
        for r in bad[:30]:
            lines.append(f"| {r['peptide_id']} | {r['status']} | {r.get('failed_stage') or ''} |")
        if len(bad) > 30:
            lines.append(f"| … | … | … (+{len(bad) - 30}) |")

    with open(report_md, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    # CSV: 상세 메트릭
    detailed_csv = os.path.join(report_dir, "detailed_metrics.csv")
    with open(detailed_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["peptide_id", "target_name", "topology", "kd_nm", "pkd",
                    "label", "delta_g_kcal", "delta_g_corrected", "qmmm_interact_kcal",
                    "status", "failed_stage"])
        for r in records_all:
            w.writerow([
                r["peptide_id"], r["target_name"], r["topology"],
                r["kd_nm"] if r["kd_nm"] is not None else "",
                f"{r['pkd']:.3f}" if r.get("pkd") and math.isfinite(r["pkd"]) else "",
                r.get("label") if r.get("label") is not None else "",
                r.get("delta_g_kcal") if r.get("delta_g_kcal") is not None else "",
                r.get("delta_g_corrected") if r.get("delta_g_corrected") is not None else "",
                r.get("qmmm_interact_kcal") if r.get("qmmm_interact_kcal") is not None else "",
                r.get("status") or "",
                r.get("failed_stage") or "",
            ])

    failed_csv = os.path.join(report_dir, "failed_cases.csv")
    with open(failed_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["peptide_id", "status", "failed_stage"])
        for r in records_all:
            if r.get("status") not in ("SUCCESS", "PARTIAL"):
                w.writerow([r["peptide_id"], r.get("status") or "", r.get("failed_stage") or ""])

    return report_md


def evaluate_and_report(records: list[dict],
                        results_dir: str,
                        report_dir: str,
                        figures_dir: str,
                        logger: logging.Logger,
                        pipeline_label: str = "af2") -> int:
    """파이프라인 결과를 수집해 메트릭을 계산하고 리포트를 생성한다.

    Returns:
        0 = 성공, 2 = 유효 결과 없음
    """
    sci = _import_sci()

    merged = collect_results(records, results_dir)

    # 유효한 평가 대상: ΔG 산출 + 실험 pKd 존재
    valid = [r for r in merged
             if r.get("delta_g_kcal") is not None
             and math.isfinite(r.get("pkd", float("nan")))]
    logger.info(f"[evaluate] 전체={len(merged)}, 유효={len(valid)}")

    if not valid:
        logger.error("[evaluate] 유효한 결과가 없음 — run 단계가 완료되었는지 확인")
        return 2

    pkd = [r["pkd"]          for r in valid]
    sco = [r["delta_g_kcal"] for r in valid]
    lab = [r.get("label") if r.get("label") is not None else -1 for r in valid]

    metrics      = compute_metrics(sci, pkd, sco, lab)
    cv_results   = run_kfold_cv(sci, valid)
    loto_results = run_loto_cv(sci, valid)
    subgroup     = run_subgroup_analysis(sci, valid)

    os.makedirs(figures_dir, exist_ok=True)
    figures = generate_figures(sci, valid, metrics, subgroup, figures_dir)
    report_md = write_report(
        report_dir, merged, valid,
        metrics, cv_results, loto_results, subgroup, figures,
        pipeline_label=pipeline_label,
    )
    logger.info(f"[evaluate] 리포트 작성: {report_md}")

    # 콘솔 요약
    print("\n=== UPDD Benchmark Results ===")
    print(f"Overall (n={metrics.get('n_valid')}):")
    print(f"  Spearman ρ = {metrics.get('spearman_rho', float('nan')):+.3f} "
          f"(p={metrics.get('spearman_p', float('nan')):.3g})")
    print(f"  Kendall τ  = {metrics.get('kendall_tau',  float('nan')):+.3f}")
    print(f"  Pearson r  = {metrics.get('pearson_r',    float('nan')):+.3f}")
    auc = metrics.get("roc_auc")
    if auc is not None:
        print(f"  ROC-AUC    = {auc:.3f}")
    print(f"  EF 1/5/10% = {metrics.get('ef_1pct', 0):.2f} / "
          f"{metrics.get('ef_5pct', 0):.2f} / {metrics.get('ef_10pct', 0):.2f}")
    print(f"Tier: {_tier_label(metrics.get('spearman_rho'))}")
    return 0


# ==========================================
# 8. CLI
# ==========================================
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="UPDD_benchmark.py",
        description="UPDD 파이프라인 벤치마킹 — 실험 구조 기반 스코어링 검증",
    )
    p.add_argument(
        "--mode",
        choices=list(MODE_PRESETS.keys()),
        default="standard",
        help="벤치마크 모드: quick(10분), standard(수시간), full(수일). 기본: standard",
    )
    p.add_argument(
        "--topology",
        choices=["all", "linear", "cyclic_htc", "cyclic_ss"],
        default="all",
        help="평가할 토폴로지 필터. 기본: all",
    )
    p.add_argument(
        "--max-samples",
        type=int,
        default=0,
        help="최대 처리 샘플 수 (0=전체). 디버깅용.",
    )
    p.add_argument(
        "--skip-download",
        action="store_true",
        help="PDB 다운로드 단계를 건너뛴다 (이미 다운로드된 경우).",
    )
    p.add_argument(
        "--skip-run",
        action="store_true",
        help="파이프라인 실행 단계를 건너뛴다 (기존 결과로 평가만 수행).",
    )
    p.add_argument(
        "--use-experimental",
        action="store_true",
        help="AF2 대신 실험 PDB 복합체 구조를 직접 사용한다 (스코어링 단독 평가용). "
             "체인 식별 실패 시 해당 데이터 포인트를 건너뛴다.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    # 디렉토리 자동 생성 (AF2 파이프라인 포함)
    for d in (BENCHMARK_DIR, PDB_DIR, PREPROCESS_DIR, AF2_INPUT_DIR, AF2_RESULT_DIR,
              RESULTS_DIR, REPORT_DIR, FIGURES_DIR):
        os.makedirs(d, exist_ok=True)

    # 데이터셋 존재 확인
    if not os.path.exists(DATASET_CSV):
        print(f"[!] 데이터셋을 찾을 수 없습니다: {DATASET_CSV}")
        print(f"    UPDD_benchmark_dataset.csv 를 {BENCHMARK_DIR}/ 에 배치하세요.")
        return 1

    logger = setup_benchmark_logger(BENCHMARK_DIR)
    pipeline_label = "experimental" if args.use_experimental else "af2"
    logger.info(f"=== UPDD Benchmark 시작 (mode={args.mode}, topology={args.topology}, "
                f"pipeline={pipeline_label}) ===")

    # 1. 데이터셋 로드 + 필터링
    records = load_dataset(DATASET_CSV)
    if args.topology != "all":
        records = [r for r in records if r["topology"] == args.topology]
    if args.max_samples > 0:
        records = records[:args.max_samples]
    assign_labels(records)
    logger.info(f"데이터셋: {len(records)} 개 로드 (topology={args.topology})")

    # 2. PDB 다운로드
    if not args.skip_download:
        n_dl, n_sk, n_fl = download_all_pdbs(records, PDB_DIR, logger)
        logger.info(f"PDB 다운로드: 성공={n_dl}, 스킵={n_sk}, 실패={n_fl}")
    else:
        logger.info("PDB 다운로드 스킵 (--skip-download)")

    # AF2 바이너리 자동 탐지 (experimental 모드에서는 선택적이므로 실패해도 경고만)
    af2_bin = None
    if not args.use_experimental:
        af2_bin = find_colabfold_bin()
        if af2_bin:
            logger.info(f"[AF2] colabfold_batch: {af2_bin}")
        else:
            logger.warning("[AF2] colabfold_batch 바이너리를 찾지 못함 — AF2 단계가 실패합니다")

    # 3. 파이프라인 실행
    if not args.skip_run:
        try:
            run_all_benchmarks(
                records, PDB_DIR, RESULTS_DIR, args.mode, PROGRESS_FILE, logger,
                preprocess_dir=PREPROCESS_DIR,
                af2_input_dir=AF2_INPUT_DIR,
                af2_result_dir=AF2_RESULT_DIR,
                af2_bin=af2_bin,
                use_experimental=args.use_experimental,
            )
        except KeyboardInterrupt:
            logger.warning("[main] 사용자 중단 — 평가 없이 종료")
            return 130
    else:
        logger.info("파이프라인 실행 스킵 (--skip-run)")

    # 4. 평가 + 리포트 (항상 실행)
    rc = evaluate_and_report(
        records, RESULTS_DIR, REPORT_DIR, FIGURES_DIR, logger,
        pipeline_label=pipeline_label,
    )

    logger.info("=== UPDD Benchmark 완료 ===")
    return rc


if __name__ == "__main__":
    sys.exit(main())
