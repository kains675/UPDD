#!/usr/bin/env python
"""
UPDD.py — Universal Peptide Drug Discovery Pipeline
RFdiffusion (RFpeptides) → ProteinMPNN → AF2 → ncAA 치환 → MD → QM/MM → MM-GBSA → 랭킹

[Architecture Zenith Update]
- 단일 진실 공급원(ncaa_registry.py) 연동을 통한 ncAA 메타데이터 중앙 집중화
- [CRITICAL FIX] _run_parameterization 호출 시 ncaa_def 객체가 아닌 ncaa_def.code 문자열을 넘기도록 오타 수정
"""

import os
import sys
import time
import json
import subprocess
import glob
import shutil
import csv
import logging
import argparse
import concurrent.futures
import builtins as _builtins
import logging as _logging
import datetime as _datetime
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any

# ==========================================
# 1. 환경 및 경로 설정 (환경 변수 기반)
# ==========================================
# [Refactor] 하드코딩된 절대 경로를 환경 변수로 대체하여 이식성 확보.
# UPDD_HOME / UPDD_DIR / UPDD_TARGET_DIR / UPDD_PROJECT_DIR 환경 변수가
# 설정되어 있으면 우선 사용하고, 부재 시 기존 기본값으로 폴백한다.
HOME_DIR    = os.environ.get("UPDD_HOME", os.path.expanduser("~"))
UPDD_DIR    = os.environ.get("UPDD_DIR", os.path.join(HOME_DIR, "UPDD_proj"))
TARGET_DIR  = os.environ.get("UPDD_TARGET_DIR", os.path.join(UPDD_DIR, "target"))
UTILS_DIR   = os.environ.get("UPDD_UTILS_DIR", os.path.join(UPDD_DIR, "utils"))
LOG_DIR     = os.environ.get("UPDD_LOG_DIR", os.path.join(UPDD_DIR, "log"))
PROJECT_DIR = os.environ.get("UPDD_PROJECT_DIR", os.path.join(HOME_DIR, "ai_projects"))

# ==========================================
# Default tuning-knob environment
# ==========================================
# UPDD 가 내부 모듈 (utils/run_qmmm.py 등) 이나 subprocess (MD, Step 6.6 G2,
# Step 8 G3 antechamber) 에 전달하는 tuning 상수들은 env 변수로 gate 되어
# 있다 (SciVal-approved bundle + R-17). 개별 export 를 잊으면 조용히 기본
# 동작이 달라지는 버그 지뢰가 되므로, **os.environ.setdefault** 로 안전한
# default 를 주입한다. 사용자가 쉘에서 명시적으로 export 한 값은 그대로
# 우선한다 (setdefault 가 덮어쓰지 않음).
_UPDD_DEFAULT_ENV = {
    # --- R-17 Stage 4 gates ---
    "UPDD_ALLOW_ANTECHAMBER": "1",      # G3: cofactor FF param 자동 regen 허용
    # --- SCF bundle (SciVal B-series, wB97X-D direct-SCF) ---
    "UPDD_VRAM_POOL_FRACTION": "0.75",  # gpu4pyscf cupy mempool ceiling
    "UPDD_MIN_AO_BLKSIZE":     "64",    # B8: 16 GB VRAM boundary
    "UPDD_SCF_DIRECT_TOL":     "1e-12", # B5: direct-SCF tolerance
    "UPDD_SCF_CONV_TOL":       "1e-7",  # B3: ranking-safe convergence
    "UPDD_SCF_INIT_GUESS":     "minao", # B6: closed-shell safe init
}
for _k, _v in _UPDD_DEFAULT_ENV.items():
    os.environ.setdefault(_k, _v)

# The interim cofactor-injection diagnostic flag must be an **explicit**
# opt-in per run — never inherit from a prior shell that tested it.
# Clearing it here prevents stale leakage into production pipeline runs.
os.environ.pop("UPDD_COFACTOR_INJECT_DIAGNOSTIC", None)

# Registry Import (단일 계약)
sys.path.append(UTILS_DIR)
try:
    from ncaa_registry import resolve_ncaa_definition, NCAA_REGISTRY_DATA
except ImportError:
    sys.exit("[!] utils/ncaa_registry.py 파일을 찾을 수 없습니다. 경로를 확인하세요.")

# [v0.5] Per-design MD status SSOT
try:
    from utils_common import (load_md_status, write_md_status,  # type: ignore
                              classify_md_outcome, get_design_tier,
                              EXPLOSION_MARKER_PATTERNS, PARTIAL_PLACEHOLDERS)
    _V05_MD_STATUS_AVAILABLE = True
except ImportError:
    _V05_MD_STATUS_AVAILABLE = False

# [v0.4] Dashboard / State / Config — UI 탈바꿈 (과학 로직 불변)
# 세 유틸리티는 utils/ 에 있으며 UTILS_DIR 이 이미 sys.path 에 추가되어 있다.
# ImportError 발생 시 대시보드를 비활성화하고 v0.3.3 레거시 모드로 자동 폴백한다.
try:
    from updd_state import PipelineState, check_outputs_exist  # type: ignore
    from updd_dashboard import run_dashboard as _run_dashboard_standalone  # type: ignore
    from updd_dashboard import run_dashboard_in_thread  # type: ignore
    from updd_config import get_workers  # type: ignore
    _V04_DASHBOARD_AVAILABLE = True
except ImportError as _v04_imp_err:
    _V04_DASHBOARD_AVAILABLE = False
    _V04_IMPORT_ERROR = str(_v04_imp_err)
    PipelineState = None  # type: ignore
    run_dashboard_in_thread = None  # type: ignore
    def check_outputs_exist(step, inputs=None, outputs_dir=""):  # type: ignore
        return False
    def get_workers(step, inputs=None):  # type: ignore
        """v0.4 모듈 불가 시 fallback — parallel_workers 또는 1 반환."""
        if inputs and "parallel_workers" in inputs:
            try:
                return max(1, int(inputs["parallel_workers"]))
            except (ValueError, TypeError):
                return 1
        return 1

# [v0.4.1] 모듈 글로벌 PipelineState 참조 — main() 에서 설정되며,
# run_rfdiffusion / execute_qmmm 등 장시간 실행 함수가 진행도를 업데이트할 때 사용.
# main() 밖에서는 None 이 기본이므로 UPDD 를 라이브러리로 import 해도 영향 없음.
_pipeline_state: Optional["PipelineState"] = None

# [v42 FIX] run_cyclic_fold.py 는 colabdesign (JAX) 의존성을 가지므로 직접
# import 하면 md_simulation conda 환경에서 ImportError 가 발생한다.
# subprocess 호출로 변경하여 colabfold/pixi 환경에서 실행한다.
# 직접 import 블록은 제거하고 SCRIPT_PATHS["run_cyclic_fold"] 로 대체한다.
HAS_CYCLIC_FOLD = os.path.exists(os.path.join(UTILS_DIR, "run_cyclic_fold.py"))

# ==========================================
# 2. 모듈 전역 상태 컨테이너 (find_colabfold_bin 호출 이전에 정의)
# ==========================================
@dataclass
class _UPDDConfig:
    """모듈 전역 가변 상태를 캡슐화한 Config 컨테이너.

    파이프라인 실행 중 변경되는 상태(현재 logger, 현재 로그 파일 경로,
    ColabFold 바이너리 경로 등)를 단일 객체로 응집하여 `global` 사용을
    제거하고 테스트 격리성과 유지보수성을 확보한다.
    """
    logger: Optional[_logging.Logger] = None
    current_log_file: Optional[str] = None
    af2_bin: str = ""

_config = _UPDDConfig()

# [Backwards-Compat Wrappers]
# 기존 코드가 모듈 전역 변수(_updd_logger, _current_log_file, AF2_BIN)에 직접
# 접근할 가능성에 대비하여 동일 이름의 모듈 속성을 유지한다. 내부 함수는
# _config.* 를 권장 진입점으로 사용한다.
_updd_logger = None
_current_log_file = None


def prefer_fixed_script(path: Optional[str]) -> Optional[str]:
    if not path:
        return path
    root, ext = os.path.splitext(path)
    fixed_path = root + "_fixed" + ext
    return fixed_path if os.path.exists(fixed_path) else path

def find_colabfold_bin() -> Optional[str]:
    candidates: List[str] = []
    w = shutil.which("colabfold_batch")
    if w:
        candidates.append(w)
    candidates += glob.glob(os.path.join(PROJECT_DIR, "localcolabfold", ".pixi", "envs", "*", "bin", "colabfold_batch"))
    candidates += glob.glob(os.path.join(PROJECT_DIR, ".venv*", "bin", "colabfold_batch"))
    candidates += glob.glob(os.path.join(PROJECT_DIR, "venv*", "bin", "colabfold_batch"))
    candidates += glob.glob(os.path.expanduser("~/.pixi/envs/*/bin/colabfold_batch"))
    for base in [UPDD_DIR, HOME_DIR]:
        candidates += glob.glob(os.path.join(base, ".pixi", "envs", "*", "bin", "colabfold_batch"))
    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    return None

# [Refactor] AF2_BIN은 _config.af2_bin으로 위임. 모듈 전역 동기화를 위해
# 동일 이름의 변수를 유지한다(하위 호환). main() 진입 시 _config가 갱신된다.
_config.af2_bin = find_colabfold_bin() or ""
AF2_BIN = _config.af2_bin

RF_MODELS_DIR = os.path.join(PROJECT_DIR, "RFdiffusion/models")
RF_MODEL = {
    "linear":     os.path.join(RF_MODELS_DIR, "Complex_base_ckpt.pt"),
    "cyclic_htc": os.path.join(RF_MODELS_DIR, "Complex_base_ckpt.pt"),
    "cyclic_ss":  os.path.join(RF_MODELS_DIR, "Complex_base_ckpt.pt"),
    "cyclic_nm":  os.path.join(RF_MODELS_DIR, "Complex_base_ckpt.pt"),
    "bicyclic":   os.path.join(RF_MODELS_DIR, "Complex_base_ckpt.pt"),
}

ENV_NAMES = {
    "md":   "md_simulation",
    "rf":   "rfdiffusion",
    "mpnn": "proteinmpnn",
    "qmmm": "qmmm",
    "ncaa": "ncaa",
}

# [v0.3 #2] 자식 프로세스 stdout 즉시 flush — PYTHONUNBUFFERED=1.
# subprocess.Popen 이 PIPE 를 사용하면 자식 Python 의 stdout 이 블록 버퍼링
# (~8KB) 으로 전환되어 DFT cycle 로그가 실시간으로 보이지 않는다. 이 환경
# 변수를 자식에 주입하면 print 호출이 줄 단위로 즉시 flush 되어 tee/리다이렉션
# 없이도 실시간 로그가 보장된다. 모듈 레벨에서 한 번 생성하여 모든 호출부에서
# 재사용한다 (CLAUDE.md v0.3 #2).
_CHILD_ENV_UNBUFFERED = os.environ.copy()
_CHILD_ENV_UNBUFFERED["PYTHONUNBUFFERED"] = "1"

SCRIPT_PATHS = {
    "rf_inference": os.path.join(PROJECT_DIR, "RFdiffusion/scripts/run_inference.py"),
    "mpnn_parse":   os.path.join(PROJECT_DIR, "ProteinMPNN/helper_scripts/parse_multiple_chains.py"),
    "mpnn_run":     os.path.join(PROJECT_DIR, "ProteinMPNN/protein_mpnn_run.py"),
    "grafting":     os.path.join(UTILS_DIR, "graft_ligand_generic.py"),
    "process_af2":  os.path.join(UTILS_DIR, "prepare_af2_input.py"),
    "preprocess_target": os.path.join(UTILS_DIR, "preprocess_target.py"),
    "rank_af2":     os.path.join(UTILS_DIR, "rank_results.py"),
    "ncaa_mutate":  os.path.join(UTILS_DIR, "ncaa_mutate.py"),
    "parameterize": os.path.join(UTILS_DIR, "parameterize_ncaa.py"),
    "run_md":       os.path.join(UTILS_DIR, "run_restrained_md.py"),
    "extract_snap": os.path.join(UTILS_DIR, "extract_snapshots.py"),
    "run_qmmm":     os.path.join(UTILS_DIR, "run_qmmm.py"),
    "run_mmgbsa":   os.path.join(UTILS_DIR, "run_mmgbsa.py"),
    "rank_qmmm":    os.path.join(UTILS_DIR, "rank_results_qmmm.py"),
    "run_cyclic_fold": os.path.join(UTILS_DIR, "run_cyclic_fold.py"),
}

TOPOLOGY_TYPES = {
    1: {"name": "Linear",           "abbr": "linear",     "rf_T": 200, "cyclic": False, "canonical": "linear",     "prediction": "colabfold",   "note": "일반 선형 펩타이드 (기본값)"},
    2: {"name": "Head-to-Tail Cyclic","abbr": "cyclic_htc","rf_T": 50,  "cyclic": True,  "canonical": "cyclic_htc", "prediction": "afcycdesign", "note": "N-C 말단 아미드 결합 고리화"},
    3: {"name": "Disulfide Cyclic", "abbr": "cyclic_ss",  "rf_T": 50,  "cyclic": True,  "canonical": "cyclic_ss",  "prediction": "rfpeptides",  "note": "두 Cys S-S 결합 고리화"},
    4: {"name": "N-Methyl Cyclic",  "abbr": "cyclic_nm",  "rf_T": 50,  "cyclic": True,  "canonical": "cyclic_htc", "prediction": "afcycdesign", "note": "N-메틸 잔기 포함 고리화 (내부적으로 HTC 처리)"},
    5: {"name": "Bicyclic",         "abbr": "bicyclic",   "rf_T": 50,  "cyclic": True,  "canonical": "bicyclic",   "prediction": "rfpeptides",  "note": "thioether + disulfide 이중 고리"},
}

# ==========================================
# 3. 헬퍼 함수
# ==========================================
_original_input = _builtins.input

def setup_logger(target_name: str, append_from: Optional[str] = None) -> logging.Logger:
    global _updd_logger, _current_log_file
    os.makedirs(LOG_DIR, exist_ok=True)
    target_log_dir = os.path.join(LOG_DIR, target_name)
    os.makedirs(target_log_dir, exist_ok=True)
    timestamp = _datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(target_log_dir, f"updd_{timestamp}.log")
    if append_from and os.path.exists(append_from):
        with open(log_file, 'w', encoding='utf-8') as out_f, open(append_from, 'r', encoding='utf-8', errors='ignore') as in_f:
            out_f.write(in_f.read())
    logger = logging.getLogger('UPDD')
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    for h in logger.handlers[:]:
        try:
            h.close()
        except (OSError, ValueError) as _close_err:
            print(f"[LOG] 기존 핸들러 정리 실패(무시): {_close_err}")
        logger.removeHandler(h)
    fh = logging.FileHandler(log_file, encoding='utf-8')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
    logger.addHandler(fh)
    _config.logger = logger
    _config.current_log_file = log_file
    _updd_logger = logger
    _current_log_file = log_file
    print(f'[LOG] 로그 파일: {log_file}')
    return logger

def log(msg: str, level: str = 'info') -> None:
    if _config.logger:
        getattr(_config.logger, level, _config.logger.info)(msg)

def logged_input(prompt: str = '') -> str:
    if prompt:
        log(f'[PROMPT] {prompt}')
    ans = _original_input(prompt)
    log(f'[INPUT] {ans}')
    return ans

def print_step(msg: str) -> None:
    line = '=' * 60
    print('\n' + line + '\n ' + msg + '\n' + line + '\n')
    log(line)
    log(' ' + msg)
    log(line)

def check_error(returncode: int, msg: str) -> None:
    if returncode != 0:
        print(f"\n[!] {msg} 실패!\n Error Code: {returncode}\n (15초 후 자동 종료)")
        time.sleep(15)
        sys.exit(1)

def run_conda_command(script_path: str, args: List[str], description: str, env_key: str) -> None:
    print(f"\n>>> {description} 실행 중... (Env: {env_key})")
    if not os.path.exists(script_path):
        sys.exit(f"[!] 스크립트 없음: {script_path}")
    env_name = ENV_NAMES.get(env_key)
    cmd = ["conda", "run", "-n", env_name, "--no-capture-output", "python", script_path] + args
    try:
        # [v0.3 #2] PYTHONUNBUFFERED=1 주입 — 자식 Python 의 stdout 즉시 flush.
        # 콘솔/리다이렉션/tee 어느 경로든 실시간 로그가 보장된다.
        subprocess.run(cmd, check=True, env=_CHILD_ENV_UNBUFFERED)
        print(f" {description} 완료.")
    except subprocess.CalledProcessError as e:
        check_error(e.returncode, description)
    except KeyboardInterrupt:
        sys.exit("\n[중단됨]")


def _find_colabfold_python() -> Optional[str]:
    """localcolabfold pixi 환경의 Python 바이너리를 탐색한다.

    ColabDesign (AfCycDesign) 은 JAX 가 설치된 pixi 환경에서만 실행 가능하다.
    conda 환경과 별개이므로 직접 바이너리 경로를 탐색한다.
    """
    for base in (PROJECT_DIR, HOME_DIR):
        hits = glob.glob(os.path.join(base, "**", "localcolabfold",
                                      ".pixi", "envs", "*", "bin", "python3"), recursive=False)
        # non-recursive 먼저, 그래도 없으면 직접 경로
        if not hits:
            hits = glob.glob(os.path.join(base, "localcolabfold",
                                          ".pixi", "envs", "*", "bin", "python3"))
        for h in hits:
            if os.path.isfile(h) and os.access(h, os.X_OK):
                return h
    return None


def run_pixi_command(script_path: str, args: List[str], description: str) -> None:
    """localcolabfold pixi 환경의 Python 으로 스크립트를 실행한다.

    ColabDesign (AfCycDesign) 전용. conda 환경에는 colabdesign 이 없으므로
    pixi 환경의 Python 바이너리를 직접 호출한다.
    """
    pybin = _find_colabfold_python()
    if pybin is None:
        sys.exit(f"[!] localcolabfold pixi Python 을 찾을 수 없습니다. "
                 f"colabdesign 설치 환경을 확인하세요.")
    print(f"\n>>> {description} 실행 중... (pixi: {os.path.dirname(os.path.dirname(pybin))})")
    if not os.path.exists(script_path):
        sys.exit(f"[!] 스크립트 없음: {script_path}")
    cmd = [pybin, script_path] + args
    try:
        subprocess.run(cmd, check=True)
        print(f" {description} 완료.")
    except subprocess.CalledProcessError as e:
        check_error(e.returncode, description)
    except KeyboardInterrupt:
        sys.exit("\n[중단됨]")


# [v3 1-1] arrow_menu / select_topology / select_ncaa / get_target_preprocess_options /
# gather_all_inputs 는 utils/updd_cli.py 로 분리되었다. 본 모듈은 UI 책임을 갖지 않는다.
from updd_cli import (  # noqa: E402
    arrow_menu,
    select_topology as _select_topology_ui,
    select_ncaa as _select_ncaa_ui,
    get_target_preprocess_options as _get_target_preprocess_options_ui,
    gather_all_inputs as _gather_all_inputs_ui,
)

def get_continuous_contigs(pdb_path: str, chain_id: str) -> str:
    residues: List[int] = []
    try:
        with open(pdb_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("ATOM") and len(line) > 26 and line[21] == chain_id:
                    try:
                        res_num = int(line[22:26].strip())
                        if not residues or residues[-1] != res_num:
                            residues.append(res_num)
                    except ValueError as _parse_err:
                        log(f"[get_continuous_contigs] residue 번호 파싱 실패: {_parse_err}", level="warning")
    except (IOError, OSError) as _io_err:
        log(f"[get_continuous_contigs] PDB 파일 접근 실패: {pdb_path} ({_io_err})", level="warning")
        return chain_id + "1-100"
    if not residues:
        return chain_id + "1-100"
    contigs: List[str] = []
    start, prev = residues[0], residues[0]
    for res in residues[1:]:
        if res != prev + 1:
            contigs.append(chain_id + str(start) + "-" + str(prev))
            start = res
        prev = res
    contigs.append(chain_id + str(start) + "-" + str(prev))
    return "/".join(contigs)

def check_rf_model(model_path: str) -> bool:
    if os.path.exists(model_path):
        return True
    print(f"[!] 모델 가중치 없음: {model_path}")
    return False

def infer_binder_chains_from_entry(entry, target_chain):
    all_chains = sorted([k.replace("seq_chain_", "") for k in entry.keys() if k.startswith("seq_chain_")])
    target_c = (target_chain or "A").upper()
    binder_chains = [c for c in all_chains if c != target_c]
    return binder_chains, all_chains

def make_tied_positions_jsonl(jsonl_path: str, output_path: str, target_chain: str = "A") -> None:
    tied_dict: Dict[str, Any] = {}
    with open(jsonl_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            entry = json.loads(line)
            raw_name = entry.get("name", "unknown")
            stem_name = os.path.splitext(os.path.basename(raw_name))[0]
            binder_chains, _ = infer_binder_chains_from_entry(entry, target_chain)
            tied_val = [{chain_id: [1, len(entry.get(f"seq_chain_{chain_id}", ""))]} for chain_id in binder_chains if entry.get(f"seq_chain_{chain_id}", "")]
            tied_dict[stem_name] = tied_val
            if raw_name != stem_name:
                tied_dict[raw_name] = tied_val
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(json.dumps(tied_dict) + "\n")

def load_json_manifest(path: str) -> Optional[Dict[str, Any]]:
    # [v35 AUDIT] Principle 5: bare 'except: return None' eliminated — now catches only JSONDecodeError.
    # [v35 AUDIT] Principle 4/v27: encoding="utf-8" enforced to prevent UnicodeDecodeError on manifests containing 한글 metadata.
    if not os.path.exists(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except json.JSONDecodeError as _err:
        print(f"  [!] Manifest 파싱 실패 (손상된 JSON): {path} — {_err}")
        return None

# ==========================================
# 모듈 실행 함수들
# ==========================================
# [v3 1-1] UI 함수는 updd_cli.py 에 위임. 본 모듈에는 UPDD.py 고유의 도메인
# 데이터(TOPOLOGY_TYPES, NCAA_REGISTRY_DATA) 를 주입한 얇은 wrapper 만 남긴다.
def get_target_preprocess_options(target_pdb):
    return _get_target_preprocess_options_ui(target_pdb)

def run_target_preprocessing(target_pdb: str, keep_mode: str, chain_in: str, keep_res: str,
                              target_id: str = "") -> str:
    """Step 3 orchestration wrapper for ``utils/preprocess_target.py``.

    [R-17 G1] When ``target_id`` is provided and ``target_cards/<id>.json``
    exists, the subprocess is launched with ``--target_id`` so the cofactor
    preservation gate auto-injects declared required resnames into
    ``--keep_hetatms`` and validates their survival after HETATM filtering.
    For cofactor-absent targets the gate is a transparent no-op (returns
    before raising).
    """
    output_clean = os.path.splitext(target_pdb)[0] + "_clean.pdb"
    args = ["--input", target_pdb, "--output", output_clean, "--chains", chain_in, "--hetatm_mode", keep_mode]
    if keep_res:
        args += ["--keep_hetatms", keep_res]
    # [R-17 G1] Pass target_id when the card exists so preprocess_target.py
    # can load target_card.cofactor_residues and enforce keep-list + presence
    # check. Missing card → legacy behavior (the subprocess emits a [WARN]).
    if target_id:
        _card_path = os.path.join(UPDD_DIR, "target_cards", f"{target_id}.json")
        if os.path.exists(_card_path):
            args += ["--target_id", target_id]
            print(f"  [R-17 G1] target_card '{target_id}' 감지 — cofactor keep-list 자동 주입")
    run_conda_command(SCRIPT_PATHS["preprocess_target"], args, "Target Preprocessing", env_key="md")
    return output_clean

def select_topology():
    return _select_topology_ui(TOPOLOGY_TYPES)

def select_ncaa(binder_len_str: str, binder_chain: str = "B", binder_chains: Optional[List[str]] = None):
    return _select_ncaa_ui(binder_len_str, NCAA_REGISTRY_DATA,
                           binder_chain=binder_chain, binder_chains=binder_chains)

def gather_all_inputs(base_name: str, target_pdb: str) -> Dict[str, Any]:
    return _gather_all_inputs_ui(base_name, target_pdb, TOPOLOGY_TYPES, NCAA_REGISTRY_DATA)


# ==========================================
# 파이프라인 코어 실행 함수 (RFdiffusion ~ QMMM)
# ==========================================
def run_ligand_grafting(work_pdb, ref_pdb, target_chain="", ref_chain=""):
    # [v35 AUDIT] Principle 7: --target_chain / --ref_chain now threaded through
    # to graft_ligand_generic.py so CA superposition aligns the correct chain
    # rather than hardcoded 'chainid 0'. Empty strings preserve legacy behavior.
    # [v3 1-4] TECH DEBT 해소: main() 의 preprocessing 직후 단계에서 본 함수가
    # 호출되도록 wiring 완료. inputs["grafting_needed"] 와 inputs["ref_pdb"] 가
    # 모두 truthy 인 경우 work_pdb 가 grafted 결과로 교체된다.
    grafted_pdb = os.path.join(TARGET_DIR, "grafted_target.pdb")
    graft_args = ["--target", work_pdb, "--ref", ref_pdb, "--output", grafted_pdb]
    if target_chain:
        graft_args += ["--target_chain", target_chain]
    if ref_chain:
        graft_args += ["--ref_chain", ref_chain]
    run_conda_command(SCRIPT_PATHS["grafting"], graft_args, "Ligand Grafting", env_key="md")
    return grafted_pdb

def run_rfdiffusion(work_pdb: str, topology: Dict[str, Any], hotspots: str, final_contig: str,
                    design_dir: str, total_designs: int, parallel_workers: int, target_chain: str) -> None:
    print_step("Step 4: RFdiffusion")
    model_path = RF_MODEL.get(topology["abbr"], RF_MODEL["linear"])
    check_rf_model(model_path)

    def run_rfd_worker(worker_id: int, worker_designs: int) -> None:
        w_prefix = os.path.join(design_dir, f"design_w{worker_id}")
        existing_count = len(glob.glob(f"{w_prefix}*.pdb"))
        if existing_count >= worker_designs:
            return
        remaining_designs = worker_designs - existing_count
        if existing_count > 0:
            w_prefix = f"{w_prefix}_resume{existing_count}"

        rf_log_dir = os.path.join(os.path.dirname(design_dir), "rf_logs", f"worker_{worker_id}")
        w_rf_args = [
            "inference.input_pdb=" + work_pdb, "inference.output_prefix=" + w_prefix,
            "inference.num_designs=" + str(remaining_designs), "contigmap.contigs=[" + final_contig + "]",
            "diffuser.T=" + str(topology["rf_T"]), "inference.ckpt_override_path=" + model_path,
            f"hydra.run.dir={rf_log_dir}",
        ]
        if topology["cyclic"]:
            w_rf_args.append("+inference.cyclic_peptide=True")
        if topology["abbr"] == "cyclic_ss":
            w_rf_args.append("+inference.add_cyclic_disulfide=True")
        if hotspots:
            w_rf_args.append("ppi.hotspot_res=[" + ",".join([target_chain + h.strip() for h in hotspots.split(",")]) + "]")
        run_conda_command(SCRIPT_PATHS["rf_inference"], w_rf_args, f"RFdiffusion (Worker {worker_id})", env_key="rf")

    if parallel_workers > 1:
        # [v0.4.1 B] as_completed 로 변환하여 worker 완료 개수를 state 에 반영.
        # 예외 동작은 legacy (wait) 와 동일하게 — worker 실패는 warn 만 하고 계속 진행.
        _ps = _pipeline_state
        with concurrent.futures.ThreadPoolExecutor(max_workers=parallel_workers) as executor:
            futures = []
            for i in range(parallel_workers):
                w_designs = total_designs // parallel_workers + (1 if i < total_designs % parallel_workers else 0)
                if w_designs > 0:
                    futures.append(executor.submit(run_rfd_worker, i + 1, w_designs))
            n_workers_actual = len(futures)
            if _ps is not None and n_workers_actual > 0:
                _ps.update_step("rfdiff", completed=0,
                                detail=f"0/{n_workers_actual} workers done")
            _done_workers = 0
            for fut in concurrent.futures.as_completed(futures):
                try:
                    fut.result()
                except Exception as _we:  # legacy: worker 실패 무시
                    print(f"  [!] RFdiffusion worker 실패: {_we}")
                _done_workers += 1
                if _ps is not None:
                    _ps.update_step("rfdiff", completed=_done_workers,
                                    detail=f"{_done_workers}/{n_workers_actual} workers done")
    else:
        run_rfd_worker(1, total_designs)
        _ps = _pipeline_state
        if _ps is not None:
            _ps.update_step("rfdiff", completed=1, detail="1/1 workers done")

def run_proteinmpnn(design_dir: str, target_chain: str, topology: Dict[str, Any], valid_only: bool = False):
    print_step("Step 5: ProteinMPNN — Sequence Design")
    mpnn_out_dir = design_dir.replace("designs", "mpnn_results")
    os.makedirs(mpnn_out_dir, exist_ok=True)
    jsonl_path = os.path.join(mpnn_out_dir, "parsed_pdbs.jsonl")

    if not valid_only:
        run_conda_command(SCRIPT_PATHS["mpnn_parse"], ["--input_path=" + design_dir, "--output_path=" + jsonl_path], "MPNN Parsing", env_key="mpnn")

    chain_id_jsonl = os.path.join(mpnn_out_dir, "chain_id.jsonl")
    if not valid_only:
        chain_id_dict: Dict[str, Any] = {}
        with open(jsonl_path, "r", encoding="utf-8") as f:
            for line in f:
                if not line.strip():
                    continue
                entry = json.loads(line.strip())
                target_c = target_chain.upper()
                binder_chains, _ = infer_binder_chains_from_entry(entry, target_c)
                chain_id_dict[entry["name"]] = [binder_chains, [target_c]]
        with open(chain_id_jsonl, "w", encoding="utf-8") as f:
            f.write(json.dumps(chain_id_dict) + "\n")

    seqs_dir = os.path.join(mpnn_out_dir, "seqs")
    os.makedirs(seqs_dir, exist_ok=True)
    if len(glob.glob(os.path.join(seqs_dir, "*.fa"))) > 0:
        return mpnn_out_dir, jsonl_path

    mpnn_args = [
        "--jsonl_path", jsonl_path, "--chain_id_jsonl", chain_id_jsonl,
        "--out_folder", mpnn_out_dir, "--num_seq_per_target", "2",
        "--batch_size", "1", "--sampling_temp", "0.1"
    ]
    if topology["cyclic"]:
        tied_jsonl = os.path.join(mpnn_out_dir, "tied_positions.jsonl")
        make_tied_positions_jsonl(jsonl_path, tied_jsonl, target_chain=target_chain)
        mpnn_args += ["--tied_positions_jsonl", tied_jsonl]

    if not valid_only:
        run_conda_command(SCRIPT_PATHS["mpnn_run"], mpnn_args, "ProteinMPNN", env_key="mpnn")
    return mpnn_out_dir, jsonl_path

# ==========================================
# [v3 1-2] run_alphafold 4분할 — 헬퍼 함수 (R-1: run_alphafold 정의보다 위에 배치)
# ==========================================
def _fetch_target_msa(target_csv: str, target_msa_dir: str, target_a3m_path: str,
                      rep_target_seq: str, af2_bin: str) -> None:
    """[v3 1-2] 타겟 단일 서열에 대한 MSA 를 ColabFold 로 1회 요청한다.

    PENDING 상태가 30분(MAX_PENDING_SEC) 이상 지속되면 warp-cli 로 IP rollover 를
    시도한 뒤 재시도한다. 정상 종료되면 함수 종료, 비정상 종료 시 check_error 가
    파이프라인 전체를 즉시 중단시킨다.

    Args:
        target_csv:      ColabFold 에 전달할 단일 서열 CSV 파일 경로 (입력)
        target_msa_dir:  ColabFold 출력 디렉토리
        target_a3m_path: 생성될 target_only.a3m 경로 (사전 존재 시 본 함수는 호출되지 않음)
        rep_target_seq:  대표 타겟 서열 (CSV 파일에 기록할 sequence)
        af2_bin:         colabfold_batch 바이너리 절대 경로
    """
    with open(target_csv, "w", encoding="utf-8") as out_csv:
        out_csv.write("id,sequence\n")
        out_csv.write(f"target_only,{rep_target_seq}\n")
    print("\n  [🛡️] 타겟 MSA 데이터를 1회 요청합니다.")

    MAX_PENDING_SEC = 1800
    while True:
        needs_reroll = False
        pending_start = None
        process = subprocess.Popen(
            [af2_bin, target_csv, target_msa_dir, "--msa-only"],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
        )
        buffer = ""
        while True:
            char = process.stdout.read(1)
            if not char and process.poll() is not None:
                break
            if char:
                sys.stdout.write(char)
                sys.stdout.flush()
                if char in ('\r', '\n'):
                    line = buffer.strip()
                    buffer = ""
                    if not line:
                        continue
                    if "PENDING" in line or "Sleeping for" in line:
                        if pending_start is None:
                            pending_start = time.time()
                        elif time.time() - pending_start > MAX_PENDING_SEC:
                            print(f"\n\n[!] 🚨 PENDING 지연으로 IP 교체를 시도합니다.")
                            process.terminate()
                            needs_reroll = True
                            break
                    else:
                        pending_start = None
                else:
                    buffer += char
        process.wait()
        if needs_reroll:
            subprocess.run(["warp-cli", "disconnect"], capture_output=True)
            time.sleep(3)
            subprocess.run(["warp-cli", "connect"], capture_output=True)
            time.sleep(10)
            continue
        if process.returncode == 0:
            return
        check_error(process.returncode, "Target MSA Generation")
        return


def _build_a3m_inputs(mpnn_out_dir: str, af2_out_dir: str, csv_path: str,
                     target_seqs: Dict[str, str], target_a3m_path: str,
                     use_local_msa: bool) -> int:
    """[v3 1-2] MPNN 시퀀스와 타겟 MSA 를 결합하여 ColabFold 입력 (CSV + A3M) 을 생성한다.

    Args:
        mpnn_out_dir:    MPNN 결과 디렉토리 (seqs/ 하위 *.fa 파일 사용)
        af2_out_dir:     AF2 출력 디렉토리 (a3m_inputs/ 및 *.a3m 생성 위치)
        csv_path:        ColabFold 가 읽을 input CSV 경로 (입력+출력)
        target_seqs:     {design_name: target_chain_sequence} 사전
        target_a3m_path: _fetch_target_msa 가 생성한 a3m 파일 경로
        use_local_msa:   True 이면 design 별 a3m 파일을 함께 생성 (local MSA mode)

    Returns:
        int: 생성된 고유 서열 조합 개수
    """
    target_msa_headers: List[str] = []
    target_msa_seqs: List[str] = []
    with open(target_a3m_path, "r", encoding="utf-8") as f:
        curr_h, curr_s = "", []
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if curr_h:
                    target_msa_headers.append(curr_h)
                    target_msa_seqs.append("".join(curr_s))
                curr_h = line
                curr_s = []
            else:
                curr_s.append(line)
        if curr_h:
            target_msa_headers.append(curr_h)
            target_msa_seqs.append("".join(curr_s))

    seqs_dir = os.path.join(mpnn_out_dir, "seqs")
    seen_sequences = set()
    unique_count = 0
    a3m_dir = os.path.join(af2_out_dir, "a3m_inputs")
    os.makedirs(a3m_dir, exist_ok=True)

    with open(csv_path, "w", encoding="utf-8") as out_csv:
        out_csv.write("id,sequence\n")
        for fa_file in glob.glob(os.path.join(seqs_dir, "*.fa")):
            base_name = os.path.basename(fa_file).replace(".fa", "")
            t_seq_clean = target_seqs.get(base_name, "").replace("-", "").strip().upper()
            if not t_seq_clean:
                continue

            with open(fa_file, "r", encoding="utf-8") as fa:
                lines = fa.read().splitlines()

            sample_idx = 1
            for i in range(len(lines)):
                if lines[i].startswith(">T="):
                    binder_seq = lines[i + 1].strip().replace("-", "").strip().upper()
                    if binder_seq in seen_sequences:
                        continue
                    seen_sequences.add(binder_seq)
                    unique_count += 1

                    job_id = f"{base_name}_s{sample_idx}"
                    out_csv.write(f"{job_id},{t_seq_clean}:{binder_seq}\n")

                    if use_local_msa:
                        job_a3m = os.path.join(af2_out_dir, f"{job_id}.a3m")
                        a3m_lines = [f">101\t102", f"{t_seq_clean}{binder_seq}"]
                        for h, s in zip(target_msa_headers[1:], target_msa_seqs[1:]):
                            a3m_lines.append(f">{h.lstrip('>').strip()}")
                            a3m_lines.append(f"{s}{'-' * len(binder_seq)}")
                        with open(job_a3m, "w", encoding="utf-8") as f:
                            f.write("\n".join(a3m_lines) + "\n")
                    sample_idx += 1

    return unique_count


def _run_colabfold(csv_path: str, af2_out_dir: str, af2_bin: str) -> None:
    """[v3 1-2] ColabFold 를 실행한다 (mmseqs2_uniref_env, num-recycle=3, unpaired_paired)."""
    print("\n ColabFold 예측 시작...")
    try:
        subprocess.run(
            [af2_bin, csv_path, af2_out_dir,
             "--num-recycle", "3",
             "--pair-mode", "unpaired_paired",
             "--msa-mode", "mmseqs2_uniref_env"],
            check=True,
        )
        print(" AlphaFold2 예측 완료.")
    except subprocess.CalledProcessError as e:
        check_error(e.returncode, "ColabFold 예측")


def run_alphafold(design_dir: str, mpnn_out_dir: str, jsonl_path: str, target_chain: str,
                  msa_mode: str, use_local_msa: bool, valid_only: bool = False,
                  topology_mode: Optional[Dict[str, Any]] = None,
                  work_pdb: Optional[str] = None,
                  binder_chain: str = "B") -> str:
    """[v3 1-2] _fetch_target_msa → _build_a3m_inputs → (_run_colabfold | run_cyclic_fold)
    의 얇은 오케스트레이터.

    [v5 해법 C] topology_mode 가 cyclic 토폴로지를 가리키는 경우, 전용 cyclic
    fold 도구 (AfCycDesign / RFpeptides) 또는 ColabFold 폴백 + N-C 거리 검증
    경로로 분기한다. CSV 빌드 / ranking 단계는 두 경로에서 동일하므로 본
    함수의 외층에서 한 번만 수행된다 (중복 코드 제거 + 데이터 계약 일관성).
    """
    is_cyclic = bool(topology_mode and topology_mode.get("cyclic", False))
    if is_cyclic:
        print_step("Step 6: Stability Check — Cyclic Peptide Fold (AfCycDesign / RFpeptides)")
    else:
        print_step("Step 6: Stability Check — AlphaFold2 (ColabFold)")
    af2_out_dir = design_dir.replace("designs", "af2_results")
    os.makedirs(af2_out_dir, exist_ok=True)
    csv_path = os.path.join(af2_out_dir, "input.csv")

    if not valid_only and len(glob.glob(os.path.join(af2_out_dir, "*_scores.json"))) == 0:
        print(" [AF2] 타겟 서열과 바인더 서열을 결합하여 복합체 입력(CSV)을 생성합니다.")
        target_seqs: Dict[str, str] = {}
        with open(jsonl_path, "r", encoding="utf-8") as f:
            for line in f:
                if not line.strip():
                    continue
                entry = json.loads(line.strip())
                target_seqs[entry["name"]] = entry.get(f"seq_chain_{target_chain.upper()}", "")

        target_msa_dir = os.path.join(af2_out_dir, "target_msa")
        os.makedirs(target_msa_dir, exist_ok=True)
        target_csv = os.path.join(target_msa_dir, "target_only.csv")
        target_a3m_path = os.path.join(target_msa_dir, "target_only.a3m")

        rep_target_name = list(target_seqs.keys())[0]
        rep_target_seq = target_seqs[rep_target_name].replace("-", "").strip().upper()

        if not os.path.exists(target_a3m_path) or os.path.getsize(target_a3m_path) == 0:
            _fetch_target_msa(target_csv, target_msa_dir, target_a3m_path,
                              rep_target_seq, _config.af2_bin)

        unique_count = _build_a3m_inputs(mpnn_out_dir, af2_out_dir, csv_path,
                                         target_seqs, target_a3m_path, use_local_msa)
        print(f" [AF2] {unique_count}개의 고유 서열 조합 준비 완료.")

        if is_cyclic:
            # [v44] cyclic 분기: ColabFold + geometric gap closure.
            # run_cyclic_fold.py 를 pixi 환경에서 실행한다 (ColabFold CLI 사용).
            canonical_topo = topology_mode.get("canonical", topology_mode.get("abbr", "cyclic_htc"))
            if HAS_CYCLIC_FOLD:
                cyc_args = [
                    "--input_csv",     csv_path,
                    "--output_dir",    af2_out_dir,
                    "--topology_type", canonical_topo,
                    "--af2_bin",       _config.af2_bin or "",
                ]
                try:
                    run_pixi_command(
                        SCRIPT_PATHS["run_cyclic_fold"],
                        cyc_args,
                        "Cyclic 구조 예측 (ColabFold + gap closure)",
                    )
                except SystemExit:
                    print(" [Warning] run_cyclic_fold subprocess 실패 — ColabFold 폴백 (gap closure 없음)")
                    _run_colabfold(csv_path, af2_out_dir, _config.af2_bin)
            else:
                print(" [Warning] run_cyclic_fold.py 없음 — ColabFold 직접 호출 (gap closure 없음)")
                _run_colabfold(csv_path, af2_out_dir, _config.af2_bin)
        else:
            _run_colabfold(csv_path, af2_out_dir, _config.af2_bin)

    if not valid_only:
        run_conda_command(SCRIPT_PATHS["rank_af2"],
                          ["--inputdir", af2_out_dir,
                           "--outputcsv", os.path.join(af2_out_dir, "af2_ranking.csv")],
                          "AF2 랭킹", env_key="mpnn")
    return af2_out_dir

def _reinsert_cofactors_into_af2_outputs(af2_dir: str, crystal_pdb: str,
                                         target_id: str,
                                         target_chain: str = "A") -> dict:
    """[R-17 G2] Post-AF2 cofactor reinsertion dispatch per rank_001 design.

    AF2 / ColabFold predictions do not retain HETATM records. For
    cofactor-present targets (e.g. 6WGN with GNP / Mg²⁺), declared cofactors
    must be transplanted from the crystal PDB into each AF2 output before
    the pipeline continues into ncAA mutation (Step 7) and MD setup (Step 9),
    otherwise the downstream stages silently drop the cofactors (R-17 SSOT
    violation). The per-design rewrite is done in-place with a backup copy
    preserved under ``af2_dir/_pre_reinsert_backup/``.

    Behavior:
        - target_id empty OR target_cards/<id>.json missing → log + no-op.
        - target_card loaded but ``cofactor_residues`` empty / no required
          entries → log + no-op (cofactor-absent target).
        - Idempotent: if ``_pre_reinsert_backup/`` already has all rank_001
          basenames, skip (resume mode).
        - Any ``CofactorReinsertionError`` propagates → pipeline halts with
          the original exception message (low Cα count, RMSD too high,
          missing crystal HETATM, missing crystal PDB).

    Args:
        af2_dir: Output directory of the AF2 stage (rank_001 *.pdb already
            filtered by ``_filter_top_rank`` / ``_filter_by_iptm``).
        crystal_pdb: Path to the preprocessed crystal PDB (``work_pdb``) —
            still contains HETATM because G1 preserved declared cofactors.
        target_id: target_card stem (e.g. "6WGN"). Empty → no-op.
        target_chain: Target chain identifier for Kabsch alignment. Passed
            through to ``reinsert_cofactors`` as ``target_chain_override``.

    Returns:
        Dict with ``status``, ``target_id``, ``n_processed``, ``n_skipped``,
        ``details`` (list of per-design diagnostics). Status values:
        ``no_target_id``, ``no_card``, ``no_required_cofactors``,
        ``ok`` (all designs processed), ``skip_cached`` (already processed).

    Raises:
        CofactorReinsertionError: Bubbled from ``reinsert_cofactors`` when a
            per-design transplant is geometrically unsafe or inputs are
            inconsistent. Halts the pipeline — R-17 fail-fast.
    """
    # [R-17 G2] Early short-circuits for cofactor-absent path.
    if not target_id:
        print("  [R-17 G2] target_id 미지정 — cofactor 재삽입 건너뜀")
        return {"status": "no_target_id", "target_id": "", "n_processed": 0,
                "n_skipped": 0, "details": []}

    _card_path = os.path.join(UPDD_DIR, "target_cards", f"{target_id}.json")
    if not os.path.exists(_card_path):
        print(f"  [R-17 G2] target_card 없음 ({_card_path}) — 재삽입 건너뜀")
        return {"status": "no_card", "target_id": target_id, "n_processed": 0,
                "n_skipped": 0, "details": []}

    # Load target_card (auto-upgrade to v0.6.5 handled by utils_common).
    if UTILS_DIR not in sys.path:
        sys.path.append(UTILS_DIR)
    try:
        from utils_common import load_target_card  # type: ignore
        from reinsert_cofactors_post_af2 import reinsert_cofactors  # type: ignore
    except ImportError as _imp_err:
        print(f"  [R-17 G2][WARN] module import 실패 ({_imp_err}) — 재삽입 건너뜀")
        return {"status": "import_error", "target_id": target_id,
                "n_processed": 0, "n_skipped": 0, "details": []}

    _cards_dir = os.path.join(UPDD_DIR, "target_cards")
    try:
        target_card = load_target_card(target_id, cards_dir=_cards_dir)
    except (FileNotFoundError, ValueError) as _load_err:
        print(f"  [R-17 G2][WARN] target_card 로드 실패 ({_load_err}) — 재삽입 건너뜀")
        return {"status": "card_load_error", "target_id": target_id,
                "n_processed": 0, "n_skipped": 0, "details": []}

    required_entries = [e for e in (target_card.get("cofactor_residues") or [])
                        if isinstance(e, dict) and e.get("required")]
    if not required_entries:
        print(f"  [R-17 G2] '{target_id}' 선언된 required cofactor 없음 — no-op")
        return {"status": "no_required_cofactors", "target_id": target_id,
                "n_processed": 0, "n_skipped": 0, "details": []}

    # Crystal PDB sanity.
    if not crystal_pdb or not os.path.exists(crystal_pdb):
        print(f"  [R-17 G2][WARN] crystal PDB 없음 ({crystal_pdb}) — 재삽입 건너뜀")
        return {"status": "no_crystal", "target_id": target_id,
                "n_processed": 0, "n_skipped": 0, "details": []}

    rank_pdbs = sorted([p for p in glob.glob(os.path.join(af2_dir, "*.pdb"))
                        if "rank_001" in os.path.basename(p)])
    if not rank_pdbs:
        print(f"  [R-17 G2][Info] rank_001 PDB 없음 ({af2_dir}) — no-op")
        return {"status": "no_inputs", "target_id": target_id,
                "n_processed": 0, "n_skipped": 0, "details": []}

    backup_dir = os.path.join(af2_dir, "_pre_reinsert_backup")
    os.makedirs(backup_dir, exist_ok=True)
    # Rejected designs (e.g., G2 Kabsch RMSD > RED threshold) are moved here
    # so downstream globs on ``af2_dir/*.pdb`` exclude them automatically.
    rejected_dir = os.path.join(af2_dir, "_rejected_g2_reinsertion")

    cof_target_chain = target_card.get("target_chain") or target_chain or "A"
    print(f"  [R-17 G2] cofactor 재삽입 — target_id={target_id}, chain={cof_target_chain}, "
          f"declared={[e.get('resname') for e in required_entries]}")

    # Lazy-import so the soft-reject path compiles even if upstream logic
    # changes the import surface. `CofactorReinsertionError` is the ONLY
    # exception we treat as fail-soft; other failures (IO, env, bugs)
    # still propagate.
    try:
        from utils.cofactor_errors import CofactorReinsertionError as _G2_REJECT_EXC
    except Exception:  # pragma: no cover — fall back to a sentinel we never raise
        class _G2_REJECT_EXC(Exception):  # type: ignore[no-redef]
            pass

    details: List[Dict[str, Any]] = []
    rejected: List[Dict[str, Any]] = []
    n_processed = 0
    n_skipped = 0
    n_rejected = 0
    for pdb_path in rank_pdbs:
        base = os.path.basename(pdb_path)
        backup_path = os.path.join(backup_dir, base)
        # Idempotent skip — backup 존재 = 이미 재삽입 완료.
        if os.path.exists(backup_path):
            details.append({"design": base, "status": "skipped_cached"})
            n_skipped += 1
            continue
        # Preserve the pristine AF2 output before overwriting.
        shutil.copy2(pdb_path, backup_path)
        # Reinsert HETATM into in-place file (atomic: write to .tmp → replace).
        tmp_out = pdb_path + ".reinsert.tmp"
        try:
            report = reinsert_cofactors(
                af2_pdb=backup_path,
                crystal_pdb=crystal_pdb,
                target_card=target_card,
                out_pdb=tmp_out,
                target_chain_override=cof_target_chain,
            )
            os.replace(tmp_out, pdb_path)
            details.append({
                "design": base,
                "status": "ok",
                "rmsd_angstrom": report.get("rmsd_angstrom"),
                "n_ca_pairs": report.get("n_ca_pairs"),
                "inserted": [m.get("resname") for m in report.get("inserted_cofactors", [])],
            })
            n_processed += 1
        except _G2_REJECT_EXC as _g2_err:
            # G2 filter rejection — soft-exclude this design so the pipeline
            # continues with the remaining valid candidates. The failing
            # design is moved out of the active af2_dir glob and logged.
            if os.path.exists(tmp_out):
                try:
                    os.remove(tmp_out)
                except OSError:
                    pass
            os.makedirs(rejected_dir, exist_ok=True)
            reject_target = os.path.join(rejected_dir, base)
            try:
                # pdb_path was overwritten by shutil.copy2; move the original
                # AF2 PDB (still in backup_path) into the rejected bin and
                # remove the active file so downstream globs skip it.
                shutil.copy2(backup_path, reject_target)
            except (OSError, shutil.Error):
                pass
            try:
                os.remove(pdb_path)
            except OSError:
                pass
            reason = str(_g2_err)
            # Record the reason + any structural signal we have (RMSD, pairs)
            # from the exception message. The full diagnostic lives in
            # `rejected_designs.json` written below.
            rejected.append({
                "design": base,
                "status": "rejected_g2",
                "reason": reason,
                "rejected_pdb": reject_target,
                "backup_pdb": backup_path,
            })
            details.append({"design": base, "status": "rejected_g2", "reason": reason})
            n_rejected += 1
            print(f"  [R-17 G2][REJECT] {base} — {reason}")
            continue
        except Exception:
            # Non-G2 failure (IO error, unexpected bug): restore the original
            # AF2 PDB so a resume can retry cleanly, then re-raise.
            if os.path.exists(tmp_out):
                try:
                    os.remove(tmp_out)
                except OSError:
                    pass
            try:
                shutil.copy2(backup_path, pdb_path)
            except (OSError, shutil.Error):
                pass
            raise

    # Persist a small per-run log alongside the backups for traceability.
    report_path = os.path.join(backup_dir, "reinsert_report.json")
    try:
        with open(report_path, "w", encoding="utf-8") as f:
            json.dump({
                "target_id": target_id,
                "target_chain": cof_target_chain,
                "crystal_pdb": os.path.abspath(crystal_pdb),
                "n_processed": n_processed,
                "n_skipped": n_skipped,
                "n_rejected": n_rejected,
                "details": details,
            }, f, indent=2)
    except (IOError, OSError) as _rep_err:
        log(f"[R-17 G2] report 기록 실패 (무시): {_rep_err}", level="warning")

    # Separate rejected-designs manifest (downstream consumers can read this
    # to know which design_ids were filtered out and why).
    if rejected:
        try:
            reject_manifest = os.path.join(rejected_dir, "rejected_designs.json")
            with open(reject_manifest, "w", encoding="utf-8") as f:
                json.dump({
                    "target_id": target_id,
                    "gate": "R-17 G2 cofactor reinsertion",
                    "rejected_at": _datetime.datetime.now().isoformat(timespec="seconds"),
                    "rejections": rejected,
                }, f, indent=2)
            print(f"  [R-17 G2][REJECT-MANIFEST] {reject_manifest}")
        except (IOError, OSError, AttributeError) as _rm_err:
            log(f"[R-17 G2] rejection manifest 기록 실패 (무시): {_rm_err}", level="warning")

    status_tag = "skip_cached" if (n_processed == 0 and n_skipped > 0) else "ok"
    print(f"  [R-17 G2] 완료 — processed={n_processed}, "
          f"skipped_cached={n_skipped}, rejected={n_rejected}")
    return {"status": status_tag, "target_id": target_id,
            "n_processed": n_processed, "n_skipped": n_skipped,
            "n_rejected": n_rejected, "details": details,
            "rejected": rejected}


def _filter_top_rank(af2_dir):
    """AF2 결과에서 rank_001 PDB만 남기고 나머지를 _lower_ranks/로 이동한다.

    ColabFold는 각 서열에 대해 5개 모델(rank_001~005)을 생성한다.
    MD 이하 단계에서는 rank_001(ipTM+pTM 최고)만 사용하면 충분하다.
    rank_002~005는 삭제하지 않고 보관하여 필요 시 복원 가능.
    """
    lower_dir = os.path.join(af2_dir, "_lower_ranks")

    if os.path.exists(lower_dir):
        print(f"  [Skip] rank 필터링 이미 완료됨")
        return

    all_pdbs = glob.glob(os.path.join(af2_dir, "*.pdb"))

    lower = [p for p in all_pdbs
             if "rank_00" in os.path.basename(p)
             and "rank_001" not in os.path.basename(p)]

    if not lower:
        print(f"  [Info] 하위 rank PDB 없음 (이미 rank_001만 존재)")
        return

    os.makedirs(lower_dir, exist_ok=True)
    moved = 0
    for pdb_path in lower:
        dest = os.path.join(lower_dir, os.path.basename(pdb_path))
        if not os.path.exists(dest):
            shutil.move(pdb_path, dest)
            moved += 1

    for pdb_path in lower:
        json_candidates = [
            pdb_path.replace(".pdb", ".json"),
            pdb_path.replace("_unrelaxed_", "_scores_").replace(".pdb", ".json"),
        ]
        for jpath in json_candidates:
            if os.path.exists(jpath):
                shutil.move(jpath, os.path.join(lower_dir, os.path.basename(jpath)))

    remaining = len(glob.glob(os.path.join(af2_dir, "*.pdb")))
    rank_001_count = len([p for p in glob.glob(os.path.join(af2_dir, "*.pdb"))
                          if "rank_001" in os.path.basename(p)])
    print(f"  [Rank Filter] {moved}개 하위 rank PDB → _lower_ranks/ 이동. "
          f"rank_001 {rank_001_count}개 유지. (MD 대상: {remaining}개)")


def _filter_by_iptm(af2_dir, min_iptm=0.3):
    # type: (str, float) -> None
    """AF2 랭킹 결과에서 ipTM 임계값 이상인 디자인만 남긴다.

    af2_ranking.csv를 읽어 ipTM < min_iptm인 PDB를 _low_iptm/ 로 이동한다.
    ipTM 필터 후 후보가 0개이면 임계값을 0.2로 자동 완화한다.

    과학적 근거:
    - ipTM < 0.3인 구조는 결합 모드 자체가 불확실 (AF2 자신 없음)
    - 이런 구조의 MD + MM-GBSA 결과는 과학적 신뢰도가 낮음
    - Watson et al. 2023 (RFdiffusion) 에서도 ipTM 기반 필터링 사용

    비파괴적: PDB 삭제 없이 _low_iptm/ 디렉토리로 이동.
    멱등: _low_iptm/ 이미 존재 시 스킵.
    """
    low_dir = os.path.join(af2_dir, "_low_iptm")

    if os.path.exists(low_dir):
        print(f"  [Skip] ipTM 필터링 이미 완료됨")
        return

    ranking_csv = os.path.join(af2_dir, "af2_ranking.csv")
    if not os.path.exists(ranking_csv):
        print(f"  [Warning] af2_ranking.csv 없음 — ipTM 필터 건너뜀")
        return

    # CSV 읽기
    with open(ranking_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        print(f"  [Warning] af2_ranking.csv 비어있음")
        return

    # ipTM 칼럼명 탐색 (iptm, ipTM, i_ptm 등)
    iptm_col = None
    for col_name in ["ipTM", "iptm", "i_ptm", "interface_ptm"]:
        if col_name in rows[0]:
            iptm_col = col_name
            break

    if iptm_col is None:
        print(f"  [Warning] af2_ranking.csv에 ipTM 칼럼 없음 — 필터 건너뜀")
        return

    # 통과/탈락 분류
    pass_ids = []   # type: list
    fail_ids = []   # type: list
    for row in rows:
        try:
            iptm_val = float(row[iptm_col])
        except (ValueError, TypeError):
            continue
        design_id = row.get("id", row.get("design_id", row.get("name", "")))
        if iptm_val >= min_iptm:
            pass_ids.append((design_id, iptm_val))
        else:
            fail_ids.append((design_id, iptm_val))

    # 후보 0개이면 임계값 완화
    if not pass_ids and fail_ids:
        fallback_iptm = 0.2
        print(f"  [ipTM Filter] 임계값 {min_iptm}에서 통과 후보 0개 "
              f"→ {fallback_iptm}으로 완화")
        pass_ids = [(did, v) for did, v in fail_ids if v >= fallback_iptm]
        fail_ids = [(did, v) for did, v in fail_ids if v < fallback_iptm]

    if not pass_ids:
        print(f"  [Warning] ipTM 완화 후에도 후보 0개 — 필터 건너뜀")
        return

    # 탈락 PDB를 _low_iptm/ 로 이동
    os.makedirs(low_dir, exist_ok=True)

    fail_design_ids = {did for did, _ in fail_ids}
    moved = 0
    for pdb_path in glob.glob(os.path.join(af2_dir, "*.pdb")):
        basename = os.path.basename(pdb_path)
        for did in fail_design_ids:
            if basename.startswith(did):
                dest = os.path.join(low_dir, basename)
                if not os.path.exists(dest):
                    shutil.move(pdb_path, dest)
                    moved += 1
                break

    # 결과 리포트
    print(f"  [ipTM Filter] 임계값: {min_iptm}")
    print(f"    통과: {len(pass_ids)}개 (MD 진행)")
    for did, v in sorted(pass_ids, key=lambda x: -x[1])[:5]:
        print(f"      {did}: ipTM={v:.3f}")
    if len(pass_ids) > 5:
        print(f"      ... 외 {len(pass_ids) - 5}개")
    print(f"    탈락: {len(fail_ids)}개 → _low_iptm/ 이동 ({moved}개 PDB)")

    # 필터 리포트 CSV 저장
    filter_report = os.path.join(af2_dir, "iptm_filter_report.csv")
    with open(filter_report, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["design_id", "ipTM", "status"])
        for did, v in sorted(pass_ids, key=lambda x: -x[1]):
            w.writerow([did, f"{v:.4f}", "PASS"])
        for did, v in sorted(fail_ids, key=lambda x: -x[1]):
            w.writerow([did, f"{v:.4f}", "FAIL"])
    print(f"    리포트: {filter_report}")


def _filter_by_rmsd(af2_dir, crystal_pdb, target_chain="A", rmsd_threshold_a=None):
    # type: (str, str, str, float) -> None
    """AF2 rank_001 PDB 의 target chain Cα RMSD (vs crystal) 로 필터링한다.

    Fail-soft: 임계값 초과 design 은 ``_high_rmsd/`` 로 이동하여 downstream
    glob 에서 제외된다. 이미 rank + ipTM 필터로 좁혀진 후보들에 대해
    예측 구조적 타당성 (target backbone 정합도) 을 추가 필터로 적용한다.

    과학적 근거:
    - AF2-multimer chain A RMSD 2.5 Å 초과 = AF2 모델의 target fold 예측이
      crystal 대비 크게 벗어남 → 이 상태로 MD/QM/MM 진행 시 cofactor
      transplant (G2) 나 charge topology 분석이 불안정
    - 2.0-2.5 Å 구간은 normal AF2 noise (Bryant 2022 Nat Commun IQR
      1.5-3.2); 그 이상은 outlier → pre-filter 단계에서 제외하여
      downstream compute (MD, SCF) 낭비 방지

    비파괴적: PDB 삭제 없이 _high_rmsd/ 디렉토리로 이동.
    멱등: _high_rmsd/ 이미 존재 시 스킵.
    """
    high_dir = os.path.join(af2_dir, "_high_rmsd")
    if os.path.exists(high_dir):
        print(f"  [Skip] RMSD 필터링 이미 완료됨")
        return
    if not crystal_pdb or not os.path.exists(crystal_pdb):
        print(f"  [Info] crystal PDB 부재 ({crystal_pdb}) — RMSD 필터 건너뜀")
        return

    # Reuse Kabsch + CA extraction from the G2 module (DRY).
    try:
        from utils.reinsert_cofactors_post_af2 import (
            DEFAULT_RMSD_THRESHOLD_A as _G2_DEFAULT_THRESHOLD,
            MIN_CA_PAIRS as _G2_MIN_CA,
            _extract_ca_by_chain_resnum as _extract_ca,
            _kabsch as _kabsch_fn,
        )
    except ImportError as _ie:
        print(f"  [Warning] reinsert module import 실패 ({_ie}) — RMSD 필터 건너뜀")
        return

    threshold = float(rmsd_threshold_a) if rmsd_threshold_a is not None else float(_G2_DEFAULT_THRESHOLD)
    ca_cry = _extract_ca(crystal_pdb, target_chain)
    if not ca_cry:
        print(f"  [Warning] crystal chain '{target_chain}' Cα 없음 — RMSD 필터 건너뜀")
        return

    rank_pdbs = sorted([p for p in glob.glob(os.path.join(af2_dir, "*.pdb"))
                        if "rank_001" in os.path.basename(p)])
    if not rank_pdbs:
        print(f"  [Info] rank_001 PDB 없음 — RMSD 필터 no-op")
        return

    import numpy as _np
    pass_designs = []   # type: list
    fail_designs = []   # type: list
    for pdb_path in rank_pdbs:
        base = os.path.basename(pdb_path)
        design_id = base.split("_unrelaxed")[0]
        ca_af2 = _extract_ca(pdb_path, target_chain)
        shared = sorted(set(ca_cry) & set(ca_af2))
        if len(shared) < _G2_MIN_CA:
            # Too few pairs — cannot reliably RMSD-filter this design.
            # Err on the side of keeping it (G2 itself will enforce again).
            pass_designs.append((design_id, float("nan"), len(shared), "under_min_pairs"))
            continue
        P = _np.array([ca_cry[r] for r in shared])
        Q = _np.array([ca_af2[r] for r in shared])
        _, _, rmsd = _kabsch_fn(P, Q)
        if rmsd <= threshold:
            pass_designs.append((design_id, float(rmsd), len(shared), "pass"))
        else:
            fail_designs.append((design_id, float(rmsd), len(shared), pdb_path))

    if not fail_designs:
        print(f"  [RMSD Filter] 임계값: {threshold:.2f} Å | 통과: {len(pass_designs)}개 "
              f"(전부 통과, no reject)")
    else:
        os.makedirs(high_dir, exist_ok=True)
        moved = 0
        for design_id, rmsd, n_pairs, pdb_path in fail_designs:
            basename = os.path.basename(pdb_path)
            dest = os.path.join(high_dir, basename)
            if not os.path.exists(dest):
                try:
                    shutil.move(pdb_path, dest)
                    moved += 1
                except (OSError, shutil.Error):
                    pass
            # companion score json (if exists) goes with the PDB
            json_candidates = [
                pdb_path.replace(".pdb", ".json"),
                pdb_path.replace("_unrelaxed_", "_scores_").replace(".pdb", ".json"),
            ]
            for jpath in json_candidates:
                if os.path.exists(jpath):
                    try:
                        shutil.move(jpath, os.path.join(high_dir, os.path.basename(jpath)))
                    except (OSError, shutil.Error):
                        pass
        print(f"  [RMSD Filter] 임계값: {threshold:.2f} Å")
        print(f"    통과: {len(pass_designs)}개")
        print(f"    탈락: {len(fail_designs)}개 → _high_rmsd/ 이동 ({moved}개 PDB)")
        for did, rmsd, n_pairs, _p in sorted(fail_designs, key=lambda x: -x[1])[:5]:
            print(f"      {did}: RMSD={rmsd:.3f} Å (n={n_pairs})")

    # Unified report (pass + fail).
    filter_report = os.path.join(af2_dir, "rmsd_filter_report.csv")
    try:
        with open(filter_report, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["design_id", "rmsd_angstrom", "n_ca_pairs", "status"])
            for did, rmsd, n_pairs, note in sorted(pass_designs, key=lambda x: x[1] if x[1] == x[1] else 0.0):
                status = "PASS" if note != "under_min_pairs" else "PASS_WARN_FEW_PAIRS"
                rmsd_str = f"{rmsd:.4f}" if rmsd == rmsd else "nan"
                w.writerow([did, rmsd_str, n_pairs, status])
            for did, rmsd, n_pairs, _p in sorted(fail_designs, key=lambda x: -x[1]):
                w.writerow([did, f"{rmsd:.4f}", n_pairs, "FAIL"])
        print(f"    리포트: {filter_report}")
    except (IOError, OSError) as _rep_err:
        log(f"[RMSD Filter] report 기록 실패 (무시): {_rep_err}", level="warning")


def _run_pre_filter(af2_dir, min_iptm=0.3, crystal_pdb=None, target_chain="A",
                     rmsd_threshold_a=None):
    # type: (str, float, str, str, float) -> None
    """통합 AF2 pre-filter: Rank → ipTM → RMSD (이 순서).

    각 단계는 fail-soft (임계값 미달 design 은 전용 서브디렉토리로 이동;
    downstream glob 에서 자동 제외). Pipeline 은 전체적으로 멈추지 않으며
    통과한 design 들로 계속 진행한다. 과거 분리돼 있던 rank 필터 +
    ipTM 필터를 하나로 묶고, RMSD 필터를 추가한다 — AF2 prediction 이
    crystal 대비 구조적으로 크게 벗어난 design 을 downstream compute
    (MD, QM/MM) 투입 이전에 솎아낸다.
    """
    _filter_top_rank(af2_dir)
    _filter_by_iptm(af2_dir, min_iptm=min_iptm)
    _filter_by_rmsd(
        af2_dir,
        crystal_pdb=crystal_pdb,
        target_chain=target_chain or "A",
        rmsd_threshold_a=rmsd_threshold_a,
    )


def _run_ncaa_mutation(af2_dir, ncaa_code, ncaa_label, positions, target_out_dir, binder_chain="B", target_chain="A"):
    if not ncaa_code or ncaa_code == "none":
        print("  [Info] Canonical Aminoacid 모드이므로 ncAA 치환 단계를 건너뜁니다.")
        return af2_dir
        
    ncaadir = os.path.join(target_out_dir, "ncaa_mutants")
    manifest_path = os.path.join(ncaadir, "mutation_manifest.json")
    
    manifest = load_json_manifest(manifest_path)
    if manifest and manifest.get("ncaa_code") == ncaa_code and manifest.get("input_positions") == positions:
        print(f"  [Skip] ncAA 치환 이미 완료됨 (Manifest 일치)")
        return ncaadir

    os.makedirs(ncaadir, exist_ok=True)
    args = [
        "--inputdir", af2_dir,
        "--outputdir", ncaadir,
        "--ncaa_code", ncaa_code,
        "--residues", positions,
        "--binder_chain", binder_chain,
        "--target_chain", target_chain,
    ]
    run_conda_command(SCRIPT_PATHS["ncaa_mutate"], args, "ncAA Mutation", env_key="ncaa")
    return ncaadir

def _extract_cofactor_source_pdbs(
    crystal_pdb: str,
    resnames,
    out_dir: str,
) -> dict:
    """Extract per-resname HETATM/ATOM blocks from ``crystal_pdb`` and write
    one minimal PDB per residue into ``out_dir``.

    Returns a ``{resname_upper: abs_pdb_path}`` dict (resnames with zero
    matching atoms are omitted).
    """
    if not crystal_pdb or not os.path.exists(crystal_pdb):
        return {}
    targets = {str(r).strip().upper() for r in resnames if r}
    if not targets:
        return {}

    os.makedirs(out_dir, exist_ok=True)
    buckets: dict = {r: [] for r in targets}
    with open(crystal_pdb, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not (line.startswith("HETATM") or line.startswith("ATOM")):
                continue
            if len(line) < 20:
                continue
            rn = line[17:20].strip().upper()
            if rn in buckets:
                buckets[rn].append(line)

    out_map: dict = {}
    for rn, lines in buckets.items():
        if not lines:
            continue
        out_path = os.path.join(out_dir, f"{rn}_source.pdb")
        with open(out_path, "w", encoding="utf-8") as f:
            f.writelines(lines)
            f.write("END\n")
        out_map[rn] = os.path.abspath(out_path)
    return out_map


def _run_cofactor_parameterization(target_id: str, crystal_pdb: str = "") -> dict:
    """[R-17 G3] Ensure cofactor force-field parameters are cached.

    Dispatches ``utils/parameterize_cofactor.ensure_cofactor_params`` for
    cofactor-present targets. Behavior:

    - ``target_id`` empty / card missing / no required cofactors → no-op,
      returns ``{"status": "noop"}``.
    - All declared cofactors already cached (frcmod + mol2 exist) or are
      AMBER-shipped ions (Li-Merz 12-6-4) → logs "cached/ion_builtin" per
      resname and returns ``{"status": "ok", ...}``.
    - Cache miss and ``UPDD_ALLOW_ANTECHAMBER=1`` → regeneration attempt
      (opt-in). On success, files land under ``params/<resname>/``.
      Per-cofactor source PDBs are extracted from ``crystal_pdb`` into
      ``<project>/params/_sources/`` and passed to
      ``ensure_cofactor_params(source_pdbs=...)``.
    - Cache miss and opt-in off → emit a WARNING (no halt). G4 (MD setup)
      will raise the hard error if the parameters are actually needed.

    Args:
        target_id: target_card stem (e.g. "6WGN"). Empty → no-op.
        crystal_pdb: Crystal/preprocessed PDB path used to extract per-
            cofactor source coordinates for antechamber. Required when
            regeneration is needed; ignored when every cofactor is already
            cached or is an ion_builtin.

    Returns:
        Dict with ``status`` in {``noop``, ``ok``, ``warn_missing``,
        ``error_unexpected``} plus the per-resname report under ``report``.
    """
    if not target_id:
        return {"status": "noop", "reason": "no_target_id"}

    _card_path = os.path.join(UPDD_DIR, "target_cards", f"{target_id}.json")
    if not os.path.exists(_card_path):
        return {"status": "noop", "reason": "no_card"}

    if UTILS_DIR not in sys.path:
        sys.path.append(UTILS_DIR)
    try:
        from utils_common import load_target_card  # type: ignore
        from parameterize_cofactor import ensure_cofactor_params  # type: ignore
        from cofactor_errors import CofactorParamMissingError  # type: ignore
    except ImportError as _imp_err:
        print(f"  [R-17 G3][WARN] module import 실패 ({_imp_err}) — skip")
        return {"status": "warn_missing", "reason": f"import_error: {_imp_err}"}

    _cards_dir = os.path.join(UPDD_DIR, "target_cards")
    try:
        target_card = load_target_card(target_id, cards_dir=_cards_dir)
    except (FileNotFoundError, ValueError) as _load_err:
        print(f"  [R-17 G3][WARN] target_card 로드 실패 ({_load_err}) — skip")
        return {"status": "warn_missing", "reason": f"card_load_error: {_load_err}"}

    required = [e for e in (target_card.get("cofactor_residues") or [])
                if isinstance(e, dict) and e.get("required")]
    if not required:
        return {"status": "noop", "reason": "no_required_cofactors"}

    # Build source_pdbs dict from crystal PDB for any required cofactor that
    # needs antechamber regeneration. Safe to compute unconditionally: if
    # every cofactor is already cached, ensure_cofactor_params skips the
    # regen branch and the extracted stubs are simply unused.
    source_pdbs: dict = {}
    if crystal_pdb and os.path.exists(crystal_pdb):
        _src_dir = os.path.join(UPDD_DIR, "params", "_sources")
        _resnames = [str(e.get("resname", "")).strip().upper() for e in required]
        source_pdbs = _extract_cofactor_source_pdbs(crystal_pdb, _resnames, _src_dir)

    # antechamber/parmchk2 ship inside conda envs (ncaa/md_simulation/qmmm),
    # not on the orchestrator's base PATH. Inject the ncaa env bin dir so
    # ensure_cofactor_params() can find them when UPDD_ALLOW_ANTECHAMBER=1.
    import shutil as _shutil
    if not _shutil.which("antechamber"):
        for _env in ("ncaa", "md_simulation", "qmmm"):
            _bin = f"/home/san/miniconda3/envs/{_env}/bin"
            if os.path.isfile(os.path.join(_bin, "antechamber")):
                os.environ["PATH"] = _bin + os.pathsep + os.environ.get("PATH", "")
                break

    try:
        report = ensure_cofactor_params(target_card, source_pdbs=source_pdbs)
    except CofactorParamMissingError as _param_err:
        # Do NOT halt — G4 will raise the hard error if the params are
        # actually consumed. Here we only emit an operator-visible WARN so
        # the user can pre-supply params/<resname>/{*.frcmod,*.mol2} or set
        # UPDD_ALLOW_ANTECHAMBER=1 before MD setup.
        print(f"  [R-17 G3][WARN] GNP/cofactor params missing: {_param_err}. "
              f"Set UPDD_ALLOW_ANTECHAMBER=1 or ship params/<resname>/* "
              f"before Step 9 MD, otherwise G4 will halt.")
        return {"status": "warn_missing", "reason": str(_param_err)}

    for entry in report:
        resname = entry.get("resname", "?")
        status = entry.get("status", "?")
        method = entry.get("method", "?")
        if status == "cache_hit":
            print(f"  [R-17 G3] {resname}: cached, skip ({method})")
        elif status == "ion_builtin":
            print(f"  [R-17 G3] {resname}: ion_builtin ({method}) — AMBER-shipped")
        elif status == "regenerated":
            print(f"  [R-17 G3] {resname}: regenerated via antechamber AM1-BCC "
                  f"(mol2={entry.get('mol2')})")
        elif status == "regenerated_bcc_fallback":
            print(f"  [R-17 G3] {resname}: regenerated via antechamber AM1-BCC "
                  f"FALLBACK ({method} requested external RESP .mol2; supply "
                  f"ff_parameters.resp_mol2 before v0.7 calibration). "
                  f"mol2={entry.get('mol2')}")
        elif status == "regenerated_resp_external":
            print(f"  [R-17 G3] {resname}: derived frcmod from user RESP .mol2 "
                  f"via parmchk2 (mol2={entry.get('mol2')})")
        else:
            print(f"  [R-17 G3] {resname}: status={status}")
    return {"status": "ok", "report": report}


def _run_parameterization(struct_dir, param_dir, ncaa_def, target_id: str = "",
                          crystal_pdb: str = ""):
    # [R-17 G3] Before ncAA-specific parameterization, verify cofactor FF
    # parameters are cached for cofactor-present targets. No-op for
    # cofactor-absent targets (target_id empty, card missing, or empty
    # cofactor_residues). Non-blocking — WARNs only when cache is missing
    # and antechamber opt-in is off; G4 raises later if params are needed.
    # crystal_pdb enables on-demand antechamber regeneration by supplying
    # per-cofactor source coordinates extracted from the preprocessed PDB.
    if target_id:
        _run_cofactor_parameterization(target_id, crystal_pdb=crystal_pdb)

    if ncaa_def is None:
        print("  [Info] 야생형(Wild-Type) 모드이므로 파라미터화 단계를 건너뜁니다.")
        return ""

    manifest_path = os.path.join(param_dir, f"{ncaa_def.code}_params_manifest.json")
    manifest = load_json_manifest(manifest_path)
    if manifest and manifest.get("ncaa_code") == ncaa_def.code:
        print(f"  [Skip] 파라미터화 이미 완료됨 (Manifest 일치)")
        return manifest_path

    os.makedirs(param_dir, exist_ok=True)
    # [CRITICAL FIX] ncaa_def 객체가 아닌 ncaa_def.code 문자열을 넘기도록 수정
    run_conda_command(SCRIPT_PATHS["parameterize"], ["--inputdir", struct_dir, "--outputdir", param_dir, "--ncaa_code", ncaa_def.code], "ncAA Parameterization", env_key="ncaa")
    return manifest_path

def _run_admet(struct_dir, admet_mode, target_out_dir, status, status_file):
    if admet_mode in ["cell", "bbb"]:
        print_step("Step 8.5: ADMET (투과성) 톨게이트 필터링 진행")
        if "admet" in status.get("completed_steps", []):
            print("  ✅ [Skip] ADMET 필터링 작업은 이미 완료되어 건너뜁니다.")
            return

        if UTILS_DIR not in sys.path:
            sys.path.append(UTILS_DIR)
        # [Refactor] 선택적 의존성(RDKit / admet_filter)은 환경별 설치 여부가
        # 다르므로 함수 내부에 try/except ImportError 로 래핑한 채 유지한다.
        try:
            from admet_filter import check_admet_rules
            from rdkit import Chem
        except ImportError as e:
            print(f"  [!] RDKit 모듈 오류로 인해 필터를 스킵합니다: {e}")
            return

        valid_pdbs: List[str] = []
        pdb_list = glob.glob(os.path.join(struct_dir, "*.pdb"))
        for pdb_file in pdb_list:
            temp_pdb = pdb_file.replace(".pdb", "_temp_binder.pdb")
            with open(pdb_file, "r", encoding="utf-8") as f_in, open(temp_pdb, "w", encoding="utf-8") as f_out:
                for line in f_in:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        if line[21] == "B":
                            f_out.write(line)
                    elif line.startswith("CONECT") or line.startswith("END"):
                        f_out.write(line)

            mol = Chem.MolFromPDBFile(temp_pdb, sanitize=False)
            if os.path.exists(temp_pdb):
                os.remove(temp_pdb)

            if mol:
                try:
                    mol.UpdatePropertyCache(strict=False)
                    Chem.FastFindRings(mol)
                except Exception as e:
                    print(f"  [!] mol 전처리 실패, 스킵: {os.path.basename(pdb_file)} ({e})")
                    continue
                is_pass, msg = check_admet_rules(mol, mode=admet_mode)
                if is_pass:
                    valid_pdbs.append(pdb_file)
                    print(f"  [PASS] {os.path.basename(pdb_file)[:30]}... : {msg}")
                else:
                    print(f"  [FAIL] {os.path.basename(pdb_file)[:30]}... : {msg}")
                    os.rename(pdb_file, pdb_file + ".fail")
        if not valid_pdbs:
            sys.exit("\n  [!] 모든 바인더가 ADMET 필터링을 통과하지 못했습니다! 파이프라인을 종료합니다.")
        else:
            print(f"\n  -> 총 {len(valid_pdbs)}/{len(pdb_list)}개의 바인더가 통과되어 MD로 이동합니다.")
            if "admet" not in status.setdefault("completed_steps", []):
                status["completed_steps"].append("admet")
            with open(status_file, "w", encoding="utf-8") as f:
                json.dump(status, f, indent=4)

def _run_md_and_snapshots(struct_dir, params_manifest, target_out_dir, ncaa_label, ncaa_code, canonical_topo, binder_chain, graph_policy="strict", dcd_backup_path="", target_id: str = ""):
    md_out_dir = os.path.join(target_out_dir, "mdresult")
    snap_dir = os.path.join(target_out_dir, "snapshots")
    manifest_path = os.path.join(md_out_dir, "md_manifest.json")

    manifest = load_json_manifest(manifest_path)
    ncaa_code_check = ncaa_code if ncaa_code else "none"

    if manifest and manifest.get("ncaa_code", "none") == ncaa_code_check and manifest.get("topology") == canonical_topo:
         print(f"  [Skip] MD 수행 이미 완료됨 (Manifest 일치, Status: {manifest.get('status')})")
    else:
        os.makedirs(md_out_dir, exist_ok=True)
        # [v4 H-3] MD step 수 기본값을 100ps(50,000) → 10ns(5,000,000) 로 상향한다.
        # 과학적 근거: 펩타이드 결합 시뮬레이션에서 신뢰할 만한 MM-GBSA 값을 얻으려면
        # 최소 10ns 가 필요하다 (Sun et al. JCIM 2014). 이 변경은 H-1 에서 활성화한
        # CUDA mixed-precision 가속 (RTX 5070 Ti 기준 50-100x 가속) 을 전제로 한다.
        # 환경 변수 `UPDD_MD_STEPS` 로 운영자가 override 할 수 있다 (CPU-only 환경
        # 또는 빠른 반복 실험 시 50000 으로 회귀 가능).
        md_steps = os.environ.get("UPDD_MD_STEPS", "5000000")
        md_args = [
            "--inputdir", struct_dir,
            "--outputdir", md_out_dir,
            "--steps", md_steps,
            "--topology", canonical_topo,
            "--binder_chain", binder_chain,
            "--graph_policy", graph_policy,
        ]
        if dcd_backup_path and os.path.exists(dcd_backup_path):
            md_args.extend(["--hdd_path", dcd_backup_path])

        # [R-17 G4] Forward target_id so run_restrained_md registers cofactor
        # GAFF templates (GNP/ATP/NAD/…) on the ForceField before addHydrogens
        # + createSystem. Silent no-op when the target_card is absent.
        if target_id:
            _card_path = os.path.join(UPDD_DIR, "target_cards", f"{target_id}.json")
            if os.path.exists(_card_path):
                md_args.extend(["--target_id", target_id])

        if not ncaa_code or ncaa_code == "none":
            md_args += ["--ncaa_label", "none", "--ncaa_code", ""]
        else:
            md_args += ["--params_manifest", params_manifest, "--ncaa_label", ncaa_label, "--ncaa_code", ncaa_code]
            
        md_args += ["--dt_fs", "2.0"]
        # Extract current seed from md_args if present (default 42)
        try:
            _seed_idx = md_args.index("--seed")
            _seed_val = int(md_args[_seed_idx + 1])
        except (ValueError, IndexError):
            _seed_val = 42
            _seed_idx = None
        run_conda_command(SCRIPT_PATHS["run_md"], md_args, "Restrained MD (Pass 1: 2fs)", env_key="md")

        # [v0.6.7] Probabilistic-NaN retry: before falling back to dt=1fs,
        # try dt=2fs with a different seed. Seed changes affect Langevin
        # random-number stream; some seed+state combinations trigger NaN
        # on the first integration step even after a clean minimization.
        # Empirically: 1 in 3-5 random draws; a single re-seed usually recovers.
        exploded_after_pass1 = glob.glob(os.path.join(md_out_dir, "*_EXPLODED_dt2fs.log"))
        all_input_pdbs_p1 = glob.glob(os.path.join(struct_dir, "*.pdb"))
        failed_no_final_p1 = [
            os.path.splitext(os.path.basename(p))[0] for p in all_input_pdbs_p1
            if not os.path.exists(os.path.join(md_out_dir,
                                               os.path.splitext(os.path.basename(p))[0] + "_final.pdb"))
        ]
        reseed_targets = set()
        for elog in exploded_after_pass1:
            reseed_targets.add(os.path.basename(elog).replace("_EXPLODED_dt2fs.log", ""))
        for stem in failed_no_final_p1:
            reseed_targets.add(stem)

        if reseed_targets:
            print(f"\n  [Pass 1b] {len(reseed_targets)}개 디자인 2fs 폭발 — seed 변경 재시도 먼저 시도...")
            _archive_ts_1b = _datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            _archive_dir_1b = os.path.join(md_out_dir, "_archive", f"{_archive_ts_1b}_pre_reseed")
            os.makedirs(_archive_dir_1b, exist_ok=True)
            _archive_patterns_1b = [
                "_EXPLODED_dt2fs.log", "_EXPLODED_MD.log", "_EXPLODED_NVT.log",
                "_EXPLODED_COLD_NVT.log", "_EXPLODED_MIN.log",
                "_EXPLODED_CYCLIC_MIN.log", "_EXPLODED_CYCLIC_RELAX.log",
                "_partial_md.pdb", "_partial_nvt.pdb",
                "_final.pdb", "_md.log", "_restrained.dcd",
            ]
            for stem in reseed_targets:
                for pat in _archive_patterns_1b:
                    src = os.path.join(md_out_dir, stem + pat)
                    if os.path.exists(src) or os.path.islink(src):
                        shutil.move(src, os.path.join(_archive_dir_1b, stem + pat))
                if dcd_backup_path and os.path.isdir(dcd_backup_path):
                    _hdd_dcd = os.path.join(dcd_backup_path, stem + "_restrained.dcd")
                    if os.path.exists(_hdd_dcd):
                        shutil.move(_hdd_dcd,
                                    os.path.join(_archive_dir_1b, stem + "_restrained_HDD.dcd"))
            reseed_args = list(md_args)
            # Prime-offset for deterministic second attempt (not original seed).
            new_seed = _seed_val + 13 if _seed_val != -1 else 7
            if _seed_idx is not None:
                reseed_args[_seed_idx + 1] = str(new_seed)
            else:
                reseed_args += ["--seed", str(new_seed)]
            run_conda_command(SCRIPT_PATHS["run_md"], reseed_args,
                              f"Restrained MD (Pass 1b: 2fs seed={new_seed})", env_key="md")

        # EXPLODED 마커가 있는 것 + _final.pdb가 없는 것 = 재시도 대상 전수
        exploded_2fs = glob.glob(os.path.join(md_out_dir, "*_EXPLODED_dt2fs.log"))
        # _final.pdb가 없는 디자인도 재시도 대상에 포함 (NVT 폭발 등)
        all_input_pdbs = glob.glob(os.path.join(struct_dir, "*.pdb"))
        failed_no_final = []
        for inp_pdb in all_input_pdbs:
            stem = os.path.splitext(os.path.basename(inp_pdb))[0]
            final_candidate = os.path.join(md_out_dir, stem + "_final.pdb")
            if not os.path.exists(final_candidate):
                failed_no_final.append(stem)

        # 재시도 대상: EXPLODED 마커 + final 없는 것 (중복 제거)
        retry_targets = set()
        for elog in exploded_2fs:
            retry_targets.add(os.path.basename(elog).replace("_EXPLODED_dt2fs.log", ""))
        for stem in failed_no_final:
            retry_targets.add(stem)

        n_retry = len(retry_targets)
        if n_retry > 0 and canonical_topo in ("cyclic_htc", "cyclic_nm"):
            print(f"\n  [2-Pass] {n_retry}개 디자인이 2fs에서 실패. 1fs로 재시도합니다.")
            print(f"    (EXPLODED: {len(exploded_2fs)}개, NVT/기타 실패: {len(failed_no_final)}개)")
            # [v0.5 Bug 1 fix] 모든 explosion marker + partial placeholder를 archive로 이동
            _archive_ts = _datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            _archive_dir = os.path.join(md_out_dir, "_archive", _archive_ts)
            os.makedirs(_archive_dir, exist_ok=True)
            _marker_patterns = EXPLOSION_MARKER_PATTERNS if _V05_MD_STATUS_AVAILABLE else [
                "_EXPLODED_dt2fs.log", "_EXPLODED_dt1fs.log",
                "_EXPLODED_MD.log", "_EXPLODED_NVT.log", "_EXPLODED_COLD_NVT.log",
                "_EXPLODED_MIN.log", "_EXPLODED_CYCLIC_MIN.log", "_EXPLODED_CYCLIC_RELAX.log",
            ]
            _partial_patterns = PARTIAL_PLACEHOLDERS if _V05_MD_STATUS_AVAILABLE else [
                "_partial_md.pdb", "_partial_nvt.pdb",
            ]
            # BUGFIX (2026-04-20): Pass 1 explosion path with
            # completion_ratio >= _MIN_COMPLETION_RATIO writes a checkpoint
            # snapshot to ``<stem>_final.pdb`` (PARTIAL_SUCCESS tier). Without
            # archiving that file, the Pass 2 (1fs) run_restrained_md SKIPs
            # the design via ``os.path.exists(out_pdb) and not
            # exploded_markers`` (markers already archived above).  Archive
            # ``_final.pdb`` + ``_md.log`` + the HDD ``_restrained.dcd`` so
            # Pass 2 retry actually re-enters MD.
            _stale_final_patterns = ["_final.pdb", "_md.log",
                                     "_restrained.dcd"]
            for stem in retry_targets:
                for pat in _marker_patterns + _partial_patterns + _stale_final_patterns:
                    src = os.path.join(md_out_dir, stem + pat)
                    if os.path.exists(src) or os.path.islink(src):
                        shutil.move(src, os.path.join(_archive_dir, stem + pat))
                # Also clear the HDD-side partial DCD if present (kept on
                # the external drive by the direct-stream fix).
                if dcd_backup_path and os.path.isdir(dcd_backup_path):
                    _hdd_dcd = os.path.join(dcd_backup_path,
                                            stem + "_restrained.dcd")
                    if os.path.exists(_hdd_dcd):
                        shutil.move(_hdd_dcd, os.path.join(
                            _archive_dir, stem + "_restrained_HDD.dcd"))

            retry_args = list(md_args)
            dt_idx = retry_args.index("--dt_fs")
            retry_args[dt_idx + 1] = "1.0"
            run_conda_command(SCRIPT_PATHS["run_md"], retry_args,
                              f"Restrained MD (Pass 2: 1fs, {n_retry}개 재시도)", env_key="md")

            still_exploded = glob.glob(os.path.join(md_out_dir, "*_EXPLODED_dt1fs.log"))
            if still_exploded:
                print(f"  [2-Pass] 1fs에서도 {len(still_exploded)}개 폭발. "
                      f"체크포인트 복구된 부분 궤적으로 진행합니다.")

        md_manifest = {"status": "SUCCESS_OR_WARNING", "topology": canonical_topo, "graph_policy": graph_policy, "ncaa_code": ncaa_code_check}
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(md_manifest, f, indent=2)

        # [v0.5] _md_status.json 생성/갱신 — per-design MD tier SSOT
        if _V05_MD_STATUS_AVAILABLE:
            md_status = load_md_status(target_out_dir) or {
                "schema_version": "0.5.0",
                "project": os.path.basename(target_out_dir),
                "topology": canonical_topo,
                "created_by": "UPDD.py",
                "designs": {},
            }
            for inp_pdb in all_input_pdbs:
                stem = os.path.splitext(os.path.basename(inp_pdb))[0]
                outcome = classify_md_outcome(md_out_dir, stem)
                md_status["designs"][stem] = outcome
            write_md_status(target_out_dir, md_status)
            _n_tiers = {}  # type: Dict[str, int]
            for d in md_status["designs"].values():
                t = d.get("tier", "UNKNOWN")
                _n_tiers[t] = _n_tiers.get(t, 0) + 1
            print(f"  [v0.5] _md_status.json 갱신: {_n_tiers}")

    if dcd_backup_path and os.path.exists(dcd_backup_path):
        print_step("Step 9.5: Symlink Setup for MDtraj")
        for hdd_dcd in glob.glob(os.path.join(dcd_backup_path, "*.dcd")):
            ssd_dcd = os.path.join(md_out_dir, os.path.basename(hdd_dcd))
            if not os.path.exists(ssd_dcd) and not os.path.islink(ssd_dcd):
                try:
                    os.symlink(hdd_dcd, ssd_dcd)
                except OSError as _ln_err:
                    log(f"[symlink] DCD 링크 생성 실패 (무시): {_ln_err}", level="warning")
        print("  ✅ 스냅샷용 심볼릭 링크 설정이 완료되었습니다.")

    print_step("Step 10: Snapshot Extraction (MDtraj)")
    if len(glob.glob(os.path.join(snap_dir, "**", "*.pdb"), recursive=True)) > 0:
        print(f"  [Skip] Snapshot 추출 이미 완료됨")
    else:
        os.makedirs(snap_dir, exist_ok=True)
        # [v35 AUDIT] Principle 7: argv contract realigned with extract_snapshots.py argparse.
        # Previously passed --inputdir/--num_snapshots/--topology — none of which were declared
        # by the receiving script, causing argparse.ArgumentError on every invocation.
        # Canonical flags: --md_dir, --outputdir, --n_snapshots, --binder_chain.
        snap_args = [
            "--md_dir", md_out_dir,
            "--outputdir", snap_dir,
            "--n_snapshots", "5",
            "--binder_chain", binder_chain,
        ]
        # [R-17 G5] Pass target_id when the card exists so extract_snapshots.py
        # honors cofactor_residues (strip_solvent keeps declared cofactors +
        # per-snapshot presence validation raises CofactorMissingError when
        # the cofactor atoms drop out of the trajectory).
        if target_id:
            _card_path = os.path.join(UPDD_DIR, "target_cards", f"{target_id}.json")
            if os.path.exists(_card_path):
                snap_args += ["--target_id", target_id]
                print(f"  [R-17 G5] target_card '{target_id}' 감지 — cofactor-aware snapshot 추출")
        run_conda_command(SCRIPT_PATHS["extract_snap"], snap_args, "Extract Snapshots", env_key="md")

    if dcd_backup_path and os.path.exists(dcd_backup_path):
        for link in glob.glob(os.path.join(md_out_dir, "*.dcd")):
            if os.path.islink(link):
                os.unlink(link)
        print("  ✅ SSD의 DCD 심볼릭 링크 정리를 완료하였습니다.")

    return md_out_dir, snap_dir

# ==========================================
# [v0.3 #6] VRAM Watchdog 헬퍼 (자동 시작/정리)
# ==========================================
# proactive_oom watchdog 을 별도 프로세스로 띄워 cycle 로그를 모니터링한다.
# v0.3 에서는 --dry-run 으로만 실행하여 kill 행위는 하지 않고 데이터 수집만 한다.
# v0.4 에서 --dry-run 제거 + 실제 kill 활성화 예정.
def _start_qmmm_watchdog(parallel_workers: int) -> Optional[subprocess.Popen]:
    """qmmm_proactive_oom.py 를 dry-run 모드로 시작한다.

    조건: 워커가 2개 이상이고 watchdog 스크립트가 존재할 때만 가동.
    실패 경로(예: pynvml/nvidia-smi 부재)는 모두 무시하며, 어떠한 예외도
    파이프라인 본체에 전파하지 않는다.
    """
    watchdog_script = os.path.join(UPDD_DIR, "utils", "qmmm_proactive_oom.py")
    if not (os.path.exists(watchdog_script) and parallel_workers > 1):
        return None
    log_base_dir = (
        os.path.dirname(_config.current_log_file)
        if _config.current_log_file else LOG_DIR
    )
    watchdog_log = os.path.join(log_base_dir, "qmmm_live.log")
    try:
        cmd = [
            sys.executable, watchdog_script,
            "--log", watchdog_log,
            "--vram-threshold", "90",
            "--interval", "30",
            "--dry-run",  # v0.3: 모니터링 only — kill 비활성화 (v0.4 에서 해제)
        ]
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
            env=_CHILD_ENV_UNBUFFERED,
        )
        print(f"  🛡️  VRAM Watchdog 시작 (PID: {proc.pid}, dry-run, log={watchdog_log})")
        return proc
    except (OSError, FileNotFoundError) as _wd_err:
        print(f"  [!] Watchdog 시작 실패 (무시하고 진행): {_wd_err}")
        return None


def _stop_qmmm_watchdog(proc: Optional[subprocess.Popen]) -> None:
    """Watchdog 프로세스를 안전하게 종료한다 (정상/예외 경로 모두 호출 가능)."""
    if proc is None:
        return
    try:
        proc.terminate()
        proc.wait(timeout=5)
        print(f"  🛡️  VRAM Watchdog 종료")
    except subprocess.TimeoutExpired:
        proc.kill()
        print(f"  🛡️  VRAM Watchdog 강제 종료")
    except OSError:
        pass


def execute_qmmm(snap_dir: str, qmmm_mode: str, ncaa_element: str, af2_out_dir: str,
                 parallel_workers: int, target_out_dir: str, binder_chain: str = "B",
                 target_id: Optional[str] = None) -> str:
    """QM/MM 단계 orchestrator.

    [v0.6.1] ``target_id`` 가 주어지고 ``target_cards/{target_id}.json`` 이 존재하면
    ``run_qmmm.py`` 에 ``--target-id`` 를 전달하여 topology mode (snapshot-invariant
    QM region) 로 실행된다. 카드가 없으면 legacy distance-based mode (fast/plaid/full)
    로 폴백하여 이전 워크플로우와의 호환을 유지한다.
    """
    print_step("Step 11: QM/MM — wB97X-D / PySCF")
    qmmm_out = os.path.join(target_out_dir, "qmmm_results")
    os.makedirs(qmmm_out, exist_ok=True)
    elem_arg = ncaa_element if ncaa_element else "none"

    # [v0.6.1] target_card 존재 여부 확인 → topology mode 활성화 가드.
    # target_id 가 None 이거나 카드 JSON 이 없으면 legacy mode 유지 (backward compat).
    _effective_target_id: Optional[str] = None
    if target_id:
        _card_path = os.path.join(UPDD_DIR, "target_cards", f"{target_id}.json")
        if os.path.exists(_card_path):
            _effective_target_id = target_id
            print(f"  [v0.6.1] Topology mode 활성화 — target_card: {_card_path}")
        else:
            print(f"  [v0.6.1] target_card 미발견 ({_card_path}) — legacy mode ({qmmm_mode}) 유지")

    # ── Phase 1: 스냅샷 존재 디자인만 필터링 ──────────────────
    af2_ranking_path = os.path.join(af2_out_dir, "af2_ranking.csv")
    if not os.path.exists(af2_ranking_path):
        print("  [!] af2_ranking.csv 없음 — QM/MM 건너뜀")
        return qmmm_out

    ranking_data: List[Dict[str, Any]] = []
    with open(af2_ranking_path, 'r', encoding="utf-8") as f:
        for row in csv.DictReader(f):
            ranking_data.append({'id': row['id'], 'plddt': float(row['plddt'])})
    ranking_data.sort(key=lambda x: x['plddt'], reverse=True)

    valid_ids = []
    skipped_ids = []
    for item in ranking_data:
        snaps = glob.glob(os.path.join(snap_dir, f"{item['id']}*.pdb"))
        if snaps:
            valid_ids.append(item['id'])
        else:
            skipped_ids.append(item['id'])

    if skipped_ids:
        print(f"  [Filter] 스냅샷 없는 디자인 {len(skipped_ids)}개 제외")
    if not valid_ids:
        print("  [!] 유효한 스냅샷이 없습니다 — QM/MM 건너뜀")
        return qmmm_out

    # ── Phase 2: Resume — 이미 정상 완료된 디자인 스킵 ────────
    def _is_qmmm_done(design_id: str) -> bool:
        """해당 디자인의 모든 스냅샷 QM/MM이 정상 완료되었는지 확인.

        [v0.5 Bug 5 fix] FAIL tier 디자인은 QM/MM 계산 불필요 (ranking_include=false).
        """
        if _V05_MD_STATUS_AVAILABLE:
            _tier = get_design_tier(target_out_dir, design_id)
            if _tier == "FAIL":
                return True
        snap_pdbs = sorted(glob.glob(os.path.join(snap_dir, f"{design_id}*.pdb")))
        if not snap_pdbs:
            return False
        for snap_pdb in snap_pdbs:
            stem = os.path.splitext(os.path.basename(snap_pdb))[0]
            json_path = os.path.join(qmmm_out, f"{stem}_qmmm_{qmmm_mode}.json")
            if not os.path.exists(json_path):
                return False
            try:
                with open(json_path, "r", encoding="utf-8") as jf:
                    data = json.load(jf)
                # run_qmmm.py 출력 키: energy_total_hartree, energy_qm_hartree
                e_val = data.get("energy_total_hartree", data.get("energy_qm_hartree", 0))
                if not e_val or e_val == 0:
                    return False  # OOM/에러로 생성된 빈 JSON
            except (json.JSONDecodeError, KeyError, TypeError):
                return False
        return True

    already_done = []
    need_calc = []
    for d_id in valid_ids:
        if _is_qmmm_done(d_id):
            already_done.append(d_id)
        else:
            need_calc.append(d_id)

    # [v0.4.1 B] 정확한 total (valid design 수) 로 state 갱신 →
    # dashboard 가 "running 3/8" 같이 진행도를 정확히 표시.
    _ps = _pipeline_state
    if _ps is not None:
        _ps.update_step(
            "qmmm",
            completed=len(already_done),
            detail=f"{len(already_done)}/{len(valid_ids)} designs",
        )

    if already_done:
        print(f"  [Resume] 이미 완료된 디자인 {len(already_done)}개 스킵")
    if not need_calc:
        print(f"  [Resume] 모든 디자인 완료 — QM/MM 재계산 불필요")
    else:
        print(f"  [우선순위] {len(need_calc)}개 디자인 계산 예정 "
              f"(전체 {len(valid_ids)}개 중 {len(already_done)}개 완료)")

        # [v0.3 #6] VRAM Watchdog 자동 시작 (Phase 3 직전, dry-run 모드).
        # Phase 3 종료 후 (정상/예외 경로 모두) finally 블록에서 정리한다.
        watchdog_proc = _start_qmmm_watchdog(parallel_workers)

        # ── QM/MM 워커 정의 (태그 + OOM 감지) ────────────────
        base_cmd = [
            "conda", "run", "-n", ENV_NAMES.get("qmmm", "qmmm"),
            "--no-capture-output", "python", "-u", SCRIPT_PATHS["run_qmmm"],
        ]

        def qmmm_worker(design_id: str) -> dict:
            qmmm_args = [
                "--snapdir",      snap_dir,
                "--outputdir",    qmmm_out,
                "--ncaa_elem",    elem_arg,
                "--mode",         qmmm_mode,
                "--qm_xc",        "wb97xd",
                "--filter",       design_id,
                "--binder_chain", binder_chain,
            ]
            # [v0.6.1] Topology mode wiring — target card 존재 시에만 주입.
            # run_qmmm.py 의 --target-id 가 설정되면 --mode 는 무시되고
            # target_cards/{id}.json 의 QM partitioning 정책이 적용된다.
            if _effective_target_id is not None:
                qmmm_args += ["--target-id", _effective_target_id]
            cmd = base_cmd + qmmm_args

            # 짧은 태그: "design_w4_1_s2" → "w4_1_s2"
            tag = design_id.replace("design_", "").split("_unrelaxed")[0]
            if len(tag) > 10:
                tag = tag[:10]

            try:
                # [v0.3 #2] PYTHONUNBUFFERED=1 — PIPE 통과 시 블록 버퍼링 방지.
                # [v0.3.1] fix(pipeline): write QM/MM subprocess output to qmmm_live.log for watchdog
                # [v0.3.2] fix(pipeline): align qmmm_live.log path with watchdog 
                # DFT cycle 로그가 즉시 [tag] prefix 와 함께 표시되어 watchdog 파싱과
                # 운영자 모니터링이 동시 가능.
                proc = subprocess.Popen(
                    cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, bufsize=1, env=_CHILD_ENV_UNBUFFERED,
                )
                oom_detected = False
                _log_dir = os.path.dirname(_config.current_log_file) if _config.current_log_file else LOG_DIR
                qmmm_log_path = os.path.join(_log_dir, "qmmm_live.log")
                for line in proc.stdout:
                    stripped = line.rstrip()
                    if stripped:
                        print(f"  [{tag}] {stripped}")
                        # watchdog이 파싱할 수 있도록 파일에도 기록
                        with open(qmmm_log_path, "a") as _lf:
                            _lf.write(f"  [{tag}] {stripped}\n")
                            _lf.flush()
                        if "out of memory" in stripped.lower() or "cudaErrorMemoryAllocation" in stripped:
                            oom_detected = True
                proc.wait()

                if oom_detected:
                    return {"id": design_id, "status": "OOM",
                            "msg": f"CUDA OOM (returncode={proc.returncode})"}
                elif proc.returncode != 0:
                    return {"id": design_id, "status": "ERROR",
                            "msg": f"비정상 종료 (returncode={proc.returncode})"}
                else:
                    return {"id": design_id, "status": "SUCCESS", "msg": ""}
            except Exception as e:
                return {"id": design_id, "status": "ERROR", "msg": str(e)}

        # ── Phase 3: Adaptive 병렬 실행 (4→3→2→1) ────────────
        pending_ids = list(need_calc)
        _qmmm_cap = int(os.environ.get("UPDD_QMMM_MAX_WORKERS", 1))
        current_workers = min(_qmmm_cap, parallel_workers, len(pending_ids))
        all_success_ids = list(already_done)
        max_passes = 4

        for pass_num in range(1, max_passes + 1):
            if not pending_ids:
                break

            print(f"\n  🔥 [QM/MM Pass {pass_num}] {len(pending_ids)}개 디자인, "
                  f"{current_workers}개 GPU 워커")

            pass_results = []
            try:
                with concurrent.futures.ThreadPoolExecutor(max_workers=current_workers) as executor:
                    future_map = {
                        executor.submit(qmmm_worker, d_id): d_id
                        for d_id in pending_ids
                    }
                    # [v0.4.1 B] 현재 pass 의 worker 수 / pending 수를 이벤트로 안내.
                    # (Step.workers 를 직접 갱신하려면 update_step API 확장 필요 —
                    #  detail 문자열로 pass 번호를 담고 이벤트 로그에 남기는 것으로 충분.)
                    if _ps is not None:
                        _ps.emit_event(
                            "info",
                            f"🔥 QM/MM Pass {pass_num} — workers={current_workers}, "
                            f"pending={len(pending_ids)}",
                        )
                    for future in concurrent.futures.as_completed(future_map):
                        result = future.result()
                        pass_results.append(result)
                        if result["status"] == "SUCCESS":
                            print(f"  ✅ {result['id']} 완료")
                        elif result["status"] == "OOM":
                            print(f"  ⚠️  {result['id']} OOM — 다음 pass에서 재시도")
                        else:
                            print(f"  ❌ {result['id']} 에러: {result['msg']}")
                        # [v0.4.1 B] 한 design 완료 시점에 dashboard 진행도 갱신.
                        if _ps is not None:
                            _ok_count = (
                                len(already_done)
                                + sum(1 for r in pass_results if r["status"] == "SUCCESS")
                            )
                            _ps.update_step(
                                "qmmm", completed=_ok_count,
                                detail=f"{_ok_count}/{len(valid_ids)} designs "
                                       f"(pass {pass_num})",
                            )
            except KeyboardInterrupt:
                print(f"\n  [!] 🛑 수동 개입! 병렬 작업을 강제 중단합니다.\n")
                break

            succeeded = [r for r in pass_results if r["status"] == "SUCCESS"]
            oom_failed = [r for r in pass_results if r["status"] == "OOM"]

            all_success_ids.extend([r["id"] for r in succeeded])

            if oom_failed:
                pending_ids = [r["id"] for r in oom_failed]
                # [#14 v0.2] OOM 실패 디자인의 기존 JSON 정리 — 단, 유효 결과는 보존.
                # 유효성 기준: converged=True AND interaction_kcal not in (None, 0.0).
                # 무효(0.0/None/converged=False) JSON 만 제거하여 재계산 대상으로 돌린다.
                for r in oom_failed:
                    for stale in glob.glob(os.path.join(qmmm_out, f"{r['id']}*_qmmm_*.json")):
                        try:
                            with open(stale, "r", encoding="utf-8") as _sf:
                                _sdata = json.load(_sf)
                            _ie = _sdata.get("interaction_kcal")
                            _cv = _sdata.get("converged", False)
                            if _cv and _ie not in (None, 0.0):
                                print(f"  [Resume] 유효 JSON 보존: {stale}")
                                continue
                        except (json.JSONDecodeError, IOError, KeyError, TypeError):
                            pass
                        os.remove(stale)
                new_workers = max(1, current_workers - 1)
                print(f"\n  [Adaptive] OOM {len(oom_failed)}개 감지. "
                      f"워커 축소: {current_workers} → {new_workers}")
                current_workers = new_workers
            else:
                pending_ids = []

        if pending_ids:
            print(f"  [!] {len(pending_ids)}개 디자인이 단일 워커에서도 OOM. "
                  f"VRAM 부족으로 건너뜁니다.")

        total = len(valid_ids)
        n_success = len(all_success_ids)
        n_oom = len(pending_ids)
        n_error = total - n_success - n_oom
        print(f"\n  [QM/MM 요약] 총 {total}개 | "
              f"✅ 성공 {n_success} | ⚠️ OOM {n_oom} | ❌ 에러 {n_error}")

        # [v0.3 #6] Phase 3 종료 — Watchdog 정리.
        # KeyboardInterrupt 는 위 inner try 에서 break 으로 처리되어 정상 흐름으로 도달한다.
        # 비정상 예외 시는 부모 프로세스가 종료되며 OS 가 자식(watchdog) 을 회수한다.
        _stop_qmmm_watchdog(watchdog_proc)

    # ── [v55 유지] 개별 JSON → 통합 summary 병합 ─────────────
    all_results = []
    failed = []
    for jf in sorted(glob.glob(os.path.join(qmmm_out, f"*_qmmm_{qmmm_mode}.json"))):
        try:
            with open(jf, "r", encoding="utf-8") as f:
                result = json.load(f)
            all_results.append(result)
        except (json.JSONDecodeError, IOError) as _je:
            failed.append(os.path.basename(jf))
    summary_path = os.path.join(qmmm_out, "qmmm_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump({
            "method":    f"wb97xd/6-31G*",
            "n_calc":    len(all_results),
            "n_failed":  len(failed),
            "failed":    failed,
            "results":   all_results,
        }, f, indent=2)
    converged = [r for r in all_results if r.get("converged")]
    print(f"  [QM/MM 통합] {len(converged)}/{len(all_results)} 수렴, 요약: {summary_path}")

    return qmmm_out

# ==========================================
# [v35 AUDIT] Principle 3/4: Structural validation entrypoint
# ==========================================
# Rationale: `--validate` mode is the autonomous-verification contract declared in
# UPDD.md Principle 3. Prior to this audit, main() referenced `validate_pipeline(...)`
# but no such function was defined — attempting `python UPDD.py --validate` crashed
# with NameError, silently breaking Principle 3. This function performs a read-only
# structural audit of all expected artifacts produced by the pipeline stages and
# emits a per-stage OK/MISSING report. It does NOT re-execute any stage.
def validate_pipeline(base_name, design_dir, mpnn_out_dir, af2_out_dir,
                      ncaa_mut_dir, params_dir, md_dir, snap_dir,
                      qmmm_dir, mmgbsa_dir, final_csv, ncaa_code):
    print_step(f"[Validate] Structural audit of '{base_name}' artifact tree")

    def _stage(label, path, predicate):
        status = "OK" if path and predicate(path) else "MISSING"
        sigil = "✅" if status == "OK" else "⚠️"
        print(f"  {sigil} {label:<30} {path or '(not applicable)'} [{status}]")
        return status == "OK"

    def _has_any(pattern):
        return len(glob.glob(pattern, recursive=True)) > 0

    results = {
        "designs":      _stage("Step 1 RFdiffusion",    design_dir,  lambda p: _has_any(os.path.join(p, "**", "*.pdb"))),
        "mpnn":         _stage("Step 3 ProteinMPNN",    mpnn_out_dir, lambda p: _has_any(os.path.join(p, "**", "*.fa*"))),
        "af2":          _stage("Step 5 AlphaFold2",     af2_out_dir,  lambda p: _has_any(os.path.join(p, "**", "*.pdb"))),
        "ncaa_mutants": _stage("Step 7 ncAA Mutation",  ncaa_mut_dir if ncaa_code else "", lambda p: os.path.isdir(p) and _has_any(os.path.join(p, "*.pdb"))),
        "params":       _stage("Step 8 Parameterization", params_dir if ncaa_code else "", lambda p: os.path.isdir(p) and _has_any(os.path.join(p, "*_params_manifest.json"))),
        # DCD trajectories live on the external drive when --hdd_path was
        # supplied (run_restrained_md streams DCDReporter there directly).
        # Validate MD completion via _final.pdb / _md.log which always stay
        # on SSD md_dir regardless of the HDD backup toggle.
        "md":           _stage("Step 9 Restrained MD",  md_dir,      lambda p: os.path.isdir(p) and (_has_any(os.path.join(p, "*_final.pdb")) or _has_any(os.path.join(p, "*_md.log")) or _has_any(os.path.join(p, "*.dcd")))),
        "snapshots":    _stage("Step 10 Snapshots",     snap_dir,    lambda p: _has_any(os.path.join(p, "**", "*.pdb"))),
        "qmmm":         _stage("Step 11 QM/MM",         qmmm_dir,    lambda p: _has_any(os.path.join(p, "**", "qmmm_summary.json"))),
        "mmgbsa":       _stage("Step 12 MM-GBSA",       mmgbsa_dir,  lambda p: _has_any(os.path.join(p, "**", "mmgbsa_summary.json"))),
        "ranking":      _stage("Step 13 Final Ranking", final_csv,   lambda p: os.path.isfile(p)),
    }

    n_ok = sum(1 for v in results.values() if v)
    n_total = len(results)
    print(f"\n  [Validate] {n_ok}/{n_total} stages produced expected artifacts.")
    if n_ok < n_total:
        print(f"  [Validate] Missing stages: {[k for k, v in results.items() if not v]}")
    print(f"  [Validate] Structural audit complete — no re-execution performed.")
    return results


# ==========================================
# 메인 실행부
# ==========================================
def main() -> None:
    # [v0.4] argparse — v0.3.3 호환 + dashboard/observer 플래그
    _parser = argparse.ArgumentParser(add_help=False,
        description="UPDD Pipeline v0.4 (btop-style dashboard)")
    _parser.add_argument("--validate", action="store_true",
        help="실행 없이 결과 파일 검증만 수행")
    _parser.add_argument("--no-dashboard", action="store_true",
        help="Dashboard UI 비활성화 (v0.3.3 레거시 print 모드)")
    _parser.add_argument("--dashboard-only", action="store_true",
        help="Pipeline 실행 없이 기존 state.json 을 관찰 (SSH 재접속/모니터링 용)")
    _parser.add_argument("--log-dir", type=str, default=None,
        help="--dashboard-only 모드에서 관찰할 로그 디렉토리 (생략 시 최신 자동 탐색)")
    _args, _ = _parser.parse_known_args()
    VALIDATE_ONLY = _args.validate
    NO_DASHBOARD = _args.no_dashboard or (not _V04_DASHBOARD_AVAILABLE)
    DASHBOARD_ONLY = _args.dashboard_only

    # ── [v0.4] --dashboard-only: 파이프라인 우회, 관찰 전용 ──
    if DASHBOARD_ONLY:
        if not _V04_DASHBOARD_AVAILABLE:
            sys.exit(f"[!] --dashboard-only 사용 불가 (v0.4 모듈 로드 실패: {_V04_IMPORT_ERROR})")
        target_log = _args.log_dir
        if not target_log:
            _cands = [c for c in glob.glob(os.path.join(LOG_DIR, "*"))
                      if os.path.isdir(c)
                      and os.path.exists(os.path.join(c, "state.json"))]
            if not _cands:
                sys.exit(f"[!] state.json 이 있는 로그 디렉토리 없음. --log-dir 지정 필요.\n    검색 대상: {LOG_DIR}")
            target_log = max(_cands, key=os.path.getmtime)
            print(f"[dashboard-only] 자동 탐색 결과: {target_log}")
        _run_dashboard_standalone(log_dir=target_log)
        return

    setup_logger("session_start")
    _builtins.input = logged_input

    _header_label = "UPDD Pipeline v0.4 (btop dashboard)" if not NO_DASHBOARD else "UPDD Pipeline v0.4 (legacy mode)"
    print_step(_header_label)

    if not _config.af2_bin:
        print(" [!] colabfold_batch 바이너리를 찾을 수 없습니다.")
        sys.exit(1)

    os.makedirs(TARGET_DIR, exist_ok=True)
    pdb_files = sorted(glob.glob(os.path.join(TARGET_DIR, "*.pdb")))
    if not pdb_files:
        sys.exit(f"[!] {TARGET_DIR} 에 PDB 파일이 없습니다.")

    choice = arrow_menu(f"target/ 디렉토리 PDB 타겟 선택", [os.path.basename(p) for p in pdb_files])
    target_pdb = pdb_files[choice]

    _target_name = os.path.splitext(os.path.basename(target_pdb))[0]
    setup_logger(_target_name, append_from=_config.current_log_file)
    base_name = os.path.basename(target_pdb).replace(".pdb", "")
    outputs_base = os.path.join(UPDD_DIR, "outputs")
    os.makedirs(outputs_base, exist_ok=True)

    status = {}
    resume_mode = False
    target_out_dir = ""

    possible_folders = sorted(glob.glob(os.path.join(outputs_base, f"{base_name}*")))
    valid_projects = [(f, os.path.join(f, "updd_status.json")) for f in possible_folders if os.path.exists(os.path.join(f, "updd_status.json"))]

    if valid_projects:
        options = [f"{os.path.basename(f)} (상태 파일 존재)" for f, _ in valid_projects]
        options.append("✨ 새로운 프로젝트 생성")
        choice_idx = arrow_menu(f"'{base_name}' 타겟에 대한 이전 프로젝트들이 발견되었습니다.", options)
        
        if choice_idx < len(valid_projects):
            target_out_dir, status_file = valid_projects[choice_idx]
            with open(status_file, "r", encoding="utf-8") as f:
                inputs = json.load(f)
            resume_mode = True
            print(f"\n  [▶️ Resume] '{os.path.basename(target_out_dir)}' 프로젝트를 이어서 진행합니다!\n")
            
            ncaa_def = resolve_ncaa_definition(inputs.get("ncaa_code", "none"))
            topology_mode = inputs.get("topology", TOPOLOGY_TYPES[1])
            binder_chain = inputs.get("binder_chain", "B")
            target_chain = inputs.get("target_chain", "A")
            ncaa_positions = inputs.get("ncaa_positions", "")
            qmmm_mode = inputs.get("qmmm_mode", "fast")
            parallel_workers = int(inputs.get("parallel_workers", 4))
            admet_mode = inputs.get("admet_mode", "none")
            msa_choice = inputs.get("msa_choice", "local")
            USE_LOCAL_MSA = False if msa_choice == "server" else True
            msa_mode = "single_sequence" if USE_LOCAL_MSA else "MMseqs2"
            dcd_backup_path = inputs.get("dcd_backup_path", "")
            
            work_pdb = os.path.join(TARGET_DIR, f"{base_name}_clean.pdb")
            if inputs.get("grafted"):
                work_pdb = os.path.join(TARGET_DIR, "grafted_target.pdb")

    if not resume_mode:
        print("\n  [New] 새로운 프로젝트를 시작합니다.")
        inputs = gather_all_inputs(base_name, target_pdb)
        work_pdb = inputs["preprocess_opts"]["target_pdb"]

        # [R-17 Stage 4 orchestrator wiring] Step 3 Target Preprocessing — executes
        # utils/preprocess_target.py when the user opted in. Previously this
        # function existed (``run_target_preprocessing``) but was never invoked
        # from main(); preprocess was effectively a no-op. Wiring it here with
        # ``target_id=base_name`` activates G1 (cofactor keep-list + presence
        # validation) when target_cards/<base_name>.json is present. For cofactor
        # -absent targets G1 is a transparent no-op.
        _pre_opts = inputs.get("preprocess_opts", {})
        if _pre_opts.get("do_preprocess") and not VALIDATE_ONLY:
            print_step("Step 3: Target Preprocessing (HETATM filter + PDBFixer)")
            work_pdb = run_target_preprocessing(
                _pre_opts.get("target_pdb", target_pdb),
                _pre_opts.get("keep_mode", "auto"),
                _pre_opts.get("chain_in", ""),
                _pre_opts.get("keep_res", ""),
                target_id=base_name,
            )

        topology_mode = inputs["topology"]
        ncaa_def = resolve_ncaa_definition(inputs["ncaa_code"])
        binder_chain = inputs["binder_chain"]
        target_chain = inputs["target_chain"]
        ncaa_positions = inputs["ncaa_positions"]
        qmmm_mode = inputs["qmmm_mode"]
        admet_mode = inputs["admet_mode"]
        parallel_workers = inputs["parallel_workers"]
        msa_choice = inputs["msa_choice"]
        USE_LOCAL_MSA = False if msa_choice == "server" else True
        msa_mode = "single_sequence" if USE_LOCAL_MSA else "MMseqs2"
        topo_name = topology_mode.get("abbr", "linear")
        folder_name = f"{base_name}_{topo_name}_{inputs['ncaa_code']}_{inputs['n_designs']}_{inputs['binder_len']}" if ncaa_def else f"{base_name}_{topo_name}_{inputs['n_designs']}_{inputs['binder_len']}"
        target_out_dir = os.path.join(outputs_base, folder_name)
        os.makedirs(target_out_dir, exist_ok=True)

        # ── [RESTORED] DCD 백업 경로 구성 (이전 버전에서 누락됨) ──
        hdd_name = inputs.get("hdd_name", "")
        dcd_backup_path = ""
        if hdd_name:
            if hdd_name.startswith("/"):
                dcd_backup_path = os.path.join(hdd_name, "UPDD_DCD_Backup", folder_name)
            else:
                dcd_backup_path = os.path.join("/media/san", hdd_name, "UPDD_DCD_Backup", folder_name)
            os.makedirs(dcd_backup_path, exist_ok=True)
            print(f"  [DCD Backup] 궤적 저장 경로: {dcd_backup_path}")

        inputs["dcd_backup_path"] = dcd_backup_path

        status_file = os.path.join(target_out_dir, "updd_status.json")
        with open(status_file, "w", encoding="utf-8") as f:
            json.dump(inputs, f, indent=4)

    design_dir = os.path.join(target_out_dir, "designs")
    os.makedirs(design_dir, exist_ok=True)

    # [v3 1-4] Ligand Grafting 단계 연결: gather_all_inputs 가 grafting_needed=True 이고
    # ref_pdb 경로가 제공된 경우, 새로운 프로젝트 진입 시 work_pdb 를 grafted_target.pdb 로
    # 교체한다. resume_mode 에서는 이전 실행의 결정(grafted=True/False) 을 그대로 따른다.
    if not VALIDATE_ONLY and not resume_mode and inputs.get("grafting_needed") and inputs.get("ref_pdb"):
        print_step("Step 2: Ligand Grafting (결합 부위 모방)")
        work_pdb = run_ligand_grafting(work_pdb, inputs["ref_pdb"], target_chain=target_chain)
        # status_file 에 grafted 플래그를 기록하여 resume_mode 가 동일 경로를 재현하도록 함
        try:
            with open(status_file, "r", encoding="utf-8") as f:
                _persisted = json.load(f)
            _persisted["grafted"] = True
            with open(status_file, "w", encoding="utf-8") as f:
                json.dump(_persisted, f, indent=4)
        except (IOError, OSError, json.JSONDecodeError) as _persist_err:
            log(f"[grafting] status_file grafted 플래그 갱신 실패(무시): {_persist_err}", level="warning")

    # ── [v0.4] 완료 step 이름 → PipelineState step 이름 매핑 ──
    # updd_status.json 의 legacy key → updd_state 의 11단계 key 매핑.
    # 일부 save_step 호출 하나가 state 2개를 완료시키는 경우도 있음 (md_snapshots).
    _STATE_STEP_MAP: Dict[str, List[str]] = {
        "rfdiffusion":   ["rfdiff"],
        "proteinmpnn":   ["proteinmpnn"],
        "alphafold":     ["af2"],
        "ncaa_mutation": ["ncaa_mutation"],
        "md_snapshots":  ["md", "snapshots"],
        "qmmm":          ["qmmm"],
        "mmgbsa":        ["mmgbsa"],
        "final_ranking": ["rank"],
    }

    # ── [v0.4] Dashboard 시작 (사용자 입력 완료 후, pipeline 실행 직전) ──
    # rich.Live(screen=True) 가 alternate screen 을 점유하므로, arrow_menu /
    # gather_all_inputs 등 대화형 입력이 끝난 뒤에만 dashboard 를 띄운다.
    # stdout/stderr 은 <log_dir>/pipeline.log 로 리다이렉트 → print() 가 dashboard
    # 아래 화면을 오염시키지 않고 로그 파일에 구조화되어 저장된다.
    state = None
    dashboard_thread = None
    dashboard_stop = None
    pipeline_log_fh = None
    _orig_stdout = sys.stdout
    _orig_stderr = sys.stderr
    use_dashboard = (not NO_DASHBOARD) and (not VALIDATE_ONLY) and _V04_DASHBOARD_AVAILABLE
    # [FIX 2026-04-20] PipelineState(state.json) 는 `--no-dashboard` 여부와
    # 무관하게 항상 기록한다. --dashboard-only 옵저버나 외부 모니터가 현재
    # 파이프라인의 진행 상태를 볼 수 있도록 한다. --no-dashboard 는 오직
    # "live 화면 렌더링"만 끈다.
    enable_state = (not VALIDATE_ONLY) and _V04_DASHBOARD_AVAILABLE

    log_dir_for_state = os.path.join(LOG_DIR, _target_name)

    if enable_state:
        try:
            os.makedirs(log_dir_for_state, exist_ok=True)
            state = PipelineState(
                log_dir=log_dir_for_state,
                project_name=os.path.basename(target_out_dir),
            )
            # [v0.4.1] Resume 초기 상태: "completed_steps 플래그" AND "실제 산출물 존재"
            # 두 조건 모두 만족해야 SKIPPED 로 표시. 한쪽이라도 false 면 PENDING 유지하여
            # 아래 step wrapper 가 start_step → 실행 → complete_step 으로 자연스럽게 흐름.
            _already_done = inputs.get("completed_steps", []) if isinstance(inputs, dict) else []
            for _cs in _already_done:
                for _sn in _STATE_STEP_MAP.get(_cs, []):
                    if check_outputs_exist(_sn, inputs, target_out_dir):
                        state.complete_step(_sn, skipped=True)
                    # 산출물 없으면 SKIPPED 마킹하지 않음 → wrapper 가 재실행 분기 진입.
            # preprocess 는 사용자 입력 단계에서 이미 완료되었다 (gather_all_inputs).
            state.complete_step("preprocess", skipped=True)
            state.emit_event("info", f"📦 Project: {os.path.basename(target_out_dir)}")

            # [v0.4.1] 모듈 글로벌 state 참조 공개 → run_rfdiffusion/execute_qmmm
            # 같은 장시간 실행 함수가 진행도를 업데이트할 때 읽는다.
            global _pipeline_state
            _pipeline_state = state
        except (OSError, RuntimeError, ImportError) as _state_err:
            print(f"[!] state.json 초기화 실패 (비-dashboard 모드에서도 필요): "
                  f"{_state_err}", file=_orig_stdout)
            state = None

    if use_dashboard and state is not None:
        try:
            # stdout 리다이렉트 "전에" dashboard 용 Console 을 캡처한다.
            # (Console 이 sys.stdout 을 나중에 읽으면 pipeline.log 로 빨려 들어감.)
            from rich.console import Console as _RichConsole
            _dashboard_console = _RichConsole(file=sys.__stdout__, force_terminal=True)

            pipeline_log_path = os.path.join(log_dir_for_state, "pipeline.log")
            pipeline_log_fh = open(pipeline_log_path, "a", encoding="utf-8", buffering=1)
            pipeline_log_fh.write(f"\n{'='*60}\n")
            pipeline_log_fh.write(f"UPDD v0.4 started at {_datetime.datetime.now().isoformat()}\n")
            pipeline_log_fh.write(f"Project: {os.path.basename(target_out_dir)}\n")
            pipeline_log_fh.write(f"{'='*60}\n")
            sys.stdout = pipeline_log_fh
            sys.stderr = pipeline_log_fh

            dashboard_thread, dashboard_stop = run_dashboard_in_thread(
                log_dir_for_state, console=_dashboard_console,
            )

            # atexit: 예외/정상 종료 모두에서 stdout 복원 보장 (finally 대체).
            import atexit as _atexit
            def _v04_final_cleanup():
                try:
                    if dashboard_stop is not None:
                        dashboard_stop.set()
                    if dashboard_thread is not None:
                        dashboard_thread.join(timeout=3)
                except (OSError, RuntimeError):
                    pass
                try:
                    sys.stdout = _orig_stdout
                    sys.stderr = _orig_stderr
                except (OSError, ValueError):
                    pass
                try:
                    if pipeline_log_fh is not None and not pipeline_log_fh.closed:
                        pipeline_log_fh.close()
                except (OSError, ValueError):
                    pass
            _atexit.register(_v04_final_cleanup)

            # signal handler: Ctrl+C → dashboard 정리 후 exit 130.
            import signal as _signal
            def _on_interrupt(signum=None, frame=None):
                if state is not None:
                    state.emit_event("warn", f"⚠️  Signal {signum} — shutting down")
                _v04_final_cleanup()
                sys.exit(130)
            try:
                _signal.signal(_signal.SIGINT, _on_interrupt)
                _signal.signal(_signal.SIGTERM, _on_interrupt)
            except (ValueError, OSError):
                # non-main thread 등에서는 signal 설정 불가 — 무시하고 atexit 에 의존.
                pass
        except (OSError, RuntimeError, ImportError) as _dash_err:
            # 실패 시 v0.3.3 legacy 모드로 안전 폴백.
            print(f"[!] Dashboard 시작 실패, legacy 모드로 전환: {_dash_err}",
                  file=_orig_stdout)
            use_dashboard = False
            state = None
            dashboard_thread = None
            dashboard_stop = None
            try:
                if pipeline_log_fh is not None and not pipeline_log_fh.closed:
                    pipeline_log_fh.close()
            except (OSError, ValueError):
                pass
            sys.stdout = _orig_stdout
            sys.stderr = _orig_stderr

    # ── [v0.4.1] save_step 은 JSON 기록 전용 (state 관리는 wrapper 에서) ──
    def save_step(step_name):
        """완료된 단계를 updd_status.json 에 기록한다 (resume skip 판단용)."""
        try:
            with open(status_file, "r", encoding="utf-8") as f:
                _st = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            _st = {}
        if step_name not in _st.setdefault("completed_steps", []):
            _st["completed_steps"].append(step_name)
        with open(status_file, "w", encoding="utf-8") as f:
            json.dump(_st, f, indent=4)

    # ── [v0.4.1] step wrapper: completed_steps ∩ outputs_exist → skipped ──
    def _step_really_done(legacy_name: str) -> bool:
        """completed_steps 플래그와 실제 산출물 둘 다 있을 때만 True."""
        if not isinstance(inputs, dict):
            return False
        if legacy_name not in inputs.get("completed_steps", []):
            return False
        # legacy key 에 매핑된 state key 중 하나라도 산출물이 있어야 완료로 간주.
        state_names = _STATE_STEP_MAP.get(legacy_name, [legacy_name])
        return any(check_outputs_exist(sn, inputs, target_out_dir) for sn in state_names)

    def _mark_skipped(legacy_name: str) -> None:
        if state is None:
            return
        for _sn in _STATE_STEP_MAP.get(legacy_name, [legacy_name]):
            # 이미 resume 초기화에서 SKIPPED 로 표시됐을 수 있음. complete_step 은
            # 멱등이므로 재호출해도 안전 (같은 SKIPPED 유지).
            state.complete_step(_sn, skipped=True)

    def _mark_running(legacy_name: str, total: int = 0, workers: int = 1,
                      detail: str = "") -> None:
        if state is None:
            return
        for _sn in _STATE_STEP_MAP.get(legacy_name, [legacy_name]):
            state.start_step(_sn, total=total, workers=workers, detail=detail)

    def _mark_complete(legacy_name: str) -> None:
        if state is None:
            return
        for _sn in _STATE_STEP_MAP.get(legacy_name, [legacy_name]):
            state.complete_step(_sn)

    def _mark_failed(legacy_name: str, err: str) -> None:
        if state is None:
            return
        for _sn in _STATE_STEP_MAP.get(legacy_name, [legacy_name]):
            state.fail_step(_sn, err)

    # ── [v0.4] step-independent worker 수 resolve (legacy parallel_workers 유지) ──
    # backward-compat: parallel_workers 가 각 run_* 에 그대로 전달되어 기존 동작 보존.
    # 추가로 rfdiff / qmmm / af2 / mmgbsa / md 에 대해 step 별 worker 를 resolve 한다.
    rfdiff_workers = get_workers("rfdiff", inputs)
    qmmm_workers = get_workers("qmmm", inputs)
    mmgbsa_workers = get_workers("mmgbsa", inputs)
    af2_workers = get_workers("af2", inputs)

    # ── [코어 파 파이프라인 가동] ──
    if not VALIDATE_ONLY:
        safe_contig = get_continuous_contigs(work_pdb, target_chain)
        final_contig = "'" + safe_contig + "/0 " + inputs.get("binder_len", "10-40") + "'"
        # Step: RFdiffusion
        if _step_really_done("rfdiffusion"):
            _mark_skipped("rfdiffusion")
        else:
            # [v0.4.1 B] total = worker 수 (run_rfdiffusion 내부에서 worker 완료 때마다 update)
            _mark_running("rfdiffusion",
                          total=max(1, min(rfdiff_workers, int(inputs.get("n_designs", 5)))),
                          workers=rfdiff_workers)
            try:
                run_rfdiffusion(work_pdb, topology_mode, inputs.get("hotspots", ""), final_contig, design_dir, int(inputs.get("n_designs", 5)), parallel_workers, target_chain)
                _mark_complete("rfdiffusion")
            except Exception as _e:
                _mark_failed("rfdiffusion", str(_e))
                raise
        save_step("rfdiffusion")

    # Step: ProteinMPNN
    if _step_really_done("proteinmpnn"):
        _mark_skipped("proteinmpnn")
        mpnn_out_dir, jsonl_path = run_proteinmpnn(design_dir, target_chain, topology_mode, valid_only=VALIDATE_ONLY)
    else:
        _mark_running("proteinmpnn", total=0, workers=1)
        try:
            mpnn_out_dir, jsonl_path = run_proteinmpnn(design_dir, target_chain, topology_mode, valid_only=VALIDATE_ONLY)
            _mark_complete("proteinmpnn")
        except Exception as _e:
            _mark_failed("proteinmpnn", str(_e))
            raise
    save_step("proteinmpnn")

    # Step: AF2 (alphafold) — [v5 해법 C] topology_mode/work_pdb/binder_chain 전달
    if _step_really_done("alphafold"):
        _mark_skipped("alphafold")
        af2_out_dir = run_alphafold(
            design_dir, mpnn_out_dir, jsonl_path, target_chain, msa_mode, USE_LOCAL_MSA,
            valid_only=VALIDATE_ONLY,
            topology_mode=topology_mode,
            work_pdb=work_pdb,
            binder_chain=binder_chain,
        )
    else:
        _mark_running("alphafold", total=0, workers=af2_workers,
                      detail="ColabFold / CyclicFold")
        try:
            af2_out_dir = run_alphafold(
                design_dir, mpnn_out_dir, jsonl_path, target_chain, msa_mode, USE_LOCAL_MSA,
                valid_only=VALIDATE_ONLY,
                topology_mode=topology_mode,
                work_pdb=work_pdb,
                binder_chain=binder_chain,
            )
            _mark_complete("alphafold")
        except Exception as _e:
            _mark_failed("alphafold", str(_e))
            raise
    save_step("alphafold")

    if VALIDATE_ONLY:
        validate_pipeline(base_name, design_dir, mpnn_out_dir, af2_out_dir, target_out_dir + "/ncaa_mutants" if ncaa_def else "", target_out_dir + "/params", target_out_dir + "/mdresult", target_out_dir + "/snapshots", target_out_dir + "/qmmm_results", target_out_dir + "/mmgbsa_results", os.path.join(target_out_dir, "final_ranking.csv"), ncaa_def.code if ncaa_def else None)
        return

    print_step("Step 6.5: Pre-Filter (Rank, ipTM, RMSD)")
    _iptm_cutoff = float(inputs.get("iptm_cutoff", 0.3))
    # RMSD threshold: user-configurable via gather_all_inputs (default 2.5 Å
    # per AF2-multimer noise envelope, Bryant 2022 Nat Commun IQR 1.5-3.2).
    # Missing or legacy inputs fall back to the G2 module default
    # (DEFAULT_RMSD_THRESHOLD_A = 2.5) inside _run_pre_filter.
    _rmsd_cutoff_in = inputs.get("rmsd_cutoff")
    try:
        _rmsd_cutoff = float(_rmsd_cutoff_in) if _rmsd_cutoff_in is not None else None
    except (TypeError, ValueError):
        _rmsd_cutoff = None
    _run_pre_filter(
        af2_out_dir,
        min_iptm=_iptm_cutoff,
        crystal_pdb=work_pdb,
        target_chain=target_chain,
        rmsd_threshold_a=_rmsd_cutoff,
    )

    # [R-17 G2] Post-AF2 cofactor reinsertion — runs before Step 7 so
    # ncAA mutation and MD setup see HETATM-augmented rank_001 PDBs. For
    # cofactor-absent targets (target_card missing or cofactor_residues=[])
    # this is a logged no-op. For cofactor-present targets (e.g. 6WGN), a
    # CofactorReinsertionError halts the pipeline (fail-fast per R-17).
    print_step("Step 6.6: R-17 G2 Cofactor Reinsertion (pre-ncAA/MD)")
    _reinsert_cofactors_into_af2_outputs(
        af2_out_dir, crystal_pdb=work_pdb,
        target_id=base_name, target_chain=target_chain,
    )

    # Step: ncAA Mutation
    print_step("Step 7: ncAA Configuration & Mutation")
    if _step_really_done("ncaa_mutation"):
        _mark_skipped("ncaa_mutation")
        struct_dir = _run_ncaa_mutation(af2_out_dir, ncaa_def.code if ncaa_def else "none", ncaa_def.label if ncaa_def else "none", ncaa_positions, target_out_dir, binder_chain, target_chain)
    else:
        _mark_running("ncaa_mutation", total=0, workers=1)
        try:
            struct_dir = _run_ncaa_mutation(af2_out_dir, ncaa_def.code if ncaa_def else "none", ncaa_def.label if ncaa_def else "none", ncaa_positions, target_out_dir, binder_chain, target_chain)
            _mark_complete("ncaa_mutation")
        except Exception as _e:
            _mark_failed("ncaa_mutation", str(_e))
            raise
    save_step("ncaa_mutation")

    # Step: Parameterize (save_step 레거시 엔트리 없음 — state 훅만 직접 관리)
    # [R-17 G3] _run_parameterization accepts target_id so cofactor FF
    # parameters (frcmod/mol2) are verified / regenerated for cofactor-present
    # targets before MD setup consumes them in Step 9. No-op when target_card
    # is absent or declares no required cofactors.
    print_step("Step 8: ncAA Parameterization")
    param_dir = os.path.join(target_out_dir, "params")
    _param_really_done = check_outputs_exist("parameterize", inputs, target_out_dir)
    if _param_really_done:
        _mark_skipped("parameterize")
        params_manifest = _run_parameterization(
            struct_dir, param_dir, ncaa_def, target_id=base_name,
            crystal_pdb=work_pdb,
        )
    else:
        _mark_running("parameterize", total=0, workers=1)
        try:
            params_manifest = _run_parameterization(
                struct_dir, param_dir, ncaa_def, target_id=base_name,
                crystal_pdb=work_pdb,
            )
            _mark_complete("parameterize")
        except Exception as _e:
            _mark_failed("parameterize", str(_e))
            raise

    _run_admet(struct_dir, admet_mode, target_out_dir, status if not resume_mode else inputs, status_file)

    # Step: MD + Snapshots (save_step 한 번에 state 2 개 complete 처리)
    # [R-17 G5] target_id=base_name is threaded into _run_md_and_snapshots so
    # the Step 10 extract_snapshots.py subprocess receives --target_id and
    # enforces the cofactor preservation gate during snapshot extraction
    # (strip_solvent exclusion + per-snapshot presence validation).
    print_step("Step 9: Restrained MD")
    canonical_topo = topology_mode.get("canonical", topology_mode.get("abbr", "linear"))
    if _step_really_done("md_snapshots"):
        _mark_skipped("md_snapshots")
        md_out_dir, snap_dir = _run_md_and_snapshots(struct_dir, params_manifest, target_out_dir, ncaa_def.label if ncaa_def else "none", ncaa_def.code if ncaa_def else "", canonical_topo, binder_chain, graph_policy="strict", dcd_backup_path=dcd_backup_path, target_id=base_name)
    else:
        _mark_running("md_snapshots", total=0, workers=get_workers("md", inputs),
                      detail="OpenMM restrained MD")
        try:
            md_out_dir, snap_dir = _run_md_and_snapshots(struct_dir, params_manifest, target_out_dir, ncaa_def.label if ncaa_def else "none", ncaa_def.code if ncaa_def else "", canonical_topo, binder_chain, graph_policy="strict", dcd_backup_path=dcd_backup_path, target_id=base_name)
            _mark_complete("md_snapshots")
        except Exception as _e:
            _mark_failed("md_snapshots", str(_e))
            raise
    save_step("md_snapshots")

    # Step: QM/MM — [v0.4.1 B] total = snap_dir 의 PDB 개수
    _snap_count = 0
    if os.path.isdir(snap_dir):
        _snap_count = len([f for f in os.listdir(snap_dir) if f.endswith(".pdb")])
    # [v0.6.1] target_id 배선 — PDB stem (base_name) 이 target_cards/{id}.json 과 매칭.
    # 카드가 있으면 execute_qmmm 내부에서 topology mode 로 전환, 없으면 legacy mode 유지.
    _qmmm_target_id = inputs.get("target_id") or inputs.get("pdb_id") or base_name
    if _step_really_done("qmmm"):
        _mark_skipped("qmmm")
        qmmm_out = execute_qmmm(snap_dir, qmmm_mode, ncaa_def.element if ncaa_def else "", af2_out_dir, parallel_workers, target_out_dir, binder_chain, target_id=_qmmm_target_id)
    else:
        _mark_running("qmmm", total=_snap_count, workers=qmmm_workers,
                      detail="DFT QM/MM (gpu4pyscf)")
        try:
            qmmm_out = execute_qmmm(snap_dir, qmmm_mode, ncaa_def.element if ncaa_def else "", af2_out_dir, parallel_workers, target_out_dir, binder_chain, target_id=_qmmm_target_id)
            _mark_complete("qmmm")
        except Exception as _e:
            _mark_failed("qmmm", str(_e))
            raise
    save_step("qmmm")

    # ── [RESTORED] Step 12: MM-GBSA 결합 자유에너지 계산 ──
    # [v0.4.1 B] polling thread 로 5s 마다 mmgbsa_results/*.json 개수 → state.update_step.
    print_step("Step 12: MM-GBSA Binding Free Energy")
    mmgbsa_out = os.path.join(target_out_dir, "mmgbsa_results")
    os.makedirs(mmgbsa_out, exist_ok=True)
    mmgbsa_args = [
        # [#18 v0.2] MD 최종 PDB 1개 → MD 스냅샷 다수로 변경 (ensemble ΔG).
        "--md_dir",          snap_dir,
        "--outputdir",       mmgbsa_out,
        "--ncaa_elem",       ncaa_def.element if ncaa_def else "none",
        "--receptor_chain",  target_chain,
        "--binder_chain",    binder_chain,
    ]
    # [R-17 G6] forward target_id so run_mmgbsa registers GNP / ATP cofactor
    # GAFF templates on the ForceField before createSystem.
    if base_name:
        _card_path = os.path.join(UPDD_DIR, "target_cards", f"{base_name}.json")
        if os.path.exists(_card_path):
            mmgbsa_args.extend(["--target_id", base_name])
    if _step_really_done("mmgbsa"):
        _mark_skipped("mmgbsa")
        run_conda_command(SCRIPT_PATHS["run_mmgbsa"], mmgbsa_args, "MM-GBSA ΔG 계산", env_key="md")
    else:
        _mark_running("mmgbsa", total=_snap_count, workers=mmgbsa_workers,
                      detail="AMBER/OpenMM ensemble ΔG")
        # [v0.4.1] polling thread — blocking subprocess 인 run_mmgbsa 의 진행도를
        # mmgbsa_results/*.json 개수로 추정 (과학 로직 건드리지 않음).
        import threading as _thr
        _mmgbsa_stop = _thr.Event()
        def _mmgbsa_poll():
            while not _mmgbsa_stop.wait(5.0):
                try:
                    _n = sum(1 for f in os.listdir(mmgbsa_out)
                             if f.endswith(".json") or f.endswith(".csv"))
                except OSError:
                    _n = 0
                if state is not None:
                    state.update_step("mmgbsa", completed=_n)
        _mmgbsa_poller = _thr.Thread(target=_mmgbsa_poll, daemon=True)
        _mmgbsa_poller.start()
        try:
            run_conda_command(SCRIPT_PATHS["run_mmgbsa"], mmgbsa_args, "MM-GBSA ΔG 계산", env_key="md")
            _mark_complete("mmgbsa")
        except Exception as _e:
            _mark_failed("mmgbsa", str(_e))
            raise
        finally:
            _mmgbsa_stop.set()
            _mmgbsa_poller.join(timeout=2)
    save_step("mmgbsa")

    # ── [RESTORED] Step 13: Final Ranking (pLDDT + QM/MM + ΔG 종합) ──
    print_step("Step 13: Final Ranking (pLDDT + QM/MM + ΔG)")
    final_csv = os.path.join(target_out_dir, "final_ranking.csv")
    rank_args = [
        "--af2_dir",    af2_out_dir,
        "--qmmm_dir",   qmmm_out,
        "--mmgbsa_dir",  mmgbsa_out,
        "--topology",   topology_mode.get("abbr", "linear"),
        "--outputcsv",  final_csv,
    ]
    if _step_really_done("final_ranking"):
        _mark_skipped("final_ranking")
        run_conda_command(SCRIPT_PATHS["rank_qmmm"], rank_args, "최종 후보 Ranking", env_key="qmmm")
    else:
        _mark_running("final_ranking", total=0, workers=1)
        try:
            run_conda_command(SCRIPT_PATHS["rank_qmmm"], rank_args, "최종 후보 Ranking", env_key="qmmm")
            _mark_complete("final_ranking")
        except Exception as _e:
            _mark_failed("final_ranking", str(_e))
            raise
    save_step("final_ranking")

    print_step("UPDD Pipeline 완료")
    topo_display = topology_mode.get("name", "Unknown")
    ncaa_display = ncaa_def.label if ncaa_def else "Wild-Type (없음)"
    print(f"  Topology : {topo_display}")
    print(f"  ncAA     : {ncaa_display}")
    print(f"  최종 산출: {final_csv}")
    print(f"  QM/MM    : {qmmm_out}")
    print(f"  MM-GBSA  : {mmgbsa_out}")
    print(f"  MD 결과  : {os.path.join(target_out_dir, 'mdresult')}")
    if dcd_backup_path:
        print(f"  DCD 백업 : {dcd_backup_path}")
    print()

    # ── [v0.4] 정상 종료: dashboard 마지막 상태 반영 후 깔끔히 정리 ──
    if use_dashboard and state is not None:
        try:
            state.emit_event("pipeline_complete", "🎉 Pipeline finished")
            time.sleep(2)  # 마지막 상태가 dashboard 에 렌더되도록
            if dashboard_stop is not None:
                dashboard_stop.set()
            if dashboard_thread is not None:
                dashboard_thread.join(timeout=5)
        except (OSError, RuntimeError):
            pass
        finally:
            sys.stdout = _orig_stdout
            sys.stderr = _orig_stderr
            try:
                if pipeline_log_fh is not None and not pipeline_log_fh.closed:
                    pipeline_log_fh.close()
            except (OSError, ValueError):
                pass
            # 복원된 터미널에 요약 1줄 출력 (사용자 피드백용)
            print(f"✅ UPDD v0.4 완료 — 최종 산출: {final_csv}")

if __name__ == "__main__":
    main()
