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

# Registry Import (단일 계약)
sys.path.append(UTILS_DIR)
try:
    from ncaa_registry import resolve_ncaa_definition, NCAA_REGISTRY_DATA
except ImportError:
    sys.exit("[!] utils/ncaa_registry.py 파일을 찾을 수 없습니다. 경로를 확인하세요.")

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

def run_target_preprocessing(target_pdb: str, keep_mode: str, chain_in: str, keep_res: str) -> str:
    output_clean = os.path.splitext(target_pdb)[0] + "_clean.pdb"
    args = ["--input", target_pdb, "--output", output_clean, "--chains", chain_in, "--hetatm_mode", keep_mode]
    if keep_res:
        args += ["--keep_hetatms", keep_res]
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
        with concurrent.futures.ThreadPoolExecutor(max_workers=parallel_workers) as executor:
            futures = []
            for i in range(parallel_workers):
                w_designs = total_designs // parallel_workers + (1 if i < total_designs % parallel_workers else 0)
                if w_designs > 0:
                    futures.append(executor.submit(run_rfd_worker, i + 1, w_designs))
            concurrent.futures.wait(futures)
    else:
        run_rfd_worker(1, total_designs)

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

def _run_parameterization(struct_dir, param_dir, ncaa_def):
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

def _run_md_and_snapshots(struct_dir, params_manifest, target_out_dir, ncaa_label, ncaa_code, canonical_topo, binder_chain, graph_policy="strict", dcd_backup_path=""):
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
        
        if not ncaa_code or ncaa_code == "none":
            md_args += ["--ncaa_label", "none", "--ncaa_code", ""]
        else:
            md_args += ["--params_manifest", params_manifest, "--ncaa_label", ncaa_label, "--ncaa_code", ncaa_code]
            
        md_args += ["--dt_fs", "2.0"]
        run_conda_command(SCRIPT_PATHS["run_md"], md_args, "Restrained MD (Pass 1: 2fs)", env_key="md")

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
            # 기존 EXPLODED 마커 정리 + _final.pdb 삭제 로직은 유지
            for elog in exploded_2fs:
                base = os.path.basename(elog).replace("_EXPLODED_dt2fs.log", "")
                for ext in ["_restrained.dcd", "_final.pdb", "_md.log"]:
                    target_f = os.path.join(md_out_dir, base + ext)
                    if os.path.exists(target_f):
                        os.remove(target_f)
                os.remove(elog)

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
                 parallel_workers: int, target_out_dir: str, binder_chain: str = "B") -> str:
    print_step("Step 11: QM/MM — wB97X-D / PySCF")
    qmmm_out = os.path.join(target_out_dir, "qmmm_results")
    os.makedirs(qmmm_out, exist_ok=True)
    elem_arg = ncaa_element if ncaa_element else "none"

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
        """해당 디자인의 모든 스냅샷 QM/MM이 정상 완료되었는지 확인."""
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
            "--no-capture-output", "python", SCRIPT_PATHS["run_qmmm"],
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
            cmd = base_cmd + qmmm_args

            # 짧은 태그: "design_w4_1_s2" → "w4_1_s2"
            tag = design_id.replace("design_", "").split("_unrelaxed")[0]
            if len(tag) > 10:
                tag = tag[:10]

            try:
                # [v0.3 #2] PYTHONUNBUFFERED=1 — PIPE 통과 시 블록 버퍼링 방지.
                # DFT cycle 로그가 즉시 [tag] prefix 와 함께 표시되어 watchdog 파싱과
                # 운영자 모니터링이 동시 가능.
                proc = subprocess.Popen(
                    cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                    text=True, bufsize=1, env=_CHILD_ENV_UNBUFFERED,
                )
                oom_detected = False
                for line in proc.stdout:
                    stripped = line.rstrip()
                    if stripped:
                        print(f"  [{tag}] {stripped}")
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
        current_workers = min(parallel_workers, len(pending_ids))
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
                    for future in concurrent.futures.as_completed(future_map):
                        result = future.result()
                        pass_results.append(result)
                        if result["status"] == "SUCCESS":
                            print(f"  ✅ {result['id']} 완료")
                        elif result["status"] == "OOM":
                            print(f"  ⚠️  {result['id']} OOM — 다음 pass에서 재시도")
                        else:
                            print(f"  ❌ {result['id']} 에러: {result['msg']}")
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
        "md":           _stage("Step 9 Restrained MD",  md_dir,      lambda p: os.path.isdir(p) and _has_any(os.path.join(p, "*.dcd"))),
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
    _parser = argparse.ArgumentParser(add_help=False)
    _parser.add_argument("--validate", action="store_true", help="실행 없이 결과 파일 검증만 수행")
    _args, _ = _parser.parse_known_args()
    VALIDATE_ONLY = _args.validate

    setup_logger("session_start")
    _builtins.input = logged_input

    print_step("UPDD Pipeline v28: Architectural & Logic Fusion")

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

    # ── [RESTORED] save_step 헬퍼 (resume 시 skip 판단용) ──
    def save_step(step_name):
        """완료된 단계를 updd_status.json에 기록한다 (resume 시 skip 판단용)."""
        try:
            with open(status_file, "r", encoding="utf-8") as f:
                _st = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            _st = {}
        if step_name not in _st.setdefault("completed_steps", []):
            _st["completed_steps"].append(step_name)
        with open(status_file, "w", encoding="utf-8") as f:
            json.dump(_st, f, indent=4)

    # ── [코어 파 파이프라인 가동] ──
    if not VALIDATE_ONLY:
        safe_contig = get_continuous_contigs(work_pdb, target_chain)
        final_contig = "'" + safe_contig + "/0 " + inputs.get("binder_len", "10-40") + "'"
        run_rfdiffusion(work_pdb, topology_mode, inputs.get("hotspots", ""), final_contig, design_dir, int(inputs.get("n_designs", 5)), parallel_workers, target_chain)
        save_step("rfdiffusion")

    mpnn_out_dir, jsonl_path = run_proteinmpnn(design_dir, target_chain, topology_mode, valid_only=VALIDATE_ONLY)
    save_step("proteinmpnn")
    # [v5 해법 C] topology_mode/work_pdb/binder_chain 을 전달하여 cyclic 분기를 활성화한다.
    af2_out_dir = run_alphafold(
        design_dir, mpnn_out_dir, jsonl_path, target_chain, msa_mode, USE_LOCAL_MSA,
        valid_only=VALIDATE_ONLY,
        topology_mode=topology_mode,
        work_pdb=work_pdb,
        binder_chain=binder_chain,
    )
    save_step("alphafold")

    if VALIDATE_ONLY:
        validate_pipeline(base_name, design_dir, mpnn_out_dir, af2_out_dir, target_out_dir + "/ncaa_mutants" if ncaa_def else "", target_out_dir + "/params", target_out_dir + "/mdresult", target_out_dir + "/snapshots", target_out_dir + "/qmmm_results", target_out_dir + "/mmgbsa_results", os.path.join(target_out_dir, "final_ranking.csv"), ncaa_def.code if ncaa_def else None)
        return

    print_step("Step 6.5: AF2 Top Rank Filter (rank_001)")
    _filter_top_rank(af2_out_dir)

    print_step("Step 6.6: ipTM Pre-filter")
    _iptm_cutoff = float(inputs.get("iptm_cutoff", 0.3))
    _filter_by_iptm(af2_out_dir, min_iptm=_iptm_cutoff)

    print_step("Step 7: ncAA Configuration & Mutation")
    struct_dir = _run_ncaa_mutation(af2_out_dir, ncaa_def.code if ncaa_def else "none", ncaa_def.label if ncaa_def else "none", ncaa_positions, target_out_dir, binder_chain, target_chain)
    save_step("ncaa_mutation")

    print_step("Step 8: ncAA Parameterization")
    param_dir = os.path.join(target_out_dir, "params")
    params_manifest = _run_parameterization(struct_dir, param_dir, ncaa_def)

    _run_admet(struct_dir, admet_mode, target_out_dir, status if not resume_mode else inputs, status_file)

    print_step("Step 9: Restrained MD")
    canonical_topo = topology_mode.get("canonical", topology_mode.get("abbr", "linear"))
    md_out_dir, snap_dir = _run_md_and_snapshots(struct_dir, params_manifest, target_out_dir, ncaa_def.label if ncaa_def else "none", ncaa_def.code if ncaa_def else "", canonical_topo, binder_chain, graph_policy="strict", dcd_backup_path=dcd_backup_path)
    save_step("md_snapshots")

    qmmm_out = execute_qmmm(snap_dir, qmmm_mode, ncaa_def.element if ncaa_def else "", af2_out_dir, parallel_workers, target_out_dir, binder_chain)
    save_step("qmmm")

    # ── [RESTORED] Step 12: MM-GBSA 결합 자유에너지 계산 ──
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
    run_conda_command(SCRIPT_PATHS["run_mmgbsa"], mmgbsa_args, "MM-GBSA ΔG 계산", env_key="md")
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
    run_conda_command(SCRIPT_PATHS["rank_qmmm"], rank_args, "최종 후보 Ranking", env_key="qmmm")
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

if __name__ == "__main__":
    main()
