#!/usr/bin/env python
"""
utils/updd_cli.py
-----------------
UPDD 파이프라인의 터미널 UI(인터랙티브 메뉴 / 사용자 입력 수집) 레이어.

[v3 1-1] UPDD.py 의 입출력 헬퍼와 메뉴 함수들을 본 모듈로 분리하여,
오케스트레이션(파이프라인 단계 호출 흐름) 과 UI(사용자 상호작용) 의 책임을
명확히 분리한다.

본 모듈은 다음 함수를 노출한다:
    - arrow_menu(title, options) -> int
    - select_topology(topology_types) -> dict
    - select_ncaa(binder_len_str, registry_data, binder_chain="B", binder_chains=None)
        -> Tuple[Optional[NCAADef], str]
    - get_target_preprocess_options(target_pdb) -> dict
    - gather_all_inputs(base_name, target_pdb, topology_types, registry_data)
        -> Dict[str, Any]

주의:
    - UI 함수는 외부 의존성(도메인 데이터, 토폴로지 정의 등) 을 인자로 받아
      순수 UI 책임만 수행한다. 이로써 UPDD.py 의 import 사이클이 단절된다.
"""

import sys
import tty
import termios
import random
from typing import Any, Dict, List, Optional, Tuple


# ==========================================
# 인터랙티브 화살표 메뉴
# ==========================================
def arrow_menu(title: str, options: List[str]) -> int:
    """위/아래 화살표로 옵션을 선택하는 터미널 인터랙티브 메뉴.

    Args:
        title: 메뉴 상단에 출력할 타이틀
        options: 선택지 문자열 리스트

    Returns:
        int: 선택된 옵션의 0-based 인덱스
    """
    print(f"\n==================================================\n 💾 {title}\n==================================================")
    current = 0
    fd = sys.stdin.fileno()
    old_settings = termios.tcgetattr(fd)
    try:
        for _ in options:
            print()
        sys.stdout.write(f"\033[{len(options)}A")
        tty.setraw(sys.stdin.fileno())
        while True:
            for i, opt in enumerate(options):
                if i == current:
                    sys.stdout.write(f"\r\033[K  ▶ \033[92m{opt}\033[0m\n")
                else:
                    sys.stdout.write(f"\r\033[K    {opt}\n")
            sys.stdout.write(f"\033[{len(options)}A")
            sys.stdout.flush()
            ch = sys.stdin.read(1)
            if ch == '\x1b':
                sys.stdin.read(1)
                key = sys.stdin.read(1)
                if key == 'A':
                    current = max(0, current - 1)
                elif key == 'B':
                    current = min(len(options) - 1, current + 1)
            elif ch in ('\r', '\n'):
                break
            elif ch == '\x03':
                raise KeyboardInterrupt
    finally:
        termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        sys.stdout.write(f"\033[{len(options)}B\n")
        sys.stdout.flush()
    return current


def _print_step(msg: str) -> None:
    """[Local Helper] UI 모듈 내부에서 사용하는 step 헤더 출력 (UPDD.py 의 print_step 미러).

    UPDD.py 의 print_step 은 logger 와도 결합되어 있어, UI 모듈은 stdout 출력만
    수행하는 단순화된 사본을 사용한다. 인터랙티브 메뉴는 logger 가 아직 초기화
    되지 않은 시점에도 호출될 수 있으므로 logger 의존을 차단한다.
    """
    line = '=' * 60
    print('\n' + line + '\n ' + msg + '\n' + line + '\n')


# ==========================================
# Step 0: 타겟 전처리 옵션
# ==========================================
def get_target_preprocess_options(target_pdb: str) -> Dict[str, Any]:
    """타겟 PDB 전처리 옵션 (체인 / HETATM 유지 모드 / 강제 유지 잔기) 을 인터랙티브하게 수집한다."""
    _print_step("Step 0: Target Preprocessing 설정")
    if arrow_menu("원본 타겟 PDB 전처리를 자동화합니다.", ["[y] 전처리 실행 (물 분자 제거, END 레코드 추가 등)", "[n] 원본 그대로 유지"]) == 1:
        return {"target_pdb": target_pdb, "do_preprocess": False, "keep_mode": "none", "chain_in": "", "keep_res": ""}

    keep_mode = ["auto", "manual", "none"][arrow_menu("특수 HETATM 유지 모드", ["auto (기본값)", "manual", "none"])]
    c_opts = ["A", "B", "C", "D", "직접 입력"]
    c_idx = arrow_menu("유지할 타겟 체인 ID를 선택해 주십시오", c_opts)
    chain_in = input(" 타겟 체인을 직접 입력해 주십시오 (쉼표로 구분): ").strip() if c_idx == 4 else c_opts[c_idx]
    keep_res = ""
    if keep_mode in ["auto", "manual"]:
        keep_res = input(" 추가로 강제 유지할 잔기(residue) 이름들을 입력해 주십시오 (예: GDP,MG) [Enter=없음]: ").strip()

    return {"target_pdb": target_pdb, "do_preprocess": True, "keep_mode": keep_mode, "chain_in": chain_in, "keep_res": keep_res}


# ==========================================
# 토폴로지 / ncAA 선택
# ==========================================
def select_topology(topology_types: Dict[int, Dict[str, Any]]) -> Dict[str, Any]:
    """펩타이드 토폴로지 사전(linear/cyclic_*/bicyclic) 에서 하나를 선택한다.

    Args:
        topology_types: UPDD.py 의 TOPOLOGY_TYPES dict (1-based key → 메타데이터)
    """
    keys = list(topology_types.keys())
    opts = [f"{topology_types[k]['name']:<22} | {topology_types[k]['note']}" for k in keys]
    choice_idx = arrow_menu("펩타이드 토폴로지를 선택해 주십시오", opts)
    topo = topology_types[keys[choice_idx]]
    print(f" 선택됨: {topo['name']} (RFdiffusion T={topo['rf_T']} / Cyclic={topo['cyclic']})")
    return topo


def select_ncaa(binder_len_str: str,
                registry_data: List[Any],
                binder_chain: str = "B",
                binder_chains: Optional[List[str]] = None) -> Tuple[Optional[Any], str]:
    """ncAA 사용 여부 / 종류 / 삽입 위치를 인터랙티브하게 선택한다.

    Args:
        binder_len_str: "10-40" 형식의 길이 범위 문자열
        registry_data:  ncaa_registry.NCAA_REGISTRY_DATA (NCAADef 리스트)
        binder_chain:   기본 binder chain 문자
        binder_chains:  복수 binder chain 후보 (없으면 [binder_chain])

    Returns:
        Tuple[Optional[NCAADef], str]: (선택된 NCAADef 또는 None, 위치 문자열)
    """
    if arrow_menu("ncAA 사용 여부를 선택해 주십시오", ["0: 천연 아미노산만 사용", "1: ncAA 사용 (목록에서 선택)"]) == 0:
        print(" ncAA 없이 표준 아미노산으로 진행합니다.")
        return None, ""

    opts = [f"{spec.code:<5} | {spec.label:<25} | {spec.charge_model_note}" for spec in registry_data]
    choice_idx = arrow_menu("사용할 ncAA를 선택해 주십시오 (Registry 기반)", opts)
    selected_spec = registry_data[choice_idx]
    print(f" 선택됨: {selected_spec.label} (Code: {selected_spec.code})")

    max_len = 40
    try:
        max_len = int(binder_len_str.split("-")[-1])
    except (ValueError, IndexError) as _binder_err:
        print(f"  [WARN] binder_len 파싱 실패, 기본값(40) 사용: {_binder_err}")

    pos_opts = [
        "1: 특정 위치 직접 입력 (예: 3,7,12)",
        "2: 랜덤 삽입 (개수 지정)",
        "3: N-말단 (처음 N개 잔기)",
        "4: C-말단 (마지막 N개 잔기)",
        "5: 균등 분포 (N개, 간격 자동 계산)",
        "6: 지능형 자동 선택 (AUTO_SASA) — BSA + pLDDT 기반"
    ]
    pos_mode = arrow_menu(f"삽입 위치를 선택해 주십시오 (바인더 길이: {binder_len_str})", pos_opts)

    if pos_mode == 0:
        raw = input(" 위치 번호를 직접 입력해 주십시오 (예: 3,7): ").strip()
        positions = [int(x.strip()) for x in raw.split(",") if x.strip().isdigit()]
    elif pos_mode == 1:
        n = int(input(f" 삽입할 개수를 입력해 주십시오 (1~{max(1, max_len//3)}): ").strip() or "1")
        positions = sorted(random.sample(range(1, max_len + 1), min(n, max_len)))
    elif pos_mode == 2:
        n = int(input(" N-말단에서 몇 개를 삽입하시겠습니까?: ").strip() or "1")
        positions = list(range(1, n + 1))
    elif pos_mode == 3:
        n = int(input(" C-말단에서 몇 개를 삽입하시겠습니까?: ").strip() or "1")
        positions = list(range(max_len - n + 1, max_len + 1))
    elif pos_mode == 4:
        n = int(input(" 균등 배치할 개수를 입력해 주십시오: ").strip() or "2")
        step = max_len // (n + 1)
        positions = [step * (i + 1) for i in range(n)]
    else:
        n = int(input(" AUTO_SASA 선택 잔기 개수를 입력해 주십시오 (예: 3): ").strip() or "3")
        auto_sasa_str = f"AUTO_SASA:{n}"
        print(f" 지능형 선택 모드 활성화: {auto_sasa_str}")
        return selected_spec, auto_sasa_str

    binder_chains = binder_chains or [binder_chain]
    if len(binder_chains) > 1:
        bc_idx = arrow_menu("ncAA를 적용할 바인더 체인을 선택해 주십시오", [f"Chain {c}" for c in binder_chains])
        binder_chain = binder_chains[bc_idx]

    positions_str = ",".join([binder_chain + str(p) for p in positions])
    print(f" 최종 치환 잔기: {positions_str}")
    return selected_spec, positions_str


# ==========================================
# 통합 입력 수집 (gather_all_inputs)
# ==========================================
def gather_all_inputs(base_name: str,
                      target_pdb: str,
                      topology_types: Dict[int, Dict[str, Any]],
                      registry_data: List[Any]) -> Dict[str, Any]:
    """파이프라인 1회 실행에 필요한 모든 사용자 입력을 인터랙티브하게 수집한다.

    Args:
        base_name:      타겟 PDB 의 base name (확장자 제거)
        target_pdb:     원본 타겟 PDB 경로
        topology_types: TOPOLOGY_TYPES dict (UPDD.py 에서 주입)
        registry_data:  NCAA_REGISTRY_DATA (UPDD.py 에서 주입)

    Returns:
        Dict[str, Any]: UPDD.py main() 가 status_file 에 저장할 inputs dict
    """
    preprocess_opts = get_target_preprocess_options(target_pdb)
    _print_step("Step 1: Hotspot Selection")
    hotspots = input(" 공략할 아미노산 번호를 쉼표로 구별하여 입력해 주십시오 (예: 68,71,72) [Enter=없음]: ").strip()
    grafting_needed = arrow_menu("Step 2: Ligand Grafting (결합 부위 모방)", ["[n] 필요 없음 (기본값)", "[y] Grafting 수행"]) == 1
    ref_pdb = ""
    if grafting_needed:
        ref_pdb = input(" Reference PDB 경로를 입력해 주십시오 (예: 4G0N.pdb): ").strip()

    _print_step("Step 3: Peptide Topology & Shape")
    topology = select_topology(topology_types)
    binder_len = input(" 바인더 길이를 입력해 주십시오 (예: 8-12) [기본값: 10-40]: ").strip() or "10-40"

    c_opts = ["A", "B", "C", "D", "직접 입력"]
    c_idx = arrow_menu("타겟 체인 ID를 선택해 주십시오", c_opts)
    target_chain = input(" 타겟 체인을 직접 입력해 주십시오: ").strip() if c_idx == 4 else c_opts[c_idx]

    n_designs = input(" 디자인 생성 개수를 입력해 주십시오 [기본값: 5]: ").strip() or "5"
    # [v43 FIX — Bug 2] MPNN 이 각 backbone 에 2 개 서열을 생성하므로 최종 AF2 후보는 2 배.
    mpnn_samples = 2
    total_candidates = int(n_designs) * mpnn_samples
    print(f"  → RFdiffusion {n_designs}개 backbone × ProteinMPNN {mpnn_samples}개 서열 = "
          f"총 {total_candidates}개 AF2 후보 (의도된 다양성 확보)")

    print("\n==================================================")
    print("  그 외 파라미터 사전 설정")
    print("==================================================")
    admet_mode = ["none", "cell", "bbb"][arrow_menu("[ADMET 투과성 필터] 선택", ["none (기본값)", "cell (세포막 투과성)", "bbb (뇌혈관 장벽)"])]

    default_binder = chr(ord(target_chain[0].upper()) + 1) if target_chain else 'B'
    ncaa_def, ncaa_positions = select_ncaa(binder_len, registry_data, binder_chain=default_binder, binder_chains=[default_binder])

    qmmm_mode = ["fast", "full", "plaid"][arrow_menu("[QM/MM 계산 모드] 선택", ["fast (기본값)", "full", "plaid"])]

    gpu_opts = ["1개(기본값)", "2개", "3개", "4개", "8개", "직접 입력"]
    gpu_idx = arrow_menu("[GPU 병렬화] 동시 가동 수 선택", gpu_opts)
    if gpu_idx == 5:
        parallel_workers = int(input(" 동시 가동 수를 입력해 주십시오: ").strip() or 1)
    else:
        parallel_workers = [1, 2, 3, 4, 8][gpu_idx]

    msa_choice = ["local", "server"][arrow_menu("[AF2 MSA 탐색] 모드 선택", ["local (기본값)", "server"])]
    hdd_name = input("\n 궤적(DCD) 백업용 외장저장장치 이름을 입력해 주십시오 (예: HDD, 없으면 Enter): ").strip()

    iptm_cutoff = input("\n ipTM 필터 임계값 (기본 0.3, Enter로 스킵): ").strip()
    iptm_cutoff = float(iptm_cutoff) if iptm_cutoff else 0.3

    # Pre-filter RMSD threshold (Step 6.5 RMSD sub-filter).
    # Kabsch Cα RMSD of AF2 rank_001 chain A vs crystal — rejects designs
    # where AF2 predicted the target fold too far from the reference.
    # Default 2.5 Å accepts the full AF2-multimer noise envelope (Bryant
    # 2022 Nat Commun IQR 1.5-3.2 Å); set 2.0 Å for stricter pre-screen.
    rmsd_cutoff_in = input(
        "\n 구조 RMSD 필터 임계값 (Å, 기본 2.5, Enter로 스킵; "
        "더 엄격히 하려면 2.0 입력): "
    ).strip()
    try:
        rmsd_cutoff = float(rmsd_cutoff_in) if rmsd_cutoff_in else 2.5
    except ValueError:
        rmsd_cutoff = 2.5

    return {
        "preprocess_opts": preprocess_opts,
        "hotspots": hotspots,
        "grafting_needed": grafting_needed,
        "ref_pdb": ref_pdb,
        "topology": topology,
        "binder_len": binder_len,
        "target_chain": target_chain,
        "n_designs": n_designs,
        "admet_mode": admet_mode,
        "ncaa_label": ncaa_def.label if ncaa_def else "none",
        "ncaa_code": ncaa_def.code if ncaa_def else "none",
        "ncaa_positions": ncaa_positions,
        "binder_chain": default_binder,
        "qmmm_mode": qmmm_mode,
        "parallel_workers": parallel_workers,
        "msa_choice": msa_choice,
        "hdd_name": hdd_name,
        "iptm_cutoff": iptm_cutoff,
        "rmsd_cutoff": rmsd_cutoff,
    }
