#!/usr/bin/env bash
#
# scripts/setup_vm.sh — v0.9 dual-GPU 통합: VM (V100) 환경 자동 셋업
#
# 동작:
#   1. VM SSH 연결 + V100 + disk + cuda toolkit pre-flight
#   2. Miniconda 자동 설치 (없으면)
#   3. mamba 설치 (faster solver)
#   4. host 의 qmmm + md_simulation conda env 미러 생성
#   5. UPDD source rsync (utils/scripts/UPDD.py/target_cards/params)
#   6. 모든 import sanity check (openmm CUDA, gpu4pyscf, cupy, mdtraj, parmed, ambertools 등)
#   7. 에러 자동 해결 (mamba retry, channel fallback, pip install fallback)
#
# 사용법:
#   bash scripts/setup_vm.sh                      # 자동 진행
#   UPDD_VM_SSH_TARGET=san@... bash scripts/setup_vm.sh
#   bash scripts/setup_vm.sh --recreate           # 기존 env 삭제 후 재생성
#   bash scripts/setup_vm.sh --skip-source-sync   # source rsync 건너뛰기
#
# Conda env 이름: qmmm + md_simulation (host 와 동일)
# Conda prefix on VM: /home/san/miniconda3 (host 와 동일 경로)

set -euo pipefail

# ──────────────────────────────────────────────────────────────
# Config
# ──────────────────────────────────────────────────────────────
VM_TARGET="${UPDD_VM_SSH_TARGET:-san@192.168.122.155}"
HOST_PROJECT_ROOT="${UPDD_HOST_PROJECT_ROOT:-/home/san/UPDD_proj}"
VM_PROJECT_ROOT="${UPDD_VM_PROJECT_ROOT:-/home/san/UPDD_proj}"
VM_CONDA_PREFIX="${UPDD_VM_CONDA_PREFIX:-/home/san/miniconda3}"
VM_SCRATCH_ROOT="${UPDD_VM_SCRATCH_ROOT:-/var/scratch}"   # VM 내 PYSCF scratch + tmp
LOG="/tmp/setup_vm_$(date +%Y%m%d_%H%M%S).log"
RECREATE_ENVS=0
SKIP_SOURCE_SYNC=0
MINICONDA_INSTALLER_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"

# argument parsing
for arg in "$@"; do
    case "$arg" in
        --recreate) RECREATE_ENVS=1 ;;
        --skip-source-sync) SKIP_SOURCE_SYNC=1 ;;
        --help|-h)
            sed -n '1,30p' "$0"
            exit 0
            ;;
    esac
done

# ──────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────
log() { printf "[%(%H:%M:%S)T] %s\n" -1 "$*" | tee -a "$LOG"; }
err() { printf "[%(%H:%M:%S)T] ERROR: %s\n" -1 "$*" | tee -a "$LOG" >&2; }
ok()  { printf "[%(%H:%M:%S)T] ✓ %s\n" -1 "$*" | tee -a "$LOG"; }

vm_run() {
    # 단일 명령 SSH 실행. stdout/stderr → log. pipefail 로 ssh 실패 propagate.
    local cmd="$1"
    set -o pipefail
    ssh -o ConnectTimeout=10 -o StrictHostKeyChecking=no -o LogLevel=ERROR \
        "$VM_TARGET" "$cmd" 2>&1 | tee -a "$LOG"
    local rc=${PIPESTATUS[0]}
    return "$rc"
}

vm_run_quiet() {
    # 결과만 반환 (log 에 안 남김). exit code 그대로.
    ssh -o ConnectTimeout=10 -o StrictHostKeyChecking=no -o LogLevel=ERROR \
        "$VM_TARGET" "$1" 2>/dev/null
}

vm_run_exit_only() {
    ssh -o ConnectTimeout=10 -o StrictHostKeyChecking=no -o LogLevel=ERROR \
        "$VM_TARGET" "$1" >/dev/null 2>&1
}

retry_vm() {
    # 명령 N회 재시도 (간격 5s/30s/120s)
    local cmd="$1"; local label="$2"
    local delays=(5 30 120)
    for delay in "${delays[@]}"; do
        if vm_run "$cmd"; then return 0; fi
        err "$label 실패 — ${delay}s 후 재시도"
        sleep "$delay"
    done
    err "$label 3회 재시도 모두 실패"
    return 1
}

# ──────────────────────────────────────────────────────────────
# Phase 1: Pre-flight
# ──────────────────────────────────────────────────────────────
log "========================================"
log "UPDD VM (V100) 환경 자동 셋업 시작"
log "VM target: $VM_TARGET"
log "Log file:  $LOG"
log "========================================"

log ""
log "[1/6] Pre-flight check"

# SSH 연결
if ! vm_run_exit_only 'echo ok'; then
    err "VM SSH 연결 실패 — VM 가 실행 중인지 확인 (virsh list --all)"
    exit 1
fi
ok "SSH 연결 OK"

# V100 확인
if ! vm_run_exit_only 'nvidia-smi --query-gpu=name --format=csv,noheader | grep -qi V100'; then
    err "V100 미감지 — VFIO 패스스루 확인"
    vm_run 'nvidia-smi -L || true'
    exit 1
fi
GPU_INFO=$(vm_run_quiet 'nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader')
ok "GPU 감지: $GPU_INFO"

# CUDA toolkit 확인
CUDA_VERSION=$(vm_run_quiet 'ls -d /usr/local/cuda-12* 2>/dev/null | head -1 | xargs basename 2>/dev/null || echo none')
if [[ "$CUDA_VERSION" == "none" ]]; then
    err "CUDA 12.x toolkit 미설치 (/usr/local/cuda-12* 부재)"
    log "→ 자동 설치 필요. cudatoolkit 은 conda env 가 자체 포함하므로 진행 가능, 단 cupy-cuda12x 는 system CUDA 필요 가능"
fi
ok "CUDA: $CUDA_VERSION"

# Disk space (최소 30GB 권장)
DISK_AVAIL_GB=$(vm_run_quiet "df -BG --output=avail / | tail -1 | tr -d 'G '")
DISK_AVAIL_GB="${DISK_AVAIL_GB:-0}"
if [[ "$DISK_AVAIL_GB" -lt 30 ]]; then
    err "Disk 부족 (${DISK_AVAIL_GB} GB free, 30 GB 권장)"
    vm_run 'df -h /'
    log "→ 진행 중단. VM disk 확장 후 재실행."
    exit 1
fi
ok "Disk free: ${DISK_AVAIL_GB} GB"

# ──────────────────────────────────────────────────────────────
# Phase 2: Miniconda 설치
# ──────────────────────────────────────────────────────────────
log ""
log "[2/6] Miniconda 셋업"

if vm_run_exit_only "test -x $VM_CONDA_PREFIX/bin/conda"; then
    ok "Miniconda 이미 설치됨 ($VM_CONDA_PREFIX)"
else
    log "Miniconda 미설치 — 자동 설치 진행"
    vm_run "cd /tmp && wget -nv -O miniconda.sh $MINICONDA_INSTALLER_URL && bash miniconda.sh -b -p $VM_CONDA_PREFIX && rm miniconda.sh" \
        || { err "Miniconda 설치 실패"; exit 1; }
    vm_run "$VM_CONDA_PREFIX/bin/conda init bash"
    ok "Miniconda 설치 완료"
fi

# conda config (--set 와 --add 를 한 invocation 에 섞을 수 없음 → 분리)
vm_run "$VM_CONDA_PREFIX/bin/conda config --set always_yes true || true"
vm_run "$VM_CONDA_PREFIX/bin/conda config --set channel_priority flexible || true"
vm_run "$VM_CONDA_PREFIX/bin/conda config --add channels conda-forge || true"

# Anaconda Terms of Service 자동 수락 (2026 정책 — 미수락 시 env create 차단)
log "Anaconda channel ToS 자동 수락"
vm_run "$VM_CONDA_PREFIX/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main || true"
vm_run "$VM_CONDA_PREFIX/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r || true"
vm_run "$VM_CONDA_PREFIX/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/msys2 || true"

# defaults 채널을 제거하고 conda-forge 만 사용 — Anaconda ToS 회피 + 일관성
log "defaults 채널 제거 (conda-forge only)"
vm_run "$VM_CONDA_PREFIX/bin/conda config --remove channels defaults 2>/dev/null || true"

# mamba (faster solver)
if ! vm_run_exit_only "test -x $VM_CONDA_PREFIX/bin/mamba"; then
    log "mamba 설치 (faster solver)"
    vm_run "$VM_CONDA_PREFIX/bin/conda install -n base -c conda-forge mamba -y" \
        || log "mamba 설치 실패 — conda 로 fallback"
fi
if vm_run_exit_only "test -x $VM_CONDA_PREFIX/bin/mamba"; then
    SOLVER="$VM_CONDA_PREFIX/bin/mamba"
    ok "Solver: mamba"
else
    SOLVER="$VM_CONDA_PREFIX/bin/conda"
    ok "Solver: conda (mamba unavailable)"
fi

# ──────────────────────────────────────────────────────────────
# Phase 3: Source code rsync
# ──────────────────────────────────────────────────────────────
log ""
log "[3/6] UPDD source rsync"

if [[ "$SKIP_SOURCE_SYNC" -eq 1 ]]; then
    log "--skip-source-sync 지정 — source rsync 건너뛰기"
else
    vm_run "mkdir -p $VM_PROJECT_ROOT"
    # rsync: utils + scripts + UPDD.py + target_cards + environment.yml. outputs/archive/log 등 제외.
    rsync -avz --delete-after \
        --exclude='outputs/' --exclude='archive/' --exclude='log/' \
        --exclude='__pycache__/' --exclude='.git/' --exclude='_archive/' \
        --exclude='_audit/' --exclude='analysis/' --exclude='chkfile_archive/' \
        --exclude='*.dcd' --exclude='*.nc' --exclude='*.h5' --exclude='*.chk' \
        --exclude='paper1*' --exclude='inputs/' --exclude='plan/' \
        -e 'ssh -o StrictHostKeyChecking=no' \
        "$HOST_PROJECT_ROOT/" \
        "$VM_TARGET:$VM_PROJECT_ROOT/" 2>&1 | tail -20 | tee -a "$LOG"
    ok "Source rsync 완료"
fi

# ──────────────────────────────────────────────────────────────
# Phase 4: qmmm env 생성 (host 의 environment.yml 미러)
# ──────────────────────────────────────────────────────────────
log ""
log "[4/6] qmmm conda env 생성"

if vm_run_exit_only "test -d $VM_CONDA_PREFIX/envs/qmmm" && [[ "$RECREATE_ENVS" -eq 0 ]]; then
    ok "qmmm env 이미 존재 — skip (재생성하려면 --recreate)"
else
    if [[ "$RECREATE_ENVS" -eq 1 ]] && vm_run_exit_only "test -d $VM_CONDA_PREFIX/envs/qmmm"; then
        log "기존 qmmm env 삭제"
        vm_run "$VM_CONDA_PREFIX/bin/conda env remove -n qmmm -y"
    fi

    # qmmm env 생성 (host 의 environment.yml 그대로 사용 — conda-forge only)
    log "qmmm env 생성 중 (3-10 분 예상)"
    if ! retry_vm "$SOLVER env create -f $VM_PROJECT_ROOT/environment.yml --quiet" "qmmm env 생성"; then
        err "qmmm env 생성 3회 실패 — pip-only fallback 시도 (가장 안전)"
        # Fallback: 최소 spec 으로 재시도 (conda-forge only, no defaults)
        vm_run "$VM_CONDA_PREFIX/bin/conda env create -n qmmm --channel conda-forge --file $VM_PROJECT_ROOT/environment.yml" \
            || { err "qmmm env 생성 최종 실패"; exit 1; }
    fi
    ok "qmmm env 생성 완료"
fi

# ──────────────────────────────────────────────────────────────
# Phase 5: md_simulation env 생성 (heredoc spec)
# ──────────────────────────────────────────────────────────────
log ""
log "[5/6] md_simulation conda env 생성"

# heredoc 으로 VM 안에 spec 전달
MD_SIM_SPEC=$(cat <<'YAML'
# UPDD md_simulation env spec (VM mirror, derived from host conda env export --from-history)
# Used by Stage 9 (Restrained MD), Stage 10 (Snapshot Extract), Stage 12 (MM-GBSA via run_mmgbsa.py)
name: md_simulation
channels:
  - conda-forge
dependencies:
  - python=3.11
  - openmm=8.2
  - openmmforcefields
  - openmmtools
  - mdtraj=1.11
  - parmed=4.3
  - pdbfixer
  - ambertools=24.8
  - cudatoolkit
  - networkx
  - scipy
  - scikit-learn
  - numpy
  - pandas
  - matplotlib-base
  - pyyaml
  - psutil
  - tqdm
  - gfortran
  - gcc
  - gxx
  - pip
  - pip:
      - pytest
YAML
)

if vm_run_exit_only "test -d $VM_CONDA_PREFIX/envs/md_simulation" && [[ "$RECREATE_ENVS" -eq 0 ]]; then
    ok "md_simulation env 이미 존재 — skip"
else
    if [[ "$RECREATE_ENVS" -eq 1 ]] && vm_run_exit_only "test -d $VM_CONDA_PREFIX/envs/md_simulation"; then
        log "기존 md_simulation env 삭제"
        vm_run "$VM_CONDA_PREFIX/bin/conda env remove -n md_simulation -y"
    fi
    # spec 을 VM 으로 전송 후 생성
    echo "$MD_SIM_SPEC" | ssh -o StrictHostKeyChecking=no "$VM_TARGET" "cat > /tmp/md_simulation_vm.yml"
    log "md_simulation env 생성 중 (5-15 분 예상, ambertools 가 크다)"
    if ! retry_vm "$SOLVER env create -f /tmp/md_simulation_vm.yml --quiet" "md_simulation env 생성"; then
        err "md_simulation env 생성 3회 실패 — conda-forge only fallback"
        vm_run "$VM_CONDA_PREFIX/bin/conda env create -n md_simulation --channel conda-forge --file /tmp/md_simulation_vm.yml" \
            || { err "md_simulation env 생성 최종 실패"; exit 1; }
    fi
    ok "md_simulation env 생성 완료"
fi

# ──────────────────────────────────────────────────────────────
# Phase 6: Sanity check + scratch 디렉토리
# ──────────────────────────────────────────────────────────────
log ""
log "[6/6] Sanity check"

# Scratch + chkfile 디렉토리 (VM 측, SSD 가 없으니 root partition 사용)
vm_run "mkdir -p $VM_SCRATCH_ROOT/pyscf_scratch $VM_SCRATCH_ROOT/tmp $VM_SCRATCH_ROOT/mmpbsa_scratch && chmod 700 $VM_SCRATCH_ROOT/*"
ok "VM scratch 디렉토리 ($VM_SCRATCH_ROOT/*) 준비"

# qmmm env sanity
log "qmmm env import check"
qmmm_imports=$(vm_run_quiet "$VM_CONDA_PREFIX/envs/qmmm/bin/python -c '
import sys
ok = True
try:
    import openmm; print(\"openmm\", openmm.__version__)
except Exception as e: print(\"openmm FAIL\", e); ok=False
try:
    import pyscf; print(\"pyscf\", pyscf.__version__)
except Exception as e: print(\"pyscf FAIL\", e); ok=False
try:
    import gpu4pyscf; print(\"gpu4pyscf\", gpu4pyscf.__version__)
except Exception as e: print(\"gpu4pyscf FAIL\", e); ok=False
try:
    import cupy; print(\"cupy\", cupy.__version__, \"cuda_count=\", cupy.cuda.runtime.getDeviceCount())
except Exception as e: print(\"cupy FAIL\", e); ok=False
try:
    import mdtraj; print(\"mdtraj\", mdtraj.__version__)
except Exception as e: print(\"mdtraj FAIL\", e); ok=False
try:
    import parmed; print(\"parmed\", parmed.__version__)
except Exception as e: print(\"parmed FAIL\", e); ok=False
sys.exit(0 if ok else 1)
'")
echo "$qmmm_imports" | tee -a "$LOG"
if echo "$qmmm_imports" | grep -q FAIL; then
    err "qmmm env import 실패 detected"
    QMMM_OK=0
else
    ok "qmmm env imports OK"
    QMMM_OK=1
fi

# md_simulation env sanity
log ""
log "md_simulation env import check"
md_imports=$(vm_run_quiet "$VM_CONDA_PREFIX/envs/md_simulation/bin/python -c '
import sys
ok = True
try:
    import openmm; print(\"openmm\", openmm.__version__)
except Exception as e: print(\"openmm FAIL\", e); ok=False
try:
    import mdtraj; print(\"mdtraj\", mdtraj.__version__)
except Exception as e: print(\"mdtraj FAIL\", e); ok=False
try:
    import parmed; print(\"parmed\", parmed.__version__)
except Exception as e: print(\"parmed FAIL\", e); ok=False
try:
    import openmmforcefields; print(\"openmmforcefields\", openmmforcefields.__version__)
except Exception as e: print(\"openmmforcefields FAIL\", e); ok=False
sys.exit(0 if ok else 1)
'")
echo "$md_imports" | tee -a "$LOG"
if echo "$md_imports" | grep -q FAIL; then
    err "md_simulation env import 실패 detected"
    MD_OK=0
else
    ok "md_simulation env imports OK"
    MD_OK=1
fi

# ambertools CLI (in md_simulation env)
log "ambertools CLI sanity (md_simulation)"
ambertools_check=$(vm_run_quiet "$VM_CONDA_PREFIX/envs/md_simulation/bin/antechamber -h 2>&1 | head -3 || $VM_CONDA_PREFIX/envs/md_simulation/bin/tleap -h 2>&1 | head -3" || true)
if [[ -n "$ambertools_check" ]]; then
    ok "ambertools CLI 사용 가능"
else
    err "ambertools CLI 미감지"
fi

# OpenMM CUDA platform check (qmmm env)
log "OpenMM CUDA platform check (qmmm env)"
cuda_check=$(vm_run_quiet "$VM_CONDA_PREFIX/envs/qmmm/bin/python -c '
import openmm
plat = openmm.Platform.getPlatformByName(\"CUDA\")
print(\"OpenMM CUDA platform OK, num_devices=\", plat.getPropertyDefaultValue(\"DeviceIndex\"))
'" 2>&1)
echo "$cuda_check" | tee -a "$LOG"
if echo "$cuda_check" | grep -q "OK"; then
    ok "OpenMM CUDA platform OK"
else
    err "OpenMM CUDA platform check 실패"
fi

# gpu4pyscf SCF dry-run (very small system, ~10 atoms)
log "gpu4pyscf SCF dry-run (H2O minimal basis)"
gpu4pyscf_dryrun=$(vm_run_quiet "$VM_CONDA_PREFIX/envs/qmmm/bin/python -c '
from pyscf import gto
import gpu4pyscf.dft as dft
mol = gto.M(atom=\"O 0 0 0; H 0 0 1; H 0 1 0\", basis=\"sto-3g\", verbose=0)
mf = dft.RKS(mol, xc=\"b3lyp\").to_gpu()
e = mf.kernel()
print(\"SCF energy:\", e, \"converged:\", mf.converged)
' 2>&1" || echo "DRYRUN FAIL")
echo "$gpu4pyscf_dryrun" | tee -a "$LOG"
if echo "$gpu4pyscf_dryrun" | grep -q "converged: True"; then
    ok "gpu4pyscf V100 dry-run SUCCESS"
else
    err "gpu4pyscf V100 dry-run 실패 — full output 위 참조"
fi

# ──────────────────────────────────────────────────────────────
# Final Summary
# ──────────────────────────────────────────────────────────────
log ""
log "========================================"
log "셋업 완료"
log "  qmmm env:          $([ $QMMM_OK -eq 1 ] && echo OK || echo FAIL)"
log "  md_simulation env: $([ $MD_OK   -eq 1 ] && echo OK || echo FAIL)"
log "  GPU:               $GPU_INFO"
log "  Scratch:           $VM_SCRATCH_ROOT/{pyscf_scratch,tmp,mmpbsa_scratch}"
log "  Log:               $LOG"
log "========================================"

if [[ $QMMM_OK -eq 1 && $MD_OK -eq 1 ]]; then
    log "✓ Phase 2 Step 2 완료 — 다음 단계 진행 OK"
    exit 0
else
    err "✗ 일부 env 가 실패 — log 확인 후 --recreate 로 재시도"
    exit 1
fi
