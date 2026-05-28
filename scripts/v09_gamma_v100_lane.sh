#!/usr/bin/env bash
# ============================================================================
# v09_gamma_v100_lane.sh — V100 parallel lane for γ pilot (Stage A only)
# ----------------------------------------------------------------------------
# Companion to scripts/v09_gamma_production.sh which runs on host 5070 Ti
# sequentially. This script runs the REVERSE-order half of Cp4 seeds on VM
# V100 (Tesla V100-PCIE-32GB via VFIO passthrough, ssh san@192.168.122.155).
#
# Work-splitting strategy:
#   Host launcher (already running): s7 → s23 → s42 → s83 → ... (forward)
#   V100 lane (this script):         s19_reseed55 → s251 → s163 → s101 (reverse)
# 두 lane 이 만나는 지점 (s101 이후) 에서 host 가 "final.pdb exists, skip" → 자연 splitter
#
# Per-seed flow:
#   1. Host 측 seedir 준비 (_md_input + γ XML + manifest 복사, host launcher 와 동일)
#   2. seedir → VM rsync (mdresult 제외)
#   3. VM 에서 conda run -n md_simulation python utils/run_restrained_md.py ...
#   4. VM mdresult → host rsync (final.pdb 가 host launcher's skip 신호 됨)
#
# Output logs co-locate w/ host pilot: outputs/analysis/v09_gamma_production_<TS>/
#
# Pre-req:
#   - setup_vm.sh 가 실행되어 VM 에 conda envs (qmmm + md_simulation) 존재
#   - VM 에 /home/san/UPDD_proj/utils/ 가 rsync 되어 있음 (setup_vm.sh Phase 3)
#   - host γ pilot 이 이미 시작되어 logdir 존재 (이 script 가 logdir 자동 감지)
# ============================================================================
set -u

PROJ=/home/san/UPDD_proj
VM_TARGET="${UPDD_VM_SSH_TARGET:-san@192.168.122.155}"
VM_PROJ="${UPDD_VM_PROJECT_ROOT:-/home/san/UPDD_proj}"
VM_CONDA_PREFIX=/home/san/miniconda3
PY_MD_VM="${VM_CONDA_PREFIX}/envs/md_simulation/bin/python"

CP4_REF_DIR="${PROJ}/outputs/2QKI_Cp4_calib_s7"
GAMMA_XML="${PROJ}/params/MTR_gaff2_layer3d_gamma.xml"
HYDROGENS_SRC="${CP4_REF_DIR}/params/MTR_hydrogens.xml"
MANIFEST_SRC_TEMPLATE="${CP4_REF_DIR}/params/MTR_params_manifest.json"

# V100 lane: host 의 latter-half (5th-8th seed) 를 forward order 로 처리.
# Host order: s7 → s23 → s42 → s83 → s101 → s163 → s251 → s19_reseed55
# V100 takes 5th-8th (s101/s163/s251/s19_reseed55). V100 가 항상 host 보다 빠르게
# 해당 seed 완료 → host launcher 의 skip-if-exists 가 자연 join (collision-free).
V100_SEEDS_INT=(101 163 251 55)
V100_SEEDS_LABEL=(s101 s163 s251 s19_reseed55)

# Auto-detect host launcher's logdir (latest v09_gamma_production_*)
LOGDIR=$(ls -td "${PROJ}"/outputs/analysis/v09_gamma_production_* 2>/dev/null | head -1)
if [ -z "$LOGDIR" ]; then
    LOGDIR="${PROJ}/outputs/analysis/v09_gamma_v100_lane_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$LOGDIR"
fi

V100_DISPATCH_LOG="${LOGDIR}/v100_lane_dispatch.log"

log() { printf '[v09_v100_lane %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$V100_DISPATCH_LOG"; }
err() { printf '[v09_v100_lane %s] ERROR: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$V100_DISPATCH_LOG" >&2; }

log "================================================================"
log "v0.9 γ V100 lane START (parallel to host pilot)"
log "VM target: $VM_TARGET"
log "Reverse-order seeds: ${V100_SEEDS_LABEL[*]}"
log "Logdir: $LOGDIR"
log "================================================================"

# Pre-flight VM connectivity
if ! ssh -o ConnectTimeout=5 -o StrictHostKeyChecking=no "$VM_TARGET" 'nvidia-smi --query-gpu=name --format=csv,noheader | grep -qi V100'; then
    err "VM SSH 연결 또는 V100 미감지 — abort"
    exit 2
fi
log "VM V100 connectivity OK"

# Pre-flight: VM outputs/ + reference dir 생성 (setup_vm.sh 가 outputs/ 를 exclude 했으므로
# 처음 V100 lane 실행 시 부재). reference _md_input 도 미리 sync (각 seed 처리 시 또
# sync 하지만 빠른 첫-시도 보장).
log "Pre-flight: VM 측 outputs/ + reference dir 준비"
ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR "$VM_TARGET" \
    "mkdir -p ${VM_PROJ}/outputs/$(basename ${CP4_REF_DIR})/_md_input ${VM_PROJ}/outputs/$(basename ${CP4_REF_DIR})/params" \
    2>&1 | tee -a "$V100_DISPATCH_LOG"
rsync -az --exclude=mdresult --exclude=snapshots --exclude='snapshots_*' \
    -e 'ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR' \
    "${CP4_REF_DIR}/_md_input/" "${VM_TARGET}:${VM_PROJ}/outputs/$(basename ${CP4_REF_DIR})/_md_input/" 2>&1 \
    | tail -3 >> "$V100_DISPATCH_LOG"
rsync -az \
    -e 'ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR' \
    "${CP4_REF_DIR}/params/" "${VM_TARGET}:${VM_PROJ}/outputs/$(basename ${CP4_REF_DIR})/params/" 2>&1 \
    | tail -3 >> "$V100_DISPATCH_LOG"
# 글로벌 params (γ XML 포함)
ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR "$VM_TARGET" \
    "mkdir -p ${VM_PROJ}/params" 2>&1 | tee -a "$V100_DISPATCH_LOG"
rsync -az "${PROJ}/params/MTR_gaff2_layer3d_gamma.xml" \
    -e 'ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR' \
    "${VM_TARGET}:${VM_PROJ}/params/" 2>&1 | tail -3 >> "$V100_DISPATCH_LOG"
log "Pre-flight: VM 측 준비 완료"

# Per-seed loop
for i in "${!V100_SEEDS_INT[@]}"; do
    seed_int="${V100_SEEDS_INT[$i]}"
    seed_label="${V100_SEEDS_LABEL[$i]}"
    seedir="${PROJ}/outputs/2QKI_Cp4_gamma_calib_${seed_label}"

    # Skip if already done (host launcher 가 이미 했거나 prior V100 run)
    if [ -f "${seedir}/mdresult/2QKI_Cp4_final.pdb" ]; then
        log "  $seed_label: SKIP (final.pdb exists)"
        continue
    fi

    log "  $seed_label: setup local dir"
    mkdir -p "${seedir}/_md_input" "${seedir}/params" "${seedir}/mdresult"

    # Copy reference _md_input PDBs (same as host launcher)
    for f in 2QKI_Cp4.pdb 2QKI_Cp4_renum.pdb; do
        [ -f "${CP4_REF_DIR}/_md_input/${f}" ] && \
            cp -n "${CP4_REF_DIR}/_md_input/${f}" "${seedir}/_md_input/"
    done

    cp -f "$GAMMA_XML" "${seedir}/params/MTR_gaff2.xml"
    cp -n "$HYDROGENS_SRC" "${seedir}/params/"
    cp -n "$MANIFEST_SRC_TEMPLATE" "${seedir}/params/"

    # Update manifest paths (host paths — VM 측에서 동일 path 사용)
    /home/san/miniconda3/envs/qmmm/bin/python -c "
import json
mf = '${seedir}/params/MTR_params_manifest.json'
with open(mf) as f: m = json.load(f)
m['xml_path'] = '${seedir}/params/MTR_gaff2.xml'
m['hydrogens_path'] = '${seedir}/params/MTR_hydrogens.xml'
m['charge_source'] = 'v0.9 gamma V100 lane parallel: same charge model as host γ production'
with open(mf, 'w') as f: json.dump(m, f, indent=2)
"

    log "  $seed_label: VM 측 seedir 생성 + rsync (without mdresult)"
    # mkdir on VM first (rsync 가 destination parent 부재 시 silent fail)
    ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR "$VM_TARGET" \
        "mkdir -p ${VM_PROJ}/outputs/2QKI_Cp4_gamma_calib_${seed_label}/{_md_input,params,mdresult}" \
        2>&1 | tee -a "$V100_DISPATCH_LOG"
    rsync -az --exclude=mdresult \
        -e 'ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR' \
        "${seedir}/" "${VM_TARGET}:${VM_PROJ}/outputs/2QKI_Cp4_gamma_calib_${seed_label}/" 2>&1 \
        | tail -5 >> "$V100_DISPATCH_LOG"

    # rsync rc check
    if ! ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR "$VM_TARGET" \
        "test -f ${VM_PROJ}/outputs/2QKI_Cp4_gamma_calib_${seed_label}/_md_input/2QKI_Cp4.pdb"; then
        err "  $seed_label: VM 측 _md_input PDB 부재 — rsync 실패 추정, skip"
        continue
    fi
    if ! ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR "$VM_TARGET" \
        "test -f ${VM_PROJ}/outputs/2QKI_Cp4_gamma_calib_${seed_label}/params/MTR_gaff2.xml"; then
        err "  $seed_label: VM 측 γ XML 부재 — rsync 실패 추정, skip"
        continue
    fi
    log "  $seed_label: VM 측 rsync 검증 OK"

    log "  $seed_label: V100 MD launch (5 ns, cyclic_ss, γ XML, seed=${seed_int})"

    REMOTE_LOG="${LOGDIR}/v100_md_gamma_${seed_label}.log"
    REMOTE_CMD="cd ${VM_PROJ} && \
        UPDD_NCAA_AMBER14_PATCH=0 \
        UPDD_MD_CUDA_DEVICE=0 \
        UPDD_MMGBSA_PLATFORM=CUDA \
        ${PY_MD_VM} utils/run_restrained_md.py \
            --inputdir outputs/2QKI_Cp4_gamma_calib_${seed_label}/_md_input \
            --outputdir outputs/2QKI_Cp4_gamma_calib_${seed_label}/mdresult \
            --params_manifest outputs/2QKI_Cp4_gamma_calib_${seed_label}/params/MTR_params_manifest.json \
            --steps 2500000 --topology cyclic_ss --binder_chain B --graph_policy strict \
            --target_id 2QKI --dt_fs 2.0 --platform CUDA --seed ${seed_int} \
            --ncaa_label MTR --ncaa_code MTR"

    ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR "$VM_TARGET" "$REMOTE_CMD" \
        > "$REMOTE_LOG" 2>&1
    rc=$?

    if [ "$rc" -ne 0 ]; then
        err "  $seed_label: V100 MD FAIL rc=$rc (see ${REMOTE_LOG})"
        continue
    fi

    log "  $seed_label: V100 MD complete — rsync mdresult → host"
    rsync -az --delete-after \
        -e 'ssh -o StrictHostKeyChecking=no -o LogLevel=ERROR' \
        "${VM_TARGET}:${VM_PROJ}/outputs/2QKI_Cp4_gamma_calib_${seed_label}/mdresult/" \
        "${seedir}/mdresult/" 2>&1 | tail -5 >> "$V100_DISPATCH_LOG"

    if [ -f "${seedir}/mdresult/2QKI_Cp4_final.pdb" ]; then
        log "  $seed_label: V100 lane DONE + synced to host (final.pdb present)"
    else
        err "  $seed_label: rsync 후 final.pdb 부재 (확인 필요)"
    fi
done

log "================================================================"
log "V100 lane complete. Host launcher 의 Stage A skip-if-exists 로 자연 join."
log "================================================================"
