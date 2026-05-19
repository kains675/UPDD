#!/bin/bash
# ============================================================================
# T1 Phase 1.5 Stage D — CUDA mode rescue (2026-05-19)
# ----------------------------------------------------------------------------
# CPU 1traj minimization 이 OpenMM CPU platform 에서 stuck (3h+ 0 progress)
# → MD Stage A 완료로 GPU 자유 → CUDA platform 으로 Stage D 재실행.
#
# 처리 대상: Cp4 fresh 2 (s53, s89) + WT 8 = 10 work units.
#            Cp4 s113/s173 는 reextract LOAD_FAILED (snapshot 0개) — MISS.
#
# 예상 wall: 2-lane × ~30분 per seed × 5 batches = ~150min 또는
#            (snapshot 당 ~30-60s CUDA) × 25 snap × 5 batches = ~60-90min
# ============================================================================

set -u
PY=/home/san/miniconda3/envs/qmmm/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh
conda activate qmmm

# CUDA mode — MD 끝났으니 GPU 사용 가능
unset CUDA_VISIBLE_DEVICES
export UPDD_MMGBSA_PLATFORM=CUDA

TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/t1_phase1_5_stage_d_cuda_${TS_TAG}"
mkdir -p "$LOGDIR"
LOG="${LOGDIR}/stage_d_cuda.log"

log () {
    printf '[stage_d_cuda %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG"
}

LANES=2

log "================================================================"
log "Stage D CUDA rescue START"
log "logdir=$LOGDIR  lanes=$LANES  platform=CUDA"
log "================================================================"

run_one_pbsa () {
    local label="$1" snap_dir="$2" out_dir="$3" target_id="$4" ncaa_elem="$5"
    if [ -f "${out_dir}/mmpbsa_summary.json" ]; then
        log "  $label: SKIP (summary exists)"
        return 0
    fi
    if ! ls "${snap_dir}"/*.pdb >/dev/null 2>&1; then
        log "  $label: MISS no pdb"
        return 1
    fi
    log "  $label: START"
    local START=$(date +%s)
    "$PY" -u scripts/run_mmpbsa.py \
        --md_dir "$snap_dir" --outputdir "$out_dir" \
        --ncaa_elem "$ncaa_elem" --receptor_chain A --binder_chain B \
        --target_id "$target_id" --protocol 1traj \
        > "${LOGDIR}/mmpbsa_${label}.log" 2>&1
    local RC=$?
    local END=$(date +%s)
    log "  $label: DONE rc=$RC elapsed=$((END-START))s"
}

WORK=(
    "Cp4_s53|outputs/2QKI_Cp4_calib_s53/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s53/mmpbsa_results_postl387_v2|2QKI|MTR"
    "Cp4_s89|outputs/2QKI_Cp4_calib_s89/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s89/mmpbsa_results_postl387_v2|2QKI|MTR"
    "WT_s19|outputs/2QKI_WT_calib_s19/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s19/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s42|outputs/2QKI_WT_calib_s42/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s42/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s83|outputs/2QKI_WT_calib_s83/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s83/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s101|outputs/2QKI_WT_calib_s101/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s101/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s127|outputs/2QKI_WT_calib_s127/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s127/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s163|outputs/2QKI_WT_calib_s163/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s163/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s199|outputs/2QKI_WT_calib_s199/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s199/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s251|outputs/2QKI_WT_calib_s251/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s251/mmpbsa_results_postl387_v2|2QKI|none"
)

log "----------------------------------------------------------------"
log "Phase D: ${#WORK[@]} work units, ${LANES} lanes"
log "----------------------------------------------------------------"

running=0
for entry in "${WORK[@]}"; do
    IFS='|' read -r label snap out tgt ncaa <<< "$entry"
    run_one_pbsa "$label" "$snap" "$out" "$tgt" "$ncaa" &
    running=$((running + 1))
    if [ "$running" -ge "$LANES" ]; then
        wait -n
        running=$((running - 1))
    fi
done
wait

log "Phase D complete."

# ============================================================================
# Phase E — Branched ΔΔG aggregation
# ============================================================================
log "----------------------------------------------------------------"
log "Phase E: branched_ddg aggregation"
log "----------------------------------------------------------------"
"$PY" -m utils.branched_ddg \
    --wt "outputs/2QKI_WT_calib_*/mmpbsa_results_postl387_v2" \
    --variant "outputs/2QKI_Cp4_calib_*/mmpbsa_results_postl387_v2" \
    --output "${LOGDIR}/branched_ddg/" \
    >> "$LOG" 2>&1 && log "Phase E DONE" || log "Phase E WARN non-zero exit"

# Summary
log "----------------------------------------------------------------"
log "Stage D CUDA rescue COMPLETE"
log "Final state:"
for entry in "${WORK[@]}"; do
    IFS='|' read -r label snap out tgt ncaa <<< "$entry"
    if [ -f "${out}/mmpbsa_summary.json" ]; then
        log "  $label: ✓ done"
    else
        log "  $label: ✗ failed"
    fi
done
log "branched_ddg.json: ${LOGDIR}/branched_ddg/branched_ddg.json"
log "================================================================"
