#!/bin/bash
# ============================================================================
# T1 Phase 1.5 Stage D CUDA — Supplementary Lanes (Plan C, 2026-05-19)
# ----------------------------------------------------------------------------
# 기존 stage_d_cuda.sh (PID 472540, 2-lane: WT_s19/s42/s83/s101) 와 병행 진행.
# 이 launcher 는 추가 2 lane 으로 WT_s127/s163/s199/s251 처리.
# 총 효과: 4 lanes 동시 진행 → wall time 약 절반 단축.
#
# 손실: 0 (기존 lane 진행 중 작업 보존)
# Final Stage E aggregation: 본 script 가 자기 work 끝난 후 호출
#   (기존 launcher 도 자기 끝날 때 Stage E 호출하지만 결과 partial — 무시)
# ============================================================================

set -u
PY=/home/san/miniconda3/envs/qmmm/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh
conda activate qmmm

unset CUDA_VISIBLE_DEVICES
export UPDD_MMGBSA_PLATFORM=CUDA

TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/t1_phase1_5_stage_d_cuda_supp_${TS_TAG}"
mkdir -p "$LOGDIR"
LOG="${LOGDIR}/stage_d_cuda_supp.log"

log () {
    printf '[supp %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG"
}

LANES=2

log "================================================================"
log "Stage D CUDA Supplementary lanes START (Plan C)"
log "logdir=$LOGDIR  lanes=$LANES (기존 2 lane + 본 2 lane = 총 4 lanes)"
log "분담: WT_s127, WT_s163, WT_s199, WT_s251"
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
    "WT_s127|outputs/2QKI_WT_calib_s127/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s127/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s163|outputs/2QKI_WT_calib_s163/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s163/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s199|outputs/2QKI_WT_calib_s199/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s199/mmpbsa_results_postl387_v2|2QKI|none"
    "WT_s251|outputs/2QKI_WT_calib_s251/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s251/mmpbsa_results_postl387_v2|2QKI|none"
)

log "----------------------------------------------------------------"
log "Phase D supp: ${#WORK[@]} work units, ${LANES} lanes"
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

log "Phase D supp complete (4 seeds processed)."

# Wait briefly to ensure 기존 launcher 도 자기 work 끝나도록 (timing margin)
log "Waiting 60s for main launcher (PID 472540) to finalize before Stage E..."
sleep 60

# Check main launcher status
if pgrep -af "t1_phase1_5_stage_d_cuda.sh" | grep -v supp > /dev/null 2>&1; then
    log "  Main launcher still running — waiting up to 30 min more..."
    for i in 1 2 3 4 5 6; do
        sleep 300  # 5min intervals, max 30min wait
        if ! pgrep -af "t1_phase1_5_stage_d_cuda.sh" | grep -v supp > /dev/null 2>&1; then
            log "  Main launcher finished after $((i*5)) min wait."
            break
        fi
    done
fi

# ============================================================================
# FINAL Phase E — Branched ΔΔG aggregation (모든 lane 완료 후)
# ============================================================================
log "----------------------------------------------------------------"
log "FINAL Phase E: branched_ddg aggregation (모든 10 seeds)"
log "----------------------------------------------------------------"

# Verify all expected done
DONE_COUNT=0
for sys in 2QKI_Cp4_calib_s53 2QKI_Cp4_calib_s89 2QKI_WT_calib_s19 2QKI_WT_calib_s42 2QKI_WT_calib_s83 2QKI_WT_calib_s101 2QKI_WT_calib_s127 2QKI_WT_calib_s163 2QKI_WT_calib_s199 2QKI_WT_calib_s251; do
    [ -f "outputs/$sys/mmpbsa_results_postl387_v2/mmpbsa_summary.json" ] && DONE_COUNT=$((DONE_COUNT+1))
done
log "  pre-aggregation done count: $DONE_COUNT / 10"

"$PY" -m utils.branched_ddg \
    --wt "outputs/2QKI_WT_calib_*/mmpbsa_results_postl387_v2" \
    --variant "outputs/2QKI_Cp4_calib_*/mmpbsa_results_postl387_v2" \
    --output "${LOGDIR}/branched_ddg/" \
    >> "$LOG" 2>&1 && log "FINAL Phase E DONE — branched_ddg.json at ${LOGDIR}/branched_ddg/" \
                   || log "FINAL Phase E WARN non-zero exit"

# Summary
log "================================================================"
log "Stage D CUDA Plan C 통합 완료"
log "Final state (10 seeds):"
for entry in \
    "Cp4_s53|outputs/2QKI_Cp4_calib_s53" \
    "Cp4_s89|outputs/2QKI_Cp4_calib_s89" \
    "WT_s19|outputs/2QKI_WT_calib_s19" \
    "WT_s42|outputs/2QKI_WT_calib_s42" \
    "WT_s83|outputs/2QKI_WT_calib_s83" \
    "WT_s101|outputs/2QKI_WT_calib_s101" \
    "WT_s127|outputs/2QKI_WT_calib_s127" \
    "WT_s163|outputs/2QKI_WT_calib_s163" \
    "WT_s199|outputs/2QKI_WT_calib_s199" \
    "WT_s251|outputs/2QKI_WT_calib_s251"; do
    IFS='|' read -r label outdir <<< "$entry"
    if [ -f "${outdir}/mmpbsa_results_postl387_v2/mmpbsa_summary.json" ]; then
        log "  $label: ✓ done"
    else
        log "  $label: ✗ failed"
    fi
done
log "FINAL branched_ddg.json: ${LOGDIR}/branched_ddg/branched_ddg.json"
log "================================================================"
