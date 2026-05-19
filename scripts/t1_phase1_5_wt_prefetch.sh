#!/bin/bash
# ============================================================================
# T1 Phase 1.5 — WT Prefetch (parallel to Stage A Cp4 MD on GPU)
# ============================================================================
# Purpose: Run WT 8 seeds' reextract + MMPBSA postl387_v2 on CPU while the
# main t1_phase1_5_fresh.sh script is doing Cp4 fresh MD on GPU. The main
# script's Stage C (reextract) and Stage D (MMPBSA) are idempotent — if
# snapshots_n25_postl387_patch_v2/ has >= 25 PDBs (Stage C SKIP), and if
# mmpbsa_summary.json exists (Stage D SKIP), the main script auto-skips the
# WT portion when it gets there.
#
# Resource separation:
#   - Reextract: mdtraj+numpy CPU only (no openmm Simulation, no GPU)
#   - MMPBSA: UPDD_MMGBSA_PLATFORM=CPU forced (OpenMM CPU platform, no CUDA)
#   - 4-lane parallel CPU; MD on GPU stays uncontested
#
# Pre-conditions:
#   - WT 8 seeds (s19, s42, s83, s101, s127, s163, s199, s251) have
#     mdresult/2QKI_WT_restrained.dcd in place
#   - cyclic_ss baseline cohort (verified by s127/s199 SG-SG std<0.05)
#
# Estimated wall: reextract 4-lane ~20-30 min + MMPBSA 4-lane ~25-35 min
# = TOTAL ~45-65 min
# ============================================================================

set -u
PY=/home/san/miniconda3/envs/qmmm/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh

# ============================================================================
# CPU affinity 자동 분배 (main script 와 동일 logic — 2026-05-19 hardening)
# ----------------------------------------------------------------------------
# 8 phys core 시스템에서 main MD (core 0) 와 충돌 회피: WT prefetch 는 cores 1-7.
# 8 phys 초과 → free scheduling fallback.
# ============================================================================
PHYS_CORES=$(lscpu -p 2>/dev/null | grep -v "^#" | awk -F, '{print $2}' | sort -nu | wc -l)
LOGICAL_CORES=$(nproc)
if [ "${PHYS_CORES:-0}" -ge 4 ] && [ "${PHYS_CORES:-0}" -le 8 ]; then
    MMPBSA_AFFINITY="1-$((PHYS_CORES - 1)),$((PHYS_CORES + 1))-$((LOGICAL_CORES - 1))"
    MMPBSA_PREFIX="taskset -c ${MMPBSA_AFFINITY}"
    USE_AFFINITY=1
else
    MMPBSA_AFFINITY=""
    MMPBSA_PREFIX=""
    USE_AFFINITY=0
fi

WT_SEEDS=(s19 s42 s83 s101 s127 s163 s199 s251)
TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/t1_phase1_5_wt_prefetch_${TS_TAG}"
mkdir -p "$LOGDIR"
LOG="${LOGDIR}/prefetch.log"

log () {
    printf '[wt_prefetch %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" | tee -a "$LOG"
}

log "================================================================"
log "T1 Phase 1.5 WT Prefetch START"
log "host=$(hostname)"
log "logdir=$LOGDIR"
log "WT seeds: ${WT_SEEDS[*]}"
log "phys_cores=$PHYS_CORES  logical=$LOGICAL_CORES  affinity=${USE_AFFINITY}"
if [ "$USE_AFFINITY" -eq 1 ]; then
    log "  MMPBSA affinity: ${MMPBSA_AFFINITY} (avoids core 0 reserved for MD)"
else
    log "  affinity off (free scheduling)"
fi
log "================================================================"

# ============================================================================
# Phase 1 — Reextract (CPU mdtraj, 4-lane parallel)
# ============================================================================
log "----------------------------------------------------------------"
log "Phase 1: reextract 8 WT seeds (4-lane CPU)"
log "----------------------------------------------------------------"

conda activate qmmm

LANES=4
running=0
for seed in "${WT_SEEDS[@]}"; do
    SYS="2QKI_WT_calib_${seed}"
    TARGET="${PROJ}/outputs/${SYS}/snapshots_n25_postl387_patch_v2"
    if [ -d "$TARGET" ] && [ "$(ls $TARGET/*.pdb 2>/dev/null | wc -l)" -ge 25 ]; then
        log "  $seed: reextract SKIP (already has >=25 PDBs)"
        continue
    fi
    (
        log "  $seed: reextract START  affinity=[${MMPBSA_AFFINITY:-free}]"
        CUDA_VISIBLE_DEVICES="" $MMPBSA_PREFIX $PY scripts/phase_beta_repbsa_v2_reextract.py \
            --filter "${SYS}" \
            --out-json "${LOGDIR}/reextract_${SYS}.json" \
            > "${LOGDIR}/reextract_${SYS}.log" 2>&1
        log "  $seed: reextract DONE rc=$?"
    ) &
    running=$((running + 1))
    if [ "$running" -ge "$LANES" ]; then
        wait -n
        running=$((running - 1))
    fi
done
wait

log "Phase 1 (reextract) complete."

# ============================================================================
# Phase 2 — MMPBSA postl387_v2 (CPU forced, 4-lane parallel)
# ============================================================================
log "----------------------------------------------------------------"
log "Phase 2: MMPBSA 8 WT seeds (4-lane CPU forced)"
log "----------------------------------------------------------------"

run_one_pbsa () {
    local seed="$1"
    local SYS="2QKI_WT_calib_${seed}"
    local SNAP="${PROJ}/outputs/${SYS}/snapshots_n25_postl387_patch_v2"
    local OUT="${PROJ}/outputs/${SYS}/mmpbsa_results_postl387_v2"

    if [ -f "${OUT}/mmpbsa_summary.json" ]; then
        log "  $seed: MMPBSA SKIP (summary already exists)"
        return 0
    fi
    if ! ls "${SNAP}"/*.pdb >/dev/null 2>&1; then
        log "  $seed: MMPBSA MISS (no snapshots)"
        return 1
    fi
    mkdir -p "$OUT"
    local START
    START=$(date +%s)
    log "  $seed: MMPBSA START  affinity=[${MMPBSA_AFFINITY:-free}]"
    CUDA_VISIBLE_DEVICES="" UPDD_MMGBSA_PLATFORM=CPU $MMPBSA_PREFIX "$PY" -u scripts/run_mmpbsa.py \
        --md_dir "$SNAP" \
        --outputdir "$OUT" \
        --ncaa_elem none \
        --receptor_chain A \
        --binder_chain B \
        --target_id 2QKI \
        --protocol 1traj \
        > "${LOGDIR}/mmpbsa_${SYS}.log" 2>&1
    local RC=$?
    local END
    END=$(date +%s)
    log "  $seed: MMPBSA DONE rc=$RC elapsed=$((END-START))s"
    return $RC
}

LANES=4
running=0
for seed in "${WT_SEEDS[@]}"; do
    run_one_pbsa "$seed" &
    running=$((running + 1))
    if [ "$running" -ge "$LANES" ]; then
        wait -n
        running=$((running - 1))
    fi
done
wait

log "Phase 2 (MMPBSA) complete."

# ============================================================================
# Summary
# ============================================================================
log "----------------------------------------------------------------"
log "WT Prefetch COMPLETE"
log "Reextract logs: ${LOGDIR}/reextract_*.log"
log "MMPBSA logs: ${LOGDIR}/mmpbsa_*.log"
log "----------------------------------------------------------------"

# Print per-seed completion summary
log "Final state:"
for seed in "${WT_SEEDS[@]}"; do
    SYS="2QKI_WT_calib_${seed}"
    SNAP="${PROJ}/outputs/${SYS}/snapshots_n25_postl387_patch_v2"
    OUT="${PROJ}/outputs/${SYS}/mmpbsa_results_postl387_v2"
    N_PDB=$(ls $SNAP/*.pdb 2>/dev/null | wc -l)
    if [ -f "${OUT}/mmpbsa_summary.json" ]; then
        log "  $seed: snapshots=$N_PDB  mmpbsa=DONE"
    else
        log "  $seed: snapshots=$N_PDB  mmpbsa=MISS"
    fi
done
