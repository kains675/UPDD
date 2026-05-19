#!/bin/bash
# ============================================================================
# T1 N→30 Sampling Expansion — Phase 1 orchestrator (2026-05-18)
# ============================================================================
# Thin wrapper around three existing helpers (no new pipeline logic):
#   1. scripts/phase_beta_repbsa_v2_reextract.py — snapshot re-extraction
#      under CONECT v55, parameterised on --out-subdir for Case A vs B.
#   2. scripts/phase_beta_repbsa_v2.sh — MMPBSA driver (idempotent, skip
#      on existing summary).
#   3. scripts/phase_beta_repbsa_v2_aggregate.py — family-level σ_btwn /
#      SE / z_SE / Tier-1/2/3 classification.
#
# Phase 1 work units (per plan/t1_n30_expansion_plan_20260518.md):
#   1A. Cp4 (2QKI) postl387_v2 — 7 MD-ready seeds → snap extract + MMPBSA
#   1B. 2QKH_MTR25 postl387 v1 — 9 in-flight restart + 3 MD-ready snap+PBSA
#   1C. 2QKI_WT v2 — 2 pilot seeds (s7, s23) for protocol validation
#
# Pre-launch conditional approval conditions C1–C4 satisfied per plan §0.
# Phase 2 (Cp4 fresh MD × 15) + Phase 3 (WT full N=30) DEFERRED — gated on
# Phase 1 aggregate results + user re-approval.
# ============================================================================

set -u

PY=/home/san/miniconda3/envs/qmmm/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh
conda activate qmmm

# Production env. Cohort consistency note (2026-05-18):
#
# Cp4 baseline n=8 published in Paper 1 v2 §3.3 (DOI 10.26434/chemrxiv.15002948/v2)
# was generated with MD-time `params/MTR_gaff2.xml` carrying GAFF2-unpatched
# charges (q_N = -0.8938, type=n8) — confirmed by per-seed `mdresult/*_universal.xml`
# audit. To preserve cohort integrity for the N→30 expansion, this orchestrator
# launches new MD with the amber14SB patch explicitly DISABLED so the MD sampling
# landscape matches the published baseline. PBSA scoring still reads the current
# `params/MTR_gaff2.xml` (patched), reproducing the Paper §3.5 "patch ON at scoring
# time" framing without re-perturbing the MD ensemble. 2QKH_MTR25 cohort is
# homogeneously patched (all 27 seeds q_N = -0.4157) so this OFF default is
# Cp4-specific; MTR25 already-prepared seeds are unaffected.
export UPDD_MTR_AMBER14_PATCH=0
export UPDD_MMGBSA_PLATFORM=CUDA

LAUNCH_TS=$(date +%s)
TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/t1_n30_expansion_${TS_TAG}"
mkdir -p "$LOGDIR"
DISPATCH_LOG="${LOGDIR}/dispatch.log"

# Hand the same log path to phase_beta_repbsa_v2.sh via T1_LOG env var (the
# v2 script honors T1_LOG override when set; falls back to historical path).
export T1_LOG="${LOGDIR}/mmpbsa.log"

# Stage status JSONs (per-stage state for restart / aggregation)
REEXTRACT_V2_JSON="${LOGDIR}/reextract_v2_status.json"
REEXTRACT_V1_JSON="${LOGDIR}/reextract_v1_status.json"

log () {
    local msg="$1"
    printf '[t1_n30_expansion %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$msg" \
        | tee -a "$DISPATCH_LOG"
}

# ----------------------------------------------------------------------------
# Pre-flight checks
# ----------------------------------------------------------------------------
log "================================================================"
log "T1 N→30 expansion — Phase 1 START"
log "host=$(hostname)  gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
log "logdir=$LOGDIR"
log "================================================================"

# GPU exclusive — bail if another user/compute app holds the GPU.
EXT_APPS=$(nvidia-smi --query-compute-apps=pid --format=csv,noheader 2>/dev/null \
           | grep -v '^$' | wc -l)
SELF_PGID=$(ps -o pgid= -p $$ | tr -d ' ')
if [ "$EXT_APPS" -gt 0 ]; then
    # Tolerate self-launched MMPBSA children; check if any external PIDs hold GPU.
    EXTERNAL=0
    while read -r pid; do
        [ -z "$pid" ] && continue
        pid_pgid=$(ps -o pgid= -p "$pid" 2>/dev/null | tr -d ' ')
        if [ -n "$pid_pgid" ] && [ "$pid_pgid" != "$SELF_PGID" ]; then
            EXTERNAL=$((EXTERNAL + 1))
        fi
    done < <(nvidia-smi --query-compute-apps=pid --format=csv,noheader)
    if [ "$EXTERNAL" -gt 0 ]; then
        log "ABORT: GPU has $EXTERNAL external compute apps. Free GPU first."
        exit 2
    fi
fi
log "GPU exclusive check: PASS"

# Thermal sanity (post thermal-paste refresh — should be cool at idle)
T_GPU=$(nvidia-smi --query-gpu=temperature.gpu --format=csv,noheader,nounits | head -1)
log "GPU idle temperature: ${T_GPU}°C"
if [ "$T_GPU" -gt 60 ]; then
    log "WARN: idle T_GPU=${T_GPU}°C > 60°C — investigate cooling"
fi

# ----------------------------------------------------------------------------
# Stage 1: snapshot re-extraction
# ----------------------------------------------------------------------------
log "----------------------------------------------------------------"
log "Stage 1a: snapshot reextract (Cp4 + 2QKI_WT) → v2 path"
log "----------------------------------------------------------------"
$PY scripts/phase_beta_repbsa_v2_reextract.py \
    --filter 2QKI_ \
    --out-json "$REEXTRACT_V2_JSON" \
    >> "$DISPATCH_LOG" 2>&1
RC_V2=$?
log "Stage 1a exit code: $RC_V2  (status JSON: $REEXTRACT_V2_JSON)"
if [ "$RC_V2" -ne 0 ]; then
    log "WARN: Stage 1a non-zero exit — review $REEXTRACT_V2_JSON before proceeding"
fi

log "----------------------------------------------------------------"
log "Stage 1b: snapshot reextract (2QKH_MTR25) → v1 path (Case A)"
log "----------------------------------------------------------------"
$PY scripts/phase_beta_repbsa_v2_reextract.py \
    --filter 2QKH_MTR25 \
    --out-subdir snapshots_n25_postl387_patch \
    --out-json "$REEXTRACT_V1_JSON" \
    >> "$DISPATCH_LOG" 2>&1
RC_V1=$?
log "Stage 1b exit code: $RC_V1  (status JSON: $REEXTRACT_V1_JSON)"
if [ "$RC_V1" -ne 0 ]; then
    log "WARN: Stage 1b non-zero exit — review $REEXTRACT_V1_JSON before proceeding"
fi

# ----------------------------------------------------------------------------
# Stage 2: MMPBSA — N-lane parallel (T1_LANES env, default 4), idempotent
# ----------------------------------------------------------------------------
# CPU profiling (2026-05-18) showed AmberTools `sander` is single-thread
# CPU-bound (~1 core at 99%). The 9800X3D 8C/16T has 15 cores idle when the
# sequential `phase_beta_repbsa_v2.sh` driver runs. We replace the Stage 2
# sequential call with an inline 4-lane parallel dispatcher: each lane spawns
# its own `run_mmpbsa.py` → independent `sander` process on a separate core.
# updd_config.py declares `mmgbsa=4` as the canonical concurrency for this
# workload. RAM cost per lane ≈ 2% (~640 MB) → 4 lanes ≈ 2.5 GB < 32 GB free.
# Each lane writes its own seed log under $LOGDIR/mmpbsa_<label>.log; the
# aggregated dispatch line stays in $T1_LOG. GPU unused by sander.
#
# Sequential reference: scripts/phase_beta_repbsa_v2.sh (preserved for
# single-lane runs and as the historical Phase β #84 record).
LANES="${T1_LANES:-4}"

log "----------------------------------------------------------------"
log "Stage 2: MMPBSA parallel sweep (${LANES} lanes, idempotent)"
log "  per-lane log dir: $LOGDIR"
log "  dispatch log:     $T1_LOG"
log "----------------------------------------------------------------"

run_one_pbsa () {
    local label="$1"
    local snap_dir="$2"
    local out_dir="$3"
    local target_id="$4"
    local ncaa_elem="$5"

    if [ -f "${out_dir}/mmpbsa_summary.json" ]; then
        printf '[%s SKIP %s]\n' "$(date +%H:%M:%S)" "$label" >> "$T1_LOG"
        return 0
    fi
    if ! ls "${snap_dir}"/*.pdb >/dev/null 2>&1; then
        printf '[%s MISS %s] no pdb in %s\n' "$(date +%H:%M:%S)" "$label" "$snap_dir" >> "$T1_LOG"
        return 1
    fi
    mkdir -p "$out_dir"
    local seed_log="${LOGDIR}/mmpbsa_${label}.log"
    local start; start=$(date +%s)
    printf '[%s START %s]\n' "$(date +%H:%M:%S)" "$label" >> "$T1_LOG"
    UPDD_MMGBSA_PLATFORM=CUDA "$PY" scripts/run_mmpbsa.py \
        --md_dir "$snap_dir" --outputdir "$out_dir" \
        --ncaa_elem "$ncaa_elem" --receptor_chain A --binder_chain B \
        --target_id "$target_id" --protocol 1traj \
        > "$seed_log" 2>&1
    local rc=$?
    local end; end=$(date +%s)
    printf '[%s DONE  %s] rc=%d elapsed=%ds\n' "$(date +%H:%M:%S)" "$label" "$rc" "$((end-start))" >> "$T1_LOG"
}

# Work units: "label|snap_dir|out_dir|target_id|ncaa_elem"
# Order: SKIP-fast items first (baseline DONE → 0s each), then new work.
WORK_UNITS=(
    # Cp4 baseline 8 (ALREADY_DONE, fast skip)
    "2QKI_Cp4_s7|outputs/2QKI_Cp4_calib_s7/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s7/mmpbsa_results_postl387_v2|2QKI|MTR"
    "2QKI_Cp4_s19_reseed55|outputs/2QKI_Cp4_calib_s19_reseed55/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s19_reseed55/mmpbsa_results_postl387_v2|2QKI|MTR"
    "2QKI_Cp4_s23|outputs/2QKI_Cp4_calib_s23/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s23/mmpbsa_results_postl387_v2|2QKI|MTR"
    "2QKI_Cp4_s42|outputs/2QKI_Cp4_calib_s42/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s42/mmpbsa_results_postl387_v2|2QKI|MTR"
    "2QKI_Cp4_s83|outputs/2QKI_Cp4_calib_s83/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s83/mmpbsa_results_postl387_v2|2QKI|MTR"
    "2QKI_Cp4_s101|outputs/2QKI_Cp4_calib_s101/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s101/mmpbsa_results_postl387_v2|2QKI|MTR"
    "2QKI_Cp4_s163|outputs/2QKI_Cp4_calib_s163/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s163/mmpbsa_results_postl387_v2|2QKI|MTR"
    "2QKI_Cp4_s251|outputs/2QKI_Cp4_calib_s251/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s251/mmpbsa_results_postl387_v2|2QKI|MTR"
    # 7TL8 2 (ALREADY_DONE)
    "7TL8_MTR6_s19|outputs/7TL8_MTR6_calib_s19/snapshots_n25_postl387_patch_v2|outputs/7TL8_MTR6_calib_s19/mmpbsa_results_postl387_v2|7TL8|MTR"
    "7TL8_MTR6_s42|outputs/7TL8_MTR6_calib_s42/snapshots_n25_postl387_patch_v2|outputs/7TL8_MTR6_calib_s42/mmpbsa_results_postl387_v2|7TL8|MTR"
    # T1 Phase 1A — Cp4 GAFF2 new 3 (active work)
    "2QKI_Cp4_s127|outputs/2QKI_Cp4_calib_s127/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s127/mmpbsa_results_postl387_v2|2QKI|MTR"
    "2QKI_Cp4_s19|outputs/2QKI_Cp4_calib_s19/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s19/mmpbsa_results_postl387_v2|2QKI|MTR"
    "2QKI_Cp4_s199|outputs/2QKI_Cp4_calib_s199/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_calib_s199/mmpbsa_results_postl387_v2|2QKI|MTR"
    # T1 Phase 1B — 2QKH_MTR25 in-flight restart 9 (active work, v1 path)
    "2QKH_MTR25_s317|outputs/2QKH_MTR25_calib_s317/snapshots_n25_postl387_patch|outputs/2QKH_MTR25_calib_s317/mmpbsa_results_postl387|2QKH|MTR"
    "2QKH_MTR25_s367|outputs/2QKH_MTR25_calib_s367/snapshots_n25_postl387_patch|outputs/2QKH_MTR25_calib_s367/mmpbsa_results_postl387|2QKH|MTR"
    "2QKH_MTR25_s421|outputs/2QKH_MTR25_calib_s421/snapshots_n25_postl387_patch|outputs/2QKH_MTR25_calib_s421/mmpbsa_results_postl387|2QKH|MTR"
    "2QKH_MTR25_s827|outputs/2QKH_MTR25_calib_s827/snapshots_n25_postl387_patch|outputs/2QKH_MTR25_calib_s827/mmpbsa_results_postl387|2QKH|MTR"
    "2QKH_MTR25_s881|outputs/2QKH_MTR25_calib_s881/snapshots_n25_postl387_patch|outputs/2QKH_MTR25_calib_s881/mmpbsa_results_postl387|2QKH|MTR"
    "2QKH_MTR25_s937|outputs/2QKH_MTR25_calib_s937/snapshots_n25_postl387_patch|outputs/2QKH_MTR25_calib_s937/mmpbsa_results_postl387|2QKH|MTR"
    "2QKH_MTR25_s1361|outputs/2QKH_MTR25_calib_s1361/snapshots_n25_postl387_patch|outputs/2QKH_MTR25_calib_s1361/mmpbsa_results_postl387|2QKH|MTR"
    "2QKH_MTR25_s1423|outputs/2QKH_MTR25_calib_s1423/snapshots_n25_postl387_patch|outputs/2QKH_MTR25_calib_s1423/mmpbsa_results_postl387|2QKH|MTR"
    "2QKH_MTR25_s1483|outputs/2QKH_MTR25_calib_s1483/snapshots_n25_postl387_patch|outputs/2QKH_MTR25_calib_s1483/mmpbsa_results_postl387|2QKH|MTR"
    # T1 Phase 1C — 2QKI_WT pilot 2 (active work)
    "2QKI_WT_s7|outputs/2QKI_WT_calib_s7/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s7/mmpbsa_results_postl387_v2|2QKI|none"
    "2QKI_WT_s23|outputs/2QKI_WT_calib_s23/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s23/mmpbsa_results_postl387_v2|2QKI|none"
)

# N-lane dispatch via bash background + wait -n throttle.
running=0
for entry in "${WORK_UNITS[@]}"; do
    IFS='|' read -r label snap out tgt ncaa <<< "$entry"
    run_one_pbsa "$label" "$snap" "$out" "$tgt" "$ncaa" &
    running=$((running + 1))
    if [ "$running" -ge "$LANES" ]; then
        wait -n
        running=$((running - 1))
    fi
done
wait
RC_MM=$?
log "Stage 2 exit code: $RC_MM  (parallel ${LANES}-lane sweep complete)"

# ----------------------------------------------------------------------------
# Stage 3: aggregate
# ----------------------------------------------------------------------------
log "----------------------------------------------------------------"
log "Stage 3: aggregate σ_btwn / SE / z_SE / Tier classification"
log "----------------------------------------------------------------"
$PY scripts/phase_beta_repbsa_v2_aggregate.py \
    >> "$DISPATCH_LOG" 2>&1
RC_AGG=$?
log "Stage 3 exit code: $RC_AGG"

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
ELAPSED=$(( $(date +%s) - LAUNCH_TS ))
ELAPSED_H=$(echo "scale=2; $ELAPSED / 3600" | bc)
log "================================================================"
log "T1 N→30 expansion — Phase 1 COMPLETE  elapsed=${ELAPSED}s (~${ELAPSED_H}h)"
log "  Stage 1a (Cp4 + WT v2 reextract):   rc=$RC_V2"
log "  Stage 1b (MTR25 v1 reextract):      rc=$RC_V1"
log "  Stage 2  (MMPBSA sweep):            rc=$RC_MM"
log "  Stage 3  (aggregate):               rc=$RC_AGG"
log ""
log "Next: review outputs/analysis/phase_beta_repbsa_v2_aggregate.json,"
log "      then decide Phase 2 (Cp4 fresh MD × 15) per plan §8."
log "================================================================"

# Exit code: 0 if all stages succeeded, else non-zero summary
if [ "$RC_V2" -eq 0 ] && [ "$RC_V1" -eq 0 ] && [ "$RC_MM" -eq 0 ] && [ "$RC_AGG" -eq 0 ]; then
    exit 0
else
    exit 1
fi
