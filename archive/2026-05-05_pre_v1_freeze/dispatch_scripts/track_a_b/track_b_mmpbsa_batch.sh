#!/bin/bash
# Track B: MM-PBSA batch on 10 Fix 4 trajectories (X.B Step B2)
# Rev. 5 hardening: CPU platform for GPU-sharing with Track A
# Generated 2026-04-23 by pipeline-runner

set -u

PY=/home/san/miniconda3/envs/qmmm/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

# AmberTools availability via conda env activation for MMPBSA.py subprocess
source /home/san/miniconda3/etc/profile.d/conda.sh
conda activate qmmm

# 2026-04-23 switch: CPU platform gave 70h ETA on 21k-atom complex minimize
# Track A uses <1 GB VRAM, 14.7 GB free → CUDA safe for Track B minimize (~5-10s/snap)
# Reverts to B1 validated architecture (B1 single-snap 90s/snap on CUDA)
export UPDD_MMGBSA_PLATFORM=CUDA

run_mmpbsa() {
    local variant="$1"
    local seed="$2"
    local snap_dir="$3"
    local outdir="$4"
    local ncaa_elem="$5"  # 'none' for WT, 'MTR' for Cp4

    if [ -f "${outdir}/mmpbsa_summary.json" ]; then
        echo "[SKIP 2QKI_${variant}_s${seed}] already complete: ${outdir}/mmpbsa_summary.json"
        return 0
    fi

    local start_ts=$(date +%s)
    local log="/tmp/track_b_mmpbsa_${variant}_s${seed}.log"
    echo "[START 2QKI_${variant}_s${seed}] $(date -Is) log=${log}"

    mkdir -p "${outdir}"

    $PY scripts/run_mmpbsa.py \
      --md_dir "${snap_dir}" \
      --outputdir "${outdir}" \
      --ncaa_elem "${ncaa_elem}" \
      --receptor_chain A \
      --binder_chain B \
      --target_id 2QKI \
      --protocol 1traj \
      > "${log}" 2>&1
    local rc=$?
    local end_ts=$(date +%s)
    local elapsed=$((end_ts - start_ts))
    echo "[DONE  2QKI_${variant}_s${seed}] $(date -Is) rc=${rc} elapsed=${elapsed}s"
    return $rc
}

echo "=========================================="
echo "TRACK B MM-PBSA BATCH START: $(date -Is)"
echo "Platform: CPU (UPDD_MMGBSA_PLATFORM=${UPDD_MMGBSA_PLATFORM})"
echo "=========================================="

# WT (5 reps)
run_mmpbsa WT 7   outputs/2QKI_WT_calib_s7/snapshots_n25    outputs/2QKI_WT_calib_s7/mmpbsa_results_fix4    none
run_mmpbsa WT 19  outputs/2QKI_WT_calib_s19/snapshots_n25   outputs/2QKI_WT_calib_s19/mmpbsa_results_fix4   none
run_mmpbsa WT 42  outputs/2QKI_WT_calib_s42/snapshots_n25   outputs/2QKI_WT_calib_s42/mmpbsa_results_fix4   none
run_mmpbsa WT 23  outputs/2QKI_WT_calib_s23/snapshots       outputs/2QKI_WT_calib_s23/mmpbsa_results_fix4   none
run_mmpbsa WT 101 outputs/2QKI_WT_calib_s101/snapshots      outputs/2QKI_WT_calib_s101/mmpbsa_results_fix4  none

# Cp4 (5 reps, incl s19_reseed55)
run_mmpbsa Cp4 7    outputs/2QKI_Cp4_calib_s7/snapshots_n25          outputs/2QKI_Cp4_calib_s7/mmpbsa_results_fix4           MTR
run_mmpbsa Cp4 19r  outputs/2QKI_Cp4_calib_s19_reseed55/snapshots    outputs/2QKI_Cp4_calib_s19_reseed55/mmpbsa_results_fix4 MTR
run_mmpbsa Cp4 42   outputs/2QKI_Cp4_calib_s42/snapshots_n25         outputs/2QKI_Cp4_calib_s42/mmpbsa_results_fix4          MTR
run_mmpbsa Cp4 23   outputs/2QKI_Cp4_calib_s23/snapshots             outputs/2QKI_Cp4_calib_s23/mmpbsa_results_fix4          MTR
run_mmpbsa Cp4 101  outputs/2QKI_Cp4_calib_s101/snapshots            outputs/2QKI_Cp4_calib_s101/mmpbsa_results_fix4         MTR

echo "=========================================="
echo "TRACK B MM-PBSA BATCH END:   $(date -Is)"
echo "=========================================="
