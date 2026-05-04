#!/bin/bash
# ============================================================================
# Phase β MM-PBSA Re-run — Full Sweep (Mandatory + Optional)
# ============================================================================
# Inputs : outputs/<system>/snapshots_n25_postl387_patch/  (post L387 v54 patch)
# Outputs: outputs/<system>/mmpbsa_results_postl387/       (NEW; preserves
#          mmpbsa_results_fix4_postpatch/ per R-7).
#
# Mandatory tier (4 systems, ~4 h GPU): BROKEN→INTACT verifications.
# Optional tier  (21 systems, ~21 h GPU): Phase α regression + multi-system σ.
#
# Idempotent: skips entries whose mmpbsa_summary.json already exists.
# Pattern: scripts/track_b_cp4_postpatch.sh + track_b_1ebp_mtr13_postpatch.sh.
# ============================================================================

set -u

PY=/home/san/miniconda3/envs/qmmm/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh
conda activate qmmm

export UPDD_MMGBSA_PLATFORM=CUDA

LOG="${PROJ}/outputs/analysis/phase_beta_mmpbsa_rerun_20260429.log"
mkdir -p "$(dirname "$LOG")"

echo "==========================================" | tee -a "$LOG"
echo "PHASE β MM-PBSA RE-RUN START: $(date --iso-8601=seconds)" | tee -a "$LOG"
echo "host=$(hostname) gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"

run_pbsa () {
    local label="$1"
    local snap_dir="$2"
    local outdir="$3"
    local target_id="$4"
    local ncaa_elem="$5"

    if [ -f "${outdir}/mmpbsa_summary.json" ]; then
        echo "SKIP ${label} (already done: ${outdir}/mmpbsa_summary.json)" | tee -a "$LOG"
        return 0
    fi
    if ! ls "${snap_dir}"/*.pdb >/dev/null 2>&1; then
        echo "MISS ${label} no pdb in ${snap_dir}" | tee -a "$LOG"
        return 1
    fi
    mkdir -p "${outdir}"
    local start=$(date +%s)
    echo "[START ${label}] $(date --iso-8601=seconds) snap=${snap_dir} out=${outdir} target=${target_id} ncaa=${ncaa_elem}" | tee -a "$LOG"
    UPDD_MMGBSA_PLATFORM=CUDA $PY scripts/run_mmpbsa.py \
        --md_dir "${snap_dir}" \
        --outputdir "${outdir}" \
        --ncaa_elem "${ncaa_elem}" \
        --receptor_chain A --binder_chain B \
        --target_id "${target_id}" \
        --protocol 1traj \
        >> "$LOG" 2>&1
    local rc=$?
    local end=$(date +%s)
    local elapsed=$((end-start))
    echo "[DONE ${label}] rc=${rc} elapsed=${elapsed}s" | tee -a "$LOG"
}

# --------------------------------------------------------------------------
# Convenience wrapper: outputs/<sys>/snapshots_n25_postl387_patch/
#                   →  outputs/<sys>/mmpbsa_results_postl387/
# --------------------------------------------------------------------------
run_sys () {
    local sys="$1"
    local target_id="$2"
    local ncaa_elem="$3"
    run_pbsa "${sys}" \
             "outputs/${sys}/snapshots_n25_postl387_patch" \
             "outputs/${sys}/mmpbsa_results_postl387" \
             "${target_id}" "${ncaa_elem}"
}

# ==========================================================================
# Tier 1 — Mandatory (4 systems, BROKEN→INTACT verification)
# ==========================================================================
echo "------------------------------------------" | tee -a "$LOG"
echo "TIER 1 — Mandatory (4 systems)" | tee -a "$LOG"
echo "------------------------------------------" | tee -a "$LOG"

run_sys 7TL8_MTR6_calib_s7   7TL8 MTR
run_sys 3IOL_NML20_calib_s19 3IOL MLE
run_sys 3IOL_NML20_calib_s23 3IOL MLE
run_sys 3IOL_NML20_calib_s42 3IOL MLE

echo "[MANDATORY-COMPLETE] $(date --iso-8601=seconds)" | tee -a "$LOG"

# ==========================================================================
# Tier 2 — Optional (21 systems, Phase α regression + multi-system σ_btwn)
# ==========================================================================
echo "------------------------------------------" | tee -a "$LOG"
echo "TIER 2 — Optional (21 systems)" | tee -a "$LOG"
echo "------------------------------------------" | tee -a "$LOG"

# 1EBP_MTR13 (5 seeds)
run_sys 1EBP_MTR13_calib_s7   1EBP MTR
run_sys 1EBP_MTR13_calib_s19  1EBP MTR
run_sys 1EBP_MTR13_calib_s23  1EBP MTR
run_sys 1EBP_MTR13_calib_s42  1EBP MTR
run_sys 1EBP_MTR13_calib_s101 1EBP MTR

# 2QKI_Cp4 (8 seeds; Cp4 binder uses MTR for 4(1MeW))
run_sys 2QKI_Cp4_calib_s7            2QKI MTR
run_sys 2QKI_Cp4_calib_s19_reseed55  2QKI MTR
run_sys 2QKI_Cp4_calib_s23           2QKI MTR
run_sys 2QKI_Cp4_calib_s42           2QKI MTR
run_sys 2QKI_Cp4_calib_s83           2QKI MTR
run_sys 2QKI_Cp4_calib_s101          2QKI MTR
run_sys 2QKI_Cp4_calib_s163          2QKI MTR
run_sys 2QKI_Cp4_calib_s251          2QKI MTR

# 2QKH_MTR25 (3 seeds)
run_sys 2QKH_MTR25_calib_s7  2QKH MTR
run_sys 2QKH_MTR25_calib_s19 2QKH MTR
run_sys 2QKH_MTR25_calib_s42 2QKH MTR

# 7TL8_MTR6 (2 INTACT seeds; s7 is mandatory above)
run_sys 7TL8_MTR6_calib_s19 7TL8 MTR
run_sys 7TL8_MTR6_calib_s42 7TL8 MTR

# 1YCR_NML22 (3 INTACT seeds; s19/s42 excluded — real binder detachment)
run_sys 1YCR_NML22_calib_s7   1YCR MLE
run_sys 1YCR_NML22_calib_s23  1YCR MLE
run_sys 1YCR_NML22_calib_s101 1YCR MLE

echo "==========================================" | tee -a "$LOG"
echo "PHASE β MM-PBSA RE-RUN END: $(date --iso-8601=seconds)" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"
