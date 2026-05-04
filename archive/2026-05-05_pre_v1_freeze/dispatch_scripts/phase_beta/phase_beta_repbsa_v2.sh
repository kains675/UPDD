#!/bin/bash
# ============================================================================
# Phase β re-PBSA v2 — 10 systems with CONECT v55 patch
# ============================================================================
# Inputs : outputs/<sys>/snapshots_n25_postl387_patch_v2/  (CONECT v55)
# Outputs: outputs/<sys>/mmpbsa_results_postl387_v2/       (NEW; preserves
#          mmpbsa_results_postl387/ from #84 per R-7).
#
# Idempotent: skips entries whose mmpbsa_summary.json already exists.
# Pattern: mirrors scripts/phase_beta_mmpbsa_rerun.sh.
# ============================================================================

set -u

PY=/home/san/miniconda3/envs/qmmm/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh
conda activate qmmm

export UPDD_MMGBSA_PLATFORM=CUDA

LOG="${PROJ}/outputs/analysis/phase_beta_repbsa_v2_20260429.log"
mkdir -p "$(dirname "$LOG")"

echo "==========================================" | tee -a "$LOG"
echo "PHASE β RE-PBSA v2 START: $(date --iso-8601=seconds)" | tee -a "$LOG"
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
    local start
    start=$(date +%s)
    echo "[START ${label}] $(date --iso-8601=seconds) snap=${snap_dir} out=${outdir} target=${target_id} ncaa=${ncaa_elem}" | tee -a "$LOG"
    UPDD_MMGBSA_PLATFORM=CUDA "$PY" scripts/run_mmpbsa.py \
        --md_dir "${snap_dir}" \
        --outputdir "${outdir}" \
        --ncaa_elem "${ncaa_elem}" \
        --receptor_chain A --binder_chain B \
        --target_id "${target_id}" \
        --protocol 1traj \
        >> "$LOG" 2>&1
    local rc=$?
    local end
    end=$(date +%s)
    local elapsed=$((end - start))
    echo "[DONE ${label}] rc=${rc} elapsed=${elapsed}s" | tee -a "$LOG"
}

# Convenience wrapper: outputs/<sys>/snapshots_n25_postl387_patch_v2/
#                  →  outputs/<sys>/mmpbsa_results_postl387_v2/
run_sys () {
    local sys="$1"
    local target_id="$2"
    local ncaa_elem="$3"
    run_pbsa "${sys}" \
             "outputs/${sys}/snapshots_n25_postl387_patch_v2" \
             "outputs/${sys}/mmpbsa_results_postl387_v2" \
             "${target_id}" "${ncaa_elem}"
}

# ==========================================================================
# 10 systems — Phase β re-PBSA v2 (CONECT v55 unblocks #84 failures)
# ==========================================================================
echo "------------------------------------------" | tee -a "$LOG"
echo "Phase β re-PBSA v2 — 10 systems" | tee -a "$LOG"
echo "------------------------------------------" | tee -a "$LOG"

# 2QKI_Cp4 (8 seeds; Cp4 binder uses MTR for 4(1MeW))
run_sys 2QKI_Cp4_calib_s7            2QKI MTR
run_sys 2QKI_Cp4_calib_s19_reseed55  2QKI MTR
run_sys 2QKI_Cp4_calib_s23           2QKI MTR
run_sys 2QKI_Cp4_calib_s42           2QKI MTR
run_sys 2QKI_Cp4_calib_s83           2QKI MTR
run_sys 2QKI_Cp4_calib_s101          2QKI MTR
run_sys 2QKI_Cp4_calib_s163          2QKI MTR
run_sys 2QKI_Cp4_calib_s251          2QKI MTR

# 7TL8_MTR6 (2 INTACT seeds; s7 already done in #84 with v55-incompatible snap dir;
#            v55 path proved CONECT INTACT for s7 Cp4 too — but s7 in #84 used
#            the prior dir, not v2; aggregation will use #84's s7 unchanged.)
run_sys 7TL8_MTR6_calib_s19 7TL8 MTR
run_sys 7TL8_MTR6_calib_s42 7TL8 MTR

echo "==========================================" | tee -a "$LOG"
echo "PHASE β RE-PBSA v2 END: $(date --iso-8601=seconds)" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"
