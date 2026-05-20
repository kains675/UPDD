#!/bin/bash
# ============================================================================
# v0.8 RESP replicate pilot — Layer 3C (full Khoury) ×2 + Layer 3B (β hybrid) ×2
# per existing seed. Purpose: tighten Wilson 95 % CI on the per-layer crash
# rate and separate stochastic noise from systemic instability via per-seed
# 3-trajectory consistency check.
#
# Existing single-trajectory observations to be replicated:
#   Layer 3C (full Khoury, 8 seeds): s17/29/47/61/73/79/97/131
#     - 6 stable (s17/47/61/73/79/97) + 2 explode (s29 58 %, s131 8 %)
#   Layer 3B (β hybrid, 4 seeds): s17/47/73/29
#     - 3 stable (s17/47/29) + 1 explode (s73 8 %)
#
# Replication: 2 additional trajectories per seed under different OpenMM
# seed integers (original × 1000 + replicate_index), so each original seed
# becomes a 3-trajectory triplet.
#
# Cohort dir naming: outputs/2QKI_Cp4_{resp|hybrid}_calib_s{N}_r{R}/
#   (suffix _r2, _r3 distinguishes from original which has no suffix.)
#
# Stage A (MD stability) only — Stage B/C/D not in scope for replication test.
#
# Compute envelope: 24 trajectories × ~40 min sequential single GPU ≈ 16 h.
# ============================================================================

set -u

PY=/home/san/miniconda3/envs/qmmm/bin/python
PY_MD=/home/san/miniconda3/envs/md_simulation/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh

PHYS_CORES=$(lscpu -p 2>/dev/null | grep -v "^#" | awk -F, '{print $2}' | sort -nu | wc -l)
LOGICAL_CORES=$(nproc)
if [ "${PHYS_CORES:-0}" -ge 4 ] && [ "${PHYS_CORES:-0}" -le 8 ]; then
    MD_AFFINITY="0,${PHYS_CORES}"
    MD_PREFIX="taskset -c ${MD_AFFINITY}"
else
    MD_PREFIX=""
fi

# Charge model XMLs (built earlier)
RESP_XML="${PROJ}/params/MTR_gaff2_resp.xml"
HYBRID_XML="${PROJ}/params/MTR_gaff2_hybrid.xml"
HYDROGENS_SRC="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_hydrogens.xml"
MANIFEST_SRC_TEMPLATE="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_params_manifest.json"
CP4_REF_DIR="${PROJ}/outputs/2QKI_Cp4_calib_s7"

LAUNCH_TS=$(date +%s)
TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/v08_resp_replicate_${TS_TAG}"
mkdir -p "$LOGDIR"
DISPATCH_LOG="${LOGDIR}/dispatch.log"

log () {
    printf '[v08_replicate %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$DISPATCH_LOG"
}

abort () {
    log "ABORT: $*"
    exit 2
}

log "================================================================"
log "v0.8 RESP replicate pilot START"
log "host=$(hostname)  gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
log "logdir=$LOGDIR"
log "Compute: 24 trajectories sequential, ~16 h envelope"
log "Purpose: Wilson 95% CI tightening + per-seed stochastic vs systemic separation"
log "================================================================"

# Pre-flight: verify both XMLs exist
[ -f "$RESP_XML" ] || abort "Full-RESP XML not found: $RESP_XML"
[ -f "$HYBRID_XML" ] || abort "Hybrid XML not found: $HYBRID_XML"

LAYER_3C_SEEDS=(17 29 47 61 73 79 97 131)
LAYER_3B_SEEDS=(17 47 73 29)
REPLICATE_INDICES=(2 3)  # _r2 and _r3 (original is unsuffixed, effectively _r1)

run_md () {
    local layer="$1"           # "resp" or "hybrid"
    local xml_path="$2"
    local seed_int="$3"
    local rep_idx="$4"
    local charge_source_note="$5"

    local effective_seed=$(( seed_int * 1000 + rep_idx ))
    local seed_tag="s${seed_int}_r${rep_idx}"
    local seedir="${PROJ}/outputs/2QKI_Cp4_${layer}_calib_${seed_tag}"

    if [ -f "${seedir}/mdresult/2QKI_Cp4_restrained.dcd" ]; then
        log "  ${layer} ${seed_tag}: SKIP (DCD exists)"
        return 0
    fi

    log "  ${layer} ${seed_tag}: setup dir (effective_seed=${effective_seed})"
    mkdir -p "${seedir}/_md_input" "${seedir}/params" "${seedir}/mdresult"

    for f in 2QKI_Cp4.pdb 2QKI_Cp4_renum.pdb; do
        [ -f "${CP4_REF_DIR}/_md_input/${f}" ] && cp -n "${CP4_REF_DIR}/_md_input/${f}" "${seedir}/_md_input/"
    done

    cp -f "$xml_path" "${seedir}/params/MTR_gaff2.xml"
    cp -n "$HYDROGENS_SRC" "${seedir}/params/"
    cp -n "$MANIFEST_SRC_TEMPLATE" "${seedir}/params/"

    $PY -c "
import json
mf = '${seedir}/params/MTR_params_manifest.json'
with open(mf) as f: m = json.load(f)
m['xml_path'] = '${seedir}/params/MTR_gaff2.xml'
m['hydrogens_path'] = '${seedir}/params/MTR_hydrogens.xml'
m['charge_source'] = '${charge_source_note}'
m['replicate_index'] = ${rep_idx}
m['effective_seed'] = ${effective_seed}
m['origin_seed'] = ${seed_int}
with open(mf, 'w') as f: json.dump(m, f, indent=2)
"

    log "  ${layer} ${seed_tag}: MD launch (5ns, cyclic_ss, seed=${effective_seed})"
    conda activate md_simulation
    UPDD_NCAA_AMBER14_PATCH=0 UPDD_MMGBSA_PLATFORM=CUDA \
        $MD_PREFIX "$PY_MD" utils/run_restrained_md.py \
            --inputdir "${seedir}/_md_input" \
            --outputdir "${seedir}/mdresult" \
            --params_manifest "${seedir}/params/MTR_params_manifest.json" \
            --steps 2500000 \
            --topology cyclic_ss \
            --binder_chain B \
            --graph_policy strict \
            --target_id 2QKI \
            --dt_fs 2.0 \
            --platform CUDA \
            --seed "$effective_seed" \
            --ncaa_label MTR \
            --ncaa_code MTR \
            > "${LOGDIR}/md_${layer}_${seed_tag}.log" 2>&1
    local rc=$?
    conda deactivate

    if [ "$rc" -ne 0 ] || [ ! -f "${seedir}/mdresult/2QKI_Cp4_restrained.dcd" ]; then
        log "  ${layer} ${seed_tag}: MD FAIL rc=$rc"
    else
        log "  ${layer} ${seed_tag}: MD DONE"
    fi
}

# Stage 1: Layer 3C replicates (16 trajectories)
log "----------------------------------------------------------------"
log "Stage 1: Layer 3C (full Khoury RESP-A2) replicates"
log "----------------------------------------------------------------"
for seed_int in "${LAYER_3C_SEEDS[@]}"; do
    for rep_idx in "${REPLICATE_INDICES[@]}"; do
        run_md "resp" "$RESP_XML" "$seed_int" "$rep_idx" \
            "Khoury 2014 OMW RESP-A2 (DOI 10.1021/sb400168u) replicate r${rep_idx}"
    done
done

# Stage 2: Layer 3B replicates (8 trajectories)
log "----------------------------------------------------------------"
log "Stage 2: Layer 3B (β hybrid sidechain-only Khoury) replicates"
log "----------------------------------------------------------------"
for seed_int in "${LAYER_3B_SEEDS[@]}"; do
    for rep_idx in "${REPLICATE_INDICES[@]}"; do
        run_md "hybrid" "$HYBRID_XML" "$seed_int" "$rep_idx" \
            "Khoury 2014 OMW RESP-A2 hybrid (sidechain-only, N/H/NE1 preserved) replicate r${rep_idx}"
    done
done

# Final acceptance evaluation
log "================================================================"
log "Replicate pilot COMPLETE"
log "----------------------------------------------------------------"

LAYER_3C_PASS=0; LAYER_3C_FAIL=0
for seed_int in "${LAYER_3C_SEEDS[@]}"; do
    for rep_idx in "${REPLICATE_INDICES[@]}"; do
        F="${PROJ}/outputs/2QKI_Cp4_resp_calib_s${seed_int}_r${rep_idx}/mdresult/2QKI_Cp4_final.pdb"
        if [ -f "$F" ]; then
            LAYER_3C_PASS=$((LAYER_3C_PASS+1))
        else
            LAYER_3C_FAIL=$((LAYER_3C_FAIL+1))
        fi
    done
done

LAYER_3B_PASS=0; LAYER_3B_FAIL=0
for seed_int in "${LAYER_3B_SEEDS[@]}"; do
    for rep_idx in "${REPLICATE_INDICES[@]}"; do
        F="${PROJ}/outputs/2QKI_Cp4_hybrid_calib_s${seed_int}_r${rep_idx}/mdresult/2QKI_Cp4_final.pdb"
        if [ -f "$F" ]; then
            LAYER_3B_PASS=$((LAYER_3B_PASS+1))
        else
            LAYER_3B_FAIL=$((LAYER_3B_FAIL+1))
        fi
    done
done

log "Layer 3C replicates: ${LAYER_3C_PASS} PASS / ${LAYER_3C_FAIL} FAIL (of 16)"
log "Layer 3B replicates: ${LAYER_3B_PASS} PASS / ${LAYER_3B_FAIL} FAIL (of 8)"
log "Combined with original observations:"
log "  Layer 3C total: $((LAYER_3C_PASS + 6)) / $((LAYER_3C_PASS + LAYER_3C_FAIL + 8)) stable"
log "  Layer 3B total: $((LAYER_3B_PASS + 3)) / $((LAYER_3B_PASS + LAYER_3B_FAIL + 4)) stable"
log "================================================================"
