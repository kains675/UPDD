#!/bin/bash
# ============================================================================
# v0.9 γ (Layer 3D fallback) smoke test
#
# Tests whether the γ XML (amber14SB Trp bonded + Khoury OMW charges) yields
# stable MD on 2QKI_Cp4 cohort. Same seed set as the v0.8 β hybrid smoke
# (s17 / s47 / s73 stable-baseline + s29 explosion-prone) to allow direct
# comparison.
#
# Acceptance criterion (per v0.9 D3 verdict §3):
#   4/4 trajectories complete 5 ns MD without NaN
#
# Wall-clock envelope: 4 × ~40 min = ~3 h sequential single-GPU.
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

GAMMA_XML="${PROJ}/params/MTR_gaff2_layer3d_gamma.xml"
HYDROGENS_SRC="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_hydrogens.xml"
MANIFEST_SRC_TEMPLATE="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_params_manifest.json"
CP4_REF_DIR="${PROJ}/outputs/2QKI_Cp4_calib_s7"

TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/v09_gamma_smoke_${TS_TAG}"
mkdir -p "$LOGDIR"
DISPATCH_LOG="${LOGDIR}/dispatch.log"

log () {
    printf '[v09_gamma %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$DISPATCH_LOG"
}

abort () {
    log "ABORT: $*"
    exit 2
}

log "================================================================"
log "v0.9 γ smoke test START"
log "host=$(hostname)  gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
log "logdir=$LOGDIR"
log "γ XML: $GAMMA_XML"
log "Cohort: 4 seeds — s17 + s47 + s73 + s29 (matches v0.8 β hybrid smoke for direct comparison)"
log "================================================================"

[ -f "$GAMMA_XML" ] || abort "γ XML not found: $GAMMA_XML"

$PY -c "
import xml.etree.ElementTree as ET
t = ET.parse('$GAMMA_XML')
res = next(r for r in t.findall('.//Residue') if r.get('name') == 'MTR')
total = sum(float(a.get('charge')) for a in res.findall('Atom'))
print(f'gamma XML Sigma_q: {total:+.2e}')
assert abs(total) < 1e-5, f'Sigma_q not neutral: {total}'
print('PASS')
" || abort "γ XML Σq verification failed"

SMOKE_SEEDS=(17 47 73 29)

log "----------------------------------------------------------------"
log "Stage A: 4 γ-MD seeds (sequential, single-GPU)"
log "----------------------------------------------------------------"

for seed_int in "${SMOKE_SEEDS[@]}"; do
    seed="s${seed_int}"
    seedir="${PROJ}/outputs/2QKI_Cp4_gamma_calib_${seed}"

    if [ -f "${seedir}/mdresult/2QKI_Cp4_restrained.dcd" ]; then
        log "  $seed: SKIP (DCD exists)"
        continue
    fi

    log "  $seed: setup dir"
    mkdir -p "${seedir}/_md_input" "${seedir}/params" "${seedir}/mdresult"

    for f in 2QKI_Cp4.pdb 2QKI_Cp4_renum.pdb; do
        [ -f "${CP4_REF_DIR}/_md_input/${f}" ] && cp -n "${CP4_REF_DIR}/_md_input/${f}" "${seedir}/_md_input/"
    done

    cp -f "$GAMMA_XML" "${seedir}/params/MTR_gaff2.xml"
    cp -n "$HYDROGENS_SRC" "${seedir}/params/"
    cp -n "$MANIFEST_SRC_TEMPLATE" "${seedir}/params/"

    $PY -c "
import json
mf = '${seedir}/params/MTR_params_manifest.json'
with open(mf) as f: m = json.load(f)
m['xml_path'] = '${seedir}/params/MTR_gaff2.xml'
m['hydrogens_path'] = '${seedir}/params/MTR_hydrogens.xml'
m['charge_source'] = 'v0.9 gamma: amber14SB Trp bonded + Khoury OMW charges + GAFF2 N-methyl'
with open(mf, 'w') as f: json.dump(m, f, indent=2)
"

    log "  $seed: MD launch (5ns, cyclic_ss, gamma XML, seed=${seed_int})"
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
            --seed "$seed_int" \
            --ncaa_label MTR \
            --ncaa_code MTR \
            > "${LOGDIR}/md_gamma_${seed}.log" 2>&1
    rc=$?
    conda deactivate

    if [ "$rc" -ne 0 ] || [ ! -f "${seedir}/mdresult/2QKI_Cp4_restrained.dcd" ]; then
        log "  $seed: MD FAIL rc=$rc"
    else
        log "  $seed: MD DONE"
    fi
done

SUCCESS=0
for seed_int in "${SMOKE_SEEDS[@]}"; do
    seed="s${seed_int}"
    if [ -f "${PROJ}/outputs/2QKI_Cp4_gamma_calib_${seed}/mdresult/2QKI_Cp4_final.pdb" ]; then
        SUCCESS=$((SUCCESS+1))
        log "  ✓ ${seed}: Stage A success"
    else
        log "  ✗ ${seed}: Stage A FAIL"
    fi
done

log "================================================================"
log "γ smoke complete. ${SUCCESS}/4 seeds passed."
if [ "$SUCCESS" -eq 4 ]; then
    log "🟢 ACCEPTANCE PASS — γ is a viable v0.9 Layer 3D fallback path."
elif [ "$SUCCESS" -ge 3 ]; then
    log "🟡 MARGINAL — review failing seed pattern (compare to β hybrid s73 regression at 8%)."
else
    log "🔴 ACCEPTANCE FAIL — γ insufficient; escalate to v1.0+ Layer 3D joint refit (δ Track D Phase γ)."
fi
log "================================================================"
