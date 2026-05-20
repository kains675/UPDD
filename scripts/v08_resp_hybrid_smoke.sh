#!/bin/bash
# ============================================================================
# v0.8 Khoury RESP β hybrid pilot SMOKE TEST — sidechain-only Khoury + backbone/NE1
# baseline. Tests whether preserving N + H + NE1 amber14SB-patched values
# (with the other 24 atoms overridden by Khoury OMW RESP-A2) eliminates the
# MD-propagation-level instability observed in the v0.8 full-RESP pilot
# (NaN cascade in 2/8 seeds at 8% and 58% MD completion).
#
# Acceptance criterion: 4/4 seeds complete 5 ns MD without NaN.
#
# Per the v0.8 full-RESP pilot diagnosis (Q4 Option β rationale).
#
# Cohort: 4 seeds — 3 historically stable (s17, s47, s73) + 1 explosion-prone
# (s29, the 58% NaN seed; if hybrid stabilizes s29, this is the strongest
# signal of NE1 / backbone-N junction as causal).
#
# Stage A only (MD). Stage B/C/D deferred until 4/4 acceptance + Stage B
# orchestrator path fix in scripts/phase_beta_repbsa_v2_reextract.py (TARGETS
# list addition for the new cohort).
#
# Cohort dir: outputs/2QKI_Cp4_hybrid_calib_s{seed}/ (distinct from full-RESP
# 2QKI_Cp4_resp_calib_* and baseline 2QKI_Cp4_calib_*).
#
# Compute envelope: 4 seeds × 50min MD ≈ 3.3h sequential single-GPU.
# ============================================================================

set -u

PY=/home/san/miniconda3/envs/qmmm/bin/python
PY_MD=/home/san/miniconda3/envs/md_simulation/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"

source /home/san/miniconda3/etc/profile.d/conda.sh

# CPU affinity 자동 분배 (v08_resp_pilot.sh 와 동일 logic)
PHYS_CORES=$(lscpu -p 2>/dev/null | grep -v "^#" | awk -F, '{print $2}' | sort -nu | wc -l)
LOGICAL_CORES=$(nproc)
if [ "${PHYS_CORES:-0}" -ge 4 ] && [ "${PHYS_CORES:-0}" -le 8 ]; then
    MD_AFFINITY="0,${PHYS_CORES}"
    MD_PREFIX="taskset -c ${MD_AFFINITY}"
else
    MD_PREFIX=""
fi

# β hybrid charges path
HYBRID_XML="${PROJ}/params/MTR_gaff2_hybrid.xml"
HYDROGENS_SRC="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_hydrogens.xml"
MANIFEST_SRC_TEMPLATE="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_params_manifest.json"
CP4_REF_DIR="${PROJ}/outputs/2QKI_Cp4_calib_s7"

LAUNCH_TS=$(date +%s)
TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/v08_resp_hybrid_smoke_${TS_TAG}"
mkdir -p "$LOGDIR"
DISPATCH_LOG="${LOGDIR}/dispatch.log"

log () {
    printf '[v08_hybrid_smoke %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$DISPATCH_LOG"
}

abort () {
    log "ABORT: $*"
    exit 2
}

log "================================================================"
log "v0.8 Khoury RESP β hybrid SMOKE TEST START"
log "host=$(hostname)  gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
log "logdir=$LOGDIR"
log "Hybrid XML: $HYBRID_XML"
log "Cohort: 4 seeds — s17 + s47 + s73 (historically stable) + s29 (full-RESP explosion-prone)"
log "Acceptance criterion: 4/4 complete 5 ns MD"
log "================================================================"

# Pre-flight: verify hybrid XML exists + Σq close to zero + preserved atoms intact
if [ ! -f "$HYBRID_XML" ]; then
    abort "Hybrid XML not found: $HYBRID_XML — run scripts/build_khoury_mtr_hybrid_xml.py first."
fi

$PY -c "
import xml.etree.ElementTree as ET
t = ET.parse('$HYBRID_XML')
res = next(r for r in t.findall('.//Residue') if r.get('name') == 'MTR')
charges = {a.get('name'): float(a.get('charge')) for a in res.findall('Atom')}
total = sum(charges.values())
print(f'Sigma_q verification: {total:+.2e}')
assert abs(total) < 1e-5, f'Sigma_q not neutral: {total}'
# Verify preserved baseline values
expected = {'N': -0.4157, 'H': 0.2719, 'NE1': -0.3418}
for name, exp in expected.items():
    got = charges.get(name)
    assert got is not None and abs(got - exp) < 1e-4, f'{name} drift: got {got}, expected {exp}'
print(f'Preserved baseline N={charges[\"N\"]:+.4f} H={charges[\"H\"]:+.4f} NE1={charges[\"NE1\"]:+.4f}')
print('PASS')
" || abort "Hybrid XML verification failed"

# Smoke cohort: 3 stable + 1 explosion-prone
SMOKE_SEEDS=(17 47 73 29)

log "----------------------------------------------------------------"
log "Stage A: 4 hybrid-RESP MD seeds (sequential, single-GPU)"
log "----------------------------------------------------------------"

for seed_int in "${SMOKE_SEEDS[@]}"; do
    seed="s${seed_int}"
    seedir="${PROJ}/outputs/2QKI_Cp4_hybrid_calib_${seed}"

    if [ -f "${seedir}/mdresult/2QKI_Cp4_restrained.dcd" ]; then
        log "  $seed: SKIP (DCD exists)"
        continue
    fi

    log "  $seed: setup dir"
    mkdir -p "${seedir}/_md_input" "${seedir}/params" "${seedir}/mdresult"

    for f in 2QKI_Cp4.pdb 2QKI_Cp4_renum.pdb; do
        [ -f "${CP4_REF_DIR}/_md_input/${f}" ] && cp -n "${CP4_REF_DIR}/_md_input/${f}" "${seedir}/_md_input/"
    done

    cp -f "$HYBRID_XML" "${seedir}/params/MTR_gaff2.xml"
    cp -n "$HYDROGENS_SRC" "${seedir}/params/"
    cp -n "$MANIFEST_SRC_TEMPLATE" "${seedir}/params/"

    $PY -c "
import json
mf = '${seedir}/params/MTR_params_manifest.json'
with open(mf) as f: m = json.load(f)
m['xml_path'] = '${seedir}/params/MTR_gaff2.xml'
m['hydrogens_path'] = '${seedir}/params/MTR_hydrogens.xml'
m['charge_source'] = 'Khoury 2014 OMW RESP-A2 hybrid (sidechain-only; backbone N/H + NE1 preserved from amber14SB-patched baseline)'
with open(mf, 'w') as f: json.dump(m, f, indent=2)
"

    qN=$(grep -m1 '<Atom name="N" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    qNE1=$(grep -m1 '<Atom name="NE1" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    log "  $seed: hybrid q_N=$qN (preserved amber14SB-patched) / q_NE1=$qNE1 (preserved baseline)"

    log "  $seed: MD launch (5ns, cyclic_ss, hybrid charges, seed=${seed_int})"
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
            > "${LOGDIR}/md_hybrid_${seed}.log" 2>&1
    rc=$?
    conda deactivate

    if [ "$rc" -ne 0 ] || [ ! -f "${seedir}/mdresult/2QKI_Cp4_restrained.dcd" ]; then
        log "  $seed: MD FAIL rc=$rc (see ${LOGDIR}/md_hybrid_${seed}.log)"
    else
        log "  $seed: MD DONE"
    fi
done

# Final acceptance evaluation
SUCCESS=0
for seed_int in "${SMOKE_SEEDS[@]}"; do
    seed="s${seed_int}"
    if [ -f "${PROJ}/outputs/2QKI_Cp4_hybrid_calib_${seed}/mdresult/2QKI_Cp4_final.pdb" ]; then
        SUCCESS=$((SUCCESS+1))
        log "  ✓ ${seed}: Stage A success"
    else
        log "  ✗ ${seed}: Stage A FAIL"
    fi
done

log "================================================================"
log "Smoke test complete. ${SUCCESS}/4 seeds passed acceptance criterion."
if [ "$SUCCESS" -eq 4 ]; then
    log "🟢 ACCEPTANCE PASS — proceed to full n=8 hybrid pilot + Stage B/C/D path fix."
elif [ "$SUCCESS" -ge 3 ]; then
    log "🟡 MARGINAL — review failing seed pattern (3/4 stable suggests partial mitigation)."
else
    log "🔴 ACCEPTANCE FAIL — hybrid approach insufficient; escalate to Track D Phase γ (joint charges + bonds refit) or revert to baseline."
fi
log "FINAL state per seed: see Stage A summary above."
log "================================================================"
