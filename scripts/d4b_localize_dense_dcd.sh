#!/bin/bash
# ============================================================================
# D4b localize dense-DCD — root-cause discrimination for the residual D3-refit
# MTR MD crashes (P3 #106, follow-up to D4 stability SMOKE FAIL).
#
# Context: D4 stability SMOKE (5-seed 25ns, D3-refit MTR FF) ended 2/2 crash
# (s211 ~11.7ns, s277 ~8.7ns; abrupt local divergence). The residual instability
# must be discriminated between two hypotheses:
#   H1 (MTR-improper)  -- indole improper barrier (k=1.1) deferred in D3 refit;
#                         a force-field-term defect, timestep-INDEPENDENT.
#   H2 (generic dt=2fs) -- generic integrator over-energization at dt=2fs in the
#                         cyclic context; timestep-DEPENDENT.
#
# Method: re-run the two crashing seeds with a dense DCD (1ps/frame instead of
# the 50ps default) to LOCALIZE which internal coordinate diverges first
# (NE1-region vs backbone/SS), PLUS a dt=1fs control arm on s211 that converts
# the correlational localization into a near-causal discriminator:
#   dt=2fs crash + dt=1fs stable -> H2 (timestep instability)
#   dt=2fs crash + dt=1fs crash  -> H1 (FF-term defect, timestep-independent)
#
# Dense DCD is enabled via UPDD_MD_DCD_INTERVAL (env-var override of the single
# DCDReporter interval in run_restrained_md.py; default behaviour unchanged when
# unset). Frame capture during a crashing 100ps checkpoint chunk IS flushed to
# disk before the NaN exception, so the divergence TREND over the final dense
# frames is recoverable (the literal NaN frame is not captured -- read the trend
# over the last ~10-50 frames, not a single frame).
#
# Runs (3 total, sequential on a single GPU; concurrent execution forbidden =
# contention):
#   run  seed  dt_fs  --steps      DCD_INTERVAL  purpose            outdir suffix
#   1    211   2.0    12500000     500  (=1ps)   localize           _localize_s211_dt2fs
#   2    277   2.0    12500000     500  (=1ps)   localize backup    _localize_s277_dt2fs
#   3    211   1.0    25000000     2000 (=2ps)   dt control         _localize_s211_dt1fs
# (run 3 keeps 25ns physical time at dt=1fs = 25M steps; DCD 2000 step = 2ps.)
#
# FF path: the FROZEN refit XML params/MTR_gaff2_layer3d_d3refit.xml is installed
# read-only as the per-seed cohort params/MTR_gaff2.xml exactly like
# d4_stability_smoke.sh (anti-fragmentation). Only dt differs across runs; the FF
# is identical (the precondition that makes the dt arm a valid discriminator).
#
# Dirs: NEW per-run cohort dirs outputs/2QKI_Cp4_d3refit_localize_s{seed}_dt{2,1}fs.
# The previously crashed dirs outputs/2QKI_Cp4_d3refit_calib_s211/_s277 are
# PRESERVED and never clobbered (R-7).
#
# Honest tally: each run is judged by the presence of an _EXPLODED_ marker (NOT
# by exit code "DONE"). A crash here is EXPECTED and is the diagnostic signal --
# do NOT claim "fixed"/"stable" from this experiment (it is a root-cause
# diagnostic, not a fix). Acceptance/geometry analysis (NE1-CM bond, NE1
# improper out-of-plane, backbone N-CA/N-H, cyclic-SS Cys-Cys) runs separately
# on the dense DCDs after completion.
#
# Usage:
#   scripts/d4b_localize_dense_dcd.sh                       # 3 runs, sequential
#   CUDA_VISIBLE_DEVICES=0 scripts/d4b_localize_dense_dcd.sh   # pin GPU
#   DRY_RUN=1 scripts/d4b_localize_dense_dcd.sh             # print plan, no MD
# ============================================================================

set -u

PY=/home/san/miniconda3/envs/qmmm/bin/python
PY_MD=/home/san/miniconda3/envs/md_simulation/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ" || { echo "cannot cd $PROJ"; exit 2; }

source /home/san/miniconda3/etc/profile.d/conda.sh

# ----- knobs (env-overridable) -----
DRY_RUN="${DRY_RUN:-0}"

# ----- CPU affinity (d4_stability_smoke.sh logic) -----
PHYS_CORES=$(lscpu -p 2>/dev/null | grep -v "^#" | awk -F, '{print $2}' | sort -nu | wc -l)
LOGICAL_CORES=$(nproc)
if [ "${PHYS_CORES:-0}" -ge 4 ] && [ "${PHYS_CORES:-0}" -le 8 ]; then
    MD_AFFINITY="0,${PHYS_CORES}"
    MD_PREFIX="taskset -c ${MD_AFFINITY}"
else
    MD_PREFIX=""
fi

# ----- refit FF path (frozen, read-only) -----
REFIT_XML="${PROJ}/params/MTR_gaff2_layer3d_d3refit.xml"
HYDROGENS_SRC="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_hydrogens.xml"
MANIFEST_SRC_TEMPLATE="${PROJ}/outputs/2QKI_Cp4_calib_s7/params/MTR_params_manifest.json"
CP4_REF_DIR="${PROJ}/outputs/2QKI_Cp4_calib_s7"

TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/d4b_localize_dense_dcd_${TS_TAG}"
mkdir -p "$LOGDIR"
mkdir -p "${PROJ}/logs"
DISPATCH_LOG="${LOGDIR}/dispatch.log"
DRIVER_LOG="${PROJ}/logs/d4b_localize_dense_dcd.log"

log () {
    printf '[d4b_localize %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$DISPATCH_LOG" "$DRIVER_LOG"
}

abort () {
    log "ABORT: $*"
    exit 2
}

log "================================================================"
log "D4b localize dense-DCD START (D3 refit MTR FF, 2QKI_Cp4 cyclic_ss)"
log "host=$(hostname)  gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
log "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<all>}"
log "logdir=$LOGDIR  driver_log=$DRIVER_LOG"
log "refit XML (frozen): $REFIT_XML"
log "3 runs sequential: s211 dt2fs (1ps DCD) | s277 dt2fs (1ps DCD) | s211 dt1fs (2ps DCD)"
log "Goal: localize divergent coordinate (NE1 vs backbone/SS) + discriminate H1(improper)/H2(dt2fs)"
log "FORBIDDEN (R-11/R-18): 'fixed'/'stable'/lambda-rate claims -- this is a root-cause diagnostic"
log "================================================================"

# ----- Pre-flight: refit XML must exist + Sigma q == 0 + frozen anchors intact -----
# (identical verification to d4_stability_smoke.sh)
if [ ! -f "$REFIT_XML" ]; then
    abort "refit XML not found: $REFIT_XML -- run scripts/d4_build_refit_xml.py first."
fi

$PY -c "
import sys
import xml.etree.ElementTree as ET
t = ET.parse('$REFIT_XML')
res = next(r for r in t.findall('.//Residue') if r.get('name') == 'MTR')
charges = {a.get('name'): float(a.get('charge')) for a in res.findall('Atom')}
total = sum(charges.values())
print('Sigma_q (MTR) = %+.2e' % total)
assert abs(total) < 1e-5, 'Sigma_q not neutral: %r' % total
# frozen anchors (amber14SB Maier)
expected = {'N': -0.4157, 'H': 0.2719, 'CA': -0.0275, 'HA': 0.1123,
            'C': 0.5973, 'O': -0.5679, 'NE1': -0.3418}
for name, exp in expected.items():
    got = charges.get(name)
    assert got is not None and abs(got - exp) < 1e-4, \
        '%s drift: got %r expected %r' % (name, got, exp)
# NE1-CM bond r_eq + new rotor Proper present
bond = next((b for b in t.findall('.//HarmonicBondForce/Bond')
             if b.get('type1') == 'protein-CT' and b.get('type2') == 'protein-NA'), None)
assert bond is not None and abs(float(bond.get('length')) - 0.14482) < 1e-6, \
    'NE1-CM bond r_eq != 0.14482: %r' % (bond.get('length') if bond is not None else None)
rotor = next((p for p in t.findall('.//PeriodicTorsionForce/Proper')
              if (p.get('type1'), p.get('type2'), p.get('type3'), p.get('type4'))
              == ('protein-CW', 'protein-NA', 'protein-CT', 'protein-HC')), None)
assert rotor is not None, 'N-methyl rotor Proper missing'
assert abs(float(rotor.get('k1')) - 0.461282) < 1e-3, \
    'rotor k1 != 0.461282: %r' % rotor.get('k1')
assert int(rotor.get('periodicity1')) == 3 and abs(float(rotor.get('phase1'))) < 1e-9
print('Pre-flight PASS: Sigma q=%+.2e, 7 frozen anchors intact, NE1-CM r_eq=0.14482, rotor k1=%.6f'
      % (total, float(rotor.get('k1'))))
" || abort "refit XML pre-flight verification failed"

# ----- one-arm launcher helper: (seed_int, dt_fs, steps, dcd_interval, suffix) -----
# Mirrors d4_stability_smoke.sh::run_one_seed per-seed setup (fresh dir, frozen
# XML install, manifest patch, charge echo) but with a per-run dt / steps / dense
# DCD interval and a distinct localize cohort dir.
run_one_arm () {
    local seed_int="$1"
    local dt_fs="$2"
    local steps="$3"
    local dcd_interval="$4"
    local suffix="$5"
    local seedir="${PROJ}/outputs/2QKI_Cp4_d3refit_${suffix}"

    log "  [${suffix}] setup dir $seedir (seed=${seed_int}, dt=${dt_fs}fs, steps=${steps}, DCD=${dcd_interval} step)"
    mkdir -p "${seedir}/_md_input" "${seedir}/params" "${seedir}/mdresult"

    for f in 2QKI_Cp4.pdb 2QKI_Cp4_renum.pdb; do
        [ -f "${CP4_REF_DIR}/_md_input/${f}" ] && [ ! -f "${seedir}/_md_input/${f}" ] && \
            cp "${CP4_REF_DIR}/_md_input/${f}" "${seedir}/_md_input/"
    done

    # install frozen refit XML as cohort MTR_gaff2.xml (read-only source; the
    # FROZEN template + refit source are never modified -- cp only).
    cp -f "$REFIT_XML" "${seedir}/params/MTR_gaff2.xml"
    [ -f "$HYDROGENS_SRC" ] && [ ! -f "${seedir}/params/$(basename "$HYDROGENS_SRC")" ] && \
        cp "$HYDROGENS_SRC" "${seedir}/params/"
    [ -f "$MANIFEST_SRC_TEMPLATE" ] && [ ! -f "${seedir}/params/$(basename "$MANIFEST_SRC_TEMPLATE")" ] && \
        cp "$MANIFEST_SRC_TEMPLATE" "${seedir}/params/"

    $PY -c "
import json
mf = '${seedir}/params/MTR_params_manifest.json'
with open(mf) as f: m = json.load(f)
m['xml_path'] = '${seedir}/params/MTR_gaff2.xml'
m['hydrogens_path'] = '${seedir}/params/MTR_hydrogens.xml'
m['charge_source'] = 'D3 RESP-A2 (NE1/backbone frozen amber14SB) + NA-CT r_eq 1.4482A + N-methyl rotor V3; refit MTR_gaff2_layer3d_d3refit.xml'
with open(mf, 'w') as f: json.dump(m, f, indent=2)
" || { log "  [${suffix}] manifest patch FAIL"; return 1; }

    local qN qNE1 qCM
    qN=$(grep -m1 '<Atom name="N" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    qNE1=$(grep -m1 '<Atom name="NE1" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    qCM=$(grep -m1 '<Atom name="CM" ' "${seedir}/params/MTR_gaff2.xml" | grep -oP 'charge="\K[^"]+')
    log "  [${suffix}] refit q_N=$qN (frozen) q_NE1=$qNE1 (frozen) q_CM=$qCM (RESP-A2)"

    local ns
    ns=$(awk "BEGIN{print ${steps}*${dt_fs}/1e6}")
    if [ "$DRY_RUN" = "1" ]; then
        log "  [${suffix}] DRY_RUN -- would launch ${ns}ns MD (steps=${steps}, dt=${dt_fs}fs, DCD=${dcd_interval} step), not launching"
        return 0
    fi

    log "  [${suffix}] MD launch (${ns}ns physical, cyclic_ss, dt=${dt_fs}fs, dense DCD=${dcd_interval} step, seed=${seed_int})"
    conda activate md_simulation
    UPDD_NCAA_AMBER14_PATCH=0 UPDD_MMGBSA_PLATFORM=CUDA UPDD_MD_DCD_INTERVAL="$dcd_interval" \
        $MD_PREFIX "$PY_MD" utils/run_restrained_md.py \
            --inputdir "${seedir}/_md_input" \
            --outputdir "${seedir}/mdresult" \
            --params_manifest "${seedir}/params/MTR_params_manifest.json" \
            --steps "$steps" \
            --topology cyclic_ss \
            --binder_chain B \
            --graph_policy strict \
            --target_id 2QKI \
            --dt_fs "$dt_fs" \
            --platform CUDA \
            --seed "$seed_int" \
            --ncaa_label MTR \
            --ncaa_code MTR \
            > "${LOGDIR}/md_${suffix}.log" 2>&1
    local rc=$?
    conda deactivate

    log "  [${suffix}] run finished rc=$rc (see ${LOGDIR}/md_${suffix}.log)"
    return 0
}

log "----------------------------------------------------------------"
log "Stage A: 3 localize runs (sequential -- single GPU, concurrency forbidden)"
log "----------------------------------------------------------------"

# run 1: s211 dt=2fs, 25ns (12.5M step), DCD 500 step = 1ps   -- localize
run_one_arm 211 2.0 12500000 500  "localize_s211_dt2fs"
# run 2: s277 dt=2fs, 25ns (12.5M step), DCD 500 step = 1ps   -- localize backup
run_one_arm 277 2.0 12500000 500  "localize_s277_dt2fs"
# run 3: s211 dt=1fs, 25ns physical (25M step), DCD 2000 step = 2ps -- dt control
run_one_arm 211 1.0 25000000 2000 "localize_s211_dt1fs"

# ----- honest tally: crash judged by _EXPLODED_ marker, NOT exit code -----
# A crash here is the EXPECTED diagnostic signal; rc-based "DONE" would be a
# false-green (a partial-recovery run exits rc=0 yet exploded).
log "================================================================"
log "D4b localize Stage A done -- crash tally (by _EXPLODED_ marker):"
ARMS="localize_s211_dt2fs localize_s277_dt2fs localize_s211_dt1fs"
N_CRASH=0
N_COMPLETE=0
N_INCONCLUSIVE=0
for suffix in $ARMS; do
    seedir="${PROJ}/outputs/2QKI_Cp4_d3refit_${suffix}"
    expl=$(ls "${seedir}/mdresult/"_EXPLODED_* 2>/dev/null | head -1)
    final="${seedir}/mdresult/2QKI_Cp4_final.pdb"
    if [ -n "$expl" ]; then
        N_CRASH=$((N_CRASH + 1))
        log "  CRASH  ${suffix}: explosion marker $(basename "$expl")"
    elif [ -f "$final" ]; then
        N_COMPLETE=$((N_COMPLETE + 1))
        log "  COMPLETE ${suffix}: 25ns reached, no explosion marker (final.pdb present)"
    else
        N_INCONCLUSIVE=$((N_INCONCLUSIVE + 1))
        log "  INCONCLUSIVE ${suffix}: neither marker nor final.pdb (check ${LOGDIR}/md_${suffix}.log)"
    fi
done

log "----------------------------------------------------------------"
log "Tally: ${N_CRASH} crashed / ${N_COMPLETE} completed / ${N_INCONCLUSIVE} inconclusive (of 3)"
log "Discrimination (read dense-DCD divergence trend on crashed arms):"
log "  s211 dt2fs crash + s211 dt1fs STABLE  -> H2 (generic dt=2fs integrator instability)"
log "  s211 dt2fs crash + s211 dt1fs CRASH   -> H1 (FF-term defect, timestep-independent)"
log "  divergence at NE1-region (improper / NE1-CM bond) AND dt1fs crash -> H1 strongly favoured"
log "  divergence at backbone / cyclic-SS                                -> H1 excluded (H2)"
log "NEXT: read final ~10-50 dense frames of each crashed DCD (NE1-CM bond, NE1 improper"
log "  out-of-plane, backbone N-CA/N-H, cyclic-SS Cys-Cys). _final.pdb is recovery-only"
log "  (crash-100ps grid) -- do NOT use it for localization; read the DCD trend."
log "Honest framing (R-18): root LOCALIZE to H1/H2 only. NO 'fixed'/'stable'/lambda claims."
log "================================================================"
