#!/bin/bash
# ============================================================================
# v0.9 γ (Layer 3D) production pilot — Cp4 ↔ WT paired ΔΔG re-test
# ----------------------------------------------------------------------------
# Goal: deploy γ XML (params/MTR_gaff2_layer3d_gamma.xml) on V100 + 5070 Ti
# dual-GPU to re-measure the Phase 1.5 paired Cp4 ↔ WT ΔΔG (= +12.26 kcal/mol
# under v0.7.1 patch-on XML, vs Magotti SSOT [-3.0, -1.4] kcal/mol) using γ
# XML on the same Cp4 cohort. WT cohort is bit-identical under γ (γ alters
# only the MTR residue), so WT MD is re-used from the existing cohort.
#
# Cohort (matches manuscript §6.1 + Appendix S-X table):
#   Cp4 n=8: s7, s23, s42, s83, s101, s163, s251, s19_reseed55 (Phase 1.5)
#   WT  n=10: s7, s19, s23, s42, s83, s101, s127, s163, s199, s251 (existing)
#
# Stages:
#   A: γ MD launch on Cp4 n=8 (V100 device pin, ~30-60 min per seed,
#      total ~4-8 h sequential)
#   B: snapshot reextract postl387_v2 for the new γ-MD trajectories (CPU)
#   C: MMPBSA postl387_v2 for Cp4 n=8 γ + verify WT n=10 existing summaries
#      (parallel lanes, CUDA)
#   D: branched_ddg paired aggregation (γ Cp4 vs existing WT)
#
# Acceptance criteria (manuscript §6.1):
#   1. Sign recovery: paired ΔΔG sign matches Magotti SSOT [-3.0, -1.4]
#   2. Magnitude convergence: |ΔΔG_gamma − Magotti| ≤ 3 kcal/mol
#   3. σ_btwn stability: |σ_btwn_gamma − 6.76| / 6.76 ≤ 1× (i.e., 0 - 13.52)
# Informative failure modes (manuscript §6.1):
#   a. Sign retained X.A → methodology limit, v1.0+ R&D
#   b. Sign recovered, magnitude > 3 kcal/mol off → Δq budget refinement
#   c. σ_btwn departure > 1.5× → sampling adequacy under Definition 3
#
# Device policy:
#   - V100 (GPU 1): primary γ MD lane (Stage A)
#   - 5070 Ti (GPU 0): MMPBSA CUDA lane (Stage C) — can run concurrent w/ A
#     once the first MD seed completes
#
# Wall-clock envelope (estimated, V100 + 5070 Ti):
#   Stage A: 8 × ~30-60 min = 4-8 h
#   Stage B: 8 × ~10 min   = ~80 min
#   Stage C: 18 lanes / 4 parallel × ~5 min = ~25 min (existing WT skipped)
#   Stage D: ~1 min
#   Total: ~6-10 h
# ============================================================================

set -u

PY=/home/san/miniconda3/envs/qmmm/bin/python
PY_MD=/home/san/miniconda3/envs/md_simulation/bin/python
PROJ=/home/san/UPDD_proj
cd "$PROJ"
source /home/san/miniconda3/etc/profile.d/conda.sh

# ============================================================================
# Configuration
# ============================================================================
CP4_SEEDS=(7 23 42 83 101 163 251)
CP4_RESEED_LABEL="s19_reseed55"  # special-case: dir name s19_reseed55, integer seed 55
CP4_RESEED_INT=55

WT_SEEDS=(7 19 23 42 83 101 127 163 199 251)

GAMMA_XML="${PROJ}/params/MTR_gaff2_layer3d_gamma.xml"
CP4_REF_DIR="${PROJ}/outputs/2QKI_Cp4_calib_s7"
HYDROGENS_SRC="${CP4_REF_DIR}/params/MTR_hydrogens.xml"
MANIFEST_SRC_TEMPLATE="${CP4_REF_DIR}/params/MTR_params_manifest.json"

TS_TAG=$(date +%Y%m%d_%H%M%S)
LOGDIR="${PROJ}/outputs/analysis/v09_gamma_production_${TS_TAG}"
mkdir -p "$LOGDIR"
DISPATCH_LOG="${LOGDIR}/dispatch.log"

# Device pin: V100 on GPU 1 (default per plan/v100_dual_gpu_prep_20260519 §3.A);
# override at launch time with `UPDD_GAMMA_MD_GPU=N bash scripts/v09_gamma_production.sh`
GAMMA_MD_GPU="${UPDD_GAMMA_MD_GPU:-1}"
PBSA_GPU="${UPDD_GAMMA_PBSA_GPU:-0}"

# CPU affinity for MD (8-core preferred)
PHYS_CORES=$(lscpu -p 2>/dev/null | grep -v "^#" | awk -F, '{print $2}' | sort -nu | wc -l)
if [ "${PHYS_CORES:-0}" -ge 4 ] && [ "${PHYS_CORES:-0}" -le 8 ]; then
    MD_PREFIX="taskset -c 0,${PHYS_CORES}"
else
    MD_PREFIX=""
fi

PBSA_LANES="${UPDD_GAMMA_PBSA_LANES:-4}"

log () {
    printf '[v09_gamma_prod %s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" \
        | tee -a "$DISPATCH_LOG"
}

abort () {
    log "ABORT: $*"
    exit 2
}

log "================================================================"
log "v0.9 γ production pilot START"
log "host=$(hostname)  GPUs detected:"
nvidia-smi --query-gpu=index,name,memory.total,memory.used --format=csv,noheader 2>&1 \
    | tee -a "$DISPATCH_LOG"
log "logdir=$LOGDIR"
log "γ XML: $GAMMA_XML"
log "MD device pin:    GPU ${GAMMA_MD_GPU} (set UPDD_GAMMA_MD_GPU to override)"
log "PBSA device pin:  GPU ${PBSA_GPU}    (set UPDD_GAMMA_PBSA_GPU to override)"
log "PBSA lanes:       ${PBSA_LANES}     (set UPDD_GAMMA_PBSA_LANES to override)"
log "Cohort: Cp4 n=$((${#CP4_SEEDS[@]} + 1)) (${CP4_SEEDS[*]} + ${CP4_RESEED_LABEL})"
log "        WT  n=${#WT_SEEDS[@]} (${WT_SEEDS[*]})"
log "================================================================"

[ -f "$GAMMA_XML" ] || abort "γ XML not found: $GAMMA_XML"

# Σq sanity check
$PY -c "
import xml.etree.ElementTree as ET
t = ET.parse('$GAMMA_XML')
res = next(r for r in t.findall('.//Residue') if r.get('name') == 'MTR')
total = sum(float(a.get('charge')) for a in res.findall('Atom'))
print(f'gamma XML Sigma_q: {total:+.2e}')
assert abs(total) < 1e-5, f'Sigma_q not neutral: {total}'
print('PASS')
" || abort "γ XML Σq verification failed"

# ============================================================================
# Stage A — γ MD on Cp4 n=8 (V100 lane)
# ============================================================================
log "----------------------------------------------------------------"
log "Stage A: γ MD on Cp4 n=8 (sequential, GPU ${GAMMA_MD_GPU})"
log "----------------------------------------------------------------"

# Sequence: standard 7 seeds + s19_reseed55 special case
A_SEEDS_INT=("${CP4_SEEDS[@]}" "$CP4_RESEED_INT")
A_SEEDS_LABEL=()
for s in "${CP4_SEEDS[@]}"; do A_SEEDS_LABEL+=("s${s}"); done
A_SEEDS_LABEL+=("$CP4_RESEED_LABEL")

for i in "${!A_SEEDS_INT[@]}"; do
    seed_int="${A_SEEDS_INT[$i]}"
    seed_label="${A_SEEDS_LABEL[$i]}"
    seedir="${PROJ}/outputs/2QKI_Cp4_gamma_calib_${seed_label}"

    if [ -f "${seedir}/mdresult/2QKI_Cp4_final.pdb" ]; then
        log "  $seed_label: SKIP Stage A (final.pdb exists)"
        continue
    fi

    log "  $seed_label: setup dir"
    mkdir -p "${seedir}/_md_input" "${seedir}/params" "${seedir}/mdresult"

    for f in 2QKI_Cp4.pdb 2QKI_Cp4_renum.pdb; do
        [ -f "${CP4_REF_DIR}/_md_input/${f}" ] && \
            cp -n "${CP4_REF_DIR}/_md_input/${f}" "${seedir}/_md_input/"
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
m['charge_source'] = 'v0.9 gamma production: amber14SB Trp bonded + Khoury OMW full RESP-A2 + NA-*-*-CT improper'
with open(mf, 'w') as f: json.dump(m, f, indent=2)
"

    log "  $seed_label: MD launch (5ns, cyclic_ss, gamma XML, seed=${seed_int}, gpu=${GAMMA_MD_GPU})"
    conda activate md_simulation
    CUDA_VISIBLE_DEVICES="${GAMMA_MD_GPU}" \
        UPDD_NCAA_AMBER14_PATCH=0 \
        UPDD_MMGBSA_PLATFORM=CUDA \
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
            > "${LOGDIR}/md_gamma_${seed_label}.log" 2>&1
    rc=$?
    conda deactivate

    if [ "$rc" -ne 0 ] || [ ! -f "${seedir}/mdresult/2QKI_Cp4_final.pdb" ]; then
        log "  $seed_label: MD FAIL rc=$rc (see ${LOGDIR}/md_gamma_${seed_label}.log)"
    else
        log "  $seed_label: MD DONE"
    fi
done

# Stage A summary
A_PASS=0
for seed_label in "${A_SEEDS_LABEL[@]}"; do
    if [ -f "${PROJ}/outputs/2QKI_Cp4_gamma_calib_${seed_label}/mdresult/2QKI_Cp4_final.pdb" ]; then
        A_PASS=$((A_PASS+1))
    fi
done
log "Stage A: ${A_PASS}/${#A_SEEDS_LABEL[@]} Cp4 γ trajectories complete"

[ "$A_PASS" -ge 6 ] || abort "Stage A: only ${A_PASS}/8 trajectories completed — paired analysis power insufficient (need ≥6 for z_SE > 2 at σ_btwn~6.76)"

# ============================================================================
# Stage B — snapshot reextract postl387_v2 for Cp4 γ trajectories
# ============================================================================
log "----------------------------------------------------------------"
log "Stage B: snapshot reextract postl387_v2 (Cp4 γ n=${A_PASS}, CPU)"
log "----------------------------------------------------------------"
conda activate qmmm

for seed_label in "${A_SEEDS_LABEL[@]}"; do
    system_seed="2QKI_Cp4_gamma_calib_${seed_label}"
    target="${PROJ}/outputs/${system_seed}/snapshots_n25_postl387_patch_v2"
    if [ -d "$target" ] && [ "$(ls "$target"/*.pdb 2>/dev/null | wc -l)" -ge 25 ]; then
        log "  ${system_seed}: SKIP (25 snapshots exist)"
        continue
    fi
    if [ ! -f "${PROJ}/outputs/${system_seed}/mdresult/2QKI_Cp4_final.pdb" ]; then
        log "  ${system_seed}: SKIP (no final.pdb — Stage A failed)"
        continue
    fi
    log "  ${system_seed}: reextract"
    CUDA_VISIBLE_DEVICES="" $PY scripts/phase_beta_repbsa_v2_reextract.py \
        --filter "${system_seed}" \
        --out-json "${LOGDIR}/reextract_${system_seed}.json" \
        >> "${LOGDIR}/reextract.log" 2>&1
done

log "Stage B complete."

# ============================================================================
# Stage C — MMPBSA postl387_v2 (Cp4 γ n=8 + WT existing n=10)
# ============================================================================
log "----------------------------------------------------------------"
log "Stage C: MMPBSA postl387_v2 sweep (${PBSA_LANES} lanes, GPU ${PBSA_GPU})"
log "----------------------------------------------------------------"

export CUDA_VISIBLE_DEVICES="${PBSA_GPU}"
export UPDD_MMGBSA_PLATFORM=CUDA

PBSA_LOG="${LOGDIR}/mmpbsa.log"

run_one_pbsa () {
    local label="$1" snap_dir="$2" out_dir="$3" target_id="$4" ncaa_elem="$5"
    if [ -f "${out_dir}/mmpbsa_summary.json" ]; then
        printf '[%s SKIP %s] summary exists\n' "$(date +%H:%M:%S)" "$label" >> "$PBSA_LOG"
        return 0
    fi
    if ! ls "${snap_dir}"/*.pdb >/dev/null 2>&1; then
        printf '[%s MISS %s] no pdb\n' "$(date +%H:%M:%S)" "$label" >> "$PBSA_LOG"
        return 1
    fi
    mkdir -p "$out_dir"
    local seed_log="${LOGDIR}/mmpbsa_${label}.log"
    local start; start=$(date +%s)
    printf '[%s START %s]\n' "$(date +%H:%M:%S)" "$label" >> "$PBSA_LOG"
    "$PY" -u scripts/run_mmpbsa.py \
        --md_dir "$snap_dir" --outputdir "$out_dir" \
        --ncaa_elem "$ncaa_elem" --receptor_chain A --binder_chain B \
        --target_id "$target_id" --protocol 1traj \
        > "$seed_log" 2>&1
    local rc=$?
    local end; end=$(date +%s)
    printf '[%s DONE %s] rc=%d elapsed=%ds\n' "$(date +%H:%M:%S)" "$label" "$rc" "$((end-start))" >> "$PBSA_LOG"
    return $rc
}

# Build WORK list: Cp4 γ + WT (existing)
WORK=()
for seed_label in "${A_SEEDS_LABEL[@]}"; do
    WORK+=("Cp4_gamma_${seed_label}|outputs/2QKI_Cp4_gamma_calib_${seed_label}/snapshots_n25_postl387_patch_v2|outputs/2QKI_Cp4_gamma_calib_${seed_label}/mmpbsa_results_postl387_v2|2QKI|MTR")
done
for s in "${WT_SEEDS[@]}"; do
    WORK+=("WT_s${s}|outputs/2QKI_WT_calib_s${s}/snapshots_n25_postl387_patch_v2|outputs/2QKI_WT_calib_s${s}/mmpbsa_results_postl387_v2|2QKI|none")
done

log "Phase C: ${#WORK[@]} work units, ${PBSA_LANES} lanes"
running=0
for entry in "${WORK[@]}"; do
    IFS='|' read -r label snap out tgt ncaa <<< "$entry"
    run_one_pbsa "$label" "$snap" "$out" "$tgt" "$ncaa" &
    running=$((running + 1))
    if [ "$running" -ge "$PBSA_LANES" ]; then
        wait -n
        running=$((running - 1))
    fi
done
wait
log "Stage C complete."

# ============================================================================
# Stage D — branched_ddg paired aggregation (γ Cp4 ↔ existing WT)
# ============================================================================
log "----------------------------------------------------------------"
log "Stage D: branched_ddg paired (γ Cp4 ↔ existing WT)"
log "----------------------------------------------------------------"

mkdir -p "${LOGDIR}/branched_ddg"
"$PY" -m utils.branched_ddg \
    --wt "outputs/2QKI_WT_calib_*/mmpbsa_results_postl387_v2" \
    --variant "outputs/2QKI_Cp4_gamma_calib_*/mmpbsa_results_postl387_v2" \
    --output "${LOGDIR}/branched_ddg/" \
    >> "$DISPATCH_LOG" 2>&1 \
    && log "Stage D DONE" \
    || log "Stage D WARN non-zero exit"

# ============================================================================
# Final summary
# ============================================================================
log "================================================================"
log "v0.9 γ production pilot COMPLETE"
log "----------------------------------------------------------------"
log "Cp4 γ trajectories  (Stage A): ${A_PASS}/${#A_SEEDS_LABEL[@]}"

C_PASS=0
C_TOTAL="${#WORK[@]}"
for entry in "${WORK[@]}"; do
    IFS='|' read -r label snap out tgt ncaa <<< "$entry"
    if [ -f "${out}/mmpbsa_summary.json" ]; then
        C_PASS=$((C_PASS+1))
    fi
done
log "MMPBSA summaries   (Stage C): ${C_PASS}/${C_TOTAL}"

if [ -f "${LOGDIR}/branched_ddg/branched_ddg.json" ]; then
    log "branched_ddg.json: ${LOGDIR}/branched_ddg/branched_ddg.json"
    log "----------------------------------------------------------------"
    log "Acceptance criteria evaluation (manuscript §6.1):"
    "$PY" -c "
import json, math
with open('${LOGDIR}/branched_ddg/branched_ddg.json') as f:
    d = json.load(f)
v = d.get('variant', {}) if 'variant' in d else d
mu_v = v.get('mean')
sig_v = v.get('sigma_btwn')
ddg = d.get('paired_ddg', {}).get('mean') if isinstance(d.get('paired_ddg'), dict) else d.get('paired_ddg')
ci95 = d.get('paired_ddg', {}).get('ci95') if isinstance(d.get('paired_ddg'), dict) else None
print(f'  Cp4 γ ⟨⟨ΔG⟩⟩ = {mu_v}')
print(f'  Cp4 γ σ_btwn = {sig_v}')
print(f'  Paired ΔΔG   = {ddg}')
print(f'  CI95         = {ci95}')
# Magotti SSOT [-3.0, -1.4]; baseline σ_btwn 6.76
mag_lo, mag_hi = -3.0, -1.4
sig_baseline = 6.76
if ddg is None:
    print('  WARN: paired ΔΔG missing from branched_ddg.json — manual review required')
else:
    # Criterion 1: sign recovery
    sign_recovered = (ddg < 0)
    print(f'  Criterion 1 (sign recovery vs Magotti): {\"PASS\" if sign_recovered else \"FAIL\"}')
    # Criterion 2: magnitude convergence within 3 kcal/mol
    closest_anchor = max(mag_lo, min(mag_hi, ddg))
    mag_off = abs(ddg - closest_anchor)
    print(f'  Criterion 2 (|ΔΔG - Magotti| ≤ 3): {mag_off:.2f} kcal/mol — {\"PASS\" if mag_off <= 3.0 else \"FAIL\"}')
    # Criterion 3: σ_btwn stability within 1×
    if sig_v is not None:
        sig_ratio = abs(sig_v - sig_baseline) / sig_baseline
        print(f'  Criterion 3 (σ_btwn within 1× of 6.76): ratio={sig_ratio:.2f} — {\"PASS\" if sig_ratio <= 1.0 else \"FAIL\"}')
    else:
        print('  Criterion 3: σ_btwn missing — manual review required')
"
fi

log "================================================================"
log "Next steps:"
log "  1. Review ${LOGDIR}/branched_ddg/branched_ddg.json"
log "  2. Fill in outputs/paper1/v3_pending/v3_1_addendum_template.md"
log "  3. If all 3 criteria PASS → ChemRxiv v3.1 addendum submit"
log "  4. If any criterion FAIL → consult §6.1 informative failure mode mapping"
log "================================================================"
