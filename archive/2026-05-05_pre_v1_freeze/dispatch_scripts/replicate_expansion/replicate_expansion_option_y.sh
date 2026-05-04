#!/bin/bash
# UPDD Option Y — Replicate expansion (seeds 23, 101)
# Rev. 5 Gate 2 resolution pivot: snap count axis exhausted at 5ns MD
# (25/29 runs: snap-spacing/tau < 1.0). New axis = replicate (3 -> 5 reps).
#
# Targets (user priority):
#   Tier 1 (WT, clean baselines): 2QKI, 3IOL, 2QKH, 1YCR, 1EBP   -> 5 WT
#   Tier 2 (variants, GPU permitting): 2QKI_Cp4, 3IOL_NML20,
#                                       1YCR_NML22, 1EBP_MTR13   -> 4 variants
#   Not in scope: 7TL8_MTR6 (bimodality), 2QKH_MTR25 (prior 0/3 valid), KRAS
#
# 2 new seeds per target -> 10 (tier 1) + 8 (tier 2 if included) = 10-18 runs
# Est ~35-45 min/run -> 7-14h total GPU sequential
# MM-GBSA: n=25 snapshots, v0.6.7 3-traj protocol (no UPDD_MMGBSA_PROTOCOL env)
#
# Safety: skip if MM-GBSA JSON already exists (idempotent)
#         Pass 1b reseed on MD failure (seed+13)
#         GPU budget halt: abort if >14h elapsed

set -u
PY=/home/san/miniconda3/envs/md_simulation/bin/python
export PATH=/home/san/miniconda3/envs/md_simulation/bin:$PATH
cd /home/san/UPDD_proj

SEEDS=(23 101)
N_SNAP=25
LAUNCH_TS=$(date +%s)
GPU_BUDGET_H=14

# Priority tier inclusion (toggle variants off if GPU tight)
INCLUDE_TIER2=${OPTION_Y_INCLUDE_TIER2:-1}

log_status() {
    echo "[optionY $(date +%H:%M)] $*"
}

check_budget() {
    local now=$(date +%s)
    local elapsed_h=$(( (now - LAUNCH_TS) / 3600 ))
    if [ "$elapsed_h" -ge "$GPU_BUDGET_H" ]; then
        log_status "GPU budget ${GPU_BUDGET_H}h exceeded (elapsed=${elapsed_h}h) - HALT"
        return 1
    fi
    return 0
}

run_md_full() {
    local outdir="$1"; local pdb_in="$2"; local topo="$3"; local ncaa_label="$4"; local ncaa_code="$5"
    local target_id="$6"; local seed="$7"; local manifest="$8"
    local steps=2500000  # 5ns x 2fs
    local tag="${outdir##*/}"

    if [ -f "${outdir}/mmgbsa_results/mmgbsa_summary.json" ]; then
        log_status "SKIP ${tag} (already done)"
        return 0
    fi
    check_budget || return 2

    log_status "MD ${tag} seed=${seed} 5ns..."
    local args=(--inputdir "$(dirname "$pdb_in")" --outputdir "${outdir}/mdresult"
        --steps "$steps" --topology "$topo" --binder_chain B --graph_policy strict
        --target_id "$target_id" --dt_fs 2.0 --platform CUDA --seed "$seed")
    if [ "$ncaa_code" != "" ]; then
        args+=(--params_manifest "$manifest" --ncaa_label "$ncaa_label" --ncaa_code "$ncaa_code")
    else
        args+=(--ncaa_label none --ncaa_code "")
    fi
    UPDD_MMGBSA_PLATFORM=CUDA $PY utils/run_restrained_md.py "${args[@]}" \
        > "/tmp/optionY_${tag}.log" 2>&1

    # Pass 1b-style seed retry
    if ! grep -qE "성공 \((ncAA|야생형) MD\) : 1 개|부분 성공 \(Recovery\) : 1 개|PARTIAL_SUCCESS" "/tmp/optionY_${tag}.log" 2>/dev/null; then
        local reseed=$((seed+13))
        log_status "WARN ${tag} seed=${seed} fail -> seed=${reseed} retry"
        rm -f "${outdir}/mdresult/"*EXPLODED* "${outdir}/mdresult/"*_partial* \
              "${outdir}/mdresult/"*_final.pdb "${outdir}/mdresult/"*_md.log \
              "${outdir}/mdresult/"*_restrained.dcd 2>/dev/null
        local args2=(--inputdir "$(dirname "$pdb_in")" --outputdir "${outdir}/mdresult"
            --steps "$steps" --topology "$topo" --binder_chain B --graph_policy strict
            --target_id "$target_id" --dt_fs 2.0 --platform CUDA --seed "$reseed")
        if [ "$ncaa_code" != "" ]; then
            args2+=(--params_manifest "$manifest" --ncaa_label "$ncaa_label" --ncaa_code "$ncaa_code")
        else
            args2+=(--ncaa_label none --ncaa_code "")
        fi
        UPDD_MMGBSA_PLATFORM=CUDA $PY utils/run_restrained_md.py "${args2[@]}" \
            > "/tmp/optionY_${tag}_reseed.log" 2>&1
    fi
    if ! grep -qE "성공 \((ncAA|야생형) MD\) : 1 개|최종 구조 저장" "/tmp/optionY_${tag}.log" "/tmp/optionY_${tag}_reseed.log" 2>/dev/null; then
        log_status "FAIL ${tag} 2 seeds both exploded - skipping snap/MM-GBSA"
        return 1
    fi

    # Snap (n=25) + MM-GBSA (v0.6.7 3-traj default)
    log_status "snap+MMGBSA ${tag} (n=${N_SNAP})..."
    $PY utils/extract_snapshots.py \
        --md_dir "${outdir}/mdresult" --outputdir "${outdir}/snapshots" \
        --n_snapshots "$N_SNAP" --binder_chain B --target_id "$target_id" \
        >> "/tmp/optionY_${tag}.log" 2>&1
    UPDD_MMGBSA_PLATFORM=CUDA $PY utils/run_mmgbsa.py \
        --md_dir "${outdir}/snapshots" --outputdir "${outdir}/mmgbsa_results" \
        --ncaa_elem none --receptor_chain A --binder_chain B --target_id "$target_id" \
        >> "/tmp/optionY_${tag}.log" 2>&1
    if [ -f "${outdir}/mmgbsa_results/mmgbsa_summary.json" ]; then
        $PY -c "
import json
d = json.load(open('${outdir}/mmgbsa_results/mmgbsa_summary.json'))
ie = list(d.get('interaction_entropy_per_design', {}).values())
s = ie[0] if ie else {'std_delta_g_kcal': 0}
n = d.get('n_calc', 0)
print(f'[DONE ${tag}] <dG>={d[\"mean_dg\"]:+.2f} sigma={s[\"std_delta_g_kcal\"]:.2f} n={n}')"
    else
        log_status "FAIL ${tag} no MM-GBSA summary"
        return 1
    fi
    return 0
}

mutate_and_prep() {
    local outdir="$1"; local src_pdb="$2"; local ncaa_code="$3"; local resid="$4"; local dst_stem="$5"
    mkdir -p "${outdir}/_md_input" "${outdir}/mdresult" "${outdir}/snapshots" "${outdir}/mmgbsa_results" "${outdir}/params"
    $PY -c "
import sys; sys.path.insert(0, '/home/san/UPDD_proj')
from utils.ncaa_mutate import mutate_pdb
from utils import ncaa_registry as r
d = next(x for x in r.NCAA_REGISTRY_DATA if x.code == '${ncaa_code}')
mutate_pdb('${src_pdb}', [('B', ${resid})], d, '${outdir}/_md_input/${dst_stem}.pdb')
"
    local xml_resname=$($PY -c "
import sys; sys.path.insert(0, '/home/san/UPDD_proj')
from utils import ncaa_registry as r
print(next(x for x in r.NCAA_REGISTRY_DATA if x.code == '${ncaa_code}').xml_resname)
")
    cp /tmp/calib_params/${xml_resname}_* "${outdir}/_md_input/"
    cp /tmp/calib_params/${ncaa_code}_params_manifest.json "${outdir}/_md_input/"
    cp /tmp/calib_params/${xml_resname}_* "${outdir}/params/"
    cp /tmp/calib_params/${ncaa_code}_params_manifest.json "${outdir}/params/"
}

# Verify ncAA params cache
if [ ! -f "/tmp/calib_params/MTR_params_manifest.json" ] \
   || [ ! -f "/tmp/calib_params/NML_params_manifest.json" ]; then
    log_status "ERR ncAA params cache missing at /tmp/calib_params/ - regenerate first"
    exit 1
fi
log_status "params cache OK: $(ls /tmp/calib_params/*.xml 2>/dev/null | wc -l) XML files"

# ============================================================
# Priority order configs:
#   target_id | topology | wt_pdb | ncaa_code | xml_resname | resid | variant_label
# ============================================================
# Tier 1 (WT, clean baselines) - strict priority order per user
declare -a TIER1_CONFIGS=(
    "2QKI      cyclic_ss   outputs/2QKI_replicate_seed42/_md_input/2QKI_compstatin_4W9A.pdb  MTR  MTR  4   Cp4"
    "3IOL      linear      target/3IOL_clean.pdb                                               NML  MLE  20  NML20"
    "2QKH      linear      target/2QKH_clean.pdb                                               MTR  MTR  25  MTR25"
    "1YCR      linear      target/1YCR_clean.pdb                                               NML  MLE  22  NML22"
    "1EBP      cyclic_ss   target/1EBP_clean.pdb                                               MTR  MTR  13  MTR13"
)
# Tier 2 (variants) - 2QKH_MTR25 excluded (prior 0/3 valid); 7TL8_MTR6 excluded (bimodality)
declare -a TIER2_CONFIGS=(
    "2QKI      cyclic_ss   outputs/2QKI_replicate_seed42/_md_input/2QKI_compstatin_4W9A.pdb  MTR  MTR  4   Cp4"
    "3IOL      linear      target/3IOL_clean.pdb                                               NML  MLE  20  NML20"
    "1YCR      linear      target/1YCR_clean.pdb                                               NML  MLE  22  NML22"
    "1EBP      cyclic_ss   target/1EBP_clean.pdb                                               MTR  MTR  13  MTR13"
)

# --- TIER 1: WT runs ---
log_status "=== TIER 1 (WT) launch ==="
for cfg in "${TIER1_CONFIGS[@]}"; do
    read -r target_id topo wt_pdb ncaa_code xml_resname resid variant <<< "$cfg"
    for seed in "${SEEDS[@]}"; do
        check_budget || { log_status "BUDGET HALT before ${target_id}_WT_s${seed}"; break 2; }
        D="outputs/${target_id}_WT_calib_s${seed}"
        if [ -f "${D}/mmgbsa_results/mmgbsa_summary.json" ]; then
            log_status "SKIP ${target_id}_WT_s${seed} (summary exists)"
            continue
        fi
        mkdir -p "${D}/_md_input" "${D}/mdresult" "${D}/snapshots" "${D}/mmgbsa_results"
        cp "$wt_pdb" "${D}/_md_input/${target_id}_WT.pdb"
        run_md_full "$D" "${D}/_md_input/${target_id}_WT.pdb" "$topo" "none" "" "$target_id" "$seed" ""
    done
done

# --- TIER 2: variants (GPU budget permitting) ---
if [ "$INCLUDE_TIER2" = "1" ]; then
    log_status "=== TIER 2 (variants) launch ==="
    for cfg in "${TIER2_CONFIGS[@]}"; do
        read -r target_id topo wt_pdb ncaa_code xml_resname resid variant <<< "$cfg"
        for seed in "${SEEDS[@]}"; do
            check_budget || { log_status "BUDGET HALT before ${target_id}_${variant}_s${seed}"; break 2; }
            D="outputs/${target_id}_${variant}_calib_s${seed}"
            if [ -f "${D}/mmgbsa_results/mmgbsa_summary.json" ]; then
                log_status "SKIP ${target_id}_${variant}_s${seed} (summary exists)"
                continue
            fi
            mutate_and_prep "$D" "$wt_pdb" "$ncaa_code" "$resid" "${target_id}_${variant}"
            run_md_full "$D" "${D}/_md_input/${target_id}_${variant}.pdb" "$topo" \
                        "$ncaa_code" "$xml_resname" "$target_id" "$seed" \
                        "${D}/_md_input/${ncaa_code}_params_manifest.json"
        done
    done
else
    log_status "=== TIER 2 SKIPPED (OPTION_Y_INCLUDE_TIER2=0) ==="
fi

log_status "========================================================"
log_status "OPTION Y ALL RUNS COMPLETE (or budget halted)"
log_status "========================================================"

# Aggregate status JSON
mkdir -p outputs/analysis
$PY << 'PYEOF'
import json, glob, os
results = {}
for seed in (23, 101):
    for jf in sorted(glob.glob(f'/home/san/UPDD_proj/outputs/*_calib_s{seed}/mmgbsa_results/mmgbsa_summary.json')):
        tag = os.path.basename(os.path.dirname(os.path.dirname(jf)))
        d = json.load(open(jf))
        ie = list(d.get('interaction_entropy_per_design', {}).values())
        s = ie[0] if ie else {'std_delta_g_kcal': 0}
        results[tag] = {
            'mean_dg': d['mean_dg'],
            'sigma_wi': s['std_delta_g_kcal'],
            'n_calc': d.get('n_calc', 0),
            'seed': seed,
        }
out = {
    'generated_ts': __import__('time').time(),
    'n_completed': len(results),
    'runs': results,
}
with open('/home/san/UPDD_proj/outputs/analysis/option_y_replicate_status_20260422.json', 'w') as f:
    json.dump(out, f, indent=2)
print(f'[option_y] status: {len(results)} new runs completed')
for tag in sorted(results):
    r = results[tag]
    print(f'  {tag:<35s} <dG>={r["mean_dg"]:+.2f} sigma={r["sigma_wi"]:.2f} n={r["n_calc"]}')
PYEOF
