#!/bin/bash
# UPDD P1.1 n=10 expansion — 2QKI_Cp4 + 2QKI_WT only
# Rev. 5 dispatch 2026-04-22: X.A' -> X.A promotion via replicate stabilization
#
# Scope:
#   2QKI_WT  (existing n=5: seeds 7, 19, 23, 42, 101)
#   2QKI_Cp4 (existing n=5: seeds 7, 19, 23, 42, 101; plus s19_reseed55)
#   Add 5 new seeds per variant -> n=10 total per variant
#
# MD: 5ns (2.5M steps x 2fs), restrained, CUDA
# Snap: n=25 per replicate
# MM-GBSA: v0.6.7 3-traj (NO UPDD_MMGBSA_PROTOCOL env — default 3-traj)
# Baseline comparison: no Fix 4
#
# Safety: skip if mmgbsa_summary.json exists (idempotent)
#         Pass 1b reseed on MD failure (seed+13)
#         GPU budget halt: abort if >10h elapsed
#         PARTIAL_SUCCESS: preserve artifacts (do NOT rm mmgbsa_summary.json or DCD)

set -u
PY=/home/san/miniconda3/envs/md_simulation/bin/python
export PATH=/home/san/miniconda3/envs/md_simulation/bin:$PATH
cd /home/san/UPDD_proj

# New seeds non-colliding with existing {7, 19, 23, 42, 55, 101} (Cp4 reseed55)
# and documented used set {42, 7, 19, 23, 101, 55, 114, 36, 47}
SEEDS=(83 127 163 199 251)
N_SNAP=25
LAUNCH_TS=$(date +%s)
GPU_BUDGET_H=10

log_status() {
    echo "[p1_1_n10 $(date +%H:%M)] $*"
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
        > "/tmp/p1_1_n10_${tag}.log" 2>&1

    # Pass 1b-style seed retry on MD failure
    if ! grep -qE "성공 \((ncAA|야생형) MD\) : 1 개|부분 성공 \(Recovery\) : 1 개|PARTIAL_SUCCESS" "/tmp/p1_1_n10_${tag}.log" 2>/dev/null; then
        local reseed=$((seed+13))
        log_status "WARN ${tag} seed=${seed} fail -> seed=${reseed} retry"
        # Preserve forensic trace but clear incomplete MD artifacts
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
            > "/tmp/p1_1_n10_${tag}_reseed.log" 2>&1
    fi
    if ! grep -qE "성공 \((ncAA|야생형) MD\) : 1 개|최종 구조 저장" "/tmp/p1_1_n10_${tag}.log" "/tmp/p1_1_n10_${tag}_reseed.log" 2>/dev/null; then
        log_status "FAIL ${tag} 2 seeds both exploded - skipping snap/MM-GBSA (expected noise for Cp4)"
        return 1
    fi

    # Snap (n=25) + MM-GBSA (v0.6.7 3-traj default — NO UPDD_MMGBSA_PROTOCOL)
    log_status "snap+MMGBSA ${tag} (n=${N_SNAP}, 3-traj)..."
    $PY utils/extract_snapshots.py \
        --md_dir "${outdir}/mdresult" --outputdir "${outdir}/snapshots" \
        --n_snapshots "$N_SNAP" --binder_chain B --target_id "$target_id" \
        >> "/tmp/p1_1_n10_${tag}.log" 2>&1
    UPDD_MMGBSA_PLATFORM=CUDA $PY utils/run_mmgbsa.py \
        --md_dir "${outdir}/snapshots" --outputdir "${outdir}/mmgbsa_results" \
        --ncaa_elem none --receptor_chain A --binder_chain B --target_id "$target_id" \
        >> "/tmp/p1_1_n10_${tag}.log" 2>&1
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
if [ ! -f "/tmp/calib_params/MTR_params_manifest.json" ]; then
    log_status "ERR MTR ncAA params cache missing at /tmp/calib_params/ - regenerate first"
    exit 1
fi
log_status "params cache OK: MTR_params_manifest.json present"
log_status "launch: seeds=${SEEDS[*]} n_snap=${N_SNAP} budget=${GPU_BUDGET_H}h"

# ============================================================
# Scope: 2QKI_WT + 2QKI_Cp4 only (per Rev. 5 P1.1 dispatch)
#   target_id | topology | wt_pdb | ncaa_code | xml_resname | resid | variant_label
# ============================================================
# Source WT PDB verified md5-identical across existing s42/s101 runs and
# outputs/2QKI_replicate_seed42/_md_input/2QKI_compstatin_4W9A.pdb
WT_PDB="outputs/2QKI_replicate_seed42/_md_input/2QKI_compstatin_4W9A.pdb"

# --- 2QKI_WT runs (5 new seeds) ---
log_status "=== 2QKI_WT launch (5 seeds) ==="
for seed in "${SEEDS[@]}"; do
    check_budget || { log_status "BUDGET HALT before 2QKI_WT_s${seed}"; break; }
    D="outputs/2QKI_WT_calib_s${seed}"
    if [ -f "${D}/mmgbsa_results/mmgbsa_summary.json" ]; then
        log_status "SKIP 2QKI_WT_s${seed} (summary exists)"
        continue
    fi
    mkdir -p "${D}/_md_input" "${D}/mdresult" "${D}/snapshots" "${D}/mmgbsa_results"
    cp "$WT_PDB" "${D}/_md_input/2QKI_WT.pdb"
    run_md_full "$D" "${D}/_md_input/2QKI_WT.pdb" "cyclic_ss" "none" "" "2QKI" "$seed" ""
done

# --- 2QKI_Cp4 runs (5 new seeds) ---
log_status "=== 2QKI_Cp4 launch (5 seeds) ==="
for seed in "${SEEDS[@]}"; do
    check_budget || { log_status "BUDGET HALT before 2QKI_Cp4_s${seed}"; break; }
    D="outputs/2QKI_Cp4_calib_s${seed}"
    if [ -f "${D}/mmgbsa_results/mmgbsa_summary.json" ]; then
        log_status "SKIP 2QKI_Cp4_s${seed} (summary exists)"
        continue
    fi
    mutate_and_prep "$D" "$WT_PDB" "MTR" "4" "2QKI_Cp4"
    run_md_full "$D" "${D}/_md_input/2QKI_Cp4.pdb" "cyclic_ss" \
                "MTR" "MTR" "2QKI" "$seed" \
                "${D}/_md_input/MTR_params_manifest.json"
done

log_status "========================================================"
log_status "P1.1 n=10 EXPANSION COMPLETE (or budget halted)"
log_status "========================================================"

# Aggregate status JSON across all Cp4/WT seeds (existing + new)
mkdir -p outputs/analysis
$PY << 'PYEOF'
import json, glob, os, math, time

variants = {'2QKI_WT': [], '2QKI_Cp4': []}
for jf in sorted(glob.glob('/home/san/UPDD_proj/outputs/2QKI_WT_calib_s*/mmgbsa_results/mmgbsa_summary.json')):
    tag = os.path.basename(os.path.dirname(os.path.dirname(jf)))
    d = json.load(open(jf))
    ie = list(d.get('interaction_entropy_per_design', {}).values())
    s = ie[0] if ie else {'std_delta_g_kcal': 0.0}
    variants['2QKI_WT'].append({
        'tag': tag, 'mean_dg': d['mean_dg'],
        'sigma_wi': s['std_delta_g_kcal'],
        'n_calc': d.get('n_calc', 0),
    })
for jf in sorted(glob.glob('/home/san/UPDD_proj/outputs/2QKI_Cp4_calib_s*/mmgbsa_results/mmgbsa_summary.json')):
    tag = os.path.basename(os.path.dirname(os.path.dirname(jf)))
    if '_reseed' in tag:
        continue  # exclude reseed sibling from primary aggregation
    d = json.load(open(jf))
    ie = list(d.get('interaction_entropy_per_design', {}).values())
    s = ie[0] if ie else {'std_delta_g_kcal': 0.0}
    variants['2QKI_Cp4'].append({
        'tag': tag, 'mean_dg': d['mean_dg'],
        'sigma_wi': s['std_delta_g_kcal'],
        'n_calc': d.get('n_calc', 0),
    })

summary = {'generated_ts': time.time(), 'variants': {}}
for vname, entries in variants.items():
    if not entries:
        continue
    dgs = [e['mean_dg'] for e in entries]
    sigmas = [e['sigma_wi'] for e in entries]
    n = len(dgs)
    mu = sum(dgs) / n
    var = sum((x - mu) ** 2 for x in dgs) / max(n - 1, 1)
    sigma_btwn = math.sqrt(var)
    sigma_wi_avg = sum(sigmas) / n
    summary['variants'][vname] = {
        'n_reps': n,
        'entries': entries,
        'mean_mean_dg': mu,
        'sigma_between': sigma_btwn,
        'sigma_within_avg': sigma_wi_avg,
    }

# ΔΔG aggregation (Cp4 - WT)
if '2QKI_WT' in summary['variants'] and '2QKI_Cp4' in summary['variants']:
    wt = summary['variants']['2QKI_WT']
    cp4 = summary['variants']['2QKI_Cp4']
    ddg = cp4['mean_mean_dg'] - wt['mean_mean_dg']
    sigma_ddg = math.sqrt(cp4['sigma_between']**2 + wt['sigma_between']**2)
    n_min = min(cp4['n_reps'], wt['n_reps'])
    sem_ddg = sigma_ddg / math.sqrt(n_min) if n_min > 0 else float('nan')
    t_crit_df9 = 2.262  # df=9, 95% two-tailed
    ci95_half = t_crit_df9 * sem_ddg
    lit_ddg = -3.96
    z_lit = abs(ddg - lit_ddg) / sigma_ddg if sigma_ddg > 0 else float('nan')
    # Classification: X.A (z<1), X.A' (1<=z<2), X.B (z>=2)
    if z_lit < 1.0:
        tier = 'X.A'
    elif z_lit < 2.0:
        tier = "X.A'"
    else:
        tier = 'X.B'
    summary['ddg_aggregate'] = {
        'ddg_n10': ddg,
        'sigma_ddg_n10': sigma_ddg,
        'sem_ddg_n10': sem_ddg,
        'ci95_half': ci95_half,
        'ci95_lower': ddg - ci95_half,
        'ci95_upper': ddg + ci95_half,
        'n_min': n_min,
        'literature_ddg': lit_ddg,
        'z_to_literature': z_lit,
        'ci95_contains_literature': (ddg - ci95_half) <= lit_ddg <= (ddg + ci95_half),
        'tier_classification': tier,
    }

out_path = '/home/san/UPDD_proj/outputs/analysis/p1_1_n10_expansion_status_20260422.json'
with open(out_path, 'w') as f:
    json.dump(summary, f, indent=2)

print(f'[p1_1_n10] aggregation saved: {out_path}')
for vname, v in summary['variants'].items():
    print(f'  {vname}: n={v["n_reps"]} <<dG>>={v["mean_mean_dg"]:+.2f} sigma_btwn={v["sigma_between"]:.2f} sigma_wi_avg={v["sigma_within_avg"]:.2f}')
if 'ddg_aggregate' in summary:
    d = summary['ddg_aggregate']
    print(f'  DDG_n{d["n_min"]} = {d["ddg_n10"]:+.2f} +/- {d["sigma_ddg_n10"]:.2f} (SEM={d["sem_ddg_n10"]:.2f})')
    print(f'  CI95 = [{d["ci95_lower"]:+.2f}, {d["ci95_upper"]:+.2f}]  contains lit(-3.96)={d["ci95_contains_literature"]}')
    print(f'  z_to_lit = {d["z_to_literature"]:.2f}  tier = {d["tier_classification"]}')
PYEOF
