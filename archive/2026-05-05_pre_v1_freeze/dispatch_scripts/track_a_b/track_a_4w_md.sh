#!/bin/bash
# UPDD X.B Track A Step A2 — 4W (compstatin-native H9) MD 5-rep + Fix 4 MM-GBSA
# Rev. 5 hardening; Route A' v0.6.8 calibration
#
# Input: inputs/4W_from_2QKI_h9_restored.pdb (Coder A1 artifact, md5 fbdb74ee9b...)
# Target: 2QKI, cyclic_ss (Cys2-Cys12), 5ns restrained MD, 25-snap, 1-traj Fix 4
# Seeds: 37, 89, 151, 211, 283 (Pass 1b retry: seed+13)
# Platform: CUDA
#
# Output: outputs/4W_calib_s{seed}/{mdresult,snapshots,mmgbsa_results}
# Log:    /tmp/track_a_4w_md_s{seed}.log

set -u
PY=/home/san/miniconda3/envs/md_simulation/bin/python
export PATH=/home/san/miniconda3/envs/md_simulation/bin:$PATH
cd /home/san/UPDD_proj

SEEDS=(37 89 151 211 283)
N_SNAP=25
TARGET_ID=2QKI
TOPO=cyclic_ss
STEPS=2500000          # 5ns x 2fs
LAUNCH_TS=$(date +%s)
GPU_BUDGET_H=6         # hard cap

log_status() {
    echo "[trackA $(date +%H:%M)] $*"
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
    local seed="$1"
    local D="outputs/4W_calib_s${seed}"
    local tag="4W_s${seed}"
    local logf="/tmp/track_a_4w_md_s${seed}.log"

    if [ -f "${D}/mmgbsa_results/mmgbsa_summary.json" ]; then
        log_status "SKIP ${tag} (summary exists)"
        return 0
    fi
    check_budget || return 2

    log_status "MD ${tag} launch 5ns..."
    UPDD_MMGBSA_PLATFORM=CUDA $PY utils/run_restrained_md.py \
        --inputdir "${D}/_md_input" \
        --outputdir "${D}/mdresult" \
        --steps "$STEPS" --topology "$TOPO" \
        --binder_chain B --graph_policy strict \
        --target_id "$TARGET_ID" --dt_fs 2.0 --platform CUDA --seed "$seed" \
        --ncaa_label none --ncaa_code "" \
        > "$logf" 2>&1

    # Pass 1b-style reseed on failure
    if ! grep -qE "성공 \((ncAA|야생형) MD\) : 1 개|부분 성공 \(Recovery\) : 1 개|PARTIAL_SUCCESS" "$logf" 2>/dev/null; then
        local reseed=$((seed+13))
        log_status "WARN ${tag} seed=${seed} fail -> seed=${reseed} retry"
        rm -f "${D}/mdresult/"*EXPLODED* "${D}/mdresult/"*_partial* \
              "${D}/mdresult/"*_final.pdb "${D}/mdresult/"*_md.log \
              "${D}/mdresult/"*_restrained.dcd 2>/dev/null
        local logf2="/tmp/track_a_4w_md_s${seed}_reseed.log"
        UPDD_MMGBSA_PLATFORM=CUDA $PY utils/run_restrained_md.py \
            --inputdir "${D}/_md_input" \
            --outputdir "${D}/mdresult" \
            --steps "$STEPS" --topology "$TOPO" \
            --binder_chain B --graph_policy strict \
            --target_id "$TARGET_ID" --dt_fs 2.0 --platform CUDA --seed "$reseed" \
            --ncaa_label none --ncaa_code "" \
            > "$logf2" 2>&1
        if ! grep -qE "성공 \((ncAA|야생형) MD\) : 1 개|최종 구조 저장" "$logf" "$logf2" 2>/dev/null; then
            log_status "FAIL ${tag} 2 seeds both failed - skip downstream"
            return 1
        fi
    fi

    # Snap extract (n=25)
    log_status "snap ${tag} (n=${N_SNAP})..."
    UPDD_MMGBSA_PLATFORM=CUDA $PY utils/extract_snapshots.py \
        --md_dir "${D}/mdresult" --outputdir "${D}/snapshots" \
        --n_snapshots "$N_SNAP" --binder_chain B --target_id "$TARGET_ID" \
        >> "$logf" 2>&1

    # MM-GBSA Fix 4 (1-traj protocol)
    log_status "MMGBSA Fix 4 (1-traj) ${tag}..."
    UPDD_MMGBSA_PROTOCOL=1traj UPDD_MMGBSA_PLATFORM=CUDA $PY utils/run_mmgbsa.py \
        --md_dir "${D}/snapshots" --outputdir "${D}/mmgbsa_results" \
        --ncaa_elem none --receptor_chain A --binder_chain B --target_id "$TARGET_ID" \
        >> "$logf" 2>&1

    if [ -f "${D}/mmgbsa_results/mmgbsa_summary.json" ]; then
        $PY -c "
import json
d = json.load(open('${D}/mmgbsa_results/mmgbsa_summary.json'))
ie = list(d.get('interaction_entropy_per_design', {}).values())
s = ie[0] if ie else {'std_delta_g_kcal': 0}
n = d.get('n_calc', 0)
print(f'[DONE ${tag}] <dG>={d[\"mean_dg\"]:+.2f} sigma={s[\"std_delta_g_kcal\"]:.2f} n={n}')"
    else
        log_status "FAIL ${tag} no MM-GBSA summary generated"
        return 1
    fi
    return 0
}

log_status "========================================================"
log_status "Track A Step A2: 4W (H9 restored) 5-rep MD + Fix 4 MM-GBSA"
log_status "Seeds: ${SEEDS[*]}  Target: ${TARGET_ID}  Topo: ${TOPO}"
log_status "========================================================"

for seed in "${SEEDS[@]}"; do
    check_budget || { log_status "BUDGET HALT before s${seed}"; break; }
    run_md_full "$seed"
done

log_status "========================================================"
log_status "Track A 4W complete (or halted)"
log_status "========================================================"

# Aggregate
$PY << 'PYEOF'
import json, glob, os, statistics
results = {}
for seed in (37, 89, 151, 211, 283):
    jf = f'/home/san/UPDD_proj/outputs/4W_calib_s{seed}/mmgbsa_results/mmgbsa_summary.json'
    if not os.path.exists(jf):
        continue
    d = json.load(open(jf))
    ie = list(d.get('interaction_entropy_per_design', {}).values())
    s = ie[0] if ie else {'std_delta_g_kcal': 0}
    results[seed] = {
        'mean_dg': d['mean_dg'],
        'sigma_wi': s['std_delta_g_kcal'],
        'n_calc': d.get('n_calc', 0),
    }

means = [r['mean_dg'] for r in results.values()]
sigs  = [r['sigma_wi'] for r in results.values()]
agg = {
    'generated_ts': __import__('time').time(),
    'n_completed': len(results),
    'per_seed': results,
}
if means:
    agg['mean_of_means'] = sum(means)/len(means)
    agg['sigma_btwn']   = statistics.stdev(means) if len(means) >= 2 else 0.0
    agg['sigma_wi_avg'] = sum(sigs)/len(sigs)
os.makedirs('/home/san/UPDD_proj/outputs/analysis', exist_ok=True)
with open('/home/san/UPDD_proj/outputs/analysis/track_a_4w_fix4_20260423.json', 'w') as f:
    json.dump(agg, f, indent=2)

print(f'[track_a_4w] {len(results)} seeds done')
for seed in sorted(results):
    r = results[seed]
    print(f'  s{seed}  <dG>={r["mean_dg"]:+.3f}  sigma_wi={r["sigma_wi"]:.3f}  n={r["n_calc"]}')
if means:
    print(f'  <<dG>>     = {agg["mean_of_means"]:+.3f}')
    print(f'  sigma_btwn = {agg["sigma_btwn"]:.3f}')
    # 4W9A baseline (WT n=10, P0.3+P1.1) from cp4_fix4_symmetric_n10_20260423.json
    wt_mean = -62.06
    delta = agg["mean_of_means"] - wt_mean
    print(f'  4W - 4W9A  = {delta:+.3f} kcal/mol (H9A expected ~0.3)')
PYEOF

log_status "Aggregate saved to outputs/analysis/track_a_4w_fix4_20260423.json"
