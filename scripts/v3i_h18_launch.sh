#!/bin/bash
# V3I (Cp01 Val3->Ile3) H18-grade densified re-validation — task #115 / ADR-0014 addendum.
# WHY: existing n=6 V3I (12-state/400cyc, O2_resolved_discordant +1.333) FAILS H18 hist-OVL>=0.45
#      on ALL 12 bound legs AND all 12 free legs (analysis/v3i_h18_regate_20260629/). The discordant
#      +sign may be a 12-state mid-ladder sampling artifact. Re-run BOTH legs at matched H18 grade.
# SCHEDULE (V3I-specific equipartition, curve-decided, 15 states/dir; SciVal GO_WITH_COND):
#   leg-up  lambda2 = 0,0.05,0.1,0.15,0.2,0.3,0.4,0.5   (splits (0,1) leg-up cliff)
#   leg-down lambda1 = 0.025,0.05,0.1,0.2,0.3,0.4,0.5   (splits (5,6) apex cliff)
# C7 soft-core canon FROZEN (untouched). Endpoints (0,0)->(0.5,0.5) invariant => dG-unbiased, ranking-only (R-11).
# wt-only single-endpoint (copy1 Ile vs copy2 Val, one 2QKI WT scaffold). seeds fixed order = rep-index pairing.
# 2-way MAX (--max-concurrent 2): VM 23 GiB RAM + V100 32GB. free (cheap, no DCD) THEN bound (DCD stride2).
set -u
cd /home/san/UPDD_proj || exit 1
PY=/home/san/miniconda3/envs/atm/bin/python
SCRIPT=/home/san/UPDD_proj/scripts/trackb_inplace_rbfe_production.py
OR=/home/san/UPDD_proj/outputs/_trackb/twocopy_v3i_h18_n6
SEEDS=s7,s101,s127,s23,s163,s199
L2='0,0.05,0.1,0.15,0.2,0.3,0.4,0.5'
L1='0.025,0.05,0.1,0.2,0.3,0.4,0.5'
mkdir -p "$OR"

echo "=== V3I H18 re-validation START $(date -u +%FT%TZ) host=$(hostname) ==="

echo "=== FREE leg start $(date -u +%FT%TZ) ==="
$PY $SCRIPT \
  --endpoints wt --seeds $SEEDS --leg free --twocopy --mutation v3i_val_ile_res3 \
  --auto-search-displacement --accept-sep-nm 1.5 \
  --directions dplus,dminus --n-windows-half 6 \
  --lambda2-rampdown "$L2" --lambda1-rampdown "$L1" \
  --n-cycles 800 --md-steps-per-cycle 250 \
  --pool --max-concurrent 2 --device-index 0 --mintimeid 600 \
  --out-root "$OR" > "$OR/driver_free.log" 2>&1
echo "=== FREE leg done rc=$? $(date -u +%FT%TZ) ==="

echo "=== BOUND leg start $(date -u +%FT%TZ) ==="
$PY $SCRIPT \
  --endpoints wt --seeds $SEEDS --leg bound --twocopy --mutation v3i_val_ile_res3 \
  --auto-search-displacement --accept-sep-nm 1.5 \
  --directions dplus,dminus --n-windows-half 6 \
  --lambda2-rampdown "$L2" --lambda1-rampdown "$L1" \
  --dcd --dcd-stride-cycles 2 \
  --n-cycles 800 --md-steps-per-cycle 250 \
  --pool --max-concurrent 2 --device-index 0 --mintimeid 600 \
  --out-root "$OR" > "$OR/driver_bound.log" 2>&1
echo "=== BOUND leg done rc=$? $(date -u +%FT%TZ) ==="
echo "=== ALL V3I H18 DONE $(date -u +%FT%TZ) ==="
