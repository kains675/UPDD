#!/bin/bash
# Relaunch SAPT gold (FNO-CCSD(T)) + DFT (tier-1) and tier-2 SAPT with scratch on
# the root NVMe (334G free, fastest SSD) after the WT-gold disk-full SIGSEGV on the
# 92G SATA SSD. qm.json is now merge-safe, so the tier-1 re-run of only gold+DFT
# preserves the already-computed silver (sapt0/sapt2p3). Host-side scratch run.
cd /home/san/UPDD_proj || exit 1
export SAPT_PSI_SCRATCH=/home/san/psi4_scratch_gold
mkdir -p "$SAPT_PSI_SCRATCH"
LOG=logs/sapt_gold_relaunch.log
PY=/home/san/miniconda3/envs/psi4/bin/python
echo "[$(date -u +%FT%TZ)] START tier-1 gold+DFT (scratch=$SAPT_PSI_SCRATCH, free=$(df -h / | awk 'NR==2{print $4}'))" >> "$LOG"
taskset -c 0-15 $PY scripts/sapt_dimer_benchmark.py \
  --variants WT MTR W4A --tiers 1 --qm-methods fno_ccsd_t wb97xd --qm-only >> "$LOG" 2>&1
echo "[$(date -u +%FT%TZ)] tier-1 done rc=$?; START tier-2 SAPT" >> "$LOG"
taskset -c 0-15 $PY scripts/sapt_dimer_benchmark.py \
  --variants WT MTR W4A --tiers 2 --qm-only >> "$LOG" 2>&1
echo "[$(date -u +%FT%TZ)] tier-2 done rc=$?; ALL COMPLETE" >> "$LOG"
