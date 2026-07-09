#!/bin/bash
# Complete the FF / DFT / SAPT2+(3) triangulation WITHOUT the disk-prohibitive
# gold FNO-CCSD(T) (deferred to a 2TB external SSD run):
#   - tier-1 DFT (wb97xd): small scratch, defchararray-shim-fixed, safe on NVMe
#   - tier-2 SAPT (sapt0 + sapt2p3): moderate scratch, no CCSD(T) (vv|vv) blowup
# qm.json merge preserves the already-computed tier-1 silver (sapt0/sapt2p3).
# Host-side scratch on root NVMe (351G free); neither method explodes it.
cd /home/san/UPDD_proj || exit 1
export SAPT_PSI_SCRATCH=/home/san/psi4_scratch_gold
mkdir -p "$SAPT_PSI_SCRATCH"
LOG=logs/sapt_triangulation.log
PY=/home/san/miniconda3/envs/psi4/bin/python
echo "[$(date -u +%FT%TZ)] START tier-1 DFT wb97xd (skip gold; scratch=$SAPT_PSI_SCRATCH, /free=$(df -h / | awk 'NR==2{print $4}'))" >> "$LOG"
taskset -c 0-15 $PY scripts/sapt_dimer_benchmark.py \
  --variants WT MTR W4A --tiers 1 --qm-methods wb97xd --qm-only >> "$LOG" 2>&1
echo "[$(date -u +%FT%TZ)] tier-1 DFT done rc=$?; START tier-2 SAPT" >> "$LOG"
taskset -c 0-15 $PY scripts/sapt_dimer_benchmark.py \
  --variants WT MTR W4A --tiers 2 --qm-only >> "$LOG" 2>&1
echo "[$(date -u +%FT%TZ)] tier-2 done rc=$?; TRIANGULATION COMPLETE (gold deferred to external SSD)" >> "$LOG"
