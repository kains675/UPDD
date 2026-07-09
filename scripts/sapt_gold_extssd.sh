#!/bin/bash
# Gold FNO-CCSD(T) tier-1 on the 2TB external SSD (/media/san/ExtSSD, 1.8T free,
# USB-C SSD) — huge enough for the ~400-600G CCSD(T) scratch that overflowed the
# internal disks. Runs SEQUENTIALLY after the tier-2 SAPT triangulation finishes
# (a concurrent run risks a 46+46 GB memory collision on 60 GB RAM), so it uses
# all cores (0-15). qm.json merge preserves the silver + DFT; adds gold only.
cd /home/san/UPDD_proj || exit 1
export PYTHONPATH=/home/san/UPDD_proj/_env_shim:$PYTHONPATH
export SAPT_PSI_SCRATCH=/media/san/ExtSSD/psi4_scratch_gold
mkdir -p "$SAPT_PSI_SCRATCH"
LOG=logs/sapt_gold_extssd.log
PY=/home/san/miniconda3/envs/psi4/bin/python
echo "[$(date -u +%FT%TZ)] START gold tier-1 fno_ccsd_t on ExtSSD (free=$(df -h /media/san/ExtSSD | awk 'NR==2{print $4}'), cores 0-15, sequential)" >> "$LOG"
taskset -c 0-15 $PY scripts/sapt_dimer_benchmark.py \
  --variants WT MTR W4A --tiers 1 --qm-methods fno_ccsd_t --qm-only >> "$LOG" 2>&1
echo "[$(date -u +%FT%TZ)] gold done rc=$?; GOLD COMPLETE (ExtSSD)" >> "$LOG"
