#!/bin/bash
# Re-run WT + MTR gold. Memory RESTORED 30GB -> 50GB (2026-07-05): the 30GB
# "headroom" reduction was a failed hypothesis; the canonical-orthogonalization fix
# got the light (monomer-type) fragments through cleanly (CCSD(T) totals), but the
# heavy dimer fragment segfaults (NOT OOM per dmesg) at the FNO (OV|OV) integral
# presort after writing ~11e9 two-electron integrals. Larger psi4 memory gives the
# presort bigger in-core sort buffers -> may avoid the pathological path. 50GB on
# 60GB RAM keeps ~10GB OS headroom (sequential single job, no 46+46 collision).
# canonical fix (sapt_dimer_benchmark.py) + defchararray sitecustomize both active.
# NOTE (R-18): if it segfaults again at the same FNO presort, memory was NOT the
# cause -> stop, do not blind-retry, escalate to DF-CCSD(T) path.
cd /home/san/UPDD_proj || exit 1
export PYTHONPATH=/home/san/UPDD_proj/_env_shim:$PYTHONPATH
export SAPT_PSI_SCRATCH=/media/san/ExtSSD/psi4_scratch_gold
export SAPT_PSI_MEM="50 GB"
mkdir -p "$SAPT_PSI_SCRATCH"
LOG=logs/sapt_gold_wtmtr.log
PY=/home/san/miniconda3/envs/psi4/bin/python
echo "[$(date -u +%FT%TZ)] START gold WT+MTR re-run (DF-CCSD(T): cc_type/mp2_type=df -> IWL sort 회피 [psi4 #35 32-bit overflow]; canonical+sitecustomize 유지, mem=50GB, ExtSSD)" >> "$LOG"
taskset -c 0-15 $PY scripts/sapt_dimer_benchmark.py \
  --variants WT MTR --tiers 1 --qm-methods fno_ccsd_t --qm-only >> "$LOG" 2>&1
echo "[$(date -u +%FT%TZ)] gold WT+MTR done rc=$?; DONE" >> "$LOG"
