#!/bin/bash
# Host-side auto-launch: wait for STEP-2 Gate-A MD (PID 2662427) to finish,
# then launch the SAPT retry (50GB, aug-cc-pVDZ gold). Runs independent of the
# interactive session. Created 2026-07-03 to avoid tying the retry to a terminal.
cd /home/san/UPDD_proj || exit 1
LOG=logs/sapt_autolaunch.log
echo "[$(date -u +%FT%TZ)] watcher started; waiting for STEP-2 PID 2662427" >> "$LOG"
while kill -0 2662427 2>/dev/null; do sleep 300; done
echo "[$(date -u +%FT%TZ)] STEP-2 PID 2662427 exited; settling 60s then launching SAPT retry" >> "$LOG"
sleep 60   # let GPU MD RAM fully release
PSI_SCRATCH=/media/san/San taskset -c 0-15 nohup /home/san/miniconda3/envs/psi4/bin/python \
  scripts/sapt_dimer_benchmark.py --variants WT MTR W4A --tiers 1 2 --qm-only \
  > logs/sapt_retry_50gb.log 2>&1 </dev/null &
echo "[$(date -u +%FT%TZ)] SAPT retry launched PID $! (50GB / aug-cc-pVDZ gold / cores 0-15)" >> "$LOG"
