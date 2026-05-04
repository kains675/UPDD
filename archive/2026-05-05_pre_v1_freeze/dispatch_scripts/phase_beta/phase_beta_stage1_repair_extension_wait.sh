#!/bin/bash
# Wait for the MD reseed dispatch to finish, then run the 2 ns extension MD.
# This script is intended to be launched in parallel; it idles until the master
# log shows "MD reseed phase complete", then proceeds.
set -u
LOGDIR=/home/san/UPDD_proj/outputs/analysis/phase_beta_stage1_repair_run_logs
EXT_LOG="$LOGDIR/extension_2ns.log"
DISPATCH_LOG="$LOGDIR/dispatch.log"

cd /home/san/UPDD_proj

echo "[ext_wait $(date '+%H:%M:%S')] watching $DISPATCH_LOG for completion ..." | tee -a "$LOGDIR/ext_wait.log"

# Poll until completion marker appears (or 14 h cap)
START=$(date +%s)
HARD_CAP=$((14*3600))
while true; do
    if grep -q "MD reseed phase complete" "$DISPATCH_LOG" 2>/dev/null; then
        echo "[ext_wait $(date '+%H:%M:%S')] reseed phase done — launching extension" | tee -a "$LOGDIR/ext_wait.log"
        break
    fi
    NOW=$(date +%s)
    if [ $((NOW - START)) -gt $HARD_CAP ]; then
        echo "[ext_wait $(date '+%H:%M:%S')] hard cap exceeded — bail" | tee -a "$LOGDIR/ext_wait.log"
        exit 1
    fi
    sleep 60
done

# Launch the extension
echo "[ext_wait $(date '+%H:%M:%S')] starting 2 ns extension MD ..." | tee -a "$LOGDIR/ext_wait.log"
/home/san/miniconda3/envs/qmmm/bin/python \
    /home/san/UPDD_proj/scripts/phase_beta_stage1_repair_extension.py \
    > "$EXT_LOG" 2>&1
RC=$?
echo "[ext_wait $(date '+%H:%M:%S')] extension finished, rc=$RC" | tee -a "$LOGDIR/ext_wait.log"
exit $RC
