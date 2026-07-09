#!/bin/bash
cd /home/san/UPDD_proj || exit 1
LOG=logs/sapt_gold_extssd.log
echo "[$(date -u +%FT%TZ)] watcher: tier-2 삼각표(PID 2785733) 종료 대기" >> "$LOG"
while kill -0 2785733 2>/dev/null; do sleep 120; done
echo "[$(date -u +%FT%TZ)] tier-2 종료 감지; 30s RAM 해제 후 gold(ExtSSD, full cores) 시작" >> "$LOG"
sleep 30
bash scripts/sapt_gold_extssd.sh
