#!/usr/bin/env bash
set -uo pipefail

ROOT=/home/san/UPDD_proj/outputs/_trackb/w4a_postdensify_dcd_20260716
OUTDIR=${ROOT}/hourly_monitor
LATEST=${OUTDIR}/latest.md
ALL=${OUTDIR}/reports.md
RUN_SESSION=w4a_postdcd_20260716
ATM_PYTHON=/home/san/miniconda3/envs/atm/bin/python
GATE_SCRIPT=/home/san/UPDD_proj/analysis/w4a_postdensify_dcd_20260716/postdensify_gate.py
GATE_LOG=${ROOT}/postdensify_gate.log
GATE_JSON=/home/san/UPDD_proj/analysis/w4a_postdensify_dcd_20260716/postdensify_gate.json
GATE_LOCK=${ROOT}/postdensify_gate.lock

run_once() {
  local ts stamp snap tmp
  ts=$(date -Is)
  stamp=$(date +%Y%m%d_%H%M%S)
  mkdir -p "${OUTDIR}"
  snap=${OUTDIR}/report_${stamp}.md
  tmp=${snap}.tmp

  {
    echo "# W4A Post-densify DCD Monitor"
    echo
    echo "- timestamp: ${ts}"
    if tmux has-session -t "${RUN_SESSION}" 2>/dev/null; then
      echo "- runner tmux: alive"
    else
      echo "- runner tmux: missing"
    fi
    echo
    echo "## Run State"
    if [ -f "${ROOT}/run_state.json" ]; then
      jq . "${ROOT}/run_state.json"
    else
      echo "run_state.json not written yet"
    fi
    echo
    echo "## Progress"
    if [ -f "${ROOT}/progress.json" ]; then
      jq '{status,n_completed,n_expected,median_elapsed_s,eta_s,pending_cell_ids}' "${ROOT}/progress.json"
    else
      echo "progress.json not written yet"
    fi
    echo
    echo "## Process"
    pgrep -af '[p]ostdensify_dcd.py' || echo "no sidecar process found"
    echo
    echo "## GPU"
    nvidia-smi --query-gpu=timestamp,name,memory.used,utilization.gpu,temperature.gpu \
      --format=csv,noheader
    nvidia-smi --query-compute-apps=pid,process_name,used_memory \
      --format=csv,noheader,nounits 2>/dev/null || true
    echo
    echo "## Host Memory"
    free -h
    echo
    echo "## Failure Scan"
    if [ -f "${ROOT}/postdensify_dcd.log" ]; then
      rg -n 'Traceback|RuntimeError|ValueError|NaN|CUDA_ERROR|FAILED|Killed|Out Of Memory' \
        "${ROOT}/postdensify_dcd.log" || echo "no failure pattern found"
      echo
      echo "## Log Tail"
      tail -n 80 "${ROOT}/postdensify_dcd.log"
    else
      echo "runner log not written yet"
    fi
  } > "${tmp}"

  mv "${tmp}" "${snap}"
  cp "${snap}" "${LATEST}"
  {
    echo
    echo "<!-- ${ts} -->"
    cat "${snap}"
  } >> "${ALL}"
}

if [ "${1:-}" = "--loop" ]; then
  while true; do
    run_once
    state=$(jq -r '.status // "UNKNOWN"' "${ROOT}/run_state.json" 2>/dev/null || echo UNKNOWN)
    if [ "${state}" = "COMPLETE" ]; then
      exec 9>"${GATE_LOCK}"
      flock 9
      if [ ! -f "${GATE_JSON}" ]; then
        "${ATM_PYTHON}" "${GATE_SCRIPT}" > "${GATE_LOG}" 2>&1
      fi
      flock -u 9
      exit 0
    fi
    if [ "${state}" = "FAILED" ]; then
      exit 0
    fi
    sleep 3600
  done
else
  run_once
fi
