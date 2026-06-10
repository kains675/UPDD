#!/usr/bin/env bash
# DEPRECATED 2026-05-31 — Track B v1 launcher architecturally invalid
# (single-context direction-flip NaN). v1 archived under
# outputs/_trackb/_archive/v1_architectural_failed_20260531/.
#
# Resume against v1 .out files is forbidden: forward-only data violates
# the UWHAM bidirectional requirement.
#
# USE INSTEAD: scripts/trackb_production_v2_asyncre.py (upstream
# atom_openmm rbfe_structprep + rbfe_production, separate replica walkers
# per direction = no Context flip).
#
# ================================================================
# Original header preserved for historical context only:
#
# Track B production v1 resume launcher — after r11 NaN crash at State 11.
#
# Crash 2026-05-30T23:57: cp4/bound/r11 (state 11, λ=0.5, d=-1) NaN at first
# integrator step. The State 10→11 transition flips the ATMForce Direction
# parameter sign — the potential-gradient discontinuity overwhelms the 2fs
# timestep before the first sample writes.
#
# Engineering fix (in scripts/trackb_production.py): per-state direction-flip
# warmup (0.1ps @ 1fs) before sampling. Resume support: --resume-from-state N
# fast-forwards through completed states with no .out writes.
#
# Strategy:
#   1) cp4/bound: resume State 11..21 (preserves the 11 good r0..r10 files)
#   2) cp4/free + wt/bound + wt/free: full 22-state schedule from scratch
#
# Both invocations use the same v1 output root so uwham can analyze
# everything together when all legs finish.
#
# Hardware: host 5070Ti (GPU 0). V100 is occupied by Track A on the VM.
#
# Crash-state evidence archived (preserved, not deleted) to
# outputs/_trackb/production_v1/_archive/crash_state11_20260530T2357/

set -euo pipefail

cd /home/san/UPDD_proj

ROOT_OUT="outputs/_trackb/production_v1"
LOG_DIR="${ROOT_OUT}"
TS="$(date +%Y%m%dT%H%M%S)"
META_LOG="${LOG_DIR}/launch_resume_${TS}.log"
PY=/home/san/miniconda3/envs/atm/bin/python

# Common args (identical to v1 launch on the surface — only resume flags differ).
COMMON_ARGS=(
  --seed-tag s7
  --out-root "${ROOT_OUT}"
  --endpoints cp4,wt
  --legs bound,free
  --n-replicates 1
  --prod-ns 1.0
  --equil-ps 500
  --sample-ps 50
  --displacement-nm 5.0,0,0
  --timestep-fs 2.0
  --base-seed 20260530
  --direction-flip-warmup-ps 1.0
  --direction-flip-warmup-fs 1.0
  --prod-timestep-fs-after-flip 1.0
)

{
  echo "Resume launch ${TS}"
  echo "=================================================="
  echo "Step 1/2: cp4/bound resume from State 11 (skip State 0 equilibration —"
  echo "          system was equilibrated for ~6 min + 10 states of MD in v1)"
  echo
  CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=8 "${PY}" scripts/trackb_production.py \
    "${COMMON_ARGS[@]}" \
    --resume-endpoints cp4 --resume-legs bound \
    --resume-from-state 11 --resume-skip-equilibration

  echo
  echo "=================================================="
  echo "Step 2/2: cp4/free + wt/bound + wt/free from State 0 (full schedule)"
  echo
  CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=8 "${PY}" scripts/trackb_production.py \
    "${COMMON_ARGS[@]}" \
    --resume-endpoints cp4 --resume-legs free

  CUDA_VISIBLE_DEVICES=0 OMP_NUM_THREADS=8 "${PY}" scripts/trackb_production.py \
    "${COMMON_ARGS[@]}" \
    --resume-endpoints wt --resume-legs bound,free

  echo
  echo "=================================================="
  echo "Resume launch ${TS} DONE"
} 2>&1 | tee -a "${META_LOG}"
