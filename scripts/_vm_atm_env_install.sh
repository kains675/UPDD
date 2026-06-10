#!/bin/bash
# scripts/_vm_atm_env_install.sh
#
# VM-side `atm` conda env provisioning script (F3 fix, 2026-05-31).
#
# Purpose: Track B per-direction production launcher
# (`scripts/trackb_per_direction_production.py`) defaults
# `--gpu-host vm` and `--abfe-bin
# /home/san/miniconda3/envs/atm/bin/abfe_production`. The VM (192.168.122.155)
# initially had only `qmmm` + `md_simulation` envs, so the default abfe_bin
# path failed with rc=127. This script provisions a minimal `atm` env on
# the VM matching the host atm spec (openmm 8.5.1 + atom-openmm 8.4.0).
#
# Usage (from host):
#
#     scp scripts/_vm_atm_env_install.sh san@192.168.122.155:/tmp/
#     ssh san@192.168.122.155 'bash /tmp/_vm_atm_env_install.sh \
#         > /tmp/vm_atm_install.log 2>&1 &'
#
# Verification (host-side):
#
#     ssh san@192.168.122.155 \
#         "/home/san/miniconda3/envs/atm/bin/abfe_production --help"
#     # rc=0, prints abfe_production banner
#
# Safety: does NOT touch the existing `qmmm` env (Track A QM batch
# uses it; PID 1717987-1896877 series). Creates `atm` from scratch.
# Use mamba (already on VM at conda 26.3.2 / mamba 2.6.2) for fast solve.

set -e
source /home/san/miniconda3/etc/profile.d/conda.sh

echo "=== STEP 0: Pre-check existing atm env ==="
if conda env list | awk '{print $1}' | grep -qx atm; then
    echo "WARN: atm env already exists. Listing contents:"
    conda list -n atm 2>&1 | head -20
    echo "Will NOT re-create; will try install into existing env."
else
    echo "atm env does not exist. Will create."
    echo ""
    echo "=== STEP 1: Create atm env (python 3.11) ==="
    mamba create -n atm -c conda-forge python=3.11 -y
fi

echo ""
echo "=== STEP 2: Install openmm 8.5.1 + scientific stack via mamba ==="
# Host atm spec (reference):
#   openmm                         8.5.1     conda-forge
#   mdtraj                         1.11.1    conda-forge
#   numpy                          2.4.6     conda-forge
#   scipy                          1.17.1    conda-forge
#   openmmtools                    0.26.0    conda-forge
mamba install -n atm -c conda-forge \
    openmm=8.5.1 \
    mdtraj \
    numpy \
    scipy \
    openmmtools \
    -y

echo ""
echo "=== STEP 3: Install atom-openmm via pip ==="
# atom-openmm is pip-only; pin to 8.4.0 to match host.
conda run -n atm pip install --upgrade pip
conda run -n atm pip install atom-openmm==8.4.0

echo ""
echo "=== STEP 4: Verify ==="
echo "--- python version ---"
conda run -n atm python --version
echo "--- openmm ---"
conda run -n atm python -c 'import openmm; print(openmm.version.full_version)'
echo "--- atom_openmm import ---"
conda run -n atm python -c 'import atom_openmm; print(getattr(atom_openmm, "__version__", "no __version__ attr"))'
echo "--- abfe_production --help ---"
/home/san/miniconda3/envs/atm/bin/abfe_production --help 2>&1 | head -5 || \
    echo "WARN abfe_production --help nonzero rc"
echo ""
echo "=== INSTALL COMPLETE ==="
