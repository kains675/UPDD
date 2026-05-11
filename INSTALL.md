# UPDD — Installation & Reproducibility Guide

This document is the **step-by-step setup guide** for external users who want to reproduce, validate, or extend the UPDD pipeline. For a quick orientation see `README.md`; for the release notes see `CHANGELOG.md`; for the full methodology see Paper 1 v1 (ChemRxiv link in `README.md` Citation section).

The setup proceeds in three tiers depending on what you intend to do:

- **Tier A — Code review / unit tests / schema regression** (no GPU, ~15 min, ~3 GB disk). Sufficient for inspecting the R-15 / R-16 charge guards, R-17 cofactor preservation gates, branched ΔΔG schema, σ_btwn / σ_w decomposition statistics, and pytest suite. Does NOT cover MD, QM/MM, or MM-PBSA scoring.
- **Tier B — Full pipeline reproduction on a single workstation** (1× NVIDIA GPU with ≥ 16 GB VRAM, CUDA 12.x, ~50 GB disk). Reproduces every stage of the published Paper 1 v1 results from the six-target ncAA benchmark.
- **Tier C — Docker / containerised CI path** (CPU-only, runs on GitHub Actions or any Docker host). Mirrors Tier A inside a container.

---

## 1. System requirements

### 1.1 Tier A — CPU-only (code review, schema regression)

- **OS**: Linux (Ubuntu 22.04 / 24.04 tested), macOS 13+ should work
- **CPU**: any x86-64 (no GPU required)
- **RAM**: 8 GB minimum, 16 GB recommended
- **Disk**: ~3 GB for the conda environment + repo
- **Python**: 3.10 (managed via conda; do not use system Python)

### 1.2 Tier B — Full pipeline (MD + QM/MM + MM-PBSA)

Adds on top of Tier A:

- **GPU**: NVIDIA RTX 5070 Ti 16 GB (tested) or equivalent (RTX 4090, RTX 5090 OK; H100 OK; older Pascal-era GPUs without `cuda_arch >= 70` may fail on `gpu4pyscf`)
- **NVIDIA driver**: ≥ 525.x
- **CUDA Toolkit**: 12.x (system-installed; `cupy-cuda12x` and `gpu4pyscf` link against it)
- **VRAM**: 16 GB minimum (matches the production target; less may force smaller QM partitions and prevent reproduction of the published n_qm = 466 → 600 budget regime)
- **RAM**: 32 GB recommended (one MD replicate at the published box size uses ~6 GB; running 4 in parallel uses ~24 GB)
- **Disk**: ~50 GB for the env + repo + the six-target snapshot caches that ship via Zenodo

### 1.3 Hardware accessibility note

Most existing peptide-design pipelines (RFdiffusion, ProteinDJ, RFpeptides) are benchmarked on A100 / A30 HPC clusters. UPDD is engineered to run on hardware accessible to independent researchers — a single 16 GB NVIDIA GPU and a consumer-class workstation. See `README.md` § "Why UPDD?" for the hardware-accessibility rationale.

---

## 2. Tier A setup — conda environment (≈ 15 min)

### 2.1 Install Miniconda (skip if you already have conda or mamba)

```bash
# Linux x86-64
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
source $HOME/miniconda3/bin/activate
conda init bash    # or `conda init zsh` if you use zsh
# reopen your shell, or `source ~/.bashrc`
```

### 2.2 Clone the repo

```bash
git clone https://github.com/kains675/UPDD.git
cd UPDD
```

### 2.3 Create the `qmmm` environment from `environment.yml`

```bash
conda env create -f environment.yml
conda activate qmmm
```

This installs the Tier-1 packages pinned to the production versions documented in `environment.yml`:

- `python 3.10` (the canonical interpreter for the codebase)
- `openmm 8.4` (MD engine, `CustomGBForce` GBn2 in `utils/run_mmgbsa.py`)
- `pyscf 2.12` + `gpu4pyscf 1.6` (QM/MM SCF, `utils/run_qmmm.py` with ωB97X-D3 / def2-SVP)
- `mdtraj 1.10` (snapshot extraction + PBC unwrap, `utils/extract_snapshots.py`)
- `parmed 4.3` (AmberTools prmtop / mbondi2 radius override)
- `pdbfixer 1.12` (target preprocessing)
- `rdkit 2025.09` (ncAA SMILES → 3D, `utils/parameterize_ncaa.py`)
- `ambertools 24.8` (`MMPBSA.py` invocation from `scripts/run_mmpbsa.py`)
- `cupy-cuda12x 14.0` (gpu4pyscf dependency; can be skipped for CPU-only Tier A)

### 2.4 Verify the environment

```bash
# 1. Confirm the canonical interpreter path
which python
# Expect: /home/<user>/miniconda3/envs/qmmm/bin/python

# 2. Confirm the key Tier-1 imports resolve
python -c "import openmm, mdtraj, pyscf, parmed, rdkit; print('OK')"

# 3. Run the pytest baseline
python -m pytest tests/ -q
# Expected: 294 pass / 9 skip / 1 fail
# (The 1 fail is `test_admet_filter::test_lipinski_constants_are_module_level` —
#  a pre-existing MW = 500 vs 1200 semantic mismatch, separate triage.
#  All other 294 tests passing confirms a clean Tier A install.)
```

> ⚠️ **Common pitfall**: do NOT invoke `pytest` directly from outside the env. The tests import `openmm`, `mdtraj`, and `numpy` from the `qmmm` env. Running `pytest` via system Python or via `~/miniconda3/bin/python` (root miniconda, no `openmm`) silently fails 13+ tests with `ModuleNotFoundError`. The canonical invocation is `python -m pytest tests/` with the `qmmm` env active, as shown above.

### 2.5 Verify charge-guard infrastructure (R-15 / R-16)

```bash
# Audit production ncAA force-field XMLs (read-only, no compute)
python -m pytest tests/test_charge_topology.py tests/test_audit_charge_consistency.py -v
# Expected: 23/23 pass (R-15 magnitude + R-16 binder chemistry SSOT)
```

---

## 3. Tier B setup — full pipeline (additional ≈ 30 min for GPU stack)

### 3.1 Install NVIDIA driver + CUDA Toolkit 12.x

Follow the official NVIDIA driver installation for your distribution. Verify:

```bash
nvidia-smi
# Expect: driver version ≥ 525.x, CUDA version ≥ 12.x in the top-right
```

### 3.2 Verify gpu4pyscf + cupy resolve and detect the GPU

```bash
conda activate qmmm
python -c "import cupy; print('cupy:', cupy.__version__); print('GPU:', cupy.cuda.runtime.getDeviceCount())"
# Expect: cupy: 14.0.x, GPU: 1 (or more)

python -c "import gpu4pyscf; print('gpu4pyscf:', gpu4pyscf.__version__)"
# Expect: gpu4pyscf: 1.6.x
```

### 3.3 Configure AmberTools paths (auto-detected from the qmmm env)

`scripts/run_mmpbsa.py` requires `tleap` and `MMPBSA.py` on PATH. The conda env's bin directory is the canonical location:

```bash
# Set AMBERHOME and ensure PATH covers tleap / MMPBSA.py
export AMBERHOME=$CONDA_PREFIX
export PATH=$AMBERHOME/bin:$PATH

# Verify
which tleap
# Expect: /home/<user>/miniconda3/envs/qmmm/bin/tleap

which MMPBSA.py
# Expect: /home/<user>/miniconda3/envs/qmmm/bin/MMPBSA.py
```

You can persist these by adding the two `export` lines to your `~/.bashrc` (or a project-specific `.envrc` if you use direnv).

### 3.4 Smoke test — end-to-end on a single snapshot

```bash
# This runs one MM-PBSA snapshot under the published Cp4 protocol.
# Takes ~2 min on RTX 5070 Ti, ~5 min on slower GPUs.
python scripts/run_mmpbsa.py \
    --md_dir outputs/2QKI_Cp4_calib_s7/snapshots_n25_postl387_patch_v2 \
    --outputdir /tmp/updd_smoke \
    --target_id 2QKI \
    --binder_chain B \
    --receptor_chain A \
    --max_snaps 1

# Expected last line:
#   ΔG_bind(PBSA) = -10.X kcal/mol  ✅  (matches the published Cp4 s7 / snap01 anchor)
```

### 3.5 Download the published benchmark snapshots (Zenodo)

The full snapshot caches for the six target families (Cp4, MTR13, NML20, NML22, MTR25, MTR6 — 13 seeds × 25 snapshots each) are archived at Zenodo `10.5281/zenodo.20067323`. Download and place under `outputs/` to reproduce the published Table 1 numbers.

---

## 4. Tier C setup — Docker (CPU-only, CI / reviewer path)

### 4.1 Build the image

```bash
git clone https://github.com/kains675/UPDD.git
cd UPDD
docker build -t updd:v0.7.1 .
```

### 4.2 Run the test suite inside the container

```bash
docker run --rm updd:v0.7.1 python -m pytest tests/ -q
# Expected: 294 pass / 9 skip / 1 fail
```

### 4.3 Inspect the image interactively

```bash
docker run --rm -it updd:v0.7.1 bash
# Inside container:
#   python -m pytest tests/ -v -k "charge or audit"
#   cat utils/parameterize_ncaa.py | grep -A 5 "UPDD_MTR_AMBER14_PATCH"
```

The image is CPU-only and does NOT include `gpu4pyscf` / `cupy` (the GPU stack is too large for general-purpose Docker distribution and requires NVIDIA-Container-Toolkit at host-side). The GPU production path uses the local conda env from Tier B above.

---

## 5. Project layout (where to find what)

```
UPDD/
├── README.md                    # Quick orientation + citation
├── INSTALL.md                   # This file
├── CHANGELOG.md                 # Release notes
├── LICENSE                      # MIT
├── environment.yml              # conda env spec (Tier-1 packages)
├── Dockerfile                   # Tier C container build
├── .github/workflows/ci.yml     # GitHub Actions CI
├── utils/                       # Stage-4 evaluation code (read-only protected files)
│   ├── parameterize_ncaa.py     # ncAA FF parameterisation + amb14SB Trp overlay
│   ├── run_qmmm.py              # QM/MM SCF (gpu4pyscf, wB97X-D3 / def2-SVP)
│   ├── run_mmgbsa.py            # MM-GBSA (1-traj, GBn2)
│   ├── run_restrained_md.py     # Restrained MD (OpenMM, ff14SB on protein)
│   ├── ncaa_registry.py         # 25-entry ncAA registry (MTR, NML, MLE, OMW, …)
│   ├── amber_charges.py         # AMBER charge utilities
│   ├── branched_ddg.py          # Paired (variant, WT) ΔΔG engine + Light D Hybrid
│   ├── extract_snapshots.py     # K-Means n=25 + L387 v54 four-layer PBC repair
│   └── audit_charge_consistency.py  # R-15 / R-16 / R-17 audit module
├── scripts/                     # Stage-4 runners + analysis scripts
│   ├── run_mmpbsa.py            # MM-PBSA (AmberTools MMPBSA.py)
│   ├── simulated_al_cycle.py    # Simulated active-learning cycle (synthetic oracle)
│   ├── verify_pbc_integrity.py  # 70-system PBC scan
│   └── ...
├── target_cards/                # Per-target metadata (binder/receptor chain IDs, cofactors, K_d anchors)
├── outputs/                     # Benchmark seeds + MD trajectories + MM-PBSA results
│   ├── 2QKI_Cp4_calib_s*/       # Compstatin Cp4 (8 seeds)
│   ├── 1EBP_MTR13_calib_s*/     # EPO 1EBP MTR13 variants (5 seeds)
│   ├── ...                      # 6 target families total
│   └── _archive/pre_amb14_patch_20260427/   # Pre-patch (legacy) audit baseline
└── tests/                       # 294 pass / 9 skip / 1 pre-existing fail
```

---

## 6. Common issues & troubleshooting

### 6.1 `ModuleNotFoundError: No module named 'openmm'` on pytest

You are running `pytest` from outside the `qmmm` env. Fix:

```bash
conda activate qmmm
python -m pytest tests/
```

### 6.2 `RuntimeError: AMBERHOME not set` on `scripts/run_mmpbsa.py`

You did not export `AMBERHOME` and `PATH`. Fix per §3.3 above. To make this permanent inside the qmmm env, you can use:

```bash
conda env config vars set AMBERHOME=$CONDA_PREFIX -n qmmm
conda activate qmmm   # re-activate to apply
```

### 6.3 `cupy.cuda.runtime.CUDARuntimeError: cudaErrorNoDevice`

Either the NVIDIA driver is not installed, or you are inside a container without GPU passthrough. Verify with `nvidia-smi` on the host first. For Docker, add `--gpus all` to the `docker run` invocation. The CPU-only Tier C path skips the GPU stack entirely.

### 6.4 `pytest test_verify_pbc_integrity.py` fails with `ImportError: arrow_menu`

This is a known regression caused by `sys.path` pollution between `scripts/updd_cli.py` and `utils/updd_tui.py` (legacy `updd_cli.py` shadow). The fix is in v0.7.1 — make sure your local copy is at `fe4309b` or later. See CHANGELOG `[v0.7.1]` → "Code unification (Phase 1+2)".

### 6.5 `gpu4pyscf` SCF non-convergence on small VRAM

UPDD targets 16 GB VRAM (RTX 5070 Ti). If you have less, the QM/MM partition budget `max_n_qm` (default 600 atoms) may exceed your VRAM. See `target_cards/<system>.json` → `"max_n_qm"` and reduce as needed; expect slightly different n_qm partitioning per system.

### 6.6 MM-PBSA scaling on multiple seeds

Multiple `scripts/run_mmpbsa.py` invocations on the same GPU share GPU compute and CUDA context. Empirically, 4-way parallel on RTX 5070 Ti runs at ~25 % per-process throughput (per-snap timing increases ~4× under 4-way contention). For batch sweeps, allocate one GPU per process or use a job scheduler.

### 6.7 R-15 / R-16 audit blocking your run

The charge consistency audit (`utils/audit_charge_consistency.py`) raises `ChargeDeclarationMismatch` or `BinderChargeMismatch` when the runtime-computed charge differs from `target_card.target_iso_net_charge` or `target_card.binder_net_charge`. This is intentional — silent magnitude mismatches were the v0.6 incident class that motivated R-15 / R-16. Fix the target card to match the chemistry-true charge from the snapshot PDB. See `target_cards/2QKI.json` and `utils/charge_topology.py::compute_binder_chem_charge` for the canonical examples.

---

## 7. Verifying you match the published numbers

The Paper 1 v1 / v2 Table 1 numbers (`outputs/paper1/paper1_v2_full_manuscript.md` §3.3) were generated under the v0.7.1 production state described in this document. To verify your installation matches:

```bash
# Cp4 family-level baseline:  ⟨⟨ΔG⟩⟩ = -1.823, σ_btwn = 6.656, CI = 3.65, Tier-3
python scripts/aggregate_mmpbsa_family.py --family 2QKI_Cp4 --subdir mmpbsa_results_postl387_v2

# 1EBP_MTR13 family-level baseline:  ⟨⟨ΔG⟩⟩ = -20.442, σ_btwn = 3.188, CI = 0.156, Tier-1
python scripts/aggregate_mmpbsa_family.py --family 1EBP_MTR13 --subdir mmpbsa_results_postl387
```

If your numbers differ by more than ±0.1 kcal/mol on the family-level ⟨⟨ΔG⟩⟩, check:

1. Are you on `commit fe4309b` or later? (`git log -1 --oneline`)
2. Is `UPDD_MTR_AMBER14_PATCH` unset? (must default to ON per `utils/parameterize_ncaa.py:1493`)
3. Did you regenerate snapshots, or are you using the Zenodo cache? (regenerating MD requires identical RNG seeds and OpenMM versions to be reproducible at floating-point precision)
4. Are AmberTools / OpenMM / mdtraj versions pinned per `environment.yml`?

The v2 Appendix S-X patch-on / patch-off sensitivity envelope (`paper1_v2_full_manuscript.md` Appendix S-X) is the **methodology-disclosure artifact**, not the production baseline. The frozen production numbers are the patch-on values reported in Table 1.

---

## 8. Citation

If you use UPDD or this guide in your research, please cite the Zenodo software DOI and the Paper 1 preprint per `README.md` § "Citation".

---

## 9. Support

This is an independent-researcher project; there is no commercial support. Issues and discussion: https://github.com/kains675/UPDD/issues
