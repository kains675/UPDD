# UPDD — Universal Peptide Drug Discovery

> AI-driven computational pipeline for designing cyclic peptide binders 
> with non-canonical amino acid (ncAA) integration.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Release](https://img.shields.io/badge/release-v0.7.1-blue.svg)](https://github.com/kains675/UPDD/releases/tag/v0.7.1-paper1-v1)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20067323.svg)](https://doi.org/10.5281/zenodo.20067323)

## Overview

UPDD is an end-to-end *in silico* peptide drug design system that bridges generative protein design with quantum-level binding validation. The pipeline generates, evaluates, and ranks cyclic peptide candidates targeting user-specified protein interfaces — from backbone generation to binding free energy estimation.

This release (v0.7.1) accompanies **Paper 1 v1** — a Capability Level 1 methodology release documenting the σ_btwn / σ_w ensemble-quality decomposition, the Convergence Index dual interpretation, the four-layer onion-peeling verification protocol, and the simulated active-learning cycle on the six-target ncAA benchmark. See the [Citation](#citation) section below for the ChemRxiv preprint reference.

## Pipeline
RFdiffusion → ProteinMPNN → AlphaFold2 → ncAA Substitution → MD → QM/MM → MM-PBSA → Ranking

| Stage | Method | Purpose |
|-------|--------|---------|
| Backbone Generation | RFdiffusion (RFpeptides) | *De novo* cyclic peptide backbone design |
| Sequence Design | ProteinMPNN | Sequence optimization for designed backbones |
| Structure Prediction | ColabFold (AF2 Multimer) | Complex structure prediction + ipTM/pLDDT scoring |
| ncAA Integration | Custom mutation engine | Non-canonical amino acid substitution + GAFF2 parametrization |
| Molecular Dynamics | OpenMM (AMBER ff14SB) | Restrained MD with **3-pass NaN recovery** |
| Electronic Structure | gpu4pyscf (ωB97X-D/6-31G*) | QM/MM single-point energy with adaptive GPU parallelism |
| Binding Energy | OpenMM (PBSA / GBn2) | MM-PBSA / MM-GBSA ΔG estimation |

## Why UPDD?

Current AI-driven peptide design pipelines typically stop at structure prediction, relying on experimental validation to assess binding. Furthermore, state-of-the-art tools like RFdiffusion and ProteinMPNN were trained exclusively on canonical amino acids and **cannot model ncAAs** [[1]](#references). UPDD addresses both gaps:

| Capability | RFdiffusion | NCFlow | CyclicChamp | **UPDD** |
|------------|:---:|:---:|:---:|:---:|
| *De novo* backbone generation | ✅ | — | ✅ | ✅ |
| Protein binder design | ✅ | — | — | ✅ |
| ncAA integration | — | ✅ | ✅ | ✅ |
| Force field parametrization | — | — | Rosetta | GAFF2 |
| Molecular dynamics | — | — | validation only | ✅ |
| QM/MM (DFT) | — | — | — | **✅** |
| MM-PBSA / MM-GBSA ΔG | — | — | — | **✅** |
| End-to-end automation | — | — | — | **✅** |

**UPDD is, to our knowledge, the first pipeline that integrates AI-driven generative design with ncAA substitution, MD simulation, and quantum-mechanical binding validation in a single automated workflow.**

## Capability Level

UPDD's scientific-capability progression is stratified into three levels:

- **Level 1 — Decision support** (current scope, v0.7.1): architectural integrity verified, sign-significance demonstrated, σ_btwn / σ_w decomposition + Convergence Index framework, layered verification protocol with empirical bug-catching evidence. ΔΔG outputs serve as decision support for ranking ncAA variants in iterative design cycles, *not* as quantitative experimental predictions.
- **Level 2 — Cross-system reproducibility** (v0.8 target): sign-significance reproduces on three+ systems spanning multiple ncAA classes; bit-identical WT control as routine reproducibility benchmark.
- **Level 3 — Quantitative match** (v1.0+ target): ΔΔG predictions within field-comparable distance of experimental anchors; Khoury R²=0.388 ceiling acknowledged.

The v0.7.1 release described in Paper 1 v1 is explicitly a Level 1 release.

## Key Features

- **ncAA Support**: 25 registered ncAAs at v0.7 (MTR-class verified, N-methyl class deferred to N-methyl overlay)
- **Cyclic Peptide Topologies**: Head-to-tail, disulfide, thioether cyclization with AF2 gap closure
- **Adaptive GPU Scheduling**: 4→3→2→1 automatic worker reduction on VRAM OOM
- **Robust MD**: 3-pass timestep strategy (Pass 1 dt=2fs → Pass 1b reseed at seed+13 → Pass 2 dt=1fs cyclic-only) with checkpoint recovery
- **PBC-aware Snapshots**: Four-layer L387 v54 PBC repair + CONECT v55 atom-index disambiguation patch + intra-residue bond integrity guard for HETATM ncAAs
- **Branched ΔΔG**: Bit-identical wild-type control with intervention isolation verified at floating-point precision
- **σ Decomposition**: Orthogonal-axis ensemble-noise statistics (between-replicate σ_btwn vs within-replicate σ_w) with the Convergence Index dual interpretation
- **Resume Support**: Interrupted QM/MM and MM-GBSA calculations resume from last completed snapshot

## Hardware Requirements

UPDD is designed to run on a **single consumer-grade workstation** — no HPC cluster or cloud GPU required.

| Component | Spec | Note |
|-----------|------|------|
| OS | Ubuntu 24.04.4 LTS | |
| GPU | NVIDIA RTX 5070 Ti (16GB) | Single consumer GPU — adaptive parallelism handles VRAM constraints |
| CPU | AMD Ryzen 9800X3D (8C/16T) | MM-PBSA runs on CPU in parallel while GPU handles DFT |
| RAM | 32GB | Sufficient for all pipeline stages |

Most existing pipelines (RFdiffusion [[2]](#references), ProteinDJ) are benchmarked on A100/A30 HPC clusters. UPDD achieves the same workflow — including QM/MM DFT — on hardware accessible to independent researchers and small labs.

## Reproducibility

See **[`INSTALL.md`](INSTALL.md)** for the step-by-step setup guide covering three reproducibility tiers:

- **Tier A** — CPU-only conda env for code review / pytest / schema regression (~15 min, ~3 GB)
- **Tier B** — Full pipeline with GPU for MD + QM/MM + MM-PBSA reproduction (~50 GB, 1× NVIDIA GPU ≥ 16 GB VRAM)
- **Tier C** — Docker container for CI / external reviewers (CPU-only)

Reproducibility primitives shipped in the repo:

- `environment.yml` — conda environment specification (`qmmm` env, Tier-1 packages pinned)
- `Dockerfile` — CPU-only validation container (Tier C)
- `.github/workflows/ci.yml` — continuous integration (lint, unit tests, Docker smoke test)
- `tests/` — regression suite (294 pass / 9 skip / 1 pre-existing fail)

The release is archived at Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20067323.svg)](https://doi.org/10.5281/zenodo.20067323)

## Citation

If you use UPDD in your research, please cite:

**Paper:** Kang, I. (2026). *An Iterative Evaluation Framework for Non-Canonical Amino Acid Peptide Drug Discovery: Sign-Direction Reproducible Multi-Layer Verification*. ChemRxiv. DOI: [10.26434/chemrxiv.15002948/v2](https://doi.org/10.26434/chemrxiv.15002948/v2).

**Software:** Kang, I. (2026). UPDD: Universal Peptide Drug Discovery — v0.7.1 [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.20067323. ORCID: [0009-0007-0753-0636](https://orcid.org/0009-0007-0753-0636).

## References

1. Lee & Kim, "Design of peptides with non-canonical amino acids using flow matching," *bioRxiv* (2025). — Documents that RFdiffusion / BindCraft cannot model ncAAs.
2. Watson et al., "De novo design of protein structure and function with RFdiffusion," *Nature* **620**, 1089–1100 (2023).
3. Dauparas et al., "Robust deep learning–based protein sequence design using ProteinMPNN," *Science* **378**, 49–56 (2022).
4. Zhu et al., "Heuristic energy-based cyclic peptide design," *PLoS Comput Biol* **21**(4), e1012290 (2025). — Rosetta-based ncAA cyclic peptide design (no binder design, no DFT).
5. Rettie et al., "Cyclic peptide design with RFpeptides," *Nat Chem Biol* (2025).

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Author

PotionMaker ([Insan Kang](https://orcid.org/0009-0007-0753-0636)) — Independent Researcher  
Bio-Organic Chemistry & Computational Peptide Chemistry & AI-driven Drug Design

For correspondence, please use the email address listed on the [ORCID profile](https://orcid.org/0009-0007-0753-0636) or the ChemRxiv preprint linked in the Citation section.
