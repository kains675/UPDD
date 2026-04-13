# UPDD — Universal Peptide Drug Discovery

> AI-driven computational pipeline for designing cyclic peptide binders with non-canonical amino acid (ncAA) integration.

## Overview

UPDD is an end-to-end *in silico* peptide drug design system that bridges generative protein design with quantum-level binding validation. The pipeline generates, evaluates, and ranks cyclic peptide candidates targeting user-specified protein interfaces — from backbone generation to binding free energy estimation.

## Pipeline

```
RFdiffusion → ProteinMPNN → AlphaFold2 → ncAA Substitution → MD → QM/MM → MM-GBSA → Ranking
```

| Stage | Method | Purpose |
|-------|--------|---------|
| Backbone Generation | RFdiffusion (RFpeptides) | *De novo* cyclic peptide backbone design |
| Sequence Design | ProteinMPNN | Sequence optimization for designed backbones |
| Structure Prediction | ColabFold (AF2 Multimer) | Complex structure prediction + ipTM/pLDDT scoring |
| ncAA Integration | Custom mutation engine | Non-canonical amino acid substitution + GAFF2 parametrization |
| Molecular Dynamics | OpenMM (AMBER ff14SB) | Restrained MD with 2-pass NaN recovery |
| Electronic Structure | gpu4pyscf (ωB97X-D/6-31G*) | QM/MM single-point energy with adaptive GPU parallelism |
| Binding Energy | OpenMM (GBn2) | MM-GBSA ΔG estimation |

## Key Features

- **ncAA Support**: 13+ non-canonical amino acids (N-methylalanine, silaproline, halogenated phenylalanines, D-amino acids, etc.)
- **Cyclic Peptide Topologies**: Head-to-tail, disulfide, thioether cyclization with AF2 gap closure
- **Adaptive GPU Scheduling**: 4→3→2→1 automatic worker reduction on VRAM OOM
- **Robust MD**: 2-pass timestep strategy (2fs→1fs) with checkpoint recovery
- **PBC-aware Snapshots**: Minimum-image convention unwrapping for multi-chain complexes
- **Resume Support**: Interrupted QM/MM and MM-GBSA calculations resume from last completed snapshot

## Hardware Requirements

Developed and tested on:
- **GPU**: NVIDIA RTX 5070 Ti (16GB VRAM)
- **CPU**: AMD Ryzen 9800X3D (8C/16T)
- **RAM**: 32GB

## License

This project is proprietary. All rights reserved.

## Author

PotionMaker — Computational Peptide Chemistry & AI-driven Drug Design
