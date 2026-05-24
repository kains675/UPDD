---
date: <% tp.date.now("YYYY-MM-DD") %>
type: experiment
cycle: 0              # 0 = baseline, 1+ = AL cycle
system: ""            # e.g. 2QKI_Cp4, 1EBP_MTR13
ncaa: ""              # e.g. MTR, Cp4
status: planned       # planned | running | done | failed | invalidated
tier: ""              # 1 | 2 | 3 per UPDD Tier classification
ddg_kcal_mol: null
sigma_btwn: null
sigma_w: null
ci: null              # Convergence Index
z_se: null
n_replicates: null
seed: null
schema: ""            # e.g. branched_ddg/0.3
workspace: ""         # absolute path on /media/san/...
tags:
  - experiment
---

# <% tp.date.now("YYYY-MM-DD") %> — <system>_<ncaa> cycle <cycle>

## Hypothesis
What we expect to learn or rank.

## Protocol
- Pipeline stage(s): RFdiffusion / ProteinMPNN / AF2 / MD / QM/MM / MM-PBSA
- Parameters: dt, force field, n_qm, λ-schedule, etc.
- Patches active: e.g. L387-v54 + CONECT-v55 + AMBER14 q_N=-0.4157

## Results
Fill in once `status: done`. Reference the aggregator output:
- `mmpbsa_summary.json`: `<path>`
- σ-decomposition: σ_btwn = , σ_w = , CI = 

## Interpretation
Sign vs WT? Tier-preserving? Inside Magotti SSOT envelope?

## Anomalies / detachment
Any frames flagged by `intra_residue_bond_check` or `detachment_metric > 0.5`?

## Links
- Decision: [[ADR-XXXX]]
- System note: [[<system>]]
- ncAA note: [[<ncaa>]]
