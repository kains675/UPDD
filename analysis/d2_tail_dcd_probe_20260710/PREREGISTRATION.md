# D2 Tail DCD Probe Preregistration

Date: 2026-07-10

## Purpose

D2 follows the read-only tail-mechanism probe. Its goal is structural pathology
diagnosis, not a new free-energy estimate. It records short per-walker DCDs for
the high-leverage W23A and W4A carved cells, including a local water shell around
the alchemical site.

## Output Isolation

- Source production artifacts are read-only.
- Serialized source boxes are copied into
  `outputs/_trackb/d2_tail_dcd_probe_20260710/`.
- D2 never writes inside the source production roots.
- D2 does not run UWHAM and does not update any production `dgb` result.

## Run Length

- `n_cycles = 50`
- `md_steps_per_cycle = 250`
- `timestep_fs = 1.0`
- physical time per walker = 12.5 ps
- `dcd_stride_cycles = 1`

Rationale: this is a short structural visibility probe. It preserves the
production per-cycle MD step size and samples enough cycles to expose immediate
site-water or basin-lock behavior without attempting convergence.

## DCD Atom Selection

For each cell:

- all non-solvent atoms;
- all atoms in water residues initially within `0.60 nm` of the alchemical-site
  heavy atoms;
- site heavy atoms are chain `B`, residue id `7` for W23A and chain `B`,
  residue id `4` for W4A, across both two-copy residues.

Water diagnostics use two preregistered thresholds:

- intrusion: any selected water oxygen within `0.26 nm` of a site heavy atom;
- contact: any selected water oxygen within `0.35 nm` of a site heavy atom.

## Cells

W23A:

- `w23a_s101_bound_dplus`
- `w23a_s101_bound_dminus`
- `w23a_s101_free_dplus`
- `w23a_ctrl_s42_bound_dminus`

W4A C1 carved:

- `w4a_carved_s127_free_dplus`
- `w4a_carved_s127_bound_dminus`
- `w4a_carved_ctrl_s101_free_dplus`
- `w4a_carved_ctrl_s163_bound_dminus`

## Stop Gates

- source serialized XML/PDB missing;
- source seed/replicate mismatch;
- D2 output root collision with a different source box;
- selected alchemical-site heavy atoms absent;
- selected local-water shell absent when source PDB contains waters near the site;
- OpenMM NaN fail-fast;
- GPU memory/runtime failure;
- DCD frame count not matching `.out` row count at stride 1.

## Decision Use

- If outlier cells show water intrusion or persistent local-water site occupancy
  that controls lack, promote carve redesign `A` for that mechanism.
- If outlier cells show stable no-intrusion trajectories but poor endpoint/apex
  exchange, prefer targeted `B-local` densify plus seed expansion.
- If both outliers and controls show the same structural signature, treat it as
  cohort-wide setup/sampling pathology, not seed-specific proof.

## Non-Goals

- No calibrated ddG claim.
- No carve-neutrality proof.
- No production-engine change.
