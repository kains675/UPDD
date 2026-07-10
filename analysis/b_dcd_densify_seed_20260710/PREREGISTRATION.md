# B-DCD Densify + Seed Expansion Preregistration

Date: 2026-07-10

## Decision

Keep the B-first path, but attach a DCD structural gate. Promote A carve
redesign only if, after densify and seed expansion, either:

1. the bound/dminus control baseline persists, or
2. the high-leverage cells remain amplified above matched controls.

Forward-only C remains a diagnostic cross-check only.

## Scientific Regime

- Track B RBFE values remain ranking/SIGN-only.
- `dgb = dgbind1 - dgbind2` is the primary axis.
- `dgbind1` remains diagnostic only.
- DCD structural readouts are pathology gates, not calibrated structural
  distances or FE estimators.

## D2 Basis

D2 confirmed a local water/site close-contact signal, but not a seed-unique one.

Key D2 fractions, using `frac<0.26 nm`:

- W23A `s101 bound/dminus`: `0.293`
- W23A control `s42 bound/dminus`: `0.180`
- W4A carved `s127 bound/dminus`: `0.515`
- W4A carved control `s163 bound/dminus`: `0.275`
- W4A carved `s127 free/dplus`: `0.167`
- W4A carved control `s101 free/dplus`: `0.175`

Interpretation: bound/dminus has a broader structural baseline that high-
leverage seeds can amplify. Densify/seed expansion must therefore carry a
structural gate, not just a scalar `.out` gate.

## Campaigns

### Campaign 1: W4A Carved/Uncarved H18 Densified Control

Purpose: resolve whether C1 remains INDETERMINATE after a denser, same-protocol
paired cohort.

Scope:

- system: 2QKI W4A Trp4->Ala, WT endpoint
- compare: carved vs uncarved
- legs: free and bound
- directions: dplus and dminus
- seeds: `s7,s19,s23,s42,s83,s101,s127,s163,s199,s251`
- matched cohort size: `n=10`

Protocol:

- two-copy ATS
- fixed displacement: `4.0 nm`
- mutation: `w4a_trp_ala_res4`
- densified H18 ladder:
  - `lambda1_rampdown = 0.025,0.05,0.1,0.2,0.3,0.4,0.5`
  - `lambda2_rampup = 0,0.05,0.1,0.15,0.2,0.3,0.4,0.5`
- `n_cycles = 800`
- `md_steps_per_cycle = 250`
- `timestep_fs = 1.0`
- `staged_min = true`
- `reseed_perm_seed = 20260618`
- `reseed_endpoint = true`
- `reseed_endpoint_band_lambda2_max = 0.25`
- `reseed_endpoint_equil_steps = 2000`
- `mintimeid = 600`

Carve arm:

- `--carve-void-waters`
- `--carve-cutoff-nm 0.26`

Uncarved arm:

- no carve flags

Rationale:

- Rerun carved and uncarved under the same H18 dense schedule instead of mixing
  old 400-cycle C1 values with new values.
- Use `n=10` to reduce single-seed leverage.
- Keep W4A fixed-displacement `4.0 nm` for continuity with the original C1 box
  family, while explicitly declaring it to satisfy the launcher displacement
  gate.

### Campaign 2: W23A Seed Expansion

Purpose: keep W23A as the in-place construction large-effect anchor while
reducing leverage of `s101`.

Current state:

- available 1YCR scaffolds: `s7,s19,s23,s42,s101`
- current cohort size: `n=5`
- additional 1YCR scaffold seeds are not currently present in
  `outputs/1YCR_WT_calib_s*/mdresult/1YCR_WT_final.pdb`.

Protocol when additional 1YCR scaffolds are available:

- use the existing W23A orchestrator protocol unchanged:
  - H18 densified ladder
  - `n_cycles = 800`
  - `md_steps_per_cycle = 250`
  - `mintimeid = 600`
  - `--staged-min`
  - `--carve-void-waters --carve-cutoff-nm 0.26`
  - auto-search displacement, `accept_sep_nm = 1.5`
- append new seeds in a fresh out-root or rerun a full same-protocol cohort;
  never mix rep-indexed directories with different seed order.

Stop condition:

- do not claim W23A seed expansion until additional 1YCR scaffold final PDBs
  exist and pass the scaffold availability gate.

## DCD Structural Gate

The production launcher's `--dcd` excludes solvent, so it is insufficient for
the water-site gate. The structural gate is therefore a sidecar local-water DCD
probe, based on the D2 method:

- run from serialized production boxes after the FE run;
- copy source XML/PDB into an isolated DCD-gate output root;
- record all non-solvent atoms plus initial waters within `0.60 nm` of the
  alchemical-site heavy atoms;
- use one DCD frame per cycle for the short structural probe;
- compare high-leverage bound/dminus cells to matched controls.

Primary structural metrics:

- `frac_lt_0p26`: fraction of frames with a selected water oxygen within
  `0.26 nm` of a site heavy atom;
- `frac_lt_0p35`: contact fraction within `0.35 nm`;
- site heavy-atom RMSD p95;
- DCD frame count vs `.out` row count.

Post-densify A-promotion gate:

- bound/dminus control baseline persists if matched control median
  `frac_lt_0p26 >= 0.15`;
- high-leverage amplification persists if
  `high_leverage - control_median >= 0.10` and
  `high_leverage / control_median >= 1.5`;
- either condition after densify+seed expansion promotes A carve redesign.

## Statistical Gate

For FE results:

- paired same-seed delta is primary;
- robust drop-max-|z| seed recheck is mandatory;
- bootstrap CI is preferred where implemented;
- point-threshold-only C1 classification is rejected;
- C1 remains INDETERMINATE unless paired/robust evidence resolves it.

## Stop Gates

- missing seed scaffold or final PDB;
- out-root seed collision or rep-index/seed mismatch;
- mixed carved/uncarved protocol between arms;
- UWHAM schedule mismatch or fallback to hardcoded state count;
- OpenMM NaN;
- `.out` row count mismatch;
- DCD frame count mismatch;
- adjacent overlap collapse below the existing Track B hard floor;
- DCD post-densify structural gate triggers A-promotion criteria;
- sign flip or new heterogeneity against the current W23A/W4A anchors.

## Immediate Runner Plan

1. Run dry-run pre-registration for W4A H18 carved and uncarved free/bound arms.
2. Keep W23A seed expansion blocked until additional 1YCR scaffolds exist.
3. Use `bdcd_gate.py` to evaluate the current D2 summary as the precheck
   baseline and later the post-densify sidecar summaries.
4. Launch GPU production only after dry-run artifacts pass the declaration-to-
   artifact check.
