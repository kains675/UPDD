# W4A Post-densify Bound/Dminus DCD Structural Gate

Status: `FROZEN BEFORE SIDECAR OUTPUT`

Date: 2026-07-16

Regime: Track B ranking/SIGN-only; structural mechanism diagnostic, not a free-
energy estimate.

## Question

After the H18 densify plus seed-expansion campaign, does the W4A carved
`bound/dminus` water-intrusion baseline or the preregistered `s127`
high-leverage amplification persist strongly enough to promote an A-branch
carve redesign?

The completed FE campaign is not converged at the hardened post-warmup
transport gate. Its free-energy values must not select cells or tune this
probe. The apparent local `4-5` and `9-10` dplus/dminus labels refer to the same
physical lambda2 `0.20 -> 0.30` interval; they are not two independent physical
bottlenecks.

## Frozen Source Cohort

Source root:

`outputs/_trackb/w4a_bdcd_h18_20260710_rerun1`

Only the following two completed arms are inputs:

- `carved_bound/wt/bound/rep0..9`
- `uncarved_bound/wt/bound/rep0..9`

The positional seed order is frozen as:

`s7,s19,s23,s42,s83,s101,s127,s163,s199,s251`

For every seed, only endpoint `wt`, leg `bound`, and direction `dminus` are
probed. This yields 20 cells: 10 carved and 10 uncarved. No cell is selected or
removed using the completed campaign's estimated free energy.

The source preregistration, `run_manifest.json`, serialized XML/PDB, and source
`dminus` control file are hashed into `source_inventory.json` before launch.
Each source direction must contain 15 walkers x 800 finite rows with 11 numeric
columns, states `0..14`, and no recorded NaN.

## Frozen Sidecar Protocol

The sidecar reuses the tested D2 implementation without modifying the Track B
production engine:

- isolated output root:
  `outputs/_trackb/w4a_postdensify_dcd_20260716`
- one short ladder per cell, serial execution on one GPU;
- H18 schedule copied from the source declaration:
  - lambda1 ramp-down: `0.025,0.05,0.10,0.20,0.30,0.40,0.50`
  - lambda2 ramp-up: `0.00,0.05,0.10,0.15,0.20,0.30,0.40,0.50`
- `50` cycles x `250` MD steps/cycle, `1 fs` timestep;
- staged minimization, 500 minimization iterations, 500 backward-equilibration
  steps;
- RNG seed `20260613 + 1000 * (replicate_index + 1)`, paired between carved
  and uncarved arms;
- one DCD frame per cycle;
- DCD selection: all non-solvent atoms plus every atom of each water residue
  initially within `0.60 nm` of W4A site heavy atoms (`chain B`, residue `4`);
- no UWHAM and no free-energy estimate from sidecar output.

The serialized source box is an initial production construction, not a
checkpoint from the end of the 800-cycle FE trajectory. This probe therefore
tests short-run structural propensity under a fixed paired protocol; it does
not claim to replay the FE tail event.

## Primary Metrics

For every cell:

- `frac_lt_0p26`: fraction of DCD frames whose minimum selected-water oxygen to
  site-heavy distance is below `0.26 nm`;
- `frac_lt_0p35`: analogous contact fraction below `0.35 nm`;
- site-heavy RMSD p95;
- DCD frame count and `.out` row count.

The required DCD count is 15 walkers x 50 cycles = 750 frames per cell. Frames
from walkers and adjacent cycles are correlated; these fractions are
descriptive mechanism metrics, not independent Bernoulli observations.

## Frozen A/B Gate

The high-leverage value `H` is the carved `s127 bound/dminus`
`frac_lt_0p26`, fixed from the earlier D2 diagnosis.

The matched control baseline `C` is the median `frac_lt_0p26` over the other
nine carved seeds. The full non-high carved cohort is used so that the control
cannot be selected post hoc.

Two preregistered A-promotion conditions are evaluated:

1. control baseline persists: `C >= 0.15`;
2. high-leverage amplification persists:
   `H - C >= 0.10` and `H / C >= 1.5`.

If `C == 0`, the ratio is undefined and condition 2 does not trigger. If either
condition triggers, verdict is `PROMOTE_A_CARVE_REDESIGN`. Otherwise verdict is
`CONTINUE_B_TRANSPORT_PILOT`.

The uncarved 10-seed median and paired carved-minus-uncarved seed differences
are secondary context only. They may distinguish construction-wide from
carve-associated water propensity, but cannot independently trigger or cancel
the frozen A-promotion rule.

## Stop Gates

Stop without a scientific verdict on any of the following:

- source seed, replicate, arm carve declaration, schedule, or hash mismatch;
- missing serialized source artifact;
- source `.out` count, state range, column count, or finite-value failure;
- staged destination collision with different content;
- OpenMM NaN;
- missing local site-heavy selection or selected water oxygens;
- DCD count other than 750 per completed cell;
- DCD/.out frame mismatch;
- incomplete 20-cell cohort at gate evaluation.

Interrupted partial sidecar cells are preserved under the output root's
`_archive/` before a clean cell restart. Completed, integrity-valid cells are
resume-skipped.

## Branch After This Gate

- `PROMOTE_A_CARVE_REDESIGN`: stop B transport tuning and preregister an
  explicitly dG-faithful carve redesign comparison.
- `CONTINUE_B_TRANSPORT_PILOT`: run overlap-only schedule/cadence pilots. Pilot
  selection must not use free-energy agreement, and full matched n=10 reruns
  remain blocked until post-burn-in round-trip and adjacent-overlap gates pass.
