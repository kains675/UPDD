# ATM-Nested Transformed-Site Dynamic Ghost P0 Preregistration

Registered: 2026-07-17 KST, before implementation or P0 output

Status: `FROZEN_BEFORE_IMPLEMENTATION_OR_OUTPUT`

SciVal verdict: `CONDITIONAL APPROVE`

## Trigger

The P1SD sampled bridge failed at `w4a_union_s101_free`, `xi=0.00`. P1CAP
showed that all 50 failed frames were dplus cap-active, with raw endpoint gaps
far beyond the ATM soft-core boundary. The previous dynamic ghost was a
top-level force evaluated only at stored particle coordinates. It therefore
could not protect the transformed W4A ring image used by the nested ATM
Hamiltonian.

## Frozen Inputs

Use the exact normalized P1SD artifacts for:

- `w4a_union_s101_bound`
- `w4a_union_s127_bound`
- `w4a_union_s163_bound`
- `w4a_union_s101_free`
- `w4a_union_s127_free`
- `w4a_union_s163_free`

For each cell, consume only:

- `normalized_source_system.xml`
- `normalized_source.pdb`
- `positions_float64.npy`
- `normalization_result.json`

The source System must have one top-level `ATMForce`, no dynamic ghost, and the
same particle count and parent box declared by the frozen P1SD inventory.

## Frozen Hamiltonian Change

Append one `CustomNonbondedForce` to the existing `ATMForce`, not to the
top-level `System`. The force:

- uses the exact nine heavy atoms of chain B TRP4;
- interacts only with all explicit-water oxygen atoms;
- obtains sigma and epsilon from the canonical nested `NonbondedForce`;
- uses Lorentz-Berthelot mixing and the previously frozen soft-core WCA
  expression;
- has global parameters `g` and `alpha_sc=0.5`;
- has no charge term, attraction, added particle, deleted water, switching
  function, or long-range correction;
- uses `CutoffPeriodic` with a `0.40 nm` neighbor cutoff;
- declares the `g` energy-parameter derivative.

All selected ring atoms must use one common `ParticleOffsetDisplacement` for
`u1`, with no `u0` displacement. Every selected water oxygen and every atom in
the selected probe water must have zero `FixedDisplacement` in both endpoints.

## Context-Free KEEPER Gates

Before any Context exists, all six cells must pass:

1. exactly one ATM force and one canonical nested nonbonded force;
2. exactly nine declared ring heavy atoms and one oxygen per explicit water;
3. no pre-existing top-level or nested dynamic ghost;
4. top-level force count unchanged and ATM nested-force count increased by one;
5. particle count, masses, constraints, virtual sites, and net charge
   unchanged;
6. all pre-existing top-level and ATM-nested force serializations unchanged;
7. ATM energy expression, global parameters, and all particle transformations
   unchanged;
8. nested ghost expression, parameters, interaction group, cutoff, and
   derivative match the frozen protocol;
9. XML serialization round-trip reproduces the same declaration report.

Any failure is `KEEPER_BLOCK`.

## Deterministic Probe Construction

The target ring atom is `NE1`. Define the outward vector from the stored ring
heavy-atom centroid to stored `NE1`. The `u1` probe uses the same vector because
the nine ring atoms share one rigid `4.0 nm` particle-offset displacement.

Set the probe oxygen at `0.85` times the NE1-water WCA cutoff from target `NE1`.
Choose one three-atom water whose normalized H-centroid vector has the largest
dot product with the outward vector, breaking ties by the lowest oxygen index.
Translate the whole water rigidly; do not rotate it or alter its internal
coordinates.

Build two diagnostic coordinate sets:

- `stored_probe`: oxygen placed outward from stored `NE1`;
- `transformed_probe`: oxygen placed outward from transformed `NE1`.

These coordinates are diagnostics only and must not enter sampling.

## Reference RUNNER Gates

Run one isolated process and one OpenMM Context at a time per System variant.
Use OpenMM `Reference`; CUDA is forbidden. MD steps and minimization steps are
both zero.

For the unmodified source and serialized nested-ghost System:

- compare `u0_dplus` at `g=0` to unmodified `u0_dplus`;
- compare `u1_dminus` at `g=0` to unmodified `u1_dminus`;
- require absolute energy delta at most `1e-5 kJ/mol`;
- require maximum force-component delta at most
  `1e-5 kJ/mol/nm`.

At source coordinates with `g=1`, require finite total/raw endpoint energies,
forces, and `dV/dg`. Reject gross behavior above `1e3 kJ/mol` isolated raw
ghost energy or `1e5 kJ/mol/nm` isolated force-component magnitude.

For each probe, isolate the ghost by `g=1 - g=0` at both endpoint states:

- intended endpoint isolated energy must exceed `0.1 kJ/mol`;
- opposite endpoint absolute isolated energy must be at most
  `1e-6 kJ/mol`;
- intended endpoint oxygen-force projection along the outward vector must
  exceed `1.0 kJ/mol/nm`;
- opposite endpoint maximum isolated force component must be at most
  `1e-5 kJ/mol/nm`;
- the vector sum of isolated forces must be at most
  `1e-5 kJ/mol/nm` per Cartesian component.

The raw endpoint energy difference returned by `ATMForce` is recorded at every
state. The `g` derivative must exist and be finite, but no equality between a
derivative and a finite `g` difference is assumed because the WCA coupling is
nonlinear.

## Decisions

1. `P0_PASS`: every KEEPER and RUNNER gate passes in all six cells. Proceed
   only to a separately preregistered state-specific DCD structural replay.
2. `REJECT_GHOST`: any finite/parity/serialization/probe gate fails. Do not
   sample; diagnose force placement or ATM transformation semantics.
3. `INDETERMINATE`: execution or source integrity prevents a complete result.
   Preserve partial output and rerun only after a documented repair.

Even `P0_PASS` does not authorize production and does not establish
free-energy neutrality.

