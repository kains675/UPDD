# P0v3r1 PATH Diagnosis

Date: 2026-07-17 KST

## Formal outcome

`REJECT_UNION_SOURCE` at `w4a_union_s101_bound` (`1/6` source cells built;
the remaining five were not launched). No OpenMM Context, MD, minimization,
GPU, or free-energy estimator was used.

Frozen inventory digest:
`ed472c770ad439fcff0b748c91efc65fd0f82076cbebe7552ff51767d8b18d48`.

## Passed source evidence

The union construction itself passed every preregistered physical and artifact
gate:

- particles `300424`, water/Na/Cl `96534/263/263`, total solvent `97060`;
- source and parent net charge both `-1.7053025658242404e-13 e`;
- nine placeholders, exact source parameters, target max error
  `8.881784197001252e-16 nm`, final placeholder count zero;
- stored/transformed ring-water minima `0.38331698/0.39665273 nm`, with zero
  pairs below `0.26 nm`;
- cubic right-handed box, canonical solute ordering and two-copy separation
  PASS;
- the frozen solvent RNG was consumed, and the global RNG state was restored
  across the `addSolvent` call.

The static total-box solvent-density delta was `-1.3923%`; this is report-only
under P0v3r1 and is not an equilibrium-density or free-energy gate.

## Failed check and pathology

The sole false check was `caller_rng_state_unchanged`. It compared the global
Python RNG before and after the entire canonical two-copy build. The frozen
protocol requires save/restore only around the cell-seeded `addSolvent` call,
which passed. Therefore this whole-builder check exceeded the preregistered
scope and cannot by itself reject union solvation.

The changed RNG state nevertheless exposed a real reproducibility gap. Relative
to the frozen uncarved parent, the rebuilt source retained identical solute atom
identity/order but changed nine hydrogen coordinates. The largest differences
were ALA4 `HB1/HB2/HB3` at `0.2155/0.2178/0.2137 nm`; six additional residue-4
hydrogens differed by `0.00024-0.00750 nm`. This is consistent with the existing
unseeded PDBFixer/OpenMM appearing-hydrogen placement path.

PATH classification:
`AUDIT_SCOPE_ERROR_WITH_REAL_APPEARING_H_REPRODUCIBILITY_SIGNAL`.

## Required correction

Do not reclassify or overwrite P0v3r1. A new preregistered revision must:

1. remove only the non-preregistered whole-builder RNG-invariance check;
2. retain the hard checks that the frozen solvent RNG is consumed and restored
   specifically around `addSolvent`;
3. use the existing deterministic per-unit appearing-H placement and bounded
   R2 retry with `K=5`, recording the full derived seed trail;
4. rebuild all six sources fresh rather than reuse the rejected cell;
5. retain every count, charge, placeholder, geometry, bridge, CPU-only, and
   stop-on-first-failure gate unchanged;
6. continue to prohibit MD, GPU, free-energy estimation, and automatic next-stage
   launch.

No threshold or scientific output from the rejected source is used to select an
appearing-H seed.

