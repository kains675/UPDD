# ATM-Nested Dynamic Ghost P0 Preregistration Amendment R2

Registered: 2026-07-17 KST, after context-free KEEPER development smoke and
before official P0 inventory, KEEPER output, or RUNNER output

Status: `FROZEN_AFTER_KEEPER_SMOKE_COORDINATE_PRECISION_DISCOVERY_BEFORE_OFFICIAL_OUTPUT`

SciVal verdict: `CONDITIONAL APPROVE`

## Preserved R0 And R1

R0 and R1 remain unchanged. R2 changes no Hamiltonian, particle selection,
probe metric, scientific threshold, or runtime.

## KEEPER-Smoke Finding

The context-free `w4a_union_s101_free` KEEPER smoke stopped before force
construction because the frozen normalized coordinates give an ATM anchor
offset norm of `4.000026749911 nm`, outside the R0 declaration tolerance of
`1e-8 nm`.

A read-only six-cell audit found offset norms from `3.999977757438` to
`4.000026749911 nm`. All nine selected ring atoms still share the exact same
declared `ParticleOffsetDisplacement` indices in each cell, and all have no
`u0` displacement. The deviation is consistent with the coordinate precision
of the normalized PDB-derived position arrays, not transformation drift.

## Frozen R2 Correction

Keep the expected offset norm at `4.0 nm` and change only its declaration
tolerance from `1e-8` to `1e-4 nm`.

The `1e-4 nm` tolerance is tied to source coordinate serialization precision,
not selected from the maximum observed deviation. Exact transformation type,
common destination/origin indices, absence of `u0` displacement, and zero
water displacements remain hard gates.

This amendment makes no free-energy or production claim.

