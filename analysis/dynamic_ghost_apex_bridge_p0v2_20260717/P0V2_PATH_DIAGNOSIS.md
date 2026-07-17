# Explicit Apex Bridge P0v2 PATH Diagnosis

Date: 2026-07-17 KST

## Formal Verdict

`REJECT_BRIDGE` under the frozen P0v2 protocol. The first two cells,
`w4a_uncarved_s101_bound` and `w4a_uncarved_s127_bound`, passed. The third
cell, `w4a_uncarved_s163_bound`, failed the absolute interior force-identity
tolerance, and the runner stopped before all three free cells as preregistered.

PATH cause: **`SOURCE_INPUT_PATHOLOGY`**.

The bridge algebra itself reproduced both source endpoints exactly and
returned the exact interior derivative. The scientific blocker is a severe
water overlap in the uncarved source state-1 coordinate image, not an observed
free-energy result and not evidence that the bridge term is zero.

## Run Record

- inventory digest:
  `8a5c2eabd6125bf70a8dbf6b334d29ded9d00091c6aab597e0c4af1c22cd3615`;
- KEEPER: 6/6 context-free PASS, zero MD/GPU, peak RSS about 2.03 GiB;
- RUNNER: 3/6 cells evaluated, zero MD/minimization/GPU, peak RSS about
  2.98 GiB, wall time 13 minutes;
- no free energy was estimated;
- all worker processes exited and returned their OpenMM memory.

| cell | verdict | Hminus-Hplus (kJ/mol) | max force gap (kJ/mol/nm) |
|---|---|---:|---:|
| s101 bound | PASS | 1599.9153046575375 | 84923.01350708821 |
| s127 bound | PASS | 2042.9575286796317 | 132518.30820455504 |
| s163 bound | REJECT | 1331583646.915236 | 276209950501.09717 |

For s101 and s127, endpoint energy/force error was exactly zero, interior
energy error was exactly zero, maximum interior force error was at most
`5.82e-11 kJ/mol/nm`, and `dV/dBridgeXi` error was exactly zero.

For s163, endpoint energy/force error and derivative error were still exactly
zero. The linear-force errors were `1.52587890625e-5` at `xi=0.25` and
`3.0517578125e-5 kJ/mol/nm` at `xi=0.75`, just above the frozen `1e-5`
absolute tolerance. These are machine-rounding residuals at force scales of
`6.9e10-2.1e11 kJ/mol/nm`; relaxing the gate would not make the source state
physically usable.

## Source Pathology

OpenMM `ATMForce.getPerturbationEnergy()` at the stored s163 bound coordinates
returned:

```text
u0 =      -3386060.307565201 kJ/mol
u1 =    2659781997.345677000 kJ/mol
Hminus = 1328197968.519055800 kJ/mol
```

The largest force is on transformed W4A TRP `B:4:CZ2` (index 10154) and the
opposing atom is water oxygen `4:-15244:O` (index 161169). Applying the frozen
state-1 `ParticleOffsetDisplacement` places them only
`0.05695665018239689 nm` apart under the triclinic minimum-image convention.
The large opposite forces on the ALA/TRP anchor carbons are the same collision
propagated through the particle-offset reference atoms.

The dynamic ghost cannot repair this initial clash. It acts at the stored
coordinates of the disappearing TRP site, while this water overlaps the TRP
after it is transformed into the appearing site.

## Six-Cell Geometry Audit

The same context-free state-1 ring/water audit found:

| cell | uncarved min (nm) | pairs <0.20 nm | carved min (nm) | carved pairs <0.26 nm |
|---|---:|---:|---:|---:|
| s101 bound | 0.25748 | 0 | 0.33929 | 0 |
| s127 bound | 0.24472 | 0 | 0.33300 | 0 |
| s163 bound | 0.05696 | 5 | 0.27641 | 0 |
| s101 free | 0.23996 | 0 | 0.26273 | 0 |
| s127 free | 0.10183 | 1 | 0.32097 | 0 |
| s163 free | 0.04850 | 5 | 0.26025 | 0 |

Existing carved sources remove the immediate state-1 overlaps, but they delete
2-6 waters depending on seed and leg. Their thermodynamic neutrality is not
demonstrated, and the completed DCD gate already showed that a static initial
carve does not prevent dynamic water re-entry. Therefore the carved sources
are a useful mechanics comparator, not an automatically valid FE route.

## Scientific Interpretation

P0v2 establishes three points:

1. An explicit `Hplus -> Hminus` bridge can reproduce the existing direction-
   specific apex Hamiltonians exactly and exposes a consumed, exact
   `dV/dBridgeXi` observable.
2. The endpoint gap is large, seed dependent, and cannot be assumed to cancel.
3. The uncarved source arm is not a valid starting ensemble for this bridge;
   at least s163 bound/free and likely s127 free contain appearing-site water
   clashes before sampling begins.

This is not a quantitative Delta G verdict. P0v2 sampled no coordinates and
estimated no free energy.

## Ranked Next Design

Return to SCIVAL before implementation and preregister a source-generation
P0v3 with all three components:

1. **Water-count-preserving union solvation.** Build the initial solvent
   exclusion against both ATM coordinate images, then place any displaced
   waters in bulk so target water count/density is preserved. Gate every source
   on transformed ring-water minimum distance before Context creation.
2. **Dynamic disappearing-site ghost.** Retain the water-only ghost to prevent
   later cavity refill; a static union-solvated start does not replace it.
3. **Explicit sampled apex bridge.** Retain the exact endpoint bridge and its
   separately reported bound/free MBAR term.

Existing carved inputs may be used only as a CPU mechanics comparator or as a
construction reference. They must not be promoted to an FE route without a
declared water-count/density policy and correction or closure evidence.

No gate relaxation, free-cell resume, sampled bridge launch, P1/P2, B tuning,
forward-only `dgbind1`, `base=u1` capping, W23A launch, commit, or push is
authorized by this result.
