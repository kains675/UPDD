# P1SD PATH Diagnosis

Status: `FAIL_P1SD_APEX_DEGENERACY_AND_TRANSFORMED_VOID_PENETRATION`

The official P1 schedule-discovery run stopped at the preregistered first
failure in `w4a_union_s101_free`, xi `0.00`. The six-cell normalization,
KEEPER, and CPU Reference stages passed. Fifteen bound windows passed their
execution gates; the first free window returned `WINDOW_FAIL`. The remaining
14 windows were not launched.

No production free energy, schedule selection, densification authorization,
round-trip claim, mutation verdict, or reusable production sample exists.

## Formal Runner Failure

The final explicit cross-state readback had maximum absolute error
`1.531239e-3 kJ/mol`, above the frozen `1e-3 kJ/mol` ceiling. This is a formal
RUNNER failure and is not relabeled as a pass.

The five signed target-state errors were nearly common-mode. Their mean was
`-1.490705e-3 kJ/mol`; after subtracting that mean, the largest residual was
`4.616871e-5 kJ/mol`. The relative linear bridge relation therefore remained
numerically intact. A future implementation should re-evaluate the restored
generating state and gate relative `Delta U`, while retaining absolute
total-energy drift as a diagnostic. That correction does not authorize reuse
or resumption of this run.

## Root Finding 1: The Apex Bridge Is Algebraically Flat

Let `Delta = u1-u0` with `UOffset=0`. In the ordinary soft-core linear branch,
`abs(Delta) < Ubcore = 418.4 kJ/mol`:

```text
uscplus  = Delta
uscminus = -Delta
Hplus    = u0 + 0.5*Delta
Hminus   = u1 - 0.5*Delta
         = (u0+u1)/2
```

Thus `Hplus == Hminus` exactly and the bridge slope is zero. The two
directional apex Hamiltonians differ only while the soft-core cap is active.

This is observed directly in the official bound samples:

- `s101_bound`: 0/250 nonzero slopes;
- `s127_bound`: 5/250 nonzero slopes, all in the initial `xi=0` transient;
- `s163_bound`: 0/250 nonzero slopes.

Overall, 745/750 bound samples were exactly flat. The apparent overlap of
`0.2` is therefore state identity, not evidence that a meaningful bridge is
well mixed. Midpoint insertion cannot recover a signal from an algebraically
flat coordinate. The exploratory bound MBAR values are not physical mutation
free energies.

## Root Finding 2: Transformed-Site Water Penetration

The failed `s101_free/xi=0` window remained finite and temperature-stable, but
the transformed ring-water distance reached `0.1867466 nm`. Thirty-nine of 50
frames were below `0.26 nm`, involving nine distinct waters. `HOH340` contacted
the transformed ring from the first frame through the final frame; the global
minimum was to transformed `NE1` at frame 2.

The stored ring remained protected (`0.2796548 nm` minimum). This state
specificity matches the implementation: the top-level dynamic ghost interacts
with water using the stored coordinates of the disappearing ring atoms. It has
no force site at `r_ring + r_destination - r_origin`, so it cannot prevent
water entry at the ATM-transformed image.

This is structural, not a single-water outlier and not repairable by xi
densification. The smaller bound warning in `s163` (`0.25960/0.25909 nm`
stored/transformed minima) points in the same direction.

## Decision

Do not resume the remaining free windows and do not densify this apex bridge.
Preserve all partial artifacts as rejected pilot evidence.

The next preregistration should separate two questions:

1. Replace apex-to-apex free-energy sampling with a cap-activation diagnostic
   based on the raw `u1-u0` distribution and the fraction of configurations
   outside `+/-Ubcore`. Directional apex disagreement is a deterministic
   soft-core-tail diagnostic in this formulation.
2. Promote the transformed-site carve/ghost redesign. The exclusion must act
   at the ATM-transformed ring coordinates and pass per-state DCD penetration
   gates before any new free-energy pilot.

If an explicit thermodynamic bridge is still required later, bridge physical
endpoint Hamiltonians rather than the direction labels at lambda `0.5`, and
freeze a fresh protocol after transformed-site structural validation.

Evidence:

- `P1SD_PARTIAL_SUMMARY.json`
- `P1SD_RUN_STATE.json`
- `sampled_schedule_discovery/w4a_union_s101_free/xi_0p00/window_result.json`
- `sampled_schedule_discovery/w4a_union_s101_free/xi_0p00/trajectory.dcd`
