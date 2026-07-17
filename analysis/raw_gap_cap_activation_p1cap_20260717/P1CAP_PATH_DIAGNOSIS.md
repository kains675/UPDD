# P1CAP PATH Diagnosis

Status: `CAP_TAIL_PRESENT_DIAGNOSTIC_ONLY`

P1CAP R1 passed all formal gates for 4 cells, 16 completed P1SD windows, and
800 saved frames. It used one CUDA Context per cell, no MD steps, and explicit
low-level DCD Angstrom-to-nanometer conversion. No frame was inside the frozen
5 kJ/mol ambiguity band.

## Bound: Apex Degeneracy Confirmed

Of 750 bound frames, 745 (`99.33%`) were in the ordinary linear branch and
therefore had exactly identical dplus/dminus apex Hamiltonians. The only five
cap-active frames were `s127_bound`, xi `0.00`, cycles 1-5, with positive raw
gaps from `1178.81` to `1769.38 kJ/mol`. Cycle 6 fell to `21.48 kJ/mol`, and
all remaining frames were linear.

This is a short startup tail, not persistent sampling along a useful bridge.
The 14 other bound windows were 50/50 linear. Densifying xi cannot recover a
coordinate that is algebraically flat over nearly all sampled configurations.

## Free: Persistent Positive Cap Tail

The failed `s101_free`, xi `0.00` window was cap-active in all 50 frames, all
in the dplus direction:

```text
raw u1-u0 min     906.11 kJ/mol
raw u1-u0 median 2065.53 kJ/mol
raw u1-u0 max    5211.05 kJ/mol
```

Bound and free samples therefore occupy different algebraic regimes: bound is
almost entirely flat while this free window is persistently soft-core capped.
Directional apex disagreement is a leg-specific tail diagnostic, not a common
physical free-energy coordinate.

## Structural Association

Within the 50-frame free trajectory, raw gap and transformed ring-water
distance had Pearson `r=-0.793` and Spearman `rho=-0.775`. The 39 frames below
`0.26 nm` had mean raw gap `2331.08 kJ/mol`, versus `1556.12 kJ/mol` for the
other 11 frames.

This is exploratory, temporally autocorrelated evidence from one rejected
window, so it does not prove causality. It is nevertheless consistent with the
existing transformed-site penetration diagnosis: closer water contact tracks
a larger positive raw perturbation tail.

## Numerical Integrity

- analytic versus same-Context slope maximum error: `0.001375 kJ/mol`, below
  the frozen `0.05 kJ/mol` floor;
- original versus replay cap-class mismatches: `0/800`;
- original versus replay slope maximum difference: `0.031302 kJ/mol`;
- final DCD versus float64 raw-gap maximum difference: `0.023251 kJ/mol`;
- final DCD versus float64 class mismatches: `0/16`.

R0 is separately preserved as a formal implementation failure: the low-level
DCD reader's Angstrom values were initially passed as nanometers. R0 has no
scientific interpretation and was not reused.

## Decision

Do not resume or densify P1SD, estimate a free energy from its samples, or
reuse them in production.

Promote the transformed-site excluded-volume redesign. The next preregistered
stage should be a zero-MD P0 audit that places the ghost WCA force inside
`ATMForce`, so the same declared interaction is evaluated at stored and
ATM-transformed ring coordinates. It must pass endpoint force/energy,
serialization, declaration, and explicit probe-water exclusion gates before
any short dynamics or free-energy pilot.

