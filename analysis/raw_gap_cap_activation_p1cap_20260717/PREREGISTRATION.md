# P1CAP Raw-Gap Cap-Activation Diagnostic Preregistration

Status: `FROZEN_BEFORE_IMPLEMENTATION_OR_P1CAP_OUTPUT`

Date: 2026-07-17

## Purpose

P1SD established that the directional apex bridge is algebraically flat when
the raw ATM perturbation gap is inside the ordinary soft-core branch. P1CAP
will directly measure `Delta_raw = u1-u0` on every saved frame from the
completed P1SD trajectories and determine when either directional soft-core
cap was active.

This is a pathology diagnostic on rejected pilot samples. It is not a
free-energy calculation, does not rehabilitate P1SD, and cannot authorize
production sampling, densification, mutation ranking, or a sign claim.

## Frozen Cohort

Only the 16 completed P1SD windows are eligible:

- `w4a_union_s101_bound`: xi `0.00, 0.25, 0.50, 0.75, 1.00`
- `w4a_union_s127_bound`: xi `0.00, 0.25, 0.50, 0.75, 1.00`
- `w4a_union_s163_bound`: xi `0.00, 0.25, 0.50, 0.75, 1.00`
- `w4a_union_s101_free`: xi `0.00`

The expected cohort is 4 cells, 16 windows, and 800 DCD frames. The 14 P1SD
windows that were not launched are outside scope and must not be imputed or
described as sampled.

## Computation

For each DCD frame, the frozen normalized bridge system is evaluated with
OpenMM CUDA mixed precision and no integration or minimization. The existing
`ATMForce.getPerturbationEnergy(context)` return values define:

```text
Delta_raw = u1 - u0
```

The actual bridge globals must be read back from the Context. The frozen
expected values are:

```text
UOffset = 0 kJ/mol
Ubcore  = 418.4 kJ/mol
Umax    = 836.8 kJ/mol
Acore   = 0.0625
```

For directional input `x`, the declared soft-core map is:

```text
softcore(x) = x                                      if x <= Ubcore
softcore(x) = (Umax-Ubcore)*f(x)+Ubcore             if x > Ubcore
```

using the exact ATM `y`, `z`, and `f` formula. The analytic apex-bridge slope
is:

```text
Hminus-Hplus =
    Delta_raw + 0.5*(softcore(-Delta_raw)-softcore(Delta_raw))
```

It is exactly zero in the ordinary branch.

## Numerical Classification

DCD coordinates are float32 trajectory records. To avoid converting
serialization noise at the branch boundary into a scientific claim, cap
classification uses a frozen `5.0 kJ/mol` guard band around `|Ubcore|`:

- `linear`: `|Delta_raw| <= Ubcore - 5.0`
- `ambiguous`: `Ubcore - 5.0 < |Delta_raw| < Ubcore + 5.0`
- `cap_active`: `|Delta_raw| >= Ubcore + 5.0`

Positive active gaps identify the dplus cap; negative active gaps identify the
dminus cap. Ambiguous frames are reported separately and never counted as
strict active or strict linear frames.

Each replayed frame must compare the analytic slope with the same Context's
`dV/dBridgeXi`. The frozen tolerance is:

```text
allowed = max(0.05 kJ/mol, 1e-5*max(1, |analytic|, |observed|))
```

The original P1SD slope is retained as a serialization diagnostic. It is not
used to reconstruct `Delta_raw`. The final DCD frame is also re-evaluated from
the corresponding float64 `final_positions.npy`; strict classifications must
agree unless either value is ambiguous.

## Execution Contract

- Absolute interpreter: `/home/san/miniconda3/envs/atm/bin/python`.
- Platform: CUDA device 0, mixed precision.
- Exactly one live OpenMM Context per worker.
- One separate worker process per completed cell; process exit is the memory
  release boundary.
- No MD steps, minimization, velocity assignment, checkpoint continuation, or
  coordinate writing.
- Inputs and implementation are hashed in a prebuild inventory after coding.
- Existing P1SD files are read-only.

## Gates

KEEPER passes only if:

1. this preregistration, protocol, frozen manifest, P1SD parent evidence, and
   every used P1SD artifact match their declared hashes;
2. the exact 4-cell/16-window/800-frame cohort is present;
3. normalized topology, System, DCD, sample table, and float64 final positions
   agree on particle/frame counts;
4. all cells declare the frozen ATM soft-core parameters; and
5. no OpenMM Context is created during KEEPER.

RUNNER passes only if:

1. all 800 frame evaluations are finite;
2. all same-Context analytic slope checks pass;
3. every window retains exactly 50 rows in source and replay output;
4. fixed P1SD box vectors are preserved during replay;
5. final DCD/float64 strict classifications agree outside the ambiguity band;
   and
6. all four cell workers exit successfully.

Any artifact drift, non-finite energy, particle/frame mismatch, parameter
drift, algebra mismatch, or strict final-frame classification mismatch is a
formal stop.

## Required Reporting

Report every frame and summarize, without pooling away cell/window identity:

- raw-gap min, median, max, selected quantiles, and sign;
- strict linear, ambiguous, and strict cap-active counts/fractions;
- dplus versus dminus active counts;
- original nonzero-slope versus replay cap-class contingency;
- analytic/Context slope error;
- DCD/float64 final-frame sensitivity; and
- the highest-absolute-gap frames.

PATH must issue one of:

- `CONFIRM_APEX_DEGENERACY`: cap activation is absent or confined to a sparse
  transient and cannot define a useful sampled bridge;
- `CAP_TAIL_PRESENT_DIAGNOSTIC_ONLY`: persistent cap activation exists but is
  only a soft-core-tail diagnostic, not a free-energy coordinate;
- `INDETERMINATE_NUMERICAL`: ambiguity or replay validation prevents a physical
  diagnosis; or
- `FAIL_INTEGRITY`: a formal KEEPER/RUNNER gate failed.

Regardless of the cap result, the already observed transformed-site water
penetration remains an independent structural blocker. A subsequent
excluded-volume redesign must act at ATM-transformed ring coordinates and pass
state-specific DCD penetration gates before any new free-energy pilot.

