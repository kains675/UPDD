# P1SD Preregistration Amendment R2

Frozen locally on 2026-07-17 after pre-official code review and before any
official P1SD output.

## Scope

This amendment supersedes only how the linear bridge slope is evaluated for
cross-state energies. R1 remains authoritative for the mixed-precision
readback tolerance. All earlier preregistrations, manifests, and development
outputs remain unchanged.

No sampled Hamiltonian, xi grid, sampling length, seed, thermodynamic state,
parent box, physical gate, MBAR setting, overlap threshold, or reuse policy
changes.

## Trigger Evidence

The exact endpoint-preserving bridge expression uses outer `select` branches.
OpenMM therefore reports `dV/dBridgeXi = 0` at exact xi `0` and `1`, even
though the interior branch has the required constant slope. The frozen
P0v3r2 `w4a_union_s101_bound` Reference result reports:

```text
xi       0.00       0.25       0.50       0.75       1.00
dV/dxi   0.0        357.661    357.661    357.661    0.0  kJ/mol
```

Using the endpoint derivative as the cross-state slope would therefore make
endpoint samples contribute an incorrect constant energy row. This was found
during code review before official inventory creation; no P1SD sample has been
used for a scientific decision.

## Revised Slope Contract

For every saved configuration from every generating xi:

1. retain the potential energy evaluated at the generating xi;
2. set `BridgeXi=0.5` without integrating or changing coordinates;
3. read the interior `dV/dBridgeXi`, which equals `Hminus-Hplus`;
4. restore the generating xi in a `finally` path before any next MD step;
5. use only that interior slope in
   `U_j = U_i + (xi_j-xi_i) * slope`;
6. preserve the derivative read at the generating xi as a diagnostic only.

Every slope must be finite, every restore check must pass, and the R1
first/final explicit cross-state readback gate remains mandatory.

## Scientific Review Verdict

`APPROVE_P1SD_R2_SLOPE_PROBE`.

The probe is an energy query at fixed coordinates, not a dynamics step, and
does not alter the ensemble sampled by a window. Any restore, finite-value, or
explicit-readback failure stops the run.
