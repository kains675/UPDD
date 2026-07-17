# P1SD Preregistration Amendment R1

Frozen locally on 2026-07-17 after the CUDA development smoke and before any
official P1SD output.

## Scope

This amendment supersedes only the mixed-precision explicit cross-state energy
readback tolerance in `PREREGISTRATION.md`. The original preregistration,
protocol, freeze manifest, and failed development output remain unchanged.

No sampling length, xi grid, seed, thermodynamic state, force construction,
parent box, physical gate, MBAR setting, or decision threshold changes.

## Trigger Evidence

The isolated `w4a_union_s101_free`, xi `0.50` CUDA development window passed
all physical, trajectory, temperature, fixed-box, and finite-value checks. It
failed only the original absolute `1e-5 kJ/mol` cross-state readback gate:

- first sample maximum error: `5.207047797739506e-5 kJ/mol`;
- final sample maximum error: `4.199275281280279e-5 kJ/mol`;
- a nonzero-gap source-minimized diagnostic reached
  `1.5489209908992052e-4 kJ/mol`;
- CUDA double precision reproduced the linear identity at the same final
  coordinates to the reported energy precision.

The failed development window result SHA256 is
`512dd71ffc4228b4a1e3915e3c4bb6e2967d17ba22318e259da842c812ef7566`.
Its samples are development-only and are prohibited from official P1SD reuse.

## Revised Readback Gate

For each observed/expected energy pair, define:

```text
energy_scale = max(1 kJ/mol, abs(observed), abs(expected))
allowed_error = min(1e-3, max(1e-5, 1e-9 * energy_scale)) kJ/mol
```

Every explicit first/final cross-state readback must satisfy its own
`abs(observed-expected) <= allowed_error` gate. The `1e-3 kJ/mol` ceiling is
about `4e-4 kT` at 300 K and remains far below a physically meaningful bridge
error while accommodating CUDA mixed-precision reduction noise in total
energies near `-6e5 kJ/mol`.

## Scientific Review Verdict

`CONDITIONAL_APPROVE_P1SD_R1_ONLY`.

Conditions remain unchanged: complete the six-cell integrity and Reference
preflights, stop on the first physical or execution failure, treat all MBAR
values as exploratory schedule diagnostics, and do not launch a densified or
production schedule automatically.
