# P1CAP Preregistration Amendment R1

Status: `FROZEN_AFTER_R0_DCD_UNIT_FAILURE_BEFORE_R1_OUTPUT`

Date: 2026-07-17

## Trigger

P1CAP R0 stopped at the first cell with a formal Runner failure. The low-level
`mdtraj.formats.DCDTrajectoryFile` reader declares `distance_unit` as
`angstroms` and returned the first bound unit-cell length as `145.795`.
The R0 worker passed those raw values to OpenMM as nanometers instead of
converting them to the expected `14.5795 nm`.

All five attempted R0 windows therefore failed the fixed-box and same-Context
algebra gates. Those values are unit-corrupted implementation output and have
no scientific interpretation. No later cell was launched and no R0 result may
be reused.

R0 inventory, KEEPER result, run state, and worker output are preserved under
their R0 names and `_archive/20260717_2037_r0_dcd_angstrom_unit_fail/`.

## R1 Correction

R1 retains the original cohort, hypotheses, soft-core values, 5 kJ/mol
classification guard, slope tolerance, zero-MD contract, and claim
prohibitions unchanged.

The only runtime correction is an explicit unit contract:

```text
DCDTrajectoryFile.distance_unit must equal "angstroms"
xyz_nm         = xyz_raw * 0.1
box_lengths_nm = box_lengths_raw * 0.1
angles_deg     = angles_raw
```

No scale is inferred from coordinate magnitude. Any other declared DCD unit is
a formal KEEPER/Runner stop.

R1 KEEPER must read the first frame of every DCD, verify the declared unit, and
confirm that converted box vectors match the frozen parent box within
`1e-4 nm`. R1 workers must independently re-check the declared unit before
replay.

R1 uses new inventory, KEEPER, run-state, replay, and summary paths. It may run
only after this amendment and `protocol_r1.json` are hash-frozen in
`FROZEN_MANIFEST_R1.json`.

