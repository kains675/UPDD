# ATM-Nested Dynamic Ghost P0 PATH Diagnosis

Completed: 2026-07-18 KST

PATH verdict: `P0_PASS_ENDPOINT_SELECTIVE_NO_PATHOLOGY`

## Scope

This diagnosis covers the six-cell, zero-MD, zero-minimization OpenMM
`Reference` audit only. It does not estimate a free energy, establish
delta-G neutrality, validate equilibrium solvent density, or authorize
production.

## Integrity

- inventory digest: `7dc94e231613aa8b5fc2ecf5d4dd04607fd539b1f477af93f618695f77c5ec46`
- KEEPER digest: `6f2ba6fd52fce0efddb725f02e4e1867146cbf5ccfe9ce09d77dfdd0672b4830`
- KEEPER: `6/6 PASS`, zero Contexts
- RUNNER: `6/6 P0_PASS`
- platform: `Reference`
- MD/minimization steps: `0/0`
- GPU used: `false`

## Physical And Numerical Checks

Across all six W4A union-source cells:

- `g=0` source energy, raw endpoint energy, and force parity were exactly zero
  within reported Reference precision;
- intended stored/transformed probe energies ranged from
  `2.3705916172` to `2.3959255067 kJ/mol`;
- the maximum within-cell stored/transformed energy mismatch was
  `4.66e-10 kJ/mol`;
- every endpoint-matched opposite-site double-difference was exactly
  `0 kJ/mol`;
- the minimum outward oxygen-force projection was
  `149.9546002635 kJ/mol/nm`;
- the maximum isolated force-balance component was
  `1.16e-11 kJ/mol/nm`;
- all `dV/dg` values were finite;
- all declared ATM ring offsets remained within
  `3.9999777574` to `4.0000267499 nm`.

The probe energy and outward force are positive, finite, and nearly
endpoint-symmetric. The opposite endpoint is inactive after the preregistered
source-baseline subtraction. Force balance is many orders of magnitude below
its gate. These results support correct nested-coordinate wiring rather than a
top-level stored-coordinate artifact.

## Source Contacts

Unmodified normalized source coordinates had endpoint-isolated ghost energies
from `0` to `1.2772110444 kJ/mol` and a maximum isolated force component of
`59.9425717199 kJ/mol/nm`. These are finite near-cutoff source contacts, remain
far below the frozen gross-behavior gates, and do not indicate a new clash
tail. R1 preserves their absolute values as diagnostics and removes them only
from the probe-specific leakage metric.

## Pathology Assessment

- no sign inversion or endpoint swap;
- no opposite-site leakage;
- no force imbalance;
- no serialization or declaration drift;
- no bound/free or seed-specific gate heterogeneity;
- no non-finite derivative, energy, or force;
- no evidence that adding the nested force globally flattens the Hamiltonian.

The P0 Hamiltonian-wiring question passes. The remaining uncertainty is
trajectory-level: this point-coordinate audit does not show whether the nested
term activates on the transformed-site penetrations observed in the failed
P1SD DCD frames or whether it is sufficient to remove the raw-gap cap tail.

## Ranked Next Step

Preregister a zero-MD, state-specific DCD replay using the preserved P1SD/P1CAP
frames. On each frame, compare source and ATM-nested systems at matched
endpoints and record:

1. stored/transformed ring-water minimum distances;
2. endpoint-specific nested-ghost energy and force activation;
3. raw `u1-u0` gap before and after the nested term;
4. soft-core cap classification before and after;
5. the failed `w4a_union_s101_free`, `xi=0.00` cohort separately from the
   bound controls.

That replay remains diagnostic and cannot be used as equilibrium sampling
under the new Hamiltonian. Only after a separate replay PASS should a fresh
sampled schedule-discovery run be preregistered.

