# ATM-Nested Dynamic Ghost P0 SciVal Verdict

Registered: 2026-07-17 KST, before implementation or P0 output

Verdict: `CONDITIONAL APPROVE`

## Scientific Question

Test whether the disappearing W4A TRP ring can carry the same water-only
soft-core WCA term inside the canonical `ATMForce`, so that the term is
evaluated at stored coordinates in `u0` and at the declared ATM-transformed
coordinates in `u1`.

This is a Hamiltonian-wiring audit. It is not a free-energy calculation and
does not establish equilibrium solvent density, thermodynamic neutrality, or
production readiness.

## Conditions

1. Use the six completed P1SD normalized union-source cells without rebuilding
   their solvent, box, topology, coordinates, or canonical nested forces.
2. Add exactly one water-only `CustomNonbondedForce` inside the existing
   `ATMForce`; do not add particles, delete waters, alter charges, or add
   electrostatics or attraction.
3. Preserve the top-level force count, all pre-existing nested-force
   serializations, ATM energy expression, and every particle transformation.
4. Run OpenMM `Reference` only, in one isolated subprocess per cell, with zero
   MD and zero minimization steps.
5. Require exact-scope selection, `g=0` source energy/force parity, finite
   derivatives, serialization round-trip integrity, and state-specific
   stored/transformed probe activation.
6. Stop with `REJECT_GHOST` on any declaration drift, non-finite value,
   opposite-endpoint leakage above the frozen tolerance, incorrect force
   direction, or failed parity.

## Prohibited Claims

- no absolute or relative free-energy estimate;
- no claim that the term is delta-G neutral;
- no reuse of probe coordinates as sampled configurations;
- no authorization of CUDA sampling or production.

