# Count-Preserving Union Solvation P0v3 Revision 1

Date: 2026-07-17 KST

Status: **FROZEN BEFORE OFFICIAL SOURCE OR P0v3r1 OUTPUT**

Scientific regime: UPDD R-11 ranking/SIGN-only mechanism validation.

## Inheritance

This revision inherits every cohort, count, RNG seed, source declaration,
placeholder definition, dynamic-ghost contract, apex-bridge equation,
execution boundary, tolerance, outcome, and stop condition from:

```text
analysis/dynamic_ghost_union_bridge_p0v3_20260717/PREREGISTRATION.md
SHA256 b768479b1c007cc8794349c39f6ed209c39ae742a58fd4ce7556f4ea2985a4bc
```

except for the single U2 density clause replaced below. The superseded design
and its preflight stop remain immutable evidence.

SCIVAL verdict: **CONDITIONAL APPROVE FOR P0v3r1 ONLY**.

No MD, minimization, GPU, FE estimate, sampled bridge, production, or W23A
launch is authorized.

## Reason For Revision

The pre-output s101-free development smoke proved that the implementation can
preserve all frozen counts exactly and remove every placeholder:

```text
particles=29992, water=9846, Na=27, Cl=27, total solvent=9900
placeholder target drift=0.0 nm, placeholders remaining=0
```

OpenMM `addSolvent(numAdded=9900)` selected a `6.8498 nm` cubic edge instead of
the padding-built parent's `6.7890 nm`. The superseded U2 gate compared solvent
molecule count to total box volume and would fail at 2.64%.

That ratio is not a valid static bulk-density estimator because total box
volume includes the solute and the deliberately enlarged union excluded
volume. Equilibrium density cannot be established without NPT sampling, which
P0v3r1 explicitly forbids. Changing a placeholder radius or molecule count to
force this malformed ratio under 1% would be post-result parameter fitting and
is prohibited.

## Replaced U2 Clause

Remove only this superseded requirement:

```text
total-solvent-molecule density differs from the parent by at most 1%
```

Replace it with:

- final box vectors and volume are finite, right-handed, cubic, and reported;
- parent and union total-box `N/V` values and relative difference are reported
  as construction diagnostics only, with no P0v3r1 pass/fail threshold;
- no claim of equilibrium or density equivalence is allowed;
- a future sampled pilot must include a separately preregistered NPT density
  and box-stability gate before any FE interpretation.

All exact count gates remain unchanged. In particular, every cell must retain
the frozen water, Na, Cl, total solvent, atom, and System particle counts. Both
stored and transformed ring images must still have zero water-oxygen pairs
below `0.26 nm`.

## P0v3r1 Outcome Boundary

`P0V3_PASS` means only that count-preserving union construction, final artifact
integrity, dynamic-ghost scope, and explicit apex-bridge identities pass on CPU
for all six cells. It does not certify equilibrium density, sampling overlap,
free energy, or Delta-G neutrality. A PASS permits preparation of a separate
NPT/sampled-bridge preregistration and no automatic launch.

