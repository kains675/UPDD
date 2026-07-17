# Dynamic Ghost P0 PATH Diagnosis

Date: 2026-07-17

## Verdict

`REJECT_GHOST` under the frozen P0 protocol. The first cell,
`w4a_uncarved_s101_bound`, failed the shared-apex Hamiltonian parity gate. The
runner stopped immediately, so the remaining five cells and all P1/P2 work
were not launched.

This is not a floating-point or serialization failure. It is a deterministic
consequence of the direction-dependent ATM soft-core Hamiltonian when
`UOffset=0` and the raw perturbation exceeds `Ubcore`.

## Raw Result

- Inventory digest:
  `e6f17322677d097e1b4c39af2e9f8ca19aaba3006a8866e37936e5653212a181`
- Cell elapsed time: `3711.924379945005 s`
- Shared-apex energy difference: `1599.9153046575375 kJ/mol`
- Shared-apex maximum force-component difference:
  `84923.01350708821 kJ/mol/nm`
- dplus apex energy: `-2984820.0572463754 kJ/mol`
- dminus apex energy: `-2983220.141941718 kJ/mol`
- Ghost force-group energy at both apexes: `0.0 kJ/mol`
- `dV/dg` at both apexes: `0.0 kJ/mol`

All other frozen checks in this cell passed. The `u0,g=0`, physical
`u1,g=0`, serialized main/correction `u1,g=1`, and one-step minimized
`u0,g=0` energy/force differences were exactly zero. The minimized position
difference was also exactly zero. Force scope, source parameters, atom count,
charge, finite values, force-group identity, and serialization round-trip all
passed.

## Deterministic Mechanism

At the stored s101 bound coordinates, Reference readback gave:

```text
u0                         = -2985094.923625607975 kJ/mol
u1                         = -2981345.360257827677 kJ/mol
Delta u = u1-u0            =     3749.563367780298 kJ/mol
softcore(+Delta u)         =      549.732758465654 kJ/mol
softcore(-Delta u)         =    -3749.563367780298 kJ/mol
```

For `Lambda1=Lambda2=0.5` and `UOffset=0`, the existing ATM expression reduces
to:

```text
H_apex,+ = u0 + 0.5*softcore(+Delta u)
H_apex,- = u1 + 0.5*softcore(-Delta u)
```

The positive perturbation is capped because it exceeds `Ubcore`, while the
negative perturbation is not capped. Therefore:

```text
H_apex,+ - H_apex,-
  = -0.5*(Delta u-softcore(+Delta u))
  = -1599.915304657537 kJ/mol
```

This reproduces the observed apex difference to numerical precision. The
dynamic ghost contributes zero at these stored coordinates, so it neither
causes nor repairs the mismatch.

## Interpretation And Boundary

The frozen assumption that the dplus and dminus `0.5/0.5` states are one
shared Hamiltonian is false for this source Hamiltonian. Removing or weakening
the apex gate would hide a coordinate-dependent thermodynamic gap and is not
an acceptable repair.

Any revision requires a new preregistration. The minimum acceptable design is
either:

1. an explicit, sampled apex bridge whose endpoints reproduce the current
   dplus and dminus apex energies and forces exactly; or
2. a direction-independent common-midpoint Hamiltonian with exact physical
   endpoint parity and independently demonstrated cycle closure.

The explicit apex bridge is the less invasive candidate because it preserves
the existing physical endpoints and measures the non-constant apex free-energy
term instead of assuming cancellation. A revised CPU-only P0 must pass before
any structural or thermodynamic GPU pilot.

The failure does not authorize B transport tuning, forward-only `dgbind1`,
`base=u1` capping, an uncorrected cavity potential, P1/P2, or W23A production.
