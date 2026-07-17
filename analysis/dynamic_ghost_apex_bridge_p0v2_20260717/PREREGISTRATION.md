# Explicit Sampled Apex Bridge P0v2 Preregistration

Date: 2026-07-17 KST

Status: **FROZEN BEFORE BRIDGE IMPLEMENTATION OR P0v2 OUTPUT**

Scientific regime: UPDD R-11 ranking/SIGN-only mechanism validation. This
revision defines an explicit thermodynamic bridge between the existing dplus
and dminus apex Hamiltonians. It does not calibrate an absolute or mutation
free energy and does not authorize a prospective claim.

## Trigger And Decision

The frozen dynamic-ghost P0 stopped at its first cell,
`w4a_uncarved_s101_bound`, because the two nominally shared apex states were
not the same Hamiltonian. On identical stored coordinates:

```text
H_apex,+ - H_apex,- = -1599.9153046575375 kJ/mol
max force-component difference = 84923.01350708821 kJ/mol/nm
```

All other tested endpoint, force-scope, charge, atom-count, finite-value, and
serialization checks passed. The failure is deterministic: with `UOffset=0`,
the current ATM soft-core caps the large positive `u1-u0` perturbation but not
the corresponding negative perturbation.

The failed P0 protocol and outputs remain immutable evidence. This revision
does not relax its shared-apex gate or reinterpret that failure as numerical
noise. Instead, it promotes the apex difference to an explicit free-energy
leg whose endpoints are the two existing apex Hamiltonians.

SCIVAL verdict: **CONDITIONAL APPROVE FOR P0v2 ONLY**.

P0v2 may test the revised Hamiltonian construction on CPU. No MD, GPU,
free-energy estimate, P1/P2 continuation, W23A launch, B-branch transport
tuning, forward-only estimate, or production launch is authorized here.

## Frozen Hamiltonians

For a configuration `x`, let `u0(x)` and `u1(x)` be the two potential values
provided by the existing serialized `ATMForce`. Keep the source constants
unchanged:

```text
UOffset = 0
Umax    = 200 kcal/mol
Ubcore  = 100 kcal/mol
Acore   = 0.0625
```

Let `S(v)` denote the existing ATM positive-tail soft-core function with those
constants. The two source apex Hamiltonians are:

```text
H_plus(x)  = u0(x) + 0.5*S( u1(x) - (u0(x)+UOffset))
H_minus(x) = u1(x) + 0.5*S(-u1(x) + (u0(x)+UOffset))
```

These definitions are not a new cap or a fitted replacement. They are exact
restatements of the current dplus and dminus `Lambda1=Lambda2=0.5` source
Hamiltonians.

Define the bridge coordinate `xi` in `[0,1]` and the bridge Hamiltonian:

```text
H_bridge(x;xi) = (1-xi)*H_plus(x) + xi*H_minus(x)
```

The implementation must use a new consumed global parameter named
`BridgeXi`, declared on the existing `ATMForce`, with an energy-parameter
derivative. It may change only the `ATMForce` energy expression and force
name, and add that parameter/derivative. It must preserve all nested forces,
particles, displacements, particle transformations, force ordering, source
parameters, atom count, and charge.

The dynamic water-only ghost from the rejected protocol is retained without
modification and held at `g=0.5` throughout this bridge. It remains a separate
top-level force after `ATMForce`, in force group 31. Therefore it contributes
identically to both endpoints and is independent of `xi`.

The serialized source globals `Lambda1`, `Lambda2`, `Alpha`, `Uh`, `W0`, and
`Direction` remain present because they belong to the source `ATMForce`, but
the bridge expression must not consume them. `Umax`, `Ubcore`, `Acore`, and
`UOffset` remain consumed. Changing any inactive legacy global at fixed
`BridgeXi` must not change the bridge energy or forces.

## Exact Thermodynamic Cycle

For environment `e` in `{bound, free}`, the revised ghost-assisted route is:

```text
DeltaG_plus,e   = G_e(H_plus,g=0.5)  - G_e(u0,g=0)
DeltaG_bridge,e = G_e(H_minus,g=0.5) - G_e(H_plus,g=0.5)
DeltaG_minus,e  = G_e(u1,g=1)        - G_e(H_minus,g=0.5)
DeltaG_off,e    = G_e(u1,g=0)        - G_e(u1,g=1)

DeltaG_phys,e = DeltaG_plus,e + DeltaG_bridge,e
              + DeltaG_minus,e + DeltaG_off,e
```

The physical binding double difference remains:

```text
DeltaDeltaG_bind,phys = DeltaG_phys,bound - DeltaG_phys,free
```

`DeltaG_bridge,e` is coordinate dependent and must be estimated from explicit
sampling. It is never zeroed, inferred from one configuration, replaced by
the endpoint potential-energy difference, or assumed to cancel between bound
and free. Existing dplus/dminus output conventions remain subordinate to this
full cycle: `dgb = dgbind1-dgbind2` is primary and `dgbind1` is diagnostic.

## Sampled Bridge Contract

The scientific route is preregistered as an explicitly sampled equilibrium
Hamiltonian path. A future sampled pilot must:

- sample both bound and free environments for every declared seed;
- use multiple declared `xi` states including exact `0` and `1` endpoints;
- propagate coordinates under each `H_bridge(x;xi)` rather than post-process
  only the source coordinates;
- evaluate each retained configuration under every declared `xi` state;
- estimate `DeltaG_bridge,e` with MBAR and report it separately by environment
  and seed;
- gate on finite reduced potentials, adjacent overlap, endpoint visitation,
  equilibration, independent-block stability, and uncertainty;
- preserve the dynamic ghost at exactly `g=0.5` for every bridge state;
- add the bridge term with the orientation shown in the frozen cycle.

The numerical `xi` schedule, sampling length, RNG seeds, equilibration cut,
and output root are intentionally not authorized by P0v2. They require a
separate frozen sampled-pilot preregistration after P0v2 passes. This prevents
choosing a schedule from an unverified implementation or silently launching a
large run. P0v2 measures no free energy.

## P0v2 Cohort And Execution Boundary

Use the same six uncarved W4A source cells as the rejected P0:

```text
s101 bound/free  control
s127 bound/free  high leverage
s163 bound/free  control
```

For every cell:

- source structure: `2QKI`;
- mutation: `w4a_trp_ala_res4`;
- endpoint: `wt`;
- construction: `twocopy`;
- platform: OpenMM `Reference` only;
- bridge audit states: `xi = 0,0.25,0.5,0.75,1`;
- ghost coupling: `g=0.5` at all bridge states;
- no MD, minimization, GPU, trajectory, or FE estimator;
- one cell per subprocess so OpenMM memory is returned on exit.

Skipping minimization is deliberate. The rejected P0 already demonstrated
exact `g=0` endpoint and one-step minimization parity. P0v2 changes only the
ATM algebra and tests exact energies and forces directly at all bridge audit
states.

## Frozen Gates

### K0: declaration-to-artifact integrity

Any failure is `REJECT_BRIDGE`:

- all trigger preregistration, diagnosis, summary, cell-result, and source
  inventory hashes match this freeze;
- all six source XML/PDB/manifests and both per-leg H18 schedules match a fresh
  frozen inventory;
- exactly one top-level `ATMForce` exists and its source expression matches the
  canonical existing Track B expression;
- the bridge changes only the permitted `ATMForce` metadata/expression and
  parameter declarations;
- nested-force XML, particle transformations, displacements, force ordering,
  atom count, and charge are invariant;
- exactly one unchanged dynamic ghost exists after `ATMForce` in group 31;
- no context is created during K0.

### K1: Reference Hamiltonian identity

Any failure is `REJECT_BRIDGE`. Tolerances are absolute:

```text
energy                <= 1e-5 kJ/mol
max force component   <= 1e-5 kJ/mol/nm
energy derivative     <= 1e-5 kJ/mol
```

At identical PDB coordinates and `g=0.5`:

- bridge `xi=0` energy and forces equal the existing dplus apex exactly;
- bridge `xi=1` energy and forces equal the existing dminus apex exactly;
- at `xi=0.25,0.5,0.75`, total energy and forces equal the corresponding
  linear combination of separately evaluated source apex values;
- at every interior state, `dV/dBridgeXi = H_minus-H_plus`;
- the reported endpoint gap equals source `H_minus-H_plus`; its magnitude is
  recorded and is not a pass/fail threshold;
- ghost force-group energy and `dV/dg` are finite and invariant across `xi`;
- changing inactive legacy globals at fixed `xi=0.5` changes neither energy
  nor forces within tolerance;
- serialization/deserialization preserves endpoint and midpoint identities;
- every energy, force, and required derivative is finite.

### K2: cycle and execution integrity

Any failure is `REJECT_BRIDGE`:

- all six cells use the same frozen bridge expression and state list;
- bridge orientation is `G(H_minus)-G(H_plus)`;
- P0v2 reports zero MD steps and no GPU use;
- completed-cell artifacts are hash checked before resume or aggregation;
- a failed cell stops the run and no later stage is launched.

## Preregistered Outcomes

1. `P0V2_PASS`: all six cells pass K0/K1/K2. The sampled bridge becomes
   scientifically eligible for a separate sampling preregistration; it does
   not start automatically.
2. `REJECT_BRIDGE`: any declaration, source, endpoint, force, derivative,
   serialization, cycle, or execution gate fails. Preserve all output and
   return to PATH diagnosis.
3. `P0V2_INTERRUPTED`: no scientific verdict. Validate completed artifacts,
   archive only invalid partial outputs, and resume with the frozen inventory.

## Stop Conditions

Stop immediately on source/freeze drift, topology/charge mismatch, endpoint or
interior identity failure, non-finite value, serialization drift, output
collision, unexpected GPU access, any MD step, R-rule concern, sign flip
against an anchor, sigma doubling, tier-boundary crossing, new cohort
heterogeneity, or attempted git push.

No P0v2 result authorizes dropping the bridge term, reusing forward-only
`dgbind1`, setting `base=u1`, weakening the source ATM cap, fitting a bridge
schedule, returning to B transport tuning, or launching W23A.

