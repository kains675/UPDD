# Count-Preserving Union Solvation And Apex Bridge P0v3 Preregistration

Date: 2026-07-17 KST

Status: **FROZEN BEFORE UNION IMPLEMENTATION OR P0v3 OUTPUT**

Scientific regime: UPDD R-11 ranking/SIGN-only mechanism validation. P0v3 is
a source-construction and Hamiltonian-integrity audit. It does not estimate a
free energy, establish equilibrium, calibrate a mutation effect, or authorize
a prospective claim.

## Trigger And Decision

Explicit apex-bridge P0v2 stopped at `w4a_uncarved_s163_bound`. The bridge
reproduced both source apex endpoints and `dV/dBridgeXi` exactly, but the
uncarved state-1 image placed a transformed TRP ring atom only
`0.05695665018239689 nm` from a water oxygen. Raw `u1` was
`2659781997.345677 kJ/mol`. Context-free geometry also found a `0.04850 nm`
contact in s163 free and a `0.10183 nm` contact in s127 free.

Existing carved sources remove these initial contacts, but they delete 2-6
waters depending on seed and leg and are not proven Delta-G neutral. Short DCD
evidence also shows dynamic water re-entry after a static carve. P0v3 therefore
separates three responsibilities:

1. count-preserving union solvation prevents initial overlap in both ATM
   coordinate images;
2. the dynamic ghost prevents later refill of the disappearing-site cavity;
3. the explicit apex bridge retains and eventually measures the nonzero
   `Hplus -> Hminus` free-energy term.

SCIVAL verdict: **CONDITIONAL APPROVE FOR P0v3 ONLY**.

No sampled bridge, MD, minimization, GPU, P1/P2, production cohort, W23A
launch, or FE verdict is authorized. Union placement is not assumed to be an
equilibrated or free-energy-neutral ensemble merely because molecule count is
preserved.

## Frozen Cohort And Parent Counts

Use the same W4A source identities as P0v2:

| cell | replicate | target solvent molecules | water | Na | Cl | RNG seed |
|---|---:|---:|---:|---:|---:|---:|
| s101 bound | 5 | 97060 | 96534 | 263 | 263 | 2026071711 |
| s127 bound | 6 | 97830 | 97300 | 265 | 265 | 2026071712 |
| s163 bound | 7 | 101585 | 101035 | 275 | 275 | 2026071713 |
| s101 free | 5 | 9900 | 9846 | 27 | 27 | 2026071721 |
| s127 free | 6 | 10066 | 10012 | 27 | 27 | 2026071722 |
| s163 free | 7 | 10600 | 10542 | 29 | 29 | 2026071723 |

`target solvent molecules` is the exact number of water plus monatomic ion
residues in the corresponding frozen uncarved source. Final water and ion
counts must each match, not only their sum.

## Union-Solvation Construction

### Canonical solute

Rebuild each source from its frozen 2QKI scaffold using the existing canonical
two-copy W4A construction:

```text
mutation_spec       = w4a_trp_ala_res4
construction        = twocopy
displacement_nm     = 4.0
auto_search         = false
padding reference   = 1.2 nm
water model         = tip3p-fb through the existing amber14 stack
ionic strength      = 0.15 M NaCl, neutralize=true
constraints         = none
```

No source charge, atom type, bonded term, ATM partition, displacement, or
soft-core constant may be changed.

### Temporary union placeholders

Before solvent placement, identify the nine disappearing TRP heavy atoms:

```text
CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2
```

For each atom, construct one temporary one-particle residue at the coordinate
obtained by moving that atom into its partner ATM image with the exact frozen
two-copy displacement vector. Each placeholder:

- has the source atom's element, mass, `sigma`, and `epsilon` read from a
  canonical unsolvated `NonbondedForce`;
- has charge exactly zero;
- has no bond, angle, torsion, constraint, virtual site, or persistent force;
- is used only by OpenMM `Modeller.addSolvent` to compute the standard
  vdW-radius exclusion;
- is recorded with its source atom index/name and target coordinate.

No placeholder radius multiplier, fitted cutoff, post-result adjustment, or
hard-coded force-field parameter is allowed.

### Exact count preservation

Call `Modeller.addSolvent` with the cell-specific frozen `numAdded` total above,
the standard ionic-strength arguments, cubic box shape, and the fixed cell RNG
seed. Save and restore Python's global RNG state around this call so no external
random stream is consumed.

After solvent and ions are placed, delete every temporary placeholder residue
before the final canonical `ForceField.createSystem` call. The final topology
and System must contain no placeholder atom, residue, type, force, mass, or
parameter. The final solvent/water/ion/particle counts must exactly equal the
corresponding uncarved parent counts.

This method preserves molecule count while allowing OpenMM to enlarge the box
slightly if the union excluded volume requires it. It does not move selected
waters by hand and does not delete a molecule from the final system.

## Final P0v3 Hamiltonian

On each count-preserving union source:

1. append the frozen dynamic water-only ghost unchanged;
2. hold ghost coupling at `g=0.5` for the apex bridge;
3. replace only the ATM energy function/name and append consumed `BridgeXi` and
   `dV/dBridgeXi`, exactly as in P0v2;
4. evaluate `xi = 0,0.25,0.5,0.75,1`.

The thermodynamic cycle remains:

```text
DeltaG_phys,e = DeltaG_plus,e + DeltaG_bridge,e
              + DeltaG_minus,e + DeltaG_off,e

DeltaG_bridge,e = G_e(Hminus,g=0.5) - G_e(Hplus,g=0.5)
```

No bridge term may be set to zero or inferred from a single-coordinate energy
gap. P0v3 estimates none of these free energies.

## Execution Stages

### U0: isolated source builds

Build all six union sources in separate subprocesses. Write serialized XML,
matching PDB, build report, input hashes, exact solvent counts, box vectors,
placeholder declarations, and final geometry. A source build failure stops
before KEEPER.

### K0: context-free KEEPER

With zero OpenMM Contexts, validate all six source artifacts, final absence of
placeholders, exact source counts, topology/System lockstep, charge, canonical
ATM expression, particle transformations, dynamic-ghost scope, bridge mutation
scope, and serialization hashes.

### P0v3: Reference Hamiltonian audit

Run one cell per subprocess on OpenMM `Reference`. No MD, minimization, GPU,
trajectory, or FE estimator is allowed. Stop at the first failed cell.

## Frozen Gates

### U1: source identity and count

Any failure is `REJECT_UNION_SOURCE`:

- exact parent seed/leg/replicate/scaffold hashes and W4A declaration;
- final water, Na, Cl, total solvent molecule, atom, and System particle counts
  equal the frozen parent exactly;
- final net charge equals the parent within `1e-8 e`;
- no placeholder name/type/residue/particle remains;
- canonical solute atom/residue ordering is unchanged;
- each placeholder parameter equals its source particle parameter exactly;
- placeholder target coordinates equal the final ATM-transformed source atom
  coordinates within `1e-6 nm` after translation alignment;
- fixed RNG seed is recorded and actually consumed;
- no static carve call or water deletion occurs.

### U2: box and geometry

Any failure is `REJECT_UNION_SOURCE`:

- final box is finite, right-handed, cubic, and periodic;
- total-solvent-molecule density differs from the parent by at most 1%;
- both stored and transformed disappearing-ring images have minimum-image
  water-oxygen distance at least `0.26 nm`;
- both images have zero ring/water-oxygen pairs below `0.26 nm`;
- no solute heavy-atom pair or periodic image violates the existing canonical
  two-copy separation gates.

The `0.26 nm` line is inherited from the pre-existing W4A carve and DCD gates;
it is not selected from P0v3 output.

### K1: final force and serialization integrity

Any failure is `REJECT_UNION_SOURCE` or `REJECT_BRIDGE` as applicable:

- exactly one canonical source `ATMForce`, then one unchanged dynamic ghost;
- source nested-force XML, atom count, charge, particle transformations, and
  force ordering are invariant under bridge construction;
- final serialized systems round-trip exactly under structured inspection;
- all energy, force, and required derivative values are finite;
- raw `|u0|` and `|u1|` are below `1e8 kJ/mol`, and maximum absolute force
  component is below `1e8 kJ/mol/nm`; these are gross-clash guards, not FE gates.

### K2: bridge identity

Absolute tolerances are unchanged from P0v2:

```text
energy              <= 1e-5 kJ/mol
max force component <= 1e-5 kJ/mol/nm
dV/dBridgeXi        <= 1e-5 kJ/mol
```

- `xi=0` equals the existing dplus apex energy and forces;
- `xi=1` equals the existing dminus apex energy and forces;
- interior energy and forces equal the declared linear combination;
- interior `dV/dBridgeXi = Hminus-Hplus`;
- ghost energy and `dV/dg` are finite and invariant across `xi`;
- inactive legacy ATM globals do not affect bridge energy or forces;
- serialization preserves endpoint and midpoint values.

### K3: execution integrity

- six source-build subprocesses and six Reference cell subprocesses are the
  maximum authorized scope;
- completed artifacts are hash checked before resume;
- MD steps, minimization steps, GPU use, and FE estimates are all exactly zero;
- any failed cell stops the stage and no sampled pilot starts automatically.

## Preregistered Outcomes

1. `P0V3_PASS`: U1/U2/K1/K2/K3 pass for all six cells. Prepare a separate
   sampled-bridge preregistration; do not launch it.
2. `REJECT_UNION_SOURCE`: count, placeholder, topology, charge, box, geometry,
   raw endpoint, or construction integrity fails.
3. `REJECT_BRIDGE`: source construction passes but final bridge Hamiltonian or
   derivative identity fails.
4. `P0V3_INTERRUPTED`: no scientific verdict; preserve valid completed output
   and resume only with the same frozen inventory.

## Stop Conditions

Stop on any frozen-file/source/code drift, placeholder residue in a final
artifact, count or ion mismatch, geometry failure, non-finite value, gross raw
endpoint clash, bridge identity failure, unexpected Context during KEEPER,
MD/minimization/GPU access, output collision, R-rule concern, sign flip against
an anchor, sigma doubling, tier-boundary crossing, new cohort heterogeneity,
or attempted git push.

No P0v3 result authorizes treating union placement as equilibrium, omitting
the ghost/deghost correction, dropping the apex bridge, resuming B tuning,
using forward-only `dgbind1`, setting `base=u1`, or launching W23A.

