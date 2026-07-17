# Dynamic Ghost / Excluded-Volume Pilot Preregistration

Date: 2026-07-16 KST

Status: **FROZEN BEFORE IMPLEMENTATION OR A-BRANCH OUTPUT**

Scientific regime: UPDD R-11 ranking/SIGN-only. This pilot is a mechanism and
thermodynamic-integrity test. It cannot produce a calibrated mutation free
energy or authorize a prospective claim.

## Decision And Scope

The completed W4A post-densify DCD gate returned
`PROMOTE_A_CARVE_REDESIGN`. B transport tuning is therefore stopped. The next
allowed scientific branch is a small W4A mechanism pilot of a dynamic,
water-only excluded-volume ghost with an explicit deghost free-energy leg.

This preregistration does not authorize:

- additional H18 window tuning or seed expansion under branch B;
- a W23A RBFE launch, despite the completed 1YCR scaffolds;
- an uncorrected permanent cavity potential;
- reuse of `dgbind1` as a physical estimator;
- a full `n=10` A-branch production cohort.

The pilot uses uncarved W4A inputs. No water is deleted. Static carve and
union-solvation remain comparators, not part of the primary A mechanism.

## Frozen Evidence

The branch change is driven by the completed 20-cell W4A structural gate:

- control median ring/site close-contact baseline: `C = 0.334667`;
- high-leverage carved `s127`: `H = 0.497333`;
- `H-C = 0.162667`, `H/C = 1.486056`;
- all nine controls individually exceeded the frozen `0.15` floor;
- formal verdict: `PROMOTE_A_CARVE_REDESIGN`.

The all-site result is supported by an exploratory ring decomposition: the
TRP indole ring, especially dminus states 9 and 10, carries the signal while
site RMSD p95 remains `0.078-0.136 nm`. This is a structural pathology signal,
not an FE estimate. Only 8/20 cells met the minimal mixing check and no cell
completed a round trip, so no FE convergence claim is imported into this
preregistration.

## Scientific Rationale

The two-copy construction has two physical endpoint configurations. In `u0`,
the outgoing TRP copy occupies its native stored coordinates. In `u1`, ATM
evaluates that TRP at the partner site, leaving its stored coordinates as a
dummy-region cavity that water can refill. Static build-time carving can remove
an initial overlap at the partner site but cannot control this later cavity
reorganization.

An auxiliary bias is acceptable only if its free-energy contribution is
explicitly added and removed through a thermodynamic cycle. This follows the
same statistical-mechanical requirement used for restraint corrections:
auxiliary potentials cannot be assumed neutral merely because they improve
sampling (Boresch et al., DOI
<https://doi.org/10.1021/jp0217839>). Dummy groups also require endpoint and
cycle checks because their partition-function contribution need not cancel
naively (Fleck et al., DOI
<https://doi.org/10.1021/acs.jctc.0c01328>). A soft-core coupling is required
to avoid insertion/deletion endpoint singularities (Beutler et al., DOI
<https://doi.org/10.1016/0009-2614(94)00397-1>).

SCIVAL verdict: **CONDITIONAL APPROVE**.

Conditions are the exact force scope, endpoint parity, separate correction
leg, state-overlap gates, and no automatic promotion defined below. An
uncorrected ghost is rejected.

## Thermodynamic Cycle

For environment `e` in `{bound, free}`, define:

- `u0`: physical stored-copy endpoint;
- `u1`: physical coordinate-swapped endpoint;
- `u1g`: the same `u1` endpoint plus the declared water-only ghost potential.

The main ghost-assisted route is

```text
DeltaG_main,e = G_e(u1g) - G_e(u0)
```

The correction removes the endpoint bias:

```text
DeltaG_off,e = G_e(u1) - G_e(u1g)
DeltaG_phys,e = DeltaG_main,e + DeltaG_off,e
```

The binding double difference is

```text
DeltaDeltaG_bind,phys = DeltaG_phys,bound - DeltaG_phys,free
```

The correction contribution must always be reported separately:

```text
DeltaDeltaG_bind,off = DeltaG_off,bound - DeltaG_off,free
```

No cancellation between bound and free may be assumed. For each ATM main or
direct route, the primary observable remains
`dgb = dgbind1 - dgbind2`. `dgbind1` is diagnostic only. The deghost leg is a
standalone endpoint-to-endpoint MBAR estimate, not a `dgbind1` half-leg. ATM
thermodynamic-cycle conventions follow the primary ATM and ATM-RBFE methods
(DOIs <https://doi.org/10.1021/acs.jctc.1c00266> and
<https://doi.org/10.1021/acs.jcim.1c01129>); MBAR follows Shirts and Chodera
(DOI <https://doi.org/10.1063/1.2978177>).

## Frozen Ghost Definition

### Interaction scope

The ghost acts only between:

- the nine W4A disappearing TRP heavy atoms
  `CG,CD1,CD2,NE1,CE2,CE3,CZ2,CZ3,CH2`; and
- every explicit-water oxygen selected by residue name `HOH,WAT,SOL` and atom
  name `O,OW,OH2`.

It does not act on water hydrogens, protein, peptide atoms outside the declared
ring, ions, cofactors, or other solvent species. It adds no charge, attraction,
bond, restraint, or water deletion. The disappearing atoms' stored coordinates
are used directly: no fixed point, translated point, virtual particle, or
ghost degree of freedom is added.

The implementation must add the ghost force after the canonical `ATMForce` is
built so that the ghost is a separate top-level force and is not evaluated
inside both ATM coordinate images.

### Force-field-derived parameters

No parameter is fit to the DCD gate or any FE result. Each atom's `sigma` and
`epsilon` are read from the nested canonical `NonbondedForce`, with
Lorentz-Berthelot mixing:

```text
sigma_ij   = (sigma_i + sigma_j) / 2
epsilon_ij = sqrt(epsilon_i * epsilon_j)
```

The frozen s127 reference readback is:

| pair class | sigma_ij (nm) | epsilon_ij (kJ/mol) | WCA cutoff (nm) |
|---|---:|---:|---:|
| aromatic C - water O | 0.328881703598677 | 0.484413968477885 | 0.369157230672847 |
| indole NE1 - water O | 0.321398154366298 | 0.681070223009209 | 0.360757230672847 |

These values were read from the nested canonical `NonbondedForce` in uncarved
W4A bound `s127` (`rep6`), serialized XML SHA256
`2cf60c082a7a22438fd0b795b9ab50d95b4c6dcad3058d2c852eff71b9255c55`.

Every pilot seed must reproduce the source per-particle parameters before a
context is created. Atom indices are topology-derived and must not be
hard-coded.

### Soft-core WCA expression

For coupling `g` in `[0,1]`, the pair potential is:

```text
rsc = (r^6 + alpha_sc*(1-g)*sigma_ij^6)^(1/6)
rc  = 2^(1/6)*sigma_ij
Vij = g*step(rc-rsc)*(
        4*epsilon_ij*((sigma_ij/rsc)^12-(sigma_ij/rsc)^6)
        + epsilon_ij
      )
Vghost = sum_ij Vij
```

Frozen constants and implementation policy:

- `alpha_sc = 0.5`;
- `g=0` is exactly zero energy and zero force;
- `g=1` is the exact force-field-derived WCA repulsive branch;
- `CustomNonbondedForce`, water/ring interaction group only;
- `CutoffPeriodic`, `cutoff = 0.40 nm`, no long-range correction;
- a dedicated force group and `dV/dg` recording are mandatory.

No amplitude multiplier, radius multiplier, contact-threshold fit, or
post-result parameter adjustment is permitted under this preregistration.

### Main-path coupling

The ghost follows a single continuous `u0 -> u1` progress coordinate while the
canonical ATM Hamiltonian remains unchanged:

```text
dplus:  g = Lambda2
dminus: g = 1 - Lambda2
```

Thus `g=0` at the dplus `u0` endpoint, `g=0.5` at both shared apexes, and
`g=1` at the dminus `u1` endpoint. `Lambda1`, the ATM soft-core constants, and
all H18 state values remain untouched. The two apex Hamiltonians must be
identical, including the ghost force.

## Pilot Cohort And Protocol

### P0: declaration and Hamiltonian audit

No MD and no GPU are allowed in P0. Use the OpenMM `Reference` platform on the
three frozen W4A source seeds:

```text
s101  control
s127  high leverage
s163  control
```

Audit both bound and free serialized systems. P0 must verify topology
selection, source parameter readback, force scope, endpoint parity, apex
parity, finite energy/forces, and serialization round-trip.

### P1: short structural mechanism pilot

- system: 2QKI W4A, `w4a_trp_ala_res4`, endpoint `wt`;
- source: uncarved inputs only;
- seeds: `s101,s127,s163`;
- leg: bound;
- directions: dplus and dminus;
- H18 15-state schedule unchanged;
- `n_cycles = 50`;
- `md_steps_per_cycle = 250`;
- `timestep_fs = 1.0`;
- `temperature_K = 300`;
- one DCD frame per cycle, with all explicit water oxygens available to the
  minimum-image analysis;
- no UWHAM or mutation FE verdict from P1.

The comparison is ghost-on versus the frozen no-ghost DCD evidence. No new
no-ghost transport-tuning run is launched in P1.

### P2: matched thermodynamic correction pilot

P2 is allowed only after P0 and P1 pass.

- seeds: `s101,s127,s163`;
- environments: bound and free;
- main route: ghost-assisted H18 dplus+dminus, `200` cycles, `250` MD
  steps/cycle, `mintimeid=100`;
- correction route: hold ATM at physical `u1`
  (`Direction=-1,Lambda1=0,Lambda2=0,Intermediate=0`) and sample `g` at
  `0,0.025,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.975,1.0`;
- correction sampling: Hamiltonian replica exchange, `200` cycles, `250` MD
  steps/cycle, one row per cycle, MBAR over the declared per-leg schedule;
- matched direct no-ghost route: the same H18 schedule and run length, used
  only as an empirical closure comparator;
- all routes use explicit per-seed RNG declarations and isolated output roots.

P2 is still a mechanism/integrity pilot. It does not authorize a W4A or W23A
scientific SIGN claim.

## Frozen Gates

### K0: declaration-to-artifact integrity

Any failure is `REJECT_GHOST`:

- exactly nine declared TRP ring atoms and all/only water oxygens selected;
- no static water deletion and no change in atom count or net charge;
- source sigma/epsilon match before context creation;
- ghost force is top-level, separately named, and assigned its own force group;
- exact protocol, source, code, and output hashes recorded per cell;
- no silent state schedule fallback.

### K1: Hamiltonian parity

Any failure is `REJECT_GHOST`:

- `(u0,g=0)` energy differs from unmodified `u0` by at most
  `1e-5 kJ/mol`, and maximum force-component difference is at most
  `1e-5 kJ/mol/nm` on Reference;
- correction `(u1,g=0)` matches unmodified physical `u1` to the same limits;
- main dminus `(u1,g=1)` matches correction `(u1,g=1)` to the same limits;
- dplus and dminus shared-apex energies at `g=0.5` agree to
  `1e-5 kJ/mol` on identical coordinates;
- ghost energy equals the dedicated force-group readback and is finite for all
  declared states;
- minimization with the ghost present does not alter the physical endpoint
  when `g=0` beyond the parity limits.

### S1: structural mechanism gate

All metrics use minimum-image distances and the ring-only primary selection.
P1 passes only if:

- control median `frac_lt_0p26 < 0.15` in bound/dminus;
- the high-leverage condition is false, where it is defined as both
  `s127-control_median >= 0.10` and `s127/control_median >= 1.5`;
- each cell has exactly `750` DCD frames and `50 x 11` finite rows per walker
  family, matching the frozen sidecar contract;
- ring/site heavy RMSD p95 is at most `0.15 nm`;
- no NaN, protein/ion force contact, gross site collapse, or new cohort
  heterogeneity appears.

The old all-site metric is retained as secondary compatibility output. It may
not replace the ring-only primary gate.

### T1: sampling and correction gate

For the ghost main and deghost correction routes:

- every adjacent pair has overlap `O >= 0.10`; any `O < 0.03` is a hard fail;
- both boundaries are visited and at least one full round trip occurs per seed
  and direction/route after burn-in;
- all state tables are finite and use the per-leg declared schedule;
- first-half versus second-half estimates differ by no more than
  `max(1.0 kcal/mol, 2*combined_SE)`;
- correction bootstrap SE is at most `1.5 kcal/mol` per environment;
- corrected total uncertainty must not exceed twice the main-route uncertainty.

Correction magnitude is never a pass/fail criterion; it must be reported even
when large. Failure of overlap, round trips, or uncertainty is
`THERMO_PILOT_INDETERMINATE`, not evidence that the correction is zero.

### T2: route closure

When the matched direct route passes its own sampling gate, compute per seed:

```text
R_e = (DeltaG_main,e + DeltaG_off,e) - DeltaG_direct,e
```

Closure supports the design only if the paired median residual is within
`1.0 kcal/mol` and each environment's residual is statistically compatible
with zero at two combined standard errors. A resolvable nonzero residual is
`REJECT_GHOST`.

If the direct route fails its already known transport gate, closure is
`UNKNOWN`, not zero and not a ghost failure. In that case no `n=10` promotion
is automatic; a new SciVal decision is required using K0/K1, S1, and T1.

## Preregistered Outcomes

1. `A_PILOT_PASS`: K0/K1, S1, T1, and interpretable T2 all pass. Prepare a new,
   separate `n=10` preregistration; do not launch automatically.
2. `A_CORRECTION_VALID_DIRECT_UNKNOWN`: K0/K1, S1, and T1 pass, but the direct
   route is non-convergent. Return to SciVal; no automatic production.
3. `A_PILOT_INDETERMINATE`: integrity passes but structural or sampling power is
   insufficient. Freeze outputs and preregister any revision before rerun.
4. `REVISE_GHOST`: structural gate fails without Hamiltonian/integrity failure.
   Parameter fitting to the observed result is forbidden; a new physical
   design and preregistration are required.
5. `REJECT_GHOST`: force scope, endpoint parity, charge/atom integrity, or
   resolvable cycle closure fails.

## Stop Conditions

Stop immediately on any R-rule violation, sign flip against an anchor, sigma
doubling, tier boundary crossing, new cohort heterogeneity, NaN, source hash
drift, endpoint/apex parity failure, non-water ghost interaction, schedule
fallback, output collision, or attempted git push.

No failure under this preregistration authorizes fallback to B tuning,
forward-only `dgbind1`, `base=u1` capping, or an uncorrected cavity bias.

## Implementation Boundary

This document freezes the scientific design only. No ghost code, GPU pilot,
W23A launch, commit, or push is part of this stage. A future implementation
must pass the six UPDD roles in order and must generate a fresh source/code
inventory before P0.
