# Union-Solvated Apex-Bridge P1SD Preregistration

Registered: 2026-07-17 KST, before P1SD implementation or output

Status: `FROZEN_BEFORE_P1SD_IMPLEMENTATION_OR_OUTPUT`

Regime: Track B mechanism and schedule discovery. Ranking/SIGN-only; no
calibrated absolute Delta G or Delta Delta G claim.

Scientific review verdict: `CONDITIONAL_APPROVE_P1SD_ONLY`.

## Trigger and correction to the P0 PATH proposal

P0v3r2 passed all six count-preserving union source and exact apex-bridge
Reference cells. It did not sample the bridge or establish overlap. Its static
union boxes are 1.14-3.44% larger by volume than the corresponding frozen parent
boxes.

The P0 PATH note proposed NPT relaxation. Existing Track B production evidence
shows that an OpenMM `MonteCarloBarostat` is unsafe with the displaced-particle
`ATMForce`: trial volume moves evaluate the shadow `u1` configuration and have
previously produced NaN. P1SD therefore does not add a barostat to an ATM
Context. It restores the exact parent NPT box deterministically, then uses the
same fixed-volume NVT regime as the validated Track B production path.

This is an R-18 correction of the next-step implementation, not a reclassification
of the completed P0v3r2 result.

## Scientific purpose

P1SD answers only:

1. whether each parent-volume-normalized union source remains finite and
   structurally safe after staged minimization and short NVT propagation;
2. whether the initial five-state bridge grid has adequate adjacent phase-space
   overlap in every seed and environment;
3. which intervals require deterministic midpoint insertion before a fresh
   production protocol is frozen.

P1SD does not estimate a production bridge correction. Its exploratory MBAR
free energies and all pilot frames are schedule-design evidence only and must
not be reused in a later production estimator. MBAR follows Shirts and Chodera,
J. Chem. Phys. 129, 124105 (2008), DOI `10.1063/1.2978177`.

## Frozen cohort

Use the completed P0v3r2 union artifacts without rebuilding or editing them:

```text
s101 bound
s127 bound
s163 bound
s101 free
s127 free
s163 free
```

All six are required because the largest bridge gap, largest fixed-coordinate
force, worst box-volume delta, and former source clash occur in different cells.

## Parent-box normalization

For each cell, use the exact target box vectors from the P0v3r2 frozen parent
contract. Do not scale solute atoms or intramolecular solvent coordinates.

For each water or monatomic ion residue:

1. use the water oxygen or ion atom as the anchor;
2. convert the anchor to wrapped fractional coordinates in the union box;
3. place the anchor at the same fractional coordinates in the parent box;
4. translate every atom in that residue by the anchor displacement.

The final topology and System box vectors must equal the parent vectors. The
following are hard gates before any CUDA Context:

- solute positions unchanged exactly in float64 memory;
- solvent intramolecular displacement vectors unchanged within `1e-12 nm`;
- particle, topology, water, Na, Cl, and net-charge contracts unchanged;
- box-vector maximum absolute error at most `1e-10 nm`;
- stored and transformed ring/water-O minima at least `0.26 nm`, with zero
  sub-threshold pairs;
- all-solute-heavy/water-O minimum at least `0.24 nm`;
- source and bridge XML force contracts unchanged except for default periodic
  box vectors;
- authoritative positions stored as float64 NumPy data; PDB is inspection-only.

The parent box is an inherited fixed-volume ensemble declaration, not a new
equilibrium-density measurement.

## CPU Reference preflight

Run all six normalized cells in isolated subprocesses before CUDA sampling.
At `xi=0,0.5,1`, require:

- finite potential energy, forces, and `dV/dBridgeXi`;
- `abs(raw u0)` and `abs(raw u1)` below `1e8 kJ/mol`;
- maximum absolute force component below `1e8 kJ/mol/nm`;
- exact endpoint, bridge-linear, derivative, ghost, inactive-global, and XML
  round-trip contracts inherited from P0v3r2;
- no MD, minimization, GPU, or free-energy estimate.

Any failure is `REJECT_PARENT_BOX_SOURCE` or `REJECT_BRIDGE` and stops P1SD.

## Frozen CUDA sampling protocol

Initial bridge grid:

```text
BridgeXi = 0.00, 0.25, 0.50, 0.75, 1.00
ghost g = 0.50 at every state
```

Execution contract:

- OpenMM platform `CUDA`, device `0`, precision `mixed`;
- one cell/window subprocess and at most one live CUDA Context;
- no barostat; exact parent box remains fixed;
- unconstrained P0v3r2 System, `1.0 fs` timestep;
- `LangevinMiddleIntegrator`, `300 K`, friction `1/ps`;
- common per-cell source preparation at physical `u0`: LocalEnergyMinimizer
  tolerance `10 kJ/mol/nm`, maximum `5000` iterations;
- each bridge window starts from the common minimized positions, then minimizes
  at its own xi with the same tolerance and maximum `1000` iterations;
- independent integrator and velocity seeds derived from the frozen SHA256 rule;
- warmup `2000` steps (`2 ps`), excluded from every matrix and trajectory;
- discovery sampling `50` cycles x `250` steps = `12.5 ps` per window;
- one DCD frame and one energy/derivative row per cycle (`50` each);
- no replica exchange and therefore no round-trip claim in P1SD.

The seed unit key is exactly:

```text
P1SD|{cell_id}|xi={xi:.2f}|{role}
```

where role is `integrator` or `velocity`. The seed is:

```text
int(SHA256(unit_key)[0:8], 16) & 0x7fffffff
```

Zero is replaced by one. Every consumed seed and unit key must be recorded.

## Exact cross-state matrix

The bridge is linear in xi. For a sample `x` generated at xi_i:

```text
U(x;xi_j) = U(x;xi_i) + (xi_j-xi_i) * dV/dBridgeXi(x)
u_j(x) = beta * U(x;xi_j)
```

This identity constructs the complete cross-state reduced-potential matrix
without additional Contexts. P1SD must verify the identity by explicit Context
readback at all five xi values for the first and final sample of every window,
with absolute energy error at most `1e-5 kJ/mol`.

Use PyMBAR 4.2 robust solver settings: maximum iterations `20000`, relative
tolerance `1e-8`, and `solver_protocol="robust"`. The overlap matrix must be
finite, every diagonal at most `1.001`, and every row sum within `1e-3` of one.

## Structural and execution gates

Every window must satisfy:

- exactly `50` finite sample rows and `50` DCD frames;
- final integrator step `14500` including warmup;
- fixed parent box throughout;
- mean sampled temperature in `[270,330] K`, and no sample above `600 K`;
- no stored or transformed ring/water-O distance below `0.10 nm`;
- finite final energy and forces, with maximum force component below `1e8`;
- XML, positions, DCD, table, and result hashes recorded;
- worker exit releases the CUDA Context before the next window starts.

The fraction of frames below `0.26 nm` is reported separately for stored and
transformed images. It is a mechanism diagnostic, not an overlap surrogate.

## Frozen overlap decision

Use the existing UPDD convention `O_i,i+1` from the PyMBAR overlap matrix.
Reduce every interval by the minimum over all six cells.

1. `GRID_CANDIDATE`: every cell and adjacent interval has `O >= 0.10`.
2. `DENSIFY_REQUIRED`: at least one interval has `0.03 <= O < 0.10`.
3. `SEVERE_BOTTLENECK`: at least one interval has `O < 0.03`.
4. `P1SD_INDETERMINATE`: MBAR unavailable/non-physical, incomplete output,
   finite-value failure, or a structural/execution gate failure.

For outcomes 2 or 3, the proposed next grid is the union of the current grid
and the exact midpoint of every interval whose six-cell minimum is below
`0.10`. P1SD does not launch that grid automatically. Outcome 1 is only a
schedule candidate; it also requires a fresh production preregistration and
fresh samples.

## Stop conditions and prohibited inferences

Stop on source/code hash drift, count/charge mismatch, geometry failure,
non-finite value, CUDA OOM/error, temperature gate failure, output collision,
first failed window, or attempted git push.

P1SD cannot authorize W23A, n=10 expansion, B transport tuning, forward-only
`dgbind1`, `base=u1` capping, a mutation SIGN verdict, or a quantitative bridge
correction. Pilot MBAR magnitude is never a pass/fail threshold.

