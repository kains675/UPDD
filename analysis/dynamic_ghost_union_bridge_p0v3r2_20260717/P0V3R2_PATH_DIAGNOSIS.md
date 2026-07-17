# P0v3r2 PATH Diagnosis

Date: 2026-07-17 KST

## Formal outcome

`P0V3R2_PASS`: all six W4A union-solvated source cells passed the frozen
CPU-only OpenMM Reference audit. No source or bridge cell was rejected.

- inventory digest:
  `0457824c66c54961d2d65fe110817b61b32e828285be5a7e2ff15fb55ba4f702`
- source summary digest:
  `d7257f100a2e90657d8f2590bf3ec461060e8822c6ece21032ee44f7d5875e08`
- platform: `Reference`
- MD steps: `0`
- minimization steps: `0`
- GPU used: `false`
- free-energy estimate produced: `false`

PATH verdict: `PASS_P0_CONSTRUCTION_ONLY`.

## Source pathology diagnosis

The count-preserving union construction removed the initial appearing-image
water clash that stopped P0v2. Every cell retained its frozen water/ion/particle
counts and parent net charge. All nine temporary placeholders matched their
ATM-transformed targets within `1.34e-15 nm` and were absent from the final
canonical system. Every deterministic appearing-H build accepted frozen attempt
1.

Across both stored and transformed images, the minimum ring/water-O distance
was `0.36778970 nm`; every cell had zero pairs below the frozen `0.26 nm`
threshold. The six source minima ranged from `0.36778970` to `0.41746486 nm`.

The critical comparison is `s163/bound`:

| metric | rejected uncarved P0v2 | union P0v3r2 |
|---|---:|---:|
| transformed ring/water-O minimum (nm) | 0.05695665 | 0.41746486 |
| raw u1 (kJ/mol) | 2659781997.345677 | -252777.525685 |
| apex force-component gap (kJ/mol/nm) | 276209950501.09717 | 0.0 |

All P0v3r2 raw endpoint magnitudes were below `3.60e6 kJ/mol`, and the largest
source force component was `2.9468e6 kJ/mol/nm`, below the frozen `1e8` gross
clash guard. The prior billion-kJ source pathology did not recur.

## Bridge and ghost diagnosis

| cell | raw u1-u0 (kJ/mol) | Hminus-Hplus (kJ/mol) | max force gap (kJ/mol/nm) |
|---|---:|---:|---:|
| s101 bound | 1231.248885 | 357.661025 | 17747.622553 |
| s127 bound | 3387.309398 | 1420.142669 | 90799.177087 |
| s163 bound | 117.975871 | 0.000000 | 0.000000 |
| s101 free | 2082.524535 | 774.636259 | 43275.863089 |
| s127 free | 2839.460955 | 1148.631021 | 70479.908060 |
| s163 free | 2735.544356 | 1097.193636 | 57734.422349 |

For all six cells:

- bridge endpoint energy and force parity error was exactly zero;
- bridge endpoint-gap identity error was exactly zero;
- interior energy and derivative identity error was exactly zero;
- maximum interior force identity error was
  `9.313225746154785e-10 kJ/mol/nm`, versus the frozen `1e-5` tolerance;
- XML round-trip energy and force error was exactly zero;
- ghost energy and `dV/dg` were invariant across bridge xi;
- inactive legacy ATM globals had exactly zero energy and force effect.

The nonzero bridge endpoint gaps are required work-landscape diagnostics, not
free-energy estimates and not failures. A fixed-coordinate gap cannot be used
as `DeltaG_bridge`, assumed to cancel between bound and free, or converted into
a mutation verdict.

## Residual unknowns

P0v3r2 changed only source construction and audited fixed coordinates. It did
not establish any of the following:

1. equilibrium NPT density or stable post-relaxation box volume;
2. stable MD after minimizing the unrelaxed source force field;
3. dynamic water re-entry behavior under union initialization plus ghost;
4. adjacent bridge-state overlap, boundary visits, or round trips;
5. bridge MBAR convergence, uncertainty, or bound-minus-free correction;
6. Delta-G neutrality, calibrated absolute Delta G, or quantitative Delta Delta G.

The static total-box solvent-count density diagnostic ranged from `-3.4384%`
to `-1.1410%` relative to the parent boxes. It remains report-only because it
does not account for solute excluded volume and is not an equilibrium
measurement. The largest fixed-coordinate source force component remains
large enough that minimization and a short stability gate are mandatory before
production sampling.

## Next authorized design step

Preregister a separate sampled-pilot revision before creating any CUDA Context.
The pilot must separate schedule discovery from a free-energy claim:

1. run count/charge/hash verification again, minimize each selected xi source,
   and perform NPT relaxation with density/volume stationarity and finite-force
   gates;
2. use the exact frozen bridge Hamiltonian with ghost `g=0.5` and collect a
   short, explicitly non-production multi-xi energy matrix;
3. apply a preregistered deterministic midpoint-insertion rule to intervals
   with inadequate adjacent overlap; pilot data may select a later schedule but
   may not contribute to its production MBAR estimate;
4. require all six seed/leg cells in schedule discovery, because the largest
   bridge gap, largest force, worst static-density diagnostic, and prior clash
   occur in different cells;
5. run each window/cell in an isolated process so a single 16-GB GPU holds at
   most one large bound Context and OpenMM memory is released at process exit;
6. freeze fresh RNG, equilibration, sampling, energy-matrix, overlap,
   autocorrelation, and stop gates before launch;
7. do not authorize W23A, n=10 expansion, or a mutation SIGN verdict from this
   P0 pass.

P0v3r2 itself is complete. Sampled work remains blocked until the new
preregistration, SciVal, and KEEPER gates pass.

