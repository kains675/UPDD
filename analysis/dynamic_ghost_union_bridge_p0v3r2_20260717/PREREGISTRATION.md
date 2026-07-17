# Dynamic-Ghost Union Apex Bridge P0v3r2 Preregistration

Registered: 2026-07-17 KST, before P0v3r2 implementation or output

Status: `FROZEN_BEFORE_P0V3R2_IMPLEMENTATION_OR_OUTPUT`

Regime: R-11 ranking/SIGN-only mechanism validation. No calibrated absolute
Delta G or Delta Delta G claim.

SciVal verdict: `CONDITIONAL_APPROVE_P0V3R2_ONLY`.

## Trigger

P0v3r1 formally stopped at `w4a_union_s101_bound` with
`REJECT_UNION_SOURCE`. Every count, charge, placeholder, box, separation, and
stored/transformed geometry gate passed. The sole false check compared Python's
global RNG before and after the entire two-copy build, although the frozen
protocol requires save/restore only around the cell-seeded `addSolvent` call.

The changed RNG state also exposed an unregistered reproducibility source:
relative to the frozen parent, nine residue-4 hydrogen coordinates changed,
with ALA4 methyl-H displacement as large as `0.21780597 nm`. P0v3r2 therefore
does not merely remove the over-broad check. It makes the existing
appearing-hydrogen placement deterministic before rebuilding every source.

P0v3r1 remains rejected and its files must not be edited, reclassified, or used
as P0v3r2 source artifacts.

## Inherited protocol

P0v3r2 inherits every P0v3r1/P0v3 clause for:

- W4A `s101,s127,s163 x bound/free` cohort and parent counts;
- canonical two-copy W4A construction and 4.0 nm displacement;
- temporary nine-heavy-atom union placeholders and exact `numAdded` solvation;
- placeholder removal, exact counts, charge, atom ordering, and box gates;
- stored/transformed ring-water minimum `0.26 nm` and zero sub-threshold pairs;
- report-only static total-box `N/V` and prohibition on equilibrium-density claims;
- unchanged dynamic ghost at `g=0.5` and explicit apex bridge;
- raw endpoint, finite-value, bridge identity, serialization, CPU Reference,
  subprocess isolation, and first-failure stop gates;
- zero MD, minimization, GPU, trajectory, FE estimation, or automatic next-stage
  launch.

No threshold, molecule count, placeholder parameter, bridge equation, or
Hamiltonian term changes in P0v3r2.

## Deterministic appearing-H construction

Resolve the mutation as the registered `MutationSpec` object
`w4a_trp_ala_res4`, then use the existing
`build_inplace_res4_twocopy_system_r2_retry` path with `retry_k=5`.

The unit key is exactly:

```text
{seed}|{leg}|w4a_trp_ala_res4
```

Attempt `i` uses the existing stable seed derivation:

```text
int(SHA256("{unit_key}#{i}")[0:8], 16) & 0x7fffffff
```

Frozen attempt seed trails:

| cell | attempts 0..4 |
|---|---|
| s101 bound | 1248013600, 821365542, 13400133, 123974878, 1391514182 |
| s127 bound | 930194524, 198058456, 1979308714, 1695657795, 643592353 |
| s163 bound | 652776389, 2100013619, 514827842, 148968126, 1852907068 |
| s101 free | 1570880720, 1238137364, 280468345, 576613175, 1311681870 |
| s127 free | 1081154822, 301540652, 1620431659, 1668652495, 550394267 |
| s163 free | 1103571852, 1876770155, 627858153, 1957885304, 2057776391 |

The first R2-valid attempt is accepted. Only the existing R2 seed-clash error
may retry; all other errors fail loud. The report must record unit key, retry
limit, attempts used, full seed trail, and accepted seed, and these must equal
the frozen derivation. No seed is selected from P0v3r1 energy or geometry output.

## RNG scope correction

The cell-specific frozen union-solvent RNG checks remain hard gates:

1. the declared solvent seed is recorded;
2. `addSolvent` consumes that seeded stream;
3. Python's global RNG state immediately before `addSolvent` is restored
   immediately afterward.

The global RNG state before versus after the entire builder is not a gate. The
entire source runs in its own subprocess, and upstream appearing-H randomness is
now replaced by the deterministic seed trail above.

## Execution and outcomes

All six P0v3r2 sources are rebuilt fresh. P0v3r1 XML/PDB files are evidence only
and must not be copied into the new source lane.

Allowed outcomes:

1. `P0V3R2_PASS`: all six source, KEEPER, raw, and bridge gates pass;
2. `REJECT_UNION_SOURCE`: any source/count/charge/placeholder/RNG/geometry/raw
   endpoint gate fails;
3. `REJECT_BRIDGE`: source gates pass but bridge identity fails;
4. `P0V3R2_INTERRUPTED`: process/integrity interruption, not a scientific pass or
   fail.

Even `P0V3R2_PASS` validates only construction and Hamiltonian identities at
fixed coordinates. It does not establish equilibrium density, Delta-G neutrality,
or authorize sampled production.

