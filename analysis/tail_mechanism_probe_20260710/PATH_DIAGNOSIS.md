# Tail-Mechanism Probe PATH Diagnosis

Date: 2026-07-10

## Runner Integrity

- Command: `/home/san/miniconda3/envs/atm/bin/python analysis/tail_mechanism_probe_20260710/tail_probe.py`
- `py_compile`: PASS
- Parsed cells: 68 total
  - `w23a_gateA`: 20 cells
  - `w4a_c1_carved`: 24 cells
  - `w4a_c1_uncarved`: 24 cells
- Parsed state rows: 876 total
- Malformed `.out` rows: 0
- Generated artifacts:
  - `tail_probe_report.md`
  - `tail_probe_summary.json`
  - `seed_summary.csv`
  - `cell_summary.csv`
  - `state_tail_summary.csv`

## PATH Read

### W23A Gate-A

- Highlight seed `s101` is also the top absolute ddG seed: `ddG=-16.005`.
- `s101` is the top seed by max global `pertE` p99: `190.829`, only `+1.221` above the non-highlight median.
- `s101` is not enriched by `frac(pertE>150)` relative to controls: `-0.008` vs control median.
- `s101` has the lowest adjacent crossing count in the cohort: `150`, `-98` below control median, but this is not a zero-crossing wall and no gate fail is present.

Interpretation: W23A `s101` has a real scalar-log tail signal plus weaker mixing than controls, but not a catastrophic single-boundary failure. This supports targeted densify/overlap repair around the endpoint/apex cells, not a conclusion that carve is globally invalid.

### W4A C1 Carved

- Highlight seed `s127` is the top absolute ddG seed: `ddG=-12.563`.
- `s127` is the top seed by max global `pertE` p99: `191.301`, `+2.681` above the non-highlight median.
- `frac(pertE>150)` is not meaningfully enriched: `+0.001` above controls.
- `s127` is not the weakest adjacent-crossing seed; its minimum crossing is `42`, while the cohort minimum is `14` at `s23`.
- Gate failures are cohort-wide in the W4A short C1 runs, not `s127`-specific.
- The seed-level effect is leg-asymmetric: bound `dgb=-6.667`, free `dgb=+5.895`.

Interpretation: W4A carved `s127` is not explained cleanly by a simple lambda underlap or single wall. The stronger read is a leg-asymmetric tail mechanism that is visible in `pertE` p99 but not in high-threshold fraction or crossing minima. This keeps C1 INDETERMINATE and argues for a short structural/DCD-on probe before treating broad densify as sufficient.

### W4A C1 Uncarved Reference

- The largest absolute uncarved seed is `s199`: `ddG=-4.408`.
- `s199` is also the top uncarved seed by max p99: `191.024`, with top `frac(pertE>150)=0.139`.
- Direct paired comparison seed `s127` has `ddG=+3.655` and p99 `189.510`, not the top uncarved outlier.

Interpretation: scalar-log tail sensitivity exists even without carve, but the carved `s127` magnitude and sign pattern are different. This does not prove carve bias; it does say the C1 carve question cannot be settled by point-threshold or scalar `.out` statistics alone.

## Decision

Best next step:

1. Run a small `D2` structural/DCD-on probe on the exact high-leverage cells before full `B`.
2. For W23A, include `s101` `bound/dplus`, `bound/dminus`, and `free/dplus` with one matched control.
3. For W4A C1 carved, include `s127` `free/dplus` and `bound/dminus`, plus one matched non-outlier carved control.
4. If D2 shows water/void penetration or a persistent structural basin, promote carve redesign `A` for that mechanism.
5. If D2 is structurally clean and the only pathology is endpoint/apex overlap, run targeted `B-local` densify plus seed expansion rather than broad blind densify.

Rejected as next immediate move:

- `C` forward-only: still risks reintroducing one-leg offset/apex bias and does not explain the two-arm C1 behavior.
- Full `A` immediately: current scalar evidence is suggestive but not sufficient to prove carve-root pathology.
- Full broad `B` immediately: useful eventually, but D2 can focus where densify should be applied.
