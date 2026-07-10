# Tail-Mechanism Probe Preregistration

Date: 2026-07-10

## Purpose

This is a read-only, 0-GPU pathology probe before launching broader seed expansion or
new carved protocols. It uses existing Track B artifacts only.

## Cohorts

- `w23a_gateA`: `outputs/_trackb/mdm2_w23a_gateA_20260708`
- `w4a_c1_carved`: `outputs/_trackb/w4a_carved_c1_20260709/*_carved`
- `w4a_c1_uncarved`: `outputs/_trackb/twocopy_w4a_*_FIXAB`

## Highlight Seeds

- W23A: `s101`, because it is the largest negative paired ddG outlier.
- W4A carved: `s127`, because it drives the carved C1 heavy-tail concern.
- W4A uncarved: `s127` for direct paired comparison and the largest absolute
  uncarved seed by the frozen reference JSON for baseline context.

## Metrics

For each cohort, seed, leg, and direction:

- `pertE` global percentiles and fractions above 120, 150, and 180 kcal/mol.
- per-state `pertE` concentration, using the state with maximum p99 as the
  primary tail-localization marker.
- manifest-declared adjacent crossing minimum, round-trips, gate pass, and
  state count.

## Decision Use

- If an outlier seed is localized to one leg/direction/state and also has weak
  adjacent crossing, prefer targeted densify around that lambda region.
- If an outlier seed has no `pertE` or mixing localization in existing `.out`
  files, treat the cause as not visible to scalar energy logs and prefer a short
  DCD-on structural probe before full seed expansion.
- If carved and uncarved behavior diverge specifically at the carved highlight
  seed, treat carve-mediated dynamics as plausible but not proven.

## Non-Goals

- No new calibrated ddG claim.
- No proof of carve neutrality.
- No production-engine change.
