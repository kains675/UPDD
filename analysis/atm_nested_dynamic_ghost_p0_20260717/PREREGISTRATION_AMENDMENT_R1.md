# ATM-Nested Dynamic Ghost P0 Preregistration Amendment R1

Registered: 2026-07-17 KST, after development smoke and before official P0
inventory, KEEPER output, or RUNNER output

Status: `FROZEN_AFTER_DEV_SMOKE_BASELINE_DISCOVERY_BEFORE_OFFICIAL_OUTPUT`

SciVal verdict: `CONDITIONAL APPROVE`

## Preserved R0

The original `PREREGISTRATION.md`, `SCIVAL_VERDICT.md`, `protocol.json`, and
`FROZEN_MANIFEST.json` remain unchanged. R0 correctly froze the Hamiltonian,
cohort, platform, probe distance, and numerical thresholds.

## Development-Smoke Finding

An unrecorded, non-official `w4a_union_s101_free` development smoke found:

- exact `g=0` source endpoint energy and force parity;
- intended stored/transformed probe activation;
- a pre-existing source `u0` ghost contribution of about `0.501 kJ/mol`;
- no corresponding source `u1` ghost contribution.

The R0 absolute opposite-endpoint metric therefore includes the unchanged
source ghost background. It cannot distinguish probe leakage from a legitimate
source contact.

## Frozen R1 Correction

For each endpoint, define the probe-specific ghost response as:

`[probe(g=1)-probe(g=0)] - [source(g=1)-source(g=0)]`

Apply the existing intended-energy, opposite-energy, radial-force,
opposite-force, and force-balance thresholds to this endpoint-matched
double-difference. Preserve the absolute `g=1-g=0` values as diagnostics.

The deterministic probe-water candidate set is further restricted to waters
whose source oxygen is at least `0.42 nm` from every stored and transformed
ring atom under minimum-image distance. Among eligible three-atom,
zero-displacement waters, retain the R0 maximum H-centroid/outward-dot rule and
lowest-oxygen-index tie break.

## Scientific Rationale

This is a measurement correction, not a relaxed threshold. It removes a
coordinate-independent source background while preserving the same force,
energy, endpoint-selectivity, serialization, and six-cell acceptance gates.
It makes no free-energy or delta-G-neutrality claim and does not authorize
sampling.

