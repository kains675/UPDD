# W4A carved-vs-uncarved — review condition C1 verdict (void-water carve dG-neutrality)

_Auto-generated 2026-07-10 01:07:40 by `scripts/w4a_c1_result.py`. Regime: R-11 SIGN/order-only, R-18 honesty._

## TL;DR

| quantity | value |
|---|---|
| **C1 VERDICT** | **INDETERMINATE** |
| **paired delta (carved − uncarved), same 6 seeds** | **-4.07 kcal**, paired \|t\|=1.55 (df=5), n=6 |
| **robust delta (drop s127 from BOTH arms)** | **-1.64 kcal**, n=5 (collapses toward noise) |
| carved-vs-uncarved per-seed correlation r | -0.43 (anti-correlated → heavy-tail artifact, not a coherent bias) |
| unpaired delta / combined SEM (secondary) | -4.07 / 2.23 → z=1.83 |
| ddG_bind CARVED (dgb, Pairing B) | -4.39 kcal (\|t\|=2.30, n=6, sign=undetermined) |
| ddG_bind UNCARVED ref (dgb, Pairing B) | -0.32 kcal (\|t\|=0.28, n=6, sign=undetermined) |
| anchor (W23A destabilizing, reporting-only) | ~4.5 kcal |
| W23A carved result (the effect C1 gates) | dgb ddG_bind = -9.81 kcal, \|t\|=6.25 (LARGE = not-flatten) |

**VERDICT — INDETERMINATE.** carve dG-neutrality cannot be resolved at this sampling: the paired same-seed test is not significant (|t|=1.55 < 3), and the robust (drop-s127) delta collapses to -1.64 kcal (40% of the point estimate), and the CI (-10.80..2.66 kcal) spans zero. The carve is NEITHER proven neutral (the CI does not exclude a several-kcal bias) NOR proven non-neutral (fails significance / collapses on outlier removal); both directions are underpowered.

## What C1 tests

W23A carved in-place two-copy RBFE returned a LARGE binding effect (dgb ddG_bind = -9.81 kcal, |t|=6.25) → the construction does NOT flatten a known large binding effect. But the void-water carve fix (needed to build W23A's crash-prone free leg) has UNPROVEN dG-neutrality. C1 runs W4A CARVED at the exact config of the existing UNCARVED W4A. Both runs share the SAME 6 seeds (s7,s23,s101,s127,s163,s199), so the carve's dG-perturbation is isolated by a **PAIRED same-seed delta** d[s] = carved_ddg[s] − uncarved_ddg[s] — NOT a difference of two means against a point-threshold.

## Paired same-seed statistics (the correct, more powerful test)

| seed | carved ddG | uncarved ddG | delta d[s] |
|---|---|---|---|
| s7 | -5.28 | -0.91 | -4.37 |
| s23 | -3.27 | -2.14 | -1.13 |
| s101 | -3.80 | +0.68 | -4.48 |
| s127 | -12.56 | +3.66 | -16.22  ← max-\|z\| outlier |
| s163 | +1.86 | +1.21 | +0.65 |
| s199 | -3.29 | -4.41 | +1.12 |

- **PAIRED mean delta = -4.07 kcal, |t| = 1.55 (df=5), n=6 → NOT significant** (|t| < 3). The CI -10.80..2.66 kcal spans zero.
- **ROBUST (drop s127, the max-|z| seed, from BOTH arms): delta = -1.64 kcal (n=5)** → the point estimate collapses toward noise once the single heavy-tail seed is removed.
- **carved-vs-uncarved per-seed correlation r = -0.43 (anti-correlated).** A systematic ~4.1 kcal carve bias would shift seeds ~coherently; instead s127 swings 16.2 kcal same-seed (carved -12.56 vs uncarved 3.66) while s23/s163 barely move → the signature of a heavy-tailed **sampling artifact**, not a reproducible carve bias.
- The paired test being LESS significant than the unpaired one (|t| 1.55 < z 1.83) confirms the two arms do NOT track a common carve signal.
- Secondary (dgbind1, forward-only, non-offset-cancelled): carved ddG_bind = -3.35 kcal (|t|=1.82), uncarved = -0.20 kcal (|t|=0.20). Reported for completeness; the primary offset-cancelled dgb is authoritative.

## Why INDETERMINATE, not the earlier "FAIL"

The script previously flagged "CARVE NON-NEUTRAL (C1 FAIL)" from a |delta|>=4.0 kcal **point-threshold** on the difference of means. That is a false-positive: it ignores the SEM (combined ~2.23 kcal) and the single-seed heavy tail (s127). Three of four proper analyses (paired |t|=1.55, robust delta=-1.64, anti-correlation r=-0.43) fail to establish a ≥4 kcal bias. Per scientific review (the formal C1 review) the correct status is INDETERMINATE: the carve is neither proven neutral nor proven non-neutral at n=6/400cyc.

## Does the carve, if real, FLATTEN or INFLATE? → INFLATE

Carved is MORE negative than uncarved (delta -4.07, SAME sign as W23A's -9.81). W23A is a Trp→Ala abolishing the deepest p53 anchor → a large-negative ddG_bind IS the signal. The carve pushes in the SAME direction, so it can only **inflate** the magnitude, never flatten it. Worst-case correction of W23A (subtract the carve delta):

| W23A | raw | minus full carve (-4.07) | minus robust carve (-1.64) |
|---|---|---|---|
| full (n=5) | -9.81 | -5.74 | -8.17 |
| ex-s101 (n=4) | -8.26 | -4.19 | -6.62 |

In every scenario W23A stays large (|ddG| ≥ ~4.2) and same-sign. The only way the carve could flatten W23A is by adding > ~8 kcal of spurious binding (~3.6σ beyond the C1 point estimate, and in the inflating not flattening direction) — not supported. So the carve **INFLATES, it does not FLATTEN**; the W23A SIGN/ORDER "large, non-flat binding effect" conclusion stands (bound-specific: free ~0 vs bound -9.8), with caveats. NOT a quantitative claim.

## Combined conclusion (funnel implication)

- **C1 = INDETERMINATE.** At n=6/400cyc the carve's effect on ddG_bind is statistically indistinguishable from zero (paired |t|=1.55, robust collapses to -1.64, anti-correlated r=-0.43). It is NOT proven neutral (CI ~-10.8..2.7 spans a several-kcal bias) and NOT proven non-neutral (fails significance). Both underpowered.
- **W23A SIGN/ORDER conclusion STANDS with caveats**: the carve inflates not flattens, so W23A's large-negative, bound-specific (free ~0 vs bound -9.8) effect survives the worst-case carve correction (stays |ddG| ≥ ~4.2, same sign). This is a SIGN/order methods finding, NOT a certified/quantitative ddG.
- This is NOT "confounded / do not publish". It IS "carve bias bounded but not proven neutral; the W23A SIGN conclusion is robust to the worst-case bias because the carve inflates, not flattens".

**OPEN DECISION (scientific-review-ranked):**
1. **(B) Densify + expand seeds** on W23A AND the W4A carved/uncarved control, with the corrected paired / robust / bootstrap statistics. Minimum-assumption, power-acquiring step; decides whether the s127/s101 tail is sampling (shrinks with n → carve salvageable) or structural (persists → A mandatory). **Top pick.**
2. **(A) Carve redesign → union-solvation / ghost-atom placeholder** (dG-faithful). The correct destination: physically occupies the displaced-copy void throughout dynamics (fixes the residual tail the t=0-only carve leaves) and dissolves the neutrality question. Becomes 1st if B shows the tail is structural.
3. **(C) Forward-only / base=u0** as a cross-check only — dropping the backward leg reintroduces the box-size soft-core apex OFFSET that cancels only in the bidirectional dgb (barnase Gate-A re-adjudication 2026-07-07). Capping base=u1 is rejected (u1 is a physical endpoint; capping masks a real clash).

## NEW RISKS (scientific review, not in the earlier summary)

1. **Irreducible transferability gap** — C1 can NEVER certify the *W23A* carve directly: uncarved W23A does not build (NaN crash — the reason the carve exists), so C1 must proxy via W4A, a DIFFERENT water set. W23A carve neutrality is structurally unprovable by the carve-vs-uncarve method.
2. **Leg-asymmetric carve does NOT cancel in the double difference** — W23A carve is bound 0 / free 5 waters (enters ddG_bind through the free leg ONLY); W4A carve is bound 3 / free 4 (near-symmetric → largely cancels). So even a conclusive C1 on W4A would not transfer to W23A.
3. **Collapsed near-zero baseline** — the W4A uncarved reference is itself O3-collapsed (min-O ~0.065, ddG -0.32, |t|=0.28). Differencing two heavy-tailed ~0 quantities inherits both tails → C1 power ~0 by construction.

## R-18 caveats (honest limits)

- The uncarved reference is O3-collapsed (min adjacent overlap 0.065–0.126, |t|=0.28, sign-undetermined). C1 can only BOUND a carve bias; it cannot certify sub-kcal neutrality — fine for a SIGN/order verdict, NOT for a calibrated ddG.
- With |z|>2 single-seed outliers at n=6, the paired-t normality assumption is violated → all |t| (W23A's 6.25 and C1's carved 2.30) are outlier-sensitive. Only conclusions robust to single-seed removal, at SIGN/ORDER, may be trusted.
- W23A magnitude (-9.81 kcal) is ~2x the typical anchor (+4–5 destabilizing). Trp23 is p53's deepest MDM2 anchor so a large effect is plausible, but FF over-estimation of buried-hydrophobic mutations is likely → R-11 order-only.
- W23A large-negative sign is robust without its s101 outlier (4-seed = -8.26 kcal); MAGNITUDE is not (already forbidden by R-11). The carve is t=0-only, so a water can still diffuse into the displaced-copy void during dynamics (the residual s101/s127 tail source the carve leaves untouched).
- SIGN convention: W23A code-negative = destabilizing is consistent with the barnase fold ADR-0024 convention → known-direction anchors 2/2 (barnase I96A + W23A) support code-neg = destabilizing.

## Provenance

- Estimator: `outputs/_trackb/w4a_carved_c1_20260709/recompute_uncarved_ref.py` (reuses frozen `trackb_inplace_rbfe_production` + `uwham` verbatim; the paired delta + r use the same `paired_difference_stats`), mintimeid=100. This script (`w4a_c1_result.py`) does NOT edit the engine or the estimator.
- Carved run: `outputs/_trackb/w4a_carved_c1_20260709/` (launch 2026-07-09 07:37, 6 seeds s7/s23/s101/s127/s163/s199, 24-state, 400 cyc, `--carve-void-waters --carve-cutoff-nm 0.26`, disp 4.0, staged-min, reseed 20260618). carve OFF = existing uncarved build byte-match; carve ON removes free 4 / bound 3 bulk waters (net-charge 0, C2 all-bulk). NB review risk #4: per-leg carved-water counts must be nailed down + fail-loud logged (W23A carve is bound 0 / free 5).
- Uncarved reference: `uncarved_reference_dgb.json` (Pairing B). Carved raw: `carved_dgb.json`. Formal review: `internal C1 scientific-review note 20260710`.
- W23A sibling result: `outputs/_trackb/mdm2_w23a_gateA_20260708/w23a_ddg.json`.
