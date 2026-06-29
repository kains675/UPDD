# #105 |ddG| Detection-Limit x GPU-day Power Table (UPDD Track B)

**Regime:** ranking-only (R-11) — sampling-design tool, NOT a magnitude/sign claim.

**Model (SciVal-adopted):** PAIRED estimator primary (SE = paired_SD/sqrt(n), Student-t df=n-1). Quadrature floor shown as `accept_with_caveat` secondary. crit = Student-t (NOT fixed 3-sigma); 3-sigma shown for transparency only.

**Cost basis:** total_h_per_seed=6.0 (measured: bound~5.1 + free~0.9, 14-state/800cyc). H18 scaling factor=1.286 (state-ratio 36/28; per-state linearity unanchored).

**ESTIMATOR-CLASS WARNING (R-18):** sigma 1.2/1.4 = paired ATS-ABFE (canonical, same engine). sigma 3.0/4.0/6.5 ~ MM-PBSA endpoint-cohort (multi-basin: MTR13~3.19, MTR25~6.36, Cp4~6.66). These are DIFFERENT physical estimators, not one comparable ladder. sigma>=3.0 rows marked `undetermined` because N-scaling cannot beat a basin-hopping floor (sigma_w ~ sigma_btwn).

| sigma_btwn | estimator class | SE_target | N (paired seeds) | df(paired) | t_crit | min_detectable_ddG (paired-t, kcal) | min_ddG (quad-t) | min_ddG (3sig) | GPU-day (measured) | GPU-day (H18) | tier |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.2 | paired-ATS-ABFE (canonical) | 0.5 | 6 | 5 | 2.571 | 1.286 | 1.575 | 2.121 | 1.500 | 1.865 | tight(detectable, N<=12) |
| 1.2 | paired-ATS-ABFE (canonical) | 1.0 | 2 | 1 | 12.706 | 12.706 | 6.085 | 4.243 | 0.500 | 0.622 | tight(detectable, N<=12) |
| 1.4 | paired-ATS-ABFE (canonical) | 0.5 | 8 | 7 | 2.365 | 1.183 | 1.517 | 2.121 | 2.000 | 2.486 | tight(detectable, N<=12) |
| 1.4 | paired-ATS-ABFE (canonical) | 1.0 | 2 | 1 | 12.706 | 12.706 | 6.085 | 4.243 | 0.500 | 0.622 | tight(detectable, N<=12) |
| 2.0 | intermediate (interpolated) | 0.5 | 16 | 15 | 2.131 | 1.065 | 1.444 | 2.121 | 4.000 | 4.972 | multi-basin-warn(N 13-40, sampling-design feasible) |
| 2.0 | intermediate (interpolated) | 1.0 | 4 | 3 | 3.182 | 3.182 | 3.461 | 4.243 | 1.000 | 1.243 | tight(detectable, N<=12) |
| 3.0 | MM-PBSA-endpoint-cohort (multi-basin) | 0.5 | 36 | 35 | 1.96 | 0.980 | 1.386 | 2.121 | 9.000 | 11.188 | undetermined(multi-basin floor; N cannot beat it) |
| 3.0 | MM-PBSA-endpoint-cohort (multi-basin) | 1.0 | 9 | 8 | 2.306 | 2.306 | 2.998 | 4.243 | 2.250 | 2.797 | undetermined(multi-basin floor; N cannot beat it) |
| 4.0 | MM-PBSA-endpoint-cohort (multi-basin) | 0.5 | 64 | 63 | 1.96 | 0.980 | 1.386 | 2.121 | 16.000 | 19.890 | undetermined(multi-basin floor; N cannot beat it) |
| 4.0 | MM-PBSA-endpoint-cohort (multi-basin) | 1.0 | 16 | 15 | 2.131 | 2.131 | 2.888 | 4.243 | 4.000 | 4.972 | undetermined(multi-basin floor; N cannot beat it) |
| 6.5 | MM-PBSA-endpoint-cohort (multi-basin) | 0.5 | 169 | 168 | 1.96 | 0.980 | 1.386 | 2.121 | 42.250 | 52.521 | undetermined(multi-basin floor; N cannot beat it) |
| 6.5 | MM-PBSA-endpoint-cohort (multi-basin) | 1.0 | 43 | 42 | 1.96 | 1.960 | 2.772 | 4.243 | 10.750 | 13.363 | undetermined(multi-basin floor; N cannot beat it) |

## #107 Cp01->V3I go/no-go

- **sigma_used:** 1.4 kcal/mol (DIRECT measured paired_SD, n=6, JSON SSOT)
- **ddg_exp:** 2.02 kcal/mol (MEDIUM-LOW confidence, Lambris ~30x, 298K)
- **clears_detection_limit:** True
- **required_n_per_endpoint:** 5
- **min_detectable_ddG at required N:** 1.7381 kcal/mol
- **gpu_day_estimate (H18-scaled, bound-leg densified, free reused):** 1.554
- **tier:** tight(detectable at N<=12, sign-resolution feasible)

### Honest caveats (R-18)
- sigma_btwn=1.40 is the DIRECT measured between-seed SD of the exact Cp01->V3I cohort (n=6, paired_SD=1.398, JSON SSOT v3i_paired_n6_verdict_20260616.json) -- NOT a borrowed proxy; W4A 1.20 (n=4, prose-only) is the canonical-class fallback floor.
- Estimator class = paired two-copy ATS ABFE (matched-seed bound-free diff SEM, Student-t). The production code uses PAIRED as primary; Welch quadrature survives only as ddint_err_quadrature_legacy. Power table built on the PAIRED model accordingly.
- ddg_exp=2.02 kcal/mol is MEDIUM-LOW confidence: back-derived from Lambris 2025 '~30x tighter' ratio (preprint, not peer-reviewed); RT=0.593 at 298K vs SPR at 310K (~4% offset -> 2.10 at 310K); true ratio plausibly 20-50x -> |ddG| 1.8-2.3.
- V3I calc sign is DISCORDANT (positive => Ile3 unfavorable in calc) with the favorable-Ile a-priori prior. Sign-RESOLUTION power != anchor AGREEMENT (R-11; magnitude forbidden).
- N=(sigma/SE)^2 is defensible for V3I ONLY because the n=6 pilot shows it is unimodal/sign-stable (sigma_btwn 1.40 not inflated by sign-flipping basins). This is known POST-pilot, not a priori.
- Ranking-only regime (R-11): table is a sampling-design tool. No magnitude/sign claim is made.

## Funnel cutoff: |ddG| floor = 1.17 kcal/mol

SciVal-prescribed: floor = t_crit(df=N-1) * sigma / sqrt(N). At N=8, sigma=1.4 (canonical paired class): t_crit(7)=2.365 * 1.4/sqrt(8) = 1.17 kcal/mol. Below this floor, sign-claims are FORBIDDEN (ranking-only). Permit a sign claim only when the Student-t CI excludes zero. The floor is SOFT (sigma is 2-sig-fig noisy at n=6, df=5).