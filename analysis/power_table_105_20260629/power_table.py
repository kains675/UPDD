#!/usr/bin/env python
"""
#105 |ddG| detection-limit x GPU-day statistical power table for UPDD Track B.
#107 Cp01->V3I go/no-go resolution.

Ranking-only regime (R-11): this is a SAMPLING-DESIGN tool, NOT a magnitude/sign claim.
Honesty (R-18): multi-basin sigma rows are marked undetermined regardless of N.

All inputs are VERIFIED grounding + SciVal-prescribed model. No value is re-guessed.

ESTIMATOR-CLASS WARNING (SciVal, R-18):
  - sigma_btwn in {1.2, 1.4} are Track B two-copy ATS ABFE PAIRED-leg between-seed SDs
    (paired_SD of per-seed ddG_bind = bound - free). Canonical class, same 2QKI engine.
  - sigma_btwn in {3.0, 4.0, 6.5} approximate the MM-PBSA endpoint-cohort SDs
    (MTR13~3.19, MTR25~6.36, Cp4~6.66). DIFFERENT estimator class, multi-basin systems.
  - These two classes are NOT one comparable ladder. The table labels them separately.

MODEL (SciVal-adopted):
  - PAIRED is mandatory for Track B ABFE: SE = paired_SD/sqrt(n), Student-t df=n-1.
  - The task's step-2 formula min_detectable_ddg = crit * sqrt(SE_cp4^2 + SE_wt^2) is the
    QUADRATURE form (MM-PBSA branched_ddg.py). SciVal verdict 'accept_with_caveat' =>
    we report BOTH the paired floor and the quadrature floor side by side.
  - crit: SciVal cutoff verdict 'accept_with_caveat' => use Student-t (NOT fixed 3.0).
    We show t-crit AND the 3-sigma fixed value for transparency.
"""
import math
import json

# ---------------------------------------------------------------------------
# VERIFIED grounding constants
# ---------------------------------------------------------------------------
BOUND_H_PER_SEED = 5.1            # measured (s7=4.92h, s163=5.28h; 14-state/800cyc, serial)
FREE_H_PER_SEED = 0.9            # measured (400cyc, 12-state, 2-concurrent wall)
TOTAL_H_PER_SEED = BOUND_H_PER_SEED + FREE_H_PER_SEED   # 6.0 = ONE seed's full RBFE (bound+free)
H18_SCALING_FACTOR = 1.286        # state-ratio 36/28 (per-state linearity unanchored); BOUND leg only
H18_FACTOR_PRESENT = True
# Per-seed H18 cost: free leg (12-state) unchanged; only the bound leg densifies (28->36 states).
H18_H_PER_SEED = BOUND_H_PER_SEED * H18_SCALING_FACTOR + FREE_H_PER_SEED

# Cp01->V3I direct measured paired SD (n=6, JSON SSOT v3i_paired_n6_verdict_20260616.json)
V3I_SIGMA = 1.40                  # 1.398 rounded to 2 sig figs
V3I_N = 6
# Experimental |ddG| target (Lambris ~30x, 298K convention). MEDIUM-LOW confidence.
DDG_EXP_KCAL = 2.02

# SE_target grid and sigma grid
SIGMA_GRID = [1.2, 1.4, 2.0, 3.0, 4.0, 6.5]
SE_TARGET_GRID = [0.5, 1.0]

# Estimator-class tag per sigma (for honest labeling, R-18)
def estimator_class(sigma):
    if sigma <= 1.6:
        return "paired-ATS-ABFE (canonical)"
    elif sigma < 3.0:
        return "intermediate (interpolated)"
    else:
        return "MM-PBSA-endpoint-cohort (multi-basin)"

# ---------------------------------------------------------------------------
# Student-t two-tailed alpha=0.05 critical values (LUT, df 1..30), z for df>30
# ---------------------------------------------------------------------------
_T_TABLE_975 = {
    1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
    8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145,
    15: 2.131, 16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093, 20: 2.086,
    21: 2.080, 22: 2.074, 23: 2.069, 24: 2.064, 25: 2.060, 26: 2.056,
    27: 2.052, 28: 2.048, 29: 2.045, 30: 2.042,
}
_Z_LARGE = 1.960

def t_crit(df):
    if df < 1:
        return None
    if df <= 30:
        return _T_TABLE_975[int(df)]
    return _Z_LARGE

# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------
def compute_row(sigma_btwn, se_target):
    # N per endpoint so that SEM <= SE_target
    n = math.ceil((sigma_btwn / se_target) ** 2)

    # crit: Student-t. For the PAIRED model df = n-1.
    # For the QUADRATURE (Welch / branched_ddg conservative) form df = 2*(n-1).
    df_paired = max(1, n - 1)
    df_quad = max(1, 2 * (n - 1))
    tcrit_paired = t_crit(df_paired)
    tcrit_quad = t_crit(df_quad)

    # Min detectable ddG (SciVal: report BOTH estimators since quadrature accept_with_caveat)
    # PAIRED floor: crit * SEM (SEM = se_target by construction of N)
    min_ddg_paired_t = tcrit_paired * se_target
    # QUADRATURE floor (task step-2 explicit formula): crit * sqrt(SE_cp4^2 + SE_wt^2)
    se_quad = math.sqrt(se_target ** 2 + se_target ** 2)  # = se_target * sqrt(2)
    min_ddg_quad_t = tcrit_quad * se_quad
    # Fixed 3-sigma reference (what SciVal said NOT to use, shown for transparency)
    min_ddg_quad_3sig = 3.0 * se_quad

    # The primary min_detectable_ddg per the production code = PAIRED Student-t floor.
    min_detectable_ddg = round(min_ddg_paired_t, 4)

    # GPU-day costs. PAIRED RBFE: one seed = one bound + one free leg = TOTAL_H_PER_SEED.
    # N = number of PAIRED seeds (NOT per-endpoint x2). Earlier x2 was a double-count (R-18 fix).
    gpu_day_measured = n * TOTAL_H_PER_SEED / 24.0
    gpu_day_h18 = n * H18_H_PER_SEED / 24.0

    # Tier classification (R-18 honest)
    eclass = estimator_class(sigma_btwn)
    if sigma_btwn >= 3.0:
        # Multi-basin: N-scaling INVALID (sigma_w ~ sigma_btwn, basin-hopping floor).
        # Marked undetermined regardless of N.
        tier = "undetermined(multi-basin floor; N cannot beat it)"
    elif n > 40:
        tier = "undetermined(N infeasible >40)"
    elif n <= 12:
        tier = "tight(detectable, N<=12)"
    else:
        tier = "multi-basin-warn(N 13-40, sampling-design feasible)"

    return {
        "sigma_btwn": sigma_btwn,
        "estimator_class": eclass,
        "se_target": se_target,
        "n_per_endpoint": int(n),
        "df_paired": int(df_paired),
        "t_crit_paired": round(tcrit_paired, 3) if tcrit_paired else None,
        "min_detectable_ddg_kcal": min_detectable_ddg,
        "min_ddg_quad_t_kcal": round(min_ddg_quad_t, 4),
        "min_ddg_quad_3sig_kcal": round(min_ddg_quad_3sig, 4),
        "gpu_day_measured_basis": round(gpu_day_measured, 3),
        "gpu_day_h18_scaled": round(gpu_day_h18, 3),
        "tier": tier,
    }

rows = []
for sigma in SIGMA_GRID:
    for se in SE_TARGET_GRID:
        rows.append(compute_row(sigma, se))

# ---------------------------------------------------------------------------
# Cp01->V3I row (#107 go/no-go)
# ---------------------------------------------------------------------------
def cp01v3i_verdict():
    sigma_used = V3I_SIGMA  # direct measured paired SD (n=6), NOT a borrowed proxy
    caveats = [
        "sigma_btwn=1.40 is the DIRECT measured between-seed SD of the exact Cp01->V3I cohort (n=6, paired_SD=1.398, JSON SSOT v3i_paired_n6_verdict_20260616.json) -- NOT a borrowed proxy; W4A 1.20 (n=4, prose-only) is the canonical-class fallback floor.",
        "Estimator class = paired two-copy ATS ABFE (matched-seed bound-free diff SEM, Student-t). The production code uses PAIRED as primary; Welch quadrature survives only as ddint_err_quadrature_legacy. Power table built on the PAIRED model accordingly.",
        "ddg_exp=2.02 kcal/mol is MEDIUM-LOW confidence: back-derived from Lambris 2025 '~30x tighter' ratio (preprint, not peer-reviewed); RT=0.593 at 298K vs SPR at 310K (~4% offset -> 2.10 at 310K); true ratio plausibly 20-50x -> |ddG| 1.8-2.3.",
        "V3I calc sign is DISCORDANT (positive => Ile3 unfavorable in calc) with the favorable-Ile a-priori prior. Sign-RESOLUTION power != anchor AGREEMENT (R-11; magnitude forbidden).",
        "N=(sigma/SE)^2 is defensible for V3I ONLY because the n=6 pilot shows it is unimodal/sign-stable (sigma_btwn 1.40 not inflated by sign-flipping basins). This is known POST-pilot, not a priori.",
        "Ranking-only regime (R-11): table is a sampling-design tool. No magnitude/sign claim is made.",
    ]

    if DDG_EXP_KCAL is None:
        return {
            "sigma_used": sigma_used,
            "tier": "blocked-need-reference-Kd",
            "honest_caveats": caveats + ["Missing: cleanly recorded parent Cp01 absolute K_D; only ~30x ratio recorded."],
            "ddg_exp_kcal": None,
            "required_n_per_endpoint": None,
            "gpu_day_estimate": None,
            "clears_detection_limit": None,
        }

    # Find min N (per endpoint) where the PAIRED Student-t min_detectable_ddg <= ddg_exp.
    # SEM = sigma_used / sqrt(N); floor = t_crit(N-1) * SEM.
    required_n = None
    for n in range(2, 201):
        sem = sigma_used / math.sqrt(n)
        tc = t_crit(max(1, n - 1))
        floor = tc * sem
        if floor <= DDG_EXP_KCAL:
            required_n = n
            break

    if required_n is None:
        return {
            "sigma_used": sigma_used,
            "tier": "undetermined(N>200 infeasible)",
            "honest_caveats": caveats,
            "ddg_exp_kcal": DDG_EXP_KCAL,
            "required_n_per_endpoint": None,
            "gpu_day_estimate": None,
            "clears_detection_limit": False,
        }

    sem = sigma_used / math.sqrt(required_n)
    tc = t_crit(max(1, required_n - 1))
    floor = tc * sem
    # GPU-day estimate (H18-scaled). Bound-leg densified; free leg reused/unchanged.
    # N = paired seeds (NOT per-endpoint x2 — R-18 double-count fix).
    gpu_day = required_n * H18_H_PER_SEED / 24.0

    # Tier for V3I: sigma 1.40 is canonical paired class (not multi-basin)
    if required_n <= 12:
        tier = "tight(detectable at N<=12, sign-resolution feasible)"
    elif required_n <= 40:
        tier = "feasible(N 13-40)"
    else:
        tier = "undetermined(N infeasible >40)"

    return {
        "sigma_used": sigma_used,
        "tier": tier,
        "honest_caveats": caveats,
        "ddg_exp_kcal": DDG_EXP_KCAL,
        "required_n_per_endpoint": int(required_n),
        "min_detectable_ddg_at_required_n": round(floor, 4),
        "gpu_day_estimate": round(gpu_day, 3),
        "clears_detection_limit": True,
    }

v3i = cp01v3i_verdict()

# ---------------------------------------------------------------------------
# funnel cutoff (SciVal: t_crit * sigma / sqrt(N); ~1.1 kcal at N=8)
# ---------------------------------------------------------------------------
# SciVal: "forbid sign claims below ~1.1 kcal at N=8". Compute concretely with the
# canonical-class sigma (V3I 1.40) at N=8: floor = t_crit(7) * 1.40 / sqrt(8).
FUNNEL_N = 8
funnel_floor = t_crit(FUNNEL_N - 1) * V3I_SIGMA / math.sqrt(FUNNEL_N)
funnel_cutoff = round(funnel_floor, 2)

# ---------------------------------------------------------------------------
# Markdown table
# ---------------------------------------------------------------------------
def make_markdown():
    lines = []
    lines.append("# #105 |ddG| Detection-Limit x GPU-day Power Table (UPDD Track B)")
    lines.append("")
    lines.append("**Regime:** ranking-only (R-11) — sampling-design tool, NOT a magnitude/sign claim.")
    lines.append("")
    lines.append("**Model (SciVal-adopted):** PAIRED estimator primary (SE = paired_SD/sqrt(n), Student-t df=n-1). "
                 "Quadrature floor shown as `accept_with_caveat` secondary. crit = Student-t (NOT fixed 3-sigma); "
                 "3-sigma shown for transparency only.")
    lines.append("")
    lines.append("**Cost basis:** total_h_per_seed=6.0 (measured: bound~5.1 + free~0.9, 14-state/800cyc). "
                 "H18 scaling factor=1.286 (state-ratio 36/28; per-state linearity unanchored).")
    lines.append("")
    lines.append("**ESTIMATOR-CLASS WARNING (R-18):** sigma 1.2/1.4 = paired ATS-ABFE (canonical, same engine). "
                 "sigma 3.0/4.0/6.5 ~ MM-PBSA endpoint-cohort (multi-basin: MTR13~3.19, MTR25~6.36, Cp4~6.66). "
                 "These are DIFFERENT physical estimators, not one comparable ladder. sigma>=3.0 rows marked "
                 "`undetermined` because N-scaling cannot beat a basin-hopping floor (sigma_w ~ sigma_btwn).")
    lines.append("")
    header = ("| sigma_btwn | estimator class | SE_target | N (paired seeds) | df(paired) | t_crit | "
              "min_detectable_ddG (paired-t, kcal) | min_ddG (quad-t) | min_ddG (3sig) | "
              "GPU-day (measured) | GPU-day (H18) | tier |")
    sep = "|" + "---|" * 13
    lines.append(header)
    lines.append(sep)
    for r in rows:
        lines.append(
            f"| {r['sigma_btwn']:.1f} | {r['estimator_class']} | {r['se_target']:.1f} | "
            f"{r['n_per_endpoint']} | {r['df_paired']} | {r['t_crit_paired']} | "
            f"{r['min_detectable_ddg_kcal']:.3f} | {r['min_ddg_quad_t_kcal']:.3f} | "
            f"{r['min_ddg_quad_3sig_kcal']:.3f} | {r['gpu_day_measured_basis']:.3f} | "
            f"{r['gpu_day_h18_scaled']:.3f} | {r['tier']} |"
        )
    lines.append("")
    lines.append("## #107 Cp01->V3I go/no-go")
    lines.append("")
    lines.append(f"- **sigma_used:** {v3i['sigma_used']} kcal/mol (DIRECT measured paired_SD, n=6, JSON SSOT)")
    lines.append(f"- **ddg_exp:** {v3i['ddg_exp_kcal']} kcal/mol (MEDIUM-LOW confidence, Lambris ~30x, 298K)")
    lines.append(f"- **clears_detection_limit:** {v3i['clears_detection_limit']}")
    lines.append(f"- **required_n_per_endpoint:** {v3i['required_n_per_endpoint']}")
    if 'min_detectable_ddg_at_required_n' in v3i:
        lines.append(f"- **min_detectable_ddG at required N:** {v3i['min_detectable_ddg_at_required_n']} kcal/mol")
    lines.append(f"- **gpu_day_estimate (H18-scaled, bound-leg densified, free reused):** {v3i['gpu_day_estimate']}")
    lines.append(f"- **tier:** {v3i['tier']}")
    lines.append("")
    lines.append("### Honest caveats (R-18)")
    for c in v3i['honest_caveats']:
        lines.append(f"- {c}")
    lines.append("")
    lines.append(f"## Funnel cutoff: |ddG| floor = {funnel_cutoff} kcal/mol")
    lines.append("")
    lines.append(f"SciVal-prescribed: floor = t_crit(df=N-1) * sigma / sqrt(N). At N={FUNNEL_N}, sigma={V3I_SIGMA} "
                 f"(canonical paired class): t_crit(7)={t_crit(7)} * {V3I_SIGMA}/sqrt({FUNNEL_N}) = {funnel_cutoff} kcal/mol. "
                 "Below this floor, sign-claims are FORBIDDEN (ranking-only). Permit a sign claim only when the "
                 "Student-t CI excludes zero. The floor is SOFT (sigma is 2-sig-fig noisy at n=6, df=5).")
    return "\n".join(lines)

md = make_markdown()
with open("/home/san/UPDD_proj/analysis/power_table_105_20260629/power_table.md", "w") as f:
    f.write(md)

# Emit JSON for the harness to consume
out = {
    "table_rows": rows,
    "cp01v3i_verdict": v3i,
    "funnel_cutoff_ddg_kcal": funnel_cutoff,
}
with open("/home/san/UPDD_proj/analysis/power_table_105_20260629/power_table_result.json", "w") as f:
    json.dump(out, f, indent=2)

print(md)
print("\n\n=== JSON ===")
print(json.dumps(out, indent=2))
