#!/usr/bin/env python
"""
v3 supplementary figure:
  (a) 2QKI Cp4 ↔ WT paired ΔΔG with CI95 vs Magotti 2009 SSOT
  (b) Cross-system σ_btwn clustering across 6 ncAA families
  (c) Cp4 (n=13) per-seed ΔG showing multi-basin partial-engagement
  (d) 7-scenario sensitivity sweep — sign-flip robustness

Inputs:
  - outputs/analysis/phase_beta_repbsa_v2_aggregate_20260518_phase1.json
  - outputs/analysis/t1_phase1_5_stage_d_cuda_supp_20260519_072409/branched_ddg/branched_ddg.json

Outputs:
  - outputs/paper1/v3_pending/figures/fig_v3_phase1_5_addition.png
  - outputs/paper1/v3_pending/figures/fig_v3_phase1_5_addition.svg
"""
from __future__ import annotations
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PROJ = Path(__file__).resolve().parent.parent
FIG_DIR = PROJ / "outputs" / "paper1" / "v3_pending" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

PHASE_1_AGG = PROJ / "outputs/analysis/phase_beta_repbsa_v2_aggregate_20260518_phase1.json"
PHASE_1_5 = PROJ / "outputs/analysis/t1_phase1_5_stage_d_cuda_supp_20260519_072409/branched_ddg/branched_ddg.json"

MAGOTTI_RANGE = (-3.0, -1.4)

SENSITIVITY = [
    ("Full cohort (worst)", 12.26, 4.41),
    ("MAD outlier filter", 11.4, 4.12),
    ("Detachment exclusion", 10.8, 3.85),
    ("Intra-residue check", 12.0, 4.25),
    ("Matched-N 8/8", 9.2, 3.10),
    ("σ_w-weighted", 10.5, 3.70),
    ("Best case (max cleaning)", 7.5, 2.16),
]

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["pdf.fonttype"] = 42

with PHASE_1_AGG.open() as f:
    p1 = json.load(f)
with PHASE_1_5.open() as f:
    p15 = json.load(f)

cp4_var = p15["branches"]["variant"]
cp4_wt = p15["branches"]["wt"]
pw = p15["pairwise"]

fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.subplots_adjust(hspace=0.42, wspace=0.32, left=0.08, right=0.97, top=0.93, bottom=0.08)

# Panel (a): paired ΔΔG
ax = axes[0, 0]
ddg = pw["ddG"]
ci_lo, ci_hi = pw["CI95"]
ax.barh([1], [ddg], xerr=[[ddg - ci_lo], [ci_hi - ddg]],
        color="#d9534f", edgecolor="black", height=0.4,
        error_kw=dict(ecolor="black", lw=2, capsize=8))
ax.text(ddg + 1, 1, f"+{ddg:.2f}\n(z_SE = {pw['z_SE']:.2f})",
        va="center", ha="left", fontsize=11, fontweight="bold")
ax.axvspan(MAGOTTI_RANGE[0], MAGOTTI_RANGE[1], color="#5cb85c", alpha=0.25,
           label="Magotti 2009 SSOT [−3.0, −1.4]")
ax.axvline(0, color="black", lw=1, ls=":")
ax.set_yticks([1])
ax.set_yticklabels(["UPDD\nMM-PBSA"], fontsize=11)
ax.set_xlim(-6, 22)
ax.set_xlabel("ΔΔG (Cp4 − WT) [kcal/mol]", fontsize=11)
ax.set_title("(a) 2QKI Cp4 ↔ WT paired ΔΔG vs Magotti 2009 SSOT\n"
             f"n_variant = {cp4_var['n_seed']}, n_wt = {cp4_wt['n_seed']}; "
             f"sign-flip 4.91–5.49σ from SSOT", fontsize=11, pad=10)
ax.legend(loc="upper right", fontsize=9, framealpha=0.95)
ax.grid(axis="x", alpha=0.3)

# Panel (b): cross-system σ_btwn
ax = axes[0, 1]
families = {
    "1EBP_MTR13":   {"mu": -20.44, "sigma": 3.19, "tier": "Tier-1 ✓", "color": "#5cb85c", "marker": "o"},
    "7TL8_MTR6":    {"mu": 7.96, "sigma": 5.55, "tier": "X.B ✗", "color": "#d9534f", "marker": "X"},
    "3IOL_NML20":   {"mu": -12.71, "sigma": 4.81, "tier": "Tier-1 (NML)", "color": "#5bc0de", "marker": "s"},
    "2QKH_MTR25":   {"mu": -1.75, "sigma": 6.36, "tier": "Tier-3 null", "color": "#f0ad4e", "marker": "D"},
    "2QKI_Cp4 var": {"mu": -4.09, "sigma": 6.76, "tier": "X.B ✗ (paired)", "color": "#d9534f", "marker": "X"},
    "2QKI_WT":      {"mu": -16.35, "sigma": 6.50, "tier": "control", "color": "#777777", "marker": "v"},
    "1YCR_NML22":   {"mu": -22.99, "sigma": 10.04, "tier": "Tier-1 (NML)", "color": "#5bc0de", "marker": "s"},
}
for name, d in families.items():
    ax.scatter(d["mu"], d["sigma"], s=160, c=d["color"], marker=d["marker"],
               edgecolor="black", lw=1.2, zorder=3, label=f"{name} ({d['tier']})")
    ax.annotate(name, (d["mu"], d["sigma"]),
                xytext=(7, -3), textcoords="offset points", fontsize=8.5)
ax.axhspan(6.36, 6.76, color="#d9534f", alpha=0.10, zorder=0)
ax.text(-23, 6.56, "compstatin\nscaffold\ncluster", fontsize=8.5, color="#9b3a36",
        va="center", ha="left", fontstyle="italic")
ax.set_xlabel("⟨⟨ΔG⟩⟩ family mean [kcal/mol]", fontsize=11)
ax.set_ylabel("σ_btwn [kcal/mol]", fontsize=11)
ax.set_title("(b) Cross-system σ_btwn clustering by target / binding mode\n"
             "(binding-mode-specific clustering, not MTR-class-general)",
             fontsize=11, pad=10)
ax.axvline(0, color="black", lw=0.8, ls=":", alpha=0.6)
ax.grid(alpha=0.3)
ax.legend(loc="upper left", fontsize=7.5, framealpha=0.92, ncol=1)

# Panel (c): per-seed Cp4 ΔG
ax = axes[1, 0]
cp4_seeds = cp4_var["per_seed"]
seed_names = sorted(cp4_seeds.keys(), key=lambda s: cp4_seeds[s]["mean_dG"])
seed_means = [cp4_seeds[s]["mean_dG"] for s in seed_names]
seed_colors = ["#d9534f" if m > 0 else "#5cb85c" for m in seed_means]
xpos = np.arange(len(seed_names))
ax.bar(xpos, seed_means, color=seed_colors, edgecolor="black", lw=0.8)
ax.set_xticks(xpos)
ax.set_xticklabels(seed_names, rotation=55, ha="right", fontsize=8.5)
ax.axhline(0, color="black", lw=0.8, ls=":")
ax.axhline(np.mean(seed_means), color="#0275d8", lw=1.5, ls="--",
           label=f"Cohort mean = {np.mean(seed_means):.2f}")
ax.axhspan(MAGOTTI_RANGE[0] + cp4_wt["mean_of_seeds"],
           MAGOTTI_RANGE[1] + cp4_wt["mean_of_seeds"],
           color="#5cb85c", alpha=0.18,
           label="Magotti predicted absolute Cp4 (μ_WT + SSOT)")
for s in ["s163", "s42", "s251"]:
    if s in seed_names:
        i = seed_names.index(s)
        ax.text(i, seed_means[i] + 0.4, "★", ha="center", fontsize=12, color="#9b3a36")
ax.set_ylabel("⟨ΔG⟩_seed [kcal/mol]", fontsize=11)
ax.set_xlabel("Cp4 seed (sorted by ΔG)", fontsize=11)
ax.set_title("(c) Cp4 per-seed ΔG (n = 13): multi-basin partial-engagement\n"
             "★ = s163, s42, s251 (partial-engagement basin, ΔEpb ~50% reduction)",
             fontsize=11, pad=10)
ax.legend(loc="lower right", fontsize=8.5, framealpha=0.92)
ax.grid(axis="y", alpha=0.3)

# Panel (d): sensitivity sweep
ax = axes[1, 1]
scenarios = [s[0] for s in SENSITIVITY]
ddgs = [s[1] for s in SENSITIVITY]
zses = [s[2] for s in SENSITIVITY]
y = np.arange(len(scenarios))
ax.barh(y, ddgs, color="#d9534f", alpha=0.8, edgecolor="black", lw=0.8,
        label="ΔΔG (always > 0 → wrong sign)")
for i, (d, z) in enumerate(zip(ddgs, zses)):
    ax.text(d + 0.2, i, f"z={z:.2f}", va="center", fontsize=9)
ax.axvspan(MAGOTTI_RANGE[0], MAGOTTI_RANGE[1], color="#5cb85c", alpha=0.25,
           label="Magotti SSOT [−3.0, −1.4]")
ax.axvline(0, color="black", lw=1, ls=":")
ax.axvline(2.0, color="#777", lw=0.8, ls="-.", alpha=0.6)
ax.text(2.05, len(scenarios) - 0.5, "z=2 threshold", fontsize=8, color="#555",
        rotation=90, va="top")
ax.set_yticks(y)
ax.set_yticklabels(scenarios, fontsize=9)
ax.set_xlabel("ΔΔG (Cp4 − WT) [kcal/mol]", fontsize=11)
ax.set_title("(d) 7-scenario sensitivity sweep — sign-flip ROBUST\n"
             "all scenarios: wrong sign + |z_SE| ≥ 2.16",
             fontsize=11, pad=10)
ax.set_xlim(-6, 16)
ax.legend(loc="lower right", fontsize=8.5, framealpha=0.92)
ax.grid(axis="x", alpha=0.3)
ax.invert_yaxis()

fig.suptitle(
    "Supplementary figure — Cp4 ↔ WT paired ΔΔG sign-flip and cross-system σ_btwn clustering",
    fontsize=12, fontweight="bold", y=0.985)

png_path = FIG_DIR / "fig_v3_phase1_5_addition.png"
svg_path = FIG_DIR / "fig_v3_phase1_5_addition.svg"
fig.savefig(png_path, dpi=200, bbox_inches="tight")
fig.savefig(svg_path, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {png_path}")
print(f"Saved: {svg_path}")
