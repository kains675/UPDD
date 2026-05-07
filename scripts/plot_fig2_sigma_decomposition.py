"""scripts/plot_fig2_sigma_decomposition.py - Paper 1 Figure 2 plotting.

Generates the 2x2 sigma_btwn / sigma_w decomposition centerpiece figure:
  - top-left  (a): Cp4 n=5 cohort pre vs post charge-fix (sigma_btwn + sigma_w)
  - top-right (b): 1EBP_MTR13 n=5 cohort pre vs post charge-fix
  - bottom-left  (c): 7TL8_MTR6_s7 Phase-beta re-PBSA refresh (sigma_w only)
  - bottom-right (d): artifact-narrow -> genuine-broad conceptual schematic

Anchor sources (hardcoded for traceability and reproducibility):
  outputs/paper1/paper1_section3_manuscript.md
    - Section 3.2 (Cp4 n=5 cohort discovery anchor; 1EBP_MTR13 paired comparator)
    - Section 3.3 verdict (post-CONECT-v55 reference, paired pre/post sigmas)
    - Section 3.3.1 (7TL8_MTR6_s7 pre-/post-L387-v54 sigma_w refresh: 38.33 -> 12.13)

These constants mirror the published Section 3 manuscript anchor and bypass any
need to re-aggregate raw data from the analysis cache. They must remain in sync
with the manuscript; if the verdict numbers change, update both files together.

Usage:
    python scripts/plot_fig2_sigma_decomposition.py
    python scripts/plot_fig2_sigma_decomposition.py --output-dir outputs/paper1/figures
    python scripts/plot_fig2_sigma_decomposition.py --format png svg pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["axes.unicode_minus"] = False
matplotlib.rcParams["mathtext.fontset"] = "dejavusans"

# ---------------------------------------------------------------------------
# Anchor data (verdict Sec 3.3 / 3.3.1 - hardcoded for traceability)
# ---------------------------------------------------------------------------
CP4_N5 = {
    "pre":  {"sigma_btwn": 7.12, "sigma_w_median": 3.62},
    "post": {"sigma_btwn": 5.50, "sigma_w_median": 6.40},
}
EBP_MTR13_N5 = {
    "pre":  {"sigma_btwn": 2.82, "sigma_w_median": 3.24},
    "post": {"sigma_btwn": 3.00, "sigma_w_median": 4.91},
}
TL8_MTR6_S7_PHASE_B_REFRESH = {
    "pre":  {"sigma_w": 38.33},
    "post": {"sigma_w": 12.13},
}

# ---------------------------------------------------------------------------
# Style palette
# ---------------------------------------------------------------------------
COLOR_BTWN_LIGHT = "#7da5c9"  # navy light shade (pre)
COLOR_BTWN_SOLID = "#1f4e79"  # navy solid (post)
COLOR_W_LIGHT    = "#f5b962"  # orange light shade (pre)
COLOR_W_SOLID    = "#d97706"  # orange solid (post)
EDGE_COLOR       = "#222222"
ANNOT_COLOR      = "#222222"

Y_LABEL = r"$\sigma$ (kcal/mol)"


# ---------------------------------------------------------------------------
# Helper: paired-bar panel for Cp4 / 1EBP_MTR13
# ---------------------------------------------------------------------------
def _paired_sigma_bars(ax, data, title, annotations, ylim_max=None):
    """Draw grouped (pre, post) x (sigma_btwn, sigma_w) bars on ax.

    Parameters
    ----------
    ax : matplotlib Axes
    data : dict with keys 'pre' and 'post', each a dict with
        'sigma_btwn' and 'sigma_w_median'.
    title : str
    annotations : dict with optional keys 'btwn_delta' and 'w_delta'
        rendered above the post bars (per-pair annot_y so each label sits
        just above the taller bar of its pair instead of at a global height).
    ylim_max : float | None
        Optional explicit y-axis upper limit; if None, defaults to ymax * 1.50.
    """
    groups = ["$\\sigma_{btwn}$", "$\\sigma_{w}$ (median)"]
    pre_vals = [data["pre"]["sigma_btwn"], data["pre"]["sigma_w_median"]]
    post_vals = [data["post"]["sigma_btwn"], data["post"]["sigma_w_median"]]

    x = np.arange(len(groups))
    width = 0.36

    bars_pre = ax.bar(
        x - width / 2, pre_vals, width,
        color=[COLOR_BTWN_LIGHT, COLOR_W_LIGHT],
        edgecolor=EDGE_COLOR, linewidth=0.8,
        label="pre-patch",
    )
    bars_post = ax.bar(
        x + width / 2, post_vals, width,
        color=[COLOR_BTWN_SOLID, COLOR_W_SOLID],
        edgecolor=EDGE_COLOR, linewidth=0.8,
        label="post-patch",
    )

    # Numeric value labels on each bar
    for bar, val in zip(list(bars_pre) + list(bars_post), pre_vals + post_vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + 0.12,
            f"{val:.2f}",
            ha="center", va="bottom",
            fontsize=8.5, color=ANNOT_COLOR,
        )

    # Delta annotations above the post bars — per-pair annot_y so each label
    # attaches to its own bar pair (pair max + 0.6) instead of floating at a
    # shared height.
    ymax = max(pre_vals + post_vals)
    btwn_pair_max = max(pre_vals[0], post_vals[0])
    w_pair_max = max(pre_vals[1], post_vals[1])
    if "btwn_delta" in annotations:
        ax.text(
            x[0] + width / 2, btwn_pair_max + 0.6,
            annotations["btwn_delta"],
            ha="center", va="bottom",
            fontsize=8.5, color=COLOR_BTWN_SOLID, fontweight="bold",
        )
    if "w_delta" in annotations:
        ax.text(
            x[1] + width / 2, w_pair_max + 0.6,
            annotations["w_delta"],
            ha="center", va="bottom",
            fontsize=8.5, color=COLOR_W_SOLID, fontweight="bold",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=10)
    ax.set_ylabel(Y_LABEL, fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_ylim(0, ylim_max if ylim_max is not None else ymax * 1.50)
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", fontsize=8.5, framealpha=0.85)


# ---------------------------------------------------------------------------
# Panel a: Cp4
# ---------------------------------------------------------------------------
def panel_a_cp4(ax):
    annotations = {
        "btwn_delta": "-23 %\n(TIGHTER)",
        "w_delta": "+77 %\n(genuine sampling)",
    }
    _paired_sigma_bars(
        ax, CP4_N5,
        title="(a) Cp4 (n=5): pre vs post charge fix",
        annotations=annotations,
        ylim_max=9.0,
    )


# ---------------------------------------------------------------------------
# Panel b: 1EBP_MTR13
# ---------------------------------------------------------------------------
def panel_b_1ebp(ax):
    annotations = {
        "btwn_delta": "+6 %",
        "w_delta": "+52 %",
    }
    _paired_sigma_bars(
        ax, EBP_MTR13_N5,
        title="(b) 1EBP_MTR13 (n=5): pre vs post charge fix",
        annotations=annotations,
    )


# ---------------------------------------------------------------------------
# Panel c: 7TL8_MTR6_s7 sigma_w refresh
# ---------------------------------------------------------------------------
def panel_c_7tl8(ax):
    pre = TL8_MTR6_S7_PHASE_B_REFRESH["pre"]["sigma_w"]
    post = TL8_MTR6_S7_PHASE_B_REFRESH["post"]["sigma_w"]

    labels = ["pre-patch\n(pre-L387 v54)", "post-patch\n(post-L387 v54)"]
    vals = [pre, post]
    colors = [COLOR_W_LIGHT, COLOR_W_SOLID]

    x = np.arange(len(labels))
    width = 0.5
    bars = ax.bar(
        x, vals, width,
        color=colors, edgecolor=EDGE_COLOR, linewidth=0.8,
    )

    for bar, val in zip(bars, vals):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            bar.get_height() + 0.6,
            f"{val:.2f}",
            ha="center", va="bottom",
            fontsize=9.5, color=ANNOT_COLOR, fontweight="bold",
        )

    # Delta annotation: replace floating arrow with a centered text annotation
    # between the two bars (Option A — anchored cleanly between pre 38.33 and
    # post 12.13 without geometric ambiguity).
    ax.text(
        (x[0] + x[1]) / 2, 25.0,
        "-68 % (PBC defect removal)",
        ha="center", va="center",
        fontsize=10.5, color=COLOR_W_SOLID, fontweight="bold",
    )
    # Residual sigma_w note (compact form, smaller font, centered to fit panel)
    ax.axhline(post, linestyle="--", color="gray", linewidth=0.8, alpha=0.8)
    ax.text(
        0.5, 17.0,
        r"residual $\sigma_w \approx 12.13$ (genuine MD heterogeneity)",
        ha="center", va="center",
        fontsize=9, color="gray", style="italic",
    )

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9.5)
    ax.set_ylabel(Y_LABEL, fontsize=10)
    ax.set_title(
        r"(c) 7TL8_MTR6_s7 Phase-$\beta$ re-PBSA $\sigma_w$ refresh",
        fontsize=11, fontweight="bold",
    )
    ax.set_ylim(0, max(vals) * 1.25)
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.6)
    ax.set_axisbelow(True)


# ---------------------------------------------------------------------------
# Panel d: artifact-narrow -> genuine-broad conceptual schematic
# ---------------------------------------------------------------------------
def panel_d_schematic(ax):
    """Conceptual diagram: artifact-narrow basins vs genuine-broad ensemble."""
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(
        "(d) artifact-narrow $\\rightarrow$ genuine-broad ensemble",
        fontsize=11, fontweight="bold",
    )

    # Left: 3 narrow Gaussians at different x-positions (high sigma_btwn, low sigma_w)
    x_left = np.linspace(0, 4, 400)
    centers_left = [0.7, 2.0, 3.3]
    sigma_narrow = 0.18
    base_y = 1.2
    amp = 1.5
    colors_narrow = ["#7da5c9", "#5b8bb6", "#3a72a3"]
    for c, col in zip(centers_left, colors_narrow):
        y = amp * np.exp(-((x_left - c) ** 2) / (2 * sigma_narrow ** 2)) + base_y
        ax.plot(x_left, y, color=col, linewidth=1.6)
        ax.fill_between(x_left, base_y, y, color=col, alpha=0.25)
    ax.text(
        2.0, 0.55,
        "artifact-narrow basins\n(high $\\sigma_{btwn}$,  low $\\sigma_w$)",
        ha="center", va="center",
        fontsize=9, color=COLOR_BTWN_SOLID, fontweight="bold",
        linespacing=1.7,
    )
    # Mark sigma_btwn on left
    ax.annotate(
        "", xy=(centers_left[2], base_y - 0.15), xytext=(centers_left[0], base_y - 0.15),
        arrowprops=dict(arrowstyle="<->", color=COLOR_BTWN_SOLID, lw=1.4),
    )
    ax.text(
        np.mean(centers_left), base_y - 0.5,
        r"$\sigma_{btwn}$ large",
        ha="center", va="top",
        fontsize=8.5, color=COLOR_BTWN_SOLID,
    )

    # Right: 3 wide Gaussians at same center (low sigma_btwn, high sigma_w)
    x_right = np.linspace(6, 10, 400)
    center_right = 8.0
    sigma_wide = 0.85
    colors_wide = ["#f5b962", "#e89a40", "#d97706"]
    for sigma_jitter, col in zip([sigma_wide, sigma_wide * 1.05, sigma_wide * 1.1], colors_wide):
        y = amp * np.exp(-((x_right - center_right) ** 2) / (2 * sigma_jitter ** 2)) + base_y
        ax.plot(x_right, y, color=col, linewidth=1.6)
        ax.fill_between(x_right, base_y, y, color=col, alpha=0.20)
    ax.text(
        8.0, 0.55,
        "genuine-broad ensemble\n(low $\\sigma_{btwn}$,  high $\\sigma_w$)",
        ha="center", va="center",
        fontsize=9, color=COLOR_W_SOLID, fontweight="bold",
        linespacing=1.7,
    )
    # Mark sigma_w on right
    ax.annotate(
        "", xy=(center_right + 1.2, base_y + 0.6), xytext=(center_right - 1.2, base_y + 0.6),
        arrowprops=dict(arrowstyle="<->", color=COLOR_W_SOLID, lw=1.4),
    )
    ax.text(
        center_right, base_y + 0.95,
        r"$\sigma_w$ large",
        ha="center", va="bottom",
        fontsize=8.5, color=COLOR_W_SOLID,
    )

    # Center arrow: charge correction
    arrow = mpatches.FancyArrowPatch(
        (4.2, 3.0), (5.8, 3.0),
        arrowstyle="-|>", mutation_scale=22,
        color="#444444", lw=2.0,
    )
    ax.add_patch(arrow)
    ax.text(
        5.0, 3.45,
        "charge\ncorrection",
        ha="center", va="bottom",
        fontsize=9.5, color="#222222", fontweight="bold",
    )

    # (no explicit "ΔG coordinate" axis label — the basin distribution shapes
    # carry the σ_btwn / σ_w contrast on their own; an x-axis label without a
    # drawn axis was over-specification.)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Plot Paper 1 Figure 2 (sigma decomposition centerpiece).",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("outputs/paper1/figures"),
        help="Directory to write fig2_sigma_decomposition.{png,svg,...}",
    )
    parser.add_argument(
        "--format", nargs="+", default=["png", "svg"],
        help="Output formats (e.g. png svg pdf).",
    )
    parser.add_argument(
        "--dpi", type=int, default=300,
        help="DPI for raster outputs (default: 300).",
    )
    args = parser.parse_args(argv)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), dpi=args.dpi)
    panel_a_cp4(axes[0, 0])
    panel_b_1ebp(axes[0, 1])
    panel_c_7tl8(axes[1, 0])
    panel_d_schematic(axes[1, 1])

    fig.suptitle(
        r"Figure 2. $\sigma_{btwn}$ / $\sigma_{w}$ decomposition: "
        r"artifact-narrow $\rightarrow$ genuine-broad transition",
        fontsize=12.5, fontweight="bold", y=0.995,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.965))

    saved = []
    for fmt in args.format:
        out = args.output_dir / f"fig2_sigma_decomposition.{fmt}"
        save_kwargs = {"bbox_inches": "tight"}
        if fmt.lower() == "png":
            save_kwargs["dpi"] = args.dpi
        fig.savefig(out, **save_kwargs)
        saved.append(out)
        print(f"saved {out}")

    plt.close(fig)
    return saved


if __name__ == "__main__":
    main()
