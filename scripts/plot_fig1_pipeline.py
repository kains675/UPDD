"""scripts/plot_fig1_pipeline.py - Paper 1 Figure 1 plotting.

Generates the 5-stage UPDD pipeline architecture diagram with branching at Stage 4:
  - Stage 1: Sequence Input Parser
  - Stage 2: Generation (RFdiffusion + ProteinMPNN, established external)
  - Stage 3: AF2 Structure Prediction (ColabFold dispatch, established external)
  - Stage 4: Evaluation (UPDD primary contribution; branched WT vs Variant)
      - sub-steps: ncAA parameterization -> MD -> QM/MM -> MM-PBSA
  - Stage 5: Wet-lab K_d feedback ingestion (cycle progression)

This is a vector architecture diagram (not a data plot).

Usage:
    python scripts/plot_fig1_pipeline.py
    python scripts/plot_fig1_pipeline.py --output-dir /tmp/fig_test
    python scripts/plot_fig1_pipeline.py --format png svg pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["axes.unicode_minus"] = False
matplotlib.rcParams["mathtext.fontset"] = "dejavusans"

# ---------------------------------------------------------------------------
# Style palette
# ---------------------------------------------------------------------------
COLOR_NAVY = "#1f4e79"          # UPDD-primary (Stage 4 sub-steps)
COLOR_NAVY_LIGHT = "#dbe7f3"
COLOR_GRAY = "#aaaaaa"          # external / boundary
COLOR_GRAY_LIGHT = "#ececec"
EDGE_COLOR = "#222222"
TEXT_COLOR = "#1a1a1a"
ANNOT_COLOR = "#444444"

# ---------------------------------------------------------------------------
# Stage definitions (top-to-bottom)
# ---------------------------------------------------------------------------
EXTERNAL_STAGES = [
    {
        "label": "Stage 1",
        "title": "Sequence Input Parser",
        "subtitle": "canonical design specification\n(ncAA registry validation)",
    },
    {
        "label": "Stage 2",
        "title": "Generation: RFdiffusion + ProteinMPNN",
        "subtitle": "established external tools\n(structure backbone + sequence design)",
    },
    {
        "label": "Stage 3",
        "title": "AF2 Structure Prediction (ColabFold)",
        "subtitle": "established external tool\n(per-residue pLDDT + PAE)",
    },
]

# Stage 4 sub-steps (UPDD primary contribution; branched WT | Variant for first three)
STAGE4_SUBSTEPS = [
    {
        "key": "ncaa_param",
        "title": "ncAA parameterization",
        "annotation": "R-15/R-16 charge guards\n(chemistry-true charge SSOT primitive)",
    },
    {
        "key": "md",
        "title": "Molecular Dynamics",
        "annotation": "L387 v54 PBC repair + CONECT v55\natom-index disambiguation\n+ 3-pass timestep recovery",
    },
    {
        "key": "qmmm",
        "title": "QM/MM single-point",
        "annotation": "ωB97X-D3 / def2-SVP\nlink-atom canonical reference\n(Senn-Thiel 2009)",
    },
    {
        "key": "mmpbsa",
        "title": "MM-PBSA scoring + ΔΔG",
        "annotation": "branched ΔΔG with bit-identical WT control\nσ_btwn / σ_w decomposition\n+ Convergence Index dual framing",
    },
]

STAGE5 = {
    "label": "Stage 5",
    "title": "Wet-lab K_d feedback ingestion",
    "subtitle": "multi-assay aggregation -> cycle progression\n(Cycle 0 -> Cycle 1+ -> plateau)",
}


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def _draw_stage_box(ax, x_center, y_center, width, height, stage,
                    is_primary=False, fontsize_title=10):
    """Draw a single rounded stage box. is_primary uses navy; otherwise gray."""
    if is_primary:
        edge = COLOR_NAVY
        face = COLOR_NAVY_LIGHT
        title_color = COLOR_NAVY
    else:
        edge = COLOR_GRAY
        face = COLOR_GRAY_LIGHT
        title_color = "#333333"

    box = FancyBboxPatch(
        (x_center - width / 2, y_center - height / 2),
        width, height,
        boxstyle="round,pad=0.02,rounding_size=0.12",
        linewidth=1.4,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(box)

    # Stage label badge (left side)
    if "label" in stage:
        ax.text(
            x_center - width / 2 + 0.12, y_center + height / 2 - 0.18,
            stage["label"],
            ha="left", va="top",
            fontsize=8.5, fontweight="bold",
            color=title_color,
        )
    # Title
    ax.text(
        x_center, y_center + height * 0.12,
        stage["title"],
        ha="center", va="center",
        fontsize=fontsize_title, fontweight="bold",
        color=title_color,
    )
    # Subtitle / annotation
    sub_key = "subtitle" if "subtitle" in stage else "annotation"
    ax.text(
        x_center, y_center - height * 0.22,
        stage[sub_key],
        ha="center", va="center",
        fontsize=7.5, color=TEXT_COLOR,
    )


def _draw_substep_box(ax, x_center, y_center, width, height, substep, branch_label=None):
    """Draw a Stage 4 sub-step box with optional branch label (Branch A / B)."""
    box = FancyBboxPatch(
        (x_center - width / 2, y_center - height / 2),
        width, height,
        boxstyle="round,pad=0.02,rounding_size=0.10",
        linewidth=1.2,
        edgecolor=COLOR_NAVY,
        facecolor=COLOR_NAVY_LIGHT,
    )
    ax.add_patch(box)

    if branch_label:
        ax.text(
            x_center, y_center + height * 0.30,
            branch_label,
            ha="center", va="center",
            fontsize=7.5, fontweight="bold", color=COLOR_NAVY,
            style="italic",
        )
        title_y = y_center + height * 0.05
        annot_y = y_center - height * 0.28
    else:
        title_y = y_center + height * 0.18
        annot_y = y_center - height * 0.18

    ax.text(
        x_center, title_y,
        substep["title"],
        ha="center", va="center",
        fontsize=9, fontweight="bold", color=COLOR_NAVY,
    )
    ax.text(
        x_center, annot_y,
        substep["annotation"],
        ha="center", va="center",
        fontsize=6.8, color=TEXT_COLOR,
    )


def _draw_arrow(ax, xy_from, xy_to, color=COLOR_NAVY, lw=1.6,
                arrowstyle="-|>", linestyle="-"):
    arrow = FancyArrowPatch(
        xy_from, xy_to,
        arrowstyle=arrowstyle, mutation_scale=14,
        color=color, lw=lw,
        shrinkA=2, shrinkB=2,
        linestyle=linestyle,
    )
    ax.add_patch(arrow)


def panel_draw_pipeline(ax):
    """Render the 5-stage vertical pipeline with branched Stage 4."""
    # Canvas: 13 wide x 18 tall
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 18)
    ax.axis("off")

    center_x = 6.5
    box_w_ext = 7.5      # external-stage box width
    box_h_ext = 1.40     # external-stage box height

    # ----- Stage 1 (top) -----
    y_s1 = 17.0
    _draw_stage_box(ax, center_x, y_s1, box_w_ext, box_h_ext, EXTERNAL_STAGES[0],
                    is_primary=False)

    # ----- Stage 2 -----
    y_s2 = 15.2
    _draw_stage_box(ax, center_x, y_s2, box_w_ext, box_h_ext, EXTERNAL_STAGES[1],
                    is_primary=False)
    _draw_arrow(ax, (center_x, y_s1 - box_h_ext / 2),
                (center_x, y_s2 + box_h_ext / 2),
                color=COLOR_GRAY)

    # ----- Stage 3 -----
    y_s3 = 13.4
    _draw_stage_box(ax, center_x, y_s3, box_w_ext, box_h_ext, EXTERNAL_STAGES[2],
                    is_primary=False)
    _draw_arrow(ax, (center_x, y_s2 - box_h_ext / 2),
                (center_x, y_s3 + box_h_ext / 2),
                color=COLOR_GRAY)

    # ----- Stage 4 header (UPDD primary contribution) -----
    y_s4_header = 11.8
    header_box = FancyBboxPatch(
        (center_x - box_w_ext / 2, y_s4_header - 0.45),
        box_w_ext, 0.9,
        boxstyle="round,pad=0.02,rounding_size=0.10",
        linewidth=1.6, edgecolor=COLOR_NAVY,
        facecolor=COLOR_NAVY,
    )
    ax.add_patch(header_box)
    ax.text(
        center_x, y_s4_header + 0.10,
        "Stage 4 — Evaluation  |  UPDD primary contribution",
        ha="center", va="center",
        fontsize=11, fontweight="bold", color="white",
    )
    ax.text(
        center_x, y_s4_header - 0.20,
        "branched WT vs Variant ensembles  |  σ-decomposition  |  Convergence Index",
        ha="center", va="center",
        fontsize=8, color="white", style="italic",
    )
    _draw_arrow(ax, (center_x, y_s3 - box_h_ext / 2),
                (center_x, y_s4_header + 0.45),
                color=COLOR_GRAY)

    # ----- Stage 4 branched sub-steps -----
    # Branches:  Branch B (Variant) on left, Branch A (WT) on right
    # 4 sub-steps stacked vertically; ncAA param + MD + QM/MM are branched;
    # MM-PBSA is the convergence (single box).
    branch_left_x = 3.4
    branch_right_x = 9.6
    sub_w = 5.2
    sub_h = 1.45
    sub_h_branched = 1.55

    # ncAA parameterization (branched)
    y_sub1 = 10.0
    _draw_substep_box(ax, branch_left_x, y_sub1, sub_w, sub_h_branched,
                      STAGE4_SUBSTEPS[0], branch_label="Branch B (Variant)")
    _draw_substep_box(ax, branch_right_x, y_sub1, sub_w, sub_h_branched,
                      STAGE4_SUBSTEPS[0], branch_label="Branch A (WT, bit-identical control)")
    # Header arrow splits to both branches
    split_y = y_s4_header - 0.45
    _draw_arrow(ax, (center_x, split_y),
                (branch_left_x, y_sub1 + sub_h_branched / 2), color=COLOR_NAVY)
    _draw_arrow(ax, (center_x, split_y),
                (branch_right_x, y_sub1 + sub_h_branched / 2), color=COLOR_NAVY)

    # MD (branched)
    y_sub2 = 8.0
    _draw_substep_box(ax, branch_left_x, y_sub2, sub_w, sub_h_branched,
                      STAGE4_SUBSTEPS[1], branch_label="Branch B (Variant)")
    _draw_substep_box(ax, branch_right_x, y_sub2, sub_w, sub_h_branched,
                      STAGE4_SUBSTEPS[1], branch_label="Branch A (WT)")
    _draw_arrow(ax, (branch_left_x, y_sub1 - sub_h_branched / 2),
                (branch_left_x, y_sub2 + sub_h_branched / 2), color=COLOR_NAVY)
    _draw_arrow(ax, (branch_right_x, y_sub1 - sub_h_branched / 2),
                (branch_right_x, y_sub2 + sub_h_branched / 2), color=COLOR_NAVY)

    # QM/MM (branched)
    y_sub3 = 6.0
    _draw_substep_box(ax, branch_left_x, y_sub3, sub_w, sub_h_branched,
                      STAGE4_SUBSTEPS[2], branch_label="Branch B (Variant)")
    _draw_substep_box(ax, branch_right_x, y_sub3, sub_w, sub_h_branched,
                      STAGE4_SUBSTEPS[2], branch_label="Branch A (WT)")
    _draw_arrow(ax, (branch_left_x, y_sub2 - sub_h_branched / 2),
                (branch_left_x, y_sub3 + sub_h_branched / 2), color=COLOR_NAVY)
    _draw_arrow(ax, (branch_right_x, y_sub2 - sub_h_branched / 2),
                (branch_right_x, y_sub3 + sub_h_branched / 2), color=COLOR_NAVY)

    # MM-PBSA convergence (single, centered)
    y_sub4 = 4.0
    _draw_substep_box(ax, center_x, y_sub4, box_w_ext, sub_h, STAGE4_SUBSTEPS[3])
    # Branches converge into the centered MM-PBSA box
    _draw_arrow(ax, (branch_left_x, y_sub3 - sub_h_branched / 2),
                (center_x - 1.5, y_sub4 + sub_h / 2), color=COLOR_NAVY)
    _draw_arrow(ax, (branch_right_x, y_sub3 - sub_h_branched / 2),
                (center_x + 1.5, y_sub4 + sub_h / 2), color=COLOR_NAVY)

    # Convergence note between the two branches at Stage 4 bottom
    ax.text(
        center_x, y_sub4 + sub_h / 2 + 0.55,
        r"$\Delta\Delta G = \langle\langle \Delta G_{\mathrm{variant}} \rangle\rangle - \langle\langle \Delta G_{\mathrm{WT}} \rangle\rangle$",
        ha="center", va="bottom",
        fontsize=9.5, color=COLOR_NAVY, fontweight="bold",
    )

    # ----- Stage 5 -----
    y_s5 = 1.9
    _draw_stage_box(ax, center_x, y_s5, box_w_ext, box_h_ext, STAGE5,
                    is_primary=False)
    _draw_arrow(ax, (center_x, y_sub4 - sub_h / 2),
                (center_x, y_s5 + box_h_ext / 2),
                color=COLOR_GRAY)

    # Cycle-back arrow: Stage 5 -> Stage 1 (left side curving up). Endpoint
    # is shifted slightly INSIDE the Stage 1 box (x = left_edge + 0.25) so
    # the arrow head visibly penetrates the box rather than floating at the
    # edge; mutation_scale bumped 16 -> 20 for clearer head.
    cycle_arrow = FancyArrowPatch(
        (center_x - box_w_ext / 2, y_s5),
        (center_x - box_w_ext / 2 + 0.25, y_s1),
        connectionstyle="arc3,rad=-0.45",
        arrowstyle="-|>", mutation_scale=20,
        color=COLOR_GRAY, lw=1.5, linestyle=(0, (4, 2)),
    )
    ax.add_patch(cycle_arrow)
    ax.text(
        0.3, (y_s1 + y_s5) / 2,
        "iterative\nrefinement\ncycle",
        ha="center", va="center",
        fontsize=8.5, color=ANNOT_COLOR, style="italic",
        rotation=90,
    )

    # ----- Legend (upper-right) -----
    legend_handles = [
        mpatches.Patch(facecolor=COLOR_NAVY_LIGHT, edgecolor=COLOR_NAVY,
                       linewidth=1.2, label="UPDD primary contribution (Stage 4)"),
        mpatches.Patch(facecolor=COLOR_GRAY_LIGHT, edgecolor=COLOR_GRAY,
                       linewidth=1.2, label="External tool / wet-lab boundary"),
    ]
    ax.legend(
        handles=legend_handles, loc="upper right",
        bbox_to_anchor=(0.99, 0.995),
        fontsize=8.5, framealpha=0.92,
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Plot Paper 1 Figure 1 (5-stage pipeline architecture).",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("outputs/paper1/figures"),
        help="Directory to write fig1_pipeline.{png,svg,...}",
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

    fig, ax = plt.subplots(figsize=(11, 14), dpi=args.dpi)
    panel_draw_pipeline(ax)

    fig.suptitle(
        "Figure 1. UPDD 5-stage pipeline — Stage 4 branched evaluation core",
        fontsize=13, fontweight="bold", y=0.985,
    )

    saved = []
    for fmt in args.format:
        out = args.output_dir / f"fig1_pipeline.{fmt}"
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
