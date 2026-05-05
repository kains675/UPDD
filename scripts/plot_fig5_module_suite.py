"""scripts/plot_fig5_module_suite.py - Paper 1 Figure 5 plotting.

Generates the Capability Level 1 module suite architecture diagram:
  - 5 modules arranged left-to-right (Sequence Parser, Branched ddG Engine,
    Iterative Cycle Manager, User CLI + Stage 1-2 wrappers, HTML Report Generator)
  - User CLI orchestration bar at the top dispatching to all 5 modules
  - Schema versions and key dataclasses annotated as edge labels
  - Test summary banner below

This is a vector architecture diagram (not a data plot).

Usage:
    python scripts/plot_fig5_module_suite.py
    python scripts/plot_fig5_module_suite.py --output-dir /tmp/fig_test
    python scripts/plot_fig5_module_suite.py --format png svg pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

# ---------------------------------------------------------------------------
# Style palette
# ---------------------------------------------------------------------------
COLOR_NAVY = "#1f4e79"          # UPDD-primary (Capability Level 1 modules)
COLOR_NAVY_LIGHT = "#dbe7f3"    # Light navy fill
COLOR_GRAY = "#aaaaaa"          # external / boundary
COLOR_GRAY_LIGHT = "#e8e8e8"    # light gray fill
EDGE_COLOR = "#222222"
TEXT_COLOR = "#1a1a1a"
ANNOT_COLOR = "#444444"

# ---------------------------------------------------------------------------
# Module definitions
# ---------------------------------------------------------------------------
MODULES = [
    {
        "title": "Module 1\nSequence Parser",
        "schema": "sequence_parser/0.1",
        "details": "ncAA registry\nvalidation\ncanonical design\nspecification",
    },
    {
        "title": "Module 2\nBranched ΔΔG\n+ Light D Hybrid",
        "schema": "branched_ddg/0.3",
        "details": "paired (variant, WT)\nbit-identical WT\ndetachment metric\nadvisory",
    },
    {
        "title": "Module 3\nIterative Cycle\nManager",
        "schema": "cycle_manager/0.2",
        "details": "Cycle 0 / Cycle 1+\nmulti-assay K_d\naggregation\nplateau detection",
    },
    {
        "title": "Module 4\nStage 1-2\nwrappers",
        "schema": "stage12_pipeline/0.1",
        "details": "RFdiffusion\nProteinMPNN\nAF2 / ColabFold\nmock CPU fallback",
    },
    {
        "title": "Module 5\nHTML Report\nGenerator",
        "schema": "html_report/0.2",
        "details": "cycle progression\ndetachment banner\nprint CSS",
    },
]

# Edge labels between adjacent modules (data contracts)
EDGE_LABELS = [
    "CandidateOutput-\nlike",      # M1 -> M2
    "BranchedDDG-\nResult",        # M2 -> M3
    "CycleEntry  (cycle_n + WetLabResult[] + PredictedDDG[])",  # M3 -> M5  (skips M4 visually)
]


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------
def _draw_module_box(ax, x_center, y_center, width, height, module, fontsize_title=8.5):
    """Draw a single rounded module box at the specified center coordinates."""
    box = FancyBboxPatch(
        (x_center - width / 2, y_center - height / 2),
        width, height,
        boxstyle="round,pad=0.02,rounding_size=0.10",
        linewidth=1.4,
        edgecolor=COLOR_NAVY,
        facecolor=COLOR_NAVY_LIGHT,
    )
    ax.add_patch(box)

    # Title (bold, navy)
    ax.text(
        x_center, y_center + height * 0.30,
        module["title"],
        ha="center", va="center",
        fontsize=fontsize_title, fontweight="bold",
        color=COLOR_NAVY,
    )
    # Schema (italic monospace)
    ax.text(
        x_center, y_center - height * 0.02,
        module["schema"],
        ha="center", va="center",
        fontsize=6.0, style="italic", family="monospace",
        color=ANNOT_COLOR,
    )
    # Details
    ax.text(
        x_center, y_center - height * 0.36,
        module["details"],
        ha="center", va="center",
        fontsize=6.2, color=TEXT_COLOR,
    )


def _draw_arrow(ax, xy_from, xy_to, color=COLOR_NAVY, lw=1.6, label=None,
                label_offset=(0, 0.32), label_fontsize=7):
    """Draw a directional arrow with optional edge label centered above the arrow."""
    arrow = FancyArrowPatch(
        xy_from, xy_to,
        arrowstyle="-|>", mutation_scale=14,
        color=color, lw=lw,
        shrinkA=2, shrinkB=2,
    )
    ax.add_patch(arrow)
    if label:
        mid_x = (xy_from[0] + xy_to[0]) / 2 + label_offset[0]
        mid_y = (xy_from[1] + xy_to[1]) / 2 + label_offset[1]
        ax.text(
            mid_x, mid_y, label,
            ha="center", va="bottom",
            fontsize=label_fontsize, color=ANNOT_COLOR, style="italic",
        )


def panel_draw_modules(ax, modules):
    """Render the 5-module flow + CLI orchestration bar + test summary."""
    # Canvas: 13 x 8
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 8)
    ax.axis("off")

    # ----- Module row geometry -----
    n = len(modules)
    box_w = 1.75
    box_h = 1.85
    row_y = 3.5
    # Centers: distribute evenly between x=1.4 and x=11.6
    # Inter-box gap = ((x_right - x_left)/(n-1)) - box_w  ~  0.80, sufficient for labels
    x_left, x_right = 1.4, 11.6
    centers = [x_left + i * (x_right - x_left) / (n - 1) for i in range(n)]

    # ----- CLI orchestration bar (top) -----
    cli_y = 6.4
    cli_h = 0.85
    cli_box = FancyBboxPatch(
        (0.6, cli_y - cli_h / 2),
        12.0, cli_h,
        boxstyle="round,pad=0.02,rounding_size=0.18",
        linewidth=1.6,
        edgecolor=COLOR_NAVY,
        facecolor=COLOR_NAVY,
    )
    ax.add_patch(cli_box)
    ax.text(
        6.5, cli_y + 0.15,
        "User CLI orchestrator (updd_cli.py)",
        ha="center", va="center",
        fontsize=12, fontweight="bold", color="white",
    )
    ax.text(
        6.5, cli_y - 0.20,
        "subcommands:  check  /  compare  /  refine  /  discover  /  cycle  /  report",
        ha="center", va="center",
        fontsize=8.5, color="white", style="italic",
    )

    # ----- Module boxes -----
    for cx, mod in zip(centers, modules):
        _draw_module_box(ax, cx, row_y, box_w, box_h, mod)

    # ----- CLI orchestration arrows (top bar -> each module top) -----
    cli_bottom_y = cli_y - cli_h / 2
    for cx in centers:
        arrow = FancyArrowPatch(
            (cx, cli_bottom_y),
            (cx, row_y + box_h / 2),
            arrowstyle="-|>", mutation_scale=10,
            color=COLOR_NAVY, lw=1.0, alpha=0.7,
            linestyle=(0, (3, 2)),
        )
        ax.add_patch(arrow)

    # Annotation for orchestration arrows
    ax.text(
        0.55, (cli_bottom_y + row_y + box_h / 2) / 2,
        "dispatch",
        ha="left", va="center",
        fontsize=7.5, color=COLOR_NAVY, style="italic", rotation=90,
    )

    # ----- Module-to-module data-contract arrows -----
    # M1 -> M2
    _draw_arrow(
        ax,
        (centers[0] + box_w / 2, row_y),
        (centers[1] - box_w / 2, row_y),
        label=EDGE_LABELS[0], label_offset=(0, 0.05), label_fontsize=6.2,
    )
    # M2 -> M3
    _draw_arrow(
        ax,
        (centers[1] + box_w / 2, row_y),
        (centers[2] - box_w / 2, row_y),
        label=EDGE_LABELS[1], label_offset=(0, 0.05), label_fontsize=6.2,
    )
    # M3 -> M5 (skip the wrappers visually, curving over Module 4)
    skip_arrow = FancyArrowPatch(
        (centers[2] + box_w / 2 - 0.10, row_y + box_h / 2 - 0.20),
        (centers[4] - box_w / 2 + 0.10, row_y + box_h / 2 - 0.20),
        connectionstyle="arc3,rad=-0.40",
        arrowstyle="-|>", mutation_scale=14,
        color=COLOR_NAVY, lw=1.6,
    )
    ax.add_patch(skip_arrow)
    # Place label above the curve apex
    ax.text(
        (centers[2] + centers[4]) / 2, row_y + box_h / 2 + 1.20,
        EDGE_LABELS[2],
        ha="center", va="bottom",
        fontsize=6.5, color=ANNOT_COLOR, style="italic",
    )

    # Side note: Module 4 also produces wrapper artifacts feeding Module 2 (mock fallback)
    ax.text(
        centers[3], row_y - box_h / 2 - 0.45,
        "(Stage 1-2 wrappers feed structural inputs; mock fallback for CPU-only environments)",
        ha="center", va="top",
        fontsize=7, color=ANNOT_COLOR, style="italic",
    )

    # ----- Test summary banner (bottom) -----
    test_y = 1.05
    test_box = FancyBboxPatch(
        (1.0, test_y - 0.45),
        11.0, 0.9,
        boxstyle="round,pad=0.02,rounding_size=0.10",
        linewidth=1.0,
        edgecolor=EDGE_COLOR,
        facecolor="#fafafa",
    )
    ax.add_patch(test_box)
    ax.text(
        6.5, test_y + 0.12,
        "Test coverage: 284 unit + 10 simulated AL cycle = 294 pass / 9 skip / 1 pre-existing fail",
        ha="center", va="center",
        fontsize=10, fontweight="bold", color=TEXT_COLOR,
    )
    ax.text(
        6.5, test_y - 0.20,
        "(pre-existing fail: test_admet_filter::test_lipinski_constants_are_module_level — separate semantic triage)",
        ha="center", va="center",
        fontsize=8, color=ANNOT_COLOR, style="italic",
    )

    # ----- Legend (upper-right) -----
    legend_handles = [
        mpatches.Patch(facecolor=COLOR_NAVY_LIGHT, edgecolor=COLOR_NAVY,
                       linewidth=1.2, label="Capability Level 1 module"),
        mpatches.Patch(facecolor=COLOR_NAVY, edgecolor=COLOR_NAVY,
                       linewidth=1.2, label="CLI orchestrator"),
    ]
    ax.legend(
        handles=legend_handles, loc="upper right",
        bbox_to_anchor=(0.99, 0.99),
        fontsize=8, framealpha=0.9,
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Plot Paper 1 Figure 5 (Capability Level 1 module suite).",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=Path("outputs/paper1/figures"),
        help="Directory to write fig5_module_suite.{png,svg,...}",
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

    fig, ax = plt.subplots(figsize=(13, 8), dpi=args.dpi)
    panel_draw_modules(ax, MODULES)

    fig.suptitle(
        "Figure 5. Capability Level 1 module suite — 5-module architecture",
        fontsize=13, fontweight="bold", y=0.97,
    )

    saved = []
    for fmt in args.format:
        out = args.output_dir / f"fig5_module_suite.{fmt}"
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
