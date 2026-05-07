"""scripts/plot_fig3_detachment_events.py - Paper 1 Figure 3 plotting.

Three-row panel composition of detachment event trajectories. Each row
shows one of the three biophysical detachment events identified in the
2026-04-29 detachment investigation audit:

  (a) 3IOL_NML20_s101 - sustained dissociation
      (fraction = 1.000, 25/25 frames)
  (b) 1YCR_NML22_s19  - transient bimodal
      (fraction = 0.760 snapshot / 0.546 trajectory)
  (c) 1YCR_NML22_s42  - sustained surface-bound
      (fraction = 1.000, 25/25 frames)

Each row is a distance-vs-time trace (left, ~70% width) paired with the
per-frame distance histogram (right, ~30% width). Distributions are
unimodal at 7-13 A (not bimodal at ~75 + ~1.3 A), distinguishing the
biophysical detachment signal from periodic-image (PBC) artifact.

The companion geometry and contact-count panels for each system live in
the same source directory and are reserved for the supplementary figure.

Source data: outputs/analysis/detachment_plots_20260429/ (12 individual
PNGs from the 2026-04-29 detachment investigation audit, anchored to
pathology/detachment_investigation_20260429.md).

Anchor data SSOT lock notes (verified 2026-05-05):
  - Source PNGs were generated 2026-04-29 22:45, AFTER all relevant
    Phase beta patches: L387 v54 (12:59), Phase beta re-extraction
    sweep of 28 systems (16:17), detachment metric introduction (18:47),
    and CONECT v55 wraparound fix (22:31). Source trajectories are
    therefore post-L387-v54 and post-CONECT-v55 clean.
  - Detachment fractions hardcoded as caption strings exactly mirror
    verdict section 3.3.3 + section 3.3.2 footnote 2 reconciled values
    (1YCR_NML22_s19 dual-resolution: 0.546 trajectory frame-by-frame /
    0.760 snapshot K-Means oversampled).
  - Figure 3 scope deliberately excludes 7TL8_MTR6 (which is the s55
    intra-residue PBC defect audit case visualized via figure 1 panel
    (c) sigma_w refresh); that scope split avoids re-using the same
    system across two figures and keeps the PR-21 intra-residue bond
    audit narrative cleanly attributed to figure 1.

Outline reference: outline section 6, lines 678-689.

Usage:
    python scripts/plot_fig3_detachment_events.py
    python scripts/plot_fig3_detachment_events.py --output-dir outputs/paper1/figures
    python scripts/plot_fig3_detachment_events.py --format png svg pdf
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["axes.unicode_minus"] = False
matplotlib.rcParams["mathtext.fontset"] = "dejavusans"

# Top crop fraction applied to source PNGs to strip baked-in matplotlib titles
# (the source detachment_plots_20260429/*.png have their own ax.set_title that
# duplicate the section header rendered by this script; cropping ~10% from the
# top removes that duplication without losing data content).
_TOP_CROP_FRAC = 0.10


# ---------------------------------------------------------------------------
# Source layout
# ---------------------------------------------------------------------------
SOURCE_DIR = Path("outputs/analysis/detachment_plots_20260429")

SYSTEMS = [
    (
        "3IOL_NML20_s101",
        "(a)",
        "Sustained dissociation (fraction = 1.000, 25/25 frames)",
    ),
    (
        "1YCR_NML22_s19",
        "(b)",
        "Transient bimodal (fraction = 0.760 snapshot / 0.546 trajectory)",
    ),
    (
        "1YCR_NML22_s42",
        "(c)",
        "Sustained surface-bound (fraction = 1.000, 25/25 frames)",
    ),
]


# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------
def _crop_top(img_array, frac=_TOP_CROP_FRAC):
    """Crop the top `frac` fraction of an image array to strip baked-in titles."""
    h = img_array.shape[0]
    return img_array[int(h * frac):, :]


def _draw_row(fig, gs, row_idx, system, panel_label, caption, source_dir):
    """Render one system row: distance-vs-time (~70%) + histogram (~30%).

    Source PNGs carry baked-in matplotlib titles (e.g., "{system}: peptide-bond
    proxy distance vs time") that duplicate the section header set on ax_left.
    We strip those by cropping the top of each source image. The right
    histogram subplot has no matplotlib-level title (the section header on the
    left subplot serves the whole row).
    """
    src_left = source_dir / f"{system}_distance_vs_time.png"
    src_right = source_dir / f"{system}_histogram.png"

    if not src_left.is_file():
        raise FileNotFoundError(f"distance_vs_time PNG missing: {src_left}")
    if not src_right.is_file():
        raise FileNotFoundError(f"histogram PNG missing: {src_right}")

    ax_left = fig.add_subplot(gs[row_idx, :6])
    ax_left.imshow(_crop_top(mpimg.imread(src_left)))
    ax_left.set_axis_off()
    ax_left.set_title(
        f"{panel_label} {system} — {caption}",
        fontsize=10.5,
        loc="left",
        pad=4,
    )

    ax_right = fig.add_subplot(gs[row_idx, 6:])
    ax_right.imshow(_crop_top(mpimg.imread(src_right)))
    ax_right.set_axis_off()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Plot Paper 1 Figure 3 (detachment event trajectories).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/paper1/figures"),
        help="Directory to write fig3_detachment_events.{png,svg,...}",
    )
    parser.add_argument(
        "--format",
        nargs="+",
        default=["png", "svg"],
        help="Output formats (e.g. png svg pdf).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="DPI for raster outputs (default: 300).",
    )
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=SOURCE_DIR,
        help="Override the source PNG directory (default: %(default)s).",
    )
    args = parser.parse_args(argv)

    source_dir = args.source_dir
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(13, 9.8), dpi=args.dpi)
    gs = fig.add_gridspec(
        nrows=3,
        ncols=10,
        hspace=0.20,
        wspace=0.18,
        left=0.03,
        right=0.985,
        top=0.93,
        bottom=0.04,
    )

    for row_idx, (system, panel_label, caption) in enumerate(SYSTEMS):
        _draw_row(fig, gs, row_idx, system, panel_label, caption, source_dir)

    fig.suptitle(
        "Figure 3. Per-frame detachment event trajectories: "
        "unimodal at 7-13 A distinguishes biophysical signal from PBC artifact",
        fontsize=13,
        fontweight="bold",
        y=0.985,
    )

    saved = []
    for fmt in args.format:
        out = args.output_dir / f"fig3_detachment_events.{fmt}"
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
