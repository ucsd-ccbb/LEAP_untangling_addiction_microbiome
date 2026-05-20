#!/usr/bin/env python3

import argparse
import sys

import pandas as pd
import matplotlib.pyplot as plt
from qiime2 import Artifact
from skbio.stats.ordination import OrdinationResults
from scipy.spatial import ConvexHull


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot a QIIME 2 PCoA (.qza) with group-colored points and convex hulls."
    )

    parser.add_argument(
        "--pcoa", required=True, help="Path to QIIME 2 PCoA artifact (.qza)"
    )
    parser.add_argument(
        "--metadata", required=True, help="Path to tab-delimited metadata file"
    )
    parser.add_argument(
        "--group-col", required=True, help="Metadata column used to color/group samples"
    )
    parser.add_argument(
        "--sample-id-col",
        default="SampleID",
        help="Metadata column containing sample IDs (default: SampleID)",
    )

    parser.add_argument(
        "--output",
        default="pcoa_plot.png",
        help="Output image path (default: pcoa_plot.png)",
    )
    parser.add_argument(
        "--title", default="Bray-curtis", help="Plot title (default: Bray-curtis)"
    )

    parser.add_argument(
        "--width", type=float, default=7.0, help="Figure width in inches (default: 7)"
    )
    parser.add_argument(
        "--height", type=float, default=7.0, help="Figure height in inches (default: 7)"
    )

    parser.add_argument(
        "--point-size", type=float, default=60.0, help="Point size (default: 60)"
    )
    parser.add_argument(
        "--alpha", type=float, default=0.18, help="Hull fill alpha (default: 0.18)"
    )
    parser.add_argument(
        "--dpi", type=int, default=300, help="Output DPI (default: 300)"
    )

    parser.add_argument("--no-legend", action="store_true", help="Suppress legend")

    parser.add_argument(
        "--group-color",
        action="append",
        default=[],
        help=(
            "Set a specific color for a group using GROUP=COLOR. Repeat for multiple groups. "
            "Example: --group-color Control=#1f77b4"
        ),
    )

    # Font sizes
    parser.add_argument(
        "--title-size", type=float, default=22.0, help="Title font size (default: 22)"
    )
    parser.add_argument(
        "--axis-label-size",
        type=float,
        default=18.0,
        help="Axis label font size (default: 18)",
    )
    parser.add_argument(
        "--tick-label-size",
        type=float,
        default=14.0,
        help="Tick label font size (default: 14)",
    )
    parser.add_argument(
        "--legend-text-size",
        type=float,
        default=14.0,
        help="Legend text size (default: 14)",
    )

    # Legend spacing controls (tighten marker/text gap + push legend closer to axes)
    parser.add_argument(
        "--legend-handletextpad",
        type=float,
        default=0.3,
        help="Gap between marker and legend text (default: 0.3)",
    )
    parser.add_argument(
        "--legend-borderaxespad",
        type=float,
        default=0.2,
        help="Padding between legend and axes (default: 0.2)",
    )
    parser.add_argument(
        "--legend-borderpad",
        type=float,
        default=0.2,
        help="Padding inside legend frame (default: 0.2)",
    )
    parser.add_argument(
        "--legend-labelspacing",
        type=float,
        default=0.3,
        help="Vertical space between legend entries (default: 0.3)",
    )
    parser.add_argument(
        "--legend-handlelength",
        type=float,
        default=1.0,
        help="Length of legend handles (default: 1.0)",
    )

    return parser.parse_args()


def parse_group_colors(group_color_args):
    color_map = {}
    for item in group_color_args:
        if "=" not in item:
            raise ValueError(
                f"Invalid --group-color value '{item}'. Expected format: GROUP=COLOR"
            )
        group, color = item.split("=", 1)
        group = group.strip()
        color = color.strip()
        if not group:
            raise ValueError(f"Invalid group name in --group-color '{item}'")
        if not color:
            raise ValueError(f"Invalid color in --group-color '{item}'")
        color_map[group] = color
    return color_map


def load_pcoa(pcoa_path: str):
    art = Artifact.load(pcoa_path)
    ordination = art.view(OrdinationResults)

    coords = ordination.samples.copy()
    coords.index.name = "SampleID"
    coords = coords.reset_index()

    ve = ordination.proportion_explained.values * 100

    if coords.shape[1] < 3:
        raise ValueError("PCoA coordinates do not contain at least two axes.")

    xcol = coords.columns[1]
    ycol = coords.columns[2]

    return coords, ve, xcol, ycol


def load_metadata(metadata_path: str, sample_id_col: str):
    # NOTE: comment="#" is convenient for plain TSVs with comment lines,
    # but it will *skip* a QIIME-style header like "#SampleID".
    # If your file is true QIIME2 metadata, consider renaming "#SampleID" to "SampleID"
    # or remove comment parsing.
    meta = pd.read_csv(metadata_path, sep="\t", comment="#")

    if sample_id_col not in meta.columns:
        raise ValueError(
            f"Metadata file must contain the sample ID column '{sample_id_col}'. "
            f"Available columns: {', '.join(meta.columns)}"
        )
    return meta


def draw_group_hull(ax, subdf, xcol, ycol, color, alpha):
    if len(subdf) < 3:
        return

    pts = subdf[[xcol, ycol]].dropna().to_numpy()
    if pts.shape[0] < 3:
        return

    try:
        hull = ConvexHull(pts)
        hull_pts = pts[hull.vertices]

        hull_x = list(hull_pts[:, 0]) + [hull_pts[0, 0]]
        hull_y = list(hull_pts[:, 1]) + [hull_pts[0, 1]]

        ax.plot(hull_x, hull_y, linewidth=1.2, color=color)
        ax.fill(hull_x, hull_y, color=color, alpha=alpha)
    except Exception:
        # Hull can fail if points are collinear, etc.
        pass


def main():
    args = parse_args()

    try:
        coords, ve, xcol, ycol = load_pcoa(args.pcoa)
        meta = load_metadata(args.metadata, args.sample_id_col)
        group_colors = parse_group_colors(args.group_color)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    if args.group_col not in meta.columns:
        print(
            f"ERROR: Group column '{args.group_col}' not found in metadata.",
            file=sys.stderr,
        )
        print(f"Metadata columns: {', '.join(meta.columns)}", file=sys.stderr)
        sys.exit(1)

    # Normalize metadata sample ID column to SampleID for merging
    meta = meta.rename(columns={args.sample_id_col: "SampleID"})

    df = coords.merge(meta, on="SampleID", how="left")

    if df[args.group_col].isna().all():
        print(
            "ERROR: After merging, all group values are missing. "
            "Check that sample IDs match between the PCoA artifact and metadata.",
            file=sys.stderr,
        )
        sys.exit(1)

    fig, ax = plt.subplots(figsize=(args.width, args.height))

    groups = [g for g in df[args.group_col].dropna().unique()]

    for group in groups:
        sub = df[df[args.group_col] == group].copy()
        group_str = str(group)

        if group_str in group_colors:
            color = group_colors[group_str]
            ax.scatter(
                sub[xcol], sub[ycol], s=args.point_size, label=group_str, color=color
            )
            draw_group_hull(ax, sub, xcol, ycol, color, args.alpha)
        else:
            sc = ax.scatter(sub[xcol], sub[ycol], s=args.point_size, label=group_str)
            auto_color = sc.get_facecolor()[0]
            draw_group_hull(ax, sub, xcol, ycol, auto_color, args.alpha)

    # Reference lines at zero
    ax.axhline(0, linestyle=(0, (5, 5)), linewidth=1.2, color="0.3")
    ax.axvline(0, linestyle=(0, (5, 5)), linewidth=1.2, color="0.3")

    # Axis labels with variance explained
    x_label = f"PCoA1 ({ve[0]:.1f}%)" if len(ve) > 0 else "PCoA1"
    y_label = f"PCoA2 ({ve[1]:.1f}%)" if len(ve) > 1 else "PCoA2"

    ax.set_xlabel(x_label, fontsize=args.axis_label_size)
    ax.set_ylabel(y_label, fontsize=args.axis_label_size)
    ax.set_title(args.title, loc="left", fontsize=args.title_size, fontweight="bold")

    # Tick sizes
    ax.tick_params(labelsize=args.tick_label_size, width=2, length=8)

    # Thicker spines
    for spine in ax.spines.values():
        spine.set_linewidth(2)

    if not args.no_legend:
        ax.legend(
            frameon=False,
            fontsize=args.legend_text_size,
            handletextpad=args.legend_handletextpad,
            borderaxespad=args.legend_borderaxespad,
            borderpad=args.legend_borderpad,
            labelspacing=args.legend_labelspacing,
            handlelength=args.legend_handlelength,
        )

    plt.tight_layout()
    plt.savefig(args.output, dpi=args.dpi)
    plt.close()

    print(f"Saved plot to: {args.output}")


if __name__ == "__main__":
    main()
