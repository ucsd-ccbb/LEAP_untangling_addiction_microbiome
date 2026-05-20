#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch


def load_abundance_table(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    return df


def load_metadata(path, sample_col=None):
    meta = pd.read_csv(path, sep="\t")
    if sample_col is None:
        sample_col = meta.columns[0]
    if sample_col not in meta.columns:
        raise ValueError(f"Sample column '{sample_col}' not found in metadata.")
    return meta.set_index(sample_col)


def relative_frequencies(count_df):
    row_sums = count_df.sum(axis=1)
    return count_df.div(row_sums.replace(0, np.nan), axis=0).fillna(0)


def short_taxon_name(name):
    parts = [p.strip() for p in str(name).split(";") if p.strip()]
    return parts[-1] if parts else str(name).strip()


def italicize_label(label):
    if label == "Other":
        return label
    escaped = label.replace(" ", r"\ ")
    return rf"$\it{{{escaped}}}$"


def desaturate_color(color, amount=0.45, lighten=0.10):
    rgb = np.array(mcolors.to_rgb(color))
    gray = np.array([rgb.mean(), rgb.mean(), rgb.mean()])
    muted = gray * (1 - amount) + rgb * amount
    muted = muted * (1 - lighten) + np.ones(3) * lighten
    return tuple(np.clip(muted, 0, 1))


def build_muted_palette(n):
    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    return [
        desaturate_color(base[i % len(base)], amount=0.42, lighten=0.08)
        for i in range(n)
    ]


def choose_pivot_taxon(row, taxa, threshold):
    for taxon in taxa:
        if row[taxon] >= threshold:
            return taxon
    return taxa[0]


def build_sample_sort_key(row, taxa, threshold):
    values = [row[t] for t in taxa]

    pivot_idx = None
    for i, v in enumerate(values):
        if v >= threshold:
            pivot_idx = i
            break
    if pivot_idx is None:
        pivot_idx = 0

    rotated = values[pivot_idx:] + values[:pivot_idx]
    return tuple([pivot_idx] + [-v for v in rotated] + [str(row.name)])


def prepare_sample_level_data(
    rel_df, meta_df, group_col, top_n, group_order=None, threshold=0.05
):
    if group_col not in meta_df.columns:
        raise ValueError(f"Grouping column '{group_col}' not found in metadata.")

    shared = rel_df.index.intersection(meta_df.index)
    if len(shared) == 0:
        raise ValueError(
            "No overlapping sample IDs between abundance table and metadata."
        )

    rel_df = rel_df.loc[shared].copy()
    meta_df = meta_df.loc[shared].copy()

    global_means = rel_df.mean(axis=0).sort_values(ascending=False)
    top_taxa = list(global_means.head(top_n).index)

    plot_df = rel_df[top_taxa].copy()
    other_cols = [c for c in rel_df.columns if c not in top_taxa]
    if other_cols:
        plot_df["Other"] = rel_df[other_cols].sum(axis=1)

    plot_df[group_col] = meta_df[group_col].astype(str)

    if group_order is None:
        ordered_groups = list(dict.fromkeys(plot_df[group_col].tolist()))
    else:
        ordered_groups = group_order
        missing = [g for g in ordered_groups if g not in set(plot_df[group_col])]
        if missing:
            raise ValueError(f"These groups were not found: {missing}")

    ordered_taxa = top_taxa.copy()
    if "Other" in plot_df.columns:
        ordered_taxa.append("Other")

    ranking_taxa = [t for t in ordered_taxa if t != "Other"]

    plot_df["_pivot_taxon"] = plot_df.apply(
        lambda row: choose_pivot_taxon(row, ranking_taxa, threshold), axis=1
    )

    pieces = []
    for g in ordered_groups:
        sub = plot_df.loc[plot_df[group_col] == g].copy()
        if sub.empty:
            continue

        sub["_sort_key"] = sub.apply(
            lambda row: build_sample_sort_key(row, ranking_taxa, threshold), axis=1
        )
        sub = sub.sort_values("_sort_key").drop(columns=["_sort_key"])
        pieces.append(sub)

    plot_df = pd.concat(pieces, axis=0)
    return plot_df, ordered_taxa, ordered_groups


def build_x_positions(group_values, ordered_groups, within_step=0.32, group_gap=0.90):
    x_positions = []
    centers = []
    boundaries = []

    x = 0.0
    idx = 0
    n = len(group_values)

    for g in ordered_groups:
        group_idx = [i for i, val in enumerate(group_values) if val == g]
        if not group_idx:
            continue

        start_x = x
        for _ in group_idx:
            x_positions.append(x)
            x += within_step

        end_x = x - within_step
        centers.append(((start_x + end_x) / 2.0, g))

        if idx + len(group_idx) < n:
            boundaries.append(end_x + (within_step / 2.0) + (group_gap / 2.0))
            x += group_gap

        idx += len(group_idx)

    return np.array(x_positions), centers, boundaries


def get_stack_order_for_sample(row, taxa):
    pivot = row["_pivot_taxon"]

    non_other = [t for t in taxa if t != "Other"]
    remaining = [t for t in non_other if t != pivot]

    global_rank = {t: i for i, t in enumerate(non_other)}
    remaining_sorted = sorted(remaining, key=lambda t: (-row[t], global_rank[t]))

    stack_order = [pivot] + remaining_sorted
    if "Other" in taxa:
        stack_order.append("Other")
    return stack_order


def apply_fontsizes(
    title_fontsize, axis_label_fontsize, tick_label_fontsize, legend_fontsize
):
    plt.rcParams.update(
        {
            "axes.titlesize": title_fontsize,
            "axes.labelsize": axis_label_fontsize,
            "xtick.labelsize": tick_label_fontsize,
            "ytick.labelsize": tick_label_fontsize,
            "legend.fontsize": legend_fontsize,
        }
    )


def make_plot(
    plot_df,
    taxa,
    group_col,
    ordered_groups,
    output_file,
    title=None,
    bar_width=0.18,
    legend_cols=1,
    show_sample_labels=False,
    within_step=0.32,
    group_gap=0.90,
    fig_height=None,
    title_fontsize=16,
    axis_label_fontsize=14,
    tick_label_fontsize=11,
    legend_fontsize=12,
    group_label_fontsize=13,
    group_label_y=-3.0,
):
    apply_fontsizes(
        title_fontsize, axis_label_fontsize, tick_label_fontsize, legend_fontsize
    )

    n_samples = plot_df.shape[0]
    group_values = plot_df[group_col].tolist()

    x, centers, boundaries = build_x_positions(
        group_values, ordered_groups, within_step=within_step, group_gap=group_gap
    )

    palette = build_muted_palette(len(taxa))
    colors = {taxon: palette[i] for i, taxon in enumerate(taxa)}
    if "Other" in colors:
        colors["Other"] = (0.82, 0.82, 0.82)

    display_labels = {}
    for taxon in taxa:
        if taxon == "Other":
            display_labels[taxon] = "Other"
        else:
            display_labels[taxon] = italicize_label(short_taxon_name(taxon))

    fig_width = max(8, n_samples * 0.18 + len(ordered_groups) * 0.6)

    if fig_height is None:
        default_height = 12
        # Estimate legend height in inches
        line_height_in = (
            legend_fontsize / 72.0
        ) * 1.45  # font pt → inches with spacing
        legend_height_in = len(taxa) * line_height_in + 0.5
        # The axes occupy only a fraction of the figure (bottom=0.12, top≈0.90)
        axes_fraction = 0.78
        fig_height = max(default_height, legend_height_in / axes_fraction)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    for i, (_, row) in enumerate(plot_df.iterrows()):
        xpos = x[i]
        bottom = 0.0
        stack_order = get_stack_order_for_sample(row, taxa)

        for taxon in stack_order:
            val = row[taxon] * 100
            if val <= 0:
                continue
            ax.bar(
                xpos,
                val,
                width=bar_width,
                bottom=bottom,
                color=colors[taxon],
                edgecolor="#444444",
                linewidth=0.28,
                align="center",
            )
            bottom += val

    for b in boundaries:
        ax.axvline(b, color="#d0d0d0", linewidth=1.4)

    ax.set_xticks(x)
    if show_sample_labels:
        ax.set_xticklabels(plot_df.index.tolist(), rotation=90)
    else:
        ax.set_xticklabels([""] * n_samples)
        ax.tick_params(axis="x", length=0)

    for center, label in centers:
        ax.text(
            center,
            group_label_y,
            label,
            ha="center",
            va="top",
            rotation=45,
            fontsize=group_label_fontsize,
            clip_on=False,
        )

    ax.set_ylabel("Relative abundance (%)")
    ax.set_ylim(0, 100)

    if len(x) > 0:
        ax.set_xlim(x.min() - bar_width * 0.7, x.max() + bar_width * 0.7)

    if title:
        ax.set_title(title, pad=12)

    legend_handles = [
        Patch(facecolor=colors[t], edgecolor="#444444", linewidth=0.28) for t in taxa
    ]
    ax.legend(
        handles=legend_handles,
        labels=[display_labels[t] for t in taxa],
        ncol=legend_cols,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        frameon=False,
        columnspacing=0.8,
        handlelength=1.0,
        handletextpad=0.4,
        fontsize=legend_fontsize,
    )

    plt.subplots_adjust(bottom=0.12, right=0.78)
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Sample-level stacked relative abundance plot grouped by metadata category (per-sample taxa ordering)."
    )
    parser.add_argument("--abundance_tsv", required=True)
    parser.add_argument("--metadata_tsv", required=True)
    parser.add_argument("--group_col", required=True)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--top_n", type=int, default=10)
    parser.add_argument("--output", required=True)
    parser.add_argument("--title", default=None)
    parser.add_argument(
        "--group_order", default=None, help="Comma-separated group order"
    )
    parser.add_argument("--show_sample_labels", action="store_true")
    parser.add_argument(
        "--switch_threshold",
        type=float,
        default=0.05,
        help="Threshold below which sample sorting falls through to the next taxon",
    )

    # Font sizes (requested naming style)
    parser.add_argument("--title_fontsize", type=float, default=16)
    parser.add_argument("--axis_label_fontsize", type=float, default=14)
    parser.add_argument("--tick_label_fontsize", type=float, default=11)
    parser.add_argument("--legend_fontsize", type=float, default=12)
    parser.add_argument("--group_label_fontsize", type=float, default=13)

    # Geometry / spacing knobs
    parser.add_argument("--fig_height", type=float, default=None)
    parser.add_argument("--bar_width", type=float, default=0.18)
    parser.add_argument("--within_step", type=float, default=0.32)
    parser.add_argument("--group_gap", type=float, default=0.90)
    parser.add_argument("--group_label_y", type=float, default=-3.0)

    args = parser.parse_args()

    abundance_df = load_abundance_table(args.abundance_tsv)
    metadata_df = load_metadata(args.metadata_tsv, sample_col=args.sample_col)
    rel_df = relative_frequencies(abundance_df)

    group_order = None
    if args.group_order:
        group_order = [x.strip() for x in args.group_order.split(",")]

    plot_df, taxa, ordered_groups = prepare_sample_level_data(
        rel_df,
        metadata_df,
        args.group_col,
        args.top_n,
        group_order=group_order,
        threshold=args.switch_threshold,
    )

    make_plot(
        plot_df,
        taxa,
        args.group_col,
        ordered_groups,
        args.output,
        title=args.title,
        bar_width=args.bar_width,
        legend_cols=1,
        show_sample_labels=args.show_sample_labels,
        within_step=args.within_step,
        group_gap=args.group_gap,
        fig_height=args.fig_height,
        title_fontsize=args.title_fontsize,
        axis_label_fontsize=args.axis_label_fontsize,
        tick_label_fontsize=args.tick_label_fontsize,
        legend_fontsize=args.legend_fontsize,
        group_label_fontsize=args.group_label_fontsize,
        group_label_y=args.group_label_y,
    )

    print(f"Saved plot to: {args.output}")


if __name__ == "__main__":
    main()
