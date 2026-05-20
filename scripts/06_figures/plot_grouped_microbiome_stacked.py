#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


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


def prepare_grouped_data(rel_df, meta_df, group_col, top_n, group_order=None):
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
    grouped = plot_df.groupby(group_col).mean()

    if group_order:
        missing = [g for g in group_order if g not in grouped.index]
        if missing:
            raise ValueError(f"These groups were not found: {missing}")
        grouped = grouped.loc[group_order]

    ordered_taxa = top_taxa.copy()
    if "Other" in grouped.columns:
        ordered_taxa.append("Other")

    return grouped[ordered_taxa]


def compute_stack_bounds(grouped_df):
    taxa = list(grouped_df.columns)
    n_groups = grouped_df.shape[0]

    lower = {taxon: np.zeros(n_groups) for taxon in taxa}
    upper = {taxon: np.zeros(n_groups) for taxon in taxa}

    cumulative = np.zeros(n_groups)
    for taxon in taxa:
        vals = grouped_df[taxon].values
        lower[taxon] = cumulative.copy()
        upper[taxon] = cumulative + vals
        cumulative = upper[taxon].copy()

    return lower, upper


def desaturate_color(color, amount=0.45, lighten=0.10):
    rgb = np.array(mcolors.to_rgb(color))
    gray = np.array([rgb.mean(), rgb.mean(), rgb.mean()])
    muted = gray * (1 - amount) + rgb * amount
    muted = muted * (1 - lighten) + np.ones(3) * lighten
    return tuple(np.clip(muted, 0, 1))


def build_muted_palette(n):
    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    colors = []
    for i in range(n):
        colors.append(desaturate_color(base[i % len(base)], amount=0.42, lighten=0.08))
    return colors


def make_plot(
    grouped_df,
    output_file,
    title=None,
    fig_width=10,
    fig_height=6,
    bar_width=0.78,
    group_spacing=1.0,
    connector_alpha=0.38,
    legend_cols=2,
    title_fontsize=14,
    axis_label_fontsize=12,
    tick_label_fontsize=10,
    legend_fontsize=9,
):
    taxa = list(grouped_df.columns)
    groups = list(grouped_df.index)

    if group_spacing <= 0:
        raise ValueError("--group_spacing must be > 0")

    x = np.arange(len(groups), dtype=float) * group_spacing

    palette = build_muted_palette(len(taxa))
    colors = {taxon: palette[i] for i, taxon in enumerate(taxa)}
    if "Other" in colors:
        colors["Other"] = (0.82, 0.82, 0.82)

    display_labels = {
        taxon: ("Other" if taxon == "Other" else short_taxon_name(taxon))
        for taxon in taxa
    }

    lower, upper = compute_stack_bounds(grouped_df)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Connector ribbons: always from bar edge to bar edge
    for i in range(len(groups) - 1):
        x_left = x[i] + bar_width / 2.0
        x_right = x[i + 1] - bar_width / 2.0

        if x_right <= x_left:
            raise ValueError(
                f"Bars overlap or touch too strongly for the chosen settings: "
                f"bar_width={bar_width}, group_spacing={group_spacing}. "
                f"Use a smaller --bar_width or larger --group_spacing."
            )

        xs = np.linspace(x_left, x_right, 160)
        t = np.linspace(0, 1, 160)
        s = 3 * t**2 - 2 * t**3

        for taxon in taxa:
            y1_low = lower[taxon][i]
            y2_low = lower[taxon][i + 1]
            y1_up = upper[taxon][i]
            y2_up = upper[taxon][i + 1]

            ys_low = y1_low + (y2_low - y1_low) * s
            ys_up = y1_up + (y2_up - y1_up) * s

            ax.fill_between(
                xs,
                ys_low * 100,
                ys_up * 100,
                color=colors[taxon],
                alpha=connector_alpha,
                linewidth=0,
                zorder=1,
            )

    bottoms = np.zeros(len(groups))
    bar_handles = []
    for taxon in taxa:
        vals = grouped_df[taxon].values * 100
        bars = ax.bar(
            x,
            vals,
            width=bar_width,
            bottom=bottoms,
            color=colors[taxon],
            edgecolor="#444444",
            linewidth=0.45,
            zorder=2,
        )
        bar_handles.append(bars[0])
        bottoms += vals

    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=45, ha="right", fontsize=tick_label_fontsize)
    ax.set_ylabel("Relative abundance (%)", fontsize=axis_label_fontsize)
    ax.tick_params(axis="y", labelsize=tick_label_fontsize)
    ax.set_ylim(0, 100)

    if len(x) > 0:
        left_margin = bar_width * 0.7
        right_margin = bar_width * 0.7
        ax.set_xlim(x[0] - left_margin, x[-1] + right_margin)

    if title is not None:
        ax.set_title(title, fontsize=title_fontsize)

    ax.legend(
        handles=bar_handles,
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

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Stacked relative abundance bar plot with connector ribbons."
    )
    parser.add_argument("--abundance_tsv", required=True)
    parser.add_argument("--metadata_tsv", required=True)
    parser.add_argument("--group_col", required=True)
    parser.add_argument("--sample_col", default=None)
    parser.add_argument("--top_n", type=int, default=10)
    parser.add_argument("--output", required=True)
    parser.add_argument("--title", default=None)
    parser.add_argument(
        "--group_order", default=None, help="Comma-separated order of groups"
    )

    parser.add_argument("--fig_width", type=float, default=10)
    parser.add_argument("--fig_height", type=float, default=6)

    parser.add_argument(
        "--bar_width",
        type=float,
        default=0.78,
        help="Width of the stacked bars in x-axis units (default: 0.78)",
    )
    parser.add_argument(
        "--group_spacing",
        type=float,
        default=1.0,
        help="Distance between centers of adjacent group columns (default: 1.0)",
    )
    parser.add_argument(
        "--connector_alpha",
        type=float,
        default=0.38,
        help="Transparency of ribbons (default: 0.38)",
    )
    parser.add_argument(
        "--legend_cols",
        type=int,
        default=2,
        help="Number of legend columns (default: 2)",
    )

    parser.add_argument(
        "--title_fontsize",
        type=float,
        default=14,
        help="Plot title font size (default: 14)",
    )
    parser.add_argument(
        "--axis_label_fontsize",
        type=float,
        default=12,
        help="Axis label font size (default: 12)",
    )
    parser.add_argument(
        "--tick_label_fontsize",
        type=float,
        default=10,
        help="Axis tick label font size (default: 10)",
    )
    parser.add_argument(
        "--legend_fontsize",
        type=float,
        default=9,
        help="Legend text font size (default: 9)",
    )

    args = parser.parse_args()

    abundance_df = load_abundance_table(args.abundance_tsv)
    metadata_df = load_metadata(args.metadata_tsv, sample_col=args.sample_col)
    rel_df = relative_frequencies(abundance_df)

    group_order = None
    if args.group_order:
        group_order = [x.strip() for x in args.group_order.split(",")]

    grouped_df = prepare_grouped_data(
        rel_df, metadata_df, args.group_col, args.top_n, group_order=group_order
    )

    make_plot(
        grouped_df=grouped_df,
        output_file=args.output,
        title=args.title,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
        bar_width=args.bar_width,
        group_spacing=args.group_spacing,
        connector_alpha=args.connector_alpha,
        legend_cols=args.legend_cols,
        title_fontsize=args.title_fontsize,
        axis_label_fontsize=args.axis_label_fontsize,
        tick_label_fontsize=args.tick_label_fontsize,
        legend_fontsize=args.legend_fontsize,
    )

    print("Genera order from bottom to top:")
    for taxon in grouped_df.columns:
        print("Other" if taxon == "Other" else short_taxon_name(taxon))

    print(f"\nSaved plot to: {args.output}")


if __name__ == "__main__":
    main()
