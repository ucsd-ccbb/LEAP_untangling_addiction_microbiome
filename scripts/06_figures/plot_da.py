#!/usr/bin/env python3

import argparse
import json
import textwrap
import zipfile
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def read_slice_csv(path):
    df = pd.read_csv(path)
    first = df.columns[0]
    df = df.rename(columns={first: "feature"})
    return df


def read_qza_jsonl(qza_path, slice_name):
    """Read a JSONL slice (e.g. 'lfc', 'q') from an ANCOM-BC2 .qza archive.

    Returns a DataFrame with 'feature' as the first column and one column
    per model term, matching the shape that read_slice_csv() produces.
    """
    target = f"data/{slice_name}.jsonl"
    with zipfile.ZipFile(qza_path, "r") as zf:
        matches = [n for n in zf.namelist() if n.endswith(target)]
        if not matches:
            available = sorted(
                {
                    n.split("data/")[-1].replace(".jsonl", "")
                    for n in zf.namelist()
                    if "/data/" in n and n.endswith(".jsonl")
                }
            )
            raise FileNotFoundError(
                f"Slice '{slice_name}' not found in {qza_path}.\n"
                f"Available slices: {', '.join(available)}"
            )
        with zf.open(matches[0]) as f:
            lines = f.read().decode("utf-8").splitlines()

    # Line 0 is the header/schema; lines 1+ are data rows
    header = json.loads(lines[0])
    id_field = header["fields"][0]["name"]  # typically 'taxon'

    rows = [json.loads(line) for line in lines[1:] if line.strip()]
    df = pd.DataFrame(rows)
    df = df.rename(columns={id_field: "feature"})
    return df


def list_qza_columns(qza_path, slice_name="lfc"):
    """Return the list of value column names from a .qza slice."""
    target = f"data/{slice_name}.jsonl"
    with zipfile.ZipFile(qza_path, "r") as zf:
        matches = [n for n in zf.namelist() if n.endswith(target)]
        if not matches:
            return []
        with zf.open(matches[0]) as f:
            header = json.loads(f.readline().decode("utf-8"))
    # Skip the first field (taxon/feature ID)
    return [field["name"] for field in header["fields"][1:]]


def choose_value_column(df, requested=None):
    cols = [c for c in df.columns if c != "feature"]
    if requested is not None:
        if requested not in df.columns:
            raise ValueError(
                f"Requested column '{requested}' not found.\n"
                f"Available columns: {', '.join(cols)}"
            )
        return requested
    if len(cols) == 1:
        return cols[0]
    raise ValueError(
        "Multiple value columns found. Please specify one with --column.\n"
        f"Available columns: {', '.join(cols)}"
    )


def read_mapping_tsv(path, key_col=None, value_col=None):
    df = pd.read_csv(path, sep="\t", dtype=str)
    if key_col is None:
        key_col = df.columns[0]
    if value_col is None:
        if len(df.columns) < 2:
            raise ValueError("Mapping TSV must contain at least two columns.")
        value_col = df.columns[1]
    return dict(zip(df[key_col], df[value_col]))


def shorten_taxonomy(label):
    if isinstance(label, str) and ";" in label:
        parts = [x.strip() for x in label.split(";") if x.strip()]
        if parts:
            return parts[-1]
    return label


def wrap_label(label, width):
    if not isinstance(label, str):
        label = str(label)
    if width is None or width <= 0:
        return label
    return "\n".join(textwrap.wrap(label, width=width, break_long_words=False))


def significance_stars(qval, t1=0.05, t2=0.01, t3=0.001):
    if pd.isna(qval):
        return ""
    try:
        q = float(qval)
    except Exception:
        return ""
    if q < t3:
        return "***"
    if q < t2:
        return "**"
    if q < t1:
        return "*"
    return ""


def build_plot_df(
    lfc_path=None,
    qval_path=None,
    qza_path=None,
    column=None,
    q_threshold=None,
    require_ss=True,
    mapping_tsv=None,
    mapping_key_col=None,
    mapping_value_col=None,
    taxonomy=False,
    wrap_width=None,
    star_threshold_1=0.05,
    star_threshold_2=0.01,
    star_threshold_3=0.001,
):
    # --- Read LFC data ---
    if qza_path is not None:
        lfc_df = read_qza_jsonl(qza_path, "lfc")
    else:
        lfc_df = read_slice_csv(lfc_path)

    value_col = choose_value_column(lfc_df, column)

    plot_df = lfc_df[["feature", value_col]].copy()
    plot_df = plot_df.rename(columns={value_col: "effect"})
    plot_df["effect"] = pd.to_numeric(plot_df["effect"], errors="coerce")
    plot_df = plot_df.dropna(subset=["effect"])

    # --- Read q-value data ---
    if qza_path is not None:
        q_df = read_qza_jsonl(qza_path, "q")
        q_col = choose_value_column(q_df, column)
        q_df = q_df[["feature", q_col]].copy()
        q_df = q_df.rename(columns={q_col: "qval"})
        q_df["qval"] = pd.to_numeric(q_df["qval"], errors="coerce")
        plot_df = plot_df.merge(q_df, on="feature", how="left")
    elif qval_path is not None:
        q_df = read_slice_csv(qval_path)
        q_col = choose_value_column(q_df, column)
        q_df = q_df[["feature", q_col]].copy()
        q_df = q_df.rename(columns={q_col: "qval"})
        q_df["qval"] = pd.to_numeric(q_df["qval"], errors="coerce")
        plot_df = plot_df.merge(q_df, on="feature", how="left")
    else:
        plot_df["qval"] = pd.NA

    if q_threshold is not None:
        plot_df = plot_df[plot_df["qval"] <= q_threshold]

    # --- Filter by sensitivity analysis (passed_ss) ---
    if require_ss and qza_path is not None:
        try:
            ss_df = read_qza_jsonl(qza_path, "passed_ss")
            ss_col = choose_value_column(ss_df, column)
            ss_df = ss_df[["feature", ss_col]].copy()
            ss_df = ss_df.rename(columns={ss_col: "passed_ss"})
            plot_df = plot_df.merge(ss_df, on="feature", how="left")
            n_before = len(plot_df)
            plot_df = plot_df[plot_df["passed_ss"] == True]
            n_after = len(plot_df)
            if n_before != n_after:
                print(
                    f"Sensitivity filter: {n_before} -> {n_after} taxa "
                    f"({n_before - n_after} removed)"
                )
            plot_df = plot_df.drop(columns=["passed_ss"])
        except FileNotFoundError:
            print(
                "Warning: passed_ss slice not found in QZA; skipping sensitivity filter."
            )

    plot_df["label"] = plot_df["feature"]

    if mapping_tsv:
        mapping = read_mapping_tsv(mapping_tsv, mapping_key_col, mapping_value_col)
        plot_df["label"] = plot_df["label"].map(lambda x: mapping.get(x, x))

    if taxonomy:
        plot_df["label"] = plot_df["label"].apply(shorten_taxonomy)

    if wrap_width:
        plot_df["label"] = plot_df["label"].apply(lambda x: wrap_label(x, wrap_width))

    plot_df["stars"] = plot_df["qval"].apply(
        lambda q: significance_stars(
            q,
            t1=star_threshold_1,
            t2=star_threshold_2,
            t3=star_threshold_3,
        )
    )

    return plot_df


def select_features(df, top_n=None, plot_all=False, sort_mode="signed"):
    if df.empty:
        return df.copy()

    if plot_all:
        selected = df.copy()
    else:
        if top_n is None:
            top_n = 15
        pos = df[df["effect"] > 0].sort_values("effect", ascending=False).head(top_n)
        neg = df[df["effect"] < 0].sort_values("effect", ascending=True).head(top_n)
        selected = pd.concat([pos, neg], axis=0)

    if selected.empty:
        return selected

    if sort_mode == "signed":
        selected = selected.sort_values("effect", ascending=True)
    elif sort_mode == "absolute":
        neg = selected[selected["effect"] < 0].copy()
        pos = selected[selected["effect"] > 0].copy()
        neg = neg.reindex(neg["effect"].abs().sort_values(ascending=False).index)
        pos = pos.reindex(pos["effect"].abs().sort_values(ascending=False).index)
        selected = pd.concat([neg, pos], axis=0)
    else:
        raise ValueError("sort_mode must be 'signed' or 'absolute'")

    return selected


def get_legend_anchor(legend_loc, x_shift=0.0, y_shift=0.0):
    base_anchor = {
        "lower right": (0.98, 0.03),
        "lower left": (0.02, 0.03),
        "upper right": (0.98, 0.97),
        "upper left": (0.02, 0.97),
        "center right": (0.98, 0.50),
        "center left": (0.02, 0.50),
        "lower center": (0.50, 0.03),
        "upper center": (0.50, 0.97),
        "center": (0.50, 0.50),
    }.get(legend_loc, (0.98, 0.03))

    return (base_anchor[0] + x_shift, base_anchor[1] + y_shift)


def make_plot(
    df,
    output,
    title,
    xlabel,
    enriched_color,
    depleted_color,
    fig_width=9,
    row_height=0.45,
    axes_padding=2.0,
    font_size=10,
    title_size=None,
    label_size=None,
    xlabel_size=None,
    tick_size=None,
    center_label_offset_frac=0.02,
    edgecolor="black",
    linewidth=1.0,
    grid=True,
    grid_alpha=0.5,
    center_linewidth=1.2,
    x_padding_frac=0.10,
    show_values=False,
    value_size=None,
    show_legend=True,
    legend_size=None,
    legend_loc="auto",
    legend_x_shift=0.0,
    legend_y_shift=0.0,
    show_stars=False,
    star_size=None,
    star_offset_frac=0.025,
):
    if df.empty:
        raise ValueError("No rows left to plot.")

    # Resolve all sizes relative to font_size when not explicitly set
    if title_size is None:
        title_size = font_size + 4
    if label_size is None:
        label_size = font_size
    if xlabel_size is None:
        xlabel_size = font_size + 1
    if tick_size is None:
        tick_size = font_size
    if value_size is None:
        value_size = font_size - 1
    if legend_size is None:
        legend_size = font_size
    if star_size is None:
        star_size = font_size + 1

    plt.rcParams["font.size"] = font_size

    n = len(df)
    fig_height = row_height * n + axes_padding

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    y = list(range(n))
    colors = [depleted_color if x < 0 else enriched_color for x in df["effect"]]

    max_abs = float(df["effect"].abs().max())
    xlim = max_abs * (1 + x_padding_frac)
    if xlim == 0:
        xlim = 1.0

    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(-0.6, n - 0.4)

    # Measure actual text height in data coordinates so bar thickness
    # and star positioning scale with the font, not with the figure.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    tmp = ax.text(0, 0, "Xg", fontsize=label_size)
    tmp_bb = tmp.get_window_extent(renderer)
    tmp.remove()
    inv = ax.transData.inverted()
    text_height_data = inv.transform((0, tmp_bb.height))[1] - inv.transform((0, 0))[1]

    bar_height = text_height_data * 1.3

    ax.barh(
        y,
        df["effect"],
        height=bar_height,
        color=colors,
        edgecolor=edgecolor,
        linewidth=linewidth,
    )

    ax.axvline(0, color="black", linewidth=center_linewidth)

    if grid:
        ax.xaxis.grid(True, linestyle=":", alpha=grid_alpha)
        ax.set_axisbelow(True)

    ax.set_yticks(y)
    ax.set_yticklabels([""] * n)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=tick_size)

    label_offset = xlim * center_label_offset_frac
    star_offset = xlim * star_offset_frac

    for yi, (_, row) in zip(y, df.iterrows()):
        effect = row["effect"]
        label = row["label"]

        if effect < 0:
            x = label_offset
            ha = "left"
        else:
            x = -label_offset
            ha = "right"

        ax.text(
            x,
            yi,
            str(label),
            ha=ha,
            va="center",
            fontsize=label_size,
            bbox=dict(facecolor="white", edgecolor="none", pad=0.5),
            zorder=3,
        )

        if show_values:
            val = f"{effect:.2f}"
            if effect < 0:
                vx = effect - (xlim * 0.02)
                vha = "right"
            else:
                vx = effect + (xlim * 0.02)
                vha = "left"

            ax.text(
                vx,
                yi,
                val,
                ha=vha,
                va="center",
                fontsize=value_size,
            )

        if show_stars:
            stars = row.get("stars", "")
            if stars:
                if effect < 0:
                    sx = effect - star_offset
                    sha = "right"
                else:
                    sx = effect + star_offset
                    sha = "left"

                ax.text(
                    sx,
                    yi - text_height_data * 0.15,
                    stars,
                    ha=sha,
                    va="center",
                    fontsize=star_size,
                    fontweight="bold",
                )

    if show_legend:
        handles = [
            Patch(
                facecolor=enriched_color, edgecolor=edgecolor, label="Enriched in OUD"
            ),
            Patch(
                facecolor=depleted_color, edgecolor=edgecolor, label="Depleted in OUD"
            ),
        ]

        legend_kwargs = dict(
            borderaxespad=0.4,
            frameon=True,
            fancybox=True,
            framealpha=1.0,
            fontsize=legend_size,
            title="Relative to reference",
            title_fontsize=legend_size + 1,
            handlelength=1.0,
            handleheight=1.0,
            labelspacing=0.4,
            borderpad=0.6,
            handletextpad=0.5,
        )

        if legend_loc == "auto":
            # Place legend outside the axes, below the xlabel, so it
            # can never overlap bar labels or stars.  bbox_inches='tight'
            # in savefig will expand the saved figure to include it.
            ax.set_xlabel(xlabel, fontsize=xlabel_size)
            ax.set_title(title, fontsize=title_size, pad=10)
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()

            # Find the bottom edge of the xlabel in axes-fraction coords
            xlabel_obj = ax.xaxis.get_label()
            xlabel_bb = xlabel_obj.get_window_extent(renderer)
            inv = ax.transAxes.inverted()
            xlabel_bottom_frac = inv.transform((0, xlabel_bb.y0))[1]

            # Small gap below the xlabel
            auto_y = xlabel_bottom_frac - 0.02
            anchor = (0.5 + legend_x_shift, auto_y + legend_y_shift)

            legend = ax.legend(
                handles=handles,
                loc="upper center",
                bbox_to_anchor=anchor,
                **legend_kwargs,
            )
        else:
            anchor_x, anchor_y = get_legend_anchor(
                legend_loc,
                x_shift=legend_x_shift,
                y_shift=legend_y_shift,
            )
            legend = ax.legend(
                handles=handles,
                loc=legend_loc,
                bbox_to_anchor=(anchor_x, anchor_y),
                **legend_kwargs,
            )

        frame = legend.get_frame()
        frame.set_facecolor("#f2f2f2")
        frame.set_edgecolor("#7a7a7a")
        frame.set_linewidth(1.0)

    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)

    # Set xlabel/title here only if auto-legend didn't already set them
    if not (show_legend and legend_loc == "auto"):
        ax.set_xlabel(xlabel, fontsize=xlabel_size)
        ax.set_title(title, fontsize=title_size, pad=10)

    # bbox_inches='tight' handles layout; tight_layout() can fail when
    # the legend sits outside the axes.
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Create a polished LEfSe-style mirrored bar plot from QIIME 2 ANCOM-BC / ANCOM-BC2 output."
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--qza",
        default=None,
        help="Path to ANCOM-BC2 .qza artifact (contains lfc/q slices as JSONL)",
    )
    input_group.add_argument(
        "--lfc", default=None, help="Path to lfc_slice.csv (ANCOM-BC CSV export)"
    )
    parser.add_argument(
        "--qval",
        default=None,
        help="Optional path to q_val_slice.csv (ignored when using --qza)",
    )
    parser.add_argument(
        "--column",
        default=None,
        help="Column name to plot (e.g. 'opioid_status::Remission')",
    )
    parser.add_argument(
        "--list-columns",
        action="store_true",
        help="Print available column names from input and exit",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="da_barplot.pdf",
        help="Output file (.png, .pdf, etc.)",
    )
    parser.add_argument("--title", default="Differential abundance", help="Plot title")
    parser.add_argument("--xlabel", default="Log fold change", help="X-axis label")

    parser.add_argument(
        "--q-threshold",
        type=float,
        default=None,
        help="Optional q-value cutoff (e.g. 0.05)",
    )
    parser.add_argument(
        "--require-ss",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Require taxa to pass ANCOM-BC2 sensitivity analysis (default: --require-ss). "
        "Use --no-require-ss to disable. Only applies to --qza input.",
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--plot-all",
        action="store_true",
        help="Plot all remaining features after filtering",
    )
    group.add_argument(
        "--top-n",
        type=int,
        default=15,
        help="Top N positive and top N negative features",
    )

    parser.add_argument(
        "--sort-mode",
        choices=["signed", "absolute"],
        default="signed",
        help="signed = order by effect from negative to positive; absolute = larger magnitudes first within each side",
    )

    parser.add_argument(
        "--enriched-color", default="red", help="Color for positive values"
    )
    parser.add_argument(
        "--depleted-color", default="blue", help="Color for negative values"
    )

    parser.add_argument(
        "--mapping-tsv",
        default=None,
        help="Optional TSV mapping file for renaming features",
    )
    parser.add_argument(
        "--mapping-key-col", default=None, help="Mapping TSV key column"
    )
    parser.add_argument(
        "--mapping-value-col", default=None, help="Mapping TSV display-name column"
    )

    parser.add_argument(
        "--taxonomy", action="store_true", help="Use last semicolon-delimited rank only"
    )
    parser.add_argument(
        "--wrap-width",
        type=int,
        default=None,
        help="Wrap label text after this many characters",
    )

    parser.add_argument("--fig-width", type=float, default=9, help="Figure width")
    parser.add_argument(
        "--row-height", type=float, default=0.45, help="Height per row in inches"
    )
    parser.add_argument(
        "--axes-padding",
        type=float,
        default=2.0,
        help="Extra figure height in inches for x-axis, labels, and legend",
    )

    parser.add_argument(
        "--font-size",
        type=float,
        default=10,
        help="Base font size (all other sizes scale from this unless overridden)",
    )
    parser.add_argument(
        "--title-size",
        type=float,
        default=None,
        help="Title font size (default: font-size + 4)",
    )
    parser.add_argument(
        "--label-size",
        type=float,
        default=None,
        help="Center label font size (default: font-size)",
    )
    parser.add_argument(
        "--xlabel-size",
        type=float,
        default=None,
        help="X-axis label font size (default: font-size + 1)",
    )
    parser.add_argument(
        "--tick-size",
        type=float,
        default=None,
        help="X-axis tick label font size (default: font-size)",
    )

    parser.add_argument(
        "--center-label-offset-frac",
        type=float,
        default=0.02,
        help="Distance of labels from zero line as fraction of x-axis half-width",
    )
    parser.add_argument(
        "--x-padding-frac",
        type=float,
        default=0.10,
        help="Extra x-axis padding as fraction of max absolute value",
    )

    parser.add_argument(
        "--show-values",
        action="store_true",
        help="Show numeric effect values at bar ends",
    )
    parser.add_argument(
        "--value-size",
        type=float,
        default=None,
        help="Font size for numeric values (default: font-size - 1)",
    )

    parser.add_argument(
        "--show-legend", action="store_true", help="Show enriched/depleted legend"
    )
    parser.add_argument(
        "--legend-size",
        type=float,
        default=None,
        help="Legend font size (default: font-size)",
    )
    parser.add_argument(
        "--legend-loc",
        default="auto",
        help="Legend anchor location (auto, lower right, upper right, lower left, upper left, etc.)",
    )
    parser.add_argument(
        "--legend-x-shift",
        type=float,
        default=0.0,
        help="Horizontal legend shift relative to anchored position (negative=left, positive=right)",
    )
    parser.add_argument(
        "--legend-y-shift",
        type=float,
        default=0.0,
        help="Vertical legend shift relative to anchored position (negative=down, positive=up)",
    )

    parser.add_argument(
        "--show-stars", action="store_true", help="Show significance stars at bar ends"
    )
    parser.add_argument(
        "--star-size",
        type=float,
        default=None,
        help="Font size for significance stars (default: font-size + 1)",
    )
    parser.add_argument(
        "--star-threshold-1",
        type=float,
        default=0.1,
        help="Threshold for one star (*)",
    )
    parser.add_argument(
        "--star-threshold-2",
        type=float,
        default=0.01,
        help="Threshold for two stars (**)",
    )
    parser.add_argument(
        "--star-threshold-3",
        type=float,
        default=0.001,
        help="Threshold for three stars (***)",
    )

    args = parser.parse_args()

    # --list-columns: print available columns and exit
    if args.list_columns:
        if args.qza:
            cols = list_qza_columns(args.qza, "lfc")
        else:
            df = read_slice_csv(args.lfc)
            cols = [c for c in df.columns if c != "feature"]
        print("Available columns:")
        for c in cols:
            print(f"  {c}")
        return

    plot_df = build_plot_df(
        lfc_path=args.lfc,
        qval_path=args.qval,
        qza_path=args.qza,
        column=args.column,
        q_threshold=args.q_threshold,
        require_ss=args.require_ss,
        mapping_tsv=args.mapping_tsv,
        mapping_key_col=args.mapping_key_col,
        mapping_value_col=args.mapping_value_col,
        taxonomy=args.taxonomy,
        wrap_width=args.wrap_width,
        star_threshold_1=args.star_threshold_1,
        star_threshold_2=args.star_threshold_2,
        star_threshold_3=args.star_threshold_3,
    )

    selected = select_features(
        plot_df,
        top_n=args.top_n,
        plot_all=args.plot_all,
        sort_mode=args.sort_mode,
    )

    make_plot(
        df=selected,
        output=args.output,
        title=args.title,
        xlabel=args.xlabel,
        enriched_color=args.enriched_color,
        depleted_color=args.depleted_color,
        fig_width=args.fig_width,
        row_height=args.row_height,
        axes_padding=args.axes_padding,
        font_size=args.font_size,
        title_size=args.title_size,
        label_size=args.label_size,
        xlabel_size=args.xlabel_size,
        tick_size=args.tick_size,
        center_label_offset_frac=args.center_label_offset_frac,
        x_padding_frac=args.x_padding_frac,
        show_values=args.show_values,
        value_size=args.value_size,
        show_legend=args.show_legend,
        legend_size=args.legend_size,
        legend_loc=args.legend_loc,
        legend_x_shift=args.legend_x_shift,
        legend_y_shift=args.legend_y_shift,
        show_stars=args.show_stars,
        star_size=args.star_size,
    )


if __name__ == "__main__":
    main()
