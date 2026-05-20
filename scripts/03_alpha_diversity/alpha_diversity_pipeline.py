#!/usr/bin/env python3
"""
All-in-one alpha diversity pipeline:
  1. Compute alpha diversity metrics from a feature table
  2. Test group significance controlling for covariates (linear model)
  3. Generate violin plots for each metric

Inputs:
  - Feature table TSV (features as rows, samples as columns; BIOM-exported format,
    with a header line starting with '#OTU ID').
  - Metadata TSV (first column = sample_name, must include group + covariate
    columns).

Outputs (written to --output-dir):
  - alpha_diversity_values.tsv         # per-sample metric values
  - alpha_diversity_significance.tsv   # OLS results per metric
  - <metric>_violin.<format>           # one violin plot per metric

Statistical model:
  metric ~ group + covariates, fit via OLS (statsmodels).
  Cohen's d is computed as 2*t / sqrt(df_resid). This is the standard
  t-to-d conversion that assumes approximately balanced group sizes; for
  highly unbalanced groups, use t * sqrt(1/n1 + 1/n2) instead.

Determinism:
  All metric and statistical computations are deterministic. The only RNG-
  dependent step is the stripplot jitter inside the violin plots; pass
  --seed to make jitter reproducible across runs.

Significance brackets:
  Brackets in the violin plots are drawn between the first two groups only
  (two-group comparison). For >2 groups, brackets are omitted.

Requirements:
  pip install pandas numpy matplotlib seaborn scipy statsmodels
"""

import argparse
import os
import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import entropy as shannon_entropy_scipy
import statsmodels.api as sm


# ---------------------------------------------------------------------------
# Alpha diversity metric functions
# ---------------------------------------------------------------------------


def shannon(counts):
    """Shannon entropy (natural log)."""
    p = counts[counts > 0] / counts.sum()
    return shannon_entropy_scipy(p)


def simpson(counts):
    """Simpson's diversity index (1 - dominance)."""
    p = counts[counts > 0] / counts.sum()
    return 1.0 - np.sum(p**2)


def observed_features(counts):
    """Number of observed features (richness)."""
    return np.sum(counts > 0)


METRICS = {
    "shannon_entropy": shannon,
    "simpson": simpson,
    "observed_features": observed_features,
}

DISPLAY_NAMES = {
    "shannon_entropy": "Shannon index",
    "simpson": "Simpson index",
    "observed_features": "Observed features",
}


# ---------------------------------------------------------------------------
# Feature table loading
# ---------------------------------------------------------------------------


def load_feature_table(path):
    """Load a BIOM-exported TSV (features as rows, samples as columns).

    Skips any leading '#' comment lines except the header line that starts
    with '#OTU ID'.
    """
    with open(path) as f:
        skip = 0
        for line in f:
            if line.startswith("#") and not line.startswith("#OTU ID"):
                skip += 1
            else:
                break
    return pd.read_csv(path, sep="\t", index_col=0, skiprows=skip)


# ---------------------------------------------------------------------------
# Significance testing
# ---------------------------------------------------------------------------


def test_significance(diversity_df, metadata, group_col, covariates):
    """
    For each metric, fit: metric ~ group + covariates via OLS.
    Returns a DataFrame with estimate, std_error, t_value, p_value, cohens_d
    for the group term(s).

    Cohen's d uses the balanced-groups approximation 2*t / sqrt(df_resid).
    """
    results = []
    for metric in diversity_df.columns:
        df = (
            diversity_df[[metric]]
            .join(metadata, how="inner")
            .dropna(subset=[metric, group_col] + covariates)
        )

        y = df[metric]
        x_vars = pd.get_dummies(
            df[[group_col] + covariates], drop_first=True, dtype=float
        )
        x_vars = sm.add_constant(x_vars)

        model = sm.OLS(y, x_vars).fit()

        # Find the group coefficient (first dummy after the constant that
        # comes from group_col)
        group_terms = [c for c in model.params.index if c.startswith(group_col + "_")]
        for term in group_terms:
            t_val = model.tvalues[term]
            df_resid = model.df_resid
            results.append(
                {
                    "metric": metric,
                    "term": term,
                    "estimate": model.params[term],
                    "std_error": model.bse[term],
                    "t_value": t_val,
                    "p_value": model.pvalues[term],
                    "cohens_d": 2 * t_val / np.sqrt(df_resid),
                }
            )

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Violin plot
# ---------------------------------------------------------------------------


def pvalue_to_stars(p):
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def plot_violin(
    df,
    metric_col,
    group_col,
    groups,
    group_labels,
    colors,
    ylabel,
    outpath,
    dot_size,
    axis_label_size,
    tick_label_size,
    x_tick_rotation,
    p_value=None,
):
    d = df.copy()
    d[group_col] = d[group_col].astype(str).str.strip()
    d = d[d[group_col].isin(groups)].copy()
    if d.empty:
        return

    d[group_col] = pd.Categorical(d[group_col], categories=groups, ordered=True)
    d = d.sort_values(group_col)

    if group_labels is not None:
        label_map = dict(zip(groups, group_labels))
        d["_display"] = d[group_col].map(label_map)
        x_col, x_order = "_display", group_labels
    else:
        x_col, x_order = group_col, groups

    palette = dict(zip(x_order, colors)) if colors else None

    plt.figure(figsize=(5, 5), dpi=200)
    ax = plt.gca()

    sns.violinplot(
        data=d,
        x=x_col,
        y=metric_col,
        order=x_order,
        palette=palette,
        inner=None,
        cut=0,
        linewidth=1,
        ax=ax,
    )
    sns.stripplot(
        data=d,
        x=x_col,
        y=metric_col,
        order=x_order,
        color="lightgray",
        edgecolor="black",
        linewidth=0.6,
        size=dot_size,
        jitter=0.20,
        alpha=0.95,
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel(ylabel, fontsize=axis_label_size)
    ax.tick_params(axis="both", which="major", labelsize=tick_label_size)
    plt.setp(ax.get_xticklabels(), rotation=x_tick_rotation, ha="right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(-0.7, len(groups) - 0.3)

    # Significance bracket (only drawn for two-group comparisons)
    if p_value is not None and len(groups) == 2:
        y_max = d[metric_col].max()
        y_min = d[metric_col].min()
        y_range = y_max - y_min

        bar_y = y_max + y_range * 0.12  # horizontal bar well above data
        leg_bottom = y_max + y_range * 0.04  # leg starts just above data
        text_y = bar_y + y_range * 0.02
        x0, x1 = 0, 1

        ax.plot(
            [x0, x0, x1, x1],
            [leg_bottom, bar_y, bar_y, leg_bottom],
            color="black",
            linewidth=1.2,
            clip_on=False,
        )
        ax.text(
            (x0 + x1) / 2,
            text_y,
            pvalue_to_stars(p_value),
            ha="center",
            va="bottom",
            fontsize=axis_label_size * 0.85,
        )

        # Expand y-axis to fit bracket + label
        ax.set_ylim(bottom=y_min - y_range * 0.05, top=text_y + y_range * 0.12)

    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compute alpha diversity, test significance, and generate violin plots."
        )
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Exported feature table TSV (features as rows, samples as columns).",
    )
    parser.add_argument("--metadata", required=True, help="Metadata TSV file.")
    parser.add_argument(
        "--output-dir", required=True, help="Output directory for all results."
    )
    parser.add_argument(
        "--group-col",
        default="opioid_status",
        help="Metadata column defining groups (default: opioid_status).",
    )
    parser.add_argument(
        "--groups",
        default="Remission,OUD",
        help=(
            "Comma-delimited ordered group values (default: Remission,OUD). "
            "Must match the values in the --group-col exactly (case-sensitive)."
        ),
    )
    parser.add_argument(
        "--group-labels",
        default=None,
        help="Comma-delimited display labels for groups (default: same as --groups).",
    )
    parser.add_argument(
        "--covariates",
        default="age,sex",
        help="Comma-delimited covariates for the linear model (default: age,sex).",
    )
    parser.add_argument(
        "--metrics",
        default=",".join(METRICS.keys()),
        help="Comma-delimited metrics to compute (default: all).",
    )
    parser.add_argument(
        "--colors",
        default="#FFC20A,#0C7BDC",
        help="Comma-delimited colors for violin plots (default: #FFC20A,#0C7BDC).",
    )
    parser.add_argument(
        "--plot-format", default="pdf", help="Plot file format (default: pdf)."
    )
    parser.add_argument("--dot-size", type=float, default=8.0)
    parser.add_argument("--axis-label-size", type=float, default=20)
    parser.add_argument("--tick-label-size", type=float, default=18)
    parser.add_argument("--x-tick-rotation", type=float, default=45)
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for stripplot jitter (default: None; jitter varies per run).",
    )
    args = parser.parse_args()

    # Seed RNGs for reproducible jitter
    if args.seed is not None:
        random.seed(args.seed)
        np.random.seed(args.seed)

    os.makedirs(args.output_dir, exist_ok=True)

    groups = [g.strip() for g in args.groups.split(",")]
    group_labels = (
        [g.strip() for g in args.group_labels.split(",")] if args.group_labels else None
    )
    covariates = [c.strip() for c in args.covariates.split(",")]
    colors = [c.strip() for c in args.colors.split(",")]
    metrics_to_run = [m.strip() for m in args.metrics.split(",")]

    # Load data
    table = load_feature_table(args.input)
    metadata = pd.read_csv(args.metadata, sep="\t")
    if "sample_name" not in metadata.columns:
        metadata = metadata.rename(columns={metadata.columns[0]: "sample_name"})
    metadata = metadata.set_index("sample_name")

    common = list(set(table.columns) & set(metadata.index))
    if not common:
        raise ValueError(
            "No overlapping sample names between feature table and metadata."
        )
    table = table[common]

    # 1. Compute alpha diversity
    print("Computing alpha diversity metrics...")
    diversity = {}
    for name in metrics_to_run:
        if name not in METRICS:
            print(f"  WARNING: unknown metric '{name}', skipping.")
            continue
        func = METRICS[name]
        diversity[name] = {
            sample: func(table[sample].values.astype(float)) for sample in common
        }
        print(f"  {name}: done")

    diversity_df = pd.DataFrame(diversity)
    diversity_df.index.name = "sample_name"

    div_path = os.path.join(args.output_dir, "alpha_diversity_values.tsv")
    diversity_df.to_csv(div_path, sep="\t")
    print(f"Diversity values saved to: {div_path}")

    # 2. Test significance
    print("\nTesting group significance...")
    results = test_significance(diversity_df, metadata, args.group_col, covariates)
    sig_path = os.path.join(args.output_dir, "alpha_diversity_significance.tsv")
    results.to_csv(sig_path, sep="\t", index=False)
    print(results.to_string(index=False))
    print(f"\nSignificance results saved to: {sig_path}")

    # 3. Generate violin plots
    # Build a lookup of p-values per metric
    pval_lookup = {}
    for _, row in results.iterrows():
        pval_lookup[row["metric"]] = row["p_value"]

    print("\nGenerating violin plots...")
    plot_df = diversity_df.join(metadata, how="inner")

    for name in diversity.keys():
        outpath = os.path.join(args.output_dir, f"{name}_violin.{args.plot_format}")
        ylabel = DISPLAY_NAMES.get(name, name)
        plot_violin(
            df=plot_df,
            metric_col=name,
            group_col=args.group_col,
            groups=groups,
            group_labels=group_labels,
            colors=colors,
            ylabel=ylabel,
            outpath=outpath,
            dot_size=args.dot_size,
            axis_label_size=args.axis_label_size,
            tick_label_size=args.tick_label_size,
            x_tick_rotation=args.x_tick_rotation,
            p_value=pval_lookup.get(name),
        )
        print(f"  {outpath}")

    print("\nDone.")


if __name__ == "__main__":
    main()
