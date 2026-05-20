#!/usr/bin/env python3
"""Plot Hill diversity profiles by group from a QIIME 2 feature table."""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def hill_number(counts, q):
    """Compute the Hill number of order q for an abundance vector."""
    p = counts[counts > 0] / counts.sum()
    if q == 1:
        return np.exp(-np.sum(p * np.log(p)))
    else:
        return np.sum(p**q) ** (1.0 / (1.0 - q))


def main():
    parser = argparse.ArgumentParser(
        description="Plot Hill diversity profiles by group."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to exported feature table TSV (samples as columns, features as rows).",
    )
    parser.add_argument("--metadata", required=True, help="Path to metadata TSV file.")
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output plot (e.g., hill_diversity.png).",
    )
    parser.add_argument(
        "--group-col",
        default="opioid_status",
        help="Metadata column defining groups (default: opioid_status).",
    )
    parser.add_argument(
        "--q-min", type=float, default=0.0, help="Minimum q value (default: 0)."
    )
    parser.add_argument(
        "--q-max", type=float, default=3.0, help="Maximum q value (default: 3)."
    )
    parser.add_argument(
        "--q-steps",
        type=int,
        default=50,
        help="Number of q values to evaluate (default: 50).",
    )
    parser.add_argument(
        "--tsv-output",
        default=None,
        help="Optional path to save per-sample Hill numbers as TSV.",
    )
    args = parser.parse_args()

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Load feature table (features as rows, samples as columns)
    # BIOM exports have a "# Constructed from biom file" comment line
    # followed by a header starting with "#OTU ID" — skip the comment,
    # keep the header.
    with open(args.input) as f:
        skip = 0
        for line in f:
            if line.startswith("#") and not line.startswith("#OTU ID"):
                skip += 1
            else:
                break
    table = pd.read_csv(args.input, sep="\t", index_col=0, skiprows=skip)
    metadata = pd.read_csv(args.metadata, sep="\t")

    # Ensure sample_name column
    if "sample_name" not in metadata.columns:
        metadata = metadata.rename(columns={metadata.columns[0]: "sample_name"})

    common_samples = list(set(table.columns) & set(metadata["sample_name"]))
    if not common_samples:
        raise ValueError(
            "No overlapping sample names between feature table and metadata."
        )

    table = table[common_samples]
    metadata = metadata[metadata["sample_name"].isin(common_samples)].set_index(
        "sample_name"
    )

    q_values = np.linspace(args.q_min, args.q_max, args.q_steps)

    # Compute Hill numbers for each sample at each q
    records = []
    for sample in common_samples:
        counts = table[sample].values.astype(float)
        group = metadata.loc[sample, args.group_col]
        for q in q_values:
            d = hill_number(counts, q)
            records.append({"sample": sample, "group": group, "q": q, "D": d})

    df = pd.DataFrame(records)

    if args.tsv_output:
        tsv_dir = os.path.dirname(args.tsv_output)
        if tsv_dir:
            os.makedirs(tsv_dir, exist_ok=True)
        df.to_csv(args.tsv_output, sep="\t", index=False)

    # Summarize: mean and SEM per group per q
    summary = df.groupby(["group", "q"])["D"].agg(["mean", "sem"]).reset_index()

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    groups = sorted(summary["group"].unique())
    colors = plt.cm.Set1.colors

    for i, group in enumerate(groups):
        g = summary[summary["group"] == group]
        color = colors[i % len(colors)]
        ax.plot(g["q"], g["mean"], label=group, color=color, linewidth=2)
        ax.fill_between(
            g["q"],
            g["mean"] - g["sem"],
            g["mean"] + g["sem"],
            alpha=0.2,
            color=color,
        )

    ax.set_xlabel("Order q", fontsize=12)
    ax.set_ylabel("Hill number D(q)", fontsize=12)
    ax.set_title("Hill Diversity Profile", fontsize=14)
    ax.legend(title=args.group_col, fontsize=10)

    # Annotate special cases
    for q_val, label in [(0, "Richness"), (1, "Shannon"), (2, "Simpson")]:
        if args.q_min <= q_val <= args.q_max:
            ax.axvline(q_val, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)
            ax.text(
                q_val + 0.05,
                ax.get_ylim()[1] * 0.95,
                f"q={q_val}\n({label})",
                fontsize=8,
                color="grey",
                va="top",
            )

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Plot saved to: {args.output}")


if __name__ == "__main__":
    main()
