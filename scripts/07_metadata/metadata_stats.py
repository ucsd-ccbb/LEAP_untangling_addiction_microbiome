#!/usr/bin/env python3
"""Test demographic balance between groups (Table 1 stats).

For continuous variables: reports mean, SD, median, IQR per group and tests
group differences with the Wilcoxon rank-sum (Mann-Whitney U) test.

For categorical variables: reports a contingency table and tests group
differences with Fisher's exact (for 2x2 tables after rare-category
filtering) or chi-squared (otherwise).

Rare categorical levels can be excluded from the test via --min-category-n
to avoid low-expected-count cells in 2x3+ tables. Excluded samples are
reported in the console output but not removed from any other analysis.

p-values are reported uncorrected; for Table 1 balance tests, multiple
testing correction is conventionally not applied.
"""

import argparse
import os
import pandas as pd
from scipy.stats import mannwhitneyu, fisher_exact, chi2_contingency


def main():
    parser = argparse.ArgumentParser(
        description="Test demographic balance between groups."
    )
    parser.add_argument("--input", required=True, help="Path to metadata TSV file.")
    parser.add_argument("--output", required=True, help="Path to output TSV file.")
    parser.add_argument(
        "--group-col",
        default="opioid_status",
        help="Column defining the two groups (default: opioid_status).",
    )
    parser.add_argument(
        "--continuous-vars",
        default="age",
        help="Comma-delimited list of continuous variables to test (default: age).",
    )
    parser.add_argument(
        "--categorical-vars",
        default="sex",
        help="Comma-delimited list of categorical variables to test (default: sex).",
    )
    parser.add_argument(
        "--min-category-n",
        type=int,
        default=2,
        help=(
            "For categorical variables, exclude levels with total n below this "
            "threshold before testing (default: 2). Excluded levels are reported "
            "in the console output."
        ),
    )
    args = parser.parse_args()

    continuous_vars = [v.strip() for v in args.continuous_vars.split(",") if v.strip()]
    categorical_vars = [
        v.strip() for v in args.categorical_vars.split(",") if v.strip()
    ]

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    meta = pd.read_csv(args.input, sep="\t")
    meta = meta.dropna(subset=[args.group_col])

    # Keep only relevant columns
    keep_cols = [args.group_col] + continuous_vars + categorical_vars
    missing = [c for c in keep_cols if c not in meta.columns]
    if missing:
        raise ValueError(f"Columns not found in metadata: {missing}")
    meta = meta[keep_cols]

    groups = meta[args.group_col].unique()
    if len(groups) != 2:
        raise ValueError(
            f"Expected 2 groups in '{args.group_col}', found {len(groups)}: {groups}"
        )

    g1_label, g2_label = sorted(groups)
    g1 = meta[meta[args.group_col] == g1_label]
    g2 = meta[meta[args.group_col] == g2_label]

    print(f"\n=== Groups ({args.group_col}) ===")
    print(f"  {g1_label}: n={len(g1)}")
    print(f"  {g2_label}: n={len(g2)}")

    results = []

    for var in continuous_vars:
        v1 = g1[var].dropna()
        v2 = g2[var].dropna()

        print(f"\n=== {var} ===")
        for label, vals in [(g1_label, v1), (g2_label, v2)]:
            print(
                f"  {label}: mean={vals.mean():.2f}, sd={vals.std():.2f}, "
                f"median={vals.median():.2f}, IQR=[{vals.quantile(0.25):.2f}, {vals.quantile(0.75):.2f}]"
            )

        stat, p = mannwhitneyu(v1, v2, alternative="two-sided")
        print(f"  Wilcoxon rank-sum: U={stat:.1f}, p={p:.4g}")

        results.append(
            {
                "variable": var,
                "test": "Wilcoxon rank-sum",
                "statistic": stat,
                "p_value": p,
                "excluded_levels": "",
            }
        )

    for var in categorical_vars:
        print(f"\n=== {var} ===")

        # Full contingency table (for display)
        ct_full = pd.crosstab(meta[args.group_col], meta[var])
        print("Full contingency table:")
        print(ct_full.to_string())

        # Identify rare levels by total n across both groups
        level_totals = ct_full.sum(axis=0)
        rare_levels = level_totals[level_totals < args.min_category_n].index.tolist()

        if rare_levels:
            print(
                f"  Excluding rare level(s) with total n < {args.min_category_n}: "
                f"{rare_levels}"
            )
            ct = ct_full.drop(columns=rare_levels)
        else:
            ct = ct_full

        if ct.shape[1] < 2:
            print(
                f"  WARNING: fewer than 2 levels remain for {var} after filtering; "
                "skipping test."
            )
            results.append(
                {
                    "variable": var,
                    "test": "skipped (insufficient levels)",
                    "statistic": float("nan"),
                    "p_value": float("nan"),
                    "excluded_levels": ",".join(map(str, rare_levels)),
                }
            )
            continue

        if ct.shape == (2, 2):
            odds_ratio, p = fisher_exact(ct)
            test_name = "Fisher's exact"
            stat = odds_ratio
            print(f"  Fisher's exact: OR={odds_ratio:.3f}, p={p:.4g}")
        else:
            chi2, p, dof, _ = chi2_contingency(ct)
            test_name = "Chi-squared"
            stat = chi2
            print(f"  Chi-squared: chi2={chi2:.3f}, df={dof}, p={p:.4g}")

        results.append(
            {
                "variable": var,
                "test": test_name,
                "statistic": stat,
                "p_value": p,
                "excluded_levels": ",".join(map(str, rare_levels)),
            }
        )

    results_df = pd.DataFrame(results)
    results_df.to_csv(args.output, sep="\t", index=False)
    print(f"\nResults written to: {args.output}")


if __name__ == "__main__":
    main()
