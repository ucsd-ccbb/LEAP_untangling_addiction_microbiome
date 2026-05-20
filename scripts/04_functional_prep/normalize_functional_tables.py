#!/usr/bin/env python3
"""
Convert a wide functional abundance table (features x samples) to
column-wise relative frequencies (each sample sums to 1).

Used to prepare One Codex MetaCyc pathway, GO term, and PFAM family wide
tables (output of merge_onecodex_functional_annotation_files.py) for
downstream import into QIIME 2 as relative frequency tables.

Samples with a total of 0 are passed through as zeros (no division by zero).

Example invocations:

  python normalize_functional_tables.py \\
    --input  pathways_wide.tsv \\
    --output pathways_wide_norm.tsv

  python normalize_functional_tables.py \\
    --input  GO_cpm_wide.tsv \\
    --output GO_cpm_wide_norm.tsv

  python normalize_functional_tables.py \\
    --input  PFAM_cpm_wide.tsv \\
    --output PFAM_cpm_wide_norm.tsv
"""

import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Normalize a wide features-by-samples table so each sample "
            "column sums to 1."
        )
    )
    ap.add_argument(
        "--input",
        required=True,
        help="Wide TSV from merge_onecodex_functional_annotation_files.py "
        "(features as rows, samples as columns).",
    )
    ap.add_argument(
        "--output",
        required=True,
        help="Output TSV with column-wise relative frequencies.",
    )
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t", index_col=0)
    col_sums = df.sum(axis=0)
    df_rel = df.div(col_sums.where(col_sums > 0), axis=1)
    df_rel.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    main()
