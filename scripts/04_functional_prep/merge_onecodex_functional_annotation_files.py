#!/usr/bin/env python3
"""
Merge One Codex per-sample functional annotation files into a wide
features-by-samples table.

Works for any One Codex functional export with the same general shape
(pathways, GO terms, PFAM families), where each per-sample file has:
  - a row-ID column (e.g., "Function ID")
  - a numeric value column (e.g., "Abundance" or "CPM")
  - optionally a row-name column (e.g., "Function Name")
  - optionally a sample-identifier column (e.g., "Sample ID" or "Sample Name")

Outputs:
  - A wide TSV with row IDs as rows and samples as columns.
  - Optionally, a 2-column TSV mapping row IDs to row names.

Example invocations:

  # MetaCyc pathways (per-sample value column is "Abundance")
  python merge_onecodex_functional_annotation_files.py \\
    --input-glob "One_Codex_Pathways_Unstratified/*.csv" \\
    --value-col Abundance \\
    --output-tsv pathways_wide.tsv \\
    --output-map pathway_function_map.tsv

  # GO terms or PFAM families (per-sample value column is "CPM")
  python merge_onecodex_functional_annotation_files.py \\
    --input-glob "One_Codex_GO_Unstratified/*.csv" \\
    --value-col CPM \\
    --output-tsv go_wide.tsv \\
    --output-map go_function_map.tsv

  python merge_onecodex_functional_annotation_files.py \\
    --input-glob "One_Codex_PFAM_Unstratified/*.csv" \\
    --value-col CPM \\
    --output-tsv pfam_wide.tsv \\
    --output-map pfam_function_map.tsv
"""

import argparse
import csv
import glob
import os

import pandas as pd


def sniff_delimiter(fp, default=","):
    """Detect CSV/TSV/semicolon delimiter for a single file."""
    with open(fp, "r", encoding="utf-8-sig", newline="") as f:
        sample = f.read(4096)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";"])
        return dialect.delimiter
    except Exception:
        return default


def pick_column_label(df, field, fp):
    """Choose a sample column label from a metadata field, with fallbacks."""
    if field.lower() == "filename":
        return os.path.splitext(os.path.basename(fp))[0]
    if field in df.columns and df[field].notna().any():
        return str(df[field].iloc[0])
    # Graceful fallback: try common fields, else filename
    for alt in ("Sample ID", "Sample Name"):
        if alt in df.columns and df[alt].notna().any():
            return str(df[alt].iloc[0])
    return os.path.splitext(os.path.basename(fp))[0]


def read_onecodex_function_file(
    fp, value_col, row_id_col, row_name_col, column_label_field
):
    """Read one One Codex per-sample functional file, return (series, name_map, label)."""
    delim = sniff_delimiter(fp)
    df = pd.read_csv(fp, sep=delim, engine="c", encoding="utf-8-sig")

    if row_id_col not in df.columns:
        raise ValueError(
            f"'{row_id_col}' not found in {fp}. Columns: {list(df.columns)}"
        )
    if value_col not in df.columns:
        raise ValueError(
            f"'{value_col}' not found in {fp}. Columns: {list(df.columns)}"
        )

    sample_label = pick_column_label(df, column_label_field, fp)
    sample_label = sample_label.replace("/", "_").replace(" ", "_")

    cols = [row_id_col, value_col]
    has_name = row_name_col in df.columns
    if has_name:
        cols.append(row_name_col)

    # Sum duplicate IDs per file (safety)
    agg_spec = {value_col: "sum"}
    if has_name:
        agg_spec[row_name_col] = "first"

    g = (
        df[cols]
        .groupby(row_id_col, as_index=True)
        .agg(agg_spec)
        .rename(columns={value_col: sample_label})
    )

    series = g[[sample_label]]
    name_map = g[row_name_col] if has_name else None
    return series, name_map, sample_label


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Merge One Codex per-sample functional files (pathways, GO, PFAM) "
            "into a wide features-by-samples table."
        )
    )
    ap.add_argument(
        "--input-glob",
        required=True,
        help="Glob for per-sample files, e.g. 'One_Codex_Pathways_Unstratified/*.csv'",
    )
    ap.add_argument(
        "--value-col",
        required=True,
        help=(
            "Numeric column to merge. Use 'Abundance' for One Codex pathway "
            "exports, 'CPM' for GO and PFAM exports."
        ),
    )
    ap.add_argument(
        "--row-id-col",
        default="Function ID",
        help="Column to use as row IDs in the wide table (default: 'Function ID')",
    )
    ap.add_argument(
        "--row-name-col",
        default="Function Name",
        help=(
            "Column used to populate the ID->name mapping in --output-map "
            "(default: 'Function Name')."
        ),
    )
    ap.add_argument(
        "--column-label-field",
        default="Sample ID",
        help=(
            "Field used to name sample columns: 'Sample ID', 'Sample Name', "
            "or 'filename' (default: 'Sample ID')."
        ),
    )
    ap.add_argument(
        "--output-tsv",
        required=True,
        help="Output wide TSV (rows = row-id-col, cols = samples).",
    )
    ap.add_argument(
        "--output-map",
        default=None,
        help="Optional 2-col TSV mapping (row-id-col, row-name-col).",
    )
    args = ap.parse_args()

    files = sorted(glob.glob(args.input_glob))
    if not files:
        raise SystemExit(f"No files matched: {args.input_glob}")

    series_list, name_maps = [], []
    for fp in files:
        s, nm, _ = read_onecodex_function_file(
            fp,
            value_col=args.value_col,
            row_id_col=args.row_id_col,
            row_name_col=args.row_name_col,
            column_label_field=args.column_label_field,
        )
        series_list.append(s)
        if nm is not None:
            name_maps.append(nm)

    wide = pd.concat(series_list, axis=1).fillna(0)
    wide = wide.sort_index()
    wide = wide[sorted(wide.columns)]
    wide.to_csv(args.output_tsv, sep="\t")

    if args.output_map:
        if name_maps:
            name_df = pd.concat(name_maps, axis=1)
            name_ser = name_df.bfill(axis=1).iloc[:, 0]
            out = pd.DataFrame(
                {args.row_id_col: name_ser.index, args.row_name_col: name_ser.values}
            )
        else:
            out = pd.DataFrame(
                {args.row_id_col: wide.index, args.row_name_col: [""] * len(wide)}
            )
        out.to_csv(args.output_map, sep="\t", index=False)


if __name__ == "__main__":
    main()
