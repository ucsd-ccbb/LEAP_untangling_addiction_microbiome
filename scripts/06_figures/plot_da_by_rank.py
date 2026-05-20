#!/usr/bin/env python3
"""
Plot differentially abundant taxa grouped by a higher taxonomic rank.

Reads ANCOMBC2 output as either:
  (a) a QIIME 2 .qza artifact containing lfc.jsonl and q.jsonl, or
  (b) separate LFC and q-value CSV files (legacy format).

Produces a stacked bar chart showing the number of significant taxa
enriched or depleted per grouping rank (e.g. species counted per family).

Examples
--------
# From a .qza artifact (QIIME 2 ANCOMBC2 output):
python plot_da_by_rank.py \
    --qza results/ancombc.qza \
    --contrast "opioid_status::Remission" \
    --top-n 15 \
    -o da_species_by_family.svg

# From legacy CSV pair:
python plot_da_by_rank.py \
    --lfc DA_species_lfc.csv \
    --qval DA_species_q_vals.csv \
    --contrast opioid_statusoud \
    -o da_species_by_family.svg
"""

import argparse
import json
import sys
import tempfile
import zipfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── Taxonomy helpers ────────────────────────────────────────────────
RANK_ORDER = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]


def parse_taxonomy(tax_str: str) -> pd.Series:
    """Split a semicolon-delimited taxonomy string into named ranks."""
    parts = str(tax_str).split(";")

    def safe(i):
        return (
            parts[i] if len(parts) > i and parts[i] not in ("", "nan") else "Unassigned"
        )

    return pd.Series({rank: safe(i) for i, rank in enumerate(RANK_ORDER)})


# ── JSONL loader ────────────────────────────────────────────────────
def read_ancombc_jsonl(path: str) -> pd.DataFrame:
    """Read an ANCOMBC2 table.jsonl file, skipping the schema header line."""
    records = []
    with open(path) as fh:
        for line in fh:
            obj = json.loads(line)
            if "doctype" in obj:
                continue
            records.append(obj)
    df = pd.DataFrame(records)
    if "taxon" in df.columns:
        df.index = df["taxon"]
        df.drop(columns=["taxon"], inplace=True)
    elif "id" in df.columns:
        df.index = df["id"]
        df.drop(columns=["id"], inplace=True)
    return df


def load_from_qza(qza_path: str, contrast: str):
    """Extract lfc.jsonl, q.jsonl, and passed_ss.jsonl from a .qza and return a merged DA DataFrame."""
    with tempfile.TemporaryDirectory() as tmpdir:
        with zipfile.ZipFile(qza_path) as zf:
            names = zf.namelist()
            lfc_member = [n for n in names if n.endswith("data/lfc.jsonl")]
            q_member = [n for n in names if n.endswith("data/q.jsonl")]
            ss_member = [n for n in names if n.endswith("data/passed_ss.jsonl")]
            if not lfc_member or not q_member:
                sys.exit(
                    "ERROR: Could not find data/lfc.jsonl and data/q.jsonl inside the .qza."
                )
            for m in lfc_member + q_member + ss_member:
                zf.extract(m, tmpdir)
            lfc_path = f"{tmpdir}/{lfc_member[0]}"
            q_path = f"{tmpdir}/{q_member[0]}"
            ss_path = f"{tmpdir}/{ss_member[0]}" if ss_member else None

        lfc_df = read_ancombc_jsonl(lfc_path)
        q_df = read_ancombc_jsonl(q_path)
        ss_df = read_ancombc_jsonl(ss_path) if ss_path else None

    if contrast not in lfc_df.columns:
        sys.exit(
            f"ERROR: contrast '{contrast}' not found in lfc.jsonl.\n"
            f"Available columns: {list(lfc_df.columns)}"
        )
    if contrast not in q_df.columns:
        sys.exit(
            f"ERROR: contrast '{contrast}' not found in q.jsonl.\n"
            f"Available columns: {list(q_df.columns)}"
        )

    da = pd.DataFrame(
        {
            "qval": q_df[contrast],
            "lfc": lfc_df[contrast],
        }
    )

    if ss_df is not None and contrast in ss_df.columns:
        passed = ss_df[contrast].astype(bool)
        n_before = len(da)
        da = da.loc[da.index.isin(passed[passed].index)]
        n_filtered = n_before - len(da)
        if n_filtered > 0:
            print(
                f"Sensitivity filter: removed {n_filtered} taxa that failed passed_ss for '{contrast}'."
            )
    elif ss_df is not None:
        print(
            f"WARNING: contrast '{contrast}' not found in passed_ss.jsonl; skipping sensitivity filter."
        )

    return da.sort_values("qval")


def load_from_csvs(lfc_path: str, qval_path: str, contrast: str):
    """Load DA results from the legacy separate-CSV format."""
    lfc_df = pd.read_csv(lfc_path)
    lfc_df.index = lfc_df["id"]
    if contrast not in lfc_df.columns:
        sys.exit(
            f"ERROR: contrast '{contrast}' not found in LFC file.\n"
            f"Available: {list(lfc_df.columns)}"
        )

    qval_df = pd.read_csv(qval_path)
    qval_df.index = qval_df["id"]
    if contrast not in qval_df.columns:
        sys.exit(
            f"ERROR: contrast '{contrast}' not found in q-value file.\n"
            f"Available: {list(qval_df.columns)}"
        )

    da = pd.DataFrame(
        {
            "qval": qval_df[contrast],
            "lfc": lfc_df[contrast],
        }
    )
    return da.sort_values("qval")


# ── CLI ─────────────────────────────────────────────────────────────
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Plot differentially abundant taxa grouped by a higher taxonomic rank.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    inp = p.add_argument_group("input (provide --qza OR both --lfc and --qval)")
    inp.add_argument(
        "--qza",
        help="QIIME 2 .qza artifact from ANCOMBC2 (contains lfc.jsonl + q.jsonl).",
    )
    inp.add_argument(
        "--lfc", help="CSV with log-fold-change values (legacy ANCOMBC export)."
    )
    inp.add_argument("--qval", help="CSV with q-values (legacy ANCOMBC export).")

    p.add_argument(
        "--contrast",
        required=True,
        help="Column name for the contrast of interest (e.g. 'opioid_status::Remission').",
    )
    p.add_argument(
        "--alpha", type=float, default=0.05, help="Significance threshold for q-values."
    )
    p.add_argument(
        "--group-rank",
        default="Family",
        choices=RANK_ORDER,
        help="Taxonomic rank to aggregate counts by.",
    )
    p.add_argument("--top-n", type=int, default=12, help="Number of groups to display.")
    p.add_argument(
        "--enriched-label",
        default="Enriched in active OUD",
        help="Legend label for LFC > 0.",
    )
    p.add_argument(
        "--depleted-label",
        default="Depleted in active OUD",
        help="Legend label for LFC < 0.",
    )
    p.add_argument(
        "--enriched-color", default="#D3290B", help="Bar color for enriched taxa."
    )
    p.add_argument(
        "--depleted-color", default="#006EB5", help="Bar color for depleted taxa."
    )
    p.add_argument(
        "--title", default=None, help="Plot title. Auto-generated if omitted."
    )
    p.add_argument(
        "--ylabel", default=None, help="Y-axis label. Auto-generated if omitted."
    )
    p.add_argument(
        "--bar-width-inches", type=float, default=0.5, help="Width per bar in inches."
    )
    p.add_argument(
        "--fig-height", type=float, default=6.0, help="Figure height in inches."
    )
    p.add_argument(
        "--dpi", type=int, default=300, help="Resolution for raster outputs."
    )
    p.add_argument(
        "-o",
        "--output",
        default="da_by_rank.svg",
        help="Output file path (svg/png/pdf).",
    )
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    # ── Validate input mode ─────────────────────────────────────────
    has_qza = args.qza is not None
    has_csv = args.lfc is not None or args.qval is not None

    if has_qza and has_csv:
        sys.exit("ERROR: Provide either --qza or --lfc/--qval, not both.")
    if not has_qza and not has_csv:
        sys.exit("ERROR: Provide either --qza or both --lfc and --qval.")
    if has_csv and (args.lfc is None or args.qval is None):
        sys.exit("ERROR: Both --lfc and --qval are required for CSV input.")

    # ── Load data ───────────────────────────────────────────────────
    if has_qza:
        da = load_from_qza(args.qza, args.contrast)
    else:
        da = load_from_csvs(args.lfc, args.qval, args.contrast)

    # ── Filter to significant ───────────────────────────────────────
    da_sig = da[da["qval"] < args.alpha].copy()
    if da_sig.empty:
        sys.exit(f"No significant taxa at alpha={args.alpha}.")

    # ── Parse taxonomy & assign direction ───────────────────────────
    tax = da_sig.index.to_series().apply(parse_taxonomy)
    da_sig = pd.concat([da_sig, tax], axis=1)
    da_sig["direction"] = np.where(
        da_sig["lfc"] > 0, args.enriched_label, args.depleted_label
    )

    # ── Count by group rank ─────────────────────────────────────────
    group_col = args.group_rank
    counts = da_sig.groupby([group_col, "direction"]).size().unstack(fill_value=0)
    sort_col = (
        args.enriched_label
        if args.enriched_label in counts.columns
        else counts.columns[0]
    )
    counts = counts.sort_values(by=sort_col, ascending=False)
    counts_top = counts.head(args.top_n)

    # ── Infer feature rank from deepest non-empty taxonomy level ────
    def _deepest_rank(row):
        for rank in reversed(RANK_ORDER):
            if row.get(rank, "Unassigned") != "Unassigned":
                return rank
        return RANK_ORDER[-1]

    feature_rank = da_sig.apply(_deepest_rank, axis=1).mode().iloc[0]

    # ── Plot ────────────────────────────────────────────────────────
    n_bars = len(counts_top)
    fig_width = max(4.0, n_bars * args.bar_width_inches + 1.5)

    fig, ax = plt.subplots(figsize=(fig_width, args.fig_height))
    y = np.arange(n_bars)

    enriched = counts_top.get(args.enriched_label, pd.Series(0, index=counts_top.index))
    depleted = counts_top.get(args.depleted_label, pd.Series(0, index=counts_top.index))

    ax.bar(y, enriched, label=args.enriched_label, color=args.enriched_color)
    ax.bar(
        y,
        depleted,
        bottom=enriched,
        label=args.depleted_label,
        color=args.depleted_color,
    )

    ax.set_xticks(y)
    ax.set_xticklabels(counts_top.index, rotation=45, ha="right")
    ax.set_ylabel(
        args.ylabel or f"Number of significant {feature_rank.lower()}", fontsize=14
    )
    ax.set_title(
        args.title
        or f"Differentially abundant {feature_rank.lower()} by {group_col.lower()}"
    )
    ax.legend(frameon=False)
    plt.tight_layout()
    fig.savefig(args.output, dpi=args.dpi)
    plt.close(fig)

    taxa_shown = int(counts_top.sum().sum())
    print(
        f"Saved → {args.output}  ({n_bars} {group_col.lower()} groups, "
        f"{taxa_shown} {feature_rank.lower()} shown, "
        f"{len(da_sig)} significant total)"
    )


if __name__ == "__main__":
    main()
