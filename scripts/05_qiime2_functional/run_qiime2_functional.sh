#!/usr/bin/env bash
#
# QIIME 2 functional analysis pipeline for One Codex pathway, GO, and PFAM
# exports that have been merged and column-normalized to relative frequencies.
#
# Performs:
#   - BIOM conversion of the normalized wide table
#   - Import as a QIIME 2 FeatureTable[RelativeFrequency] artifact
#   - Sample subsetting to those present in the metadata file
#   - Conversion to pseudo-counts (target frequency 1,000,000) via
#     qiime sapienns frequency, to satisfy ANCOM-BC2's input requirements
#   - Filtering of features present in <15 samples or with total counts <100
#   - ANCOM-BC2 differential abundance, adjusted for age and sex, at the
#     relaxed FDR threshold of q<0.10 used for functional analyses
#
# Run separately for the normalized pathway, GO, and PFAM tables.
#
# Requirements: QIIME 2 2025.4 (amplicon distribution) with the
#               q2-sapienns and q2-composition plugins, biom-format.
#
# Usage:
#   bash run_qiime2_functional.sh <normalized_table.tsv> <metadata.tsv> <output_dir>
#
# Example:
#   bash run_qiime2_functional.sh \
#     ../04_functional_prep/pathways_wide_norm.tsv \
#     ../../data/sample_metadata.tsv \
#     ../../results/qiime2_pathways

set -euo pipefail

# --- Arguments ---
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <normalized_table.tsv> <metadata.tsv> <output_dir>" >&2
    exit 1
fi

NORM_TSV="$1"
METADATA="$2"
OUTDIR="$3"

mkdir -p "${OUTDIR}"

# --- Derive a base name for outputs from the input ---
BASE="$(basename "${NORM_TSV}" .tsv)"

# --- File paths for each intermediate artifact ---
BIOM="${OUTDIR}/${BASE}.biom"
TABLE_QZA="${OUTDIR}/${BASE}.qza"
TABLE_META="${OUTDIR}/${BASE}_meta-samples-only.qza"
TABLE_COUNTS="${OUTDIR}/${BASE}_meta-samples-only_counts.qza"
TABLE_FILT="${OUTDIR}/${BASE}_meta-samples-only_counts_rm-low-prev-low-freq.qza"
ANCOMBC2="${OUTDIR}/${BASE}_meta-samples-only_counts_rm-low-prev-low-freq_ancombc2_p0-10.qza"

# --- 1. Convert normalized TSV to BIOM (HDF5) ---
echo "[1/6] Converting TSV to BIOM..."
biom convert \
    --input-fp "${NORM_TSV}" \
    --output-fp "${BIOM}" \
    --table-type "Function table" \
    --to-hdf5

# --- 2. Import BIOM as a QIIME 2 FeatureTable[RelativeFrequency] artifact ---
echo "[2/6] Importing BIOM as QIIME 2 artifact..."
qiime tools import \
    --input-path "${BIOM}" \
    --type 'FeatureTable[RelativeFrequency]' \
    --input-format BIOMV210Format \
    --output-path "${TABLE_QZA}"

# --- 3. Subset samples to those present in the metadata file ---
# Drops samples that failed One Codex QC or had incomplete metadata.
echo "[3/6] Filtering to metadata samples..."
qiime feature-table filter-samples \
    --i-table "${TABLE_QZA}" \
    --m-metadata-file "${METADATA}" \
    --o-filtered-table "${TABLE_META}"

# --- 4. Rescale relative frequencies to pseudo-counts (target 1,000,000) ---
# ANCOM-BC2 requires count-like input; sapienns frequency multiplies the
# relative frequencies by --p-target-freq to produce integer-valued counts.
echo "[4/6] Rescaling to pseudo-counts (target 1,000,000)..."
qiime sapienns frequency \
    --i-table "${TABLE_META}" \
    --p-target-freq 1000000 \
    --o-output-table "${TABLE_COUNTS}"

# --- 5. Filter features present in <15 samples or with total counts <100 ---
echo "[5/6] Filtering low-prevalence and low-abundance features..."
qiime feature-table filter-features \
    --i-table "${TABLE_COUNTS}" \
    --p-min-samples 15 \
    --p-min-frequency 100 \
    --o-filtered-table "${TABLE_FILT}"

# --- 6. ANCOM-BC2 differential abundance (relaxed FDR threshold q<0.10) ---
echo "[6/6] Running ANCOM-BC2..."
qiime composition ancombc2 \
    --i-table "${TABLE_FILT}" \
    --m-metadata-file "${METADATA}" \
    --p-fixed-effects-formula "opioid_status + age + sex" \
    --p-reference-levels opioid_status::Remission sex::Female \
    --p-alpha 0.10 \
    --p-num-processes 8 \
    --o-ancombc2-output "${ANCOMBC2}"

echo "Done. Outputs in: ${OUTDIR}"
