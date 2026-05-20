#!/usr/bin/env bash
#
# QIIME 2 taxonomy pipeline for ConQuR-corrected OneCodex counts.
#
# Performs:
#   - BIOM/QZA conversion of the ConQuR-corrected count table
#   - Feature table transposition (samples as rows)
#   - Singleton filtering (features present in <2 samples removed)
#   - Bray-Curtis beta diversity
#   - PERMANOVA (adonis) adjusted for age and sex
#   - PERMDISP for homogeneity of multivariate dispersions
#   - PCoA for visualization
#   - ANCOM-BC2 differential abundance, adjusted for age and sex
#
# Run separately for genus- and species-level corrected tables.
#
# Requirements: QIIME 2 2025.4 (amplicon distribution), biom-format
#
# Usage:
#   bash run_qiime2_taxonomy.sh <corrected_counts.tsv> <metadata.tsv> <output_dir>
#
# Example:
#   bash run_qiime2_taxonomy.sh \
#     ../01_batch_correction/species_corrected.tsv \
#     ../../data/sample_metadata.tsv \
#     ../../results/qiime2_species

set -euo pipefail

# --- Arguments ---
if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <corrected_counts.tsv> <metadata.tsv> <output_dir>" >&2
    exit 1
fi

CORRECTED_TSV="$1"
METADATA="$2"
OUTDIR="$3"

mkdir -p "${OUTDIR}"

# --- Derive a base name for outputs from the input ---
BASE="$(basename "${CORRECTED_TSV}" .tsv)"

# --- File paths for each intermediate artifact ---
BIOM="${OUTDIR}/${BASE}.biom"
TABLE_QZA="${OUTDIR}/${BASE}.qza"
TABLE_T="${OUTDIR}/${BASE}_transposed.qza"
TABLE_FILT="${OUTDIR}/${BASE}_transposed_no-singletons.qza"
DM="${OUTDIR}/${BASE}_bray-curtis.qza"
ADONIS_QZV="${OUTDIR}/${BASE}_bray-curtis_adonis.qzv"
PERMDISP_QZV="${OUTDIR}/${BASE}_bray-curtis_permdisp.qzv"
PCOA="${OUTDIR}/${BASE}_bray-curtis_pcoa.qza"
ANCOMBC2="${OUTDIR}/${BASE}_ancombc2_p0-05.qza"

# --- 1. Convert TSV to BIOM (HDF5) ---
echo "[1/9] Converting TSV to BIOM..."
biom convert \
    --input-fp "${CORRECTED_TSV}" \
    --output-fp "${BIOM}" \
    --table-type="Taxon table" \
    --to-hdf5

# --- 2. Import BIOM as a QIIME 2 FeatureTable[Frequency] artifact ---
# NOTE: this step is missing from the originally-provided commands but is
# required to produce the .qza consumed downstream. Please confirm the
# input-format flag matches what you used.
echo "[2/9] Importing BIOM as QIIME 2 artifact..."
qiime tools import \
    --input-path "${BIOM}" \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path "${TABLE_QZA}"

# --- 3. Transpose so samples are rows (ConQuR output orientation -> QIIME 2 expected) ---
echo "[3/9] Transposing feature table..."
qiime feature-table transpose \
    --i-table "${TABLE_QZA}" \
    --o-transposed-feature-table "${TABLE_T}"

# --- 4. Filter features present in <2 samples ---
echo "[4/9] Filtering singletons..."
qiime feature-table filter-features \
    --i-table "${TABLE_T}" \
    --p-min-samples 2 \
    --o-filtered-table "${TABLE_FILT}"

# --- 5. Bray-Curtis beta diversity ---
echo "[5/9] Computing Bray-Curtis distances..."
qiime diversity beta \
    --i-table "${TABLE_FILT}" \
    --p-metric braycurtis \
    --o-distance-matrix "${DM}"

# --- 6. PERMANOVA via adonis, adjusting for age and sex ---
echo "[6/9] Running adonis (PERMANOVA)..."
qiime diversity adonis \
    --i-distance-matrix "${DM}" \
    --m-metadata-file "${METADATA}" \
    --p-formula "opioid_status + age + sex" \
    --p-n-jobs auto \
    --o-visualization "${ADONIS_QZV}"

# --- 7. PERMDISP for dispersion homogeneity ---
echo "[7/9] Running PERMDISP..."
qiime diversity beta-group-significance \
    --i-distance-matrix "${DM}" \
    --m-metadata-file "${METADATA}" \
    --m-metadata-column opioid_status \
    --p-method permdisp \
    --p-pairwise \
    --o-visualization "${PERMDISP_QZV}"

# --- 8. PCoA ---
echo "[8/9] Running PCoA..."
qiime diversity pcoa \
    --i-distance-matrix "${DM}" \
    --o-pcoa "${PCOA}"

# --- 9. ANCOM-BC2 differential abundance ---
echo "[9/9] Running ANCOM-BC2..."
qiime composition ancombc2 \
    --i-table "${TABLE_FILT}" \
    --m-metadata-file "${METADATA}" \
    --p-fixed-effects-formula "opioid_status + age + sex" \
    --p-reference-levels opioid_status::Remission sex::Female \
    --p-alpha 0.05 \
    --p-num-processes 8 \
    --o-ancombc2-output "${ANCOMBC2}"

echo "Done. Outputs in: ${OUTDIR}"
