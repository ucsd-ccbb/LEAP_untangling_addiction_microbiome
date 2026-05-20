# LEAP: Untangling Addiction Microbiome

Code accompanying:

> Roy C-A, Rosenthal SB, Sathish V, Palugod B, Hernandez A, Marienfeld C, Gabriel RA.
> *Shotgun metagenomics identifies gut microbiome alterations associated with active opioid use disorder.*
> [Journal TBD] (2026). DOI: TBD.

This repository contains the analysis scripts used to produce the diversity,
differential abundance, and functional results reported in the manuscript.

## Repository contents

```
LEAP_untangling_addiction_microbiome/
├── README.md
├── CITATION.cff                # Citation metadata
├── environment.yml             # Python + QIIME 2 2025.4 conda env
├── renv.lock                   # R package versions for ConQuR step
├── scripts/
│   ├── 01_batch_correction/    # ConQuR batch correction (R)
│   ├── 02_qiime2_taxonomy/     # QIIME 2 taxonomic pipeline
│   ├── 03_alpha_diversity/     # Alpha diversity metrics + testing + plots
│   ├── 04_functional_prep/     # Merging and normalization of One Codex
│   │                           # pathway, GO, and PFAM exports
│   ├── 05_qiime2_functional/   # QIIME 2 functional pipeline
│   ├── 06_figures/             # Figure-generation scripts
│   └── 07_metadata/            # Demographic balance testing
└── figures/                    # Final figure outputs
```

## Software requirements

- **QIIME 2 2025.4** (amplicon distribution), with the `q2-composition` and
  `q2-sapienns` plugins installed. Install QIIME 2 following the official
  instructions (<https://docs.qiime2.org/>), then install `q2-sapienns` into
  the same environment.
- **Python 3.10+** with the packages listed in `environment.yml`
  (`pandas`, `numpy`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`,
  `biom-format`).
- **R 4.3+** with the packages captured in `renv.lock`. Restore the R
  environment with `renv::restore()` from the project directory.
  `ConQuR` is installed from GitHub (`wdl2459/ConQuR`).

To restore the Python environment:

```bash
conda env create -f environment.yml
conda activate qiime2-amplicon-2025.4
```

To restore the R environment:

```r
install.packages("renv")
renv::restore()
```

## Data availability

Sequencing data, sample metadata, and One Codex outputs are not included in
this repository. Access information will be provided here upon publication:

- **Raw sequencing reads:** [SRA/dbGaP accession TBD]
- **One Codex taxonomic and functional exports:** [accession or contact TBD]
- **Sample metadata:** [access mechanism TBD; subject to IRB restrictions]

The Gabriel cohort consists of human subjects enrolled under UC San Diego
IRB protocol #[number]. Data sharing is governed by the consent and IRB
terms approved for this study. Researchers interested in accessing the
data should contact the corresponding author.

## Analysis overview

The analysis pipeline proceeds in seven stages:

1. **Initial data preparation** (manual, not scripted). Raw OneCodex
   per-sample counts were downloaded in batches and combined in a
   spreadsheet, sample identifiers were harmonized, and the metadata
   columns `batch`, `opioid_status`, `age`, and `sex` were appended.
   The output is two TSVs with samples as rows, taxa as columns (one each
   at genus and species level), and trailing metadata columns.
2. **Batch correction** with ConQuR (R), adjusting for `opioid_status`,
   `age`, and `sex`.
3. **Taxonomic analysis** in QIIME 2: import, transposition, singleton
   filtering, Bray-Curtis beta diversity, PERMANOVA (adonis), PERMDISP,
   PCoA, and ANCOM-BC2 differential abundance.
4. **Alpha diversity analysis** in Python: Shannon, Simpson, and observed
   features, with OLS testing controlling for age and sex, and generation of
   violin plots.
5. **Functional preparation**: per-sample One Codex pathway, GO, and PFAM
   exports are merged into wide tables, then column-normalized to
   relative frequencies.
6. **Functional QIIME 2 pipeline**: import, sample subsetting, rescaling
   to pseudo-counts, low-prevalence/low-abundance feature filtering, and
   ANCOM-BC2 differential abundance at a relaxed FDR threshold (q<0.10).
7. **Figure generation** for all manuscript figures and supplemental
   figures.

## How to reproduce

The commands below assume the conda env is active, the R environment is
restored, and the input data files are placed in a `data/` directory at
the project root.

### Stage 2: ConQuR batch correction

```bash
Rscript scripts/01_batch_correction/conqur_batch_correction.R
```

The script reads `data/[input].tsv`, runs ConQuR with `opioid_status + age + sex`
as covariates and batch `0` as the reference, and writes a corrected TSV.
Run separately for the genus- and species-level input tables (edit the
input/output paths at the top of the script).

### Stage 3: QIIME 2 taxonomic pipeline

```bash
bash scripts/02_qiime2_taxonomy/run_qiime2_taxonomy.sh \
  results/01_batch_correction/genus_corrected.tsv \
  data/sample_metadata.tsv \
  results/02_qiime2_genus

bash scripts/02_qiime2_taxonomy/run_qiime2_taxonomy.sh \
  results/01_batch_correction/species_corrected.tsv \
  data/sample_metadata.tsv \
  results/02_qiime2_species
```

### Stage 4: Alpha diversity

The alpha diversity pipeline operates on the post-filter, TSV-exported
feature table from the QIIME 2 stage (Stage 3).

```bash
python scripts/03_alpha_diversity/alpha_diversity_pipeline.py \
  --input results/02_qiime2_genus/feature-table.tsv \
  --metadata data/sample_metadata.tsv \
  --output-dir results/03_alpha_diversity_genus \
  --seed 1

python scripts/03_alpha_diversity/alpha_diversity_pipeline.py \
  --input results/02_qiime2_species/feature-table.tsv \
  --metadata data/sample_metadata.tsv \
  --output-dir results/03_alpha_diversity_species \
  --seed 1
```

### Stage 5: Functional preparation

```bash
# Merge per-sample One Codex exports
python scripts/04_functional_prep/merge_onecodex_functional_annotation_files.py \
  --input-glob "data/One_Codex_Pathways_Unstratified/*.csv" \
  --value-col Abundance \
  --output-tsv results/04_functional_prep/pathways_wide.tsv \
  --output-map results/04_functional_prep/pathway_function_map.tsv

python scripts/04_functional_prep/merge_onecodex_functional_annotation_files.py \
  --input-glob "data/One_Codex_GO_Unstratified/*.csv" \
  --value-col CPM \
  --output-tsv results/04_functional_prep/GO_cpm_wide.tsv \
  --output-map results/04_functional_prep/GO_function_map.tsv

python scripts/04_functional_prep/merge_onecodex_functional_annotation_files.py \
  --input-glob "data/One_Codex_PFAM_Unstratified/*.csv" \
  --value-col CPM \
  --output-tsv results/04_functional_prep/PFAM_cpm_wide.tsv \
  --output-map results/04_functional_prep/PFAM_function_map.tsv

# Normalize each to column-wise relative frequencies
python scripts/04_functional_prep/normalize_functional_tables.py \
  --input results/04_functional_prep/pathways_wide.tsv \
  --output results/04_functional_prep/pathways_wide_norm.tsv

python scripts/04_functional_prep/normalize_functional_tables.py \
  --input results/04_functional_prep/GO_cpm_wide.tsv \
  --output results/04_functional_prep/GO_cpm_wide_norm.tsv

python scripts/04_functional_prep/normalize_functional_tables.py \
  --input results/04_functional_prep/PFAM_cpm_wide.tsv \
  --output results/04_functional_prep/PFAM_cpm_wide_norm.tsv
```

### Stage 6: QIIME 2 functional pipeline

```bash
bash scripts/05_qiime2_functional/run_qiime2_functional.sh \
  results/04_functional_prep/pathways_wide_norm.tsv \
  data/sample_metadata.tsv \
  results/05_qiime2_pathways

bash scripts/05_qiime2_functional/run_qiime2_functional.sh \
  results/04_functional_prep/GO_cpm_wide_norm.tsv \
  data/sample_metadata.tsv \
  results/05_qiime2_go

bash scripts/05_qiime2_functional/run_qiime2_functional.sh \
  results/04_functional_prep/PFAM_cpm_wide_norm.tsv \
  data/sample_metadata.tsv \
  results/05_qiime2_pfam
```

### Stage 7: Figures

Each command corresponds to one or more panels in the manuscript figures.
Parameter values are those used to produce the published figures.

**Figure 1A / Supplemental Figure 1A,B: Stacked relative abundance bar plots**

```bash
# Per-group, top 30 taxa (Figure 1A, genus shown)
python scripts/06_figures/plot_grouped_microbiome_stacked.py \
  --abundance_tsv data/onecodex_metagenomics_genus.tsv \
  --metadata_tsv  data/sample_metadata.tsv \
  --group_col opioid_status \
  --sample_col sample_name \
  --top_n 30 \
  --group_order Remission,OUD \
  --output figures/grouped_stacked_barplot_genus.pdf \
  --title "" \
  --axis_label_fontsize 20 --tick_label_fontsize 20 --legend_fontsize 12 \
  --bar_width 0.40 --group_spacing 0.60

# Per-sample (Supplemental Figure 1)
python scripts/06_figures/plot_grouped_samples_stacked.py \
  --abundance_tsv data/onecodex_metagenomics_genus.tsv \
  --metadata_tsv  data/sample_metadata.tsv \
  --group_col opioid_status \
  --sample_col sample_name \
  --top_n 30 \
  --group_order Remission,OUD \
  --switch_threshold 0.05 \
  --output figures/grouped_stacked_sample_barplot_genus.pdf \
  --title "" \
  --axis_label_fontsize 32 --tick_label_fontsize 28 \
  --legend_fontsize 18 --title_fontsize 26 --group_label_fontsize 28
```

**Figure 1B,C / Supplemental Figure 2: Alpha diversity violin plots**

Generated by `alpha_diversity_pipeline.py` in Stage 4.

**Figure 1D,E: Bray-Curtis PCoA**

```bash
python scripts/06_figures/plot_qiime_pcoa.py \
  --pcoa     results/02_qiime2_genus/[base]_bray-curtis_pcoa.qza \
  --metadata data/sample_metadata.tsv \
  --group-col opioid_status \
  --sample-id-col sample_name \
  --output figures/bray_pcoa_genus.pdf \
  --title "" \
  --axis-label-size 28 --tick-label-size 28 --legend-text-size 20 \
  --group-color OUD=#0C7BDC \
  --group-color Remission=#FFC20A
```

**Figure 2A: Top differentially abundant species**

```bash
python scripts/06_figures/plot_da.py \
  --lfc  results/02_qiime2_species/ancombc2_export/lfc_slice.csv \
  --qval results/02_qiime2_species/ancombc2_export/q_val_slice.csv \
  --q-threshold 0.2 \
  --top-n 10 \
  --taxonomy \
  --show-legend --show-stars \
  --column opioid_statusoud \
  --enriched-color "#c23b22" \
  --depleted-color "#2b6cb0" \
  --xlabel "ANCOM-BC log fold change" \
  --title "" --font-size 16 \
  -o figures/DA_species.pdf
```

**Figure 2B: Differentially abundant species by family**

```bash
python scripts/06_figures/plot_da_by_rank.py \
  --qza results/02_qiime2_species/[base]_ancombc2_p0-05.qza \
  --contrast "opioid_status::OUD" \
  --top-n 20 \
  --group-rank Family \
  --enriched-color "#c23b22" \
  --depleted-color "#2b6cb0" \
  -o figures/da_species_by_family.pdf
```

**Figure 3: Differentially abundant GO terms, pathways, and protein families**

`plot_da.py` is reused for each functional category (pathways, GO, PFAM)
with inputs taken from the corresponding Stage 6 ANCOM-BC2 exports.

**Supplemental Figure 3: Hill diversity profiles**

```bash
python scripts/06_figures/plot_hill_diversity_curves.py \
  --input    results/02_qiime2_genus/feature-table.tsv \
  --metadata data/sample_metadata.tsv \
  --output   figures/hill_diversity_genus.pdf

python scripts/06_figures/plot_hill_diversity_curves.py \
  --input    results/02_qiime2_species/feature-table.tsv \
  --metadata data/sample_metadata.tsv \
  --output   figures/hill_diversity_species.pdf
```

### Demographic balance testing (Table 1)

```bash
python scripts/07_metadata/metadata_stats.py \
  --input  data/sample_metadata.tsv \
  --output results/07_metadata/balance_tests.tsv \
  --group-col opioid_status \
  --continuous-vars age \
  --categorical-vars sex
```

## Use of AI-assisted tools

Several scripts in this repository were polished with the assistance of
large language models (Anthropic Claude). All AI-assisted code was reviewed,
tested against the study data, and edited by the authors, who take full
responsibility for its correctness. AI tools were not used to generate or
interpret results, design the analysis, or write the manuscript text.

## Citation

Please cite both the software (via this repository) and the associated
paper. See `CITATION.cff` for machine-readable citation metadata; GitHub
will render a "Cite this repository" button in the sidebar.

## Contact

Charles-Alexandre Roy — <c2roy@health.ucsd.edu>
Center for Computational Biology and Bioinformatics, UC San Diego.
