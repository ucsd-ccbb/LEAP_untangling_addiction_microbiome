# ConQuR batch correction for OneCodex shotgun metagenomic counts
# Covariates: opioid_status, age, sex
# Reference batch: "0"

library(ConQuR)
library(doParallel)
options(warn = -1)

# --- Inputs / outputs ---
input_file  <- "sample_taxa_counts.tsv"
output_file <- "sample_taxa_counts_corrected.tsv"
pcoa_file   <- "PCoA.pdf"

# --- Read data ---
# Input TSV has: Sample.ID column, taxonomic count columns, and 4 trailing
# metadata columns (batch, opioid_status, age, sex).
Sample_Data_real <- read.table(file = input_file, sep = "\t", header = TRUE)

# Taxa table: drop Sample.ID (col 1) and the trailing 4 metadata columns
taxa <- Sample_Data_real[, -c(1, (ncol(Sample_Data_real) - 3):ncol(Sample_Data_real))]

# Batch and covariates
batchid <- factor(Sample_Data_real[["batch"]])
covar <- data.frame(
  opioid_status = factor(Sample_Data_real[["opioid_status"]]),
  age           = Sample_Data_real[["age"]],
  sex           = factor(Sample_Data_real[["sex"]])
)

# --- Metadata balance checks ---
print(table(batchid, covar$opioid_status))
print(table(batchid, covar$sex))

# --- Filter taxa with all-zero counts ---
zero_rate    <- apply(taxa, 2, function(z) mean(z > 0))
taxa_reduced <- taxa[, which(zero_rate > 0)]

# --- ConQuR correction ---
taxa_corrected <- ConQuR(
  tax_tab    = taxa,
  batchid    = batchid,
  covariates = covar,
  batch_ref  = "0"
)

# --- Diagnostic PCoA ---
pdf(pcoa_file, width = 12, height = 6)
par(mfrow = c(2, 4))
Plot_PCoA(TAX = taxa,           factor = batchid, main = "Before Correction, Bray-Curtis")
Plot_PCoA(TAX = taxa_corrected, factor = batchid, main = "ConQuR (Default), Bray-Curtis")
dev.off()

# --- PERMANOVA R^2 (variance explained by batch and opioid status) ---
print(PERMANOVA_R2(TAX = taxa,           batchid = batchid, covariates = covar, key_index = 1))
print(PERMANOVA_R2(TAX = taxa_corrected, batchid = batchid, covariates = covar, key_index = 1))

# --- Write corrected table ---
rownames(taxa_corrected) <- Sample_Data_real[["Sample.ID"]]
write.table(
  taxa_corrected,
  file      = output_file,
  sep       = "\t",
  row.names = TRUE,
  quote     = FALSE
)
