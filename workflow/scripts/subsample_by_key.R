#!/usr/bin/env Rscript

library(magrittr)
library(SummarizedExperiment)
library(DelayedArray)
library(HDF5Array)
library(dplyr)
library(tibble)
library(BiocParallel)

################################################################################
# Snakemake Object
################################################################################

if (interactive()) { ## fake data for testing
  Snakemake <- setClass("Snakemake", slots=c(input="list", output="list",
                                             params="list", threads="numeric"))
  snakemake <- Snakemake(
    input=list(h5="data/se/se_peak10k_cell10k_48FL_assays.h5",
               rds="data/se/se_peak10k_cell10k_48FL_se.rds"),
    output=list(rds="results/se/subsamples/se_sample_300_r20241104.Rds"),
    params=list(key_subsample="sampleID",
                n_subsample="300",
                rseed="20241104"),
    threads=48L)
}

## show snakemake values
print(snakemake)

################################################################################
# Set Parameters
################################################################################
set.seed(as.integer(snakemake@params$rseed))

setAutoBPPARAM(MulticoreParam(snakemake@threads))
register(MulticoreParam(snakemake@threads))

stopifnot(dirname(snakemake@input$h5) == dirname(snakemake@input$rds))
PATH_INPUT=dirname(snakemake@input$rds)
PREFIX_INPUT=gsub("se.rds", "", basename(snakemake@input$rds))

################################################################################
# Load Data
################################################################################
se_peak_cell <- loadHDF5SummarizedExperiment(PATH_INPUT, PREFIX_INPUT)

## chromVAR will assume assay name is `counts`
assayNames(se_peak_cell) <- c("counts")

################################################################################
# Generate Subsampled Data
################################################################################

idx_cell_subsample <- colData(se_peak_cell) %>%
  as_tibble(rownames="cell_id") %>%
  group_by(across(all_of(snakemake@params$key_subsample))) %>%
  slice_sample(n=as.integer(snakemake@params$n_subsample)) %>%
  ungroup() %>%
  pull(cell_id)

se_peak_sample <- se_peak_cell[,idx_cell_subsample]

## keep in memory
assay(se_peak_sample, "counts") %<>% as("CsparseMatrix")

## subset reducedDims
metadata(se_peak_sample)$reducedDims$UMAP %<>% `[`(colnames(se_peak_sample),)
metadata(se_peak_sample)$reducedDims$LSI %<>% `[`(colnames(se_peak_sample),)

idx_peaks_nz <- rowSums(assay(se_peak_sample, "counts")) > 0
print(sprintf("Subsample has no fragments in %d peaks", sum(!idx_peaks_nz)))

## remove empty rows
se_peak_sample %<>% `[`(idx_peaks_nz,)

################################################################################
# Export SE
################################################################################
saveRDS(se_peak_sample, snakemake@output$rds)
