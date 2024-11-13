#!/usr/bin/env Rscript

message("[INFO] Loading libraries...")

library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(SummarizedExperiment)
library(magrittr)

################################################################################
# Snakemake Object
################################################################################

if (interactive()) { ## fake data for testing
  Snakemake <- setClass("Snakemake", slots=c(input="list", output="list",
                                             params="list", threads="numeric"))
  snakemake <- Snakemake(
    input=list(se="results/se/subsamples/se_sample300_rng20241104.Rds"),
    output=list(mat="results/bg/mat_bg200_rng20241104.Rds"),
    params=list(n_background="200",
                rseed="20241104"),
    threads=4L)
}

## show snakemake values
print(snakemake)

################################################################################
# Set Parameters
################################################################################
message("[INFO] Initializing parameters...")
set.seed(as.integer(snakemake@params$rseed))

register(MulticoreParam(snakemake@threads))

################################################################################
# Load Data
################################################################################
message("[INFO] Loading SummarizedExperiment...")
se_peak_cell <- readRDS(snakemake@input$se)

################################################################################
# Generate Background Peaks
################################################################################
message("[INFO] Computing GC bias...")
se_peak_cell %<>% addGCBias(genome=BSgenome.Hsapiens.UCSC.hg38)

message("[INFO] Generating background peaks...")
mat_bg <- getBackgroundPeaks(se_peak_cell,
                             bias=rowData(se_peak_cell)$bias,
                             niterations=as.integer(snakemake@params$n_background))

################################################################################
# Export Matrix
################################################################################
message("[INFO] Exporting matrix...")
saveRDS(mat_bg, snakemake@output$mat)

message("[INFO] Done!")
