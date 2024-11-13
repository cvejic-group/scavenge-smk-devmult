#!/usr/bin/env Rscript

message("[INFO] Loading libraries...")

library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(magrittr)

################################################################################
# Snakemake Object
################################################################################

if (interactive()) { ## fake data for testing
  Snakemake <- setClass("Snakemake", slots=c(input="list", output="list",
                                             params="list", threads="numeric"))
  snakemake <- Snakemake(
    input=list(se="results/se/subsamples/se_sample300_rng20241104.Rds",
               mat_bg="results/bg/mat_bg200_sample300_rng20241104.Rds",
               bed_trait="/work/DevM_analysis/06.phenotypes/data/gwas/SCAVENGE-reproducibility/data/finemappedtraits_hg19/hg38/covid19hg_B1_FMdata.hg38.bed"),
    output=list(dev_trait="results/dev/scavenge_yu_2022/dev_covid19_bg200_sample300_rng20241104.Rds"),
    params=list(rseed="20241104"),
    threads=64L)
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

message("[INFO] Loading background peaks...")
mat_bg <- readRDS(snakemake@input$mat_bg)

message("[INFO] Loading trait BED...")
trait_import <- importBedScore(rowRanges(se_peak_cell), snakemake@input$bed_trait)

################################################################################
# Computing Deviations
################################################################################
message("[INFO] Computing GC bias...")
se_peak_cell %<>% addGCBias(genome=BSgenome.Hsapiens.UCSC.hg38)

message("[INFO] Computing Deviations...")
dev_trait <- computeWeightedDeviations(object=se_peak_cell,
                                       weights=trait_import,
                                       background_peaks=mat_bg)

################################################################################
# Export Deviations Object
################################################################################
message("[INFO] Exporting deviations...")
saveRDS(dev_trait, snakemake@output$dev_trait)

message("[INFO] Done!")
