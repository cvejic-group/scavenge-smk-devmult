#!/usr/bin/env Rscript

message("[INFO] Loading libraries...")

library(chromVAR)
library(gchromVAR)
library(SummarizedExperiment)
library(DelayedArray)
library(HDF5Array)
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
    input=list(h5="data/se/se_peak10k_cell10k_48FL_assays.h5",
               rds="data/se/se_peak10k_cell10k_48FL_se.rds",
               mat_bg="results/bg/mat_bg200_rng20241104.Rds",
               bed_trait="/work/DevM_analysis/06.phenotypes/data/gwas/SCAVENGE-reproducibility/data/finemappedtraits_hg19/hg38/covid19hg_B1_FMdata.hg38.bed"),
    output=list(dev_trait="results/dev/scavenge_yu_2022/dev_covid19_bg200_rng20241104.Rds"),
    params=list(rseed="20241104"),
    threads=48L)
}

## show snakemake values
print(snakemake)

################################################################################
# Set Parameters
################################################################################
message("[INFO] Initializing parameters...")
set.seed(as.integer(snakemake@params$rseed))

setAutoBPPARAM(MulticoreParam(snakemake@threads))
register(MulticoreParam(snakemake@threads))

stopifnot(dirname(snakemake@input$h5) == dirname(snakemake@input$rds))
PATH_INPUT=dirname(snakemake@input$rds)
PREFIX_INPUT=gsub("se.rds", "", basename(snakemake@input$rds))

################################################################################
# Functions
################################################################################
## We override these so that the specialized `DelayedArray` implementation is
## used rather than trigger a coercion to a Matrix object.
#' @describeIn getFragmentsPerPeak method for SummarizedExperiment
#' @export
setMethod("getFragmentsPerPeak", c(object = "SummarizedExperiment"), 
          function(object) {
            rowSums(counts(object))
          })

#' @describeIn getFragmentsPerSample method for SummarizedExperiment
#' @export
setMethod("getFragmentsPerSample", c(object = "SummarizedExperiment"), 
          function(object) {
            colSums(counts(object))
          })

#' @describeIn getFragmentsPerPeak method for DelayedMatrix
#' @export
setMethod("getFragmentsPerPeak", c(object = "DelayedMatrix"), 
          function(object) {
            DelayedArray::rowSums(object)
          })

#' @describeIn getFragmentsPerSample method for DelayedMatrix
#' @export
setMethod("getFragmentsPerSample", c(object = "DelayedMatrix"), 
          function(object) {
            DelayedArray::colSums(object)
          })

#' @describeIn getTotalFragments method for DelayedMatrix
#' @export
setMethod("getTotalFragments", c(object = "DelayedMatrix"), 
          function(object) {
            sum(getFragmentsPerSample(object))
          })

#' @describeIn computeExpectations method for SummarizedExperiment with counts
#' slot
#' @export
setMethod("computeExpectations", c(object = "SummarizedExperiment"),
          function(object,
                   norm = FALSE,
                   group = NULL) {
            fpp <- getFragmentsPerPeak(object)
            fpp / sum(fpp)
          })

################################################################################
# Load Data
################################################################################
message("[INFO] Loading SummarizedExperiment...")
se_peak_cell <- loadHDF5SummarizedExperiment(PATH_INPUT, PREFIX_INPUT)

## chromVAR assumes assay is named "counts"
assayNames(se_peak_cell) <- c("counts")

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
