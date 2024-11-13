#!/usr/bin/env Rscript

message("[INFO] Loading libraries...")

library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(chromVAR)
library(SummarizedExperiment)
library(DelayedArray)
library(HDF5Array)
library(magrittr)

################################################################################
# Snakemake Object
################################################################################

if (interactive()) { ## fake data for testing
  Snakemake <- setClass("Snakemake", slots=c(input="list", output="list",
                                             params="list", threads="numeric"))
  snakemake <- Snakemake(
    input=list(h5="data/se/se_peak10k_cell10k_48FL_assays.h5",
               rds="data/se/se_peak10k_cell10k_48FL_se.rds"),
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

setAutoBPPARAM(MulticoreParam(snakemake@threads))
register(MulticoreParam(snakemake@threads))

stopifnot(dirname(snakemake@input$h5) == dirname(snakemake@input$rds))
PATH_INPUT=dirname(snakemake@input$rds)
PREFIX_INPUT=gsub("se.rds", "", basename(snakemake@input$rds))

###############
## FUNCTIONS ##
###############
## We override these so that the specialized `DelayedArray` implementation is used
## rather than trigger a coercion to a Matrix object.
setMethod("getFragmentsPerPeak", c(object = "SummarizedExperiment"),
          function(object) {
            rowSums(counts(object))
          })

setMethod("getFragmentsPerPeak", c(object = "DelayedMatrix"),
          function(object) {
            DelayedArray::rowSums(object)
          })

################################################################################
# Load Data
################################################################################
message("[INFO] Loading SummarizedExperiment...")
se_peak_cell <- loadHDF5SummarizedExperiment(PATH_INPUT, PREFIX_INPUT)

## chromVAR assumes assay is named "counts"
assayNames(se_peak_cell) <- c("counts")

################################################################################
# Generate Background Peaks
################################################################################
message("[INFO] Computing GC bias...")
se_peak_cell %<>% addGCBias(genome=BSgenome.Hsapiens.UCSC.hg38)


message("[INFO] Generating background peaks...")
## note that we manually call into the raw method to avoid the checks that would
## otherwise coerce the `DelayedArray` to a `Matrix`.
mat_bg <- chromVAR:::get_background_peaks_core(assay(se_peak_cell, 'counts'),
                                                     bias=rowData(se_peak_cell)$bias,
                                                     niterations=as.integer(snakemake@params$n_background))

################################################################################
# Export Matrix
################################################################################
message("[INFO] Exporting matrix...")
saveRDS(mat_bg, snakemake@output$mat)

message("[INFO] Done!")
