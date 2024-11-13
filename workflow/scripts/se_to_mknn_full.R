#!/usr/bin/env Rscript

message("[INFO] Loading libraries...")

library(SCAVENGE)
library(SummarizedExperiment)
library(HDF5Array)

################################################################################
# Snakemake Object
################################################################################

if (interactive()) { ## fake data for testing
  Snakemake <- setClass("Snakemake", slots=c(input="list", output="list",
                                             params="list", threads="numeric"))
  snakemake <- Snakemake(
    input=list(h5="data/se/se_peak10k_cell10k_48FL_assays.h5",
               rds="data/se/se_peak10k_cell10k_48FL_se.rds"),
    output=list(mat="results/knn/mat_knn30_rng20241104.Rds"),
    params=list(n_neighbors="30", rseed="20241104"),
    threads=4L)
}

## show snakemake values
print(snakemake)

################################################################################
# Set Parameters
################################################################################
set.seed(as.integer(snakemake@params$rseed))

stopifnot(dirname(snakemake@input$h5) == dirname(snakemake@input$rds))
PATH_INPUT=dirname(snakemake@input$rds)
PREFIX_INPUT=gsub("se.rds", "", basename(snakemake@input$rds))

################################################################################
# Load Data
################################################################################
message("[INFO] Loading SummarizedExperiment...")
se_peak_cell <- loadHDF5SummarizedExperiment(PATH_INPUT, PREFIX_INPUT)

message("[INFO] Computing mutual k-NN matrix...")
mat_mknn <- getmutualknn(lsimat=metadata(se_peak_cell)$reducedDims$LSI,
                         num_k=as.integer(snakemake@params$n_neighbors))

################################################################################
# Export Matrix
################################################################################
message("[INFO] Exporting matrix...")
saveRDS(mat_mknn, snakemake@output$mat)

message("[INFO] Done!")
