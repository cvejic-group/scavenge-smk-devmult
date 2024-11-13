#!/usr/bin/env Rscript

message("[INFO] Loading libraries...")

library(SCAVENGE)
library(SummarizedExperiment)

################################################################################
# Snakemake Object
################################################################################

if (interactive()) { ## fake data for testing
  Snakemake <- setClass("Snakemake", slots=c(input="list", output="list",
                                             params="list", threads="numeric"))
  snakemake <- Snakemake(
    input=list(se="results/se/subsamples/se_sample300_rng20241104.Rds"),
    output=list(mat="results/knn/mat_knn30_sample300_rng20241104.Rds"),
    params=list(n_neighbors="30", rseed="20241104"),
    threads=4L)
}

## show snakemake values
print(snakemake)

################################################################################
# Set Parameters
################################################################################
set.seed(as.integer(snakemake@params$rseed))

################################################################################
# Load Data
################################################################################
message("[INFO] Loading SummarizedExperiment...")
se_peak_cell <- readRDS(snakemake@input$se)

message("[INFO] Computing mutual k-NN matrix...")
mat_mknn <- getmutualknn(lsimat=metadata(se_peak_cell)$reducedDims$LSI,
                         num_k=as.integer(snakemake@params$n_neighbors))

################################################################################
# Export Matrix
################################################################################
message("[INFO] Exporting matrix...")
saveRDS(mat_mknn, snakemake@output$mat)

message("[INFO] Done!")
