#!/usr/bin/env Rscript

message("[INFO] Loading libraries...")

library(SCAVENGE)
library(SummarizedExperiment)
library(BiocParallel)
library(magrittr)
library(tibble)
library(dplyr)
library(readr)

################################################################################
# Snakemake Object
################################################################################

if (interactive()) { ## fake data for testing
  Snakemake <- setClass("Snakemake", slots=c(input="list", output="list",
                                             params="list", threads="numeric"))
  snakemake <- Snakemake(
    input=list(mat_knn="results/knn/mat_knn30_sample300_rng20241104.Rds",
               dev_trait="results/dev/scavenge_yu_2022/dev_covid19_bg200_sample300_rng20241104.Rds"),
    output=list(csv_trait="results/trs/scavenge_yu_2022/df_trs_covid19_knn30_bg200_sample300_rng20241104.csv.gz"),
    params=list(frac_seed_max="0.05",
                frac_scale_cut="0.01",
                gamma_rw="0.05",
                n_perms="1000",
                pvalue_sig_cell="0.05",
                q_ceiling="0.95",
                rseed="20241104",
                trait_set_id="scavenge_yu_2022",
                trait_id="covid19"),
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
message("[INFO] Loading trait deviations...")
dev_trait <- readRDS(snakemake@input$dev_trait)

message("[INFO] Loading kNN matrix...")
mat_knn <- readRDS(snakemake@input$mat_knn)

################################################################################
# Table Construction
################################################################################
message("[INFO] Reformating to tibble...")
df_dev <- assay(dev_trait, "deviations") %>% 
  t %>% `colnames<-`("deviation") %>% as_tibble(rownames="cell_id")
df_z <- assay(dev_trait, "z") %>% 
  t %>% `colnames<-`("z_score") %>% as_tibble(rownames="cell_id")
df_trs <- full_join(df_dev, df_z, by="cell_id")

head(df_trs)

################################################################################
# SCAVENGE
################################################################################
message("[INFO] Identifying SCAVENGE seeds...")
df_trs %<>%
  mutate(is_seed=seedindex(z_score=z_score, 
                           percent_cut=as.numeric(snakemake@params$frac_seed_max)))

message("[INFO] Calculating scale factor...")
scale_factor <- cal_scalefactor(z_score=df_trs$z_score,
                                percent_cut=as.numeric(snakemake@params$frac_scale_cut))

message("[INFO] Running network propagation...")
df_trs %<>%
  mutate(np_score=randomWalk_sparse(mat_knn, 
                                    queryCells=cell_id[is_seed],
                                    gamma=as.numeric(snakemake@params$gamma_rw))) %>%
  mutate(is_omitted=np_score == 0)

message(sprintf("[INFO] Omitting %d cells as unreachable.", sum(df_trs$is_omitted)))
df_trs_omitted <- filter(df_trs, is_omitted)

message("[INFO] Computing TRS...")
df_trs %<>% filter(!is_omitted) %>%
  mutate(trs_raw=capOutlierQuantile(x=np_score, q_ceiling=as.numeric(snakemake@params$q_ceiling))) %>%
  mutate(trs_scaled=max_min_scale(trs_raw) * scale_factor)

message("[INFO] Testing trait enrichment...")
df_sig <- get_sigcell_simple(knn_sparse_mat=mat_knn[df_trs$cell_id, df_trs$cell_id],
                             seed_idx=df_trs$is_seed,
                             topseed_npscore=df_trs$np_score,
                             permutation_times=as.integer(snakemake@params$n_perms),
                             true_cell_significance=as.numeric(snakemake@params$pvalue_sig_cell), 
                             rda_output=FALSE, 
                             mycores=snakemake@threads,
                             rw_gamma=as.numeric(snakemake@params$gamma_rw))

df_trs_final <- df_trs %>%
  mutate(is_enriched=df_sig$true_cell_top_idx) %>%
  bind_rows(df_trs_omitted) %>%
  mutate(trait_set_id=snakemake@params$trait_set_id,
         trait_id=snakemake@params$trait_id) %>%
  dplyr::select(cell_id, trait_set_id, trait_id, deviation, z_score,
                np_score, trs_raw, trs_scaled, is_seed, is_enriched)

################################################################################
# Export Table
################################################################################
message("[INFO] Exporting TRS table...")
write_csv(df_trs_final, snakemake@output$csv_trait)

message("[INFO] Done!")
