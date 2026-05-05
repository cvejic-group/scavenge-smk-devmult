[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20024542.svg)](https://doi.org/10.5281/zenodo.20024542)

# scavenge-smk
This is a Snakemake pipeline for running SCAVENGE on large scATAC-seq datasets.
Specifically, the pipeline is written for an `SummarizedExperiment` object
that represents peak-cell counts.

## Overview

The pipeline uses the arguments in [`config/config.yaml`](config/config.yaml) to compute chromVAR deviations and 
SCAVENGE TRS scores using peak-cell counts input as an `HDF5Array`-backed `SummarizedExperiment` object.
We recommend SNP sets (`trait_sets`) be defined modularly in separate configuration YAMLs (e.g.,
[`config/config-gwas.yaml`](config/config-gwas.yaml).
The pipeline includes both a **subsampling** mode and a **full** mode. The Snakefile `MAIN` section 
defines what outputs will be generated.

### Parameter Summary
The following SCAVENGE parameters were used at runtime:

```yaml
n_background: 200
n_neighbors: 30
n_perms: 1000
frac_seed_max: 0.05
frac_scale_cut: 0.01
gamma_rw: 0.05
q_ceiling: 0.95
pvalue_sig_cell: 0.05
```

More details available in [`config/config.yaml`](config/config.yaml).

### Inputs
- `se_peak_cell` - a peak-cell count `SummarizedExperiment` object
- BED files for fine-mapped traits - hg38 liftOver outputs from https://github.com/drewmard/t21_multiome

### Outputs
- `results/trs/` - CSV tables per trait and per trait set with chromVAR and SCAVENGE scores
- `results/reports/` - summarizations per trait and per trait set

## Running the Pipeline
The pipeline has been previously run with Snakemake v8.15 in both local and SLURM-based modes.
It is expected to run with Conda+Mamba (Miniforge installation).

### Install Snakemake from YAML
YAML representations of the Conda environment with Snakemake v8.15 are provided in minimal (".min") 
and full (".full") forms for reuse and reproduction, respectively. The minimal environment could be
created with

```bash
conda env create -n smk_8_15 -f workflow/envs/smk_8_15.min.yaml
```

### Local Execution
The pipeline can be executed locally with

```bash
snakemake --use-conda --cores 64 -s workflow/Snakefile --configfile config/config-gwas.yaml
```

Adjust the `--cores` argument accordingly.

### SLURM Execution
To run on a SLURM configuration, first configure a SLURM profile for Snakemake.
We have used

**~/.config/snakemake/slurm_basic/config.v8+.yaml**
```yaml
executor: slurm
use-conda: true
jobs: 10000
default-resources:
  mem_mb_per_cpu: 4096
```

which can be used with:

```bash
snakemake --profile slurm_basic -s workflow/Snakefile --configfile config-gwas.yaml
```

### Generating Individual Outputs
As a Snakemake pipeline, one can generate individual files *ad hoc* by specifying them in the command.
For example, one could request a specific report on trait, like `covid19` with

```bash
snakemake --profile slurm_basic \
  -s workflow/Snakefile --configfile config-gwas.yaml \
  results/reports/trs/scavenge_yu_2022/report_trs_covid19_knn30_bg200_sample300_rng20241104.html
```

Note that the parameters in file name must match those in the `config.yaml` or be overridden, e.g.,

```bash
snakemake --profile slurm_basic \
  -s workflow/Snakefile --configfile config-gwas.yaml \
  --config seed_rng=20250101 \
  results/reports/trs/scavenge_yu_2022/report_trs_covid19_knn30_bg200_sample300_rng20250101.html
```
