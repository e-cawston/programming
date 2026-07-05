---
title: "Modular RNA-seq pipeline for mst7 transcriptome"
output:
  html_document:
    toc: true
    toc_depth: 2
    theme: flatly
---




``` r
required_pkgs <- c(
  "dplyr", "readr", "purrr", "tidyr", "stringr", "tibble",
  "here", "ggplot2", "reshape2", "sleuth", "UpSetR", "RankProd"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}
```

```
## Error:
## ! Missing required packages: ggplot2, sleuth, UpSetR, RankProd
```

# Overview

This document runs the modular RNA-seq workflow for the mst7 transcriptome project. It reads the sample metadata, sources the reusable analysis scripts, runs sleuth-based DGE comparisons from the comparison table, generates DEG summary plots, creates upset plots, and generates correlation matrices.


``` r
library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(stringr)
library(tibble)
library(here)
library(ggplot2)
```

```
## Error in `library()`:
## ! there is no package called 'ggplot2'
```

``` r
library(reshape2)
```


``` r
project_dir <- here::here()
analysis_dir <- file.path(project_dir, "analysis")
scripts_dir <- file.path(analysis_dir, "scripts")
results_dir <- file.path(project_dir, "results")
fig_dir <- file.path(results_dir, "figures")
tables_dir <- file.path(results_dir, "tables")
dge_dir <- file.path(results_dir, "dge")
upset_dir <- file.path(fig_dir, "upset")
cor_dir <- file.path(fig_dir, "correlation")

for (dir_path in c(fig_dir, tables_dir, dge_dir, upset_dir, cor_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}
```


``` r
source(file.path(scripts_dir, "comparison_table.R"))
source(file.path(scripts_dir, "upset_functions.R"))
```

```
## Error:
## ! Install required packages before sourcing upset_functions.R: UpSetR
```

``` r
source(file.path(scripts_dir, "correlation_functions.R"))
```

```
## Error:
## ! Install required packages before sourcing correlation_functions.R: ggplot2
```

``` r
source(file.path(scripts_dir, "dge_functions.R"))
```

```
## Error:
## ! Install required packages before sourcing dge_functions.R: ggplot2, sleuth
```

``` r
source(file.path(scripts_dir, "rank_products_functions.R"))
```

```
## Error:
## ! Install required packages before sourcing rank_products_functions.R: ggplot2, RankProd
```


``` r
s2c <- read_csv(file.path(project_dir, "raw", "run_metadata.txt"), show_col_types = FALSE) %>%
  mutate(
    path = file.path(project_dir, "raw", path),
    name = str_remove(name, "_0h_\\d+$"),
    sample = paste0(name, "_", timepoint, "_", biorep),
    name = factor(name),
    timepoint = factor(timepoint),
    biorep = factor(biorep)
  ) %>%
  mutate(
    experiment = if_else(name %in% c("pmk1", "Guy11M"), "Exp2", "Exp1"),
    experiment = factor(experiment, levels = c("Exp1", "Exp2"))
  )
```

# Sleuth differential expression


``` r
dge_results <- run_comparisons(s2c, comparisons, out_dir = dge_dir)
```

```
## Error in `run_comparisons()`:
## ! could not find function "run_comparisons"
```


``` r
for (lab in unique(dge_results$comparison_label)) {
  if (!is.na(lab) && nzchar(lab)) {
    plot_dge_summary(dge_results, lab, out_dir = fig_dir)
  }
}
```

```
## Error:
## ! object 'dge_results' not found
```

# Upset plots


``` r
deg_sets <- load_deg_tables(
  comparison_labels = comparisons$label,
  dge_dir = dge_dir,
  qval_threshold = 0.05,
  fc_threshold = 0
)
```

```
## Error in `load_deg_tables()`:
## ! could not find function "load_deg_tables"
```

``` r
plot_upset(deg_sets, label = "all_comparisons", out_dir = upset_dir)
```

```
## Error in `plot_upset()`:
## ! could not find function "plot_upset"
```

# Correlation matrices


``` r
tpm_matrix <- load_tpm_matrix(file.path(project_dir, "raw", "all_samples_tpm_matrix.txt"))
```

```
## Error in `load_tpm_matrix()`:
## ! could not find function "load_tpm_matrix"
```

``` r
cor_matrix_all(tpm_matrix, s2c)
```

```
## Error in `cor_matrix_all()`:
## ! could not find function "cor_matrix_all"
```

``` r
cor_matrix_exp1(tpm_matrix, s2c)
```

```
## Error in `cor_matrix_exp1()`:
## ! could not find function "cor_matrix_exp1"
```

``` r
cor_matrix_exp2(tpm_matrix, s2c)
```

```
## Error in `cor_matrix_exp2()`:
## ! could not find function "cor_matrix_exp2"
```

``` r
cor_matrix_key_strains(tpm_matrix, s2c)
```

```
## Error in `cor_matrix_key_strains()`:
## ! could not find function "cor_matrix_key_strains"
```

# Rank-products plots


``` r
rp_results <- run_rank_products_comparisons(s2c, comparisons, out_dir = file.path(results_dir, "rank_products"))
```

```
## Error in `run_rank_products_comparisons()`:
## ! could not find function "run_rank_products_comparisons"
```

``` r
for (lab in unique(rp_results$comparison_label)) {
  if (!is.na(lab) && nzchar(lab)) {
    plot_rank_products_mirrored(
      rp_results,
      lab,
      out_dir = file.path(fig_dir, "rank_products")
    )
  }
}
```

```
## Error:
## ! object 'rp_results' not found
```

# Session info


``` r
sessionInfo()
```

```
## R version 4.5.3 (2026-03-11)
## Platform: x86_64-conda-linux-gnu
## Running under: Ubuntu 24.04.4 LTS
## 
## Matrix products: default
## BLAS/LAPACK: /home/Occ4m/miniforge3/envs/ds_env/lib/libopenblasp-r0.3.33.so;  LAPACK version 3.12.0
## 
## locale:
##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
## 
## time zone: Europe/London
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] reshape2_1.4.5 here_1.0.2     tibble_3.3.1   stringr_1.6.0  tidyr_1.3.2   
## [6] purrr_1.2.2    readr_2.2.0    dplyr_1.2.1   
## 
## loaded via a namespace (and not attached):
##  [1] crayon_1.5.3     vctrs_0.7.3      cli_3.6.6        knitr_1.51      
##  [5] rlang_1.2.0      xfun_0.57        stringi_1.8.7    generics_0.1.4  
##  [9] bit_4.6.0        glue_1.8.1       rprojroot_2.1.1  plyr_1.8.9      
## [13] hms_1.1.4        evaluate_1.0.5   tzdb_0.5.0       lifecycle_1.0.5 
## [17] compiler_4.5.3   Rcpp_1.1.1-1.1   pkgconfig_2.0.3  R6_2.6.1        
## [21] tidyselect_1.2.1 parallel_4.5.3   vroom_1.7.1      pillar_1.11.1   
## [25] magrittr_2.0.5   bit64_4.8.2      tools_4.5.3
```
