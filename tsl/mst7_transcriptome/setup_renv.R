#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("Installing renv...\n")
install.packages('renv', dependencies = TRUE)

cat("\nInitializing renv project...\n")
renv::init(project = '/home/Occ4m/programming/tsl/mst7_transcriptome')

cat("\nInstalling required packages...\n")
renv::install(c(
  'dplyr', 'readr', 'purrr', 'tidyr', 'stringr', 'tibble',
  'here', 'ggplot2', 'reshape2', 'UpSetR',
  'sleuth', 'RankProd'
))

cat("\nFinal check:\n")
renv::status()
