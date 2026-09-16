#!/usr/bin/env Rscript
# Install packages using pak - modern R package manager
# pak handles binary packages better than BiocManager with conda

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)

cat("=== Installing packages with pak ===\n\n")

# Install pak first if needed
if (!require('pak', quietly = TRUE)) {
  cat("Installing pak package manager...\n")
  install.packages('pak', dependencies = TRUE)
}

library(pak)

# Configure pak for better performance
cat("Configuring pak...\n")
pak::pak_config_set(list(
  "use_cran_binary" = TRUE,
  "cran_mirror" = "https://cloud.r-project.org"
))

cat("\nInstalling packages from CRAN and Bioconductor...\n\n")

# List of packages to install
packages <- c(
  # Already have these, but ensure they're up to date
  "ggplot2",
  "openxlsx",
  "dplyr", 
  "readr",
  "tibble",
  "tidyr",
  "stringr",
  "purrr",
  "here",
  "reshape2",
  "UpSetR",
  # The problematic ones - try via pak
  "sleuth",
  "RankProd"
)

cat("Installing ", length(packages), " packages...\n\n")

# Use pak to install - it handles dependencies better
pak::pkg_install(packages, ask = FALSE, upgrade = FALSE)

cat("\n\n=== Installation Complete ===\n")
cat("\nFinal package check:\n")

check_pkgs <- c('ggplot2', 'sleuth', 'UpSetR', 'RankProd', 'openxlsx', 'dplyr', 'readr', 'tidyr')
for(p in check_pkgs) {
  status <- requireNamespace(p, quietly = TRUE)
  symbol <- if(status) "✓" else "✗"
  cat(sprintf("%s %-15s %s\n", symbol, p, if(status) "" else "NOT INSTALLED"))
}
