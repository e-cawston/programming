#!/usr/bin/env Rscript
# Install packages from CRAN
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install base packages first
base_pkgs <- c('ggplot2', 'UpSetR', 'openxlsx')
for (pkg in base_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(pkg, "already installed\n")
  }
}

# Try sleuth from bioconductor
if (!require('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager')
}

if (!requireNamespace('sleuth', quietly = TRUE)) {
  cat("Installing sleuth from Bioconductor\n")
  BiocManager::install('sleuth')
}

if (!requireNamespace('RankProd', quietly = TRUE)) {
  cat("Installing RankProd from Bioconductor\n")
  BiocManager::install('RankProd')
}

# Final check
cat("\nFinal package check:\n")
pkgs <- c('ggplot2', 'sleuth', 'UpSetR', 'RankProd', 'openxlsx')
for(p in pkgs) {
  cat(p, ":", requireNamespace(p, quietly=TRUE), "\n")
}
