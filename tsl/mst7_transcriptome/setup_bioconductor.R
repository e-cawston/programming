#!/usr/bin/env Rscript
# Install BiocManager packages for new rna_analysis environment

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)

# Install BiocManager first
if (!require('BiocManager', quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages('BiocManager', dependencies = TRUE, quiet = TRUE)
}

# Install remaining dependencies
cat("Installing rhdf5...\n")
BiocManager::install('rhdf5', update = FALSE, ask = FALSE, quiet = TRUE)

cat("Installing sleuth...\n")
BiocManager::install('sleuth', update = FALSE, ask = FALSE, quiet = TRUE)

cat("Installing UpSetR from CRAN...\n")
install.packages('UpSetR', dependencies = TRUE, quiet = TRUE)

cat("Installing RankProd...\n")
BiocManager::install('RankProd', update = FALSE, ask = FALSE, quiet = TRUE)

# Final check
cat("\n=== Final Package Check ===\n")
pkgs <- c('ggplot2', 'sleuth', 'UpSetR', 'RankProd', 'dplyr', 'readr', 'ggplot2')
for(p in pkgs) {
  status <- requireNamespace(p, quietly = TRUE)
  cat(p, ":", if(status) "✓" else "✗", "\n")
}
