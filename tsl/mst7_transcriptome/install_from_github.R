#!/usr/bin/env Rscript
# Install sleuth and RankProd from GitHub/alternative sources

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)

cat("=== Installing from GitHub and alternative sources ===\n\n")

# Install devtools for remotes::install_github
if (!require('devtools', quietly = TRUE)) {
  cat("Installing devtools...\n")
  install.packages('devtools', dependencies = TRUE)
}

library(devtools)

# Try sleuth from GitHub (Caltech/sleuth original repo)
cat("\n1. Installing sleuth from GitHub...\n")
tryCatch(
  {
    devtools::install_github("pachterlab/sleuth", dependencies = TRUE, upgrade = FALSE)
    cat("   ✓ sleuth installed from GitHub\n")
  },
  error = function(e) {
    cat("   ✗ GitHub install failed:", e$message, "\n")
    cat("   Trying CRAN mirror...\n")
    tryCatch(
      install.packages('sleuth'),
      error = function(e2) {
        cat("   ✗ CRAN install also failed:", e2$message, "\n")
      }
    )
  }
)

# For RankProd, try GitHub or source
cat("\n2. Installing RankProd...\n")
tryCatch(
  {
    devtools::install_github("fredcommo/RankProd", dependencies = TRUE, upgrade = FALSE)
    cat("   ✓ RankProd installed from GitHub\n")
  },
  error = function(e) {
    cat("   ✗ GitHub install failed:", e$message, "\n")
    cat("   RankProd may need manual installation or system dependencies\n")
  }
)

# Check what we have
cat("\n\n=== Final Status ===\n")
check_pkgs <- c('sleuth', 'RankProd')
for(p in check_pkgs) {
  status <- requireNamespace(p, quietly = TRUE)
  symbol <- if(status) "✓" else "✗"
  cat(sprintf("%s %-15s\n", symbol, p))
}

cat("\nIf packages are still not installed, system dependencies may be needed.\n")
