#!/usr/bin/env Rscript
# Simple direct installation of binary packages

options(repos = c(CRAN = "https://cran.r-project.org"))
options(timeout = 600)
options(pkgType = "binary")  # Force binary packages

cat("=== Installing pre-compiled binary packages ===\n\n")

# Try installing as binaries first
packages_to_try <- list(
  list(name = "sleuth", repos = "https://bioconductor.org/packages/3.22/bioc"),
  list(name = "RankProd", repos = "https://bioconductor.org/packages/3.22/bioc"),
  list(name = "UpSetR", repos = "https://cran.r-project.org"),
  list(name = "rhdf5", repos = "https://bioconductor.org/packages/3.22/bioc")
)

cat("Attempting binary package installation...\n\n")

for (pkg_info in packages_to_try) {
  pkg_name <- pkg_info$name
  pkg_repos <- pkg_info$repos
  
  cat("Trying:", pkg_name, "\n")
  tryCatch(
    {
      install.packages(
        pkg_name,
        repos = pkg_repos,
        type = "binary",
        dependencies = FALSE
      )
      if (requireNamespace(pkg_name, quietly = TRUE)) {
        cat("  ✓ Successfully installed\n\n")
      } else {
        cat("  ✗ Package installed but cannot load\n\n")
      }
    },
    error = function(e) {
      cat("  ✗ Failed:", e$message, "\n\n")
    },
    warning = function(w) {
      cat("  ⚠ Warning:", w$message, "\n\n")
    }
  )
}

cat("=== Installation Summary ===\n")
check_pkgs <- c('sleuth', 'RankProd', 'UpSetR', 'rhdf5')
installed_count <- 0
for(p in check_pkgs) {
  status <- requireNamespace(p, quietly = TRUE)
  if (status) {
    cat("✓", p, "\n")
    installed_count <- installed_count + 1
  } else {
    cat("✗", p, "\n")
  }
}

cat("\nInstalled:", installed_count, "of", length(check_pkgs), "packages\n")
