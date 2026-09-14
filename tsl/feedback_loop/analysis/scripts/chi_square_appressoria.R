#!/usr/bin/env Rscript

script_path <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_path[grepl("^--file=", script_path)])
script_dir <- dirname(normalizePath(script_path))
source(file.path(script_dir, "appressorium_assay.R"))

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else
  file.path(script_dir, "..", "..", "raw", "appressoria_phenotype.csv")
output_file <- if (length(args) >= 2) args[[2]] else
  file.path(script_dir, "..", "..", "figures", "appressoria_chi_square_results.csv")

results <- chi_square_appressoria(input_file)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.csv(results, output_file, row.names = FALSE)
print(results)
