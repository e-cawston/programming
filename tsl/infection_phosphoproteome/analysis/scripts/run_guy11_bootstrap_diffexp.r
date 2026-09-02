## Driver script: Guy11 within-experiment time-course differential
## abundance (10,000-iteration bootstrap-t vs 0h, per experiment), then
## between-experiment kinetics correlation for the significant peptides.
## This is a multi-hour computation -- run via `Rscript` in the background,
## not inside an interactive knit. See
## analysis/scripts/phosphoproteome_diffexp_bootstrap.r for the full
## rationale (why pepdiff::compare() isn't used directly, why Modified
## Sequence/VSN-imputed data are used, the MNAR-imputation caveat, etc).

suppressMessages({
  library(pepdiff)
  library(tidyverse)
})
source("analysis/scripts/phosphoproteome_diffexp_bootstrap.r")

dir.create("results", showWarnings = FALSE, recursive = TRUE)
checkpoint_dir <- "results/guy11_bootstrap_checkpoints"

# Respect a SLURM CPU allocation if running under sbatch/srun (SLURM_CPUS_PER_TASK);
# parallel::detectCores() reports the physical node's core count, not what the
# job was actually granted, so it can over-subscribe on a shared HPC node.
n_cores <- {
  slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)
  if (!is.na(slurm_cpus) && nzchar(slurm_cpus)) {
    as.integer(slurm_cpus)
  } else {
    max(1, parallel::detectCores() - 1)
  }
}
message("Using n_cores = ", n_cores)

message("[", format(Sys.time(), "%H:%M:%S"), "] Building long-format Guy11 comparison table (VSN-normalized, imputed)...")
guy11_long_data <- build_guy11_between_experiments_data(
  path_a = "results/mst7_05_vsn_mar_mnar_imputed.xlsx",
  path_b = "results/pmk1_05_vsn_mar_mnar_imputed.xlsx",
  experiment_a = "mst7",
  experiment_b = "pmk1"
)
writexl::write_xlsx(guy11_long_data, "results/guy11_between_experiments_long_imputed.xlsx")

# ---- Stage 1: within-experiment 0h-vs-timepoint bootstrap-t, per experiment ----
message("[", format(Sys.time(), "%H:%M:%S"), "] Starting 10,000-iteration bootstrap-t, mst7 experiment (Guy11, vs 0h)...")
mst7_results <- run_bootstrap_vs_ref(
  guy11_long_data %>% dplyr::filter(experiment == "mst7"),
  compare = "timepoint", ref = "0h", n_boot = 10000, seed = 1,
  checkpoint_dir = checkpoint_dir, label = "mst7"
)

message("[", format(Sys.time(), "%H:%M:%S"), "] Starting 10,000-iteration bootstrap-t, pmk1 experiment (Guy11, vs 0h)...")
pmk1_results <- run_bootstrap_vs_ref(
  guy11_long_data %>% dplyr::filter(experiment == "pmk1"),
  compare = "timepoint", ref = "0h", n_boot = 10000, seed = 2,
  checkpoint_dir = checkpoint_dir, label = "pmk1"
)

stage1_results <- dplyr::bind_rows(
  mst7_results %>% dplyr::mutate(experiment = "mst7", .before = 1),
  pmk1_results %>% dplyr::mutate(experiment = "pmk1", .before = 1)
)
writexl::write_xlsx(stage1_results, "results/guy11_timepoint_vs_0h_bootstrap.xlsx")

ggplot2::ggsave(
  "results/guy11_timepoint_vs_0h_bootstrap_volcano.png",
  plot_diffexp_volcano(stage1_results, n_boot = 10000),
  width = 16, height = 7
)

# ---- Filter: |log2FC| >= 1 & p < 0.05 at any comparison, in either experiment ----
sig_peptides <- union(
  filter_significant_peptides(mst7_results, fc_threshold = 1, p_threshold = 0.05),
  filter_significant_peptides(pmk1_results, fc_threshold = 1, p_threshold = 0.05)
)
message(sprintf(
  "[%s] Significant peptides (|log2FC|>=1 & p<0.05, either experiment): %d of %d shared",
  format(Sys.time(), "%H:%M:%S"), length(sig_peptides), dplyr::n_distinct(guy11_long_data$peptide)
))
writeLines(sig_peptides, "results/guy11_significant_peptides.txt")

# ---- Stage 2: between-experiment kinetics correlation for significant peptides ----
message("[", format(Sys.time(), "%H:%M:%S"), "] Correlating peptide kinetics between experiments...")
kinetics_cor <- correlate_peptide_kinetics(
  guy11_long_data, sig_peptides,
  experiment_a = "mst7", experiment_b = "pmk1"
)
writexl::write_xlsx(kinetics_cor, "results/guy11_kinetics_correlation.xlsx")

cat("\nKinetics correlation summary (Pearson r per peptide, mst7 vs pmk1):\n")
print(summary(kinetics_cor$r))

ggplot2::ggsave(
  "results/guy11_kinetics_correlation_distribution.png",
  plot_kinetics_correlation_distribution(kinetics_cor),
  width = 10, height = 7
)

message("[", format(Sys.time(), "%H:%M:%S"), "] DONE. Key outputs:")
message("  results/guy11_timepoint_vs_0h_bootstrap.xlsx        -- stage 1 per-peptide log2FC/p vs 0h, both experiments")
message("  results/guy11_significant_peptides.txt              -- peptides passing |log2FC|>=1 & p<0.05")
message("  results/guy11_kinetics_correlation.xlsx             -- stage 2 per-peptide Pearson r (mst7 vs pmk1 trajectory)")
message("  results/guy11_kinetics_correlation_distribution.png -- histogram of those r values")
