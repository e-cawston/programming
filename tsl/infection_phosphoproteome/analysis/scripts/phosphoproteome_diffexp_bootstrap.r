## Guy11 within-experiment time-course differential abundance, then
## between-experiment kinetics correlation, using pepdiff's bootstrap-t test.
##
## Background
## ----------
## Guy11 (wild type) was run as the shared control genotype in both the
## mst7 experiment and the pmk1 experiment -- two separate MS submissions/
## batches. This analysis:
##
##   1. Within EACH experiment separately, tests every timepoint
##      (1h, 1.5h, 2h, 4h, 6h) against the 0h/spores reference for Guy11,
##      and keeps peptides that change at least 2-fold (|log2FC| >= 1) with
##      p < 0.05 at any of those comparisons, in either experiment.
##   2. For that significant peptide set, correlates each peptide's actual
##      mean-intensity trajectory across all 6 timepoints between the two
##      experiments (Pearson r per peptide) -- i.e. do the two experiments
##      agree not just on the size of a peptide's change, but on the SHAPE
##      of its response over time (its kinetics)?
##
## Step 2 deliberately correlates real VSN intensities (not log2FC-vs-0h
## values): fold-change-vs-0h trivially equals 0 at 0h in both experiments,
## which would artificially inflate every peptide's correlation toward
## agreement. Real per-timepoint intensities have no such shared anchor.
##
## Data used
## ---------
## VSN-normalized, IMPUTED intensities (results/<exp>_05_vsn_mar_mnar_
## imputed.xlsx), per explicit instruction. This gives every peptide a
## complete 3-replicate x 6-timepoint block in both experiments, which
## step 2 needs for a full 6-point kinetics trajectory per peptide.
##
## Caveat worth knowing: MNAR-imputed cells (a whole genotype x timepoint
## group with all 3 replicates missing) are filled with that SAMPLE's own
## low (1st-percentile) quantile value, not a real biological replicate --
## so a peptide/timepoint group that was entirely MNAR-imputed will show
## artificially low within-group variance in the bootstrap-t test, which
## can inflate its apparent significance. This isn't corrected for here;
## worth cross-checking any striking hit against the non-imputed value
## (results/<exp>_05_vsn_normalized_nonimputed.xlsx) before trusting it.
##
## Peptide identity
## -----------------
## "Peptide Sequence" is NOT unique per row (~2-6 rows can share the same
## backbone sequence at different phosphosites/modification states).
## "Modified Sequence" IS unique per row in both experiments' filtered
## data, so it is used as the peptide/site identifier fed to pepdiff.
## Analysis is restricted to phosphosites observed in both experiments,
## since step 2 needs a value in both to correlate.
##
## Why this bypasses pepdiff::compare()
## -------------------------------------
## pepdiff::compare(method = "pairwise", test = "bootstrap_t") always calls
## test_bootstrap_t() with its default n_boot = 1000 and does not forward
## '...' through to the test function (see pepdiff's compare_pairwise()),
## so there is no way to request 10,000 iterations via the high-level
## compare() API. To honour the requested 10,000-iteration bootstrap, this
## script calls pepdiff::test_bootstrap_t() directly, once per peptide per
## treatment-vs-reference comparison -- mirroring exactly the per-peptide /
## multiple-treatment-levels-vs-one-reference / BH-FDR-per-comparison logic
## pepdiff::compare(compare = "timepoint", ref = "0h") uses internally
## (see compare_pairwise()'s loop over treatment_levels), just with n_boot
## set explicitly.
##
## Fold change on VSN-scale data
## -------------------------------
## pepdiff's built-in fold_change formula is a ratio of means
## (treatment_mean / control_mean), which assumes linear-scale abundance
## input. VSN-normalized values are already on a roughly log2 scale, so a
## ratio of two VSN values is not a meaningful fold change. Here,
## log2_fc = mean(treatment) - mean(control) (a difference, since the
## inputs are already log2-like) and fold_change = 2^log2_fc.

library(pepdiff)
library(tidyverse)

# =============================================================================
# 1. Build a long Guy11-vs-Guy11 table from two VSN-normalized xlsx files
# =============================================================================

#' Read one experiment's VSN-normalized, imputed data and pivot the
#' shared-strain sample columns to long format.
#'
#' @param path Path to a `<experiment>_05_vsn_mar_mnar_imputed.xlsx` file
#' @param experiment_label Label to tag this experiment's rows with
#' @param strain Sample-name prefix identifying the shared strain (default "guy11")
#' @param peptide_col Column giving the unique peptide/site identifier
#'   (default "Modified Sequence" -- see header note; "Peptide Sequence"
#'   is NOT unique per row and must not be used here)
#' @param gene_col Column giving the gene identifier (default "Protein ID" --
#'   matches this project's convention elsewhere; the "Gene" column itself
#'   is empty/NA throughout this data)
#' @return Long tibble: peptide, gene_id, experiment, timepoint, bio_rep, value
read_guy11_long <- function(path, experiment_label, strain = "guy11",
                            peptide_col = "Modified Sequence", gene_col = "Protein ID") {
  df <- readxl::read_xlsx(path)

  # Match bare VSN intensity columns only (e.g. "guy11_0h_1"), not the
  # "guy11_0h_1 Spectral Count" / "MaxLFQ Intensity" / "Match Type" /
  # "Localization" columns that also start with the strain prefix.
  sample_pattern <- paste0("^", strain, "_([0-9]+(_[0-9]+)?h|spores)_[0-9]+$")
  sample_cols <- colnames(df)[stringr::str_detect(colnames(df), sample_pattern)]
  if (length(sample_cols) == 0) {
    stop("No bare '", strain, "_<timepoint>_<rep>' intensity columns found in ", path)
  }
  if (!peptide_col %in% colnames(df)) stop("Column not found: ", peptide_col)
  if (!gene_col %in% colnames(df)) stop("Column not found: ", gene_col)
  if (anyDuplicated(df[[peptide_col]])) {
    stop("'", peptide_col, "' is not unique per row in ", path,
         " -- pick a column that uniquely identifies each phosphosite")
  }

  df %>%
    dplyr::select(dplyr::all_of(c(peptide_col, gene_col, sample_cols))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(sample_cols), names_to = "sample", values_to = "value"
    ) %>%
    dplyr::mutate(
      # Sample names look like "guy11_0h_1" or "guy11_1_5h_2"; the 0h
      # timepoint is spelled "spores" in the pmk1 experiment's columns.
      timepoint_raw = stringr::str_extract(sample, "(?<=_)([0-9]+(_[0-9]+)?h|spores)(?=_)"),
      timepoint = dplyr::case_when(
        timepoint_raw == "spores" ~ "0h",
        TRUE ~ stringr::str_replace(timepoint_raw, "_", ".")
      ),
      bio_rep = stringr::str_extract(sample, "[0-9]+$"),
      experiment = experiment_label
    ) %>%
    dplyr::filter(!is.na(timepoint), !is.na(bio_rep)) %>%
    dplyr::rename(peptide = !!peptide_col, gene_id = !!gene_col) %>%
    dplyr::select(peptide, gene_id, experiment, timepoint, bio_rep, value)
}

#' Build the combined long-format Guy11-vs-Guy11 table, restricted to
#' phosphosites observed in both experiments.
#'
#' Peptides only seen in one experiment can never be correlated between
#' experiments and would just waste bootstrap iterations, so they're
#' dropped up front.
#'
#' @param path_a,path_b VSN-normalized, imputed xlsx paths for the two experiments
#' @param experiment_a,experiment_b Labels for the two experiments
#' @param ... Passed to [read_guy11_long()] (e.g. strain, peptide_col, gene_col)
#' @return Long tibble ready for [prepare_pepdiff_data()]
build_guy11_between_experiments_data <- function(path_a, path_b,
                                                  experiment_a, experiment_b,
                                                  ...) {
  long_a <- read_guy11_long(path_a, experiment_a, ...)
  long_b <- read_guy11_long(path_b, experiment_b, ...)

  shared_peptides <- intersect(unique(long_a$peptide), unique(long_b$peptide))
  message(sprintf(
    "%s: %d Guy11 phosphosites | %s: %d Guy11 phosphosites | shared: %d",
    experiment_a, dplyr::n_distinct(long_a$peptide),
    experiment_b, dplyr::n_distinct(long_b$peptide),
    length(shared_peptides)
  ))

  dplyr::bind_rows(long_a, long_b) %>%
    dplyr::filter(.data$peptide %in% shared_peptides)
}

# =============================================================================
# 2. Import into a pepdiff_data object
# =============================================================================

#' Write a long Guy11 comparison table to CSV and import it with
#' [pepdiff::read_pepdiff()].
#'
#' @param long_data Tibble from [build_guy11_between_experiments_data()]
#' @param csv_path Where to write the intermediate CSV (pepdiff reads from a file)
#' @return A pepdiff_data object
prepare_pepdiff_data <- function(long_data, csv_path) {
  dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(long_data, csv_path)
  pepdiff::read_pepdiff(
    csv_path,
    id = "peptide", gene = "gene_id", value = "value",
    factors = c("experiment", "timepoint"), replicate = "bio_rep"
  )
}

# =============================================================================
# 3. 10,000-iteration bootstrap-t: every timepoint vs a single reference
# =============================================================================

#' Run a bootstrap-t comparison of every non-reference level of one factor
#' against a single reference level (e.g. every timepoint vs "0h"), for
#' every peptide, calling pepdiff::test_bootstrap_t() directly so n_boot can
#' be set above pepdiff::compare()'s hardcoded default (see header note).
#'
#' No stratification: if you need to restrict to one experiment, subset
#' `data` first and call this once per subset (see run_guy11_bootstrap_diffexp.r).
#'
#' @param data A pepdiff_data object (or its `$data` tibble already subset
#'   to one experiment) -- either a pepdiff_data object or a plain tibble
#'   with the required columns is accepted
#' @param compare Factor to compare (default "timepoint")
#' @param ref Reference level of `compare` (default "0h")
#' @param n_boot Number of bootstrap iterations (default 10000)
#' @param alpha FDR significance threshold (default 0.05, matching pepdiff::compare()'s default)
#' @param fdr_method p.adjust method (default "BH", matching pepdiff::compare()'s default)
#' @param n_cores Parallel workers for [parallel::mclapply()] (default: all but one core)
#' @param seed Base RNG seed (L'Ecuyer-CMRG, for reproducible parallel bootstraps)
#' @param checkpoint_dir If not NULL, save each treatment-vs-ref comparison's
#'   results here as they finish (protects a multi-hour run against being
#'   lost partway through)
#' @param label Prefix for checkpoint filenames (e.g. an experiment name),
#'   so two calls sharing one checkpoint_dir don't collide
#' @return A tibble: peptide, gene_id, comparison, treatment, reference,
#'   n_ctrl, n_trt, log2_fc, fold_change, p_value, t_obs, fdr, significant
run_bootstrap_vs_ref <- function(data, compare = "timepoint", ref = "0h",
                                 n_boot = 10000, alpha = 0.05, fdr_method = "BH",
                                 n_cores = max(1, parallel::detectCores() - 1),
                                 seed = 1, checkpoint_dir = NULL, label = "run") {
  df <- if (inherits(data, "pepdiff_data")) data$data else data
  stopifnot(is.data.frame(df), all(c("peptide", "gene_id", compare, "value") %in% colnames(df)))

  all_levels <- unique(df[[compare]])
  trt_levels <- setdiff(all_levels, ref)
  peptides <- unique(df$peptide)

  if (!is.null(checkpoint_dir)) dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  all_results <- purrr::map_dfr(trt_levels, function(trt_level) {
    checkpoint_path <- if (!is.null(checkpoint_dir)) {
      file.path(checkpoint_dir, paste0(label, "_", trt_level, "_vs_", ref, ".rds"))
    } else NA_character_

    if (!is.na(checkpoint_path) && file.exists(checkpoint_path)) {
      message(sprintf("[%s] %s: '%s vs %s': loading existing checkpoint",
                      format(Sys.time(), "%H:%M:%S"), label, trt_level, ref))
      return(readRDS(checkpoint_path))
    }

    message(sprintf(
      "[%s] %s: '%s vs %s': %d peptides, n_boot = %d, %d core(s)",
      format(Sys.time(), "%H:%M:%S"), label, trt_level, ref, length(peptides), n_boot, n_cores
    ))
    t_start <- Sys.time()

    per_peptide <- parallel::mclapply(peptides, function(pep) {
      pep_rows <- df[df$peptide == pep, ]
      ctrl_values <- pep_rows$value[pep_rows[[compare]] == ref]
      trt_values <- pep_rows$value[pep_rows[[compare]] == trt_level]

      ctrl_values <- ctrl_values[!is.na(ctrl_values)]
      trt_values <- trt_values[!is.na(trt_values)]

      test_result <- pepdiff::test_bootstrap_t(ctrl_values, trt_values, n_boot = n_boot)

      # VSN-normalized values are already ~log2 scale: fold change is a
      # difference of means, not pepdiff's default ratio (see header note).
      log2_fc <- if (length(ctrl_values) > 0 && length(trt_values) > 0) {
        mean(trt_values) - mean(ctrl_values)
      } else {
        NA_real_
      }

      tibble::tibble(
        peptide = pep,
        gene_id = pep_rows$gene_id[1],
        comparison = paste(trt_level, "vs", ref),
        treatment = trt_level,
        reference = ref,
        n_ctrl = length(ctrl_values),
        n_trt = length(trt_values),
        log2_fc = log2_fc,
        fold_change = if (is.na(log2_fc)) NA_real_ else 2^log2_fc,
        p_value = test_result$p_value,
        t_obs = test_result$t_obs
      )
    }, mc.cores = n_cores, mc.set.seed = TRUE)

    comparison_results <- dplyr::bind_rows(per_peptide)

    elapsed <- as.numeric(Sys.time() - t_start, units = "mins")
    message(sprintf("  -> '%s vs %s' done in %.1f min", trt_level, ref, elapsed))

    if (!is.na(checkpoint_path)) saveRDS(comparison_results, checkpoint_path)

    comparison_results
  })

  # BH FDR correction within each comparison, matching how
  # pepdiff::compare_pairwise() applies FDR correction per comparison.
  all_results %>%
    dplyr::group_by(.data$comparison) %>%
    dplyr::mutate(
      fdr = stats::p.adjust(.data$p_value, method = fdr_method),
      significant = .data$fdr < alpha
    ) %>%
    dplyr::ungroup()
}

# =============================================================================
# 4. Significance filtering and between-experiment kinetics correlation
# =============================================================================

#' Peptides changing at least `fc_threshold`-fold (log2 scale) with
#' p < `p_threshold` at any comparison.
#'
#' Uses the raw (not FDR-corrected) p-value, as requested -- with ~16,000
#' peptides tested per comparison this will include some false positives;
#' the `fdr` column is still in `results` if a stricter cut is wanted later.
#'
#' @param results Output of [run_bootstrap_vs_ref()]
#' @param fc_threshold Minimum |log2 fold change| (default 1, i.e. 2-fold)
#' @param p_threshold Maximum p-value (default 0.05)
#' @return Character vector of unique peptide IDs
filter_significant_peptides <- function(results, fc_threshold = 1, p_threshold = 0.05) {
  results %>%
    dplyr::filter(abs(.data$log2_fc) >= fc_threshold, .data$p_value < p_threshold) %>%
    dplyr::pull(.data$peptide) %>%
    unique()
}

#' Correlate each peptide's mean-intensity trajectory across timepoints
#' between two experiments (Pearson r per peptide) -- i.e. do the two
#' experiments agree on the *shape* of a peptide's response over time, not
#' just the size of one fold change.
#'
#' Uses real VSN intensities (mean of bioreps at each timepoint), not
#' log2FC-vs-0h values, so there is no shared trivial zero point at 0h to
#' artificially inflate agreement (see script header note).
#'
#' @param long_data Long tibble from [build_guy11_between_experiments_data()]
#'   (columns: peptide, gene_id, experiment, timepoint, bio_rep, value)
#' @param peptides Character vector of peptide IDs to correlate (e.g. from
#'   [filter_significant_peptides()])
#' @param experiment_a,experiment_b Experiment labels to correlate against each other
#' @param timepoint_order Order of timepoints along the trajectory
#' @param min_points Minimum number of paired timepoints required to
#'   compute a correlation (default 3; Pearson r/p on fewer is unreliable)
#' @return Tibble: peptide, gene_id, n_points, r, p_value
correlate_peptide_kinetics <- function(long_data, peptides,
                                       experiment_a, experiment_b,
                                       timepoint_order = c("0h", "1h", "1.5h", "2h", "4h", "6h"),
                                       min_points = 3) {
  mean_by_tp <- long_data %>%
    dplyr::filter(.data$peptide %in% peptides) %>%
    dplyr::group_by(.data$peptide, .data$gene_id, .data$experiment, .data$timepoint) %>%
    dplyr::summarise(mean_value = mean(.data$value, na.rm = TRUE), .groups = "drop")

  wide <- mean_by_tp %>%
    tidyr::pivot_wider(
      id_cols = c(peptide, gene_id, timepoint),
      names_from = experiment, values_from = mean_value
    )

  purrr::map_dfr(peptides, function(pep) {
    pep_wide <- wide %>%
      dplyr::filter(.data$peptide == pep) %>%
      dplyr::mutate(timepoint = factor(.data$timepoint, levels = timepoint_order)) %>%
      dplyr::arrange(.data$timepoint) %>%
      dplyr::filter(!is.na(.data[[experiment_a]]), !is.na(.data[[experiment_b]]))

    n_points <- nrow(pep_wide)
    gene_id <- if (n_points > 0) pep_wide$gene_id[1] else NA_character_

    if (n_points < min_points) {
      return(tibble::tibble(
        peptide = pep, gene_id = gene_id, n_points = n_points, r = NA_real_, p_value = NA_real_
      ))
    }

    ct <- suppressWarnings(stats::cor.test(
      pep_wide[[experiment_a]], pep_wide[[experiment_b]], method = "pearson"
    ))
    tibble::tibble(
      peptide = pep, gene_id = gene_id, n_points = n_points,
      r = unname(ct$estimate), p_value = ct$p.value
    )
  })
}

# =============================================================================
# 5. Plotting helpers
# =============================================================================

#' Histogram of per-peptide kinetics correlation coefficients
#'
#' @param cor_results Output of [correlate_peptide_kinetics()]
#' @return ggplot object
plot_kinetics_correlation_distribution <- function(cor_results) {
  cor_results %>%
    dplyr::filter(!is.na(.data$r)) %>%
    ggplot(aes(x = r)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    theme_minimal() +
    labs(
      x = "Pearson r (peptide intensity trajectory, experiment A vs B)",
      y = "Number of peptides",
      title = "Between-experiment kinetics correlation of significantly-changed Guy11 peptides"
    )
}

#' Plot one peptide's intensity trajectory in both experiments side by side
#'
#' @param long_data Long tibble from [build_guy11_between_experiments_data()]
#' @param peptide_id The peptide ID to plot
#' @param timepoint_order Order of timepoints along the x-axis
#' @return ggplot object
plot_peptide_kinetics <- function(long_data, peptide_id,
                                  timepoint_order = c("0h", "1h", "1.5h", "2h", "4h", "6h")) {
  long_data %>%
    dplyr::filter(.data$peptide == peptide_id) %>%
    dplyr::mutate(timepoint = factor(.data$timepoint, levels = timepoint_order)) %>%
    ggplot(aes(x = timepoint, y = value, color = experiment, group = experiment)) +
    stat_summary(fun = mean, geom = "line", linewidth = 1) +
    stat_summary(fun = mean, geom = "point", size = 2) +
    geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
    theme_minimal() +
    labs(
      x = "Timepoint", y = "VSN-normalized intensity",
      title = paste0("Peptide kinetics: ", peptide_id), color = "Experiment"
    )
}
