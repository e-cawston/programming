## Phosphoproteome analysis helpers
## Summaries, plotting helpers, and venn set computations

library(tidyverse)
library(ggplot2)

#' Summarise phospho metrics per sample
#'
#' @param df Data frame (phospho filtered)
#' @param spec_cols Character vector of spectral count columns
#' @param lfq_cols Character vector of MaxLFQ columns
#' @param peptide_col Column name for peptide identifier (string)
#' @param protein_col Column name for protein identifier (string)
#' @return tibble with one row per sample and summary metrics
summarise_phospho_samples <- function(df, spec_cols, lfq_cols = NULL, peptide_col = "Peptide Sequence", protein_col = "Protein ID") {
  # Ensure columns exist
  missing_spec <- setdiff(spec_cols, colnames(df))
  if (length(missing_spec) > 0) stop("Spec columns missing: ", paste(missing_spec, collapse = ", "))

  # For each spectral column compute totals and counts
  res_spec <- purrr::map_dfr(spec_cols, function(col) {
    vals <- as.numeric(df[[col]])
    tibble(
      sample = col,
      sample_key = stringr::str_remove(col, "\\s+Spectral Count$"),
      total_spectral_counts = sum(vals, na.rm = TRUE),
      n_phospho_peptides = df %>% filter(replace_na(as.numeric(.data[[col]]), 0) > 0) %>% pull(!!rlang::sym(peptide_col)) %>% unique() %>% length(),
      n_phospho_proteins = df %>% filter(replace_na(as.numeric(.data[[col]]), 0) > 0) %>% pull(!!rlang::sym(protein_col)) %>% unique() %>% length()
    )
  })

  # Add MaxLFQ totals if provided
  if (!is.null(lfq_cols) && length(lfq_cols) > 0) {
    missing_lfq <- setdiff(lfq_cols, colnames(df))
    if (length(missing_lfq) > 0) stop("LFQ columns missing: ", paste(missing_lfq, collapse = ", "))
    res_lfq <- purrr::map_dfr(lfq_cols, function(col) {
      tibble(
        sample_key = stringr::str_remove(col, "\\s+MaxLFQ Intensity$"),
        total_maxlfq = sum(as.numeric(df[[col]]), na.rm = TRUE)
      )
    })
    # Spectral and MaxLFQ columns use different suffixes for the same sample.
    res <- dplyr::left_join(res_spec, res_lfq, by = "sample_key")
  } else {
    res <- res_spec %>% mutate(total_maxlfq = NA_real_)
  }

  # Clean sample labels (shorten by removing trailing ' Spectral Count' or ' MaxLFQ Intensity')
  res <- res %>%
    mutate(
      sample_short = stringr::str_remove(sample, "\\s+Spectral Count$"),
      sample_short = stringr::str_remove(sample_short, "\\s+MaxLFQ Intensity$")
    ) %>%
    select(sample, sample_short, total_spectral_counts, total_maxlfq,
           n_phospho_peptides, n_phospho_proteins)

  observed_in_any_sample <- df %>%
    filter(if_any(all_of(spec_cols), ~ replace_na(as.numeric(.x), 0) > 0))

  total_row <- res %>%
    summarise(
      sample = "Total",
      sample_short = "Total",
      total_spectral_counts = sum(total_spectral_counts, na.rm = TRUE),
      total_maxlfq = sum(total_maxlfq, na.rm = TRUE),
      n_phospho_peptides = observed_in_any_sample %>%
        pull(!!rlang::sym(peptide_col)) %>%
        unique() %>%
        length(),
      n_phospho_proteins = observed_in_any_sample %>%
        pull(!!rlang::sym(protein_col)) %>%
        unique() %>%
        length()
    )

  bind_rows(res, total_row)
}

#' Generic plotting helpers: return ggplot object
plot_total_spectral_counts <- function(summary_df) {
  ggplot(filter(summary_df, sample_short != "Total"),
         aes(x = sample_short, y = total_spectral_counts)) +
    geom_col(fill = "steelblue") +
    theme_minimal() +
    labs(x = "Sample", y = "Total spectral counts", title = "Total spectral counts per sample") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_total_maxlfq <- function(summary_df) {
  ggplot(filter(summary_df, sample_short != "Total"),
         aes(x = sample_short, y = total_maxlfq)) +
    geom_col(fill = "darkgreen") +
    theme_minimal() +
    labs(x = "Sample", y = "Total MaxLFQ intensity", title = "Total MaxLFQ per sample") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_total_phospho_peptides <- function(summary_df) {
  ggplot(filter(summary_df, sample_short != "Total"),
         aes(x = sample_short, y = n_phospho_peptides)) +
    geom_col(fill = "purple") +
    theme_minimal() +
    labs(x = "Sample", y = "Phosphorylated peptides", title = "Number of phosphorylated peptides per sample") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_total_phospho_proteins <- function(summary_df) {
  ggplot(filter(summary_df, sample_short != "Total"),
         aes(x = sample_short, y = n_phospho_proteins)) +
    geom_col(fill = "orange") +
    theme_minimal() +
    labs(x = "Sample", y = "Phosphorylated proteins", title = "Number of phosphorylated proteins per sample") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Unique counts plotting (expects columns `unique_peptides` / `unique_proteins` or similar)
plot_unique_phospho_peptides <- function(summary_df, value_col = "unique_peptides") {
  summary_df <- summary_df %>% rename(value = !!rlang::sym(value_col))
  ggplot(summary_df, aes(x = sample_short, y = value)) +
    geom_col(fill = "darkred") +
    theme_minimal() +
    labs(x = "Sample", y = "Unique phosphorylated peptides", title = "Unique phosphorylated peptides per sample") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_unique_phospho_proteins <- function(summary_df, value_col = "unique_proteins") {
  summary_df <- summary_df %>% rename(value = !!rlang::sym(value_col))
  ggplot(summary_df, aes(x = sample_short, y = value)) +
    geom_col(fill = "darkcyan") +
    theme_minimal() +
    labs(x = "Sample", y = "Unique phosphorylated proteins", title = "Unique phosphorylated proteins per sample") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Summarise one data-processing stage per sample.
summarise_phospho_stage <- function(df, experiment, stage,
                                    spec_cols = character(),
                                    lfq_cols = character(),
                                    intensity_cols = character(),
                                    peptide_col = "Peptide Sequence",
                                    protein_col = "Protein ID") {
  if (!all(c(peptide_col, protein_col) %in% colnames(df))) {
    stop("Peptide and protein identifier columns are required")
  }

  available <- function(columns) columns[columns %in% colnames(df)]
  spec_cols <- available(spec_cols)
  lfq_cols <- available(lfq_cols)
  intensity_cols <- available(intensity_cols)

  if (length(spec_cols) == 0 && length(lfq_cols) == 0 && length(intensity_cols) == 0) {
    stop("No measurement columns found for stage: ", stage)
  }

  sample_names <- unique(c(
    stringr::str_remove(spec_cols, "\\s+Spectral Count$"),
    stringr::str_remove(lfq_cols, "\\s+MaxLFQ Intensity$"),
    stringr::str_remove(intensity_cols, "\\s+Intensity$")
  ))

  result <- purrr::map_dfr(sample_names, function(sample_name) {
    spec_col <- spec_cols[stringr::str_remove(spec_cols, "\\s+Spectral Count$") == sample_name]
    lfq_col <- lfq_cols[stringr::str_remove(lfq_cols, "\\s+MaxLFQ Intensity$") == sample_name]
    intensity_col <- intensity_cols[stringr::str_remove(intensity_cols, "\\s+Intensity$") == sample_name]

    if (stage %in% c("normalized", "imputed")) {
      detected <- !is.na(as.numeric(df[[intensity_col[[1]]]]))
    } else if (length(intensity_col) > 0) {
      detected <- replace(as.numeric(df[[intensity_col[[1]]]]),
                          is.na(df[[intensity_col[[1]]]]), 0) > 0
    } else {
      detected <- replace(as.numeric(df[[spec_col[[1]]]]),
                          is.na(df[[spec_col[[1]]]]), 0) > 0
    }

    tibble::tibble(
      experiment = experiment,
      stage = stage,
      sample = sample_name,
      sample_short = sample_name,
      total_spectral_counts = if (length(spec_col) > 0) sum(as.numeric(df[[spec_col[[1]]]]), na.rm = TRUE) else NA_real_,
      total_maxlfq = if (length(lfq_col) > 0) sum(as.numeric(df[[lfq_col[[1]]]]), na.rm = TRUE) else NA_real_,
      total_intensity = if (length(intensity_col) > 0) sum(as.numeric(df[[intensity_col[[1]]]]), na.rm = TRUE) else NA_real_,
      n_phospho_peptides = dplyr::n_distinct(df[[peptide_col]][detected], na.rm = TRUE),
      n_phospho_proteins = dplyr::n_distinct(df[[protein_col]][detected], na.rm = TRUE)
    )
  })

  result
}

#' Create one overall row per experiment and processing stage.
overall_phospho_summary <- function(summary_df) {
  summary_df %>%
    group_by(experiment, stage) %>%
    summarise(
      total_spectral_counts = sum(total_spectral_counts, na.rm = TRUE),
      total_maxlfq = sum(total_maxlfq, na.rm = TRUE),
      total_intensity = sum(total_intensity, na.rm = TRUE),
      total_phospho_peptides = sum(n_phospho_peptides, na.rm = TRUE),
      total_phospho_proteins = sum(n_phospho_proteins, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(experiment, factor(stage, levels = c("filtered", "normalized", "imputed")))
}

#' Bar chart of one summary metric, faceted by stage x experiment.
#'
#' Stages get their own facet row (rather than being dodged within one
#' panel) with a free y-scale per row: metrics like total_intensity are on
#' wildly different scales before vs after VSN normalisation (raw
#' intensities vs a roughly log2 scale), so sharing one y-axis across
#' stages would make the normalized/imputed bars invisible.
plot_summary_metric <- function(summary_df, metric, y_label, title) {
  stage_levels <- intersect(c("filtered", "normalized", "imputed", "phospho_only"),
                            unique(summary_df$stage))
  summary_df <- summary_df %>%
    dplyr::mutate(stage = factor(stage, levels = stage_levels))

  ggplot(summary_df, aes(x = sample_short, y = .data[[metric]], fill = stage)) +
    geom_col() +
    facet_grid(stage ~ experiment, scales = "free") +
    theme_minimal() +
    labs(x = "Sample", y = y_label, title = title, fill = "Stage") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_stage_spectral_counts <- function(summary_df) {
  plot_summary_metric(summary_df, "total_spectral_counts", "Total spectral counts", "Spectral counts by processing stage")
}

plot_stage_maxlfq <- function(summary_df) {
  plot_summary_metric(summary_df, "total_maxlfq", "Total MaxLFQ intensity", "MaxLFQ intensity by processing stage")
}

plot_stage_intensity <- function(summary_df) {
  plot_summary_metric(summary_df, "total_intensity", "Total intensity", "Intensity by processing stage")
}

plot_stage_phospho_peptides <- function(summary_df) {
  plot_summary_metric(summary_df, "n_phospho_peptides", "Detected phosphorylated peptides", "Phosphorylated peptides by processing stage")
}

plot_stage_phospho_proteins <- function(summary_df) {
  plot_summary_metric(summary_df, "n_phospho_proteins", "Detected phosphorylated proteins", "Phosphorylated proteins by processing stage")
}

#' Boxplot of peptide intensities per sample
#'
#' One point is plotted per peptide's intensity value in that sample; the
#' boxplot summarises the distribution of those values for each sample,
#' one box per sample (consistent with the per-sample layout of the other
#' QC plots, e.g. `plot_total_spectral_counts`). Intended to be run on the
#' same dataset before and after VSN normalisation/imputation so the two
#' stages can be compared side by side.
#'
#' @param df Data frame containing intensity columns (raw/filtered or
#'   VSN-normalized/imputed).
#' @param intensity_cols Character vector of intensity column names to plot.
#' @param intensity_suffix Suffix to strip from `intensity_cols` to recover
#'   sample names, e.g. " Intensity" for the filtered stage or "" once VSN
#'   normalisation has renamed columns to bare sample names.
#' @param peptide_col Column identifying each peptide (default "Peptide Sequence").
#' @param drop_nonpositive Drop intensities <= 0 before plotting (these mark
#'   "not observed" on the raw/filtered scale, but should stay FALSE for
#'   VSN-scale data where 0/negative values are legitimate).
#' @param log_y Plot the y-axis on a log10 scale. Defaults to whatever
#'   `drop_nonpositive` is, since raw/filtered intensities span several
#'   orders of magnitude (log scale needed) while VSN-scale data is already
#'   roughly log-transformed (linear scale is appropriate).
#' @param title Plot title.
#' @param y_label Y-axis label.
#' @return ggplot object
plot_intensity_boxplot_by_sample <- function(df, intensity_cols,
                                             intensity_suffix = " Intensity",
                                             peptide_col = "Peptide Sequence",
                                             drop_nonpositive = TRUE,
                                             log_y = drop_nonpositive,
                                             title = "Peptide intensity distribution by sample",
                                             y_label = "Intensity") {
  missing_cols <- setdiff(intensity_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Intensity columns missing: ", paste(missing_cols, collapse = ", "))
  }
  if (!peptide_col %in% colnames(df)) stop("Peptide column missing: ", peptide_col)

  long <- df %>%
    dplyr::select(dplyr::all_of(c(peptide_col, intensity_cols))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(intensity_cols),
      names_to = "sample", values_to = "intensity"
    ) %>%
    dplyr::mutate(
      intensity = as.numeric(intensity),
      sample = if (nzchar(intensity_suffix)) {
        stringr::str_remove(sample, stringr::fixed(intensity_suffix))
      } else {
        sample
      }
    ) %>%
    dplyr::filter(!is.na(intensity), is.finite(intensity))

  if (drop_nonpositive) long <- dplyr::filter(long, intensity > 0)

  if (nrow(long) == 0) {
    stop("No plottable intensity values found; check intensity_suffix")
  }

  # Keep samples in the same order as intensity_cols (matches the other
  # per-sample QC plots, which use column order rather than chronological
  # timepoint order).
  sample_levels <- if (nzchar(intensity_suffix)) {
    stringr::str_remove(intensity_cols, stringr::fixed(intensity_suffix))
  } else {
    intensity_cols
  }
  long <- long %>%
    dplyr::mutate(sample = factor(sample, levels = unique(sample_levels)))

  p <- ggplot(long, aes(x = sample, y = intensity)) +
    geom_jitter(width = 0.15, alpha = 0.06, size = 0.4, color = "grey40") +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, fill = "steelblue", width = 0.5) +
    theme_minimal() +
    labs(x = "Sample", y = y_label, title = title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (log_y) {
    p <- p + scale_y_log10(labels = scales::label_comma())
  }

  p
}

#' Venn computations for peptides
#'
#' @param df_a Data frame A (filtered phospho dataset)
#' @param df_b Data frame B
#' @param peptide_col Name of peptide column
#' @return list(a_only, b_only, overlap)
venn_phospho_peptides <- function(df_a, df_b, peptide_col = "Peptide Sequence") {
  a_set <- df_a %>% pull(!!rlang::sym(peptide_col)) %>% unique() %>% na.omit()
  b_set <- df_b %>% pull(!!rlang::sym(peptide_col)) %>% unique() %>% na.omit()

  overlap <- intersect(a_set, b_set)
  a_only <- setdiff(a_set, b_set)
  b_only <- setdiff(b_set, a_set)

  list(a_only = tibble(peptide = a_only), b_only = tibble(peptide = b_only), overlap = tibble(peptide = overlap))
}

#' Venn computations for proteins
#'
#' @param df_a Data frame A
#' @param df_b Data frame B
#' @param protein_col Name of protein identifier column
#' @return list(a_only, b_only, overlap)
venn_phospho_proteins <- function(df_a, df_b, protein_col = "Protein ID") {
  a_set <- df_a %>% pull(!!rlang::sym(protein_col)) %>% unique() %>% na.omit()
  b_set <- df_b %>% pull(!!rlang::sym(protein_col)) %>% unique() %>% na.omit()

  overlap <- intersect(a_set, b_set)
  a_only <- setdiff(a_set, b_set)
  b_only <- setdiff(b_set, a_set)

  list(a_only = tibble(protein = a_only), b_only = tibble(protein = b_only), overlap = tibble(protein = overlap))
}

#' Run PCA on samples using phosphosite intensity profiles.
#'
#' Rows are features and columns are samples in `data`; the PCA treats
#' samples as observations and phosphosites as variables. Input is expected
#' to be VSN-normalized (imputed) abundance data. Features with non-finite
#' values or zero variance are excluded; features are not individually
#' scaled because the VSN values are already on a normalized scale.
#'
#' @param data Data frame with phosphosite metadata and intensity columns.
#' @param sample_cols Bare sample intensity column names.
#' @param experiment Label used for output names and plot title.
#' @param output_dir Directory for PCA figures.
#' @param workbook_path Path for the workbook containing PCA scores,
#'   variance explained, and feature loadings.
#' @return List containing the prcomp fit, scores, variance table,
#'   loadings, and score/scree plots.
run_phosphoproteome_pca <- function(data, sample_cols, experiment,
                                    output_dir = "results/figures/pca",
                                    workbook_path = file.path(
                                      "results", paste0(experiment, "_pca.xlsx")
                                    )) {
  missing_cols <- setdiff(sample_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("PCA sample columns missing: ", paste(missing_cols, collapse = ", "))
  }
  if (length(sample_cols) < 3) {
    stop("PCA requires at least three sample columns")
  }
  if (length(experiment) != 1 || is.na(experiment) || !nzchar(experiment)) {
    stop("experiment must be a non-empty string")
  }

  intensity_matrix <- vapply(
    data[sample_cols], function(x) as.numeric(x), numeric(nrow(data))
  )
  dimnames(intensity_matrix) <- list(NULL, sample_cols)
  complete_features <- apply(is.finite(intensity_matrix), 1, all)
  feature_variance <- apply(intensity_matrix, 1, stats::var, na.rm = TRUE)
  keep_features <- complete_features & is.finite(feature_variance) & feature_variance > 0
  if (sum(keep_features) < 2) {
    stop("PCA requires at least two complete, variable phosphosite features")
  }

  pca <- stats::prcomp(
    t(intensity_matrix[keep_features, , drop = FALSE]),
    center = TRUE, scale. = FALSE
  )
  if (ncol(pca$x) < 2) stop("PCA produced fewer than two components")

  sample_parts <- stringr::str_match(
    sample_cols, "^(.+?)_((?:[0-9]+(?:_[0-9]+)?h)|spores)_([0-9]+)$"
  )
  if (anyNA(sample_parts[, 2])) {
    stop("PCA sample names must contain genotype, timepoint, and replicate: ",
         paste(sample_cols[is.na(sample_parts[, 2])], collapse = ", "))
  }
  sample_timepoint <- ifelse(
    sample_parts[, 3] == "spores",
    "0h",
    stringr::str_replace(sample_parts[, 3], "_", ".")
  )
  timepoint <- factor(
    sample_timepoint, levels = c("0h", "1h", "1.5h", "2h", "4h", "6h")
  )
  if (anyNA(timepoint)) {
    stop("PCA sample names contain unsupported timepoints: ",
         paste(sample_cols[is.na(timepoint)], collapse = ", "))
  }
  score_data <- tibble::tibble(
    sample = sample_cols,
    genotype = sample_parts[, 2],
    timepoint = timepoint,
    biological_replicate = sample_parts[, 4],
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )

  eigenvalues <- pca$sdev^2
  variance_data <- tibble::tibble(
    component = paste0("PC", seq_along(eigenvalues)),
    variance = eigenvalues,
    proportion = eigenvalues / sum(eigenvalues),
    cumulative_proportion = cumsum(eigenvalues / sum(eigenvalues))
  )
  variance_label <- function(component) {
    sprintf("%s (%.1f%%)", component, 100 * variance_data$proportion[
      match(component, variance_data$component)
    ])
  }

  feature_metadata_cols <- intersect(
    c("Unique_site_ID", "Peptide Sequence", "Modified Sequence",
      "Assigned Modifications", "Protein ID"),
    colnames(data)
  )
  feature_metadata <- data[keep_features, feature_metadata_cols, drop = FALSE]
  feature_ids <- if ("Unique_site_ID" %in% colnames(feature_metadata)) {
    as.character(feature_metadata[["Unique_site_ID"]])
  } else {
    paste0("feature_", which(keep_features))
  }
  feature_ids[is.na(feature_ids) | !nzchar(feature_ids)] <-
    paste0("feature_", which(is.na(feature_ids) | !nzchar(feature_ids)))
  feature_ids <- make.unique(feature_ids)
  loadings <- data.frame(
    feature_id = feature_ids,
    feature_metadata,
    pca$rotation,
    check.names = FALSE
  )

  score_plot <- ggplot2::ggplot(
    score_data,
    ggplot2::aes(x = .data$PC1, y = .data$PC2,
                 color = .data$genotype, shape = .data$timepoint)
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::labs(
      title = paste(toupper(experiment), "PCA of phosphosite abundances"),
      x = variance_label("PC1"), y = variance_label("PC2"),
      color = "Genotype", shape = "Timepoint"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "right")

  scree_plot <- ggplot2::ggplot(
    variance_data,
    ggplot2::aes(x = factor(.data$component, levels = .data$component),
                 y = .data$proportion)
  ) +
    ggplot2::geom_col(fill = "#2F6690") +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$cumulative_proportion, group = 1),
      color = "#D55E00", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$cumulative_proportion),
      color = "#D55E00", size = 2
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(),
      sec.axis = ggplot2::sec_axis(~ ., labels = scales::label_percent(),
                                   name = "Cumulative variance explained")
    ) +
    ggplot2::labs(
      title = paste(toupper(experiment), "PCA variance explained"),
      x = "Principal component", y = "Variance explained"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(workbook_path), recursive = TRUE, showWarnings = FALSE)
  writexl::write_xlsx(
    list(scores = score_data, variance_explained = variance_data, loadings = loadings),
    workbook_path
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0(experiment, "_pca_scores.png")),
    score_plot, width = 9, height = 7, dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0(experiment, "_pca_scores.pdf")),
    score_plot, width = 9, height = 7, bg = "white"
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0(experiment, "_pca_variance.png")),
    scree_plot, width = 9, height = 6, dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0(experiment, "_pca_variance.pdf")),
    scree_plot, width = 9, height = 6, bg = "white"
  )

  message(experiment, " PCA: ", sum(keep_features), " features across ",
          length(sample_cols), " samples")
  list(
    fit = pca, scores = score_data, variance = variance_data,
    loadings = loadings, score_plot = score_plot, scree_plot = scree_plot
  )
}
