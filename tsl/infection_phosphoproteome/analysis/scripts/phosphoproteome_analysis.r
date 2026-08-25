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

  total_row <- res %>%
    summarise(
      sample = "Total",
      sample_short = "Total",
      total_spectral_counts = sum(total_spectral_counts, na.rm = TRUE),
      total_maxlfq = sum(total_maxlfq, na.rm = TRUE),
      n_phospho_peptides = sum(n_phospho_peptides, na.rm = TRUE),
      n_phospho_proteins = sum(n_phospho_proteins, na.rm = TRUE)
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
