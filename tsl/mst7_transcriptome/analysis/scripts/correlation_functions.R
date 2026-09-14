required_pkgs <- c(
  "dplyr", "readr", "tibble", "purrr", "ggplot2", "tidyr", "here", "reshape2"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  warning("Missing packages for correlation_functions.R: ", paste(missing_pkgs, collapse = ", "))
}

library(dplyr)
library(readr)
library(tibble)
library(purrr)
library(ggplot2)
library(tidyr)
library(here)
library(reshape2)

source(file.path(here("analysis", "scripts"), "label_functions.R"))

load_tpm_matrix <- function(tpm_file = here("raw", "all_samples_tpm_matrix.txt")) {
  # TPM matrix is CSV format with gene IDs as first column
  read_csv(tpm_file, show_col_types = FALSE) %>%
    column_to_rownames(var = "target") %>%
    as.matrix()
}

# Collapse biological replicates: average TPM across bioreps for each
# name x timepoint combination, returning a condensed matrix (one column
# per condition/timepoint) plus matching metadata. This lets every subset
# function below (all/Exp1/Exp2/key_strains) be reused unchanged to produce
# a "mean of bioreps per timepoint" correlation matrix alongside the
# per-sample one.
collapse_bioreps <- function(tpm_matrix, s2c) {
  meta <- s2c %>%
    distinct(name, timepoint, experiment) %>%
    arrange(name, timepoint) %>%
    mutate(sample = paste0(name, "_", timepoint))

  mean_mat <- sapply(seq_len(nrow(meta)), function(i) {
    reps <- s2c %>%
      filter(name == meta$name[i], timepoint == meta$timepoint[i]) %>%
      pull(sample) %>%
      intersect(colnames(tpm_matrix))

    rowMeans(tpm_matrix[, reps, drop = FALSE], na.rm = TRUE)
  })
  colnames(mean_mat) <- meta$sample
  rownames(mean_mat) <- rownames(tpm_matrix)

  list(tpm_matrix = mean_mat, s2c = meta)
}

# Order samples by strain then timepoint (then biorep, if present) so that
# within-strain pairs form contiguous diagonal blocks in the heatmap and
# between-strain pairs form the off-diagonal blocks.
order_metadata <- function(meta) {
  if ("biorep" %in% names(meta)) {
    meta %>% arrange(name, timepoint, biorep)
  } else {
    meta %>% arrange(name, timepoint)
  }
}

# Filter out low-expression genes and log2-transform before correlating.
# Raw-TPM correlation is dominated by a handful of very highly expressed
# genes and understates true sample similarity (e.g. biological replicate
# pairs can come out below r = 0.2). This matches the approach validated
# in Correlation_matrix_RNAseq_V2.Rmd.
prepare_log2_matrix <- function(tpm_matrix, sample_names, min_mean_tpm = 1) {
  mat <- tpm_matrix[, sample_names, drop = FALSE]
  mean_tpm <- rowMeans(mat, na.rm = TRUE)
  mat <- mat[mean_tpm > min_mean_tpm, , drop = FALSE]
  log2(mat + 1)
}

compute_correlation_matrix <- function(tpm_matrix,
                                       sample_names,
                                       method = "pearson",
                                       min_mean_tpm = 1) {
  # TPM matrix is already a matrix, subset by column (order is preserved,
  # so callers control block ordering via the order of sample_names)
  mat <- prepare_log2_matrix(tpm_matrix, sample_names, min_mean_tpm)
  cor(mat, use = "pairwise.complete.obs", method = method)
}

plot_correlation_heatmap <- function(cor_mat,
                                     label,
                                     out_dir = here("results", "figures", "correlation"),
                                     subtitle = "Pearson correlation, log2(TPM + 1), genes with mean TPM > 1",
                                     show_values = FALSE,
                                     width = 10,
                                     height = 8,
                                     relabel_strains = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  df <- melt(cor_mat)
  # melt() puts the matrix's ROW identifier first (Var1) and its COLUMN
  # identifier second (Var2) -- map row->Sample2 (y-axis) and col->Sample1
  # (x-axis) accordingly. For the square matrices used everywhere else,
  # rows and columns are the same sample set so this distinction was
  # invisible; it matters once rows and columns differ (e.g. a rectangular
  # reference-strain-vs-everything matrix).
  colnames(df) <- c("Sample2", "Sample1", "Correlation")
  col_labels <- colnames(cor_mat)
  row_labels <- rownames(cor_mat)
  if (relabel_strains) {
    col_labels <- relabel_cross_experiment_strains(col_labels)
    row_labels <- relabel_cross_experiment_strains(row_labels)
    df$Sample1 <- relabel_cross_experiment_strains(as.character(df$Sample1))
    df$Sample2 <- relabel_cross_experiment_strains(as.character(df$Sample2))
  }
  df$Sample1 <- factor(df$Sample1, levels = col_labels)
  df$Sample2 <- factor(df$Sample2, levels = rev(row_labels))

  p <- ggplot(df, aes(x = Sample1, y = Sample2, fill = Correlation)) +
    geom_tile() +
    # Correlation here is bounded [0, 1], not [-1, 1] -- there is no
    # meaningful "opposite" value, so a diverging blue/red scale implies a
    # zero-crossing that doesn't exist. A sequential single-hue ramp (white
    # -> dark green) is the correct encoding for magnitude-only data.
    scale_fill_gradient(
      low = "#f7fcf5",
      high = "#00441b",
      limits = c(0, 1)
    ) +
    # Exact nomenclature (Delta + italics + superscript) for the five
    # strains that have it -- independent of relabel_strains above, and
    # applied regardless of whether this is a cross-experiment figure.
    scale_x_discrete(labels = strain_nomenclature_labels(col_labels)) +
    scale_y_discrete(labels = strain_nomenclature_labels(rev(row_labels))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(paste("Correlation Matrix:", label), subtitle = subtitle)

  if (show_values) {
    # White text reads on the dark-green end of the ramp, black on the
    # light end; switch around the ramp's visual midpoint.
    df$label_colour <- if_else(df$Correlation > 0.55, "white", "black")
    p <- p + geom_text(
      data = df,
      aes(label = sprintf("%.2f", Correlation), colour = label_colour),
      size = 3
    ) + scale_colour_identity()
  }

  out_file <- file.path(out_dir, paste0("correlation_", label, ".png"))
  ggsave(out_file, p, width = width, height = height, dpi = 300)

  invisible(p)
}

write_correlation_table <- function(cor_mat,
                                    label,
                                    out_dir = here("results", "tables", "correlation")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  df_out <- as.data.frame(cor_mat) %>%
    tibble::rownames_to_column("sample")

  out_file <- file.path(out_dir, paste0("correlation_", label, ".csv"))
  write_csv(df_out, out_file)

  invisible(df_out)
}

# Classify every pair of columns in a correlation matrix as within-strain
# (biological replicates, or same strain at different timepoints) or
# between-strain (different strain, at the same or different timepoint),
# and summarise correlation by category. Writes both the detailed pairwise
# table and the per-category summary to out_dir.
summarise_correlation_by_type <- function(cor_mat,
                                          s2c,
                                          label,
                                          out_dir = here("results", "tables", "correlation")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  has_biorep <- "biorep" %in% names(s2c)
  meta_cols <- c("sample", "name", "timepoint", if (has_biorep) "biorep")
  meta <- s2c %>%
    filter(sample %in% colnames(cor_mat)) %>%
    distinct(across(all_of(meta_cols)))

  pairs <- as.data.frame(cor_mat) %>%
    tibble::rownames_to_column("sample1") %>%
    pivot_longer(-sample1, names_to = "sample2", values_to = "correlation") %>%
    filter(sample1 < sample2) %>%
    left_join(meta, by = c("sample1" = "sample")) %>%
    rename(name1 = name, timepoint1 = timepoint) %>%
    left_join(meta, by = c("sample2" = "sample"), suffix = c("1", "2")) %>%
    rename(name2 = name, timepoint2 = timepoint)

  pairs <- pairs %>%
    mutate(
      same_strain = name1 == name2,
      same_time = timepoint1 == timepoint2,
      pair_type = case_when(
        has_biorep & same_strain & same_time ~ "Within-strain: biological replicates",
        same_strain & !same_time ~ "Within-strain: different timepoint",
        !same_strain & same_time ~ "Between-strain: same timepoint",
        TRUE ~ "Between-strain: different timepoint"
      )
    )

  detail_file <- file.path(out_dir, paste0("correlation_pairs_", label, ".csv"))
  write_csv(pairs, detail_file)

  plot_correlation_by_type(pairs, label)

  summary_tbl <- pairs %>%
    group_by(pair_type) %>%
    summarise(
      n_pairs = n(),
      mean_cor = mean(correlation, na.rm = TRUE),
      median_cor = median(correlation, na.rm = TRUE),
      sd_cor = sd(correlation, na.rm = TRUE),
      min_cor = min(correlation, na.rm = TRUE),
      max_cor = max(correlation, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_cor))

  summary_file <- file.path(out_dir, paste0("correlation_by_type_", label, ".csv"))
  write_csv(summary_tbl, summary_file)

  invisible(summary_tbl)
}

# Visualise whether correlation is driven more by strain identity or by
# timepoint: a 2x2 grouped boxplot of every sample pair, split by whether
# the pair shares a strain (x-axis) and whether it shares a timepoint
# (fill). "Within-strain, same timepoint" is always a biological-replicate
# pair (same strain + same timepoint + different sample), so this is the
# same categorisation as summarise_correlation_by_type() above, just laid
# out as a factorial grid instead of four flat categories.
plot_correlation_by_type <- function(pairs,
                                     label,
                                     out_dir = here("results", "figures", "correlation")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot_df <- pairs %>%
    mutate(
      strain_rel = factor(
        if_else(same_strain, "Within-strain", "Between-strain"),
        levels = c("Within-strain", "Between-strain")
      ),
      time_rel = factor(
        if_else(same_time, "Same timepoint", "Different timepoint"),
        levels = c("Same timepoint", "Different timepoint")
      )
    )

  p <- ggplot(plot_df, aes(x = strain_rel, y = correlation, fill = time_rel)) +
    geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.4, width = 0.6) +
    scale_fill_manual(values = c(
      "Same timepoint" = "#2a78d6",
      "Different timepoint" = "#eb6834"
    )) +
    labs(
      title = paste("Strain vs. timepoint effect on correlation:", label),
      subtitle = "Every sample pair, split by whether it shares a strain and/or a timepoint",
      x = NULL,
      y = "Pearson correlation (log2 TPM)",
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")

  out_file <- file.path(out_dir, paste0("correlation_by_type_", label, ".png"))
  ggsave(out_file, p, width = 7, height = 6, dpi = 300)

  invisible(p)
}

# One strain's conditions (rows) against every strain's conditions
# (columns) -- a rectangular slice rather than the full square matrix, with
# the correlation value printed in each tile. Both sides share the same
# gene filter (computed once, across every sample used on either axis) so
# the two axes are on the same footing before being split into rows/cols.
cor_matrix_reference_strain <- function(tpm_matrix,
                                        s2c,
                                        reference,
                                        strains = NULL,
                                        label = paste0(reference, "_reference_mean_biorep"),
                                        min_mean_tpm = 1,
                                        relabel_strains = TRUE) {
  ordered <- order_metadata(s2c)
  if (!is.null(strains)) {
    ordered <- ordered %>% filter(name %in% union(strains, reference))
  }
  all_samples <- ordered %>% pull(sample)
  ref_samples <- ordered %>% filter(name == reference) %>% pull(sample)

  log2_mat <- prepare_log2_matrix(tpm_matrix, all_samples, min_mean_tpm)
  cor_mat <- cor(
    log2_mat[, ref_samples, drop = FALSE],
    log2_mat[, all_samples, drop = FALSE],
    use = "pairwise.complete.obs"
  )

  plot_correlation_heatmap(cor_mat, label, show_values = TRUE, width = 12, height = 4, relabel_strains = relabel_strains)
  write_correlation_table(cor_mat, label)
  cor_mat
}

cor_matrix_all <- function(tpm_matrix, s2c, label = "all_samples", min_mean_tpm = 1, show_values = FALSE) {
  sample_names <- order_metadata(s2c) %>% pull(sample)
  cor_mat <- compute_correlation_matrix(tpm_matrix, sample_names, min_mean_tpm = min_mean_tpm)

  plot_correlation_heatmap(cor_mat, label, show_values = show_values, relabel_strains = TRUE)
  write_correlation_table(cor_mat, label)
  summarise_correlation_by_type(cor_mat, s2c, label)
  cor_mat
}

cor_matrix_exp1 <- function(tpm_matrix, s2c, label = "Exp1", min_mean_tpm = 1, show_values = FALSE) {
  s2c_sub <- s2c %>% filter(experiment == "Exp1")
  sample_names <- order_metadata(s2c_sub) %>% pull(sample)

  cor_mat <- compute_correlation_matrix(tpm_matrix, sample_names, min_mean_tpm = min_mean_tpm)

  plot_correlation_heatmap(cor_mat, label, show_values = show_values)
  write_correlation_table(cor_mat, label)
  summarise_correlation_by_type(cor_mat, s2c_sub, label)
  cor_mat
}

cor_matrix_exp2 <- function(tpm_matrix, s2c, label = "Exp2", min_mean_tpm = 1, show_values = FALSE) {
  s2c_sub <- s2c %>% filter(experiment == "Exp2")
  sample_names <- order_metadata(s2c_sub) %>% pull(sample)

  cor_mat <- compute_correlation_matrix(tpm_matrix, sample_names, min_mean_tpm = min_mean_tpm)

  plot_correlation_heatmap(cor_mat, label, show_values = show_values)
  write_correlation_table(cor_mat, label)
  summarise_correlation_by_type(cor_mat, s2c_sub, label)
  cor_mat
}

cor_matrix_key_strains <- function(tpm_matrix, s2c, label = "key_strains", min_mean_tpm = 1, show_values = FALSE) {
  key <- c("Guy11", "mst7", "Guy11M", "pmk1")

  s2c_sub <- s2c %>% filter(name %in% key)
  sample_names <- order_metadata(s2c_sub) %>% pull(sample)

  cor_mat <- compute_correlation_matrix(tpm_matrix, sample_names, min_mean_tpm = min_mean_tpm)

  plot_correlation_heatmap(cor_mat, label, show_values = show_values, width = 12, height = 10, relabel_strains = TRUE)
  write_correlation_table(cor_mat, label)
  summarise_correlation_by_type(cor_mat, s2c_sub, label)
  cor_mat
}
