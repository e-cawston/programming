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

load_tpm_matrix <- function(tpm_file = here("raw", "all_samples_tpm_matrix.txt")) {
  # TPM matrix is CSV format with gene IDs as first column
  read_csv(tpm_file, show_col_types = FALSE) %>%
    column_to_rownames(var = "target") %>%
    as.matrix()
}

compute_correlation_matrix <- function(tpm_matrix, sample_names, method = "pearson") {
  # TPM matrix is already a matrix, subset by column
  mat <- tpm_matrix[, sample_names, drop = FALSE]
  cor(mat, use = "pairwise.complete.obs", method = method)
}

plot_correlation_heatmap <- function(cor_mat,
                                     label,
                                     out_dir = here("results", "figures", "correlation")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  df <- melt(cor_mat)
  colnames(df) <- c("Sample1", "Sample2", "Correlation")

  p <- ggplot(df, aes(x = Sample1, y = Sample2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue",
      high = "red",
      mid = "white",
      midpoint = 0.5,
      limit = c(0, 1)
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(paste("Correlation Matrix:", label))

  out_file <- file.path(out_dir, paste0("correlation_", label, ".png"))
  ggsave(out_file, p, width = 10, height = 8, dpi = 300)

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

cor_matrix_all <- function(tpm_matrix, s2c, label = "all_samples") {
  sample_names <- s2c$sample
  cor_mat <- compute_correlation_matrix(tpm_matrix, sample_names)

  plot_correlation_heatmap(cor_mat, label)
  write_correlation_table(cor_mat, label)
  cor_mat
}

cor_matrix_exp1 <- function(tpm_matrix, s2c, label = "Exp1") {
  sample_names <- s2c %>%
    filter(experiment == "Exp1") %>%
    pull(sample)

  cor_mat <- compute_correlation_matrix(tpm_matrix, sample_names)

  plot_correlation_heatmap(cor_mat, label)
  write_correlation_table(cor_mat, label)
  cor_mat
}

cor_matrix_exp2 <- function(tpm_matrix, s2c, label = "Exp2") {
  sample_names <- s2c %>%
    filter(experiment == "Exp2") %>%
    pull(sample)

  cor_mat <- compute_correlation_matrix(tpm_matrix, sample_names)

  plot_correlation_heatmap(cor_mat, label)
  write_correlation_table(cor_mat, label)
  cor_mat
}

cor_matrix_key_strains <- function(tpm_matrix, s2c, label = "key_strains") {
  key <- c("Guy11", "mst7", "Guy11M", "pmk1")

  sample_names <- s2c %>%
    filter(name %in% key) %>%
    pull(sample)

  cor_mat <- compute_correlation_matrix(tpm_matrix, sample_names)

  plot_correlation_heatmap(cor_mat, label)
  write_correlation_table(cor_mat, label)
  cor_mat
}
