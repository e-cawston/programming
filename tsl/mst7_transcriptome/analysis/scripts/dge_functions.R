required_pkgs <- c(
  "dplyr", "readr", "purrr", "stringr", "tibble", "here", "ggplot2", "sleuth"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  warning("Missing packages for dge_functions.R (some features may not work): ", paste(missing_pkgs, collapse = ", "))
}

suppress_library_error <- function(pkg) {
  tryCatch(library(pkg, character.only = TRUE), error = function(e) {
    warning("Could not load ", pkg, ": ", e$message)
  })
}

library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tibble)
library(here)

source(file.path(here("analysis", "scripts"), "label_functions.R"))
library(ggplot2)
suppress_library_error('sleuth')

# If sleuth didn't load, use stub functions
if (!exists('sleuth_prep')) {
  source(file.path(here('analysis'), 'scripts', 'stub_functions.R'), local = TRUE)
}

resolve_abundance_path <- function(path) {
  # sleuth_prep needs the *directory* containing a sample's kallisto output
  # (abundance.h5, abundance.tsv, run_info.json) -- not the abundance.tsv
  # file itself. run_metadata.txt stores paths pointing at abundance.tsv,
  # so those must be converted to their containing directory.
  if (is.null(path) || length(path) == 0 || is.na(path)) {
    return(NA_character_)
  }

  if (dir.exists(path)) {
    return(path)
  }

  if (file.exists(path) && basename(path) %in% c("abundance.tsv", "abundance.h5")) {
    return(dirname(path))
  }

  candidate <- file.path(here("raw"), path)
  if (dir.exists(candidate)) {
    return(candidate)
  }
  if (file.exists(candidate) && basename(candidate) %in% c("abundance.tsv", "abundance.h5")) {
    return(dirname(candidate))
  }

  path
}

prepare_comparison_subset <- function(s2c, comparison_row) {
  control <- comparison_row$control
  test <- comparison_row$test
  experiment <- comparison_row$experiment

  if (is.null(experiment) || is.na(experiment) || identical(experiment, "Cross")) {
    s2c_sub <- s2c %>%
      filter(name %in% c(control, test)) %>%
      droplevels()
  } else {
    s2c_sub <- s2c %>%
      filter(experiment == experiment) %>%
      filter(name %in% c(control, test)) %>%
      droplevels()
  }

  s2c_sub %>%
    mutate(path = map_chr(path, resolve_abundance_path)) %>%
    mutate(name = factor(name), timepoint = factor(timepoint)) %>%
    mutate(name = relevel(name, ref = control))
}

standardise_dge_result <- function(res, comparison_row, timepoint) {
  if (is.null(res) || nrow(res) == 0) {
    return(tibble())
  }

  res %>%
    mutate(
      comparison_label = comparison_row$label,
      timepoint = as.character(timepoint),
      control = comparison_row$control,
      test = comparison_row$test,
      comparison = paste(control, test, sep = "_vs_")
    )
}

run_comparison <- function(s2c,
                           comparison_row,
                           t2g = NULL,
                           aggregation_column = NULL,
                           out_dir = here("results", "dge"),
                           num_cores = 1) {
  label <- comparison_row$label
  control <- comparison_row$control
  test <- comparison_row$test

  message(">> Running comparison: ", label)

  s2c_sub <- prepare_comparison_subset(s2c, comparison_row)
  timepoints <- levels(s2c_sub$timepoint)

  results <- map(timepoints, function(tp) {
    s2c_tp <- s2c_sub %>% filter(timepoint == tp)

    if (nrow(s2c_tp) == 0) {
      message("   Timepoint ", tp, "H: no samples, skipping.")
      return(tibble())
    }

    message("   Timepoint ", tp, "H: running sleuth")

    so <- tryCatch(
      sleuth_prep(
        s2c_tp,
        ~ name,
        target_mapping = t2g,
        aggregation_column = aggregation_column,
        transformation_function = function(x) log2(x + 0.5),
        num_cores = num_cores,
        read_bootstrap_tpm = FALSE,
        extra_bootstrap_summary = FALSE
      ),
      error = function(e) {
        message("   sleuth_prep failed: ", e$message)
        NULL
      }
    )

    if (is.null(so)) {
      return(tibble())
    }

    so <- sleuth_fit(so, ~ name, "full")
    so <- sleuth_fit(so, ~ 1, "reduced")
    so <- sleuth_lrt(so, "reduced", "full")

    coef_name <- paste0("name", test)
    so <- sleuth_wt(so, coef_name)
    res <- sleuth_results(so, coef_name, "wt")
    res <- standardise_dge_result(res, comparison_row, tp)

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(out_dir, paste0("sleuth_", label, "_", tp, "H.csv"))
    write_csv(res, out_file)

    rm(so)
    gc()
    closeAllConnections()

    res
  })

  bind_rows(results)
}

run_comparisons <- function(s2c,
                            comparisons,
                            t2g = NULL,
                            aggregation_column = NULL,
                            out_dir = here("results", "dge"),
                            num_cores = 1) {
  map_dfr(seq_len(nrow(comparisons)), function(i) {
    comparison_row <- comparisons[i, , drop = FALSE]
    run_comparison(
      s2c = s2c,
      comparison_row = comparison_row,
      t2g = t2g,
      aggregation_column = aggregation_column,
      out_dir = out_dir,
      num_cores = num_cores
    )
  })
}

summarise_dge_counts <- function(dge_table,
                                 qval_threshold = 0.05,
                                 fc_threshold = 0,
                                 mirror = FALSE) {
  dge_table %>%
    filter(qval <= qval_threshold) %>%
    filter(abs(b) >= fc_threshold) %>%
    # mirror flips which side of zero "up"/"down" land on -- xor() with
    # mirror=TRUE inverts the b>0 test, giving the exact A_vs_B / B_vs_A
    # mirror image of the same underlying DEGs.
    mutate(direction = if_else(xor(b > 0, mirror), "up", "down")) %>%
    count(timepoint, direction, name = "n_deg") %>%
    mutate(
      n_plot = if_else(direction == "down", -n_deg, n_deg),
      timepoint_num = suppressWarnings(as.numeric(str_remove(timepoint, "H")))
    )
}

plot_dge_summary <- function(dge_table,
                             comparison_label,
                             out_dir = here("results", "figures", "dge"),
                             qval_threshold = 0.05,
                             fc_threshold = 0,
                             display_label = comparison_label,
                             mirror = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot_df <- dge_table %>%
    filter(comparison_label == !!comparison_label) %>%
    summarise_dge_counts(qval_threshold = qval_threshold, fc_threshold = fc_threshold, mirror = mirror) %>%
    # as.character() guards against timepoint coming back numeric when
    # dge_table has been re-read from a saved CSV rather than used live
    # from run_comparisons() -- fct_reorder() requires factor/character.
    mutate(timepoint = forcats::fct_reorder(as.character(timepoint), timepoint_num, .desc = FALSE))

  if (nrow(plot_df) == 0) {
    message("No DEGs to plot for comparison: ", comparison_label)
    return(NULL)
  }

  # A_vs_B mirrored is the B_vs_A framing -- swap the two halves of the
  # display label (not the internal comparison_label used for filtering
  # and the file name) so the title matches what's actually drawn.
  title_label <- display_label
  if (mirror) {
    parts <- strsplit(display_label, "_vs_", fixed = TRUE)[[1]]
    title_label <- if (length(parts) == 2) paste(rev(parts), collapse = "_vs_") else paste(display_label, "(mirrored)")
  }

  # Total unique DEGs across all timepoints -- a gene hitting threshold at
  # more than one timepoint is only counted once here, unlike the sum of
  # the per-bar counts above (which is a sum of per-timepoint events, not
  # a gene count).
  n_unique <- dge_table %>%
    filter(comparison_label == !!comparison_label, qval <= qval_threshold, abs(b) >= fc_threshold) %>%
    pull(target_id) %>%
    n_distinct()

  p <- ggplot(plot_df, aes(x = n_plot, y = timepoint, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    geom_text(
      aes(label = n_deg, hjust = if_else(direction == "down", 1.15, -0.15)),
      size = 3.2
    ) +
    scale_fill_manual(values = c(up = "#D62728", down = "#1F77B4")) +
    scale_x_continuous(labels = abs, expand = expansion(mult = 0.12)) +
    labs(
      title = strain_nomenclature_title(title_label, prefix = "DEG counts: "),
      subtitle = paste0("qval < ", qval_threshold, " and |log2FC| >= ", fc_threshold),
      x = "Number of DEGs",
      y = "Timepoint",
      fill = NULL,
      caption = paste0("Total unique DEGs across all timepoints: ", n_unique)
    ) +
    theme_classic(base_size = 12) +
    theme(plot.caption = element_text(hjust = 0.5, size = 10, face = "italic"))

  out_file <- file.path(out_dir, paste0("dge_summary_", comparison_label, if (mirror) "_mirrored" else "", ".png"))
  ggsave(out_file, p, width = 8, height = 5, dpi = 300)

  p
}

# Convenience: write both the forward (A_vs_B) and mirrored (B_vs_A)
# framing of one comparison's DEG summary in a single call, so both are
# available to pick from -- e.g. dge_summary_Guy11_vs_mst7.png (up = more
# in mst7) and dge_summary_Guy11_vs_mst7_mirrored.png (up = more in Guy11).
plot_dge_summary_both_directions <- function(dge_table,
                                             comparison_label,
                                             out_dir = here("results", "figures", "dge"),
                                             qval_threshold = 0.05,
                                             fc_threshold = 0,
                                             display_label = comparison_label) {
  list(
    forward = plot_dge_summary(
      dge_table, comparison_label, out_dir, qval_threshold, fc_threshold,
      display_label = display_label, mirror = FALSE
    ),
    mirrored = plot_dge_summary(
      dge_table, comparison_label, out_dir, qval_threshold, fc_threshold,
      display_label = display_label, mirror = TRUE
    )
  )
}

plot_dge_volcano <- function(dge_table,
                             comparison_label,
                             timepoint,
                             out_dir = here("results", "figures", "dge"),
                             qval_threshold = 0.05,
                             fc_threshold = 0) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot_df <- dge_table %>%
    filter(comparison_label == !!comparison_label) %>%
    filter(timepoint == !!timepoint) %>%
    mutate(
      status = case_when(
        qval <= qval_threshold & b >= fc_threshold ~ "up",
        qval <= qval_threshold & b <= -fc_threshold ~ "down",
        TRUE ~ "ns"
      )
    )

  if (nrow(plot_df) == 0) {
    return(NULL)
  }

  p <- ggplot(plot_df, aes(x = b, y = -log10(qval), colour = status)) +
    geom_point(alpha = 0.7, size = 1.2) +
    geom_hline(yintercept = -log10(qval_threshold), linetype = "dashed") +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed") +
    scale_colour_manual(values = c(up = "#D62728", down = "#1F77B4", ns = "grey70")) +
    labs(
      title = paste("Volcano plot:", comparison_label, timepoint),
      x = "log2 fold change",
      y = "-log10(q-value)",
      colour = "Status"
    ) +
    theme_classic(base_size = 12)

  out_file <- file.path(out_dir, paste0("volcano_", comparison_label, "_", timepoint, ".png"))
  ggsave(out_file, p, width = 7, height = 6, dpi = 300)

  p
}
