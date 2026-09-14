required_pkgs <- c(
  "dplyr", "readr", "purrr", "stringr", "tibble", "here", "ggplot2", "RankProd"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  warning("Missing packages for rank_products_functions.R (some features may not work): ", paste(missing_pkgs, collapse = ", "))
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
library(ggplot2)
suppress_library_error('RankProd')

source(file.path(here("analysis", "scripts"), "label_functions.R"))

# If RankProd didn't load, use stub functions
if (!exists('RP')) {
  source(file.path(here('analysis'), 'scripts', 'stub_functions.R'), local = TRUE)
}

build_tpm_matrix_from_samples <- function(s2c_subset) {
  tpm_list <- pmap(s2c_subset %>% select(sample, path), function(sample, path) {
    abundance_file <- if (file.exists(path)) {
      path
    } else {
      file.path(here("raw"), path)
    }

    read_tsv(abundance_file, show_col_types = FALSE) %>%
      select(target_id, tpm) %>%
      rename(!!sample := tpm)
  })

  tpm_mat <- reduce(tpm_list, full_join, by = "target_id")
  rownames_mat <- tpm_mat$target_id

  tpm_mat <- tpm_mat %>%
    select(-target_id) %>%
    as.matrix()

  rownames(tpm_mat) <- rownames_mat
  tpm_mat[, s2c_subset$sample, drop = FALSE]
}

run_rank_products_comparison <- function(s2c,
                                         group1,
                                         group2,
                                         comparison_label,
                                         timepoint = NULL,
                                         out_dir = here("results", "rank_products"),
                                         num_perm = 100) {
  s2c_sub <- s2c %>%
    filter(name %in% c(group1, group2)) %>%
    mutate(name = factor(name, levels = c(group1, group2))) %>%
    arrange(name) %>%
    droplevels()

  if (!is.null(timepoint)) {
    s2c_sub <- s2c_sub %>% filter(timepoint == timepoint)
  }

  if (nrow(s2c_sub) == 0) {
    return(tibble())
  }

  # NOTE: `%>% log2(. + 0.5)` would pass 2 arguments to log2() -- magrittr
  # inserts the LHS as the first argument *in addition to* substituting `.`
  # whenever `.` appears only nested inside another call (here, inside `+`).
  expr_mat <- build_tpm_matrix_from_samples(s2c_sub)
  expr_mat <- log2(expr_mat + 0.5)

  n1 <- sum(s2c_sub$name == group1)
  n2 <- sum(s2c_sub$name == group2)
  cl <- c(rep(1, n1), rep(0, n2))

  rp_result <- RP(
    data = expr_mat,
    cl = cl,
    num.perm = num_perm,
    logged = TRUE,
    gene.names = rownames(expr_mat),
    plot = FALSE
  )

  tg <- topGene(rp_result, cutoff = 1, method = "pfp", logged = TRUE, logbase = 2)

  tidy_rankprod_tbl <- function(tbl, direction, comparison_name, tp) {
    if (is.null(tbl) || nrow(tbl) == 0) {
      return(tibble())
    }

    as.data.frame(tbl) %>%
      tibble::rownames_to_column("target_id") %>%
      mutate(
        direction = direction,
        timepoint = as.character(tp),
        comparison_label = comparison_name,
        group1 = group1,
        group2 = group2
      )
  }

  res_df <- bind_rows(
    tidy_rankprod_tbl(tg$Table1, paste0("up_in_", group1), comparison_label, timepoint),
    tidy_rankprod_tbl(tg$Table2, paste0("up_in_", group2), comparison_label, timepoint)
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0("rankprod_", comparison_label, "_", timepoint, "H.csv"))
  write_csv(res_df, out_file)

  res_df
}

run_rank_products_pair <- function(s2c,
                                   comparison_row,
                                   timepoint = NULL,
                                   out_dir = here("results", "rank_products"),
                                   num_perm = 100) {
  forward <- run_rank_products_comparison(
    s2c = s2c,
    group1 = comparison_row$control,
    group2 = comparison_row$test,
    comparison_label = paste0(comparison_row$control, "_vs_", comparison_row$test),
    timepoint = timepoint,
    out_dir = out_dir,
    num_perm = num_perm
  )

  reverse <- run_rank_products_comparison(
    s2c = s2c,
    group1 = comparison_row$test,
    group2 = comparison_row$control,
    comparison_label = paste0(comparison_row$test, "_vs_", comparison_row$control),
    timepoint = timepoint,
    out_dir = out_dir,
    num_perm = num_perm
  )

  bind_rows(forward, reverse)
}

run_rank_products_comparisons <- function(s2c,
                                          comparisons,
                                          out_dir = here("results", "rank_products"),
                                          num_perm = 100) {
  timepoints <- levels(factor(s2c$timepoint))

  results <- map(seq_len(nrow(comparisons)), function(i) {
    comparison_row <- comparisons[i, , drop = FALSE]

    map(timepoints, function(tp) {
      run_rank_products_pair(
        s2c = s2c,
        comparison_row = comparison_row,
        timepoint = tp,
        out_dir = out_dir,
        num_perm = num_perm
      )
    }) %>%
      bind_rows()
  })

  bind_rows(results)
}

# Unique DEGs from a rank products result table (pools both "up in group1"
# and "up in group2" rows, i.e. both directions of one comparison).
extract_rank_product_deg_ids <- function(rp_table, pfp_threshold = 0.05) {
  rp_table %>%
    filter(pfp <= pfp_threshold) %>%
    pull(target_id) %>%
    unique()
}

# Rank-products equivalent of upset_functions.R::load_deg_tables() -- reads
# the per-timepoint rankprod_<label>_<tp>H.csv files already written by
# run_rank_products_comparisons() and returns one DEG set per comparison
# label, pooled across timepoints. Use this (instead of the sleuth-based
# loader) for comparisons that cross experiments/batches, where sleuth's
# shared-normalisation model isn't appropriate but rank products' per-
# replicate ranking still is.
load_rank_product_deg_tables <- function(comparison_labels,
                                         rp_dir = here("results", "rank_products"),
                                         pfp_threshold = 0.05) {
  deg_list <- map(comparison_labels, function(label) {
    files <- list.files(
      rp_dir,
      pattern = paste0("^rankprod_", label, "_.*H\\.csv$"),
      full.names = TRUE
    )

    if (length(files) == 0) {
      message("No rank products files found for: ", label)
      return(NULL)
    }

    deg_ids <- map(files, function(f) {
      tbl <- read_csv(f, show_col_types = FALSE)
      extract_rank_product_deg_ids(tbl, pfp_threshold)
    })

    unique(unlist(deg_ids, use.names = FALSE))
  })

  names(deg_list) <- comparison_labels
  purrr::compact(deg_list)
}

plot_rank_products_mirrored <- function(rp_table,
                                        comparison_label,
                                        out_dir = here("results", "figures", "rank_products"),
                                        pfp_threshold = 0.05,
                                        display_label = comparison_label) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot_df <- rp_table %>%
    filter(comparison_label == !!comparison_label) %>%
    filter(pfp <= pfp_threshold) %>%
    mutate(direction_label = if_else(direction == paste0("up_in_", group1), "up", "down")) %>%
    count(timepoint, direction_label, name = "n") %>%
    # timepoint_num must be derived after count() -- count() drops every
    # column that isn't a grouping var, so computing it beforehand (and
    # relying on it surviving into this mutate) silently loses the column.
    mutate(
      n_plot = if_else(direction_label == "down", -n, n),
      timepoint_num = suppressWarnings(as.numeric(str_remove(as.character(timepoint), "H"))),
      timepoint = forcats::fct_reorder(as.character(timepoint), timepoint_num, .desc = FALSE)
    )

  if (nrow(plot_df) == 0) {
    return(NULL)
  }

  # Total unique DEGs across all timepoints (a gene hitting threshold at
  # more than one timepoint is only counted once here, unlike the sum of
  # the per-bar counts, which counts per-timepoint events).
  n_unique <- rp_table %>%
    filter(comparison_label == !!comparison_label, pfp <= pfp_threshold) %>%
    pull(target_id) %>%
    n_distinct()

  p <- ggplot(plot_df, aes(x = n_plot, y = timepoint, fill = direction_label)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    geom_text(
      aes(label = n, hjust = if_else(direction_label == "down", 1.15, -0.15)),
      size = 3.2
    ) +
    scale_fill_manual(values = c(up = "#D62728", down = "#1F77B4")) +
    scale_x_continuous(labels = abs, expand = expansion(mult = 0.12)) +
    labs(
      title = strain_nomenclature_title(display_label, prefix = "Rank products: "),
      subtitle = paste0("pfp < ", pfp_threshold),
      x = "Number of DEGs",
      y = "Timepoint",
      fill = NULL,
      caption = paste0("Total unique DEGs across all timepoints: ", n_unique)
    ) +
    theme_classic(base_size = 12) +
    theme(plot.caption = element_text(hjust = 0.5, size = 10, face = "italic"))

  out_file <- file.path(out_dir, paste0("rankprod_mirrored_", comparison_label, ".png"))
  ggsave(out_file, p, width = 8, height = 5, dpi = 300)

  p
}
