## Guy11 between-experiment correlations using non-imputed VSN data filtered
## to sites significant at any timepoint (|log2FC| >= 1, raw p < 0.05) and
## observed in >=2 bioreps at >=3 timepoints in both experiments.

source("analysis/scripts/phosphoproteome_cleanup.r")

required_packages <- c("readxl", "writexl", "dplyr", "tidyr", "ggplot2", "stringr", "purrr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages in the active renv environment: ",
       paste(missing_packages, collapse = ", "))
}

input_paths <- c(
  mst7 = "results/mst7_05_vsn_normalized_nonimputed.xlsx",
  pmk1 = "results/pmk1_05_vsn_normalized_nonimputed.xlsx"
)
timepoint_order <- c("0h", "1h", "1.5h", "2h", "4h", "6h")
timepoint_hours <- c("0h" = 0, "1h" = 1, "1.5h" = 1.5, "2h" = 2, "4h" = 4, "6h" = 6)
identity_cols <- c("Peptide Sequence", "Assigned Modifications")
dir.create("results", recursive = TRUE, showWarnings = FALSE)
figure_root <- file.path("results", "figures", "guy11_vsn_filtered", "significant_sites")
abundance_figure_dir <- file.path(figure_root, "abundance_comparisons")
correlation_figure_dir <- file.path(figure_root, "correlation_distributions")
kinetics_figure_dir <- file.path(figure_root, "kinetics_examples")
for (figure_dir in c(abundance_figure_dir, correlation_figure_dir, kinetics_figure_dir)) {
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
}

diffexp_path <- "results/guy11_timepoint_vs_0h_bootstrap.xlsx"
if (!file.exists(diffexp_path)) {
  stop("Differential-expression results not found: ", diffexp_path,
       ". Run analysis/scripts/run_guy11_bootstrap_diffexp.r first.")
}
diffexp_results <- readxl::read_xlsx(diffexp_path)
required_diffexp_cols <- c("peptide", "log2_fc", "p_value")
missing_diffexp_cols <- setdiff(required_diffexp_cols, colnames(diffexp_results))
if (length(missing_diffexp_cols) > 0) {
  stop("Differential-expression results are missing columns: ",
       paste(missing_diffexp_cols, collapse = ", "))
}
significant_peptides <- diffexp_results %>%
  dplyr::filter(
    is.finite(.data$log2_fc), abs(.data$log2_fc) >= 1,
    is.finite(.data$p_value), .data$p_value < 0.05
  ) %>%
  dplyr::pull(.data$peptide) %>%
  as.character() %>%
  unique()
if (length(significant_peptides) == 0) {
  stop("No sites pass the differential-expression hit rule: |log2FC| >= 1 and p < 0.05.")
}
message("Significant Modified Sequence IDs from differential expression: ",
        length(significant_peptides))

find_guy11_columns <- function(df) {
  columns <- colnames(df)[
    stringr::str_detect(
      colnames(df),
      "^guy11_((([0-9]+)(_[0-9]+)?h)|spores)_[0-9]+$"
    )
  ]
  if (length(columns) == 0) stop("No Guy11 VSN intensity columns found")
  columns
}

peptide_identity <- function(sequence, modifications) {
  valid <- !is.na(sequence) & nzchar(sequence) &
    !is.na(modifications) & nzchar(modifications)
  id <- rep(NA_character_, length(sequence))
  id[valid] <- paste0(
    nchar(sequence[valid]), ":", sequence[valid], "|",
    nchar(modifications[valid]), ":", modifications[valid]
  )
  id
}

read_filtered_experiment <- function(path, experiment) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  filtered_path <- paste0(
    "results/guy11_", experiment, "_vsn_filtered_3tp_2bioreps.xlsx"
  )
  reuse_filtered <- file.exists(filtered_path) &&
    file.info(filtered_path)$mtime >= file.info(path)$mtime
  df <- readxl::read_xlsx(if (reuse_filtered) filtered_path else path)
  missing_identity <- setdiff(identity_cols, colnames(df))
  if (length(missing_identity) > 0) {
    stop("Peptide identity columns missing from ", path, ": ",
         paste(missing_identity, collapse = ", "))
  }
  if (!"Modified Sequence" %in% colnames(df)) {
    stop("Column 'Modified Sequence' is required to match differential-expression hits in ",
         path)
  }

  intensity_cols <- find_guy11_columns(df)
  filtered <- if (reuse_filtered) df else filter_vsn_min_timepoints_bioreps(
    df, intensity_cols = intensity_cols, min_timepoints = 3, min_bioreps = 2
  )
  filtered <- filtered %>%
    dplyr::mutate(.peptide_id = peptide_identity(
      as.character(.data[["Peptide Sequence"]]),
      as.character(.data[["Assigned Modifications"]])
    )) %>%
    dplyr::filter(!is.na(.data$.peptide_id))

  if (anyDuplicated(filtered$.peptide_id)) {
    stop("Peptide sequence + modification identity is not unique in ", path)
  }
  filtered <- filtered %>%
    dplyr::filter(as.character(.data[["Modified Sequence"]]) %in% significant_peptides)
  message(experiment, ": retained ", nrow(filtered),
          " significant peptides meeting the >=3 timepoints x >=2 bioreps rule")
  list(data = filtered, intensity_cols = intensity_cols, reused_filtered = reuse_filtered)
}

experiments <- purrr::imap(input_paths, read_filtered_experiment)
filtered_data <- purrr::map(experiments, "data")

for (experiment in names(filtered_data)) {
  if (!experiments[[experiment]]$reused_filtered) {
    writexl::write_xlsx(
      dplyr::select(filtered_data[[experiment]], -dplyr::all_of(".peptide_id")),
      paste0("results/guy11_", experiment, "_vsn_filtered_3tp_2bioreps.xlsx")
    )
  }
}

shared_ids <- Reduce(intersect, purrr::map(filtered_data, ~ .x$.peptide_id))
if (length(shared_ids) == 0) {
  stop("No differentially significant peptide sequence + modification identities pass the detection filter in both experiments")
}
message("Exact shared peptide sequence + modification identities: ", length(shared_ids))

significant_shared_data <- purrr::map(filtered_data, ~ {
  .x %>%
    dplyr::filter(.data$.peptide_id %in% shared_ids) %>%
    dplyr::select(-dplyr::all_of(".peptide_id"))
})
writexl::write_xlsx(
  c(
    purrr::map(filtered_data, ~ dplyr::select(.x, -dplyr::all_of(".peptide_id"))),
    stats::setNames(significant_shared_data, paste0(names(significant_shared_data), "_shared"))
  ),
  "results/guy11_vsn_filtered_significant_sites.xlsx"
)

to_long <- function(data, experiment) {
  intensity_cols <- experiments[[experiment]]$intensity_cols
  data %>%
    dplyr::filter(.data$.peptide_id %in% shared_ids) %>%
    dplyr::select(dplyr::all_of(c(".peptide_id", identity_cols, intensity_cols))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(intensity_cols),
      names_to = "sample", values_to = "value"
    ) %>%
    dplyr::mutate(
      timepoint_raw = stringr::str_match(
        .data$sample,
        "^guy11_((?:[0-9]+(?:_[0-9]+)?h)|spores)_[0-9]+$"
      )[, 2],
      timepoint = dplyr::case_when(
        .data$timepoint_raw == "spores" ~ "0h",
        TRUE ~ stringr::str_replace(.data$timepoint_raw, "_", ".")
      ),
      time_hours = unname(timepoint_hours[.data$timepoint]),
      bio_rep = stringr::str_extract(.data$sample, "[0-9]+$"),
      experiment = experiment,
      value = as.numeric(.data$value)
    ) %>%
    dplyr::select(
      dplyr::all_of(c(".peptide_id", identity_cols, "experiment",
                      "timepoint", "time_hours", "bio_rep", "value"))
    )
}

guy11_long <- purrr::imap_dfr(filtered_data, to_long) %>%
  dplyr::filter(is.finite(.data$value))
means <- guy11_long %>%
  dplyr::group_by(
    .data$.peptide_id, .data[["Peptide Sequence"]],
    .data[["Assigned Modifications"]], .data$experiment, .data$timepoint
  ) %>%
  dplyr::summarise(mean_intensity = mean(.data$value), .groups = "drop")

writexl::write_xlsx(
  dplyr::rename(
    means, peptide_sequence = `Peptide Sequence`,
    modifications = `Assigned Modifications`
  ),
    "results/guy11_vsn_filtered_significant_shared_timepoint_means.xlsx"
)

paired_means <- means %>%
  dplyr::mutate(timepoint = factor(.data$timepoint, levels = timepoint_order)) %>%
  tidyr::pivot_wider(
    id_cols = dplyr::all_of(c(".peptide_id", identity_cols, "timepoint")),
    names_from = "experiment", values_from = "mean_intensity"
  )

per_peptide_cor <- paired_means %>%
  dplyr::group_by(
    .data$.peptide_id, .data[["Peptide Sequence"]],
    .data[["Assigned Modifications"]]
  ) %>%
  dplyr::group_modify(function(.x, .y) {
    complete <- stats::complete.cases(.x[, c("mst7", "pmk1")])
    x <- .x$mst7[complete]
    y <- .x$pmk1[complete]
    n_points <- length(x)
    if (n_points < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) {
      return(tibble::tibble(
        n_timepoints = n_points,
        r = NA_real_, p_value = NA_real_,
        spearman_r = NA_real_, spearman_p_value = NA_real_
      ))
    }
    pearson_test <- stats::cor.test(x, y, method = "pearson")
    spearman_test <- stats::cor.test(x, y, method = "spearman", exact = FALSE)
    tibble::tibble(
      n_timepoints = n_points,
      r = unname(pearson_test$estimate),
      p_value = pearson_test$p.value,
      spearman_r = unname(spearman_test$estimate),
      spearman_p_value = spearman_test$p.value
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    peptide_sequence = `Peptide Sequence`,
    modifications = `Assigned Modifications`
  )
writexl::write_xlsx(
  per_peptide_cor,
  "results/guy11_vsn_filtered_significant_peptide_kinetics_correlation.xlsx"
)

model_predictions <- guy11_long %>%
  dplyr::group_by(
    .data$.peptide_id, .data[["Peptide Sequence"]],
    .data[["Assigned Modifications"]]
  ) %>%
  dplyr::group_modify(function(.x, .y) {
    available_mst7 <- unique(as.character(.x$timepoint[.x$experiment == "mst7"]))
    available_pmk1 <- unique(as.character(.x$timepoint[.x$experiment == "pmk1"]))
    common_timepoints <- timepoint_order[
      timepoint_order %in% available_mst7 & timepoint_order %in% available_pmk1
    ]
    if (length(common_timepoints) < 3) {
      return(tibble::tibble(
        experiment = character(), timepoint = character(),
        predicted_intensity = numeric(), n_observations = integer()
      ))
    }

    model_data <- .x %>%
      dplyr::mutate(experiment = factor(.data$experiment, levels = c("mst7", "pmk1")))
    fit <- stats::lm(
      value ~ experiment * splines::ns(time_hours, df = 2),
      data = model_data
    )
    prediction_grid <- tidyr::expand_grid(
      experiment = factor(c("mst7", "pmk1"), levels = c("mst7", "pmk1")),
      timepoint = common_timepoints
    ) %>%
      dplyr::mutate(
        time_hours = unname(timepoint_hours[.data$timepoint])
      )
    prediction_grid$predicted_intensity <- as.numeric(
      stats::predict(fit, newdata = prediction_grid)
    )
    prediction_grid %>%
      dplyr::mutate(
        n_observations = nrow(model_data),
        experiment = as.character(.data$experiment)
      )
  }) %>%
  dplyr::ungroup()

model_predicted_cor <- model_predictions %>%
  dplyr::group_by(
    .data$.peptide_id, .data[["Peptide Sequence"]],
    .data[["Assigned Modifications"]]
  ) %>%
  dplyr::group_modify(function(.x, .y) {
    paired <- .x %>%
      dplyr::select(dplyr::all_of(c("timepoint", "experiment", "predicted_intensity"))) %>%
      tidyr::pivot_wider(names_from = "experiment", values_from = "predicted_intensity")
    complete <- stats::complete.cases(paired[, c("mst7", "pmk1")])
    x <- paired$mst7[complete]
    y <- paired$pmk1[complete]
    n_points <- length(x)
    if (n_points < 3 || stats::sd(x) == 0 || stats::sd(y) == 0) {
      return(tibble::tibble(
        n_timepoints = n_points,
        r = NA_real_, p_value = NA_real_,
        spearman_r = NA_real_, spearman_p_value = NA_real_
      ))
    }
    pearson_test <- stats::cor.test(x, y, method = "pearson")
    spearman_test <- stats::cor.test(x, y, method = "spearman", exact = FALSE)
    tibble::tibble(
      n_timepoints = n_points,
      r = unname(pearson_test$estimate),
      p_value = pearson_test$p.value,
      spearman_r = unname(spearman_test$estimate),
      spearman_p_value = spearman_test$p.value
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    peptide_sequence = `Peptide Sequence`,
    modifications = `Assigned Modifications`
  )
writexl::write_xlsx(
  model_predicted_cor,
  "results/guy11_vsn_filtered_significant_model_predicted_kinetics_correlation.xlsx"
)
writexl::write_xlsx(
  dplyr::rename(
    model_predictions,
    peptide_sequence = `Peptide Sequence`,
    modifications = `Assigned Modifications`
  ),
  "results/guy11_vsn_filtered_significant_model_predicted_timepoint_means.xlsx"
)

make_example_sets <- function(data, correlation_column) {
  list(
    negative = data %>% dplyr::filter(.data[[correlation_column]] <= -0.5),
    near_zero = data %>% dplyr::filter(abs(.data[[correlation_column]]) <= 0.1),
    positive = data %>% dplyr::filter(.data[[correlation_column]] >= 0.5)
  )
}
example_sets <- list(
  pearson = make_example_sets(per_peptide_cor, "r"),
  spearman = make_example_sets(per_peptide_cor, "spearman_r")
)
set.seed(1)
example_sets <- purrr::map(example_sets, ~ purrr::map(.x, function(bucket) {
  eligible <- dplyr::filter(bucket, is.finite(.data$r), is.finite(.data$spearman_r))
  dplyr::slice_sample(eligible, n = min(10L, nrow(eligible)))
})
)

mean_long <- means %>%
  dplyr::mutate(timepoint = factor(.data$timepoint, levels = timepoint_order))
plot_kinetics_examples <- function(selection, method, bucket) {
  labels <- selection %>%
    dplyr::mutate(
      label = sprintf(
        "%s\n%s\nPearson r = %.2f; Spearman rho = %.2f",
        .data$peptide_sequence, .data$modifications,
        .data$r, .data$spearman_r
      )
    ) %>%
    dplyr::select(dplyr::all_of(c(".peptide_id", "label")))
  mean_long %>%
    dplyr::inner_join(labels, by = ".peptide_id") %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = .data$timepoint, y = .data$mean_intensity,
        color = .data$experiment, group = .data$experiment
      )
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~label, ncol = 5, scales = "free_y") +
    ggplot2::scale_color_manual(values = c(mst7 = "#0072B2", pmk1 = "#D55E00")) +
    ggplot2::labs(
      title = paste0(toupper(method), " ", bucket, "-correlation site trajectories"),
      x = "Timepoint", y = "Mean VSN-normalized intensity", color = "Experiment"
    ) +
    ggplot2::theme_classic(base_family = "sans", base_size = 10) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold", size = 11),
      axis.text = ggplot2::element_text(size = 9),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(size = 7)
    )
}

correlation_at_pair <- function(data, time_a, time_b) {
  values_a <- data %>%
    dplyr::filter(.data$experiment == "mst7", .data$timepoint == time_a) %>%
    dplyr::select(dplyr::all_of(c(".peptide_id", "mean_intensity"))) %>%
    dplyr::rename(x = mean_intensity)
  values_b <- data %>%
    dplyr::filter(.data$experiment == "pmk1", .data$timepoint == time_b) %>%
    dplyr::select(dplyr::all_of(c(".peptide_id", "mean_intensity"))) %>%
    dplyr::rename(y = mean_intensity)
  paired <- dplyr::inner_join(values_a, values_b, by = ".peptide_id")
  paired <- dplyr::filter(paired, is.finite(.data$x), is.finite(.data$y))
  n_peptides <- nrow(paired)
  r <- if (n_peptides >= 3 && stats::sd(paired$x) > 0 && stats::sd(paired$y) > 0) {
    stats::cor(paired$x, paired$y, method = "pearson")
  } else {
    NA_real_
  }
  r2 <- if (is.finite(r)) r^2 else NA_real_
  list(
    points = paired,
    summary = tibble::tibble(
      mst7_timepoint = time_a, pmk1_timepoint = time_b,
      n_peptides = n_peptides, pearson_r = r, r_squared = r2
    )
  )
}

timepoint_pairs <- tibble::tibble(
  mst7_timepoint = factor(timepoint_order, levels = timepoint_order),
  pmk1_timepoint = factor(timepoint_order, levels = timepoint_order)
)
pair_results <- purrr::pmap(
  timepoint_pairs,
  ~ correlation_at_pair(mean_long, ..1, ..2)
)
timepoint_correlations <- dplyr::bind_rows(purrr::map(pair_results, "summary")) %>%
  dplyr::mutate(
    mst7_timepoint = factor(.data$mst7_timepoint, levels = timepoint_order),
    pmk1_timepoint = factor(.data$pmk1_timepoint, levels = timepoint_order)
  )
if (!all(as.character(timepoint_correlations$mst7_timepoint) ==
         as.character(timepoint_correlations$pmk1_timepoint))) {
  stop("Abundance comparisons must use matching timepoints only")
}
writexl::write_xlsx(
  timepoint_correlations,
  "results/guy11_vsn_filtered_significant_timepoint_pair_correlations.xlsx"
)

save_figure <- function(plot, stem, width, height, directory = figure_root) {
  ggplot2::ggsave(
    file.path(directory, paste0(stem, ".png")), plot, width = width, height = height,
    units = "in", dpi = 300, bg = "white"
  )
  ggplot2::ggsave(
    file.path(directory, paste0(stem, ".pdf")), plot, width = width, height = height,
    units = "in", bg = "white"
  )
}

make_correlation_distribution <- function(data, column, title, x_label) {
  values <- data %>%
    dplyr::filter(is.finite(.data[[column]]))
  ggplot2::ggplot(values, ggplot2::aes(x = .data[[column]])) +
    ggplot2::geom_histogram(
      binwidth = 0.05, boundary = -1, color = "white", fill = "#2F6690"
    ) +
    ggplot2::geom_vline(
      xintercept = stats::median(values[[column]]),
      linetype = "dashed", color = "#D55E00", linewidth = 0.8
    ) +
    ggplot2::scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
    ggplot2::labs(
      x = x_label, y = "Number of phosphosites", title = title,
      subtitle = sprintf("n = %d; median = %.2f", nrow(values), stats::median(values[[column]]))
    ) +
    ggplot2::theme_classic(base_family = "sans", base_size = 10)
}

pearson_distribution <- make_correlation_distribution(
  per_peptide_cor, "r", "Per-site Pearson kinetics correlations (significant sites)",
  "Pearson correlation (r)"
)
spearman_distribution <- make_correlation_distribution(
  per_peptide_cor, "spearman_r", "Per-site Spearman kinetics correlations (significant sites)",
  "Spearman correlation (rho)"
)
model_distribution <- make_correlation_distribution(
  model_predicted_cor, "r", "Model-predicted Pearson kinetics correlations (significant sites)",
  "Pearson correlation (r)"
)
save_figure(
  pearson_distribution, "guy11_vsn_filtered_significant_pearson_correlation_distribution",
  9, 6, correlation_figure_dir
)
save_figure(
  spearman_distribution, "guy11_vsn_filtered_significant_spearman_correlation_distribution",
  9, 6, correlation_figure_dir
)
save_figure(
  model_distribution, "guy11_vsn_filtered_significant_model_predicted_correlation_distribution",
  9, 6, correlation_figure_dir
)

for (i in seq_len(nrow(timepoint_pairs))) {
  pair <- timepoint_pairs[i, ]
  plot_data <- pair_results[[i]]$points
  r_squared <- timepoint_correlations$r_squared[[i]]
  r_squared_label <- if (is.finite(r_squared)) {
    sprintf("R-squared = %.3f", r_squared)
  } else {
    "R-squared unavailable"
  }
  timepoint_label <- as.character(pair$mst7_timepoint)
  comparison_plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point(color = "#0072B2", alpha = 0.3, size = 0.8) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                         color = "#D55E00", linewidth = 0.8) +
    ggplot2::labs(
      title = paste0("Significant Guy11 sites: ", timepoint_label, " vs ", timepoint_label),
      subtitle = paste0("n = ", nrow(plot_data), " shared sites; ", r_squared_label),
      x = paste0("MST7 mean VSN-normalized intensity at ", timepoint_label),
      y = paste0("PMK1 mean VSN-normalized intensity at ", timepoint_label)
    ) +
    ggplot2::theme_classic(base_family = "sans", base_size = 10)
  safe_timepoint <- gsub("\\.", "_", timepoint_label)
  save_figure(
    comparison_plot,
    paste0("guy11_vsn_filtered_significant_abundance_", safe_timepoint, "_mst7_vs_pmk1"),
    7, 6, abundance_figure_dir
  )
}

for (method in names(example_sets)) {
  for (bucket in names(example_sets[[method]])) {
    examples <- example_sets[[method]][[bucket]]
    file_stem <- paste0(
      "guy11_vsn_filtered_significant_", method, "_kinetics_examples_", bucket
    )
    writexl::write_xlsx(examples, paste0("results/", file_stem, ".xlsx"))
    if (nrow(examples) > 0) {
      save_figure(
        plot_kinetics_examples(examples, method, bucket),
        file_stem,
        14, max(5, 2 + ceiling(nrow(examples) / 5) * 2.6),
        kinetics_figure_dir
      )
    }
  }
}

writeLines(
  c(
    "Significant-site Guy11 filtered VSN kinetics example figures: each facet is labeled with peptide sequence and exact assigned modifications; curves show per-timepoint means in the MST7 and PMK1 experiments, and each title identifies the Pearson- or Spearman-selected correlation bucket.",
    "Significant-site Guy11 filtered VSN abundance scatterplots: each figure compares the same timepoint between MST7 and PMK1; each point is one shared peptide sequence + modification identity; the line is a linear fit and the subtitle reports R-squared.",
    "Significant-site Guy11 filtered VSN correlation distribution figures summarize per-site Pearson, Spearman, and model-predicted Pearson trajectory correlations. Sites pass |log2FC| >= 1 and raw p < 0.05 at any timepoint in either experiment."
  ),
  "results/guy11_vsn_filtered_significant_figure_alt_text.txt"
)

message("Wrote significant-site filtered shared-peptide correlations, same-timepoint comparisons, example tables, and figures to results/.")
