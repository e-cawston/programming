## VSN normalisation and MAR/MNAR imputation helpers

find_intensity_columns <- function(df, suffix = " Intensity") {
  columns <- colnames(df)[stringr::str_ends(colnames(df), suffix)]
  columns <- columns[!stringr::str_ends(columns, " MaxLFQ Intensity")]
  if (length(columns) == 0) {
    stop("No ordinary Intensity columns found")
  }
  columns
}

make_site_ids <- function(df, id_col = "Unique_site_ID") {
  if (id_col %in% colnames(df)) {
    ids <- as.character(df[[id_col]])
    if (all(!is.na(ids) & nzchar(ids)) && !anyDuplicated(ids)) {
      return(ids)
    }
    warning("Column '", id_col, "' is missing, empty, or duplicated; using row IDs")
  }
  paste0("site_", seq_len(nrow(df)))
}

classify_missingness <- function(missing, sample_names) {
  groups <- stringr::str_remove(sample_names, "_[0-9]+$")
  group_names <- unique(groups)
  mnar <- matrix(FALSE, nrow(missing), ncol(missing), dimnames = dimnames(missing))
  mar <- mnar

  for (group in group_names) {
    columns <- which(groups == group)
    observed <- rowSums(!missing[, columns, drop = FALSE])
    mnar[, columns] <- mnar[, columns] | missing[, columns, drop = FALSE] & observed <= 1
    mar[, columns] <- mar[, columns] | missing[, columns, drop = FALSE] & observed == 2
  }

  list(mnar = mnar, mar = mar, groups = groups)
}

# Process one saved cleanup table independently and save normalized/imputed tables.
process_vsn_imputation <- function(input_path, experiment_name, output_dir = "results",
                                   id_col = "Unique_site_ID", intensity_suffix = " Intensity",
                                   k = 10, rowmax = 0.9, colmax = 0.9,
                                   quantile_prob = 0.01, stage = "05") {
  required_packages <- c("readxl", "stringr", "vsn", "impute", "writexl")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages) > 0) {
    stop("Missing required packages in the active renv environment: ",
         paste(missing_packages, collapse = ", "))
  }

  if (length(input_path) != 1 || !file.exists(input_path)) {
    stop("input_path must point to an existing Excel file")
  }
  if (length(experiment_name) != 1 || !nzchar(experiment_name)) {
    stop("experiment_name must be a non-empty string")
  }
  if (length(stage) != 1 || !grepl("^[0-9]{2}$", stage)) {
    stop("stage must be a two-digit string, such as '05'")
  }

  df <- readxl::read_xlsx(input_path)
  if (nrow(df) == 0) {
    stop("The input Excel file contains no rows: ", input_path)
  }

  intensity_cols <- find_intensity_columns(df, suffix = intensity_suffix)
  sample_names <- stringr::str_remove(intensity_cols, stringr::fixed(intensity_suffix))
  sample_parts <- stringr::str_match(sample_names, "^(.+)_([0-9]+)$")
  if (anyNA(sample_parts[, 3])) {
    stop("Intensity sample names must end with '_rep': ",
         paste(sample_names[is.na(sample_parts[, 3])], collapse = ", "))
  }
  group_counts <- table(stringr::str_remove(sample_names, "_[0-9]+$"))
  if (any(group_counts != 3)) {
    stop("Each genotype/timepoint must have exactly 3 replicates: ",
         paste(names(group_counts)[group_counts != 3], collapse = ", "))
  }

  site_ids <- make_site_ids(df, id_col = id_col)
  raw_mat <- vapply(df[intensity_cols], as.numeric, numeric(nrow(df)))
  dimnames(raw_mat) <- list(site_ids, sample_names)
  raw_mat[!is.finite(raw_mat) | raw_mat <= 0] <- NA_real_

  vsn_mat <- vsn::justvsn(raw_mat)
  dimnames(vsn_mat) <- dimnames(raw_mat)
  missingness <- classify_missingness(is.na(vsn_mat), sample_names)

  imputed_mat <- vsn_mat
  quantiles <- apply(vsn_mat, 2, function(values) {
    finite_values <- values[is.finite(values)]
    if (length(finite_values) == 0) NA_real_ else
      as.numeric(stats::quantile(finite_values, probs = quantile_prob, names = FALSE))
  })
  if (any(missingness$mnar & is.na(matrix(quantiles,
                                          nrow = nrow(vsn_mat),
                                          ncol = ncol(vsn_mat), byrow = TRUE)))) {
    stop("Cannot impute MNAR values in a sample with no finite VSN values")
  }
  for (column in seq_len(ncol(imputed_mat))) {
    imputed_mat[missingness$mnar[, column], column] <- quantiles[column]
  }

  mar_missing <- missingness$mar & is.na(imputed_mat)
  if (any(mar_missing)) {
    knn_result <- impute::impute.knn(imputed_mat, k = k,
                                     rowmax = rowmax, colmax = colmax)
    imputed_mat[mar_missing] <- knn_result$data[mar_missing]
  }

  metadata <- df[, setdiff(colnames(df), intensity_cols), drop = FALSE]
  metadata[[id_col]] <- site_ids
  metadata <- metadata[, c(id_col, setdiff(colnames(metadata), id_col)), drop = FALSE]
  normalized_df <- cbind(metadata, as.data.frame(vsn_mat, check.names = FALSE))
  imputed_df <- cbind(metadata, as.data.frame(imputed_mat, check.names = FALSE))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  normalized_path <- file.path(output_dir,
             paste0(experiment_name, "_", stage,
               "_vsn_normalized_nonimputed.xlsx"))
  imputed_path <- file.path(output_dir,
          paste0(experiment_name, "_", stage,
            "_vsn_mar_mnar_imputed.xlsx"))
  summary_path <- file.path(output_dir,
          paste0(experiment_name, "_", stage,
            "_vsn_imputation_summary.xlsx"))
  writexl::write_xlsx(normalized_df, normalized_path)
  writexl::write_xlsx(imputed_df, imputed_path)

  summary <- tibble::tibble(
    experiment = experiment_name,
    sites = nrow(vsn_mat),
    samples = ncol(vsn_mat),
    missing_before = sum(is.na(vsn_mat)),
    mnar_imputed = sum(missingness$mnar),
    mar_imputed = sum(missingness$mar),
    total_imputed = sum(missingness$mnar) + sum(missingness$mar),
    remaining_na = sum(is.na(imputed_mat))
  )
  writexl::write_xlsx(summary, summary_path)
  message(experiment_name, ": ", nrow(vsn_mat), " sites x ", ncol(vsn_mat),
          " samples; MNAR=", summary$mnar_imputed,
          ", MAR=", summary$mar_imputed,
          ", remaining NA=", summary$remaining_na)

  list(
    raw_matrix = raw_mat,
    vsn_matrix = vsn_mat,
    imputed_matrix = imputed_mat,
    normalized = normalized_df,
    imputed = imputed_df,
    missing_mnar = missingness$mnar,
    missing_mar = missingness$mar,
    summary = summary,
    normalized_path = normalized_path,
    imputed_path = imputed_path,
    summary_path = summary_path
  )
}