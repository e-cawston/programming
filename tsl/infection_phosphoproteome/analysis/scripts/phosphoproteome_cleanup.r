library(tidyverse)

# Load phosphoproteomics TSV produced by MSFragger
load_phospho_data <- function(path) {
  df <- readr::read_tsv(path, guess_max = 10000, show_col_types = FALSE)
  tibble::as_tibble(df)
}

# Helper: detect phosphorylation in a modification string
detect_phospho_in_string <- function(mod_str, mass_shift = 79.9663, tolerance = 0) {
  if (is.na(mod_str) || trimws(mod_str) == "") return(FALSE)
  s <- as.character(mod_str)

  if (stringr::str_detect(s, regex("phospho", ignore_case = TRUE))) return(TRUE)

  nums <- stringr::str_extract_all(s, "-?\\d+\\.?\\d*")[[1]]
  if (length(nums) == 0) return(FALSE)

  nums <- as.numeric(nums)
  any(abs(nums - mass_shift) <= tolerance)
}

resolve_column_name <- function(df, column) {
  column_name <- rlang::eval_tidy(rlang::enquo(column))
  if (length(column_name) != 1 || !is.character(column_name)) {
    stop("column must be a single character column name")
  }
  if (!column_name %in% colnames(df)) {
    stop("Column not found: ", column_name)
  }
  column_name
}

# Filter peptides that contain any phosphorylation (allow other mods)
filter_phospho_peptides_all <- function(df, mod_col = "Assigned Modifications",
                                        mass_shift = 79.9663, tolerance = 0) {

  mod_col_name <- resolve_column_name(df, mod_col)

  df %>%
    dplyr::filter(
      purrr::map_lgl(
        .data[[mod_col_name]],
        ~ detect_phospho_in_string(.x, mass_shift = mass_shift, tolerance = tolerance)
      )
    )
}

# Filter peptides that contain phosphorylation and NO other modifications
filter_phospho_peptides_only <- function(df, mod_col = "Assigned Modifications",
                                         mass_shift = 79.9663, tolerance = 0) {

  mod_col_name <- resolve_column_name(df, mod_col)

  df %>%
    dplyr::filter(
      purrr::map_lgl(
        .data[[mod_col_name]],
        function(mstr) {

          if (is.na(mstr) || trimws(mstr) == "") return(FALSE)

          s <- as.character(mstr)

          has_phospho <- detect_phospho_in_string(
            s,
            mass_shift = mass_shift,
            tolerance = tolerance
          )
          if (!has_phospho) return(FALSE)

          shifts <- stringr::str_extract_all(s, "\\(-?\\d+\\.?\\d*\\)")[[1]]
          shifts <- as.numeric(stringr::str_remove_all(shifts, "[()]"))

          length(shifts) > 0 && all(abs(shifts - mass_shift) <= tolerance)
        }
      )
    )
}

# Keep rows detected in at least min_bioreps for one timepoint.
filter_min_biorep_intensity <- function(df, intensity_cols = NULL, min_bioreps = 2) {
  if (length(min_bioreps) != 1 || !is.numeric(min_bioreps) || min_bioreps < 1) {
    stop("min_bioreps must be a positive number")
  }

  if (is.null(intensity_cols)) {
    intensity_cols <- colnames(df)[
      stringr::str_detect(colnames(df), " Intensity$") &
        !stringr::str_detect(colnames(df), " MaxLFQ Intensity$")
    ]
  } else {
    missing_cols <- setdiff(intensity_cols, colnames(df))
    if (length(missing_cols) > 0) {
      stop("Intensity columns not found: ", paste(missing_cols, collapse = ", "))
    }
  }

  if (length(intensity_cols) == 0) {
    stop("No intensity columns found")
  }

  sample_names <- stringr::str_remove(intensity_cols, " Intensity$")
  sample_parts <- stringr::str_match(sample_names, "^(.+)_([0-9]+)$")
  if (any(is.na(sample_parts[, 3]))) {
    stop("Intensity columns must end with a numeric biorep suffix: ",
         paste(intensity_cols[is.na(sample_parts[, 3])], collapse = ", "))
  }

  timepoint_groups <- split(seq_along(intensity_cols), sample_parts[, 2])
  intensity_values <- vapply(
    df[intensity_cols],
    function(values) as.numeric(values),
    numeric(nrow(df))
  )

  measured_bioreps <- vapply(
    timepoint_groups,
    function(indices) {
      rowSums(intensity_values[, indices, drop = FALSE] > 0, na.rm = TRUE)
    },
    numeric(nrow(df))
  )

  df[apply(measured_bioreps, 1, max) >= min_bioreps, , drop = FALSE]
}

# Keep VSN-normalized peptides observed in at least min_bioreps at
# min_timepoints or more. VSN values may be negative, so missingness is
# determined by finite/non-NA measurements rather than positivity.
filter_vsn_min_timepoints_bioreps <- function(df, intensity_cols = NULL,
                                              min_timepoints = 3,
                                              min_bioreps = 2) {
  validate_minimum <- function(value, name) {
    if (length(value) != 1 || !is.numeric(value) || is.na(value) ||
        value < 1 || value != as.integer(value)) {
      stop(name, " must be a positive integer")
    }
  }
  validate_minimum(min_timepoints, "min_timepoints")
  validate_minimum(min_bioreps, "min_bioreps")

  if (is.null(intensity_cols)) {
    intensity_cols <- colnames(df)[
      stringr::str_detect(
        colnames(df),
        "^.+_((([0-9]+)(_[0-9]+)?h)|spores)_[0-9]+$"
      )
    ]
  } else {
    missing_cols <- setdiff(intensity_cols, colnames(df))
    if (length(missing_cols) > 0) {
      stop("Intensity columns not found: ", paste(missing_cols, collapse = ", "))
    }
  }

  if (length(intensity_cols) == 0) {
    stop("No VSN intensity columns found")
  }

  sample_parts <- stringr::str_match(
    intensity_cols, "^.+_((?:[0-9]+(?:_[0-9]+)?h)|spores)_([0-9]+)$"
  )
  if (anyNA(sample_parts[, 2])) {
    stop("VSN intensity columns must end with a recognized timepoint and biorep: ",
         paste(intensity_cols[is.na(sample_parts[, 2])], collapse = ", "))
  }

  values <- vapply(df[intensity_cols], function(x) as.numeric(x), numeric(nrow(df)))
  timepoint_groups <- split(seq_along(intensity_cols), sample_parts[, 2])
  timepoint_counts <- vapply(timepoint_groups, function(indices) {
    rowSums(is.finite(values[, indices, drop = FALSE]), na.rm = TRUE)
  }, numeric(nrow(df)))

  df[rowSums(timepoint_counts >= min_bioreps) >= min_timepoints, , drop = FALSE]
}

# Full cleanup wrapper
cleanup_phospho_dataset <- function(path, mod_col = "Assigned Modifications",
                                    mass_shift = 79.9663, tolerance = 0,
                                    intensity_cols = NULL, min_bioreps = 2) {

  df <- load_phospho_data(path)

  phospho_all <- filter_phospho_peptides_all(
    df, mod_col = mod_col, mass_shift = mass_shift, tolerance = tolerance
  )

  phospho_only <- filter_phospho_peptides_only(
    df, mod_col = mod_col, mass_shift = mass_shift, tolerance = tolerance
  )

  phospho_all_biorep_filtered <- filter_min_biorep_intensity(
    phospho_all, intensity_cols = intensity_cols, min_bioreps = min_bioreps
  )
  phospho_only_biorep_filtered <- filter_min_biorep_intensity(
    phospho_only, intensity_cols = intensity_cols, min_bioreps = min_bioreps
  )

  list(
    raw = df,
    phospho_all = phospho_all,
    phospho_only = phospho_only,
    phospho_all_biorep_filtered = phospho_all_biorep_filtered,
    phospho_only_biorep_filtered = phospho_only_biorep_filtered
  )
}

