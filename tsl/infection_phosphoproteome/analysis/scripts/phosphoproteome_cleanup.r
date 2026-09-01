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


