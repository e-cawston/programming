required_pkgs <- c(
  "dplyr", "readr", "purrr", "stringr", "tibble", "here", "UpSetR"
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  warning("Missing packages for upset_functions.R: ", paste(missing_pkgs, collapse = ", "))
}

suppress_warnings <- function(expr) {
  tryCatch(expr, error = function(e) {
    warning("Failed to load package: ", e$message)
  })
}

library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tibble)
library(here)
suppress_warnings(library(UpSetR))

extract_deg_ids <- function(dge_table,
                            qval_threshold = 0.05,
                            fc_threshold = 0) {
  dge_table %>%
    filter(qval <= qval_threshold) %>%
    filter(abs(b) >= fc_threshold) %>%
    pull(target_id) %>%
    unique()
}

load_deg_tables <- function(comparison_labels,
                            dge_dir = here("results", "dge"),
                            qval_threshold = 0.05,
                            fc_threshold = 0) {
  deg_list <- map(comparison_labels, function(label) {
    files <- list.files(
      dge_dir,
      pattern = paste0("^sleuth_", label, "_.*H\\.csv$"),
      full.names = TRUE
    )

    if (length(files) == 0) {
      message("No DGE files found for: ", label)
      return(NULL)
    }

    deg_ids <- map(files, function(f) {
      tbl <- read_csv(f, show_col_types = FALSE)
      extract_deg_ids(tbl, qval_threshold, fc_threshold)
    })

    unique(unlist(deg_ids, use.names = FALSE))
  })

  names(deg_list) <- comparison_labels
  purrr::compact(deg_list)
}

plot_upset <- function(deg_list,
                       label = "deg_overlap",
                       out_dir = here("results", "figures", "upset")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (length(deg_list) < 2) {
    stop("At least two DEG sets are required for an upset plot.")
  }

  all_genes <- unique(unlist(deg_list, use.names = FALSE))

  mat <- map_dfc(names(deg_list), function(nm) {
    tibble(!!nm := all_genes %in% deg_list[[nm]])
  })

  mat <- bind_cols(tibble(gene = all_genes), mat)
  df <- as.data.frame(mat)

  png(file.path(out_dir, paste0("upset_", label, ".png")),
      width = 1800, height = 1200, res = 180)

  upset(
    df,
    nsets = length(deg_list),
    sets = names(deg_list),
    order.by = "freq",
    keep.order = TRUE
  )

  dev.off()

  invisible(df)
}
