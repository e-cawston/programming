plot_perithecia_boxplots <- function(path, pair_sep = " x ", rename_strains = NULL, font_size = 10) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)

  rename_values <- function(strains, mapping) {
    if (is.null(mapping) || length(mapping) == 0) {
      return(strains)
    }
    out <- strains
    hits <- strains %in% names(mapping)
    out[hits] <- unname(mapping[strains[hits]])
    out
  }

  default_map <- c(
    "1.95" = "ham5 1.95",
    "27.2" = "ham5 27.2",
    "ham5" = "\u0394ham5",
    "mst7" = "\u0394mst7",
    "mst11" = "\u0394mst11",
    "mst50" = "\u0394mst50",
    "pmk1" = "\u0394pmk1"
  )
  rename_map <- default_map
  if (!is.null(rename_strains) && length(rename_strains) > 0) {
    rename_map[names(rename_strains)] <- rename_strains
  }

  # ---- 1. Load CSV ----
  df <- read.csv(path, stringsAsFactors = FALSE)

  # ---- 2. Reshape to long format and parse strain pair ----
  df_long <- df %>%
    pivot_longer(
      cols = c(Top, Bottom),
      names_to = "Position",
      values_to = "Perithecia"
    ) %>%
    mutate(
      RepDate = as.character(Date),
      Rep = paste0("Rep", dense_rank(RepDate)),
      BottomStrain = str_replace(Condition, "^.*_([^_]+)_bottom_([^_]+)_top$", "\\1"),
      TopStrain = str_replace(Condition, "^.*_([^_]+)_bottom_([^_]+)_top$", "\\2"),
      Strain = if_else(Position == "Bottom", BottomStrain, TopStrain),
      StrainLabel = rename_values(Strain, rename_map),
      Pair = paste0(rename_values(BottomStrain, rename_map), pair_sep, rename_values(TopStrain, rename_map)),
      Pair = factor(Pair, levels = unique(Pair)),
      StrainLabel = factor(StrainLabel, levels = unique(StrainLabel)),
      Rep = factor(Rep, levels = unique(Rep))
    )

  # ---- 3. Boxplot grouped by pair with dotplot overlay ----
  cbPalette <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999", "#66C2A5"
  )

  plt <- ggplot(df_long, aes(x = StrainLabel, y = Perithecia)) +
    geom_boxplot(color = "black", fill = "white", outlier.shape = NA, alpha = 0.8) +
    geom_jitter(aes(color = Rep), width = 0.2, height = 0, size = 1.8, alpha = 0.8, show.legend = TRUE) +
    scale_color_manual(values = cbPalette, name = "Rep") +
    facet_wrap(~ Pair, scales = "free_x", nrow = 2) +
    theme_bw(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 8),
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    ) +
    labs(
      x = "Strain",
      y = "Perithecia count"
    )

  return(plt)
}
