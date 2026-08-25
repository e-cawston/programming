plot_genotype_boxplot <- function(filepath,
                                  genotype_order = NULL,
                                  reps = NULL) {
  library(tidyverse)
  library(readxl)

  df <- read_excel(filepath)

  # Filter by reps if provided
  if (!is.null(reps)) {
    df <- df %>% filter(Biorep %in% reps)
  }

  # Apply genotype order if provided
  if (!is.null(genotype_order)) {
    df <- df %>%
      mutate(
        Genotype = factor(Genotype, levels = genotype_order),
        Biorep   = factor(Biorep)
      )
  } else {
    df <- df %>%
      mutate(
        Genotype = factor(Genotype),
        Biorep   = factor(Biorep)
      )
  }

  ggplot(df, aes(x = Genotype, y = Count)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
      aes(color = Biorep),
      width = 0.2,
      size = 2.5,
      alpha = 0.8
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    labs(
      title = "Counts per Genotype",
      x = "Genotype",
      y = "Count",
      color = "Biorep"
    )
}

