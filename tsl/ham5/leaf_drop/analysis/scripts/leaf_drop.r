analyse_leafdrop <- function(
  path = "raw",
  annotate = c("none", "anova", "tukey"),
  rename_strains = NULL,
  base_size = 14,
  x_title = "Strain",
  y_title = expression(Lesion~area~(cm^2))
) {

  library(tidyverse)
  library(multcompView)

  annotate <- match.arg(annotate)

  # Find all CSVs
  files <- list.files(path = path, pattern = "^Results_.*\\.csv$", full.names = TRUE)

  # Read and merge
  data_list <- lapply(files, function(f) {
    df <- read.csv(f)
    strain <- basename(f) |>
      sub("^Results_", "", x = _) |>
      sub("\\.csv$", "", x = _)
    df$Strain <- strain
    df
  })

  merged <- bind_rows(data_list)

  # === Rename strains if requested ===
  if (!is.null(rename_strains)) {
    merged$Strain <- dplyr::recode(merged$Strain, !!!rename_strains)
  }

  # ANOVA
  anova_model <- aov(Mean ~ Strain, data = merged)
  anova_summary <- summary(anova_model)

  # Tukey post-hoc
  tukey_res <- TukeyHSD(anova_model)

  # Base plot (no grid, bold text, thicker axes)
  p <- ggplot(merged, aes(x = Strain, y = Area)) +
    geom_boxplot(outlier.shape = NA, fill = "grey85") +
    geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
    theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),                 # remove background grid
      axis.line = element_line(size = 1.2),         # thicker axes
      axis.text.x = element_text(
        angle = 45, hjust = 1,
        size = base_size,
        face = "bold"
      ),
      axis.text.y = element_text(
        size = base_size,
        face = "bold"
      ),
      axis.title.x = element_text(
        size = base_size + 2,
        face = "bold"
      ),
      axis.title.y = element_text(
        size = base_size + 2,
        face = "bold"
      )
    ) +
    ylab(y_title) +
    xlab(x_title)

  # === Add ANOVA p-value annotation ===
  if (annotate == "anova") {
    pval <- anova_summary[[1]]$`Pr(>F)`[1]
    p <- p + annotate(
      "text",
      x = 1,
      y = max(merged$Area),
      label = paste0("ANOVA p = ", signif(pval, 3)),
      hjust = 0,
      vjust = -0.5,
      size = base_size,
      fontface = "bold"
    )
  }

  # === Add Tukey group letters ===
  if (annotate == "tukey") {
    tukey_letters <- multcompView::multcompLetters4(anova_model, tukey_res)
    letters_vec <- tukey_letters$Strain$Letters

    letters_df <- data.frame(
      Strain = names(letters_vec),
      Letters = letters_vec
    )

    p <- p + geom_text(
      data = letters_df,
      aes(x = Strain, y = max(merged$Area) * 1.05, label = Letters),
      size = base_size,
      fontface = "bold"
    )
  }

  # Return everything
  list(
    data = merged,
    plot = p,
    anova = anova_summary,
    tukey = tukey_res
  )
}
