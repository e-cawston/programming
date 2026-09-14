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

perform_anova <- function(filepath,
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

  # Perform one-way ANOVA
  anova_result <- aov(Count ~ Genotype, data = df)
  return(anova_result)
}

perform_tukey <- function(filepath,
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

  # Perform one-way ANOVA and Tukey HSD test
  anova_result <- aov(Count ~ Genotype, data = df)
  tukey_result <- TukeyHSD(anova_result)
  return(tukey_result)
}

plot_with_tukey_annotations <- function(filepath,
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

  # Perform ANOVA and Tukey test
  anova_result <- aov(Count ~ Genotype, data = df)
  tukey_result <- TukeyHSD(anova_result)
  tukey_df <- as.data.frame(tukey_result$Genotype)
  tukey_df$comparison <- rownames(tukey_df)

  # Extract significant comparisons and prepare for plotting
  sig_comparisons <- tukey_df %>%
    filter(`p adj` < 0.05) %>%
    separate(comparison, into = c("group1", "group2"), sep = "-") %>%
    mutate(
      group1 = trimws(group1),
      group2 = trimws(group2),
      stars = ifelse(`p adj` < 0.001, "***",
                     ifelse(`p adj` < 0.01, "**",
                            ifelse(`p adj` < 0.05, "*", "ns"))),
      group1_num = as.numeric(factor(group1, levels = levels(df$Genotype))),
      group2_num = as.numeric(factor(group2, levels = levels(df$Genotype)))
    )

  # Get y-axis limits for positioning annotations
  y_max <- max(df$Count, na.rm = TRUE)
  y_min <- min(df$Count, na.rm = TRUE)
  y_range <- y_max - y_min

  # Create base plot
  p <- ggplot(df, aes(x = Genotype, y = Count)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
      aes(color = Biorep),
      width = 0.2,
      size = 2.5,
      alpha = 0.8
    ) +
    scale_y_continuous(limits = c(0, 50)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    labs(
      title = "Counts per Genotype with Tukey HSD Test",
      x = "Genotype",
      y = "Count",
      color = "Biorep"
    )

  # Add significance bars and stars if there are significant comparisons
  if (nrow(sig_comparisons) > 0) {
    for (i in seq_len(nrow(sig_comparisons))) {
      row <- sig_comparisons[i, ]
      y_pos <- y_max + (i * y_range * 0.05)
      
      # Add horizontal lines for significance bars
      p <- p +
        geom_segment(
          x = row$group1_num,
          xend = row$group2_num,
          y = y_pos,
          yend = y_pos,
          inherit.aes = FALSE,
          color = "black",
          size = 0.5
        ) +
        # Add vertical ticks at the ends
        geom_segment(
          x = row$group1_num,
          xend = row$group1_num,
          y = y_pos - (y_range * 0.01),
          yend = y_pos,
          inherit.aes = FALSE,
          color = "black",
          size = 0.5
        ) +
        geom_segment(
          x = row$group2_num,
          xend = row$group2_num,
          y = y_pos - (y_range * 0.01),
          yend = y_pos,
          inherit.aes = FALSE,
          color = "black",
          size = 0.5
        ) +
        # Add significance stars
        annotate(
          "text",
          x = (row$group1_num + row$group2_num) / 2,
          y = y_pos + (y_range * 0.02),
          label = row$stars,
          size = 5
        )
    }
  }

  return(p)
}

