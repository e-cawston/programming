read_appressoria_counts <- function(filepath) {
  required_columns <- c("Strain", "Appressoria", "Germtube_only")
  data <- read.csv(filepath, check.names = FALSE)

  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      "The input file is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(data[required_columns])) {
    stop("Required count columns contain missing values.", call. = FALSE)
  }

  if (any(data$Appressoria < 0) || any(data$Germtube_only < 0)) {
    stop("Appressoria counts cannot be negative.", call. = FALSE)
  }

  data
}

chi_square_appressoria <- function(filepath,
                                    reference_strain = "Guy11",
                                    p_adjust_method = "BH") {
  data <- read_appressoria_counts(filepath)

  if (!reference_strain %in% data$Strain) {
    stop("Reference strain not found: ", reference_strain, call. = FALSE)
  }

  strains <- unique(data$Strain)
  comparison_strains <- setdiff(strains, reference_strain)
  reference_counts <- colSums(
    data[data$Strain == reference_strain, c("Appressoria", "Germtube_only")],
    na.rm = TRUE
  )

  results <- lapply(comparison_strains, function(strain) {
    strain_counts <- colSums(
      data[data$Strain == strain, c("Appressoria", "Germtube_only")],
      na.rm = TRUE
    )
    counts <- rbind(reference_counts, strain_counts)
    test <- suppressWarnings(chisq.test(counts, correct = FALSE))

    data.frame(
      reference = reference_strain,
      strain = strain,
      reference_appressoria = unname(reference_counts["Appressoria"]),
      reference_germtube_only = unname(reference_counts["Germtube_only"]),
      strain_appressoria = unname(strain_counts["Appressoria"]),
      strain_germtube_only = unname(strain_counts["Germtube_only"]),
      statistic = unname(test$statistic),
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, results)
  results$p_adjusted <- p.adjust(results$p_value, method = p_adjust_method)
  results$significance <- ifelse(
    results$p_adjusted < 0.001, "***",
    ifelse(results$p_adjusted < 0.01, "**",
           ifelse(results$p_adjusted < 0.05, "*", "ns"))
  )
  results$significant <- results$p_adjusted < 0.05
  results
}

plot_appressoria_boxplot <- function(filepath,
                                      strain_order = NULL,
                                      annotate_significance = TRUE,
                                      reference_strain = "Guy11") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required to make the boxplot.", call. = FALSE)
  }

  data <- read_appressoria_counts(filepath)
  data$Percentage_appressoria <- with(
    data,
    100 * Appressoria / (Appressoria + Germtube_only)
  )

  if (is.null(strain_order)) {
    strain_order <- unique(data$Strain)
  }
  data$Strain <- factor(data$Strain, levels = strain_order)

  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(x = Percentage_appressoria, y = Strain)
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      fill = "gray97",
      color = "black"
    ) +
    ggplot2::geom_jitter(
      height = 0.15,
      width = 0,
      size = 1.8,
      alpha = 0.4,
      color = "darkblue"
    ) +
    ggplot2::labs(
      title = "Appressoria Formation Phenotype",
      x = "Germtubes forming appressoria (%)",
      y = "Genotype"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 110),
      expand = ggplot2::expansion(mult = c(0.02, 0.12))
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"))

  if (annotate_significance) {
    results <- chi_square_appressoria(
      filepath,
      reference_strain = reference_strain
    )
    labels <- data.frame(
      Strain = factor(results$strain, levels = strain_order),
      label = results$significance,
      x = 106,
      stringsAsFactors = FALSE
    )
    plot <- plot +
      ggplot2::geom_text(
        data = labels,
        ggplot2::aes(x = x, y = Strain, label = label),
        inherit.aes = FALSE,
        color = "red",
        fontface = "bold",
        size = 5
      )
  }

  plot
}

perform_appressoria_anova <- function(filepath) {
  data <- read_appressoria_counts(filepath)
  data$Percentage_appressoria <- with(
    data,
    100 * Appressoria / (Appressoria + Germtube_only)
  )
  stats::aov(Percentage_appressoria ~ Strain, data = data)
}

perform_appressoria_tukey <- function(filepath) {
  stats::TukeyHSD(perform_appressoria_anova(filepath))
}

plot_appressoria_boxplot_no_signif <- function(filepath,
                                               strain_order = NULL) {
  plot_appressoria_boxplot(
    filepath,
    strain_order = strain_order,
    annotate_significance = FALSE
  )
}
