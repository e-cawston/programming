required_pkgs <- c("dplyr", "readr", "purrr", "tibble", "ggplot2", "here")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  warning("Missing packages for diagnostics_functions.R: ", paste(missing_pkgs, collapse = ", "))
}

library(dplyr)
library(readr)
library(purrr)
library(tibble)
library(ggplot2)
library(here)

# Is a DEG count distinguishable from chance, or is it what you'd expect from
# pure noise?
#
# Under the null hypothesis that NONE of the tested genes are truly
# differentially expressed, raw p-values are uniformly distributed on
# [0, 1] -- so exactly alpha (e.g. 5%) of them are expected to fall below
# `alpha` purely by chance, regardless of how many genes are tested. This
# tests the observed count of p < alpha against that chance expectation
# with an exact binomial test (one-sided: more hits than chance predicts).
#
# A small p-value here means the comparison has real signal (far more
# low p-values than noise alone would produce). A large p-value means the
# observed hits are consistent with what pure noise would generate -- i.e.
# you cannot distinguish this comparison's DEGs from false discoveries.
test_pvalue_enrichment <- function(p_values, label, alpha = 0.05) {
  p_values <- p_values[!is.na(p_values)]
  n <- length(p_values)
  k <- sum(p_values < alpha)

  test <- binom.test(k, n, p = alpha, alternative = "greater")

  tibble(
    comparison_label = label,
    n_tests = n,
    n_below_alpha = k,
    observed_proportion = k / n,
    expected_proportion_by_chance = alpha,
    fold_enrichment = (k / n) / alpha,
    binom_p_value = test$p.value
  )
}

# Run test_pvalue_enrichment() over several comparisons at once and return
# one combined summary table.
compare_pvalue_enrichment <- function(p_value_list, alpha = 0.05) {
  imap_dfr(p_value_list, function(p_values, label) {
    test_pvalue_enrichment(p_values, label, alpha = alpha)
  }) %>%
    arrange(binom_p_value)
}

# The visual companion to the binomial test: a flat/uniform histogram means
# no real signal (consistent with the null); a spike near zero means real
# differential expression. Facets by timepoint if a `timepoint` column is
# present.
plot_pvalue_histogram <- function(dge_table,
                                  label,
                                  pval_col = "pval",
                                  out_dir = here("results", "figures", "diagnostics")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  df <- dge_table %>% rename(pval = all_of(pval_col)) %>% filter(!is.na(pval))

  p <- ggplot(df, aes(x = pval)) +
    geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "#2a78d6", colour = "white") +
    geom_hline(yintercept = nrow(df) / 20, linetype = "dashed", colour = "#eb6834") +
    labs(
      title = paste("p-value distribution:", label),
      subtitle = "Dashed line = count expected in each bin under the null (uniform p-values)",
      x = "raw p-value",
      y = "count"
    ) +
    theme_minimal(base_size = 12)

  if ("timepoint" %in% names(df)) {
    p <- p + facet_wrap(~timepoint, labeller = label_both)
  }

  out_file <- file.path(out_dir, paste0("pvalue_hist_", label, ".png"))
  ggsave(out_file, p, width = 9, height = 6, dpi = 300)

  invisible(p)
}

# Bar chart summarising test_pvalue_enrichment()/compare_pvalue_enrichment()
# results across comparisons: fold enrichment over the chance rate, with
# the dashed line at 1x marking "no more signal than pure noise" and each
# bar's binomial p-value annotated directly.
plot_pvalue_enrichment_summary <- function(enrichment_tbl,
                                           label,
                                           alpha = 0.05,
                                           out_dir = here("results", "figures", "diagnostics")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot_df <- enrichment_tbl %>%
    mutate(
      comparison_label = factor(comparison_label, levels = comparison_label),
      significant = binom_p_value < alpha,
      p_label = sprintf("%.2f×\np=%.1e", fold_enrichment, binom_p_value)
    )

  p <- ggplot(plot_df, aes(x = comparison_label, y = fold_enrichment, fill = significant)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_text(aes(label = p_label), vjust = -0.15, size = 3.2, lineheight = 0.9) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    scale_fill_manual(
      values = c(`TRUE` = "#2a78d6", `FALSE` = "grey70"),
      labels = c(`TRUE` = "More p<0.05 genes than chance predicts", `FALSE` = "Indistinguishable from chance"),
      drop = FALSE
    ) +
    labs(
      title = "Is the DEG count more than chance?",
      subtitle = "Fold enrichment of genes with p < 0.05 over the 5% expected by chance alone (dashed line = chance level)",
      x = NULL,
      y = "Fold enrichment over chance",
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 20, hjust = 1),
      plot.margin = margin(5.5, 5.5, 5.5, 25)
    )

  out_file <- file.path(out_dir, paste0("pvalue_enrichment_summary_", label, ".png"))
  ggsave(out_file, p, width = 8, height = 6, dpi = 300)

  invisible(p)
}

# Does a query gene set (e.g. the genes that directly differ between two
# mutants) overlap named category gene sets more than chance predicts?
# This is the hypergeometric test -- the same statistic behind GO-term /
# gene-set enrichment analysis, applied here to custom gene sets instead of
# GO terms. See ?phyper; equivalent to a one-sided Fisher's exact test on
# the 2x2 table of (in query / not) x (in category / not).
test_overlap_enrichment <- function(query_set, category_sets, background_n) {
  imap_dfr(category_sets, function(category_set, category_name) {
    overlap <- intersect(query_set, category_set)
    n_query <- length(query_set)
    n_category <- length(category_set)
    expected <- n_query * n_category / background_n

    # P(overlap >= observed) under the hypergeometric null: draw n_query
    # genes at random from background_n, n_category of which are "in the
    # category" -- phyper(k-1, ..., lower.tail=FALSE) gives P(X >= k).
    p_value <- phyper(
      length(overlap) - 1, n_category, background_n - n_category, n_query,
      lower.tail = FALSE
    )

    tibble(
      category = category_name,
      category_size = n_category,
      observed_overlap = length(overlap),
      expected_overlap = expected,
      fold_enrichment = length(overlap) / expected,
      p_value = p_value
    )
  })
}

# Grouped bar chart: observed vs. chance-expected overlap for each
# category, annotated with fold-enrichment and the hypergeometric p-value.
plot_overlap_enrichment <- function(enrichment_tbl,
                                    query_label,
                                    out_dir = here("results", "figures", "diagnostics")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot_df <- enrichment_tbl %>%
    mutate(category = factor(category, levels = category)) %>%
    tidyr::pivot_longer(
      c(observed_overlap, expected_overlap),
      names_to = "kind", values_to = "count"
    ) %>%
    mutate(kind = recode(kind, observed_overlap = "Observed", expected_overlap = "Expected by chance"))

  # Fold/p-value annotation sits above the taller of each pair's two bars,
  # offset by a FIXED amount (relative to the chart's overall max, not
  # each category's own bar height) so short-barred categories don't
  # crowd their annotation text against the per-bar count label just
  # above them -- a multiplicative per-row offset placed the annotation
  # too close whenever that category's bars were small relative to the
  # chart's shared y-axis scale.
  chart_max <- max(c(enrichment_tbl$observed_overlap, enrichment_tbl$expected_overlap), na.rm = TRUE)
  label_df <- enrichment_tbl %>%
    mutate(
      category = factor(category, levels = category),
      label = sprintf("%.1f×\np=%.1e", fold_enrichment, p_value),
      y = pmax(observed_overlap, expected_overlap) + 0.16 * chart_max
    )

  p <- ggplot(plot_df, aes(x = category, y = count, fill = kind)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(
      aes(label = ifelse(kind == "Observed", as.character(count), sprintf("%.1f", count))),
      position = position_dodge(width = 0.7), vjust = -0.4, size = 3
    ) +
    geom_text(data = label_df, aes(x = category, y = y, label = label), inherit.aes = FALSE, size = 3.2, lineheight = 0.9) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.32))) +
    scale_fill_manual(values = c("Observed" = "#2a78d6", "Expected by chance" = "#eb6834")) +
    labs(
      title = paste("Overlap enrichment:", query_label),
      subtitle = "Hypergeometric test: observed vs. chance-expected overlap with each category",
      x = NULL,
      y = "genes overlapping the query set",
      fill = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")

  out_file <- file.path(out_dir, paste0("overlap_enrichment_", query_label, ".png"))
  ggsave(out_file, p, width = 8, height = 6, dpi = 300)

  invisible(p)
}

# Compare several gene lists against the *same* set of category-enrichment
# results (e.g. curated phenotype categories tested for mst7-only,
# pmk1-only, and shared) in one figure: rows = category, columns = gene
# list, dot size = overlap count, dot colour = raw p-value (red = more
# significant, blue = less -- same convention as plot_go_enrichment()),
# with a "*" marking raw p < 0.05. Ordered by the best (smallest) p-value
# any list achieved for that category, so the most interesting rows sit
# at the top -- built for putting several test_overlap_enrichment() runs
# side by side, not tied to GO specifically.
plot_enrichment_comparison <- function(results_list,
                                       min_category_size = 5,
                                       exclude_categories = character(0),
                                       label = "comparison",
                                       title = "Enrichment comparison",
                                       out_dir = here("results", "figures", "diagnostics")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  combined <- imap_dfr(results_list, function(df, list_name) {
    df %>% mutate(query = list_name)
  }) %>%
    filter(category_size >= min_category_size, !category %in% exclude_categories, nzchar(category))

  category_order <- combined %>%
    group_by(category) %>%
    summarise(best_p = min(p_value), .groups = "drop") %>%
    arrange(best_p) %>%
    pull(category)

  plot_df <- combined %>%
    mutate(
      category = factor(category, levels = rev(category_order)),
      query = factor(query, levels = names(results_list)),
      sig_label = if_else(p_value < 0.05, "*", "")
    )

  p <- ggplot(plot_df, aes(x = query, y = category, size = observed_overlap, colour = p_value)) +
    geom_point() +
    geom_text(aes(label = sig_label), colour = "black", size = 6, vjust = 0.75, show.legend = FALSE) +
    scale_colour_gradient(low = "red", high = "blue", name = "p-value") +
    scale_size_continuous(name = "Genes\noverlapping", range = c(2, 10)) +
    labs(
      title = title,
      subtitle = "Dot size = genes overlapping; colour = raw p-value; * = p < 0.05 (uncorrected)",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  out_file <- file.path(out_dir, paste0("enrichment_comparison_", label, ".png"))
  ggsave(out_file, p, width = 7.5, height = max(4, 0.45 * length(category_order) + 1.5), dpi = 300)

  invisible(p)
}
