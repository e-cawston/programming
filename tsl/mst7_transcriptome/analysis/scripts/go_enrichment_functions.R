required_pkgs <- c("dplyr", "readr", "purrr", "tibble", "ggplot2", "here")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_pkgs) > 0) {
  warning("Missing packages for go_enrichment_functions.R: ", paste(missing_pkgs, collapse = ", "))
}

library(dplyr)
library(readr)
library(purrr)
library(tibble)
library(ggplot2)
library(here)

# Gene -> GO term associations (term2gene) and GO term id -> name/domain
# (term2name), built once in raw/go_annotation/ from UniProt + the GO
# consortium (see raw/go_annotation/README.md for provenance) since no
# annotation file was available in this project.
load_go_annotation <- function(term2gene_file = here("raw", "go_annotation", "term2gene.csv"),
                               term2name_file = here("raw", "go_annotation", "term2name.csv")) {
  list(
    term2gene = read_csv(term2gene_file, show_col_types = FALSE),
    term2name = read_csv(term2name_file, show_col_types = FALSE)
  )
}

# Over-representation analysis (ORA): for every GO term with at least
# `min_term_size` annotated genes in the universe, test whether the query
# gene set contains more of that term's genes than chance predicts, via
# the hypergeometric test (identical statistic to a one-sided Fisher's
# exact test on the term-membership x query-membership 2x2 table -- the
# same test used for GO enrichment in tools like clusterProfiler's
# enricher(), and the same one used earlier this session for the
# mst7/pmk1 gene-set overlap analysis, just applied per GO term here).
# p-values are BH-adjusted across every term actually tested.
run_go_enrichment <- function(gene_ids,
                              universe,
                              term2gene,
                              term2name = NULL,
                              min_term_size = 3,
                              qval_threshold = 0.05) {
  gene_ids <- intersect(gene_ids, universe)
  n_universe <- length(universe)
  n_query <- length(gene_ids)

  # The three GO root terms ("any biological_process/molecular_function/
  # cellular_component") are trivially carried by every annotated gene --
  # testing them produces a spurious "significant" hit reflecting nothing
  # more than "this list's genes tend to be annotated at all", not a real
  # functional signal. Always exclude them.
  go_roots <- c("GO:0008150", "GO:0003674", "GO:0005575")
  t2g <- term2gene %>% filter(target_id %in% universe, !go_id %in% go_roots)

  term_sizes <- t2g %>% count(go_id, name = "term_size") %>% filter(term_size >= min_term_size)
  t2g <- t2g %>% semi_join(term_sizes, by = "go_id")

  tested <- t2g %>%
    group_by(go_id) %>%
    summarise(
      term_size = n(),
      overlap_genes = paste(intersect(target_id, gene_ids), collapse = "/"),
      overlap = length(intersect(target_id, gene_ids)),
      .groups = "drop"
    )

  result <- tested %>%
    filter(overlap > 0) %>%
    mutate(
      expected = n_query * term_size / n_universe,
      fold_enrichment = overlap / expected,
      p_value = map2_dbl(overlap, term_size, function(k, term_n) {
        phyper(k - 1, term_n, n_universe - term_n, n_query, lower.tail = FALSE)
      })
    )

  # BH correction is across every term with >= min_term_size genes that
  # was actually tested (i.e. `tested`, before the overlap > 0 filter),
  # since untested terms (zero overlap) were never candidates for
  # significance and shouldn't shrink the correction's denominator.
  result$p_adjust <- p.adjust(
    result$p_value,
    method = "BH",
    n = nrow(tested)
  )

  if (!is.null(term2name)) {
    result <- result %>% left_join(term2name, by = "go_id")
    n_unmapped <- sum(is.na(result$go_name))
    if (n_unmapped > 0) {
      warning(n_unmapped, " GO term(s) in term2gene have no name/namespace in term2name ",
              "(e.g. an obsolete or alternate ID not resolved when the mapping was built) -- dropped.")
      result <- result %>% filter(!is.na(go_name))
    }
    result <- result %>% relocate(go_name, go_namespace, .after = go_id)
  }

  result %>%
    filter(p_adjust <= qval_threshold) %>%
    arrange(p_adjust, p_value)
}

# Standard clusterProfiler-style dotplot: GeneRatio (genes in term that
# are in the query / query size) on x, terms on y, faceted by GO domain
# (BP/CC/MF), dot size = gene count in the term, dot colour = p.adjust
# (red = most significant, blue = least, matching clusterProfiler's
# default palette). stat = "p_adjust" (default, for genuinely significant
# results) or "p_value" -- use the latter for a "borderline terms" view,
# where BH-adjusted p is often capped at exactly 1 for every term
# (collapsing the colour scale to a single value) even though the raw
# p-values still differ meaningfully from each other.
plot_go_enrichment <- function(go_result,
                               label,
                               n_query,
                               top_n_per_domain = 5,
                               stat = c("p_adjust", "p_value"),
                               out_dir = here("results", "figures", "go_enrichment")) {
  stat <- match.arg(stat)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (nrow(go_result) == 0) {
    message("No enriched GO terms to plot for: ", label)
    return(NULL)
  }

  domain_abbrev <- c(
    biological_process = "BP",
    cellular_component = "CC",
    molecular_function = "MF"
  )

  plot_df <- go_result %>%
    mutate(
      gene_ratio = overlap / n_query,
      domain = factor(domain_abbrev[go_namespace], levels = c("BP", "CC", "MF")),
      # Some GO names run to 100+ characters (e.g. certain enzyme activity
      # terms) and would otherwise blow out the left margin and push the
      # whole plot (title, legend, axis) off canvas -- wrap to multiple
      # lines instead.
      go_name_wrapped = vapply(go_name, function(x) paste(strwrap(x, width = 50), collapse = "\n"), character(1))
    ) %>%
    group_by(domain) %>%
    slice_min(.data[[stat]], n = top_n_per_domain, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(domain, gene_ratio) %>%
    # facet_grid's free y-scale is per-panel, but a shared factor's levels
    # are global -- prefixing with domain keeps each facet's terms in
    # their own gene_ratio order without terms from different facets
    # (possibly sharing a name-adjacent sort position) interleaving.
    mutate(term_key = factor(paste(domain, go_name_wrapped, sep = "___"), levels = paste(domain, go_name_wrapped, sep = "___")))

  legend_name <- if (stat == "p_adjust") "p.adjust" else "p.value"

  p <- ggplot(plot_df, aes(x = gene_ratio, y = term_key, size = overlap, colour = .data[[stat]])) +
    geom_point() +
    facet_grid(domain ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
    scale_colour_gradient(low = "red", high = "blue", name = legend_name) +
    scale_size_continuous(name = "Count") +
    labs(
      title = paste("GO enrichment:", label),
      x = "GeneRatio",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(strip.text.y = element_text(angle = 0))

  out_file <- file.path(out_dir, paste0("go_enrichment_", label, ".png"))
  ggsave(out_file, p, width = 9, height = max(4, 0.4 * nrow(plot_df) + 1.5), dpi = 300)

  invisible(p)
}
