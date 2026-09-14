# "Guy11M" is not a genetically distinct strain -- the "M" is a
# researcher's initial. Guy11 and Guy11M are the same genetic background,
# just run in different experiments/batches (Exp1 and Exp2 respectively);
# mst7 (Exp1) and pmk1 (Exp2) are labelled the same way here for the same
# reason: the pair only looks like "two different strains" when both
# members appear in the same figure. Apply these ONLY to figures that mix
# strains from both experiments (e.g. all_samples, key_strains, the
# cross-experiment comparisons) -- a figure with only Exp1 or only Exp2
# strains has no ambiguity to resolve and should keep the plain name.
#
# These relabel *display text only* (plot titles, axis ticks, legends).
# Sample IDs, file names, and the `comparisons` table are untouched, so
# nothing about how the pipeline identifies or filters data changes.
relabel_cross_experiment_strains <- function(x) {
  # Each condition is tested against the ORIGINAL string and categories
  # are mutually exclusive, so a label is only ever relabeled once --
  # chaining sequential sub() calls instead would let a later pattern
  # (e.g. "^Guy11") re-match text a prior substitution already inserted
  # (e.g. "Guy11_Exp2_2" still starts with "Guy11_"), double-relabeling it.
  is_guy11m <- grepl("^Guy11M(_|$)", x)
  is_guy11 <- !is_guy11m & grepl("^Guy11(_|$)", x)
  is_mst7 <- grepl("^mst7(_|$)", x)
  is_pmk1 <- grepl("^pmk1(_|$)", x)

  out <- x
  out[is_guy11m] <- sub("^Guy11M", "Guy11_Exp2", x[is_guy11m])
  out[is_guy11] <- sub("^Guy11", "Guy11_Exp1", x[is_guy11])
  out[is_mst7] <- sub("^mst7", "mst7_Exp1", x[is_mst7])
  out[is_pmk1] <- sub("^pmk1", "pmk1_Exp2", x[is_pmk1])
  out
}

# Same relabelling applied to a compound "<control>_vs_<test>" comparison
# label, as used in DGE/rank-products/upset plot titles and set names.
relabel_cross_experiment_comparison <- function(x) {
  vapply(x, function(lab) {
    parts <- strsplit(lab, "_vs_", fixed = TRUE)[[1]]
    paste(relabel_cross_experiment_strains(parts), collapse = "_vs_")
  }, character(1), USE.NAMES = FALSE)
}

# Exact published nomenclature for five strains: gene deletions (Delta,
# gene name in italics) and phosphosite-substitution GFP fusions (protein
# name and GFP in italics; the substitution itself a plain-text
# superscript). "PM"/"PD" are recognised too since some comparison labels
# use that short form for MST7PM/MST7PD. Independent of, and composes
# with, relabel_cross_experiment_strains() above -- applies everywhere
# these names appear, not just cross-experiment figures.
#
# Rendering this requires R's plotmath expression system (base graphics,
# no extra package): ggplot2 titles/axis labels accept a plotmath
# `expression()` in place of plain text wherever italics/superscript are
# needed. It only works where ggplot2 draws the text itself -- UpSetR's
# upset() plot draws its own set-name labels internally with no rich-text
# hook, so upset plots keep plain "mst7"/"pmk1" text.
.strain_nomenclature_dict <- list(
  "mst7"   = 'paste(Delta, italic("mst"), "7")',
  "pmk1"   = 'paste(Delta, italic("pmk"), "1")',
  "MST7WT" = 'paste(italic("MST"), "7-", italic("GFP"))',
  "MST7PM" = 'paste(italic("MST"), "7"^{"S358D"}, "-", italic("GFP"))',
  "MST7PD" = 'paste(italic("MST"), "7"^{"S358A"}, "-", italic("GFP"))',
  "PM"     = 'paste(italic("MST"), "7"^{"S358D"}, "-", italic("GFP"))',
  "PD"     = 'paste(italic("MST"), "7"^{"S358A"}, "-", italic("GFP"))'
)

# One token (e.g. a sample ID like "mst7_0_1", or a "_vs_"-split half of a
# comparison label, possibly already passed through
# relabel_cross_experiment_strains()) -> plotmath source text, meant to be
# embedded as argument(s) inside an outer paste(...) call. Anything not
# matching one of the five names above is returned as a quoted literal, so
# it still combines validly with styled terms.
strain_nomenclature_source <- function(token) {
  for (nm in names(.strain_nomenclature_dict)) {
    if (grepl(paste0("^", nm, "(_|$)"), token)) {
      suffix <- sub(paste0("^", nm), "", token)
      styled <- .strain_nomenclature_dict[[nm]]
      return(if (nzchar(suffix)) paste0(styled, ', "', suffix, '"') else styled)
    }
  }
  paste0('"', token, '"')
}

# Vectorised plotmath expressions for a set of axis-tick labels, in input
# order -- pass directly as e.g. scale_x_discrete(labels = ...).
strain_nomenclature_labels <- function(x) {
  parse(text = vapply(
    x, function(tok) paste0("paste(", strain_nomenclature_source(tok), ")"),
    character(1), USE.NAMES = FALSE
  ))
}

# A single plotmath expression for a compound "<left>_vs_<right>" plot
# title, with `prefix` (plain text, e.g. "DEG counts: ") prepended.
strain_nomenclature_title <- function(display_label, prefix = "") {
  parts <- strsplit(display_label, "_vs_", fixed = TRUE)[[1]]
  src <- if (length(parts) == 2) {
    paste0(
      'paste("', prefix, '", ',
      strain_nomenclature_source(parts[1]), ', "_vs_", ',
      strain_nomenclature_source(parts[2]),
      ")"
    )
  } else {
    paste0('paste("', prefix, display_label, '")')
  }
  parse(text = src)[[1]]
}
