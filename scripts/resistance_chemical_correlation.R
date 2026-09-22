#!/usr/bin/env Rscript
# resistance_chemical_correlation.R
# Produces two outputs:
#   (1) element_compound_correlations.csv — direct matches only (e.g. Cu MRG
#       expression vs Cu concentration), the cleanest test of co-selection.
#   (2) full_correlation_matrix.csv + correlation_heatmap.pdf — every
#       ARG/MRG category against every measured element, for broader
#       co-selection exploration (e.g. does Cu exposure also select for
#       antibiotic resistance genes?).
#
# BH correction is applied ONCE across the combined ARG+MRG resistome
# testing family, not separately per source — a single correlation is one
# test in one family, regardless of whether its category happens to be an
# ARG or an MRG. element_compound_correlations.csv is derived from that
# same corrected table, so its padj values are on the same scale as
# full_correlation_matrix.csv's, not a second, inconsistent correction.

suppressPackageStartupMessages({
  library(optparse)
  library(pheatmap)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--by_aro",      type = "character"),
  make_option("--by_compound", type = "character"),
  make_option("--chemical",    type = "character"),
  make_option("--fraction_label", type = "character", default = "bioavailable"),
  make_option("--method",      type = "character", default = "spearman"),
  make_option("--min_detected_samples", type = "integer", default = 4,
              help = paste("Minimum samples where a resistance category must have",
                            "nonzero expression before it's tested. Distinct from the",
                            "sample-overlap check below: a category can be present in",
                            "many matched samples yet detected (nonzero) in very few,",
                            "in which case a correlation would be driven by a handful",
                            "of points. [default: %default]")),
  make_option("--out_dir",     type = "character")
)))

# Resolve all input paths to absolute BEFORE changing the working directory —
# otherwise relative paths (e.g. results/resistance_contaminant/...) would
# no longer resolve once we setwd() into out_dir below.
opt$by_aro      <- normalizePath(opt$by_aro,      mustWork = FALSE)
opt$by_compound <- normalizePath(opt$by_compound, mustWork = FALSE)
opt$chemical    <- normalizePath(opt$chemical,    mustWork = FALSE)
opt$out_dir     <- normalizePath(opt$out_dir,     mustWork = FALSE)

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$out_dir)

cat("========================================\n")
cat("Resistance x Chemical Correlation (", opt$fraction_label, ")\n")
cat("========================================\n")
cat("Minimum detected samples per category:", opt$min_detected_samples, "\n")

# Common element aliases (symbol <-> a few likely full-name variants seen in
# BacMet compound fields). Extend this table if new mismatches turn up.
ALIASES <- list(
  Cu = c("Cu", "Copper"), Zn = c("Zn", "Zinc"), As = c("As", "Arsenic"),
  Cd = c("Cd", "Cadmium"), Co = c("Co", "Cobalt"), Cr = c("Cr", "Chromium"),
  Ni = c("Ni", "Nickel"), Pb = c("Pb", "Lead"), Mn = c("Mn", "Manganese"),
  Fe = c("Fe", "Iron"), Al = c("Al", "Aluminium", "Aluminum"),
  Ba = c("Ba", "Barium"), Ca = c("Ca", "Calcium"), K = c("K", "Potassium"),
  Mg = c("Mg", "Magnesium"), Na = c("Na", "Sodium"), P = c("P", "Phosphorus"),
  S = c("S", "Sulfur", "Sulphur"), Sr = c("Sr", "Strontium")
)
.exact_alias_hit <- function(name) {
  hit <- names(ALIASES)[sapply(ALIASES, function(v) tolower(name) %in% tolower(v))]
  if (length(hit)) hit[1] else NA
}

# BacMet compound names are frequently written as "Cadmium (Cd)" or
# "Zinc (Zn)" rather than the bare symbol or name alone, so an exact match
# against ALIASES fails on essentially every real row. Extract a short
# parenthetical token (1-2 letters, e.g. the "(Cd)" in "Cadmium (Cd)") and
# try that first; fall back to the original whole-name exact match for
# compounds with no parenthetical symbol (e.g. a bare "Copper").
alias_lookup <- function(name) {
  paren <- regmatches(name, regexpr("\\(([A-Za-z]{1,2})\\)", name))
  if (length(paren) && nzchar(paren)) {
    symbol <- sub("^\\(", "", sub("\\)$", "", paren))
    hit <- .exact_alias_hit(symbol)
    if (!is.na(hit)) return(hit)
  }
  .exact_alias_hit(name)
}

# ── Load chemical data (long format: SampleID, Element, Concentration) ──
chem <- read.csv(opt$chemical, stringsAsFactors = FALSE, check.names = FALSE)
meta_cols <- intersect(c("SampleID", "Country", "Season", "Replicate"), colnames(chem))
elem_cols <- setdiff(colnames(chem), meta_cols)
chem_long <- do.call(rbind, lapply(elem_cols, function(e) {
  data.frame(SampleID = chem$SampleID, Element = e,
             Concentration = as.numeric(chem[[e]]), stringsAsFactors = FALSE)
}))

# ── Helper: load a resistance abundance table into long format ─────────
load_res_long <- function(path, category_col_candidates) {
  if (!file.exists(path)) {
    cat("  NOTE:", path, "does not exist — skipping.\n")
    return(NULL)
  }
  if (file.info(path)$size == 0) {
    cat("  NOTE:", path, "is empty — skipping.\n")
    return(NULL)
  }
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  cat_col <- intersect(category_col_candidates, colnames(df))[1]
  if (is.na(cat_col)) cat_col <- colnames(df)[1]
  sample_cols <- setdiff(colnames(df), cat_col)
  long <- do.call(rbind, lapply(sample_cols, function(s) {
    data.frame(Category = df[[cat_col]], SampleID = s,
               Expression = as.numeric(df[[s]]), stringsAsFactors = FALSE)
  }))
  long
}

arg_long <- load_res_long(opt$by_aro,      c("Best_Hit_ARO", "ARO"))
mrg_long <- load_res_long(opt$by_compound, c("Compound"))

if (is.null(arg_long) && is.null(mrg_long)) {
  stop("Neither ARG nor MRG abundance table could be loaded.\n",
       "  Checked --by_aro:      ", opt$by_aro, "\n",
       "  Checked --by_compound: ", opt$by_compound)
}

# ── Correlation helper ──────────────────────────────────────────────
# Returns UNADJUSTED p-values. BH correction is applied once, after
# combining ARG and MRG results into a single resistome testing family
# (see below) — not here, per-source.
correlate_pair <- function(res_long, source_label) {
  if (is.null(res_long)) return(NULL)
  categories <- unique(res_long$Category)
  elements   <- unique(chem_long$Element)
  out <- list()
  for (cat in categories) {
    exp_sub <- res_long[res_long$Category == cat, ]
    for (el in elements) {
      chem_sub <- chem_long[chem_long$Element == el, ]
      m <- merge(exp_sub, chem_sub, by = "SampleID")
      if (nrow(m) < 4) next

      # Detection floor: require the category to be DETECTED (nonzero
      # expression) in enough of the matched samples, not merely present
      # in the merged sample set (see --min_detected_samples above).
      n_detected <- sum(m$Expression > 0, na.rm = TRUE)
      if (n_detected < opt$min_detected_samples) next

      ct <- suppressWarnings(cor.test(m$Expression, m$Concentration,
                                       method = opt$method))
      out[[length(out) + 1]] <- data.frame(
        Source = source_label, Category = cat, Element = el,
        n = nrow(m), n_detected = n_detected,
        rho = unname(ct$estimate), pvalue = ct$p.value,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

cat("\nComputing full correlation matrix (this may take a moment)...\n")
arg_corr <- correlate_pair(arg_long, "ARG")
mrg_corr <- correlate_pair(mrg_long, "MRG")
full_corr <- do.call(rbind, Filter(Negate(is.null), list(arg_corr, mrg_corr)))

if (is.null(full_corr) || nrow(full_corr) == 0) {
  stop("No correlations could be computed — check sample overlap/detection ",
       "between resistance tables and chemical data, or lower ",
       "--min_detected_samples.")
}

# Single BH correction across the WHOLE resistome-testing family (ARG + MRG
# combined) — one family, matching what is actually being tested, rather
# than two independently-corrected halves stacked together.
full_corr$padj <- p.adjust(full_corr$pvalue, method = "BH")
full_corr <- full_corr[order(full_corr$padj), ]
write.csv(full_corr, "full_correlation_matrix.csv", row.names = FALSE)
cat("Saved: full_correlation_matrix.csv (", nrow(full_corr), "pairs, ",
    "BH-corrected as one family)\n")

# ── Direct element<->compound matches (the clean co-selection test) ────
# Derived from full_corr (already correctly, uniformly BH-corrected) rather
# than re-filtering the pre-correction mrg_corr object.
mrg_rows <- full_corr[full_corr$Source == "MRG", ]
if (nrow(mrg_rows) > 0) {
  mrg_rows$MatchedElement <- sapply(mrg_rows$Category, alias_lookup)
  direct <- mrg_rows[!is.na(mrg_rows$MatchedElement) &
                      mrg_rows$MatchedElement == mrg_rows$Element, ]
  direct <- direct[order(direct$padj), ]
  write.csv(direct, "element_compound_correlations.csv", row.names = FALSE)
  cat("Saved: element_compound_correlations.csv (", nrow(direct),
      "direct element-compound pairs)\n")
  if (nrow(direct) > 0) {
    cat("\nDirect matches summary (compound expression vs matching element):\n")
    print(direct[, c("Category", "Element", "n", "n_detected", "rho", "pvalue", "padj")],
          row.names = FALSE)
  }
} else {
  cat("No MRG rows available — skipping direct element-compound matching.\n")
  write.csv(data.frame(), "element_compound_correlations.csv", row.names = FALSE)
}

# ── Heatmap of the full correlation matrix ──────────────────────────
cat("\n=== Correlation heatmap ===\n")

mat <- reshape(full_corr[, c("Category", "Element", "rho")],
               idvar = "Category", timevar = "Element", direction = "wide")
rownames(mat) <- mat$Category
mat$Category <- NULL
colnames(mat) <- sub("^rho\\.", "", colnames(mat))
mat <- as.matrix(mat)
mat[is.na(mat)] <- 0

if (nrow(mat) >= 2 && ncol(mat) >= 2) {
  pal <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(50)
  pheatmap(mat, color = pal, breaks = seq(-1, 1, length.out = 51),
           cluster_rows = nrow(mat) > 2, cluster_cols = ncol(mat) > 2,
           fontsize_row = 7, fontsize_col = 8,
           main = paste0("ARG/MRG vs chemical element correlation (",
                         opt$method, ", ", opt$fraction_label, ")"),
           filename = "correlation_heatmap.pdf", width = 10,
           height = max(6, nrow(mat) * 0.15))
  cat("Saved: correlation_heatmap.pdf\n")
} else {
  pdf("correlation_heatmap.pdf"); plot.new()
  title("Not enough categories/elements for a heatmap"); dev.off()
  cat("Skipped heatmap (matrix too small)\n")
}

cat("\n========================================\n")
cat("Resistance-chemical correlation complete.\n")
cat("========================================\n")
