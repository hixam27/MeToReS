#!/usr/bin/env Rscript
# run_maaslin2.R
# Associate gene expression with measured element concentrations using
# MaAsLin2 (Multivariate Association with Linear Models).
#
# ---------------------------------------------------------------------------
# DESIGN NOTE — why one element per model
#
# The intuitive design is:
#     expression ~ Al + As + Ba + ... + Zn + (1|country)
# With 19 elements and 14 samples that model has 19 predictors and 14
# observations. It is rank-deficient and cannot be fitted; MaAsLin2 would
# error or silently drop terms.
#
# This script therefore fits ONE element at a time:
#     expression ~ <element> + (1|country)
# looping over the elements requested. Each coefficient is interpretable, each
# model is fittable, and FDR is applied by MaAsLin2 within each run. The
# --combine_fdr option additionally re-applies BH across the pooled results of
# all elements, which is the more conservative and more honest correction when
# you intend to scan every element.
#
# ---------------------------------------------------------------------------
# TRANSFORMATION — why log10 THEN z-score
#
# Environmental concentrations are heavily right-skewed, so a few high samples
# would otherwise dominate any linear model.
#   - log10(x + 1) compresses the right tail: this is what actually fixes skew.
#   - z-score ((x - mean)/sd) is LINEAR: it shifts and rescales but does NOT
#     change distribution shape, so on its own it does not fix skew at all.
#     Its real value here is that it puts every element on a common scale, so
#     MaAsLin2 coefficients become "change per 1 SD of element" and are
#     comparable across elements with different units.
# Applying log10 then z-score gets both properties. That is the default.
#
# (For rank-based methods such as Spearman, all of this is a no-op: monotonic
#  transforms leave ranks unchanged. Transformation matters for LINEAR models
#  like MaAsLin2, which is exactly what this script runs.)
#
# ---------------------------------------------------------------------------
# SAMPLE SIZE — read before interpreting anything
#
# With ~14 samples this is an EXPLORATORY screen, not a confirmatory test.
# Adding country as a random effect over 6 levels with ~2 samples per level
# makes the random-effect variance poorly identified. The script warns when
# the design is thin and writes the diagnostics to the run log. Treat hits as
# candidates to follow up, not as established associations.
# ---------------------------------------------------------------------------
#
# Usage:
#   Rscript run_maaslin2.R \
#     --vst          results/deg/country_deg_results/vst_normalised_counts.csv \
#     --chemical     results/chemical/chemical_data_bioavailable.csv \
#     --samplesheet  samplesheet_mobiles.csv \
#     --out_dir      results/maaslin2_country \
#     --elements     "Cu Zn Pb As Cd" \
#     --random_effects "country" \
#     --transform    log10_zscore \
#     --top_var_genes 2000

suppressPackageStartupMessages({
  library(optparse)
  library(Maaslin2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--vst",            type = "character"),
  make_option("--chemical",       type = "character"),
  make_option("--samplesheet",    type = "character"),
  make_option("--out_dir",        type = "character"),
  make_option("--elements",       type = "character", default = "",
              help = "Space-separated elements to test. Empty = all in the chemical file."),
  make_option("--random_effects", type = "character", default = "country",
              help = "Space-separated samplesheet columns as random effects. Empty = none."),
  make_option("--fixed_covariates", type = "character", default = "",
              help = paste0("Space-separated samplesheet columns added as FIXED covariates ",
                            "alongside the element (use sparingly at low n).")),
  make_option("--transform",      type = "character", default = "log10_zscore",
              help = "Chemistry transform: log10_zscore (default), log10, zscore, none"),
  make_option("--top_var_genes",  type = "integer", default = 2000,
              help = "Restrict to this many most-variable genes"),
  make_option("--gene_subset",    type = "character", default = NULL,
              help = "Optional file of gene IDs (one per line) to restrict to"),
  make_option("--min_prevalence", type = "double", default = 0.1),
  make_option("--max_significance", type = "double", default = 0.05),
  make_option("--combine_fdr",    action = "store_true", default = TRUE,
              help = "Re-apply BH across the pooled results of all elements"),
  make_option("--cores",          type = "integer", default = 1)
)))

for (req in c("vst", "chemical", "samplesheet", "out_dir"))
  if (is.null(opt[[req]])) stop(paste0("--", req, " is required"))

# Resolve paths BEFORE any setwd by MaAsLin2 internals
opt$vst         <- normalizePath(opt$vst,         mustWork = FALSE)
opt$chemical    <- normalizePath(opt$chemical,    mustWork = FALSE)
opt$samplesheet <- normalizePath(opt$samplesheet, mustWork = FALSE)
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
opt$out_dir     <- normalizePath(opt$out_dir,     mustWork = FALSE)
if (!is.null(opt$gene_subset))
  opt$gene_subset <- normalizePath(opt$gene_subset, mustWork = FALSE)

split_ws <- function(x) {
  if (is.null(x)) return(character(0))
  v <- strsplit(trimws(x), "\\s+")[[1]]
  v[nchar(v) > 0]
}
RANDOM <- split_ws(opt$random_effects)
FIXCOV <- split_ws(opt$fixed_covariates)

cat("========================================\n")
cat("MaAsLin2: expression ~ element (+ covariates)\n")
cat("========================================\n")
cat("Transform      :", opt$transform, "\n")
cat("Random effects :", if (length(RANDOM)) paste(RANDOM, collapse = ", ") else "none", "\n")
cat("Fixed covars   :", if (length(FIXCOV)) paste(FIXCOV, collapse = ", ") else "none", "\n")

# ── Load expression (genes x samples -> samples x genes) ──────────────
vst <- read.csv(opt$vst, header = TRUE, check.names = FALSE)
rownames(vst) <- vst[[1]]; vst[[1]] <- NULL
expr <- as.matrix(vst)

# ── Load chemistry ───────────────────────────────────────────────────
chem <- read.csv(opt$chemical, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("SampleID" %in% colnames(chem))
rownames(chem) <- chem$SampleID
chem_meta_cols <- intersect(c("SampleID", "Country", "Season", "Replicate"),
                            colnames(chem))
elem_all  <- setdiff(colnames(chem), chem_meta_cols)
chem_mat  <- chem[, elem_all, drop = FALSE]
chem_mat[] <- lapply(chem_mat, as.numeric)

# ── Load samplesheet (for random effects / covariates) ────────────────
ss <- read.csv(opt$samplesheet, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("sample" %in% colnames(ss))
rownames(ss) <- ss$sample
need_cols <- c(RANDOM, FIXCOV)
missing_cols <- setdiff(need_cols, colnames(ss))
if (length(missing_cols))
  stop("Column(s) not in samplesheet: ", paste(missing_cols, collapse = ", "))

# ── Match samples ────────────────────────────────────────────────────
common <- Reduce(intersect, list(colnames(expr), rownames(chem_mat), rownames(ss)))
cat("Samples common to expression, chemistry and samplesheet:", length(common), "\n")
if (length(common) < 6) stop("Too few matched samples (<6).")
expr     <- expr[, common, drop = FALSE]
chem_mat <- chem_mat[common, , drop = FALSE]
ss_sub   <- ss[common, , drop = FALSE]

# ── Transform chemistry ──────────────────────────────────────────────
# See the header note: log10 fixes skew; z-score only rescales (making
# coefficients comparable across elements). log10_zscore does both.
apply_transform <- function(m, how) {
  if (how == "none") return(m)
  if (how %in% c("log10", "log10_zscore")) {
    if (any(m < 0, na.rm = TRUE))
      warning("Negative concentrations present; log10(x+1) will produce NaN.")
    m <- log10(m + 1)
  }
  if (how %in% c("zscore", "log10_zscore")) {
    m <- as.data.frame(scale(m))
  }
  m
}
if (!opt$transform %in% c("none", "log10", "zscore", "log10_zscore"))
  stop("Unknown --transform: ", opt$transform)
chem_mat <- apply_transform(chem_mat, opt$transform)

# Drop elements with no variance after transform
keep_elem <- vapply(chem_mat, function(x) {
  v <- var(x, na.rm = TRUE); is.finite(v) && v > 0
}, logical(1))
chem_mat <- chem_mat[, keep_elem, drop = FALSE]

ELEMENTS <- split_ws(opt$elements)
if (!length(ELEMENTS)) ELEMENTS <- colnames(chem_mat)
ELEMENTS <- intersect(ELEMENTS, colnames(chem_mat))
if (!length(ELEMENTS)) stop("None of the requested elements are usable.")
cat("Elements to test:", length(ELEMENTS), "->", paste(ELEMENTS, collapse = ", "), "\n")

# ── Restrict genes ───────────────────────────────────────────────────
if (!is.null(opt$gene_subset) && file.exists(opt$gene_subset)) {
  wanted <- trimws(readLines(opt$gene_subset)); wanted <- wanted[nchar(wanted) > 0]
  expr <- expr[rownames(expr) %in% wanted, , drop = FALSE]
  cat("Genes restricted to subset:", nrow(expr), "\n")
} else {
  gv   <- apply(expr, 1, var)
  keep <- order(gv, decreasing = TRUE)[seq_len(min(opt$top_var_genes, nrow(expr)))]
  expr <- expr[keep, , drop = FALSE]
  cat("Genes restricted to top", nrow(expr), "most variable\n")
}
expr <- expr[apply(expr, 1, var) > 0, , drop = FALSE]
cat("Genes after variance filter:", nrow(expr), "\n")

# MaAsLin2 wants samples as ROWS
feat <- as.data.frame(t(expr), check.names = FALSE)

# ── Degrees-of-freedom sanity check ──────────────────────────────────
n <- length(common)
n_fixed <- 1 + length(FIXCOV)      # element + fixed covariates
cat("\nDesign check: n =", n, "| fixed terms per model =", n_fixed, "\n")
if (length(RANDOM)) {
  for (rv in RANDOM) {
    lv <- length(unique(ss_sub[[rv]]))
    cat("  random effect '", rv, "': ", lv, " levels, ",
        round(n / lv, 1), " samples/level\n", sep = "")
    if (n / lv < 3)
      cat("    WARNING: <3 samples per level — the random-effect variance is",
          "poorly identified. Results are exploratory.\n")
  }
}
if (n - n_fixed < 5)
  cat("  WARNING: very few residual df. Treat results as exploratory only.\n")

# ── Run MaAsLin2, one element at a time ──────────────────────────────
all_res <- list()
for (el in ELEMENTS) {
  cat("\n=== MaAsLin2: ", el, " ===\n", sep = "")
  meta <- data.frame(chem_mat[, el, drop = FALSE], check.names = FALSE)
  colnames(meta) <- el
  for (cv in unique(c(RANDOM, FIXCOV))) meta[[cv]] <- ss_sub[[cv]]
  rownames(meta) <- common

  el_out <- file.path(opt$out_dir, paste0("element_", el))

  res <- tryCatch(
    Maaslin2(
      input_data       = feat,
      input_metadata   = meta,
      output           = el_out,
      fixed_effects    = c(el, FIXCOV),
      random_effects   = if (length(RANDOM)) RANDOM else NULL,
      # VST output is already normalised and variance-stabilised, so both of
      # MaAsLin2's own steps are switched off — applying TSS/LOG on top would
      # double-transform the data.
      normalization    = "NONE",
      transform        = "NONE",
      analysis_method  = "LM",
      min_prevalence   = opt$min_prevalence,
      max_significance = opt$max_significance,
      cores            = opt$cores,
      plot_heatmap     = TRUE,
      plot_scatter     = TRUE
    ),
    error = function(e) { cat("  ERROR for ", el, ": ", conditionMessage(e), "\n", sep = ""); NULL }
  )

  rf <- file.path(el_out, "all_results.tsv")
  if (file.exists(rf)) {
    df <- read.delim(rf, stringsAsFactors = FALSE)
    df$feature <- gsub("\\.", "-", df$feature)
    df <- df[df$metadata == el, , drop = FALSE]
    if (nrow(df)) { df$element <- el; all_res[[el]] <- df }
  }
}

if (!length(all_res)) stop("MaAsLin2 produced no results for any element.")

res <- do.call(rbind, all_res)

# ── Pooled FDR across all elements ───────────────────────────────────
# MaAsLin2's qval corrects WITHIN each element's run. Scanning every element
# means many more tests than any single run saw, so re-correct across the pool.
if (isTRUE(opt$combine_fdr)) {
  res$qval_across_elements <- p.adjust(res$pval, method = "BH")
  cat("\nApplied BH across all", nrow(res), "element-gene tests.\n")
}
ord <- if ("qval_across_elements" %in% colnames(res)) res$qval_across_elements else res$qval
res <- res[order(ord), ]

write.table(res, file.path(opt$out_dir, "maaslin2_all_elements.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: maaslin2_all_elements.tsv (", nrow(res), "tests )\n")

qcol <- if ("qval_across_elements" %in% colnames(res)) "qval_across_elements" else "qval"
sig  <- res[!is.na(res[[qcol]]) & res[[qcol]] < opt$max_significance, ]
write.table(sig, file.path(opt$out_dir, "maaslin2_significant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved: maaslin2_significant.tsv (", nrow(sig), "at", qcol, "<",
    opt$max_significance, ")\n")

if (nrow(sig)) {
  tab <- as.data.frame(table(sig$element))
  colnames(tab) <- c("element", "n_significant_genes")
  tab <- tab[order(-tab$n_significant_genes), ]
  write.table(tab, file.path(opt$out_dir, "significant_genes_per_element.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nSignificant genes per element:\n")
  print(tab, row.names = FALSE)
} else {
  cat("\nNo significant associations after correction.\n")
  cat("At n =", n, "this is an entirely plausible outcome and is worth\n")
  cat("reporting as such rather than relaxing thresholds to find hits.\n")
}

cat("\n========================================\n")
cat("MaAsLin2 complete. Exploratory at n =", n, "\n")
cat("========================================\n")
