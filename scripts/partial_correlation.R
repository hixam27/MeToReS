#!/usr/bin/env Rscript
# partial_correlation.R
# Bioindicator discovery by partial (Spearman) correlation.
#
# Finds genes whose expression correlates with element concentrations while
# statistically controlling for categorical confounders (e.g. country, season)
# via ppcor::pcor.test. This directly addresses the country/contaminant
# confounding: it isolates the chemistry-expression association from the
# variance explained by the confounders.
#
# IMPORTANT — sample size. With only ~14 samples, controlling for country
# (6 levels -> 5 dummy covariates) plus season (1 dummy) leaves very few
# residual degrees of freedom and low power. Two safeguards are built in:
#   (1) covariates are configurable (--covariates), so you can control for
#       season alone to preserve df;
#   (2) genes are restricted to a manageable, biologically-focused set
#       (--gene_subset file, or the top --top_var_genes most-variable genes),
#       both to reduce the multiple-testing burden and to keep runtime sane.
# Results should be read as exploratory bioindicator candidates, not
# definitive associations.
#
# Chemistry is transformed before correlation (--transform), because
# environmental concentrations are heavily right-skewed; log10(x + 1) is the
# default.
#
# Usage:
#   Rscript partial_correlation.R \
#     --vst        results/deg/country_deg_results/vst_normalised_counts.csv \
#     --chemical   results/chemical/chemical_data_bioavailable.csv \
#     --samplesheet samplesheet_mobiles.csv \
#     --out_dir    results/bioindicators_country \
#     --covariates "country season" \
#     --transform  log10 \
#     --top_var_genes 2000 \
#     --r_threshold 0.6 --padj_threshold 0.05

suppressPackageStartupMessages({
  library(optparse)
  library(ppcor)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--vst",            type = "character"),
  make_option("--chemical",       type = "character"),
  make_option("--samplesheet",    type = "character"),
  make_option("--out_dir",        type = "character"),
  make_option("--covariates",     type = "character", default = "country season",
              help = "Space-separated samplesheet columns to control for"),
  make_option("--transform",      type = "character", default = "log10",
              help = paste0("Chemistry transform: log10 (log10(x+1)), zscore, log10_zscore, none. ",
                            "NOTE: with --method spearman this is a NO-OP (see below).")),
  make_option("--gene_subset",    type = "character", default = NULL,
              help = "Optional file of gene IDs (one per line) to restrict to"),
  make_option("--top_var_genes",  type = "integer", default = 2000,
              help = "If no --gene_subset, use this many most-variable genes"),
  make_option("--method",         type = "character", default = "spearman"),
  make_option("--r_threshold",    type = "double", default = 0.6),
  make_option("--padj_threshold", type = "double", default = 0.05)
)))

for (req in c("vst", "chemical", "samplesheet", "out_dir"))
  if (is.null(opt[[req]])) stop(paste0("--", req, " is required"))

# Resolve paths BEFORE changing directory
opt$vst         <- normalizePath(opt$vst,         mustWork = FALSE)
opt$chemical    <- normalizePath(opt$chemical,    mustWork = FALSE)
opt$samplesheet <- normalizePath(opt$samplesheet, mustWork = FALSE)
opt$out_dir     <- normalizePath(opt$out_dir,     mustWork = FALSE)
if (!is.null(opt$gene_subset))
  opt$gene_subset <- normalizePath(opt$gene_subset, mustWork = FALSE)
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(opt$out_dir)

COVARS <- strsplit(trimws(opt$covariates), "\\s+")[[1]]
COVARS <- COVARS[nchar(COVARS) > 0]

cat("========================================\n")
cat("Partial correlation (bioindicators)\n")
cat("========================================\n")
cat("Covariates controlled for:", paste(COVARS, collapse = ", "), "\n")
cat("Chemistry transform      :", opt$transform, "\n")

# ── Load ──────────────────────────────────────────────────────────────
vst <- read.csv(opt$vst, header = TRUE, check.names = FALSE)
rownames(vst) <- vst[[1]]; vst[[1]] <- NULL
expr <- as.matrix(vst)                              # genes x samples

chem <- read.csv(opt$chemical, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("SampleID" %in% colnames(chem))
rownames(chem) <- chem$SampleID
meta_cols <- intersect(c("SampleID", "Country", "Season", "Replicate"),
                       colnames(chem))
elem_cols <- setdiff(colnames(chem), meta_cols)
chem_mat  <- chem[, elem_cols, drop = FALSE]
chem_mat[] <- lapply(chem_mat, as.numeric)

ss <- read.csv(opt$samplesheet, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("sample" %in% colnames(ss))
rownames(ss) <- ss$sample
missing_cov <- setdiff(COVARS, colnames(ss))
if (length(missing_cov))
  stop("Covariate(s) not in samplesheet: ", paste(missing_cov, collapse = ", "))

# ── Match samples across all three tables ───────────────────────────
common <- Reduce(intersect, list(colnames(expr), rownames(chem_mat), rownames(ss)))
cat("Samples common to expression, chemistry and samplesheet:", length(common), "\n")
if (length(common) < 6)
  stop("Too few matched samples (<6) for partial correlation.")
expr     <- expr[, common, drop = FALSE]
chem_mat <- chem_mat[common, , drop = FALSE]
covar_df <- ss[common, COVARS, drop = FALSE]

# ── Transform chemistry ─────────────────────────────────────────────
# Transform the chemistry.
#
# HONEST NOTE ON WHAT THIS DOES:
#   - log10(x + 1) compresses the right tail and is what actually fixes the
#     skew of environmental concentrations.
#   - z-score is LINEAR: it shifts and rescales but leaves distribution shape
#     untouched, so on its own it does NOT reduce the influence of outliers.
#     Its value is putting elements on a common scale for LINEAR models.
#   - If --method is "spearman" (the default here), BOTH are a complete no-op:
#     Spearman ranks the data first, and any monotonic transform leaves ranks
#     unchanged. Raw, logged and z-scored chemistry give IDENTICAL Spearman
#     partial correlations. The option matters only for --method pearson.
if (opt$transform %in% c("log10", "log10_zscore")) {
  if (any(chem_mat < 0, na.rm = TRUE))
    warning("Negative concentrations present; log10(x+1) will produce NaN.")
  chem_mat <- log10(chem_mat + 1)
}
if (opt$transform %in% c("zscore", "log10_zscore")) {
  chem_mat <- as.data.frame(scale(chem_mat))
}
if (!opt$transform %in% c("none", "log10", "zscore", "log10_zscore")) {
  stop("Unknown --transform: ", opt$transform)
}
if (opt$method == "spearman" && opt$transform != "none") {
  cat("Note: --transform", opt$transform, "has no effect on Spearman results",
      "(rank-based); reported for transparency.\n")
}
# Drop elements with no variance after transform
chem_mat <- chem_mat[, apply(chem_mat, 2, function(x) var(x, na.rm = TRUE) > 0),
                     drop = FALSE]
cat("Elements retained:", ncol(chem_mat), "\n")

# ── Build numeric covariate design (dummy-code factors) ─────────────
for (cv in COVARS) covar_df[[cv]] <- factor(covar_df[[cv]])
# model.matrix without intercept column -> dummy variables for each factor
covar_mm <- model.matrix(~ ., data = covar_df)
covar_mm <- covar_mm[, colnames(covar_mm) != "(Intercept)", drop = FALSE]
# Drop any zero-variance / aliased covariate columns
covar_mm <- covar_mm[, apply(covar_mm, 2, function(x) var(x) > 0), drop = FALSE]
n_resid_df <- length(common) - ncol(covar_mm) - 2   # x, y, and covariates
cat("Covariate columns:", ncol(covar_mm),
    "| approx residual df per test:", n_resid_df, "\n")
if (n_resid_df < 3)
  cat("  WARNING: very low residual df — results will be unstable. Consider",
      "controlling for fewer covariates (e.g. season only).\n")

# ── Restrict genes ──────────────────────────────────────────────────
if (!is.null(opt$gene_subset) && file.exists(opt$gene_subset)) {
  wanted <- readLines(opt$gene_subset)
  wanted <- trimws(wanted); wanted <- wanted[nchar(wanted) > 0]
  expr <- expr[rownames(expr) %in% wanted, , drop = FALSE]
  cat("Genes restricted to subset:", nrow(expr), "\n")
} else {
  gene_var <- apply(expr, 1, var)
  keep <- order(gene_var, decreasing = TRUE)[seq_len(min(opt$top_var_genes,
                                                         nrow(expr)))]
  expr <- expr[keep, , drop = FALSE]
  cat("Genes restricted to top", nrow(expr), "most variable\n")
}
# Drop zero-variance genes (constant expression breaks correlation)
expr <- expr[apply(expr, 1, var) > 0, , drop = FALSE]
cat("Genes after variance filter:", nrow(expr), "\n")

# ── Partial correlation loop ────────────────────────────────────────
cat("\nComputing partial correlations (",
    nrow(expr), "genes x", ncol(chem_mat), "elements =",
    nrow(expr) * ncol(chem_mat), "tests)...\n")

results <- vector("list", nrow(expr) * ncol(chem_mat))
idx <- 0
gene_ids <- rownames(expr)
elem_ids <- colnames(chem_mat)

for (gi in seq_len(nrow(expr))) {
  gvals <- as.numeric(expr[gi, ])
  for (ei in seq_along(elem_ids)) {
    evals <- as.numeric(chem_mat[, ei])
    ok <- is.finite(gvals) & is.finite(evals)
    if (sum(ok) - ncol(covar_mm) - 2 < 1) next
    pc <- tryCatch(
      pcor.test(gvals[ok], evals[ok], covar_mm[ok, , drop = FALSE],
                method = opt$method),
      error = function(e) NULL
    )
    if (is.null(pc)) next
    idx <- idx + 1
    results[[idx]] <- data.frame(
      gene = gene_ids[gi], element = elem_ids[ei],
      estimate = pc$estimate, p.value = pc$p.value,
      n = sum(ok), stringsAsFactors = FALSE
    )
  }
}
results <- results[seq_len(idx)]
if (idx == 0) stop("No partial correlations could be computed.")
res <- do.call(rbind, results)
res$padj <- p.adjust(res$p.value, method = "BH")
res <- res[order(res$padj, -abs(res$estimate)), ]

write.csv(res, "partial_correlations_all.csv", row.names = FALSE)
cat("Saved: partial_correlations_all.csv (", nrow(res), "pairs)\n")

sig <- res[abs(res$estimate) >= opt$r_threshold &
           res$padj < opt$padj_threshold & !is.na(res$padj), ]
write.csv(sig, "partial_correlations_significant.csv", row.names = FALSE)
cat("Saved: partial_correlations_significant.csv (", nrow(sig),
    "pairs at |r| >=", opt$r_threshold, ", padj <", opt$padj_threshold, ")\n")

# Small summary of candidate bioindicators per element
if (nrow(sig) > 0) {
  tab <- as.data.frame(table(sig$element))
  colnames(tab) <- c("element", "n_candidate_genes")
  tab <- tab[order(-tab$n_candidate_genes), ]
  write.csv(tab, "bioindicator_candidates_per_element.csv", row.names = FALSE)
  cat("\nCandidate bioindicator genes per element:\n")
  print(tab, row.names = FALSE)
}

cat("\n========================================\n")
cat("Partial correlation analysis complete.\n")
cat("========================================\n")
