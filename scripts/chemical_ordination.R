#!/usr/bin/env Rscript
# chemical_ordination.R
# Two complementary analyses linking gene expression to measured chemistry:
#   (1) envfit — fit chemical element vectors onto a PCoA of gene expression
#       to test which elements explain the observed expression variance.
#   (2) A standalone PCA of the chemistry itself, to see which elements
#       co-vary and how samples cluster on chemistry alone.

suppressPackageStartupMessages({
  library(optparse)
  library(vegan)
  library(ape)
  library(ggplot2)
  library(ggrepel)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--vst",          type = "character"),
  make_option("--chemical",     type = "character"),
  make_option("--workdir",      type = "character"),
  make_option("--fraction_label", type = "character", default = "bioavailable"),
  make_option("--permutations", type = "integer", default = 999)
)))

dir.create(opt$workdir, showWarnings = FALSE, recursive = TRUE)
opt$vst      <- normalizePath(opt$vst,      mustWork = FALSE)
opt$chemical <- normalizePath(opt$chemical, mustWork = FALSE)
opt$workdir  <- normalizePath(opt$workdir,  mustWork = FALSE)
setwd(opt$workdir)

cat("========================================\n")
cat("Chemical Ordination (", opt$fraction_label, ")\n")
cat("========================================\n")

# ── Load data ─────────────────────────────────────────────────────────
vst <- read.csv(opt$vst, header = TRUE, check.names = FALSE)
rownames(vst) <- vst[[1]]; vst[[1]] <- NULL
comm <- t(as.matrix(vst))   # samples x genes

chem <- read.csv(opt$chemical, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("SampleID" %in% colnames(chem))
rownames(chem) <- chem$SampleID
meta_cols <- intersect("SampleID", colnames(chem))
elem_cols <- setdiff(colnames(chem), meta_cols)
chem_mat  <- chem[, elem_cols, drop = FALSE]
chem_mat[] <- lapply(chem_mat, as.numeric)

common <- intersect(rownames(comm), rownames(chem_mat))
cat("Samples with both expression and chemical data:", length(common), "\n")
if (length(common) < 4) stop("Too few matched samples for ordination (<4).")

comm     <- comm[common, , drop = FALSE]
chem_mat <- chem_mat[common, , drop = FALSE]

# Drop zero-variance genes / elements
comm     <- comm[, apply(comm, 2, var) > 0, drop = FALSE]
chem_mat <- chem_mat[, apply(chem_mat, 2, function(x) var(x, na.rm = TRUE) > 0), drop = FALSE]
cat("Elements retained (non-zero variance):", ncol(chem_mat), "\n")

# ── Part 1: envfit onto a PCoA of gene expression ───────────────────
cat("\n=== envfit: chemical vectors on expression PCoA ===\n")

# Bray-Curtis requires non-negative data. VST output can be slightly
# negative for very low counts — shift to non-negative before computing
# the distance matrix (same approach used in constrained_ordination.R).
comm_for_dist <- comm
if (min(comm_for_dist) < 0) {
  shift <- abs(min(comm_for_dist))
  comm_for_dist <- comm_for_dist + shift
  cat("Note: VST had negative values; shifted by", round(shift, 3),
      "for Bray-Curtis\n")
}
dist_mat <- vegdist(comm_for_dist, method = "bray")
pcoa_res <- pcoa(dist_mat)
scores   <- as.data.frame(pcoa_res$vectors[, 1:2])
colnames(scores) <- c("Axis1", "Axis2")

# Each element is fitted on the samples that have a value for it. A single
# envfit call with na.rm = TRUE would keep only samples complete for EVERY
# element; with scattered gaps that can leave a handful of samples, too few
# for the permutation test to reach any p < 0.05.
min_n <- 5
ef_rows <- list()
for (el in colnames(chem_mat)) {
  keep <- !is.na(chem_mat[[el]])
  n_el <- sum(keep)
  if (n_el < min_n || var(chem_mat[[el]][keep]) == 0) {
    cat("  ", el, ": skipped (", n_el, " samples with a value)\n", sep = "")
    next
  }
  set.seed(42)
  f <- envfit(scores[keep, , drop = FALSE], chem_mat[keep, el, drop = FALSE],
              permutations = opt$permutations)
  ef_rows[[el]] <- data.frame(
    Element   = el,
    NMDS1     = f$vectors$arrows[1, 1],
    NMDS2     = f$vectors$arrows[1, 2],
    r2        = unname(f$vectors$r),
    pvalue    = unname(f$vectors$pvals),
    n_samples = n_el
  )
}
if (length(ef_rows) == 0) stop("No element has enough samples for envfit")
ef_df <- do.call(rbind, ef_rows)
ef_df$padj <- p.adjust(ef_df$pvalue, method = "BH")
cat("Samples per element:", paste0(ef_df$Element, "=", ef_df$n_samples, collapse = " "), "\n")
if (any(ef_df$n_samples < nrow(chem_mat)))
  cat("Note: elements with missing values are tested on fewer samples (n_samples column)\n")
ef_df <- ef_df[order(ef_df$pvalue), ]
write.csv(ef_df, "envfit_results.csv", row.names = FALSE)
cat("Saved: envfit_results.csv\n")
cat("Significant elements (p < 0.05):", sum(ef_df$pvalue < 0.05), "of", nrow(ef_df), "\n")

# Biplot: PCoA points + significant vectors only
sig <- ef_df[ef_df$pvalue < 0.05, ]
arrow_scale <- 0.5 * max(abs(range(scores)))
if (nrow(sig) > 0) {
  sig$NMDS1 <- sig$NMDS1 * sqrt(sig$r2) * arrow_scale
  sig$NMDS2 <- sig$NMDS2 * sqrt(sig$r2) * arrow_scale
}

p <- ggplot(scores, aes(Axis1, Axis2)) +
  geom_point(size = 3, colour = "grey30") +
  geom_text_repel(aes(label = rownames(scores)), size = 2.8, max.overlaps = 20) +
  { if (nrow(sig) > 0)
      geom_segment(data = sig, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                   arrow = arrow(length = unit(0.2, "cm")), colour = "#b2182b",
                   inherit.aes = FALSE)
  } +
  { if (nrow(sig) > 0)
      geom_text_repel(data = sig, aes(x = NMDS1, y = NMDS2, label = Element),
                       colour = "#b2182b", fontface = "bold", size = 3.2,
                       inherit.aes = FALSE)
  } +
  labs(title = paste0("Chemical envfit on expression PCoA (", opt$fraction_label, ")"),
       subtitle = paste0(nrow(sig), " of ", nrow(ef_df), " elements significant (p < 0.05)"),
       x = "PCoA Axis 1", y = "PCoA Axis 2") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40"))

ggsave("envfit_biplot.pdf", p, width = 8.5, height = 7, device = cairo_pdf)
cat("Saved: envfit_biplot.pdf\n")

# ── Part 2: standalone PCA of the chemistry ─────────────────────────
cat("\n=== Standalone chemical PCA ===\n")

n_na <- sum(is.na(chem_mat))
if (n_na > 0) {
  na_by_element <- colSums(is.na(chem_mat))
  na_by_element <- na_by_element[na_by_element > 0]
  cat("NOTE:", n_na, "missing value(s) across", length(na_by_element),
      "element(s) — imputing with column mean for this PCA only:\n")
  for (el in names(na_by_element)) {
    cat("  ", el, ":", na_by_element[el], "missing\n")
    chem_mat[[el]][is.na(chem_mat[[el]])] <- mean(chem_mat[[el]], na.rm = TRUE)
  }
}

chem_scaled <- scale(chem_mat)
chem_pca <- prcomp(chem_scaled, center = FALSE, scale. = FALSE)
pct_var  <- round(100 * summary(chem_pca)$importance[2, 1:2])

pca_df <- as.data.frame(chem_pca$x[, 1:2])
pca_df$SampleID <- rownames(pca_df)

p_pca <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 3, colour = "#4C8DFF") +
  geom_text_repel(aes(label = SampleID), size = 2.8, max.overlaps = 20) +
  labs(title = paste0("Chemical PCA (", opt$fraction_label, ")"),
       x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)")) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("chemical_pca.pdf", p_pca, width = 8, height = 6.5, device = cairo_pdf)
cat("Saved: chemical_pca.pdf\n")

loadings <- as.data.frame(chem_pca$rotation[, 1:2])
loadings$Element <- rownames(loadings)
write.csv(loadings, "chemical_pca_loadings.csv", row.names = FALSE)
cat("Saved: chemical_pca_loadings.csv\n")

cat("\n========================================\n")
cat("Chemical ordination complete.\n")
cat("========================================\n")
