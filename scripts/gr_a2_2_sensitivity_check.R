#!/usr/bin/env Rscript
# gr_a2_2_sensitivity_check.R
#
# Two checks on whether the 186 Greece_contaminant_vs_Greece_Natural DEGs
# depend on GR_A_2_2 specifically (size factor 4.37, the third-most-extreme
# in the dataset, sitting on the non-reference side of this contrast):
#
#   (1) Composition-robust re-normalization: refit DESeq2 with the
#       poscounts size-factor estimator (uses only positive counts per
#       gene for the ratio, standard alternative for datasets with heavy
#       zero-inflation — this dataset ranges 44.9-93% zeros by sample) and
#       compare the resulting DEG list against the original.
#   (2) Per-DEG consistency: for the original 186 DEGs, does the signal
#       hold in BOTH Greece_contaminant samples (GR_A_2_2 and GR_S_2_3),
#       or is it driven by GR_A_2_2 alone? Avoids refitting at n=1, which
#       would repeat the known-flawed methodology from the earlier
#       exclusion experiment.

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--dds",            type = "character", help = "dds_object.RData"),
  make_option("--original_degs",  type = "character", help = "Greece_contaminant_vs_Greece_Natural/differential_genes.csv"),
  make_option("--out_dir",        type = "character"),
  make_option("--pvalue",         type = "double", default = 0.05),
  make_option("--fold_change",    type = "double", default = 2)
)))

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========================================\n")
cat("GR_A_2_2 sensitivity check\n")
cat("========================================\n\n")

load(opt$dds)  # loads `dds`

# ── Part 1: composition-robust re-normalization (poscounts) ────────────
cat("=== Part 1: poscounts re-normalization ===\n\n")

cat("Original size factors:\n")
print(round(sizeFactors(dds), 3))

dds_pc <- estimateSizeFactors(dds, type = "poscounts")
cat("\nPoscounts size factors:\n")
print(round(sizeFactors(dds_pc), 3))

cat("\nGR_A_2_2 specifically:\n")
cat("  Original (median-of-ratios):", round(sizeFactors(dds)["GR_A_2_2"], 3), "\n")
cat("  Poscounts:                  ", round(sizeFactors(dds_pc)["GR_A_2_2"], 3), "\n")

cat("\nRefitting dispersions and GLM under poscounts normalization...\n")
dds_pc <- DESeq(dds_pc, quiet = TRUE)

res_pc <- results(dds_pc, contrast = c("group", "Greece_contaminant", "Greece_Natural"))
res_pc_df <- as.data.frame(res_pc)
res_pc_df$ORF_id <- rownames(res_pc_df)

log2fc_thresh <- log2(opt$fold_change)
degs_pc <- res_pc_df[!is.na(res_pc_df$padj) &
                      res_pc_df$padj < opt$pvalue &
                      abs(res_pc_df$log2FoldChange) > log2fc_thresh, ]

write.csv(degs_pc, file.path(opt$out_dir, "greece_contaminant_DEGs_poscounts.csv"), row.names = FALSE)

cat("\nOriginal DEGs (median-of-ratios): reading from", opt$original_degs, "\n")
degs_orig <- read.csv(opt$original_degs, stringsAsFactors = FALSE)
orig_ids <- degs_orig[[1]]

cat("  Original DEG count :", length(orig_ids), "\n")
cat("  Poscounts DEG count:", nrow(degs_pc), "\n")

overlap <- intersect(orig_ids, degs_pc$ORF_id)
cat("  Overlap            :", length(overlap), "\n")
cat("  Jaccard similarity :", round(length(overlap) / length(union(orig_ids, degs_pc$ORF_id)), 3), "\n")
cat("  % of original DEGs that still replicate under poscounts:",
    round(100 * length(overlap) / length(orig_ids), 1), "%\n\n")

# ── Part 2: per-DEG consistency between the two Greece_contaminant samples ──
cat("=== Part 2: per-DEG consistency (GR_A_2_2 vs GR_S_2_3) ===\n\n")

norm_counts <- counts(dds, normalized = TRUE)
gr_a2_2 <- norm_counts[, "GR_A_2_2"]
gr_s2_3 <- norm_counts[, "GR_S_2_3"]
gr_natural_mean <- rowMeans(norm_counts[, c("GR_S_1_3", "GR_A_1_3")])

consistency <- data.frame(
  ORF_id = orig_ids,
  GR_A_2_2_norm = gr_a2_2[orig_ids],
  GR_S_2_3_norm = gr_s2_3[orig_ids],
  Greece_Natural_mean_norm = gr_natural_mean[orig_ids]
)
# Does each Greece_contaminant sample individually move in the SAME
# direction away from the Greece_Natural mean as the reported DEG result?
consistency$A2_2_direction <- ifelse(consistency$GR_A_2_2_norm > consistency$Greece_Natural_mean_norm, "up", "down")
consistency$S2_3_direction <- ifelse(consistency$GR_S_2_3_norm > consistency$Greece_Natural_mean_norm, "up", "down")
consistency$both_agree <- consistency$A2_2_direction == consistency$S2_3_direction

write.csv(consistency, file.path(opt$out_dir, "greece_contaminant_per_DEG_consistency.csv"), row.names = FALSE)

n_agree <- sum(consistency$both_agree)
cat("DEGs where GR_A_2_2 and GR_S_2_3 move in the SAME direction relative to Greece_Natural:",
    n_agree, "of", nrow(consistency),
    sprintf("(%.1f%%)\n", 100 * n_agree / nrow(consistency)))
cat("DEGs where they DISAGREE (signal likely driven by one sample only):",
    nrow(consistency) - n_agree, "of", nrow(consistency), "\n\n")

cat("========================================\n")
cat("Sensitivity check complete. See:\n")
cat(" ", file.path(opt$out_dir, "greece_contaminant_DEGs_poscounts.csv"), "\n")
cat(" ", file.path(opt$out_dir, "greece_contaminant_per_DEG_consistency.csv"), "\n")
cat("========================================\n")
