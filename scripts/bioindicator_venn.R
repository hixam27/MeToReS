#!/usr/bin/env Rscript
# =============================================================================
# High-Confidence Bioindicator Discovery — Venn Diagrams
# =============================================================================
# Intersects genes independently flagged as significant by TWO different
# methods for the same chemical element:
#   - MaAsLin2 (linear models, country as random effect, season as fixed
#     covariate)
#   - ppcor partial correlations (Spearman, country + season as explicit
#     covariates)
#
# A gene flagged by BOTH methods is a substantially stronger bioindicator
# candidate than one flagged by either alone, since the two methods have
# different confounder-handling mechanisms (random effect vs. explicit
# covariate).
#
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(VennDiagram)
})

# =============================================================================
# Argument parsing
# =============================================================================
option_list <- list(
    make_option(c("--maaslin2"),     type="character",
        help="Path to maaslin2_significant.tsv [required]"),
    make_option(c("--partial_corr"), type="character",
        help="Path to partial_correlations_significant.csv [required]"),
    make_option(c("--elements"),     type="character", default="Cu Co Cd Ni Pb",
        help="Space-separated element symbols to process [default: %default]"),
    make_option(c("--out_dir"),      type="character",
        help="Output directory [required]")
)

opts <- parse_args(OptionParser(option_list=option_list))

for (arg in c("maaslin2", "partial_corr", "out_dir")) {
    if (is.null(opts[[arg]])) stop(paste("Missing required argument:", arg))
}
dir.create(opts$out_dir, recursive=TRUE, showWarnings=FALSE)

elements <- strsplit(trimws(opts$elements), "\\s+")[[1]]

cat("============================================\n")
cat("  High-Confidence Bioindicator Venn Diagrams\n")
cat("============================================\n")
cat("Elements:", paste(elements, collapse=", "), "\n\n")

# =============================================================================
# STEP 1: Load both significant-results tables
# =============================================================================
cat("=== STEP 1: Loading significant-results tables ===\n")

if (!file.exists(opts$maaslin2)) {
    stop(paste("MaAsLin2 significant-results file not found:", opts$maaslin2))
}
if (!file.exists(opts$partial_corr)) {
    stop(paste("Partial-correlation significant-results file not found:", opts$partial_corr))
}

maaslin_df <- read.delim(opts$maaslin2, stringsAsFactors=FALSE)
ppcor_df   <- read.csv(opts$partial_corr, stringsAsFactors=FALSE)

required_maaslin <- c("feature", "element", "coef", "qval")
missing_maaslin   <- setdiff(required_maaslin, colnames(maaslin_df))
if (length(missing_maaslin) > 0) {
    stop(paste0("maaslin2_significant.tsv missing expected column(s): ",
                paste(missing_maaslin, collapse=", "), "\n",
                "Available: ", paste(colnames(maaslin_df), collapse=", ")))
}

required_ppcor <- c("gene", "element", "estimate", "padj")
missing_ppcor  <- setdiff(required_ppcor, colnames(ppcor_df))
if (length(missing_ppcor) > 0) {
    stop(paste0("partial_correlations_significant.csv missing expected column(s): ",
                paste(missing_ppcor, collapse=", "), "\n",
                "Available: ", paste(colnames(ppcor_df), collapse=", ")))
}

cat("  MaAsLin2 rows           :", nrow(maaslin_df), "\n")
cat("  Partial correlation rows:", nrow(ppcor_df), "\n\n")

# =============================================================================
# STEP 2: Per-element intersection + Venn diagram
# =============================================================================
cat("=== STEP 2: Per-element intersection ===\n")

summary_rows <- list()

for (elem in elements) {

    cat("\n  --- Element:", elem, "---\n")

    maaslin_sub <- maaslin_df[maaslin_df$element == elem, , drop=FALSE]
    ppcor_sub   <- ppcor_df[ppcor_df$element == elem, , drop=FALSE]

    maaslin_genes <- unique(maaslin_sub$feature)
    ppcor_genes   <- unique(ppcor_sub$gene)

    intersect_genes <- intersect(maaslin_genes, ppcor_genes)

    n_maaslin_only <- length(setdiff(maaslin_genes, ppcor_genes))
    n_ppcor_only   <- length(setdiff(ppcor_genes, maaslin_genes))
    n_intersect    <- length(intersect_genes)

    cat("    MaAsLin2 significant genes            :", length(maaslin_genes), "\n")
    cat("    Partial correlation significant genes :", length(ppcor_genes), "\n")
    cat("    High-confidence (both methods)        :", n_intersect, "\n")

    # ── Save the intersecting gene list, with stats from BOTH methods ──
    out_csv <- file.path(opts$out_dir,
                         paste0("high_confidence_bioindicators_", elem, ".csv"))

    if (n_intersect > 0) {
        maaslin_stats <- maaslin_sub[maaslin_sub$feature %in% intersect_genes,
                                      c("feature", "coef", "qval")]
        colnames(maaslin_stats) <- c("GeneID", "maaslin2_coef", "maaslin2_qval")

        ppcor_stats <- ppcor_sub[ppcor_sub$gene %in% intersect_genes,
                                  c("gene", "estimate", "padj")]
        colnames(ppcor_stats) <- c("GeneID", "ppcor_estimate", "ppcor_padj")

        combined <- merge(maaslin_stats, ppcor_stats, by="GeneID")
        write.csv(combined, out_csv, row.names=FALSE)
    } else {
        write.csv(data.frame(GeneID=character(), maaslin2_coef=numeric(),
                              maaslin2_qval=numeric(), ppcor_estimate=numeric(),
                              ppcor_padj=numeric()),
                  out_csv, row.names=FALSE)
    }
    cat("    Saved:", basename(out_csv), "\n")

    # ── Venn diagram ─────────────────────────────────────────────────
    out_pdf <- file.path(opts$out_dir,
                         paste0("high_confidence_bioindicators_", elem, "_venn.pdf"))

    pdf(out_pdf, width=6, height=6)
    if (length(maaslin_genes) == 0 && length(ppcor_genes) == 0) {
        grid::grid.newpage()
        grid::grid.text(paste0(elem, ": no significant genes from either method"))
    } else {
        venn_plot <- draw.pairwise.venn(
            area1      = length(maaslin_genes),
            area2      = length(ppcor_genes),
            cross.area = n_intersect,
            category   = c("MaAsLin2", "Partial correlation"),
            fill       = c("#4C72B0", "#DD8452"),
            alpha      = 0.5,
            cat.pos    = c(-20, 20),
            cex        = 1.3,
            cat.cex    = 1.1,
            main       = elem,
            ind        = FALSE
        )
        grid::grid.draw(venn_plot)
    }
    dev.off()
    cat("    Saved:", basename(out_pdf), "\n")

    summary_rows[[elem]] <- data.frame(
        Element              = elem,
        MaAsLin2_significant = length(maaslin_genes),
        Ppcor_significant    = length(ppcor_genes),
        High_confidence      = n_intersect,
        MaAsLin2_only        = n_maaslin_only,
        Ppcor_only           = n_ppcor_only,
        stringsAsFactors     = FALSE
    )
}

# =============================================================================
# STEP 3: Cross-element summary
# =============================================================================
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df,
          file.path(opts$out_dir, "high_confidence_bioindicators_summary.csv"),
          row.names=FALSE)

cat("\n============================================\n")
cat("  Bioindicator intersection complete\n")
cat("============================================\n")
print(summary_df, row.names=FALSE)
cat("\nOutput directory:", opts$out_dir, "\n")
cat("============================================\n")
