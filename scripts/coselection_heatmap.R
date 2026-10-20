#!/usr/bin/env Rscript
# =============================================================================
# Co-selection Heatmap — Top ARG/MRG families vs. chemical elements
# =============================================================================
# Clustered heatmap of Spearman correlations (from full_correlation_matrix.csv)
# between the top N most abundant ARG/MRG families and the chemical elements.
#
# --element_selection decides which elements are shown:
#   all     every element (default): a descriptive figure that hides nothing
#           behind a test
#   padj    elements with envfit BH-adjusted p < --envfit_alpha
#   pvalue  elements with envfit raw p < --envfit_alpha (exploratory)
# Cells are marked "*" where the correlation's p-value (--annotate_on padj,
# the default, or pvalue) is below --corr_alpha.
#
#   full_correlation_matrix.csv        : Source, Category, Element, n, rho,
#                                         pvalue, padj
#   envfit_results.csv                 : Element, NMDS1, NMDS2, r2, pvalue
#   arg_abundance_by_aro.csv           : Best_Hit_ARO + one column per sample
#   mrg_abundance_by_compound.csv      : Compound       + one column per sample
#
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(tidyr)
    library(pheatmap)
})

# =============================================================================
# Argument parsing
# =============================================================================
option_list <- list(
    make_option(c("--full_corr"),      type="character",
        help="Path to full_correlation_matrix.csv [required]"),
    make_option(c("--envfit"),         type="character",
        help="Path to envfit_results.csv [required]"),
    make_option(c("--arg_abundance"),  type="character",
        help="Path to arg_abundance_by_aro.csv [required]"),
    make_option(c("--mrg_abundance"),  type="character",
        help="Path to mrg_abundance_by_compound.csv [required]"),
    make_option(c("--top_n"),          type="integer", default=20,
        help="Number of most abundant ARG/MRG families to show [default: %default]"),
    make_option(c("--element_selection"), type="character", default="all",
        help="all | padj | pvalue [default: %default]"),
    make_option(c("--envfit_alpha"),   type="double", default=0.05,
        help="envfit cut-off when --element_selection is padj or pvalue [default: %default]"),
    make_option(c("--corr_alpha"),     type="double", default=0.05,
        help="cut-off for asterisk annotation [default: %default]"),
    make_option(c("--annotate_on"),    type="character", default="padj",
        help="correlation column used for asterisks: padj | pvalue [default: %default]"),
    make_option(c("--out_dir"),        type="character",
        help="Output directory [required]")
)

opts <- parse_args(OptionParser(option_list=option_list))

for (arg in c("full_corr", "envfit", "arg_abundance", "mrg_abundance", "out_dir")) {
    if (is.null(opts[[arg]])) stop(paste("Missing required argument:", arg))
}
dir.create(opts$out_dir, recursive=TRUE, showWarnings=FALSE)

cat("============================================\n")
cat("  Co-selection Heatmap\n")
cat("============================================\n")

# =============================================================================
# STEP 1: Elements to show
# =============================================================================
sel <- opts$element_selection
if (!sel %in% c("all", "padj", "pvalue")) stop("--element_selection must be all, padj or pvalue")
if (!opts$annotate_on %in% c("padj", "pvalue")) stop("--annotate_on must be padj or pvalue")
cat("=== STEP 1: Elements (selection:", sel, ") ===\n")

envfit_df <- read.csv(opts$envfit, stringsAsFactors=FALSE)
required_envfit <- c("Element", "pvalue")
missing_envfit   <- setdiff(required_envfit, colnames(envfit_df))
if (length(missing_envfit) > 0) {
    stop(paste0("envfit_results.csv missing expected column(s): ",
                paste(missing_envfit, collapse=", "), "\n",
                "Available: ", paste(colnames(envfit_df), collapse=", ")))
}

if (!"padj" %in% colnames(envfit_df))
    envfit_df$padj <- p.adjust(envfit_df$pvalue, method = "BH")
sig_elements <- switch(sel,
    all    = envfit_df$Element,
    padj   = envfit_df$Element[envfit_df$padj   < opts$envfit_alpha],
    pvalue = envfit_df$Element[envfit_df$pvalue < opts$envfit_alpha])
sel_label <- switch(sel,
    all    = "all elements",
    padj   = paste0("elements with envfit q < ", opts$envfit_alpha),
    pvalue = paste0("elements with envfit p < ", opts$envfit_alpha, " (unadjusted, exploratory)"))

if (length(sig_elements) == 0) {
    cat("No elements selected (", sel_label, ") — saving a placeholder instead of a",
        "heatmap. This is an honest null result at this sample size, not a failure.\n")

    pdf(file.path(opts$out_dir, "coselection_heatmap.pdf"), width = 6, height = 4)
    grid::grid.newpage()
    grid::grid.text(paste0("No elements selected (", sel_label,
                            ")\nNo co-selection heatmap to display."))
    dev.off()

    write.csv(data.frame(Source=character(), Category=character(), Element=character(),
                          n=numeric(), rho=numeric(), pvalue=numeric(), padj=numeric()),
              file.path(opts$out_dir, "coselection_heatmap_data.csv"), row.names=FALSE)

    quit(status = 0)
}

cat("  Elements shown (", length(sig_elements), "):", paste(sig_elements, collapse=", "), "\n\n", sep="")

# =============================================================================
# STEP 2: Determine top N most abundant ARG/MRG families (combined ranking)
# =============================================================================
cat("=== STEP 2: Top", opts$top_n, "most abundant ARG/MRG families ===\n")

arg_df <- read.csv(opts$arg_abundance, stringsAsFactors=FALSE, check.names=FALSE)
mrg_df <- read.csv(opts$mrg_abundance, stringsAsFactors=FALSE, check.names=FALSE)

arg_id_col <- "Best_Hit_ARO"
mrg_id_col <- "Compound"
if (!arg_id_col %in% colnames(arg_df)) {
    stop(paste0("arg_abundance_by_aro.csv missing expected column '", arg_id_col, "'"))
}
if (!mrg_id_col %in% colnames(mrg_df)) {
    stop(paste0("mrg_abundance_by_compound.csv missing expected column '", mrg_id_col, "'"))
}

# Sample columns = numeric columns other than the ID (annotation columns
# such as gene family names are ignored).
num_cols <- function(df, id) {
    cols <- setdiff(colnames(df), id)
    cols[vapply(df[cols], is.numeric, logical(1))]
}
sample_cols_arg <- num_cols(arg_df, arg_id_col)
sample_cols_mrg <- num_cols(mrg_df, mrg_id_col)

arg_totals <- data.frame(
    Category    = arg_df[[arg_id_col]],
    Source      = "ARG",
    total_abund = rowSums(arg_df[, sample_cols_arg, drop=FALSE], na.rm=TRUE),
    stringsAsFactors = FALSE
)
mrg_totals <- data.frame(
    Category    = mrg_df[[mrg_id_col]],
    Source      = "MRG",
    total_abund = rowSums(mrg_df[, sample_cols_mrg, drop=FALSE], na.rm=TRUE),
    stringsAsFactors = FALSE
)

all_totals <- rbind(arg_totals, mrg_totals)

# Drop degenerate/empty Category entries (e.g. blank strings, stray "2")
# before ranking, so they can never occupy a top-N slot or appear as a
# heatmap row label.
n_before <- nrow(all_totals)
all_totals <- all_totals[!is.na(all_totals$Category) &
                          trimws(all_totals$Category) != "" &
                          nchar(trimws(all_totals$Category)) > 1, ]
cat("  Dropped", n_before - nrow(all_totals), "degenerate/empty Category entries\n")

top_categories_df <- all_totals %>%
    arrange(desc(total_abund)) %>%
    distinct(Category, .keep_all=TRUE) %>%
    slice_head(n=opts$top_n)

cat("  Top", nrow(top_categories_df), "families (", 
    sum(top_categories_df$Source=="ARG"), "ARG,",
    sum(top_categories_df$Source=="MRG"), "MRG):\n")
for (i in seq_len(nrow(top_categories_df))) {
    cat("    -", top_categories_df$Source[i], ":", top_categories_df$Category[i],
        "(total abundance:", round(top_categories_df$total_abund[i], 1), ")\n")
}
cat("\n")

# =============================================================================
# STEP 3: Load correlations, restrict to top families x significant elements
# =============================================================================
cat("=== STEP 3: Loading correlations ===\n")

corr_df <- read.csv(opts$full_corr, stringsAsFactors=FALSE)
required_corr <- c("Source", "Category", "Element", "rho", "pvalue")
missing_corr  <- setdiff(required_corr, colnames(corr_df))
if (length(missing_corr) > 0) {
    stop(paste0("full_correlation_matrix.csv missing expected column(s): ",
                paste(missing_corr, collapse=", "), "\n",
                "Available: ", paste(colnames(corr_df), collapse=", ")))
}

# Asterisks use the BH-adjusted p-value over the whole correlation table
# (all families x all elements), not only the cells shown.
if (!"padj" %in% colnames(corr_df))
    corr_df$padj <- p.adjust(corr_df$pvalue, method = "BH")
corr_df$p_mark <- corr_df[[opts$annotate_on]]

corr_sub <- corr_df[corr_df$Category %in% top_categories_df$Category &
                     corr_df$Element  %in% sig_elements, ,
                     drop=FALSE]
cat("  Correlation rows matching top families x significant elements:", nrow(corr_sub), "\n")

if (nrow(corr_sub) == 0) {
    cat("No correlation rows match the top families x significant elements —",
        "saving a placeholder instead of a heatmap.\n")

    pdf(file.path(opts$out_dir, "coselection_heatmap.pdf"), width = 6, height = 4)
    grid::grid.newpage()
    grid::grid.text("No correlation data available for the top families\nagainst the significant elements.")
    dev.off()

    write.csv(data.frame(Source=character(), Category=character(), Element=character(),
                          n=numeric(), rho=numeric(), pvalue=numeric(), padj=numeric()),
              file.path(opts$out_dir, "coselection_heatmap_data.csv"), row.names=FALSE)

    quit(status = 0)
}

# =============================================================================
# STEP 4: Build rho matrix + p-value matrix for asterisk annotation
# =============================================================================
cat("\n=== STEP 4: Building heatmap matrices ===\n")

rho_wide <- corr_sub %>%
    select(Category, Element, rho) %>%
    distinct(Category, Element, .keep_all=TRUE) %>%
    tidyr::pivot_wider(names_from=Element, values_from=rho) %>%
    as.data.frame()
rownames(rho_wide) <- rho_wide$Category
rho_wide$Category <- NULL

pval_wide <- corr_sub %>%
    select(Category, Element, p_mark) %>%
    distinct(Category, Element, .keep_all=TRUE) %>%
    tidyr::pivot_wider(names_from=Element, values_from=p_mark) %>%
    as.data.frame()
rownames(pval_wide) <- pval_wide$Category
pval_wide$Category <- NULL

# Ensure both matrices have identical row/column ordering before proceeding
common_rows <- intersect(rownames(rho_wide), rownames(pval_wide))
common_cols <- intersect(colnames(rho_wide), colnames(pval_wide))
rho_mat  <- as.matrix(rho_wide[common_rows, common_cols, drop=FALSE])
pval_mat <- as.matrix(pval_wide[common_rows, common_cols, drop=FALSE])

cat("  Heatmap dimensions:", nrow(rho_mat), "families x", ncol(rho_mat), "elements\n")
cat("  Missing (NA) cells:", sum(is.na(rho_mat)), "of", length(rho_mat),
    "(no correlation was computed for that family x element pair)\n\n")

# Asterisk annotation: "*" where the chosen p-value < corr_alpha, blank otherwise.
# NA cells (no computed correlation) get a blank, not an asterisk.
display_mat <- matrix("", nrow=nrow(rho_mat), ncol=ncol(rho_mat),
                       dimnames=dimnames(rho_mat))
display_mat[!is.na(pval_mat) & pval_mat < opts$corr_alpha] <- "*"

# Save the underlying data (long format, easy to inspect/reuse)
out_data_csv <- file.path(opts$out_dir, "coselection_heatmap_data.csv")
write.csv(corr_sub %>% select(-p_mark) %>% arrange(Source, Category, Element),
          out_data_csv, row.names=FALSE)
cat("Saved:", basename(out_data_csv), "\n")

# =============================================================================
# STEP 5: Clustered heatmap
# =============================================================================
cat("\n=== STEP 5: Rendering heatmap ===\n")

# pheatmap's clustering can't handle NA in the distance calculation, so
# clustering is computed on an NA->0 imputed copy (treating "no computed
# correlation" as "no relationship" for clustering purposes only), while
# the DISPLAYED matrix keeps true NA (rendered as a distinct grey via
# na_col) so missing data is never visually confused with a real rho of 0.
rho_for_clustering <- rho_mat
rho_for_clustering[is.na(rho_for_clustering)] <- 0

# Clustering needs at least 2 rows (columns); with a single significant
# element or family that dimension is shown unclustered instead of failing.
cluster_r <- nrow(rho_for_clustering) >= 2
cluster_c <- ncol(rho_for_clustering) >= 2
row_dist <- if (cluster_r) dist(rho_for_clustering) else "euclidean"
col_dist <- if (cluster_c) dist(t(rho_for_clustering)) else "euclidean"
if (!cluster_r) cat("  Only 1 family: rows not clustered\n")
if (!cluster_c) cat("  Only 1 element: columns not clustered\n")

# Row/column annotation: which Source (ARG/MRG) each family belongs to
row_annotation <- data.frame(
    Source = top_categories_df$Source[match(rownames(rho_mat), top_categories_df$Category)],
    row.names = rownames(rho_mat)
)

out_pdf <- file.path(opts$out_dir, "coselection_heatmap.pdf")
pheatmap(
    rho_mat,
    display_numbers        = display_mat,
    number_color           = "black",
    fontsize_number         = 10,
    cluster_rows            = cluster_r,
    cluster_cols            = cluster_c,
    clustering_distance_rows = row_dist,
    clustering_distance_cols = col_dist,
    clustering_method       = "average",
    color                   = colorRampPalette(c("#2166ac", "white", "#b2182b"))(101),
    breaks                  = seq(-1, 1, length.out=102),
    na_col                  = "grey85",
    annotation_row          = row_annotation,
    annotation_colors       = list(Source=c(ARG="#4C72B0", MRG="#DD8452")),
    main                    = paste0("Co-selection: top ", nrow(rho_mat),
                                     " ARG/MRG families vs ", sel_label, "\n",
                                     "(* = ", ifelse(opts$annotate_on == "padj", "BH q", "p, unadjusted"),
                                     " < ", opts$corr_alpha, ")"),
    fontsize                = 10,
    fontsize_row            = 8,
    filename                = out_pdf,
    width                   = 10,
    height                  = max(6, 0.3 * nrow(rho_mat))
)
cat("Saved:", basename(out_pdf), "\n")

cat("\n============================================\n")
cat("  Co-selection heatmap complete\n")
cat("============================================\n")
cat("Output directory:", opts$out_dir, "\n")
cat("============================================\n")
