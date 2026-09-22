#!/usr/bin/env Rscript
# resistance_deg_report.R
# (1) Annotate the per-contrast DESeq2 DEGs with ARG/MRG status (no new DESeq2),
# (2) draw per-ARO and per-compound expression heatmaps across groups.

suppressPackageStartupMessages({
  library(optparse)
  library(pheatmap)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--annotation",     type = "character"),
  make_option("--deg_dir",        type = "character"),
  make_option("--contrasts",      type = "character"),   # space-separated
  make_option("--by_aro",         type = "character"),
  make_option("--by_compound",    type = "character"),
  make_option("--samplesheet",    type = "character"),
  make_option("--group_col",      type = "character", default = "contaminant"),
  make_option("--out_table",      type = "character"),
  make_option("--out_arg_heatmap",type = "character"),
  make_option("--out_mrg_heatmap",type = "character"),
  make_option("--top_n",          type = "integer", default = 50L),
  make_option("--padj",           type = "double",  default = 0.05)
)))

contrasts <- strsplit(trimws(opt$contrasts), "\\s+")[[1]]
contrasts <- contrasts[nchar(contrasts) > 0]

## --- detect the gene-ID column in a DESeq2 results CSV (version-defensive) ---
read_deg <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  idcand <- c("ORF_id", "gene", "gene_id", "Gene", "X", "row", "Row.names")
  hit <- intersect(idcand, colnames(df))
  if (length(hit) >= 1) colnames(df)[colnames(df) == hit[1]] <- "ORF_id"
  else                  colnames(df)[1] <- "ORF_id"   # fall back to first column
  df
}

## --- 1. annotated resistant-DEG table across contrasts -----------------------
ann <- read.delim(opt$annotation, stringsAsFactors = FALSE, check.names = FALSE)

all_res <- list()
for (ct in contrasts) {
  f <- file.path(opt$deg_dir, ct, "differential_genes.csv")
  if (!file.exists(f)) { warning("missing DEG file: ", f); next }
  deg <- read_deg(f)
  m <- merge(deg, ann, by = "ORF_id")           # inner join = DE genes that are ARG/MRG
  if (nrow(m) == 0) next
  m$contrast <- ct
  all_res[[ct]] <- m
}

if (length(all_res)) {
  res <- do.call(rbind, all_res)
  front <- intersect(c("ORF_id", "contrast", "log2FoldChange", "padj", "pvalue",
                       "baseMean", "is_ARG", "Best_Hit_ARO", "Drug_Class",
                       "Resistance_Mechanism", "is_MRG", "Gene_name", "Compound"),
                     colnames(res))
  res <- res[, c(front, setdiff(colnames(res), front))]
  if ("padj" %in% colnames(res)) res <- res[order(res$contrast, res$padj), ]
  write.table(res, opt$out_table, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Resistant DEGs: ", nrow(res), " rows across ", length(all_res), " contrasts")
} else {
  writeLines("ORF_id\tcontrast\tlog2FoldChange\tpadj\tis_ARG\tis_MRG", opt$out_table)
  message("No resistant DEGs found in any contrast.")
}

## --- 2. column annotation (sample -> group) ---------------------------------
ss <- read.csv(opt$samplesheet, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("sample" %in% colnames(ss), opt$group_col %in% colnames(ss))
grp <- data.frame(group = ss[[opt$group_col]], row.names = ss$sample)

## --- helper: row z-scored heatmap, columns ordered/annotated by group -------
draw_heatmap <- function(path_in, path_out, title, top_n, ann_col) {
  if (!file.exists(path_in)) { warning("missing: ", path_in); return(invisible()) }
  m <- read.csv(path_in, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(m) == 0) { warning("empty: ", path_in); pdf(path_out); plot.new()
                      text(.5, .5, "no rows to plot"); dev.off(); return(invisible()) }
  rownames(m) <- m[[1]]; m[[1]] <- NULL
  m <- as.matrix(m)
  m <- m[rowSums(m) > 0, , drop = FALSE]
  if (nrow(m) > top_n)                                # keep most-expressed rows for readability
    m <- m[order(rowSums(m), decreasing = TRUE)[seq_len(top_n)], , drop = FALSE]
  common <- intersect(colnames(m), rownames(ann_col))
  m  <- m[, common, drop = FALSE]
  ac <- ann_col[common, , drop = FALSE]
  ord <- order(ac$group)                              # group samples together
  m <- m[, ord, drop = FALSE]; ac <- ac[ord, , drop = FALSE]
  z <- t(scale(t(m))); z[is.na(z)] <- 0               # row z-score => pattern across sites
  pheatmap(z, cluster_cols = FALSE, cluster_rows = nrow(z) > 1,
           annotation_col = ac, show_colnames = TRUE, fontsize_row = 7,
           main = title, filename = path_out,
           width = 9, height = max(4, nrow(z) * 0.18))
  message("Heatmap written: ", path_out, " (", nrow(z), " rows)")
}

draw_heatmap(opt$by_aro, opt$out_arg_heatmap,
             "ARG expression per ARO (row z-score of summed TPM)", opt$top_n, grp)
draw_heatmap(opt$by_compound, opt$out_mrg_heatmap,
             "Metal/biocide resistance per compound (row z-score)", opt$top_n, grp)
