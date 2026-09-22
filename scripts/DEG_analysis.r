#!/usr/bin/env Rscript
# =============================================================================
# DESeq2 Differential Expression Analysis for MetoReS
# Multi-contrast version: all levels vs fixed reference level
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(tximport)
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(tibble)
})

# =============================================================================
# Argument parsing
# =============================================================================
option_list <- list(
    make_option(c("-q", "--quant_dir"),       type="character"),
    make_option(c("-i", "--filtered_matrix"), type="character"),
    make_option(c("-s", "--samplesheet"),     type="character"),
    make_option(c("-o", "--output_dir"),      type="character"),
    make_option(c("--group_col"),             type="character", default="group"),
    make_option(c("--formula"),               type="character", default="~ group"),
    make_option(c("-p", "--pvalue"),          type="double",    default=0.05),
    make_option(c("-f", "--fold_change"),     type="double",    default=2.0),
    make_option(c("--ref_level"),             type="character", default=""),
    make_option(c("--ref_level_filter"),      type="character", default="",
        help="Optional 'column:value' filter to disambiguate which samples
              count as the reference level, when the reference level's name
              alone doesn't uniquely identify the intended reference samples
              (e.g. two samples share group_col=='Greece' but only one
              contaminant condition should count as the reference).
              Format: 'contaminant:Natural'. Leave empty (default) when
              the reference level is otherwise unambiguous."),
    make_option(c("--min_count"),             type="integer",   default=10),
    make_option(c("-n", "--comparison_name"), type="character", default="")
)

opts <- parse_args(OptionParser(option_list=option_list))
log2fc_cutoff <- log2(opts$fold_change)

# ref_level is now REQUIRED
if (opts$ref_level == "") {
    stop("--ref_level is required. Set it in config.yaml")
}

cat("============================================\n")
cat("  DESeq2 Multi-Contrast Analysis\n")
cat("============================================\n")
cat("Formula          :", opts$formula, "\n")
cat("Reference level  :", opts$ref_level, "\n")
cat("P-value cutoff   :", opts$pvalue, "\n")
cat("Fold change      :", opts$fold_change,
    "(log2FC:", round(log2fc_cutoff, 3), ")\n\n")

# =============================================================================
# Validate inputs
# =============================================================================
for (arg in c("quant_dir", "filtered_matrix", "samplesheet", "output_dir")) {
    if (is.null(opts[[arg]])) stop(paste("Missing required argument:", arg))
}
dir.create(opts$output_dir, recursive=TRUE, showWarnings=FALSE)

# =============================================================================
# STEP 1: Load samplesheet
# =============================================================================
cat("=== STEP 1: Loading samplesheet ===\n")
sample_info <- read.csv(opts$samplesheet, stringsAsFactors=FALSE)

required_cols <- c("sample", "sra")
missing_cols  <- setdiff(required_cols, colnames(sample_info))
if (length(missing_cols) > 0) {
    stop(paste("Samplesheet missing columns:", paste(missing_cols, collapse=", ")))
}

# Use the column specified by --group_col from Snakefile
if (opts$group_col != "group") {
    if (!opts$group_col %in% colnames(sample_info)) {
        stop(paste0("Column '", opts$group_col, "' not found in samplesheet. ",
                    "Available: ", paste(colnames(sample_info), collapse=", ")))
    }
    colnames(sample_info)[colnames(sample_info) == opts$group_col] <- "group"
    cat("  Using column '", opts$group_col, "' as grouping variable\n")
}

cat("  Samples:", nrow(sample_info), "\n")
cat("  Groups :", paste(unique(sample_info$group), collapse=", "), "\n\n")

# =============================================================================
# STEP 1b: Optional ref_level disambiguation
# =============================================================================
# If multiple samples share the group_col value used as --ref_level, but
# only some of them should actually count as "the" reference (e.g. two
# samples are both group_col=="Greece" but represent different contaminant
# conditions), --ref_level_filter identifies which ones are the REAL
# reference and excludes the rest from this analysis entirely. This is a
# general mechanism, not specific to any one dataset: it only activates
# when --ref_level_filter is explicitly set, and has no effect otherwise.
if (opts$ref_level_filter != "") {
    filter_parts <- strsplit(opts$ref_level_filter, ":", fixed=TRUE)[[1]]
    if (length(filter_parts) != 2) {
        stop(paste0("--ref_level_filter must be in 'column:value' format, got: '",
                    opts$ref_level_filter, "'"))
    }
    filter_col   <- trimws(filter_parts[1])
    filter_value <- trimws(filter_parts[2])

    if (!filter_col %in% colnames(sample_info)) {
        stop(paste0("--ref_level_filter column '", filter_col,
                    "' not found in samplesheet. Available: ",
                    paste(colnames(sample_info), collapse=", ")))
    }

    # A sample is excluded only if BOTH: it belongs to the reference group_col
    # value, AND it fails the disambiguating filter. Samples in other groups
    # are never touched by this logic.
    is_ambiguous_ref <- (sample_info$group == opts$ref_level) &
                        (sample_info[[filter_col]] != filter_value)

    if (any(is_ambiguous_ref)) {
        excluded_samples <- sample_info$sample[is_ambiguous_ref]
        cat("  ref_level_filter '", opts$ref_level_filter, "' applied:\n", sep="")
        cat("    Excluding", length(excluded_samples),
            "sample(s) that share group_col =='", opts$ref_level,
            "' but do not match ", filter_col, "=='", filter_value, "':\n", sep="")
        for (s in excluded_samples) cat("      -", s, "\n")

        sample_info <- sample_info[!is_ambiguous_ref, , drop=FALSE]
        cat("  Samples after ref_level_filter:", nrow(sample_info), "\n\n")
    } else {
        cat("  ref_level_filter '", opts$ref_level_filter,
            "' set but no ambiguous samples found — no samples excluded.\n\n", sep="")
    }
}

# =============================================================================
# STEP 2: Load filtered gene list
# =============================================================================
cat("=== STEP 2: Loading filtered gene list ===\n")
filtered_matrix <- read.csv(opts$filtered_matrix, stringsAsFactors=FALSE)
if (!"GeneID" %in% colnames(filtered_matrix)) {
    stop("Filtered matrix must have a 'GeneID' column.")
}
filtered_genes <- filtered_matrix$GeneID
cat("  Filtered genes:", length(filtered_genes), "\n\n")

# =============================================================================
# STEP 3: Import Salmon quantifications via tximport
# =============================================================================
cat("=== STEP 3: Importing Salmon counts ===\n")
quant_files <- file.path(opts$quant_dir, paste0(sample_info$sra, "_quant"), "quant.sf")
names(quant_files) <- sample_info$sample

missing_files <- quant_files[!file.exists(quant_files)]
if (length(missing_files) > 0) {
    stop(paste("Missing quant files:\n ", paste(missing_files, collapse="\n  ")))
}

txi <- tximport(quant_files, type="salmon",
                countsFromAbundance="lengthScaledTPM", txOut=TRUE)
cat("  Total genes in tximport:", nrow(txi$counts), "\n")

# =============================================================================
# STEP 4: Subset to filtered genes
# =============================================================================
cat("\n=== STEP 4: Subsetting to filtered genes ===\n")
genes_in_counts <- intersect(filtered_genes, rownames(txi$counts))
cat("  Filtered genes present:", length(genes_in_counts), "\n")

txi_filtered <- list(
    counts              = txi$counts[genes_in_counts, , drop=FALSE],
    abundance           = txi$abundance[genes_in_counts, , drop=FALSE],
    length              = txi$length[genes_in_counts, , drop=FALSE],
    countsFromAbundance = txi$countsFromAbundance
)

# =============================================================================
# STEP 5: Build colData and set reference level
# =============================================================================
cat("\n=== STEP 5: Building colData ===\n")
rownames(sample_info) <- sample_info$sample

formula_obj  <- as.formula(opts$formula)
formula_vars <- all.vars(formula_obj)

for (v in formula_vars) {
    if (!v %in% colnames(sample_info)) {
        stop(paste0("Variable '", v, "' in formula not found. Available: ",
                    paste(colnames(sample_info), collapse=", ")))
    }
    sample_info[[v]] <- factor(sample_info[[v]])
}

main_factor <- formula_vars[length(formula_vars)]
all_levels  <- levels(sample_info[[main_factor]])

cat("  Main factor:", main_factor, "\n")
cat("  Levels     :", paste(all_levels, collapse=", "), "\n")

# Validate reference level
if (!opts$ref_level %in% all_levels) {
    stop(paste0("Reference level '", opts$ref_level, "' not found. Available: ",
                paste(all_levels, collapse=", ")))
}

sample_info[[main_factor]] <- relevel(sample_info[[main_factor]], ref=opts$ref_level)
cat("  Reference  :", opts$ref_level, "\n")

non_ref_levels <- setdiff(levels(sample_info[[main_factor]]), opts$ref_level)
cat("  Contrasts to extract:", length(non_ref_levels), "\n")
for (lvl in non_ref_levels) {
    cat("    -", lvl, "vs", opts$ref_level, "\n")
}

col_data <- sample_info[colnames(txi_filtered$counts), , drop=FALSE]

# =============================================================================
# STEP 6: Create DESeqDataSet and run DESeq2 (once)
# =============================================================================
cat("\n=== STEP 6: Running DESeq2 (single fit) ===\n")
dds <- DESeqDataSetFromTximport(txi=txi_filtered, colData=col_data, design=formula_obj)

cat("  Genes before pre-filter:", nrow(dds), "\n")
keep <- rowSums(counts(dds)) >= opts$min_count
dds  <- dds[keep, ]
cat("  Genes after pre-filter :", nrow(dds), "\n")

dds <- DESeq(dds)
cat("  DESeq2 fitted. Available coefficients:\n")
for (rn in resultsNames(dds)) cat("    -", rn, "\n")

cat("\n=== STEP 7: Shared outputs (PCA + VST) ===\n")

vst_mat <- vst(dds, blind=FALSE)
vst_df  <- as.data.frame(assay(vst_mat)) %>% rownames_to_column("GeneID")
write.csv(vst_df, file.path(opts$output_dir, "vst_normalised_counts.csv"), row.names=FALSE)
cat("  Saved: vst_normalised_counts.csv\n")

save(dds, file=file.path(opts$output_dir, "dds_object.RData"))
cat("  Saved: dds_object.RData\n")

norm_counts <- counts(dds, normalized=TRUE)
write.csv(norm_counts, file.path(opts$output_dir, "deseq2_normalized_counts.csv"))
cat("  Saved: deseq2_normalized_counts.csv\n")

pca_data <- plotPCA(vst_mat, intgroup=main_factor, returnData=TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x=PC1, y=PC2, color=.data[[main_factor]])) +
    geom_point(size=4) +
    geom_text_repel(aes(label=name), size=3, max.overlaps=20) +
    xlab(paste0("PC1: ", pct_var[1], "% variance")) +
    ylab(paste0("PC2: ", pct_var[2], "% variance")) +
    theme_bw(base_size=14) +
    theme(plot.title=element_text(face="bold", hjust=0.5)) +
    ggtitle(paste0("PCA — all levels (ref: ", opts$ref_level, ")"))

ggsave(file.path(opts$output_dir, "samples_PCA.pdf"), p_pca,
       width=8, height=6, device=cairo_pdf)
cat("  Saved: samples_PCA.pdf\n")

pdf(file.path(opts$output_dir, "deseq2_dispersion.pdf"), width=8, height=6)
plotDispEsts(dds, main="Dispersion Estimates")
dev.off()
cat("  Saved: deseq2_dispersion.pdf\n\n")

# =============================================================================
# STEP 8: Extract each contrast vs reference level
# =============================================================================
cat("=== STEP 8: Extracting per-contrast results ===\n")

all_contrasts_summary <- data.frame(
    Contrast      = character(),
    Genes_tested  = integer(),
    Total_DEGs    = integer(),
    Upregulated   = integer(),
    Downregulated = integer(),
    stringsAsFactors = FALSE
)

for (lvl in non_ref_levels) {
    contrast_name <- paste0(lvl, "_vs_", opts$ref_level)
    contrast_dir  <- file.path(opts$output_dir, contrast_name)
    dir.create(contrast_dir, recursive=TRUE, showWarnings=FALSE)

    cat("\n  --- Contrast:", contrast_name, "---\n")

    res <- results(dds, contrast=c(main_factor, lvl, opts$ref_level),
                   alpha=opts$pvalue)

    n_extreme_before <- sum(abs(res$log2FoldChange) > 15, na.rm=TRUE)
    res <- lfcShrink(dds, contrast=c(main_factor, lvl, opts$ref_level),
                      type="normal", res=res)
    n_extreme_after  <- sum(abs(res$log2FoldChange) > 15, na.rm=TRUE)
    cat("  log2FC shrinkage applied (type=normal): genes with |log2FC|>15",
        "went from", n_extreme_before, "to", n_extreme_after, "\n")

    cat("  DESeq2 summary for ", contrast_name, ":\n")
    summary(res)

    res_df <- as.data.frame(res) %>%
        rownames_to_column("GeneID") %>%
        arrange(padj)

    sig_df <- res_df %>%
        filter(!is.na(padj), padj < opts$pvalue,
               abs(log2FoldChange) > log2fc_cutoff) %>%
        mutate(change = case_when(
            log2FoldChange >  log2fc_cutoff ~ "Up",
            log2FoldChange < -log2fc_cutoff ~ "Down",
            TRUE ~ "NS"
        ))

    up_genes   <- sig_df %>% filter(change == "Up")
    down_genes <- sig_df %>% filter(change == "Down")

    cat("    Total DEGs   :", nrow(sig_df), "\n")
    cat("    Upregulated  :", nrow(up_genes), "\n")
    cat("    Downregulated:", nrow(down_genes), "\n")

    write.csv(res_df,
              file.path(contrast_dir, "deseq2_results_all.csv"),
              row.names=FALSE)

    write.csv(
        sig_df %>% select(GeneID, baseMean, log2FoldChange, lfcSE, pvalue, padj, change),
        file.path(contrast_dir, "differential_genes.csv"),
        row.names=FALSE
    )

    writeLines(sig_df$GeneID,
               file.path(contrast_dir, "differential_genes_id.txt"))
    writeLines(up_genes$GeneID,
               file.path(contrast_dir, "upregulated_genes_id.txt"))
    writeLines(down_genes$GeneID,
               file.path(contrast_dir, "downregulated_genes_id.txt"))

    vol_df <- res_df %>%
        filter(!is.na(padj), !is.na(log2FoldChange)) %>%
        mutate(
            neg_log10_padj = -log10(padj),
            significance   = case_when(
                padj < opts$pvalue & log2FoldChange >  log2fc_cutoff ~ "Up",
                padj < opts$pvalue & log2FoldChange < -log2fc_cutoff ~ "Down",
                TRUE ~ "NS"
            ),
            significance = factor(significance, levels=c("Down", "NS", "Up"))
        )

    colors <- c("Down"="#2166ac", "NS"="grey60", "Up"="#b2182b")

    p_vol <- ggplot(vol_df, aes(x=log2FoldChange, y=neg_log10_padj, color=significance)) +
        geom_point(alpha=0.5, size=1.2) +
        scale_color_manual(values=colors) +
        geom_vline(xintercept=c(-log2fc_cutoff, log2fc_cutoff),
                   linetype="dashed", color="grey30") +
        geom_hline(yintercept=-log10(opts$pvalue),
                   linetype="dashed", color="grey30") +
        labs(
            title    = paste0("Volcano: ", contrast_name),
            subtitle = paste0("padj < ", opts$pvalue,
                              " | FC > ", opts$fold_change),
            x        = expression(log[2]~Fold~Change),
            y        = expression(-log[10]~(padj))
        ) +
        annotate("text", x=Inf,  y=Inf, label=paste("Up:", nrow(up_genes)),
                 hjust=1.1, vjust=1.5, color="#b2182b", fontface="bold") +
        annotate("text", x=-Inf, y=Inf, label=paste("Down:", nrow(down_genes)),
                 hjust=-0.1, vjust=1.5, color="#2166ac", fontface="bold") +
        theme_bw(base_size=14) +
        theme(plot.title=element_text(face="bold", hjust=0.5),
              plot.subtitle=element_text(hjust=0.5, size=9, color="grey40"))

    ggsave(file.path(contrast_dir, "differential_genes_volcano.pdf"),
           p_vol, width=8, height=6, device=cairo_pdf)

    pdf(file.path(contrast_dir, "differential_genes_MA.pdf"), width=8, height=6)
    plotMA(res, ylim=c(-5, 5), main=paste("MA Plot:", contrast_name))
    dev.off()

    cat("    Saved: volcano + MA plots\n")

    all_contrasts_summary <- rbind(all_contrasts_summary, data.frame(
        Contrast      = contrast_name,
        Genes_tested  = nrow(res_df),
        Total_DEGs    = nrow(sig_df),
        Upregulated   = nrow(up_genes),
        Downregulated = nrow(down_genes),
        stringsAsFactors = FALSE
    ))
}

write.csv(all_contrasts_summary,
          file.path(opts$output_dir, "all_contrasts_summary.csv"),
          row.names=FALSE)

cat("\n============================================\n")
cat("  Multi-Contrast Analysis Complete!\n")
cat("============================================\n")
cat("Reference level :", opts$ref_level, "\n")
cat("Contrasts done  :", length(non_ref_levels), "\n\n")
print(all_contrasts_summary, row.names=FALSE)
cat("\nOutput directory:", opts$output_dir, "\n")
cat("============================================\n")
