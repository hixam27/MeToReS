#!/usr/bin/env Rscript
# =============================================================================
# Resistome Summary Plots — Top ARG/MRG Gene Families
# =============================================================================
# Produces two stacked barplots: relative abundance of the top N ARG gene
# families, and the top N MRG compounds, across all samples.
#
# INPUT: abundance_with_resistance.csv 
# =============================================================================

suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(RColorBrewer)
})

# =============================================================================
# Argument parsing
# =============================================================================
option_list <- list(
    make_option(c("--abundance_with_resistance"), type="character",
        help="Path to abundance_with_resistance.csv [required]"),
    make_option(c("--samplesheet"), type="character",
        help="Samplesheet CSV with sample + country columns [required]"),
    make_option(c("--top_n"),       type="integer", default=10,
        help="Number of top families to show individually [default: %default]"),
    make_option(c("--out_dir"),     type="character",
        help="Output directory [required]")
)

opts <- parse_args(OptionParser(option_list=option_list))

for (arg in c("abundance_with_resistance", "samplesheet", "out_dir")) {
    if (is.null(opts[[arg]])) stop(paste("Missing required argument:", arg))
}
dir.create(opts$out_dir, recursive=TRUE, showWarnings=FALSE)

cat("============================================\n")
cat("  Resistome Summary Plots\n")
cat("============================================\n")
cat("Top N families per plot:", opts$top_n, "\n\n")

# =============================================================================
# STEP 1: Load inputs
# =============================================================================
cat("=== STEP 1: Loading inputs ===\n")

df          <- read.csv(opts$abundance_with_resistance, stringsAsFactors=FALSE)
sample_info <- read.csv(opts$samplesheet, stringsAsFactors=FALSE)

required_cols <- c("ORF_id", "is_ARG", "AMR_Gene_Family", "is_MRG", "Compound")
missing_cols  <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
    stop(paste0(
        "abundance_with_resistance.csv is missing expected column(s): ",
        paste(missing_cols, collapse=", "), "\n",
        "Available columns: ", paste(colnames(df), collapse=", "), "\n",
        "(This script expects the column layout confirmed against the ",
        "real pipeline output — if merge_resistance.py's output columns ",
        "have changed, update the column names in this script.)"
    ))
}

if (!"sample" %in% colnames(sample_info) || !"country" %in% colnames(sample_info)) {
    stop(paste0("Samplesheet must have 'sample' and 'country' columns. Available: ",
                paste(colnames(sample_info), collapse=", ")))
}

# Per-sample TPM columns are whatever columns in df match a real sample
# name from the samplesheet (confirmed: sample columns use the short
# 'sample' name, e.g. PL_S_1, not the full SRA identifier).
sample_cols <- intersect(colnames(df), sample_info$sample)
if (length(sample_cols) == 0) {
    stop("No columns in abundance_with_resistance.csv match any sample name in the samplesheet.")
}

cat("  Rows (genes)   :", nrow(df), "\n")
cat("  Sample columns :", length(sample_cols), "\n\n")

# is_ARG / is_MRG are written as "True"/"False" strings (from the Python
# merge step) — as.logical() correctly parses these regardless of case.
df$is_ARG <- as.logical(df$is_ARG)
df$is_MRG <- as.logical(df$is_MRG)

# =============================================================================
# STEP 2: Reusable per-category (ARG or MRG) plotting function
# =============================================================================
make_resistome_plot <- function(flag_col, family_col, category_label, file_tag) {

    cat("\n  ---", category_label, "---\n")

    cat_df <- df[!is.na(df[[flag_col]]) & df[[flag_col]] &
                 !is.na(df[[family_col]]) & df[[family_col]] != "", ,
                 drop=FALSE]
    cat("  ", category_label, " genes (", flag_col, "==TRUE, ",
        family_col, " non-empty):", nrow(cat_df), "\n", sep="")

    if (nrow(cat_df) == 0) {
        cat("  SKIPPED (", category_label, "): zero genes match.\n", sep="")
        return(invisible(NULL))
    }

    # Long format: one row per gene x sample
    long_df <- cat_df %>%
        select(all_of(family_col), all_of(sample_cols)) %>%
        rename(Family = all_of(family_col)) %>%
        tidyr::pivot_longer(cols=all_of(sample_cols), names_to="sample", values_to="tpm")

    # Sum TPM per family per sample
    family_sample_sums <- long_df %>%
        group_by(sample, Family) %>%
        summarise(tpm_sum = sum(tpm, na.rm=TRUE), .groups="drop")

    # Total resistome (this category) abundance per sample — computed
    # across ALL families, not just the top N, so proportions reflect the
    # true composition (bars sum to 100% including "Other").
    sample_totals <- family_sample_sums %>%
        group_by(sample) %>%
        summarise(sample_total = sum(tpm_sum), .groups="drop")

    # Top N families by overall abundance across all samples
    top_families <- family_sample_sums %>%
        group_by(Family) %>%
        summarise(overall_sum = sum(tpm_sum), .groups="drop") %>%
        arrange(desc(overall_sum)) %>%
        slice_head(n=opts$top_n) %>%
        pull(Family)

    cat("  Top", length(top_families), "families:", paste(top_families, collapse=", "), "\n")

    plot_df <- family_sample_sums %>%
        mutate(Family = ifelse(Family %in% top_families, Family, "Other")) %>%
        group_by(sample, Family) %>%
        summarise(tpm_sum = sum(tpm_sum), .groups="drop") %>%
        left_join(sample_totals, by="sample") %>%
        mutate(rel_abundance = ifelse(sample_total > 0, tpm_sum / sample_total, 0)) %>%
        left_join(sample_info[, c("sample", "country")], by="sample")

    family_order <- c(top_families, "Other")
    plot_df$Family <- factor(plot_df$Family, levels=family_order)

    out_csv <- file.path(opts$out_dir, paste0("resistome_", file_tag, "_relative_abundance.csv"))
    write.csv(plot_df %>% arrange(country, sample, Family), out_csv, row.names=FALSE)
    cat("  Saved:", basename(out_csv), "\n")

    n_top <- length(top_families)
    palette_colors <- if (n_top <= 12) {
        RColorBrewer::brewer.pal(max(n_top, 3), "Paired")[seq_len(n_top)]
    } else {
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_top)
    }
    fill_colors <- c(palette_colors, "grey70")
    names(fill_colors) <- family_order

    p <- ggplot(plot_df, aes(x=sample, y=rel_abundance, fill=Family)) +
        geom_col(width=0.85) +
        facet_wrap(~country, scales="free_x", nrow=1) +
        scale_fill_manual(values=fill_colors) +
        scale_y_continuous(labels=scales::percent_format()) +
        labs(
            title = paste0("Top ", n_top, " ", category_label, " Families — Relative Abundance"),
            x     = "Sample",
            y     = "Relative abundance",
            fill  = paste(category_label, "family")
        ) +
        theme_bw(base_size=12) +
        theme(
            plot.title       = element_text(face="bold", hjust=0.5),
            axis.text.x      = element_text(angle=45, hjust=1, size=8),
            strip.background = element_rect(fill="grey90"),
            strip.text       = element_text(face="bold")
        )

    out_pdf <- file.path(opts$out_dir, paste0("resistome_", file_tag, "_barplot.pdf"))
    ggsave(out_pdf, p, width=12, height=6, device=cairo_pdf)
    cat("  Saved:", basename(out_pdf), "\n")
}

# =============================================================================
# STEP 3: Generate both plots
# =============================================================================
cat("=== STEP 2: Generating plots ===\n")

make_resistome_plot("is_ARG", "AMR_Gene_Family", "ARG", "ARG")
make_resistome_plot("is_MRG", "Compound",        "MRG", "MRG")

cat("\n============================================\n")
cat("  Resistome summary plots complete\n")
cat("============================================\n")
cat("Output directory:", opts$out_dir, "\n")
cat("============================================\n")
