#!/usr/bin/env Rscript
###############################################################
## filter_genes.R — Count-aware gene filtering for DESeq2
##
## Strategy (dual-mode):
##   PRIMARY   — Count-based: keeps genes where >cutoff fraction
##               of samples in at least one group have
##               estimated counts >= min_count (default: 5)
##   SECONDARY — TPM-based:   additionally requires TPM >= min_tpm
##               in the same fraction of samples (default: 0.5)
##
## Both conditions must be satisfied so the filtered matrix is
## appropriate for both DESeq2 (count-based) and the TPM-level
## downstream analysis that consumes the filtered CSV.
##
## Output: filtered CSV with GeneID column (TPM values retained
##         for downstream compatibility; DESeq2 rule re-imports
##         raw counts via tximport using these gene IDs as a mask)
###############################################################

suppressWarnings(suppressMessages(library(optparse)))

option_list <- list(
    make_option(c("-i", "--input"), type="character",
        help="Abundance matrix CSV file (TPM, GeneID as first column) [required]"),
    make_option(c("-q", "--quant_dir"), type="character",
        help="Directory containing per-sample Salmon quant folders (for count import) [required]"),
    make_option(c("-s", "--samplesheet"), type="character",
        help="Samplesheet CSV with sample, sra, group columns [required]"),
    make_option(c("-o", "--output"), type="character",
        help="Output filtered matrix CSV file [required]"),
    make_option(c("-g", "--group_col"), type="character", default="group",
        help="Name of grouping column in samplesheet [default: group]"),
    make_option(c("-c", "--cutoff"), type="double", default=0.5,
        help="Min fraction of samples per group passing both thresholds [default: %default]"),
    make_option(c("-t", "--min_tpm"), type="double", default=1.0,
        help="Minimum TPM threshold [default: %default]"),
    make_option(c("-n", "--min_count"), type="double", default=5.0,
        help="Minimum estimated count threshold [default: %default]")
)

opts <- parse_args(OptionParser(option_list=option_list),
    args=commandArgs(trailingOnly=TRUE))

for (arg in c("input", "quant_dir", "samplesheet", "output")) {
    if (is.null(opts[[arg]])) {
        stop(paste("Missing required argument --", arg))
    }
}

# ============================================================
# Load TPM matrix
# ============================================================
cat("========================================\n")
cat("Gene Filtering — Count + TPM dual mode\n")
cat("========================================\n\n")

cat("Loading TPM matrix:", opts$input, "\n")
tpm_matrix <- read.csv(opts$input, header=TRUE, row.names=1,
                        stringsAsFactors=FALSE)
cat("  Genes:", nrow(tpm_matrix), "\n")
cat("  Samples:", ncol(tpm_matrix), "\n\n")

# ============================================================
# Load samplesheet
# ============================================================
cat("Loading samplesheet:", opts$samplesheet, "\n")
samplesheet <- read.csv(opts$samplesheet, header=TRUE,
                         stringsAsFactors=FALSE)
# use the column specified by --group_col from Snakefile
if (opts$group_col != "group") {
    if (!opts$group_col %in% colnames(samplesheet)) {
        stop(paste0("Column '", opts$group_col, "' not found in samplesheet. ",
                    "Available: ", paste(colnames(samplesheet), collapse=", ")))
    }
    colnames(samplesheet)[colnames(samplesheet) == opts$group_col] <- "group"
    cat("  Using column '", opts$group_col, "' as grouping variable\n")
}

groups <- unique(samplesheet$group)
cat("  Groups:", paste(groups, collapse=", "), "\n\n")

# ============================================================
# Load estimated counts from Salmon quant.sf files
# ============================================================
cat("Loading estimated counts from Salmon quant files...\n")

count_list <- list()
for (i in seq_len(nrow(samplesheet))) {
    sra    <- samplesheet$sra[i]
    sample <- samplesheet$sample[i]
    qfile  <- file.path(opts$quant_dir,
                         paste0(sra, "_quant"),
                         "quant.sf")

    if (!file.exists(qfile)) {
        stop(paste("quant.sf not found:", qfile))
    }

    qdata <- read.table(qfile, header=TRUE, sep="\t",
                         stringsAsFactors=FALSE)

    # NumReads column = estimated counts from Salmon
    count_vec <- setNames(qdata$NumReads, qdata$Name)
    count_list[[sample]] <- count_vec
    cat("  Loaded", sample, "(", sra, "):", length(count_vec), "genes\n")
}

# Build count matrix aligned to TPM matrix rows
all_genes    <- rownames(tpm_matrix)
count_matrix <- do.call(cbind, lapply(count_list, function(cv) {
    # Genes in TPM matrix but missing from quant → 0
    counts <- cv[all_genes]
    counts[is.na(counts)] <- 0
    counts
}))
rownames(count_matrix) <- all_genes

cat("\n  Count matrix dimensions:", nrow(count_matrix), "x", ncol(count_matrix), "\n\n")

# ============================================================
# Dual-mode filtering per group
# ============================================================
cat("Filtering parameters:\n")
cat("  Min fraction per group  :", opts$cutoff, "\n")
cat("  Min TPM                 :", opts$min_tpm, "\n")
cat("  Min estimated count     :", opts$min_count, "\n\n")

passing_genes <- c()

for (grp in groups) {
    grp_samples <- samplesheet$sample[samplesheet$group == grp]

    # Subset to samples that exist in both matrices
    grp_samples_tpm   <- intersect(grp_samples, colnames(tpm_matrix))
    grp_samples_count <- intersect(grp_samples, colnames(count_matrix))

    if (length(grp_samples_tpm) == 0) {
        cat("  WARNING: No samples found for group", grp, "in TPM matrix — skipping\n")
        next
    }

    grp_tpm   <- tpm_matrix[, grp_samples_tpm,   drop=FALSE]
    grp_count <- count_matrix[, grp_samples_count, drop=FALSE]

    n_samples <- length(grp_samples_tpm)

    # TPM filter: fraction of samples >= min_tpm
    tpm_pass_frac <- apply(grp_tpm, 1, function(x) {
        sum(x >= opts$min_tpm, na.rm=TRUE) / n_samples
    })

    # Count filter: fraction of samples >= min_count
    count_pass_frac <- apply(grp_count, 1, function(x) {
        sum(x >= opts$min_count, na.rm=TRUE) / n_samples
    })

    # Gene passes if BOTH filters exceed cutoff
    grp_passing <- names(which(tpm_pass_frac   > opts$cutoff &
                                count_pass_frac > opts$cutoff))

    cat("  Group", grp, "(n=", n_samples, "):\n")
    cat("    TPM filter passing    :", sum(tpm_pass_frac   > opts$cutoff), "\n")
    cat("    Count filter passing  :", sum(count_pass_frac > opts$cutoff), "\n")
    cat("    Both filters passing  :", length(grp_passing), "\n\n")

    passing_genes <- union(passing_genes, grp_passing)
}

# ============================================================
# Apply filter and save
# ============================================================
filtered_tpm <- tpm_matrix[passing_genes, , drop=FALSE]

cat("========================================\n")
cat("Filtering summary:\n")
cat("  Genes before filtering:", nrow(tpm_matrix), "\n")
cat("  Genes after filtering :", nrow(filtered_tpm), "\n")
cat("  Genes removed         :", nrow(tpm_matrix) - nrow(filtered_tpm), "\n")
cat("  Reduction             :", round((1 - nrow(filtered_tpm)/nrow(tpm_matrix))*100, 1), "%\n")
cat("========================================\n\n")

# Write output — GeneID as explicit column (matches downstream expectations)
out_df <- data.frame(GeneID=rownames(filtered_tpm),
                      filtered_tpm,
                      check.names=FALSE)
write.csv(out_df, opts$output, row.names=FALSE)
cat("Saved filtered matrix:", opts$output, "\n")
cat("  Dimensions:", nrow(out_df), "genes x", ncol(out_df)-1, "samples\n")
