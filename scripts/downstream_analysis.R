#!/usr/bin/env Rscript
###############################################################
## downstream_analysis.R — Multi-Contrast Version
##
## SHARED outputs (in workdir/):
##   - gene_pcoa.pdf, anosim_result.csv
##   - upset_all_DEGs.pdf, upset_upregulated.pdf
##   - all_contrast_DEG_counts.pdf, all_contrast_DEG_counts.csv
##   - kaiju_community_barplot.pdf, kaiju_community_table.csv
##
## PER-CONTRAST outputs (in workdir/per_contrast/{contrast}/):
##   - DEG_up_heatmap.pdf
##   - cog.pdf
##   - go_rich_bar.pdf, go_rich_dot.pdf
##   - ko_rich_bar.pdf, ko_rich_dot.pdf
##   - go_classification.pdf
##
###############################################################

suppressWarnings(suppressMessages(library(optparse)))

option_list <- list(
    make_option(c("-w", "--workdir"),       type="character",
        help="Output directory for downstream results [required]"),
    make_option(c("-d", "--deg_dir"),       type="character", default=NULL,
        help="Parent DESeq2 results directory containing contrast subfolders"),
    make_option(c("-k", "--kaiju_dir"),     type="character", default=NULL),
    make_option(c("-e", "--emapper_dir"),   type="character", default=NULL),
    make_option(c("-c", "--contrasts"),     type="character", default=NULL,
        help="Space-separated list of contrast names (subfolder names)"),
    make_option(c("-m", "--matrix"),        type="character",
        default="transcript_abundance_quantification_table_filter.csv"),
    make_option(c("-T", "--tax_matrix"),    type="character",
        default="transcript_abundance_quantification_table_filter_taxonomy.csv"),
    make_option(c("-a", "--emapper_all"),   type="character",
        default="all_genes_emapper.emapper.annotations"),
    make_option(c("-s", "--sample_group"),  type="character",
        default="sample_group.csv"),
    make_option(c("--samplesheet"),         type="character", default=NULL),
    make_option(c("-f", "--taxon_filter"),  type="character", default="Streptophyta"),
    make_option(c("--heatmap_min_sum"),     type="numeric",   default=50),
    make_option(c("--top_n_kaiju"),         type="integer",   default=15),
    make_option(c("--cog_funclass"),        type="character", default=NULL)
)

opts <- parse_args(OptionParser(option_list=option_list),
    args=commandArgs(trailingOnly=TRUE))

if (is.null(opts$workdir))   stop("--workdir is required")
if (is.null(opts$contrasts)) stop("--contrasts is required (space-separated list)")

# Parse contrasts
CONTRASTS <- strsplit(trimws(opts$contrasts), "\\s+")[[1]]
CONTRASTS <- CONTRASTS[nchar(CONTRASTS) > 0]
if (length(CONTRASTS) == 0) stop("No contrasts parsed from --contrasts")

cat("========================================\n")
cat("MeToReS Multi-Contrast Downstream Analysis\n")
cat("Contrasts:", length(CONTRASTS), "\n")
for (c in CONTRASTS) cat("  -", c, "\n")
cat("========================================\n\n")

# ── Resolve paths to absolute ────────────────────────────────
opts$workdir <- normalizePath(opts$workdir, mustWork=FALSE)
if (!is.null(opts$deg_dir))     opts$deg_dir     <- normalizePath(opts$deg_dir, mustWork=FALSE)
if (!is.null(opts$kaiju_dir))   opts$kaiju_dir   <- normalizePath(opts$kaiju_dir, mustWork=FALSE)
if (!is.null(opts$emapper_dir)) opts$emapper_dir <- normalizePath(opts$emapper_dir, mustWork=FALSE)
if (!is.null(opts$samplesheet)) opts$samplesheet <- normalizePath(opts$samplesheet, mustWork=FALSE)

dir.create(opts$workdir, showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(opts$workdir, "per_contrast"),
           showWarnings=FALSE, recursive=TRUE)
setwd(opts$workdir)

results_root <- dirname(opts$workdir)
if (is.null(opts$deg_dir))     opts$deg_dir     <- file.path(results_root, "deg")
if (is.null(opts$kaiju_dir))   opts$kaiju_dir   <- file.path(results_root, "qc", "kaiju")
if (is.null(opts$emapper_dir)) opts$emapper_dir <- file.path(results_root, "emapper")

# Resolve paths
matrix_path       <- file.path(opts$deg_dir, opts$matrix)
tax_matrix_path   <- file.path(opts$deg_dir, opts$tax_matrix)
emapper_all_path  <- file.path(opts$emapper_dir, opts$emapper_all)
sample_group_path <- opts$sample_group

cat("Working directory:", getwd(), "\n")
cat("DEG parent dir   :", opts$deg_dir, "\n")
cat("eggNOG directory :", opts$emapper_dir, "\n")
cat("Kaiju directory  :", opts$kaiju_dir, "\n\n")

suppressWarnings(suppressMessages({
    library(tidyverse)
    library(vegan)
    library(ape)
    library(ggplot2)
    library(dplyr)
    library(pheatmap)
    library(RColorBrewer)
    library(stringr)
    library(ComplexUpset)
}))

# ── eggNOG-mapper file reader ────────────────────────────────
read_emapper <- function(path) {
    if (!file.exists(path) || file.info(path)$size == 0) {
        cat("  WARNING: emapper file missing or empty:", path, "\n")
        return(data.frame(query=character(0), stringsAsFactors=FALSE))
    }
    lines <- readLines(path)
    lines <- lines[!grepl("^##", lines)]
    if (length(lines) == 0) return(data.frame(query=character(0), stringsAsFactors=FALSE))
    lines[1] <- sub("^#", "", lines[1])
    df <- read.table(text=paste(lines, collapse="\n"),
                     header=TRUE, sep="\t", quote="",
                     stringsAsFactors=FALSE, fill=TRUE, comment.char="")
    qcol <- grep("^query$|^X.query$|^.query$", colnames(df), value=TRUE)[1]
    if (!is.na(qcol) && qcol != "query")
        colnames(df)[colnames(df) == qcol] <- "query"
    df
}

###############################################################
## SECTION 1: Host contamination removal (Kaiju-derived IDs)
###############################################################
cat("=== SECTION 1: Removing", opts$taxon_filter, "===\n")

if (!file.exists(matrix_path))
    stop(paste("Filtered abundance matrix not found:", matrix_path))

table_filter <- read.csv(matrix_path, header=TRUE, row.names=1)
cat("  Genes before removal:", nrow(table_filter), "\n")

kaiju_filter_candidates <- c(
    file.path(opts$kaiju_dir, paste0(tolower(opts$taxon_filter), "_ids.csv")),
    file.path(opts$kaiju_dir, paste0(opts$taxon_filter, "_ids.csv")),
    file.path(opts$kaiju_dir, paste0(tolower(opts$taxon_filter), "_id.csv"))
)
kaiju_filter_path <- kaiju_filter_candidates[file.exists(kaiju_filter_candidates)][1]

if (!is.na(kaiju_filter_path)) {
    if (file.info(kaiju_filter_path)$size == 0) {
        cat("  Kaiju filter list empty — matrix unchanged\n")
        table_filter2 <- table_filter
    } else {
        contam_raw <- read.csv(kaiju_filter_path, header=FALSE, stringsAsFactors=FALSE)
        contam_ids <- contam_raw[[1]]
        contam_ids <- contam_ids[nchar(trimws(contam_ids)) > 0]
        overlap <- intersect(rownames(table_filter), contam_ids)
        cat("  Contamination genes in matrix :", length(overlap), "\n")
        table_filter2 <- if (length(overlap) > 0)
            table_filter[!rownames(table_filter) %in% contam_ids, , drop=FALSE]
        else table_filter
    }
} else {
    cat("  No contamination source found. Skipping removal.\n")
    table_filter2 <- table_filter
}

cat("  Genes after :", nrow(table_filter2), "\n")
write.csv(
    data.frame(GeneID=rownames(table_filter2), table_filter2),
    "transcript_abundance_quantification_table_filter2.csv", row.names=FALSE)
cat("\n")

###############################################################
## SECTION 2: PCoA ordination + ANOSIM
###############################################################
cat("=== SECTION 2: PCoA + ANOSIM ===\n")

if (!file.exists(sample_group_path))
    stop(paste("Sample group file not found:", sample_group_path))

group <- read.csv(sample_group_path, header=TRUE, stringsAsFactors=FALSE)
names(group) <- c("sample", "group")
cat("  Groups:", paste(unique(group$group), collapse=", "), "\n")

gene_dist <- vegdist(t(table_filter2), "bray")
res       <- pcoa(gene_dist)
vectors   <- res$vectors

x_value <- res$values[1, 2]
y_value <- res$values[2, 2]

vectors_df    <- data.frame(sampleID=rownames(vectors), vectors)
group_df      <- data.frame(sampleID=group$sample, group=group$group)
vectors_group <- left_join(vectors_df, group_df, by="sampleID")

p <- ggplot(vectors_group, aes(Axis.1, Axis.2)) +
    geom_point(aes(colour=factor(group)), size=5) +
    xlab(paste0("PCoA1 (", round(x_value * 100, 2), "%)")) +
    ylab(paste0("PCoA2 (", round(y_value * 100, 2), "%)")) +
    geom_vline(xintercept=0, linetype="dotted", linewidth=1.2) +
    geom_hline(yintercept=0, linetype="dotted", linewidth=1.2) +
    theme_bw() +
    theme(text=element_text(face="bold", size=16),
          panel.background=element_rect(fill="white", color="black", linewidth=2),
          panel.grid=element_blank(),
          axis.text=element_text(face="bold", size=14))

ggsave(p, file="gene_pcoa.pdf", width=9, height=7, device=cairo_pdf)
cat("  Saved: gene_pcoa.pdf\n")

anosim_result <- anosim(t(table_filter2), vectors_group$group,
                         permutations=999, distance="bray")
cat("  ANOSIM R:", round(anosim_result$statistic, 4),
    "  p-value:", anosim_result$signif, "\n")
write.csv(
    data.frame(Statistic=round(anosim_result$statistic, 6),
               P_value=anosim_result$signif, Permutations=999,
               Method="Bray-Curtis"),
    "anosim_result.csv", row.names=FALSE)
cat("\n")

###############################################################
## SECTION 3: Load all per-contrast DEG lists
###############################################################
cat("=== SECTION 3: Loading per-contrast DEGs ===\n")

# Build named lists: contrast → vector of gene IDs
all_DEGs_list      <- list()
upregulated_list   <- list()
downregulated_list <- list()

deg_counts <- data.frame(
    contrast      = character(),
    total_DEGs    = integer(),
    upregulated   = integer(),
    downregulated = integer(),
    stringsAsFactors = FALSE
)

for (contrast in CONTRASTS) {
    contrast_dir <- file.path(opts$deg_dir, contrast)
    all_file  <- file.path(contrast_dir, "differential_genes_id.txt")
    up_file   <- file.path(contrast_dir, "upregulated_genes_id.txt")
    down_file <- file.path(contrast_dir, "downregulated_genes_id.txt")

    if (!file.exists(all_file)) {
        cat("  WARNING: skipping", contrast, "— differential_genes_id.txt not found\n")
        next
    }

    read_ids <- function(f) {
        if (!file.exists(f) || file.info(f)$size == 0) return(character(0))
        ids <- readLines(f)
        ids[nchar(trimws(ids)) > 0]
    }

    all_ids  <- read_ids(all_file)
    up_ids   <- read_ids(up_file)
    down_ids <- read_ids(down_file)

    all_DEGs_list[[contrast]]      <- all_ids
    upregulated_list[[contrast]]   <- up_ids
    downregulated_list[[contrast]] <- down_ids

    deg_counts <- rbind(deg_counts, data.frame(
        contrast      = contrast,
        total_DEGs    = length(all_ids),
        upregulated   = length(up_ids),
        downregulated = length(down_ids),
        stringsAsFactors = FALSE
    ))

    cat("  ", contrast, ": total =", length(all_ids),
        "| up =", length(up_ids),
        "| down =", length(down_ids), "\n")
}

write.csv(deg_counts, "all_contrast_DEG_counts.csv", row.names=FALSE)

# Bar chart of DEG counts per contrast
deg_counts_long <- deg_counts %>%
    tidyr::pivot_longer(cols=c(upregulated, downregulated),
                        names_to="direction", values_to="count")

p_counts <- ggplot(deg_counts_long,
                   aes(x=reorder(contrast, -count, sum), y=count, fill=direction)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values=c(upregulated="#b2182b", downregulated="#2166ac")) +
    labs(title="DEG counts per contrast", x="", y="Number of DEGs") +
    theme_classic(base_size=12) +
    theme(plot.title=element_text(hjust=0.5, face="bold"),
          axis.text.x=element_text(angle=45, hjust=1))

ggsave("all_contrast_DEG_counts.pdf", p_counts, width=10, height=6, device=cairo_pdf)
cat("  Saved: all_contrast_DEG_counts.pdf\n\n")

###############################################################
## SECTION 4: UpSet plots across contrasts
###############################################################
cat("=== SECTION 4: UpSet plots ===\n")

# ── Helper: build presence/absence dataframe for UpSet ───────
build_upset_df <- function(gene_lists) {
    all_genes <- unique(unlist(gene_lists))
    if (length(all_genes) == 0) return(NULL)
    df <- data.frame(gene_id = all_genes, stringsAsFactors=FALSE)
    for (name in names(gene_lists)) {
        df[[name]] <- df$gene_id %in% gene_lists[[name]]
    }
    df
}

# ── UpSet plot helper ────────────────────────────────────────
make_upset <- function(df, contrast_names, title, outfile) {
    if (is.null(df) || nrow(df) == 0) {
        cat("  Skipping (no genes):", outfile, "\n")
        pdf(outfile); plot.new(); title(paste(title, "- no DEGs")); dev.off()
        return(invisible(NULL))
    }
    # ComplexUpset requires at least 1 set with ≥1 element
    set_sizes <- sapply(contrast_names, function(n) sum(df[[n]]))
    keep_sets <- names(set_sizes)[set_sizes > 0]
    if (length(keep_sets) < 2) {
        cat("  Skipping (fewer than 2 non-empty sets):", outfile, "\n")
        pdf(outfile); plot.new(); title(paste(title, "- insufficient sets")); dev.off()
        return(invisible(NULL))
    }
    p <- ComplexUpset::upset(
        df,
        intersect = keep_sets,
        name      = "",
        width_ratio = 0.18,
        min_size    = 1,
        sort_intersections_by = "cardinality",
        base_annotations = list(
            'Intersection size' = ComplexUpset::intersection_size(
                text = list(size=3),
                mapping = aes(fill='bars_color')
            ) + scale_fill_manual(values=c('bars_color'='#5B8C5A'), guide='none')
        ),
        themes = ComplexUpset::upset_modify_themes(list(
            'intersections_matrix' = theme(text=element_text(size=10))
        ))
    ) + ggtitle(title) +
        theme(plot.title=element_text(hjust=0.5, face="bold", size=13))

    ggsave(outfile, p, width=max(8, length(keep_sets)*1.4), height=6, device=cairo_pdf)
    cat("  Saved:", outfile, "\n")
}

# 1. All DEGs across contrasts
upset_all <- build_upset_df(all_DEGs_list)
make_upset(upset_all, CONTRASTS,
           "DEG overlap across contrasts (all DEGs)",
           "upset_all_DEGs.pdf")

# 2. Upregulated only
upset_up <- build_upset_df(upregulated_list)
make_upset(upset_up, CONTRASTS,
           "Upregulated gene overlap across contrasts",
           "upset_upregulated.pdf")

# Save the presence/absence matrices for reference
if (!is.null(upset_all)) {
    write.csv(upset_all, "upset_all_DEGs_matrix.csv", row.names=FALSE)
}
if (!is.null(upset_up)) {
    write.csv(upset_up, "upset_upregulated_matrix.csv", row.names=FALSE)
}
cat("\n")

###############################################################
## SECTION 5: Per-contrast DEG ID lists
###############################################################
cat("=== SECTION 5: Per-contrast gene set files ===\n")

for (contrast in names(all_DEGs_list)) {
    contrast_out <- file.path("per_contrast", contrast)
    dir.create(contrast_out, showWarnings=FALSE, recursive=TRUE)

    # Copy gene lists into per_contrast folder for reference
    writeLines(all_DEGs_list[[contrast]],
               file.path(contrast_out, "DEGs_all.txt"))
    writeLines(upregulated_list[[contrast]],
               file.path(contrast_out, "DEGs_upregulated.txt"))
    writeLines(downregulated_list[[contrast]],
               file.path(contrast_out, "DEGs_downregulated.txt"))
}
cat("  Per-contrast DEG lists saved\n\n")

###############################################################
## SECTION 6: Per-sample Kaiju community composition
###############################################################
cat("=== SECTION 6: Per-sample Kaiju Community Composition ===\n")

kaiju_sample_files <- list.files(opts$kaiju_dir,
    pattern="_kaiju_summary\\.tsv$", full.names=TRUE)

if (length(kaiju_sample_files) == 0) {
    cat("  WARNING: No per-sample Kaiju files found — skipping.\n\n")
} else {
    cat("  Found", length(kaiju_sample_files), "files\n")

    if (!is.null(opts$samplesheet) && file.exists(opts$samplesheet)) {
        ss_map <- read.csv(opts$samplesheet, stringsAsFactors=FALSE)
        sra_to_sample <- setNames(ss_map$sample, ss_map$sra)
    } else {
        sra_to_sample <- c()
    }

    kaiju_all <- lapply(kaiju_sample_files, function(f) {
        sra_id <- sub("_kaiju_summary\\.tsv$", "", basename(f))
        sample_name <- ifelse(sra_id %in% names(sra_to_sample),
                              sra_to_sample[sra_id], sra_id)
        df <- read.table(f, header=TRUE, sep="\t",
                         stringsAsFactors=FALSE, quote="")
        colnames(df) <- tolower(colnames(df))
        pct_col  <- grep("percent|fraction", colnames(df), value=TRUE)[1]
        name_col <- grep("taxon_name|name", colnames(df), value=TRUE)[1]
        if (is.na(pct_col) || is.na(name_col)) return(NULL)
        df <- df[, c(pct_col, name_col)]
        colnames(df) <- c("percent", "taxon")
        df$sample <- sample_name
        df
    })
    kaiju_all <- bind_rows(Filter(Negate(is.null), kaiju_all))

    kaiju_all <- kaiju_all %>%
        filter(!grepl("unclassified|cannot|root", taxon, ignore.case=TRUE),
               percent > 0)

    top_taxa <- kaiju_all %>% group_by(taxon) %>%
        summarise(mean_pct=mean(percent), .groups="drop") %>%
        arrange(desc(mean_pct)) %>%
        slice_head(n=opts$top_n_kaiju) %>% pull(taxon)

    kaiju_plot <- kaiju_all %>%
        filter(taxon %in% top_taxa) %>%
        mutate(taxon=factor(taxon, levels=rev(top_taxa)))

    group_map <- setNames(group$group, group$sample)
    kaiju_plot$group <- group_map[kaiju_plot$sample]

    p_kaiju <- ggplot(kaiju_plot, aes(x=sample, y=percent, fill=taxon)) +
        geom_bar(stat="identity", position="stack") +
        facet_wrap(~group, scales="free_x") +
        labs(title=paste("Per-sample Kaiju: Top", opts$top_n_kaiju, "Phyla"),
             x="Sample", y="% of classified reads", fill="Phylum") +
        theme_classic() +
        theme(plot.title=element_text(hjust=0.5, size=14, face="bold"),
              axis.text.x=element_text(angle=45, hjust=1, size=9),
              legend.text=element_text(size=8),
              legend.key.size=unit(0.8, "line"),
              strip.background=element_rect(fill="grey90", color="black")) +
        guides(fill=guide_legend(ncol=1))

    ggsave("kaiju_community_barplot.pdf", p_kaiju, width=14, height=7)
    cat("  Saved: kaiju_community_barplot.pdf\n")

    kaiju_wide <- kaiju_plot %>%
        dplyr::select(sample, taxon, percent) %>%
        pivot_wider(names_from=sample, values_from=percent, values_fill=0)
    write.csv(kaiju_wide, "kaiju_community_table.csv", row.names=FALSE)
    cat("  Saved: kaiju_community_table.csv\n\n")
}

###############################################################
## SECTION 7: Build gene-GO/KO mapping tables (universe)
###############################################################
cat("=== SECTION 7: Gene-GO/KO mappings (universe) ===\n")

if (!file.exists(emapper_all_path))
    stop(paste("eggNOG all-gene annotations not found:", emapper_all_path))

emapper_all_df <- read_emapper(emapper_all_path)
cat("  Universe genes:", nrow(emapper_all_df), "\n")

gene2go_all <- emapper_all_df %>%
    dplyr::select(GID=query, GO=GOs) %>%
    separate_rows(GO, sep=",", convert=FALSE) %>%
    dplyr::filter(!is.na(GO), GO != "", GO != "-") %>%
    mutate(GO=trimws(GO))

gene2ko_all <- emapper_all_df %>%
    dplyr::select(GID=query, KO=KEGG_ko) %>%
    separate_rows(KO, sep=",", convert=FALSE) %>%
    dplyr::filter(!is.na(KO), KO != "", KO != "-") %>%
    mutate(KO=trimws(KO))

gene2pathway_all <- emapper_all_df %>%
    dplyr::select(GID=query, Pathway=KEGG_Pathway) %>%
    separate_rows(Pathway, sep=",", convert=FALSE) %>%
    dplyr::filter(!is.na(Pathway), Pathway != "", Pathway != "-") %>%
    mutate(Pathway=trimws(Pathway))

cat("  Gene-GO mappings     :", nrow(gene2go_all), "\n")
cat("  Gene-KO mappings     :", nrow(gene2ko_all), "\n")
cat("  Gene-Pathway mappings:", nrow(gene2pathway_all), "\n")

write.csv(gene2go_all,      "gene2go_all.csv",      row.names=FALSE)
write.csv(gene2ko_all,      "gene2ko_all.csv",      row.names=FALSE)
write.csv(gene2pathway_all, "gene2pathway_all.csv", row.names=FALSE)

# Build GO term table from GO.db
go_class <- NULL
if (requireNamespace("GO.db", quietly=TRUE)) {
    suppressWarnings(suppressMessages({
        library(GO.db)
        go_df <- as.data.frame(GOTERM)
        go_class <- data.frame(
            GO_ID       = go_df$go_id,
            Description = go_df$Term,
            Ontology    = go_df$Ontology,
            stringsAsFactors = FALSE
        )
        go_class <- go_class[!is.na(go_class$GO_ID) & go_class$GO_ID != "", ]
    }))
    cat("  GO term table from GO.db:", nrow(go_class), "terms\n\n")
} else {
    cat("  WARNING: GO.db not available — GO enrichment plots will be skipped.\n\n")
}

universe_genes <- unique(emapper_all_df$query)

###############################################################
## ORA + plotting helpers
###############################################################
ora_hypergeometric <- function(gene_list, term2gene, term2name=NULL,
                                universe, p_adjust="BH") {
    gene_list <- unique(gene_list)
    universe  <- unique(universe)
    N <- length(universe); k <- length(gene_list)
    if (k == 0 || N == 0) return(NULL)

    terms <- unique(term2gene[[1]])
    results <- lapply(terms, function(term) {
        term_genes <- unique(term2gene[[2]][term2gene[[1]] == term])
        M <- length(term_genes)
        x <- length(intersect(gene_list, term_genes))
        if (x == 0) return(NULL)
        pval <- phyper(x - 1, M, N - M, k, lower.tail=FALSE)
        data.frame(
            ID=term, Count=x,
            GeneRatio=paste0(x, "/", k),
            BgRatio=paste0(M, "/", N),
            pvalue=pval,
            geneID=paste(intersect(gene_list, term_genes), collapse="/"),
            stringsAsFactors=FALSE
        )
    })
    results <- do.call(rbind, Filter(Negate(is.null), results))
    if (is.null(results) || nrow(results) == 0) return(NULL)
    results$p.adjust <- p.adjust(results$pvalue, method=p_adjust)
    results <- results[order(results$p.adjust), ]
    if (!is.null(term2name)) {
        colnames(term2name) <- c("ID", "Description")
        results <- merge(results, term2name, by="ID", all.x=TRUE)
        results$Description[is.na(results$Description)] <- results$ID[is.na(results$Description)]
    } else {
        results$Description <- results$ID
    }
    results
}

plot_ora_bar <- function(res, title, top_n=30) {
    if (is.null(res) || nrow(res) == 0)
        return(ggplot() + labs(title=paste(title, "— no enriched terms")))
    color_col <- if (diff(range(res$p.adjust, na.rm=TRUE)) > 0) "p.adjust" else "pvalue"
    df <- head(res[order(res[[color_col]]), ], top_n)
    dup <- duplicated(df$Description) | duplicated(df$Description, fromLast=TRUE)
    df$Description[dup] <- paste0(df$Description[dup], " (", df$ID[dup], ")")
    df$Description <- factor(df$Description, levels=rev(unique(df$Description)))
    df$fill_val <- df[[color_col]]
    ggplot(df, aes(x=Count, y=Description, fill=fill_val)) +
        geom_bar(stat="identity") +
        scale_fill_gradient(low="red", high="blue", name=color_col) +
        labs(title=title, x="Gene count", y="") +
        theme_classic() +
        theme(plot.title=element_text(hjust=0.5, size=11),
              axis.text.y=element_text(size=8))
}

plot_ora_dot <- function(res, title, top_n=30) {
    if (is.null(res) || nrow(res) == 0)
        return(ggplot() + labs(title=paste(title, "— no enriched terms")))
    color_col <- if (diff(range(res$p.adjust, na.rm=TRUE)) > 0) "p.adjust" else "pvalue"
    df <- head(res[order(res[[color_col]]), ], top_n)
    dup <- duplicated(df$Description) | duplicated(df$Description, fromLast=TRUE)
    df$Description[dup] <- paste0(df$Description[dup], " (", df$ID[dup], ")")
    df$Description <- factor(df$Description, levels=rev(unique(df$Description)))
    df$ratio_num <- sapply(strsplit(df$GeneRatio, "/"),
                            function(x) as.numeric(x[1]) / as.numeric(x[2]))
    df$color_val <- df[[color_col]]
    ggplot(df, aes(x=ratio_num, y=Description, size=Count, color=color_val)) +
        geom_point() +
        scale_color_gradient(low="red", high="blue", name=color_col) +
        labs(title=title, x="Gene ratio", y="", size="Count") +
        theme_classic() +
        theme(plot.title=element_text(hjust=0.5, size=11),
              axis.text.y=element_text(size=8))
}

###############################################################
## SECTION 8: Per-contrast functional analysis loop
###############################################################
cat("=== SECTION 8: Per-contrast functional analysis ===\n")

# COG funclass table (shared across contrasts)
cog_info <- NULL
if (!is.null(opts$cog_funclass) && file.exists(opts$cog_funclass)) {
    cog_info <- read.delim(opts$cog_funclass, sep="\t", header=TRUE,
                           stringsAsFactors=FALSE)
    cat("  Loaded COG funclass table:", nrow(cog_info), "categories\n")
} else {
    cat("  WARNING: cog_funclass.tab not found — COG plots will be empty placeholders\n")
}

for (contrast in CONTRASTS) {
    cat("\n  --- Contrast:", contrast, "---\n")
    contrast_out <- file.path("per_contrast", contrast)
    dir.create(contrast_out, showWarnings=FALSE, recursive=TRUE)

    up_genes   <- upregulated_list[[contrast]]
    down_genes <- downregulated_list[[contrast]]
    all_degs   <- all_DEGs_list[[contrast]]

    cat("    Genes: up =", length(up_genes),
        ", down =", length(down_genes),
        ", total =", length(all_degs), "\n")

# ── Subset the all-genes annotation by this contrast's DEG lists ────
    emapper_up_df   <- emapper_all_df[emapper_all_df$query %in% up_genes, ]
    emapper_down_df <- emapper_all_df[emapper_all_df$query %in% down_genes, ]

    # ── Reusable per-direction functional heatmap ────────────
    make_direction_heatmap <- function(direction_label, genes, emapper_df) {
        out_pdf <- file.path(contrast_out,
                             paste0("DEG_", direction_label, "_heatmap.pdf"))

        if (nrow(emapper_df) == 0 || length(genes) == 0) {
            pdf(out_pdf); plot.new()
            title(paste0("No emapper data (", direction_label, ")")); dev.off()
            cat("    Skipped DEG_", direction_label,
                "_heatmap (no emapper data)\n", sep="")
            return(invisible(NULL))
        }

        table_filter2_df <- data.frame(geneID=rownames(table_filter2), table_filter2)
        tab_join <- left_join(table_filter2_df, emapper_df,
                              by=c("geneID"="query"))
        tab_join <- tab_join[complete.cases(tab_join), ]

        if (nrow(tab_join) == 0) {
            pdf(out_pdf); plot.new()
            title(paste0("No annotated ", direction_label, "regulated genes")); dev.off()
            return(invisible(NULL))
        }

        merge_otu <- function(data, formula, FUN) {
            df <- as.data.frame(data); df[is.na(df)] <- 0
            aggregate(formula, data=df, FUN)
        }
        sample_cols  <- colnames(table_filter2)
        cols_to_keep <- intersect(c(sample_cols, "Description"),
                                   colnames(tab_join))
        if (!"Description" %in% cols_to_keep) {
            pdf(out_pdf); plot.new(); title("No Description column"); dev.off()
            return(invisible(NULL))
        }

        tab_sum <- merge_otu(data=tab_join[, cols_to_keep],
                             formula=. ~ Description, FUN=sum)
        rownames(tab_sum) <- tab_sum[, 1]
        tab_sum <- tab_sum[, -1]
        if ("-" %in% rownames(tab_sum))
            tab_sum <- tab_sum[rownames(tab_sum) != "-", ]
        tab_sum <- tab_sum[which(rowSums(tab_sum) >= opts$heatmap_min_sum), ]

        if (nrow(tab_sum) == 0) {
            pdf(out_pdf); plot.new()
            title("No categories above threshold"); dev.off()
            cat("    DEG_", direction_label,
                "_heatmap: no categories above min_sum threshold\n", sep="")
            return(invisible(NULL))
        }

        # Up = blue→red (enriched), Down = green→purple (depleted) for visual cue
        pal <- if (direction_label == "up")
            colorRampPalette(c("blue", "yellow", "red"))(50)
        else
            colorRampPalette(c("darkgreen", "yellow", "purple"))(50)

        x1 <- scale(as.data.frame(tab_sum), center=FALSE, scale=TRUE)

        # Guard: if scaled data has no variance (all rows nearly identical),
        # pheatmap fails with "breaks are not unique". Save placeholder instead.
        x1_clean <- x1[is.finite(rowSums(x1)), , drop=FALSE]
        if (nrow(x1_clean) < 2 || all(apply(x1_clean, 2, var) == 0)) {
            pdf(out_pdf); plot.new()
            title(paste0("Too few/uniform categories for heatmap (",
                         nrow(x1_clean), " rows)"))
            dev.off()
            cat("    DEG_", direction_label,
                "_heatmap: matrix too small/uniform — placeholder saved\n",
                sep="")
            return(invisible(NULL))
        }

        tryCatch({
            pheatmap(x1_clean, color=pal, border_color="black",
                     cluster_row=FALSE, cluster_col=FALSE,
                     scale="none", cellwidth=18, cellheight=18,
                     fontsize_row=10, fontsize_col=10, angle_col="315",
                     main=paste0(contrast, " — ", direction_label, "regulated"),
                     filename=out_pdf)
            while (!is.null(dev.list())) dev.off()
            cat("    Saved: DEG_", direction_label, "_heatmap.pdf (",
                nrow(x1_clean), " categories)\n", sep="")
        }, error = function(e) {
            while (!is.null(dev.list())) dev.off()
            pdf(out_pdf); plot.new()
            title(paste0("Heatmap render failed: ",
                         substr(conditionMessage(e), 1, 50)))
            dev.off()
            cat("    DEG_", direction_label,
                "_heatmap: pheatmap error (", conditionMessage(e), ") — placeholder saved\n",
                sep="")
        })
    }

    # Generate BOTH heatmaps
    make_direction_heatmap("up",   up_genes,   emapper_up_df)
    make_direction_heatmap("down", down_genes, emapper_down_df)

# ── COG classification (upregulated genes) ───────────────
    if (!is.null(cog_info) && length(up_genes) > 0 &&
        nrow(emapper_up_df) > 0) {

        emapper_contrast <- emapper_up_df
        cogs <- emapper_contrast %>%
            dplyr::select(GID=query, COG=COG_category) %>%
            dplyr::filter(!is.na(COG), COG != "-")

        if (nrow(cogs) > 0) {
            df_temp <- lapply(seq_len(nrow(cogs)), function(i) {
                the_gid  <- cogs[i, "GID"]
                the_cogs <- str_trim(
                    str_split(cogs[i, "COG"], "", simplify=FALSE)[[1]])
                tibble(GID=rep(the_gid, length(the_cogs)), COG=the_cogs)
            })
            gene2cog <- bind_rows(df_temp) %>%
                left_join(cog_info, by="COG") %>%
                filter(COG != "-", !is.na(COG_Name))

            if (nrow(gene2cog) > 0) {
                p <- ggplot(gene2cog) +
                    geom_bar(aes(x=COG, fill=COG_Name)) +
                    labs(title=paste("COG Classification —", contrast),
                         x="", y="Number of genes") +
                    theme_classic() +
                    theme(plot.title=element_text(hjust=0.5, face="bold"),
                          legend.title=element_blank(),
                          legend.key.size=unit(1, "line"),
                          legend.text=element_text(size=7.5)) +
                    guides(fill=guide_legend(ncol=1))
                ggsave(file.path(contrast_out, "cog.pdf"), p, width=16, height=7)
                cat("    Saved: cog.pdf\n")
            } else {
                pdf(file.path(contrast_out, "cog.pdf"))
                plot.new(); title("No COG hits"); dev.off()
            }
        } else {
            pdf(file.path(contrast_out, "cog.pdf"))
            plot.new(); title("No COG annotations"); dev.off()
        }
    } else {
        pdf(file.path(contrast_out, "cog.pdf"))
        plot.new(); title("COG data unavailable"); dev.off()
    }

    # ── GO ORA enrichment ────────────────────────────────────
    if (length(up_genes) > 0 && !is.null(go_class)) {
        term2gene_go <- gene2go_all %>% dplyr::select(TERM=GO, GENE=GID)
        term2name_go <- go_class[, c("GO_ID", "Description")]
        colnames(term2name_go) <- c("TERM", "NAME")

        go_rich <- ora_hypergeometric(up_genes, term2gene_go, term2name_go,
                                       universe_genes)

        if (!is.null(go_rich) && nrow(go_rich) > 0) {
            write.table(go_rich, file.path(contrast_out, "go_rich.txt"),
                        sep="\t", row.names=FALSE, quote=FALSE)
            ggsave(plot_ora_bar(go_rich, paste("GO ORA —", contrast)),
                   file=file.path(contrast_out, "go_rich_bar.pdf"),
                   width=8, height=6)
            ggsave(plot_ora_dot(go_rich, paste("GO ORA —", contrast)),
                   file=file.path(contrast_out, "go_rich_dot.pdf"),
                   width=8, height=6)
            cat("    Saved: go_rich_bar.pdf, go_rich_dot.pdf\n")
        } else {
            pdf(file.path(contrast_out, "go_rich_bar.pdf"))
            plot.new(); title("No GO enrichment"); dev.off()
            pdf(file.path(contrast_out, "go_rich_dot.pdf"))
            plot.new(); title("No GO enrichment"); dev.off()
        }
    } else {
        pdf(file.path(contrast_out, "go_rich_bar.pdf"))
        plot.new(); title("No upregulated genes / GO.db unavailable"); dev.off()
        pdf(file.path(contrast_out, "go_rich_dot.pdf"))
        plot.new(); title("No upregulated genes / GO.db unavailable"); dev.off()
    }

    # ── KEGG ORA enrichment ──────────────────────────────────
    if (length(up_genes) > 0) {
        term2gene_kegg <- gene2pathway_all %>%
            dplyr::select(TERM=Pathway, GENE=GID)
        ko_rich <- ora_hypergeometric(up_genes, term2gene_kegg, NULL,
                                       universe_genes)

        if (!is.null(ko_rich) && nrow(ko_rich) > 0) {
            write.table(ko_rich, file.path(contrast_out, "ko_rich.txt"),
                        sep="\t", row.names=FALSE, quote=FALSE)
            ggsave(plot_ora_bar(ko_rich, paste("KEGG ORA —", contrast)),
                   file=file.path(contrast_out, "ko_rich_bar.pdf"),
                   width=8, height=6)
            ggsave(plot_ora_dot(ko_rich, paste("KEGG ORA —", contrast)),
                   file=file.path(contrast_out, "ko_rich_dot.pdf"),
                   width=8, height=6)
            cat("    Saved: ko_rich_bar.pdf, ko_rich_dot.pdf\n")
        } else {
            pdf(file.path(contrast_out, "ko_rich_bar.pdf"))
            plot.new(); title("No KEGG enrichment"); dev.off()
            pdf(file.path(contrast_out, "ko_rich_dot.pdf"))
            plot.new(); title("No KEGG enrichment"); dev.off()
        }
    } else {
        pdf(file.path(contrast_out, "ko_rich_bar.pdf"))
        plot.new(); title("No upregulated genes"); dev.off()
        pdf(file.path(contrast_out, "ko_rich_dot.pdf"))
        plot.new(); title("No upregulated genes"); dev.off()
    }
}

###############################################################
## SUMMARY
###############################################################
cat("\n==========================================\n")
cat("MeToReS Multi-Contrast Downstream Complete!\n")
cat("==========================================\n")
cat("Shared outputs:\n")
cat("  gene_pcoa.pdf, anosim_result.csv\n")
cat("  upset_all_DEGs.pdf, upset_upregulated.pdf\n")
cat("  all_contrast_DEG_counts.pdf, all_contrast_DEG_counts.csv\n")
cat("  kaiju_community_barplot.pdf, kaiju_community_table.csv\n")
cat("\nPer-contrast outputs (in per_contrast/{contrast}/):\n")
cat("  DEG_up_heatmap.pdf, DEG_down_heatmap.pdf\n")
cat("  cog.pdf\n")
cat("  go_rich_bar.pdf, go_rich_dot.pdf\n")
cat("  ko_rich_bar.pdf, ko_rich_dot.pdf\n")
cat("==========================================\n")
