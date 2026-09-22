#!/usr/bin/env Rscript
###############################################################
## constrained_ordination.R — db-RDA / distance-based RDA
##
## Constrained ordination using vegan::capscale
## (distance-based RDA) on VST-normalized counts.
##
## Produces:
##   - Single-variable db-RDA per variable (one constraint each)
##   - Full model db-RDA (~ all variables)
##   - Optional conditioned models (test | Condition(control))
##   - Variance partitioning across variables (varpart)
##   - Permutation significance tests (anova.cca)
##
## Usage:
##   Rscript constrained_ordination.R \
##     --vst        results/deg/country_deg_results/vst_normalised_counts.csv \
##     --samplesheet samplesheet_mobiles.csv \
##     --workdir    results/constrained_ordination_country/ \
##     --variables  "contaminant season country" \
##     --distance   bray \
##     --conditioned "contaminant:season contaminant:country"
###############################################################

suppressWarnings(suppressMessages(library(optparse)))

option_list <- list(
    make_option(c("-v", "--vst"),        type="character",
        help="VST-normalized counts CSV (genes x samples, first col = GeneID) [required]"),
    make_option(c("-s", "--samplesheet"), type="character",
        help="Samplesheet CSV with sample + metadata columns [required]"),
    make_option(c("-w", "--workdir"),    type="character",
        help="Output directory [required]"),
    make_option(c("--variables"),        type="character", default=NULL,
        help="Space-separated metadata columns to use as constraints [required]"),
    make_option(c("--distance"),         type="character", default="bray",
        help="Distance metric (bray, jaccard, euclidean) [default: bray]"),
    make_option(c("--conditioned"),      type="character", default=NULL,
        help="Optional space-separated test:control pairs, e.g. 'contaminant:season contaminant:country'"),
    make_option(c("--permutations"),     type="integer", default=999,
        help="Number of permutations for anova.cca [default: 999]")
)

opts <- parse_args(OptionParser(option_list=option_list),
    args=commandArgs(trailingOnly=TRUE))

for (req in c("vst", "samplesheet", "workdir", "variables")) {
    if (is.null(opts[[req]])) stop(paste0("--", req, " is required"))
}

VARS <- strsplit(trimws(opts$variables), "\\s+")[[1]]
VARS <- VARS[nchar(VARS) > 0]
if (length(VARS) == 0) stop("No variables parsed from --variables")

# Parse conditioned pairs (optional)
COND_PAIRS <- list()
if (!is.null(opts$conditioned) && nchar(trimws(opts$conditioned)) > 0) {
    pairs <- strsplit(trimws(opts$conditioned), "\\s+")[[1]]
    pairs <- pairs[nchar(pairs) > 0]
    for (p in pairs) {
        parts <- strsplit(p, ":")[[1]]
        if (length(parts) == 2) {
            COND_PAIRS[[length(COND_PAIRS) + 1]] <- list(test=parts[1], cond=parts[2])
        } else {
            cat("  WARNING: ignoring malformed conditioned pair:", p, "\n")
        }
    }
}

opts$vst         <- normalizePath(opts$vst, mustWork=FALSE)
opts$samplesheet <- normalizePath(opts$samplesheet, mustWork=FALSE)
opts$workdir     <- normalizePath(opts$workdir, mustWork=FALSE)
dir.create(opts$workdir, showWarnings=FALSE, recursive=TRUE)
setwd(opts$workdir)

cat("========================================\n")
cat("MetoReS Constrained Ordination (db-RDA)\n")
cat("========================================\n")
cat("VST matrix  :", opts$vst, "\n")
cat("Variables   :", paste(VARS, collapse=", "), "\n")
cat("Distance    :", opts$distance, "\n")
cat("Conditioned :", length(COND_PAIRS), "model(s)\n")
cat("========================================\n\n")

suppressWarnings(suppressMessages({
    library(vegan)
    library(ggplot2)
    library(dplyr)
}))

# =============================================================================
# Load data
# =============================================================================
if (!file.exists(opts$vst))
    stop(paste("VST matrix not found:", opts$vst))

vst <- read.csv(opts$vst, header=TRUE, stringsAsFactors=FALSE,
                check.names=FALSE)
# First column is GeneID
rownames(vst) <- vst[[1]]
vst <- vst[, -1, drop=FALSE]

# Samples are columns; transpose so samples are rows (community matrix form)
comm <- t(as.matrix(vst))
cat("Community matrix:", nrow(comm), "samples x", ncol(comm), "genes\n")

ss <- read.csv(opts$samplesheet, stringsAsFactors=FALSE)
if (!"sample" %in% colnames(ss))
    stop("Samplesheet must have a 'sample' column")

# Match metadata rows to community matrix rows (by sample name)
rownames(ss) <- ss$sample
common_samples <- intersect(rownames(comm), rownames(ss))
if (length(common_samples) < nrow(comm)) {
    cat("  WARNING: only", length(common_samples), "of", nrow(comm),
        "samples matched the samplesheet\n")
}
comm <- comm[common_samples, , drop=FALSE]
meta <- ss[common_samples, , drop=FALSE]

# Validate requested variables
missing_vars <- setdiff(VARS, colnames(meta))
if (length(missing_vars) > 0)
    stop(paste("Variables not in samplesheet:", paste(missing_vars, collapse=", ")))

# Coerce constraint variables to factors
for (v in VARS) meta[[v]] <- factor(meta[[v]])

cat("\nMetadata levels:\n")
for (v in VARS) {
    cat("  ", v, ":", paste(levels(meta[[v]]), collapse=", "), "\n")
}
cat("\n")

# =============================================================================
# Helper: build distance matrix
# =============================================================================
# Drop zero-variance genes (cause problems with some distances)
comm <- comm[, apply(comm, 2, var) > 0, drop=FALSE]
cat("Genes after dropping zero-variance:", ncol(comm), "\n")

# Bray-Curtis works on non-negative data. VST output can be slightly negative
# for very low counts. Shift to non-negative for Bray-Curtis only.
comm_for_dist <- comm
if (opts$distance == "bray" && min(comm_for_dist) < 0) {
    shift <- abs(min(comm_for_dist))
    comm_for_dist <- comm_for_dist + shift
    cat("  Note: VST had negative values; shifted by", round(shift, 3),
        "for Bray-Curtis\n")
}

dist_mat <- vegdist(comm_for_dist, method=opts$distance)

# =============================================================================
# Helper: run one db-RDA, save ordination plot + return stats
# =============================================================================
results_summary <- data.frame(
    Model            = character(),
    Constraint       = character(),
    Constrained_var  = numeric(),
    Total_var        = numeric(),
    Pct_constrained  = numeric(),
    F_stat           = numeric(),
    P_value          = numeric(),
    stringsAsFactors = FALSE
)

run_dbrda <- function(formula_rhs, model_name, plot_title, color_var=NULL) {
    cat("\n--- Model:", model_name, "---\n")
    cat("  Formula: dist ~", formula_rhs, "\n")

    # capscale needs the formula with the distance matrix on the LHS
    form <- as.formula(paste("dist_mat ~", formula_rhs))

    mod <- tryCatch(
        capscale(form, data=meta),
        error = function(e) {
            cat("  ERROR fitting model:", conditionMessage(e), "\n")
            return(NULL)
        }
    )
    if (is.null(mod)) return(NULL)

    # Variance explained
    total_inertia       <- mod$tot.chi
    constrained_inertia <- mod$CCA$tot.chi
    pct <- round(100 * constrained_inertia / total_inertia, 2)
    cat("  Constrained variance:", pct, "%\n")

    # Permutation test (overall)
    set.seed(42)
    av <- tryCatch(
        anova.cca(mod, permutations=opts$permutations),
        error = function(e) NULL
    )
    f_stat <- if (!is.null(av)) av$F[1] else NA
    p_val  <- if (!is.null(av)) av$`Pr(>F)`[1] else NA
    cat("  Permutation test: F =", round(f_stat, 3),
        ", p =", p_val, "\n")

    # Build ordination plot (samples in constrained space)
    site_scores <- as.data.frame(scores(mod, display="sites", choices=c(1, 2)))
    colnames(site_scores) <- c("Dim1", "Dim2")
    site_scores$sample <- rownames(site_scores)

    if (!is.null(color_var) && color_var %in% colnames(meta)) {
        site_scores$grp <- meta[rownames(site_scores), color_var]
    } else {
        site_scores$grp <- meta[rownames(site_scores), VARS[1]]
    }

    # Axis labels with % variance
    eig <- mod$CCA$eig
    if (!is.null(eig) && length(eig) >= 2) {
        ax1_pct <- round(100 * eig[1] / sum(c(mod$CCA$eig, mod$CA$eig)), 1)
        ax2_pct <- round(100 * eig[2] / sum(c(mod$CCA$eig, mod$CA$eig)), 1)
        xlab_txt <- paste0("dbRDA1 (", ax1_pct, "%)")
        ylab_txt <- paste0("dbRDA2 (", ax2_pct, "%)")
    } else {
        xlab_txt <- "dbRDA1"; ylab_txt <- "dbRDA2"
    }

    p <- ggplot(site_scores, aes(Dim1, Dim2, colour=grp)) +
        geom_point(size=4) +
        geom_text(aes(label=sample), vjust=-1, size=2.8, show.legend=FALSE) +
        geom_vline(xintercept=0, linetype="dotted") +
        geom_hline(yintercept=0, linetype="dotted") +
        labs(title=plot_title,
             subtitle=paste0("Constrained variance: ", pct, "%  |  p = ", p_val),
             x=xlab_txt, y=ylab_txt, colour="") +
        theme_bw(base_size=13) +
        theme(plot.title=element_text(face="bold", hjust=0.5),
              plot.subtitle=element_text(hjust=0.5, size=9, colour="grey40"))

    outfile <- paste0("dbrda_", model_name, ".pdf")
    ggsave(outfile, p, width=8, height=6.5, device=cairo_pdf)
    cat("  Saved:", outfile, "\n")

    results_summary <<- rbind(results_summary, data.frame(
        Model           = model_name,
        Constraint      = formula_rhs,
        Constrained_var = round(constrained_inertia, 4),
        Total_var       = round(total_inertia, 4),
        Pct_constrained = pct,
        F_stat          = round(f_stat, 4),
        P_value         = p_val,
        stringsAsFactors = FALSE
    ))

    mod
}

# =============================================================================
# 1. Single-variable models
# =============================================================================
cat("=== Single-variable db-RDA models ===\n")
for (v in VARS) {
    run_dbrda(v, v,
              paste0("db-RDA constrained by ", v),
              color_var=v)
}

# =============================================================================
# 2. Full model (all variables)
# =============================================================================
cat("\n=== Full model ===\n")
full_rhs <- paste(VARS, collapse=" + ")
run_dbrda(full_rhs, "full_model",
          paste0("db-RDA — full model (", full_rhs, ")"),
          color_var=VARS[1])

# Term-wise significance for the full model (by margin)
cat("\n--- Full model: term-wise permutation test ---\n")
full_form <- as.formula(paste("dist_mat ~", full_rhs))
full_mod  <- tryCatch(capscale(full_form, data=meta), error=function(e) NULL)
if (!is.null(full_mod)) {
    set.seed(42)
    term_av <- tryCatch(
        anova.cca(full_mod, by="margin", permutations=opts$permutations),
        error = function(e) NULL
    )
    if (!is.null(term_av)) {
        term_df <- as.data.frame(term_av)
        term_df$Term <- rownames(term_df)
        write.csv(term_df, "full_model_termwise_anova.csv", row.names=FALSE)
        cat("  Saved: full_model_termwise_anova.csv\n")
        print(term_av)
    }
}

# =============================================================================
# 3. Conditioned models (optional)
# =============================================================================
if (length(COND_PAIRS) > 0) {
    cat("\n=== Conditioned models ===\n")
    for (cp in COND_PAIRS) {
        test_v <- cp$test
        cond_v <- cp$cond
        if (!test_v %in% VARS || !cond_v %in% colnames(meta)) {
            cat("  Skipping", test_v, "| Condition(", cond_v, ") — variable missing\n")
            next
        }
        rhs <- paste0(test_v, " + Condition(", cond_v, ")")
        model_name <- paste0(test_v, "_cond_", cond_v)
        run_dbrda(rhs, model_name,
                  paste0("db-RDA: ", test_v, " | controlling for ", cond_v),
                  color_var=test_v)
    }
}

# =============================================================================
# 4. Variance partitioning
# =============================================================================
if (length(VARS) >= 2 && length(VARS) <= 4) {
    cat("\n=== Variance partitioning ===\n")
    # Each explanatory table is a one-term formula EVALUATED IN meta.
    # Passing data=meta (and evaluating the formulas against it) is essential:
    # a bare ~contaminant formula otherwise carries the global environment,
    # where the factor does not exist, giving "object 'contaminant' not found".
    varpart_formulas <- lapply(VARS, function(v) {
        stats::as.formula(paste("~", v), env = environment())
    })

    vp <- tryCatch(
        do.call(varpart, c(list(Y = dist_mat), varpart_formulas,
                           list(data = meta))),
        error = function(e) {
            cat("  ERROR in varpart:", conditionMessage(e), "\n")
            NULL
        }
    )

    if (!is.null(vp)) {
        pdf("variance_partitioning.pdf", width=8, height=7)
        plot(vp, bg=c("#E64B35", "#4DBBD5", "#00A087", "#3C5488")[seq_along(VARS)],
             Xnames=VARS, digits=2)
        title(main="Variance Partitioning (db-RDA)")
        dev.off()
        cat("  Saved: variance_partitioning.pdf\n")

        # Save the fractions table
        vp_frac <- vp$part$indfract

        # vegan's varpart() labels fractions generically as X1, X2, X3, ...
        # in the order the tables were passed to varpart() — which is exactly
        # the order of VARS. Substitute the real variable names so the CSV is
        # human-readable without needing to cross-reference the config.
        frac_labels <- rownames(vp_frac)
        for (idx in seq_along(VARS)) {
            frac_labels <- gsub(paste0("X", idx, "\\b"), VARS[idx], frac_labels)
        }

        write.csv(data.frame(Fraction=frac_labels, vp_frac),
                  "variance_partitioning_fractions.csv", row.names=FALSE)
        cat("  Saved: variance_partitioning_fractions.csv\n")
    }
} else {
    cat("\n  Variance partitioning skipped (needs 2-4 variables, have ",
        length(VARS), ")\n")
}

# =============================================================================
# Save summary
# =============================================================================
write.csv(results_summary, "dbrda_anova_results.csv", row.names=FALSE)
write.csv(results_summary[, c("Model", "Pct_constrained", "P_value")],
          "variance_explained.csv", row.names=FALSE)

cat("\n========================================\n")
cat("Constrained Ordination Complete!\n")
cat("========================================\n")
print(results_summary, row.names=FALSE)
cat("\nOutput directory:", getwd(), "\n")
cat("========================================\n")
