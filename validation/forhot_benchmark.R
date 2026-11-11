#!/usr/bin/env Rscript
# forhot_benchmark.R — compare MeToReS on ForHot with Söllinger et al. 2022
# (Sci Adv 8:eabm3230; DataverseNO doi:10.18710/HZBZQE, CC0).
#
# KO level: Salmon counts of the genes that pass the MeToReS gene filter,
# summed per KEGG KO (eggNOG-mapper KEGG_ko; a gene with several KOs counts
# for each). DESeq2 ~ gradient + group (ET vs AT, paired by gradient),
# fold changes shrunk with ashr. Published contrast is A - E, so a published
# logFC > 0 means LOWER in warmed soil.
#
# Sections of the report:
#   1. Main comparison (the criteria written before the results were opened;
#      they were not frozen, so every version below is reported).
#   2. Sensitivity: size factors from the KO matrix or from Salmon mapped
#      reads, paired or unpaired design.
#   3. Ribosomal proteins relative to the other KOs (the paper reports
#      relative KO expression).
#   4. Assembly origin: each gene is named after the sample it was assembled
#      from; KO counts are rebuilt from subsets of genes by origin.
#
# Usage (from the MeToReS folder, metores environment):
#   Rscript validation/forhot_benchmark.R \
#     --results results_forhot --samplesheet examples/forhot/samplesheet_forhot.csv \
#     --published /path/to/HZBZQE --out validation/forhot_benchmark

suppressPackageStartupMessages({
  library(optparse); library(data.table); library(DESeq2)
})
o <- parse_args(OptionParser(option_list = list(
  make_option("--results"), make_option("--samplesheet"),
  make_option("--published"), make_option("--out"),
  make_option("--group_col", default = "group"),
  make_option("--min_count", type = "integer", default = 10),
  make_option("--alpha", type = "double", default = 0.05)
)))
dir.create(o$out, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(..., "\n", sep = "")
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x) formatC(x, format = "f", digits = 2)
fp <- function(p) if (is.na(p)) "n/a" else if (p < 1e-4) formatC(p, format = "e", digits = 0) else formatC(signif(p, 1), format = "fg")

# ── MeToReS inputs ──────────────────────────────────────────────────────
ss <- fread(o$samplesheet)
ss[, group := factor(get(o$group_col), levels = c("AT", "ET"))]
ss[, site := factor(site)]; ss[, gradient := factor(gradient)]
filt <- fread(file.path(o$results, "deg", paste0(o$group_col, "_deg_results"),
              "transcript_abundance_quantification_table_filter.csv"), select = 1)
genes <- filt[[1]]
say("Genes passing the MeToReS filter: ", length(genes))

cnt <- NULL
for (s in ss$sra) {
  q <- fread(file.path(o$results, "quant", paste0(s, "_quant"), "quant.sf"),
             select = c("Name", "NumReads"))
  q <- q[Name %chin% genes]
  setnames(q, "NumReads", s)
  cnt <- if (is.null(cnt)) q else merge(cnt, q, by = "Name", all = TRUE)
}
cnt[is.na(cnt)] <- 0

# Salmon mapped reads per sample (for library-size normalisation)
mapped <- sapply(ss$sra, function(s) {
  j <- readLines(file.path(o$results, "quant", paste0(s, "_quant"), "aux_info", "meta_info.json"), warn = FALSE)
  as.numeric(sub('.*"num_mapped": *([0-9]+).*', "\\1", grep('"num_mapped"', j, value = TRUE)[1]))
})
sf_mapped <- mapped / exp(mean(log(mapped)))

em_file <- Sys.glob(file.path(o$results, paste0("emapper_", o$group_col), "*.emapper.annotations"))[1]
em_lines <- readLines(em_file)
em_lines <- em_lines[!startsWith(em_lines, "##")]
em <- fread(text = paste(em_lines, collapse = "\n"), sep = "\t", quote = "")
setnames(em, 1, "query")
g2k <- em[KEGG_ko != "-" & KEGG_ko != "", .(KO = unlist(strsplit(KEGG_ko, ","))), by = query]
g2k[, KO := sub("^ko:", "", KO)]
g2k <- unique(g2k)
say("Filtered genes with a KO: ", uniqueN(g2k$query), "; distinct KOs: ", uniqueN(g2k$KO))

m <- merge(g2k, cnt, by.x = "query", by.y = "Name")

# assembly origin: gene ids start with the sample (sra) they were assembled from
sras <- ss$sra[order(-nchar(ss$sra))]
m[, origin := NA_character_]
for (s in sras) m[is.na(origin) & startsWith(query, paste0(s, "_")), origin := s]
m[, origin_site := as.character(ss$site)[match(origin, ss$sra)]]
m[, origin_group := as.character(ss$group)[match(origin, ss$sra)]]
fwrite(m[, c("query", "KO", "origin", "origin_site", "origin_group", ss$sra), with = FALSE],
       file.path(o$out, "ko_gene_counts.tsv"), sep = "\t")
say("Genes with an unknown origin: ", m[is.na(origin), uniqueN(query)])

ko_matrix <- function(rows = rep(TRUE, nrow(m))) {
  k <- m[rows, lapply(.SD, sum), by = KO, .SDcols = ss$sra]
  x <- round(as.matrix(k[, ss$sra, with = FALSE])); rownames(x) <- k$KO; x
}
kmat <- ko_matrix()
kmat <- kmat[rowSums(kmat) >= o$min_count, , drop = FALSE]
say("KOs quantified (>= ", o$min_count, " counts over all samples): ", nrow(kmat))
fwrite(data.table(KO = rownames(kmat), kmat), file.path(o$out, "ko_counts.tsv"), sep = "\t")

# DESeq2 on a KO matrix. sf: "poscounts" (from the KO matrix) or "mapped" (Salmon mapped reads).
fit_ko <- function(mat, idx, design, factor, num, den, sf = "poscounts",
                   fit = "parametric", min_rows = 1, label = NULL) {
  cd <- as.data.frame(ss[idx]); rownames(cd) <- cd$sra
  for (v in c("group", "site", "gradient")) cd[[v]] <- droplevels(cd[[v]])
  x  <- mat[, cd$sra, drop = FALSE]; x <- x[rowSums(x) >= min_rows, , drop = FALSE]
  dds <- suppressMessages(DESeqDataSetFromMatrix(x, cd, design))
  if (sf == "poscounts") dds <- estimateSizeFactors(dds, type = "poscounts")
  else sizeFactors(dds) <- sf_mapped[cd$sra]
  dds <- suppressMessages(DESeq(dds, quiet = TRUE, fitType = fit))
  res <- results(dds, contrast = c(factor, num, den), alpha = o$alpha)
  shr <- lfcShrink(dds, contrast = c(factor, num, den), type = "ashr", res = res, quiet = TRUE)
  out <- data.table(KO = rownames(shr), baseMean = shr$baseMean,
                    lfc_mle = res$log2FoldChange, lfc = shr$log2FoldChange,
                    pvalue = res$pvalue, padj = res$padj)
  if (!is.null(label)) {
    fwrite(out, file.path(o$out, paste0("ko_deseq2_", label, ".tsv")), sep = "\t")
    say("  ", label, ": ", nrow(out), " KOs tested, ", sum(out$padj < o$alpha, na.rm = TRUE),
        " at padj < ", o$alpha)
  }
  out
}
i_all <- seq_len(nrow(ss)); i_mtw <- which(ss$site == "MTW")
i_ltw <- which(ss$site == "LTW"); i_at <- which(ss$group == "AT")

say("\nKO-level DESeq2 (ashr):")
all_ <- fit_ko(kmat, i_all, ~ gradient + group, "group", "ET", "AT", label = "ET_vs_AT_both")
mtw  <- fit_ko(kmat, i_mtw, ~ gradient + group, "group", "ET", "AT", label = "ET_vs_AT_MTW")
ltw  <- fit_ko(kmat, i_ltw, ~ gradient + group, "group", "ET", "AT", label = "ET_vs_AT_LTW")
neg  <- fit_ko(kmat, i_at,  ~ site, "site", "MTW", "LTW", label = "AT_MTW_vs_LTW")

# ── published results ───────────────────────────────────────────────────
readpub <- function(f) {
  d <- suppressWarnings(fread(file.path(o$published, f)))   # files carry an extra row-index column
  if (!"K_nr" %in% names(d)) setnames(d, 2, "K_nr")      # first column is a row index
  unique(d[, .(KO = K_nr, logFC_AminusE = logFC, FDR)], by = "KO")
}
pub <- list(both = readpub("005a_DGE_Analysis_AvsE.txt"),
            MTW  = readpub("005b_DGE_Analysis_MTW_AvsE.txt"),
            LTW  = readpub("005c_DGE_Analysis_LTW_AvsE.txt"))
ribo <- fread(file.path(o$published, "006a_proteinBS_oxidP_mapping_file.txt"))
ribo <- ribo[complex == "ribosomal_proteins", unique(KO_number)]

direction <- function(p, mt) {
  p <- p[FDR < 0.05]
  j <- merge(p, mt, by = "KO")
  j[, same := sign(lfc) == -sign(logFC_AminusE)]
  list(n_pub = nrow(p), n_found = nrow(j), n_same = sum(j$same),
       pct = if (nrow(j)) 100 * mean(j$same) else NA_real_,
       rho = if (nrow(j) > 2) suppressWarnings(cor(j$lfc, -j$logFC_AminusE, method = "spearman")) else NA_real_,
       tab = j)
}
med_ribo <- function(d) { r <- d[KO %chin% ribo]; if (nrow(r)) median(r$lfc) else NA_real_ }
c3_hits  <- function(d) d[padj < o$alpha & abs(lfc) > 1, KO]

# ── 1. main comparison ──────────────────────────────────────────────────
c1  <- direction(pub$both, all_)
c1m <- direction(pub$MTW, mtw)
c1l <- direction(pub$LTW, ltw)
fwrite(c1$tab, file.path(o$out, "published_DE_KOs_vs_metores.tsv"), sep = "\t")

c2a <- med_ribo(all_) < 0; c2b <- med_ribo(mtw) < med_ribo(ltw)
neg_de <- length(c3_hits(neg)); neg_pct <- 100 * neg_de / nrow(neg)

pubDE  <- pub$both[FDR < 0.05, KO]
pub_in <- intersect(pubDE, all_$KO)
ourDE  <- all_[padj < o$alpha, KO]
recall <- if (length(pub_in)) 100 * length(intersect(pub_in, ourDE)) / length(pub_in) else NA_real_
jac    <- length(intersect(pubDE, ourDE)) / length(union(intersect(pubDE, all_$KO), ourDE))

pass <- function(x) if (isTRUE(x)) "PASS" else if (identical(x, FALSE)) "FAIL" else "n/a"
rep <- c(
 "# ForHot benchmark: MeToReS vs Söllinger et al. 2022", "",
 "The criteria below were written before the results were opened but were not frozen",
 "(not committed or tagged before the run), so every version computed is reported.", "",
 sprintf("Genes passing the MeToReS filter: %d; with a KO: %d; KOs quantified: %d.",
         length(genes), uniqueN(g2k$query), nrow(kmat)),
 sprintf("Published DE KOs (FDR < 0.05): both sites %d, MTW %d, LTW %d. Ribosomal-protein KOs in the published list: %d.",
         c1$n_pub, c1m$n_pub, c1l$n_pub, length(ribo)),
 sprintf("Share of Salmon mapped reads that end up in KO counts: %s-%s %% per sample.",
         f1(min(100 * colSums(kmat) / mapped)), f1(max(100 * colSums(kmat) / mapped))), "",
 "## 1. Main comparison", "",
 "| # | Criterion | Value | Result |", "|---|---|---|---|",
 sprintf("| C1 | Published DE KOs (both sites) with the same direction | %d of %d found (%s %%); Spearman ρ of fold changes %s | %s |",
         c1$n_same, c1$n_found, f1(c1$pct), f2(c1$rho), pass(c1$pct >= 80)),
 sprintf("| C1, MTW | same, per site | %d of %d (%s %%) | reported |", c1m$n_same, c1m$n_found, f1(c1m$pct)),
 sprintf("| C1, LTW | same, per site | %d of %d (%s %%) | reported |", c1l$n_same, c1l$n_found, f1(c1l$pct)),
 sprintf("| C2a | Ribosomal proteins, median shrunken log2FC ET vs AT | %s (n = %d) | %s |",
         f2(med_ribo(all_)), all_[KO %chin% ribo, .N], pass(c2a)),
 sprintf("| C2b | Median more negative at MTW than LTW | MTW %s (n = %d), LTW %s (n = %d) | %s |",
         f2(med_ribo(mtw)), mtw[KO %chin% ribo, .N], f2(med_ribo(ltw)), ltw[KO %chin% ribo, .N], pass(c2b)),
 sprintf("| C3 | AT: MTW vs LTW, KOs with padj < %s and |log2FC| > 1 | %d of %d (%s %%) | %s |",
         o$alpha, neg_de, nrow(neg), f1(neg_pct), pass(neg_pct <= 1)),
 sprintf("| C4 | Overlap with the published DE KOs | recall %s %% (%d of %d found); Jaccard %s | reported |",
         f1(recall), length(intersect(pub_in, ourDE)), length(pub_in), f2(jac)), "",
 sprintf("MeToReS KO-level DE (both sites, padj < %s): %d of %d KOs.", o$alpha, length(ourDE), nrow(all_)),
 "Published contrast A - E was inverted to ET vs AT before comparing directions."
)

# ── 2. sensitivity: normalisation and design ────────────────────────────
say("\nSensitivity runs:")
sens <- rbindlist(lapply(c("poscounts", "mapped"), function(sf) rbindlist(lapply(c("paired", "unpaired"), function(des) {
  if (sf == "poscounts" && des == "paired") { A <- all_; M <- mtw; L <- ltw; N <- neg } else {
    A <- fit_ko(kmat, i_all, if (des == "paired") ~ gradient + group else ~ site + group, "group", "ET", "AT", sf = sf)
    M <- fit_ko(kmat, i_mtw, if (des == "paired") ~ gradient + group else ~ group, "group", "ET", "AT", sf = sf)
    L <- fit_ko(kmat, i_ltw, if (des == "paired") ~ gradient + group else ~ group, "group", "ET", "AT", sf = sf)
    N <- if (sf == "poscounts") neg else fit_ko(kmat, i_at, ~ site, "site", "MTW", "LTW", sf = sf)
  }
  d1 <- direction(pub$both, A); dm <- direction(pub$MTW, M); dl <- direction(pub$LTW, L)
  say("  ", sf, " / ", des, " done")
  data.table(size_factors = if (sf == "poscounts") "KO matrix (poscounts)" else "Salmon mapped reads",
             design = if (des == "paired") "paired by gradient" else "unpaired",
             C1_both = sprintf("%d/%d", d1$n_same, d1$n_found),
             C1_MTW = sprintf("%d/%d", dm$n_same, dm$n_found), C1_LTW = sprintf("%d/%d", dl$n_same, dl$n_found),
             ribo_both = med_ribo(A), ribo_MTW = med_ribo(M), ribo_LTW = med_ribo(L),
             C3 = length(c3_hits(N)), C3_tested = nrow(N),
             KO_DE_both = sum(A$padj < o$alpha, na.rm = TRUE))
}))))
fwrite(sens, file.path(o$out, "sensitivity.tsv"), sep = "\t")
rep <- c(rep, "", "## 2. Sensitivity: normalisation and design", "",
 "Unpaired: ~ site + group over both sites, ~ group per site. C3 always uses ~ site on the 8 AT samples.", "",
 "| Size factors | Design | C1 both | C1 MTW | C1 LTW | Ribosomal median both / MTW / LTW | C3 | KO DE both |",
 "|---|---|---|---|---|---|---|---|",
 sens[, sprintf("| %s | %s | %s | %s | %s | %s / %s / %s | %d/%d | %d |", size_factors, design, C1_both, C1_MTW, C1_LTW,
                f2(ribo_both), f2(ribo_MTW), f2(ribo_LTW), C3, C3_tested, KO_DE_both)])

# ── 3. ribosomal proteins relative to the other KOs ─────────────────────
rel_ribo <- function(d) {
  r <- d[KO %chin% ribo & !is.na(lfc_mle), lfc_mle]; x <- d[!KO %chin% ribo & !is.na(lfc_mle), lfc_mle]
  list(n = length(r), ribo = median(r), other = median(x), diff = median(r) - median(x),
       p = if (length(r) > 1 && length(x) > 1) suppressWarnings(wilcox.test(r, x)$p.value) else NA_real_)
}
rel <- list(both = rel_ribo(all_), MTW = rel_ribo(mtw), LTW = rel_ribo(ltw))
rep <- c(rep, "", "## 3. Ribosomal proteins relative to the other KOs (main run)", "",
 "Median unshrunken log2FC ET vs AT; difference = ribosomal minus other KOs; Wilcoxon rank-sum test.", "",
 "| Sites | Ribosomal KOs | Ribosomal median | Other KOs median | Difference | P |", "|---|---|---|---|---|---|",
 sapply(names(rel), function(k) with(rel[[k]], sprintf("| %s | %d | %s | %s | %s | %s |", k, n, f2(ribo), f2(other), f2(diff), fp(p)))))

# ── 4. assembly origin ──────────────────────────────────────────────────
say("\nAssembly-origin checks:")
gpk <- m[, .(n = uniqueN(query)), by = KO][KO %chin% rownames(kmat)]
org <- function(rows, idx, design, factor, num, den)
  fit_ko(ko_matrix(rows), idx, design, factor, num, den, sf = "mapped", fit = "local", min_rows = o$min_count)
all_rows <- rep(TRUE, nrow(m))

# 4a. control sites, KO counts from genes assembled at each site
c3_main <- c3_hits(neg)
c3o <- rbindlist(lapply(c("all genes", "LTW-assembled", "MTW-assembled"), function(nm) {
  rows <- switch(nm, `all genes` = all_rows, `LTW-assembled` = m$origin_site %chin% "LTW",
                 `MTW-assembled` = m$origin_site %chin% "MTW")
  t <- org(rows, i_at, ~ site, "site", "MTW", "LTW")
  h <- c3_hits(t)
  data.table(genes = nm, KOs = nrow(t), hits = length(h), higher_MTW = t[KO %chin% h & lfc > 0, .N],
             med_c3 = median(t[KO %chin% c3_main, lfc]), n_c3 = t[KO %chin% c3_main, .N])
}))

# 4b. ET vs AT per site, KO counts from genes assembled in AT or ET samples
ro <- rbindlist(lapply(c("MTW", "LTW"), function(st) {
  idx <- if (st == "MTW") i_mtw else i_ltw
  rbindlist(lapply(c("all genes", "AT-assembled", "ET-assembled"), function(nm) {
    rows <- switch(nm, `all genes` = all_rows, `AT-assembled` = m$origin_group %chin% "AT",
                   `ET-assembled` = m$origin_group %chin% "ET")
    t <- org(rows, idx, ~ gradient + group, "group", "ET", "AT"); r <- rel_ribo(t)
    data.table(site = st, genes = nm, KOs = nrow(t), all_median = median(t$lfc_mle, na.rm = TRUE),
               ribo_n = r$n, ribo_median = r$ribo, diff = r$diff, p = r$p,
               ribo_DE = t[KO %chin% ribo & padj < o$alpha, .N])
  }))
}))
fwrite(c3o, file.path(o$out, "origin_control_sites.tsv"), sep = "\t")
fwrite(ro,  file.path(o$out, "origin_ribosomal.tsv"), sep = "\t")
gap <- ro[, .(gap = all_median[genes == "ET-assembled"] - all_median[genes == "AT-assembled"]), by = site]

rep <- c(rep, "", "## 4. Assembly origin", "",
 sprintf("Genes per KO: median %s, mean %s, maximum %d; %d of %d KOs rest on a single gene; ribosomal KOs: median %s genes.",
         f1(median(gpk$n)), f1(mean(gpk$n)), max(gpk$n), gpk[n == 1, .N], nrow(gpk), f1(median(gpk[KO %chin% ribo, n]))),
 sprintf("Genes by assembly site: LTW %d, MTW %d.", m[origin_site == "LTW", uniqueN(query)], m[origin_site == "MTW", uniqueN(query)]),
 "KO counts rebuilt from a subset of genes; size factors from Salmon mapped reads; KOs with at least",
 sprintf("%d counts in the samples of the contrast; DESeq2 fitType local.", o$min_count), "",
 sprintf("Origin gap (median log2FC ET vs AT over all KOs, ET-assembled minus AT-assembled genes): MTW %s, LTW %s.",
         f2(gap[site == "MTW", gap]), f2(gap[site == "LTW", gap])), "",
 "### 4a. Control sites (AT: MTW vs LTW)", "",
 sprintf("| Genes counted | KOs | Hits | Higher at MTW | Median shrunken log2FC of the %d main-run C3 KOs |", length(c3_main)),
 "|---|---|---|---|---|",
 c3o[, sprintf("| %s | %d | %d (%s %%) | %d | %s (n = %d) |", genes, KOs, hits, f1(100 * hits / KOs), higher_MTW, f2(med_c3), n_c3)],
 "", "### 4b. Ribosomal proteins relative to the other KOs, by assembly origin", "",
 "| Site | Genes counted | Ribosomal KOs | Ribosomal median | Difference to other KOs | P | Ribosomal KOs DE |",
 "|---|---|---|---|---|---|---|",
 ro[, sprintf("| %s | %s | %d | %s | %s | %s | %d |", site, genes, ribo_n, f2(ribo_median), f2(diff), sapply(p, fp), ribo_DE)])

writeLines(rep, file.path(o$out, "benchmark_report.md")); cat(paste(rep, collapse = "\n"), "\n")
