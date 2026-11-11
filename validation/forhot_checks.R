#!/usr/bin/env Rscript
# forhot_checks.R — checks asked for in Lorenzo's review of 04/11/2026
#   1.1  sample mapping: our per-sample KO profiles vs the published per-sample
#        KO counts (002a_KO_tax_table_raw.txt, read-based)
#   1.2  rRNA and taxonomy of the filtered genes; our reads vs the paper's
#   1.3  share of mapped reads that fall in the filtered genes
#
# Usage (from the MeToReS folder, metores environment):
#   Rscript validation/forhot_checks.R --results results_forhot \
#     --samplesheet examples/forhot/samplesheet_forhot.csv \
#     --published /path/to/HZBZQE \
#     --benchmark validation/forhot_benchmark \
#     --rrna_ids validation/forhot_checks/rrna_ids.txt \
#     --out validation/forhot_checks

suppressPackageStartupMessages({ library(optparse); library(data.table) })
o <- parse_args(OptionParser(option_list = list(
  make_option("--results"), make_option("--samplesheet"), make_option("--published"),
  make_option("--benchmark"), make_option("--rrna_ids"), make_option("--out"),
  make_option("--group_col", default = "group"),
  make_option("--fasta", default = NULL, help = "CDS of the filtered genes [default: <results>/emapper_<group_col>/filtered_cds.fasta]"),
  make_option("--min_ko_count", type = "integer", default = 5,
              help = "KOs used for sample matching: at least this many counts in every sample, both sides")
)))
dir.create(o$out, showWarnings = FALSE, recursive = TRUE)
rep <- character(); add <- function(...) { x <- unlist(list(...)); rep <<- c(rep, x); cat(x, sep = "\n"); cat("\n") }
f1 <- function(x) formatC(x, format = "f", digits = 1); f2 <- function(x) formatC(x, format = "f", digits = 2)

ss <- fread(o$samplesheet); S <- ss$sra
lab <- paste(ss$site, ss$group, sub(".*_", "", ss$gradient), sep = "_")   # e.g. MTW_AT_G1

# ── our data ────────────────────────────────────────────────────────────
filt  <- fread(file.path(o$results, "deg", paste0(o$group_col, "_deg_results"),
               "transcript_abundance_quantification_table_filter.csv"), select = 1)[[1]]
meta  <- function(s, key) { j <- readLines(file.path(o$results, "quant", paste0(s, "_quant"), "aux_info", "meta_info.json"), warn = FALSE)
  as.numeric(sub(sprintf('.*"%s": *([0-9]+).*', key), "\\1", grep(sprintf('"%s"', key), j, value = TRUE)[1])) }
mapped <- sapply(S, meta, key = "num_mapped"); processed <- sapply(S, meta, key = "num_processed")
cnt <- NULL
for (s in S) {
  q <- fread(file.path(o$results, "quant", paste0(s, "_quant"), "quant.sf"), select = c("Name", "NumReads"))
  q <- q[Name %chin% filt]; setnames(q, "NumReads", s)
  cnt <- if (is.null(cnt)) q else merge(cnt, q, by = "Name", all = TRUE)
}
cnt[is.na(cnt)] <- 0
gm <- as.matrix(cnt[, ..S]); rownames(gm) <- cnt$Name

em_file <- Sys.glob(file.path(o$results, paste0("emapper_", o$group_col), "*.emapper.annotations"))[1]
l  <- readLines(em_file); em <- fread(text = paste(l[!startsWith(l, "##")], collapse = "\n"), sep = "\t", quote = "")
setnames(em, 1, "query")
kt <- fread(file.path(o$results, "kaiju", "kaiju_gene_taxonomy.tsv"))
kt <- kt[gene_id %chin% filt]
rr <- if (!is.null(o$rrna_ids) && file.exists(o$rrna_ids)) unique(readLines(o$rrna_ids)) else character()

g <- data.table(gene = rownames(gm), total = rowSums(gm))
g[, rRNA := gene %chin% rr][, eggnog := gene %chin% em$query]
g[, KO := gene %chin% em[KEGG_ko != "-" & KEGG_ko != "", query]]
g[kt, on = .(gene = gene_id), `:=`(kaiju_status = i.status, superkingdom = i.superkingdom)]
g[is.na(superkingdom) | superkingdom == "", superkingdom := "unclassified"]
fwrite(g, file.path(o$out, "filtered_genes_checks.tsv"), sep = "\t")

add("# ForHot checks (Lorenzo's review, 04/11/2026)", "")
add("## 1.3 How much of the data the filtered genes hold", "")
add(sprintf("Filtered genes: %d. Share of Salmon mapped reads in the filtered genes: %s %% overall, %s-%s %% per sample.",
    nrow(gm), f1(100 * sum(gm) / sum(mapped)), f1(min(100 * colSums(gm) / mapped)), f1(max(100 * colSums(gm) / mapped))))

add("", "## 1.2 What the filtered genes are", "")
add(sprintf("rRNA (SortMeRNA hit on the gene sequence): %d of %d genes (%s %%), %s %% of the filtered-gene counts (per sample %s-%s %%).",
    sum(g$rRNA), nrow(g), f1(100 * mean(g$rRNA)), f1(100 * sum(g[rRNA == TRUE, total]) / sum(g$total)),
    f1(min(100 * colSums(gm[g$rRNA, , drop = FALSE]) / colSums(gm))), f1(max(100 * colSums(gm[g$rRNA, , drop = FALSE]) / colSums(gm)))))
add(sprintf("eggNOG annotation rate: all filtered genes %s %%; rRNA genes %s %%; other genes %s %%.",
    f1(100 * mean(g$eggnog)), f1(100 * mean(g[rRNA == TRUE, eggnog])), f1(100 * mean(g[rRNA == FALSE, eggnog]))))
sk <- g[, .(genes = .N, pct_genes = 100 * .N / nrow(g), pct_counts = 100 * sum(total) / sum(g$total),
            eggnog_pct = 100 * mean(eggnog), rRNA = sum(rRNA)), by = superkingdom][order(-genes)]
add("", "Kaiju superkingdom of the filtered genes:", "",
    "| Superkingdom | Genes | % genes | % counts | eggNOG annotated % | rRNA genes |", "|---|---|---|---|---|---|")
add(sk[, sprintf("| %s | %d | %s | %s | %s | %d |", superkingdom, genes, f1(pct_genes), f1(pct_counts), f1(eggnog_pct), rRNA)])
top <- g[order(-total)][1:min(50, .N)]
add("", sprintf("Top 50 filtered genes by counts: %d rRNA, %d with an eggNOG annotation, %d unclassified by Kaiju; together %s %% of the filtered-gene counts.",
    sum(top$rRNA), sum(top$eggnog), sum(top$superkingdom == "unclassified"), f1(100 * sum(top$total) / sum(g$total))))
add("", "Share of the filtered-gene counts by annotation:", "",
    sprintf("- eggNOG annotated: %s %%; with a KEGG KO: %s %%; rRNA hit: %s %%; neither eggNOG nor Kaiju nor rRNA: %s %%.",
      f1(100 * sum(g[eggnog == TRUE, total]) / sum(g$total)), f1(100 * sum(g[KO == TRUE, total]) / sum(g$total)),
      f1(100 * sum(g[rRNA == TRUE, total]) / sum(g$total)),
      f1(100 * sum(g[eggnog == FALSE & rRNA == FALSE & superkingdom == "unclassified", total]) / sum(g$total))))

# top genes by counts, with length, GC and annotation
fa <- if (!is.null(o$fasta)) o$fasta else file.path(o$results, paste0("emapper_", o$group_col), "filtered_cds.fasta")
seqlen <- gc <- setNames(numeric(), character())
if (file.exists(fa)) {
  x <- readLines(fa); h <- startsWith(x, ">"); id <- sub(" .*", "", sub("^>", "", x[h]))
  sq <- vapply(split(x[!h], cumsum(h)[!h]), paste, "", collapse = "")
  seqlen <- setNames(nchar(sq), id); gc <- setNames(100 * nchar(gsub("[^GCgc]", "", sq)) / nchar(sq), id)
}
ann <- em[, .(gene = query, Preferred_name, Description, KEGG_ko)]
tg <- merge(g[order(-total)][1:min(200, .N)], ann, by = "gene", all.x = TRUE)[order(-total)]
tg[kt, on = .(gene = gene_id), phylum := i.phylum]
tg[, `:=`(pct_counts = 100 * total / sum(g$total), length = seqlen[gene], GC = gc[gene],
          samples_ge5 = rowSums(gm[gene, , drop = FALSE] >= 5), origin = sub("_k141_.*", "", gene))]
fwrite(tg, file.path(o$out, "top_genes.tsv"), sep = "\t")
add("", "Top 30 filtered genes by counts (all 200 in top_genes.tsv):", "",
    "| Gene | % of counts | Length | GC % | Samples ≥ 5 | rRNA | Kaiju | eggNOG |", "|---|---|---|---|---|---|---|---|")
add(tg[1:min(30, .N), sprintf("| %s | %s | %s | %s | %d | %s | %s | %s |", gene, f2(pct_counts),
       ifelse(is.na(length), "", as.character(length)), ifelse(is.na(GC), "", f1(GC)), samples_ge5, ifelse(rRNA, "yes", ""),
       ifelse(superkingdom == "unclassified", "", paste(superkingdom, phylum)),
       ifelse(is.na(Preferred_name), "", paste(Preferred_name, substr(Description, 1, 40))))])

# ── published per-sample KO counts (read-based) ─────────────────────────
# written by R write.table: quoted header with one name fewer than the data (row names first)
raw <- suppressWarnings(fread(file.path(o$published, "002a_KO_tax_table_raw.txt"), sep = "\t"))
if (ncol(raw) == 27) raw[, (1) := NULL]
setnames(raw, gsub('^"|"$', "", names(raw)))
kcol <- grep("^K_?number$|^K_nr$", names(raw), value = TRUE)[1]
pcols <- names(raw)[11:26]                       # 002_data_preparation.Rmd: colSums(KO_tax_raw[,11:26])
if (is.na(kcol) || !all(sapply(raw[, ..pcols], is.numeric)))
  stop("Unexpected layout of 002a. Columns: ", paste(names(raw), collapse = ", "))
add("", sprintf("Published table: %d rows; KO column '%s'; other columns 1-10: %s; sample columns 11-26: %s.",
    nrow(raw), kcol, paste(names(raw)[1:10], collapse = ", "), paste(pcols, collapse = ", ")))
for (cc in setdiff(names(raw)[1:10], kcol)) {
  v <- raw[[cc]]
  if (is.character(v) && any(v == "d__Eukaryota", na.rm = TRUE))
    add(sprintf("Published KO-annotated reads from Eukaryota (column '%s'): %s %% overall, %s-%s %% per sample.", cc,
        f1(100 * sum(as.matrix(raw[v == "d__Eukaryota", ..pcols])) / sum(as.matrix(raw[, ..pcols]))),
        f1(min(100 * colSums(raw[v == "d__Eukaryota", ..pcols]) / colSums(raw[, ..pcols]))),
        f1(max(100 * colSums(raw[v == "d__Eukaryota", ..pcols]) / colSums(raw[, ..pcols])))))
}
# order in the published scripts: GO-A 1-4, GO-E 1-4, GN-A 1,2,4,5, GN-E 1,2,4,5 (GO = LTW, GN = MTW);
# provisional labels in our naming (GN 1,2,4,5 -> G1-G4), to be confirmed by the matching below
plab <- c(paste0("LTW_AT_G", 1:4), paste0("LTW_ET_G", 1:4), paste0("MTW_AT_G", 1:4), paste0("MTW_ET_G", 1:4))
pk <- raw[, lapply(.SD, sum), by = c(kcol), .SDcols = pcols]; setnames(pk, kcol, "KO")
pm <- as.matrix(pk[, ..pcols]); rownames(pm) <- pk$KO; colnames(pm) <- plab
ourlab <- setNames(lab, S)

add("", "### Reads per sample", "",
    "| Sample | Our reads after rRNA removal (Salmon processed) | Our mapped reads | Published KO-annotated reads | Ours / published |",
    "|---|---|---|---|---|")
for (i in seq_along(S)) { p <- sum(pm[, ourlab[S[i]]])
  add(sprintf("| %s | %s | %s | %s | %s |", ourlab[S[i]], format(processed[i], big.mark = ","), format(mapped[i], big.mark = ","),
              format(round(p), big.mark = ","), f1(processed[i] / p))) }

# ── 1.1 sample matching ─────────────────────────────────────────────────
add("", "## 1.1 Sample mapping", "")
ko <- fread(file.path(o$benchmark, "ko_counts.tsv")); om <- as.matrix(ko[, ..S]); rownames(om) <- ko$KO
colnames(om) <- ourlab[colnames(om)]
common <- intersect(rownames(om), rownames(pm))
use <- common[apply(om[common, , drop = FALSE], 1, min) >= o$min_ko_count & apply(pm[common, , drop = FALSE], 1, min) >= o$min_ko_count]
add(sprintf("KOs in both tables: %d; used for matching (>= %d counts in every sample on both sides): %d.", length(common), o$min_ko_count, length(use)))
lc <- function(m) { m <- m[use, , drop = FALSE]; log2(t(t(m) / colSums(m)) * 1e6 + 1) }
O <- lc(om); P <- lc(pm[, colnames(om)])
centre <- function(m, cols) m[, cols, drop = FALSE] - rowMeans(m[, cols, drop = FALSE])
match_tab <- function(cols, title) {
  r <- cor(centre(O, cols), centre(P, cols))
  fwrite(data.table(ours = rownames(r), round(r, 3)), file.path(o$out, paste0("sample_match_", gsub("\\W+", "_", title), ".tsv")), sep = "\t")
  add("", sprintf("**%s** (KO profiles centred within this set; Pearson r of ours vs published):", title), "",
      "| Our sample | Best published match | r | r with its own label | Rank of own label |", "|---|---|---|---|---|")
  for (s in rownames(r)) { x <- sort(r[s, ], decreasing = TRUE)
    add(sprintf("| %s | %s | %s | %s | %d |", s, names(x)[1], f2(x[1]), f2(r[s, s]), match(s, names(x)))) }
  add(sprintf("Own label is the best match for %d of %d samples.", sum(sapply(rownames(r), function(s) names(which.max(r[s, ])) == s)), nrow(r)))
}
match_tab(colnames(O), "All 16 samples")
match_tab(grep("^MTW", colnames(O), value = TRUE), "MTW only")
match_tab(grep("^LTW", colnames(O), value = TRUE), "LTW only")

# ribosomal share per sample, ours vs published (all KOs, not only the matching set)
ribo <- fread(file.path(o$published, "006a_proteinBS_oxidP_mapping_file.txt"))[complex == "ribosomal_proteins", unique(KO_number)]
rs_o <- 100 * colSums(om[rownames(om) %chin% ribo, , drop = FALSE]) / colSums(om)
rs_p <- 100 * colSums(pm[rownames(pm) %chin% ribo, , drop = FALSE]) / colSums(pm)
rs <- data.table(sample = names(rs_o), ours = rs_o, published = rs_p[names(rs_o)])
add("", "### Ribosomal-protein share of KO counts per sample (%)", "", "| Sample | Ours | Published |", "|---|---|---|")
add(rs[, sprintf("| %s | %s | %s |", sample, f1(ours), f1(published))])
add(sprintf("Spearman correlation over the 16 samples: %s; within MTW: %s; within LTW: %s.",
    f2(cor(rs$ours, rs$published, method = "spearman")),
    f2(rs[grepl("^MTW", sample), cor(ours, published, method = "spearman")]),
    f2(rs[grepl("^LTW", sample), cor(ours, published, method = "spearman")])))

writeLines(rep, file.path(o$out, "checks_report.md"))
