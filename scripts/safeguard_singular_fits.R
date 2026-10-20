#!/usr/bin/env Rscript
# safeguard_singular_fits.R — S3: singular random-effect fits behind MaAsLin2 hits.
#
# With about two samples per level of a random effect (e.g. site), lme4
# often estimates the random-effect variance at exactly zero ("singular
# fit"). The model then silently becomes an ordinary linear model without
# that grouping, and any difference between groups passes straight into the
# covariate's coefficient. If singular fits are much more common among the
# significant associations than among all tests, the hits are likely to be
# group (site) differences rather than covariate effects.
#
# This script refits every MaAsLin2 model with the same formula
#     expr ~ element + fixed covariates + (1 | random effect) ...
# records whether each fit is singular, and compares the singular fraction
# among hits with the fraction among all tests. It first checks that the
# refit reproduces MaAsLin2's p-values; if it does not, the singular-fit
# figures are reported but marked unreliable.
#
# Optional (--within_between): refit with each element split into its
# group mean (between) and the deviation from it (within), using the first
# random effect as the grouping. An association supported only by the
# between part is a group difference, whatever its cause.
#
# Report only: MaAsLin2's results are not changed.

suppressPackageStartupMessages({
  library(optparse); library(lme4); library(lmerTest); library(parallel)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--maaslin_all",      type = "character"),
  make_option("--vst",              type = "character"),
  make_option("--chemical",         type = "character"),
  make_option("--samplesheet",      type = "character"),
  make_option("--random_effects",   type = "character", default = ""),
  make_option("--fixed_covariates", type = "character", default = ""),
  make_option("--transform",        type = "character", default = "log10_zscore"),
  make_option("--alpha",            type = "double",    default = 0.05),
  make_option("--warn_ratio",       type = "double",    default = 2.0),
  make_option("--within_between",   action = "store_true", default = FALSE),
  make_option("--cores",            type = "integer",   default = 1),
  make_option("--out_tsv",          type = "character"),
  make_option("--out_json",         type = "character")
)))
for (req in c("maaslin_all", "vst", "chemical", "samplesheet", "out_tsv", "out_json"))
  if (is.null(opt[[req]])) stop("--", req, " is required")

split_ws <- function(x) { v <- strsplit(trimws(x), "\\s+")[[1]]; v[nchar(v) > 0] }
RANDOM <- split_ws(opt$random_effects)
FIXCOV <- split_ws(opt$fixed_covariates)
dir.create(dirname(opt$out_tsv), showWarnings = FALSE, recursive = TRUE)

write_json <- function(lst, path) {
  # minimal JSON writer (avoids a jsonlite dependency)
  enc <- function(v) {
    if (is.null(v) || (length(v) == 1 && is.na(v))) return("null")
    if (is.list(v)) {
      if (!is.null(names(v)))
        return(paste0("{", paste0('"', names(v), '": ', vapply(v, enc, ""), collapse = ", "), "}"))
      return(paste0("[", paste(vapply(v, enc, ""), collapse = ", "), "]"))
    }
    if (is.logical(v)) return(tolower(as.character(v)))
    if (is.numeric(v)) return(format(v, digits = 6))
    paste0('"', gsub('"', '\\\\"', as.character(v)), '"')
  }
  writeLines(enc(lst), path)
}

if (!length(RANDOM)) {
  cat("S3: MaAsLin2 was run without random effects; there is nothing to check.\n")
  write.table(data.frame(), opt$out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  write_json(list(safeguard = "S3_singular_fits", applicable = FALSE,
                  note = "no random effects in the MaAsLin2 model"), opt$out_json)
  quit(save = "no")
}

ml <- read.delim(opt$maaslin_all, stringsAsFactors = FALSE, check.names = FALSE)
qcol <- if ("qval_across_elements" %in% names(ml)) "qval_across_elements" else "qval"
ml$is_hit <- !is.na(ml[[qcol]]) & ml[[qcol]] < opt$alpha

vst <- read.csv(opt$vst, check.names = FALSE)
rownames(vst) <- vst[[1]]; vst[[1]] <- NULL
chem <- read.csv(opt$chemical, check.names = FALSE, stringsAsFactors = FALSE)
rownames(chem) <- chem$SampleID
ss <- read.csv(opt$samplesheet, stringsAsFactors = FALSE, check.names = FALSE)
rownames(ss) <- ss$sample
samples <- Reduce(intersect, list(colnames(vst), rownames(chem), ss$sample))

# Same transform as run_maaslin2.R, then MaAsLin2's own z-scoring of
# continuous metadata.
tf <- function(x) {
  x <- as.numeric(x)
  if (opt$transform %in% c("log10", "log10_zscore")) x <- log10(x + 1)
  x
}
z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

fit_one <- function(i) {
  g <- ml$feature[i]; el <- ml$element[i]
  out <- c(refit_coef = NA, refit_p = NA, singular = NA,
           within_coef = NA, within_p = NA, between_p = NA)
  if (!g %in% rownames(vst) || !el %in% colnames(chem)) return(out)
  d <- data.frame(expr = as.numeric(vst[g, samples]), x = z(tf(chem[samples, el])))
  for (cv in c(FIXCOV, RANDOM)) d[[cv]] <- ss[samples, cv]
  for (cv in FIXCOV) if (is.numeric(d[[cv]])) d[[cv]] <- z(d[[cv]])
  d <- d[complete.cases(d), , drop = FALSE]
  if (nrow(d) < 4 || var(d$expr) == 0) return(out)
  rhs <- paste(c("x", FIXCOV, sprintf("(1 | %s)", RANDOM)), collapse = " + ")
  f <- tryCatch(suppressMessages(suppressWarnings(
         lmerTest::lmer(as.formula(paste("expr ~", rhs)), data = d))), error = function(e) NULL)
  if (!is.null(f)) {
    cf <- summary(f)$coefficients
    out["refit_coef"] <- cf["x", "Estimate"]; out["refit_p"] <- cf["x", "Pr(>|t|)"]
    out["singular"] <- isSingular(f)
  }
  if (opt$within_between) {
    grp <- d[[RANDOM[1]]]
    d$x_between <- ave(d$x, grp); d$x_within <- d$x - d$x_between
    if (sum(abs(d$x_within) > 1e-12) >= 2) {
      rhs2 <- paste(c("x_within", "x_between", FIXCOV, sprintf("(1 | %s)", RANDOM)), collapse = " + ")
      f2 <- tryCatch(suppressMessages(suppressWarnings(
              lmerTest::lmer(as.formula(paste("expr ~", rhs2)), data = d))), error = function(e) NULL)
      if (!is.null(f2)) {
        cf <- summary(f2)$coefficients
        if ("x_within" %in% rownames(cf)) {
          out["within_coef"] <- cf["x_within", "Estimate"]; out["within_p"] <- cf["x_within", "Pr(>|t|)"] }
        if ("x_between" %in% rownames(cf)) out["between_p"] <- cf["x_between", "Pr(>|t|)"]
      }
    }
  }
  out
}

cat("S3: refitting", nrow(ml), "MaAsLin2 models (", length(unique(ml$element)),
    "elements ) on", opt$cores, "cores...\n")
t0 <- Sys.time()
r <- do.call(rbind, mclapply(seq_len(nrow(ml)), fit_one, mc.cores = opt$cores))
cat("  done in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

res <- cbind(ml[, intersect(c("feature", "element", "coef", "pval", qcol), names(ml))],
             is_hit = ml$is_hit, as.data.frame(r))
res$singular <- as.logical(res$singular)
if (opt$within_between) {
  res$within_q <- NA
  k <- !is.na(res$within_p); res$within_q[k] <- p.adjust(res$within_p[k], "BH")
}
write.table(res, opt$out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

# ── Fidelity check: does the refit reproduce MaAsLin2? ─────────────────
ok <- !is.na(res$refit_p) & !is.na(res$pval)
rho <- if (sum(ok) > 10) suppressWarnings(cor(res$refit_p[ok], res$pval[ok], method = "spearman")) else NA
reliable <- !is.na(rho) && rho > 0.99
cat(sprintf("  refit vs MaAsLin2 p-values: Spearman rho = %.4f over %d tests%s\n", rho, sum(ok),
            if (reliable) "" else "  -> NOT reproduced; singular-fit figures unreliable"))

# ── Singular fractions ─────────────────────────────────────────────────
fa <- mean(res$singular, na.rm = TRUE)
h  <- res[res$is_hit, ]
fh <- if (nrow(h)) mean(h$singular, na.rm = TRUE) else NA
ratio <- if (!is.na(fh) && fa > 0) fh / fa else NA
warn <- reliable && nrow(h) >= 5 && !is.na(ratio) && ratio >= opt$warn_ratio

per_el <- lapply(split(res, res$element), function(e) list(
  element = e$element[1],
  pct_singular_all = round(100 * mean(e$singular, na.rm = TRUE), 1),
  n_hits = sum(e$is_hit),
  n_hits_singular = sum(e$singular & e$is_hit, na.rm = TRUE)))

summ <- list(
  safeguard = "S3_singular_fits", applicable = TRUE,
  random_effects = paste(RANDOM, collapse = " "),
  n_tests = nrow(res), n_singular = sum(res$singular, na.rm = TRUE),
  pct_singular_all = round(100 * fa, 1),
  n_hits = nrow(h), n_hits_singular = sum(h$singular, na.rm = TRUE),
  pct_singular_hits = if (is.na(fh)) NA else round(100 * fh, 1),
  ratio_hits_vs_all = if (is.na(ratio)) NA else round(ratio, 2),
  warn_ratio = opt$warn_ratio, warning = warn,
  refit_spearman_rho = if (is.na(rho)) NA else round(rho, 4), refit_reliable = reliable,
  per_element = unname(per_el))

cat(sprintf("S3 singular fits: %.1f%% of all tests, %s of hits (%d/%d)\n", 100 * fa,
            if (is.na(fh)) "n/a" else sprintf("%.1f%%", 100 * fh),
            sum(h$singular, na.rm = TRUE), nrow(h)))
if (warn)
  cat(sprintf(paste0("SAFEGUARD_WARNING: S3 singular fits are %.1fx more common among hits than ",
                     "among all tests: the %s random effect is often estimated at zero, so ",
                     "these hits may be group differences rather than covariate effects\n"),
              ratio, paste(RANDOM, collapse = "/")))

if (opt$within_between) {
  nw <- sum(h$within_q < 0.05, na.rm = TRUE)
  nb <- sum(h$between_p < 0.05, na.rm = TRUE)
  summ$within_between <- list(grouping = RANDOM[1], n_hits = nrow(h),
                              n_hits_between_p05 = nb, n_hits_within_q05 = nw,
                              n_all_within_q05 = sum(res$within_q < 0.05, na.rm = TRUE))
  cat(sprintf("  within/between (%s): %d of %d hits have a between-group term at p < 0.05; %d have a within-group term at q < 0.05\n",
              RANDOM[1], nb, nrow(h), nw))
}
write_json(summ, opt$out_json)
