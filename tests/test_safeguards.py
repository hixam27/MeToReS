#!/usr/bin/env python3
"""test_safeguards.py — positive and negative controls for safeguards S1–S4.

Each safeguard script is run directly (no Snakemake, no real data) on small
synthetic inputs where the right answer is known by construction:

  S1  hits from genes detected in 2 of 14 samples        -> flagged
      hits from genes detected in every sample           -> not flagged
  S2  a library with a tenth of the genes of the others  -> flagged
      a declared library type that Salmon does not detect -> flagged
      a contrast whose samples mix library types         -> flagged
      the other libraries, types and contrasts           -> not flagged
  S3  hits that only follow differences between sites    -> singular fits, warning,
                                                            between-site part only
      hits that follow the variable within sites         -> no warning, within part
  S4  a hit set planted to come from one sample's assembly -> flagged
      a hit set chosen by expression level, where one deep sample assembles the
      most expressed genes (no origin effect)            -> not flagged by the
                                                            matched test

Usage (from the repository root):
    python tests/test_safeguards.py [--out_md safeguard_controls.md] [--keep DIR]

Exits with status 1 if any check fails. Needs pandas, numpy and Rscript with
optparse, lme4 and lmerTest (all in the MeToReS environment).
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")
RESULTS = []          # (safeguard, control, expected, observed, passed)


def check(sg, control, expected, observed, passed):
    RESULTS.append((sg, control, expected, observed, bool(passed)))


def run(cmd, cwd):
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if p.returncode != 0:
        sys.stderr.write(f"\nCOMMAND FAILED: {' '.join(cmd)}\n{p.stdout}\n{p.stderr}\n")
        raise SystemExit(1)
    return p.stdout


def load(path):
    with open(path) as fh:
        return json.load(fh)


# ── shared design: 7 sites x 2 samples ─────────────────────────────────────
SITES = [f"site{i}" for i in range(1, 8)]
SAMPLES = [f"{s}_{r}" for s in SITES for r in (1, 2)]
SRAS = [f"SRR{100 + i}" for i in range(len(SAMPLES))]


def samplesheet(path, extra=None):
    d = pd.DataFrame({"sample": SAMPLES, "sra": SRAS,
                      "site": [s.split("_")[0] for s in SAMPLES]})
    for k, v in (extra or {}).items():
        d[k] = v
    d.to_csv(path, index=False)
    return d


# ── S1 prevalence ──────────────────────────────────────────────────────────
def test_s1(w, rng):
    ss = os.path.join(w, "s1_samplesheet.csv")
    samplesheet(ss)
    n = 200
    mat = rng.poisson(20, (n, len(SAMPLES))).astype(float) + 1     # detected everywhere
    genes = [f"g{i}" for i in range(n)]
    rare = genes[:20]
    for i in range(20):                                             # detected in 2 samples
        keep = rng.choice(len(SAMPLES), 2, replace=False)
        row = np.zeros(len(SAMPLES)); row[keep] = rng.poisson(50, 2) + 1
        mat[i] = row
    mfile = os.path.join(w, "s1_matrix.csv")
    pd.DataFrame(mat, index=genes, columns=SAMPLES).to_csv(mfile)
    rule = ["--min_fraction", "0.5", "--samplesheet", ss, "--group_col", "site", "--min_in_group", "3"]

    out = run([sys.executable, os.path.join(SCRIPTS, "safeguard_prevalence.py"), "genes",
               "--matrix", mfile, *rule, "--out", os.path.join(w, "s1_prevalent.txt")], w)
    kept = set(open(os.path.join(w, "s1_prevalent.txt")).read().split())
    check("S1", "Filter: genes detected in 2/14 samples", "removed",
          f"{len(set(rare) - kept)}/20 removed", not (set(rare) & kept))
    check("S1", "Filter: genes detected in 14/14 samples", "kept",
          f"{len(set(genes[20:]) & kept)}/{n - 20} kept", set(genes[20:]) <= kept)

    for label, hit_genes, want in (("hits from rare genes", rare[:10], True),
                                   ("hits from common genes", genes[50:60], False)):
        hf = os.path.join(w, f"s1_hits_{want}.tsv")
        pd.DataFrame({"feature": hit_genes, "element": "X"}).to_csv(hf, sep="\t", index=False)
        js = os.path.join(w, f"s1_{want}.json")
        out = run([sys.executable, os.path.join(SCRIPTS, "safeguard_prevalence.py"), "annotate",
                   "--matrix", mfile, "--hits", hf, "--gene_col", "feature", "--label", "test",
                   *rule, "--out_tsv", os.path.join(w, f"s1_{want}.tsv"), "--out_json", js], w)
        j = load(js)
        warned = "SAFEGUARD_WARNING" in out
        if want:
            check("S1", f"Report: {label}", "10/10 flagged, warning",
                  f"{j['n_hits_low_prevalence']}/10 flagged, warning {'yes' if warned else 'no'}",
                  j["n_hits_low_prevalence"] == 10 and warned)
        else:
            check("S1", f"Report: {label}", "0/10 flagged, no warning",
                  f"{j['n_hits_low_prevalence']}/10 flagged, warning {'yes' if warned else 'no'}",
                  j["n_hits_low_prevalence"] == 0 and not warned)


# ── S2 library complexity and type ─────────────────────────────────────────
def test_s2(w, rng):
    low = SAMPLES[3]                                    # site2_2
    groups = [s.split("_")[0] for s in SAMPLES]
    # site3 is declared IU (unstranded) -> contrast site3_vs_site1 mixes types
    libtype = ["IU" if g == "site3" else "" for g in groups]
    ss = os.path.join(w, "s2_samplesheet.csv")
    samplesheet(ss, {"group": groups, "libtype": libtype})
    n = 2000
    mat = rng.poisson(10, (n, len(SAMPLES))).astype(float)
    j_low = SAMPLES.index(low)
    mat[rng.choice(n, int(0.9 * n), replace=False), j_low] = 0      # 10% of genes detected
    mfile = os.path.join(w, "s2_matrix.csv")
    pd.DataFrame(mat, index=[f"g{i}" for i in range(n)], columns=SAMPLES).to_csv(mfile)

    asm, quant, probe = (os.path.join(w, d) for d in ("s2_assembly", "s2_quant", "s2_probe"))
    wrong = SRAS[SAMPLES.index("site5_1")]               # declared ISR, reads look IU
    for s, sra, lt in zip(SAMPLES, SRAS, libtype):
        os.makedirs(os.path.join(asm, sra)); os.makedirs(os.path.join(quant, f"{sra}_quant", "aux_info"))
        os.makedirs(probe, exist_ok=True)
        with open(os.path.join(asm, sra, "orfs.fasta"), "w") as fh:
            fh.write(">o\nATG\n" * int(rng.integers(900, 1100)))
        with open(os.path.join(quant, f"{sra}_quant", "aux_info", "meta_info.json"), "w") as fh:
            json.dump({"num_mapped": int(rng.integers(9e6, 11e6))}, fh)
        with open(os.path.join(probe, f"{sra}.txt"), "w") as fh:
            fh.write("IU\n" if (lt == "IU" or sra == wrong) else "ISR\n")

    js = os.path.join(w, "s2.json")
    out = run([sys.executable, os.path.join(SCRIPTS, "safeguard_library_complexity.py"),
               "--samplesheet", ss, "--matrix", mfile, "--assembly_dir", asm, "--quant_dir", quant,
               "--warn_fraction", "0.25", "--default_libtype", "ISR", "--probe_dir", probe,
               "--deseq_group_col", "group", "--contrasts", "site2_vs_site1 site3_vs_site1",
               "--out_tsv", os.path.join(w, "s2.tsv"), "--out_json", js,
               "--out_mqc", os.path.join(w, "s2_mqc.yaml")], w)
    j = load(js)
    flagged = [f["sample"] for f in j["flagged"]]
    check("S2", "Library with 10% of the genes detected", "flagged",
          "flagged" if low in flagged else "not flagged", low in flagged)
    check("S2", "13 normal libraries", "not flagged", f"{len(set(flagged) - {low})} flagged",
          set(flagged) == {low})
    mism = [m["sample"] for m in j["libtype_mismatch"]]
    check("S2", "Declared ISR, reads unstranded", "type mismatch",
          "mismatch" if "site5_1" in mism else "no mismatch", "site5_1" in mism)
    check("S2", "13 libraries with matching types (incl. 2 declared IU)", "no mismatch",
          f"{len(set(mism) - {'site5_1'})} mismatches", set(mism) == {"site5_1"})
    mixed = [m["contrast"] for m in j["contrasts_mixing_libtypes"]]
    check("S2", "Contrast mixing ISR and IU", "flagged",
          "flagged" if "site3_vs_site1" in mixed else "not flagged", "site3_vs_site1" in mixed)
    check("S2", "Contrast with one library type", "not flagged",
          "flagged" if "site2_vs_site1" in mixed else "not flagged", "site2_vs_site1" not in mixed)


# ── S3 singular fits ───────────────────────────────────────────────────────
MIMIC_R = r"""
suppressPackageStartupMessages({library(lme4); library(lmerTest)})
a <- commandArgs(TRUE)
vst <- read.csv(a[1], check.names = FALSE); rownames(vst) <- vst[[1]]; vst[[1]] <- NULL
chem <- read.csv(a[2]); ss <- read.csv(a[3])
x <- chem$X[match(colnames(vst), chem$SampleID)]; x <- (x - mean(x)) / sd(x)
site <- ss$site[match(colnames(vst), ss$sample)]
r <- t(sapply(rownames(vst), function(g) {
  f <- suppressMessages(suppressWarnings(lmer(as.numeric(vst[g, ]) ~ x + (1 | site))))
  cf <- summary(f)$coefficients; c(cf["x", "Estimate"], cf["x", "Pr(>|t|)"]) }))
out <- data.frame(feature = rownames(vst), metadata = "X", value = "X", coef = r[, 1],
                  pval = r[, 2], element = "X")
out$qval_across_elements <- p.adjust(out$pval, "BH")
write.table(out, a[4], sep = "\t", quote = FALSE, row.names = FALSE)
"""


def test_s3(w, rng):
    ss = os.path.join(w, "s3_samplesheet.csv")
    d = samplesheet(ss)
    site_idx = np.array([SITES.index(s) for s in d["site"]])
    site_mean = rng.normal(0, 1, len(SITES))
    # pos: the variable barely varies within sites (as for soil chemistry of
    # replicate cores), so a gene that follows site differences also follows
    # the variable; neg: the variable varies within sites as well
    within_sd = {"pos": 0.1, "neg": 0.5}
    within = {k: rng.normal(0, v, len(SAMPLES)) for k, v in within_sd.items()}
    for k in within:
        pd.DataFrame({"SampleID": SAMPLES, "X": 10 + site_mean[site_idx] + within[k]}).to_csv(
            os.path.join(w, f"s3_{k}_chem.csv"), index=False)
    mimic = os.path.join(w, "s3_mimic.R")
    open(mimic, "w").write(MIMIC_R)

    def genes(kind, k, scen):
        rows = []
        for _ in range(k):
            u = rng.normal(0, 1.5, len(SITES))[site_idx]       # site effect unrelated to X
            e = rng.normal(0, 0.3, len(SAMPLES))
            if kind == "between":                              # follows site means of X only
                rows.append(3 * site_mean[site_idx] + rng.normal(0, 0.3, len(SAMPLES)))
            elif kind == "within":                             # follows X within sites
                rows.append(3 * within[scen] / within_sd[scen] * 0.5 + u + e)
            else:
                rows.append(u + e)
        return rows

    outcomes = {}
    for scen, hit_kind in (("pos", "between"), ("neg", "within")):
        rows = genes(hit_kind, 10, scen) + genes("null", 200, scen)
        chem = os.path.join(w, f"s3_{scen}_chem.csv")
        ids = [f"{hit_kind}_{i}" for i in range(10)] + [f"null_{i}" for i in range(200)]
        vst = os.path.join(w, f"s3_{scen}_vst.csv")
        pd.DataFrame(rows, index=ids, columns=SAMPLES).to_csv(vst)
        ml = os.path.join(w, f"s3_{scen}_maaslin.tsv")
        run(["Rscript", mimic, vst, chem, ss, ml], w)
        js = os.path.join(w, f"s3_{scen}.json")
        out = run(["Rscript", os.path.join(SCRIPTS, "safeguard_singular_fits.R"),
                   "--maaslin_all", ml, "--vst", vst, "--chemical", chem,
                   "--samplesheet", ss, "--random_effects", "site", "--transform", "none",
                   "--within_between", "--out_tsv", os.path.join(w, f"s3_{scen}.tsv"),
                   "--out_json", js], w)
        outcomes[scen] = (load(js), "SAFEGUARD_WARNING" in out)

    j, warned = outcomes["pos"]
    wb = j["within_between"]
    check("S3", "Hits that follow site means only", "singular fits enriched among hits, warning",
          f"singular: {j['pct_singular_hits']}% of hits vs {j['pct_singular_all']}% of tests, "
          f"warning {'yes' if warned else 'no'}",
          j["n_hits"] >= 5 and warned)
    check("S3", "Same hits, between/within refit", "between part only",
          f"between P < 0.05: {wb['n_hits_between_p05']}/{wb['n_hits']}; "
          f"within q < 0.05: {wb['n_hits_within_q05']}",
          wb["n_hits"] >= 5 and wb["n_hits_between_p05"] >= 0.8 * wb["n_hits"]
          and wb["n_hits_within_q05"] == 0)
    j, warned = outcomes["neg"]
    wb = j["within_between"]
    check("S3", "Hits that follow the variable within sites", "no warning",
          f"singular: {j['pct_singular_hits']}% of hits vs {j['pct_singular_all']}% of tests, "
          f"warning {'yes' if warned else 'no'}",
          j["n_hits"] >= 5 and not warned)
    check("S3", "Same hits, between/within refit", "within part significant",
          f"within q < 0.05: {wb['n_hits_within_q05']}/{wb['n_hits']}",
          wb["n_hits"] >= 5 and wb["n_hits_within_q05"] >= 0.8 * wb["n_hits"])
    check("S3", "Refit reproduces the association P-values", "Spearman rho > 0.99",
          f"rho = {j['refit_spearman_rho']}", j["refit_reliable"])


# ── S4 assembly origin ─────────────────────────────────────────────────────
def test_s4(w, rng):
    # 4 groups x 2 samples; contrast B_vs_A
    samples = [f"{g}_{r}" for g in "ABCD" for r in (1, 2)]
    sras = [f"ERR{200 + i}" for i in range(len(samples))]
    group = [s[0] for s in samples]
    ss = os.path.join(w, "s4_samplesheet.csv")
    pd.DataFrame({"sample": samples, "sra": sras, "group": group}).to_csv(ss, index=False)
    deep, planted = "C_1", "A_2"                    # deep sample outside the contrast
    n = 6000
    origin = rng.choice(samples, n)
    base = rng.normal(2.0, 1.0, n) + np.where(origin == deep, 1.5, 0.0)   # deep sample assembles
    ids = [f"{sras[samples.index(o)]}_k141_{i}" for i, o in enumerate(origin)]  # expressed genes
    counts = rng.poisson(np.exp(base)[:, None] * np.ones((1, len(samples))))
    mfile = os.path.join(w, "s4_matrix.csv")
    pd.DataFrame(counts, index=ids, columns=samples).to_csv(mfile)

    ids_a = np.array(ids)
    mean_expr = counts.mean(axis=1)
    # planted: 40 random genes from the planted origin + 40 random genes
    from_p = rng.choice(np.where(origin == planted)[0], 40, replace=False)
    rest = rng.choice(np.where(origin != planted)[0], 40, replace=False)
    up = ids_a[np.concatenate([from_p, rest])]
    # expression trap: 80 genes drawn from the 600 most expressed, no origin effect
    top = np.argsort(-mean_expr)[:600]
    down = ids_a[rng.choice(top, 80, replace=False)]

    cdir = os.path.join(w, "s4_deg", "B_vs_A")
    os.makedirs(cdir)
    pd.DataFrame({"GeneID": ids}).to_csv(os.path.join(cdir, "deseq2_results_all.csv"), index=False)
    pd.DataFrame({"GeneID": list(up) + list(down),
                  "change": ["Up"] * len(up) + ["Down"] * len(down)}).to_csv(
        os.path.join(cdir, "differential_genes.csv"), index=False)
    js, tsv = os.path.join(w, "s4.json"), os.path.join(w, "s4.tsv")
    run([sys.executable, os.path.join(SCRIPTS, "safeguard_assembly_origin.py"),
         "--samplesheet", ss, "--deg_dir", os.path.join(w, "s4_deg"), "--contrasts", "B_vs_A",
         "--deg_group_col", "group", "--matrix", mfile, "--out_tsv", tsv, "--out_json", js], w)
    j, t = load(js), pd.read_csv(tsv, sep="\t")
    flagged = {(f["set"], f["origin"]) for f in j["flagged"] if f["level"] == "sample"}
    check("S4", f"Hit set with half of its genes from one sample ({planted})", "flagged",
          "flagged" if ("B_vs_A:up", planted) in flagged else "not flagged",
          ("B_vs_A:up", planted) in flagged)
    down_flags = {o for s, o in flagged if s == "B_vs_A:down"}
    r = t[(t["set"] == "B_vs_A:down") & (t["level"] == "sample") & (t["origin"] == deep)].iloc[0]
    check("S4", f"Hit set chosen by expression; deep sample {deep}", "not flagged (matched test)",
          f"{'flagged' if down_flags else 'not flagged'}; matched OR {r['odds_ratio']:.1f} "
          f"vs unmatched OR {r['odds_ratio_unmatched']:.1f}",
          not down_flags)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out_md", default=None, help="write the results table as Markdown")
    ap.add_argument("--keep", default=None, help="keep the synthetic inputs and outputs here")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    w = os.path.abspath(a.keep) if a.keep else tempfile.mkdtemp(prefix="metores_sg_")
    os.makedirs(w, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    try:
        test_s1(w, rng)
        test_s2(w, rng)
        test_s3(w, rng)
        test_s4(w, rng)
        jsons = [os.path.join(w, f) for f in ("s1_True.json", "s2.json", "s3_pos.json", "s4.json")]
        run([sys.executable, os.path.join(SCRIPTS, "safeguards_summary.py"), "--jsons", *jsons,
             "--out", os.path.join(w, "SAFEGUARDS_SUMMARY.md")], w)
        check("All", "Summary report from the four JSON files", "written",
              "written" if os.path.getsize(os.path.join(w, "SAFEGUARDS_SUMMARY.md")) else "empty",
              os.path.getsize(os.path.join(w, "SAFEGUARDS_SUMMARY.md")) > 0)
    finally:
        if not a.keep:
            shutil.rmtree(w, ignore_errors=True)

    L = ["| Safeguard | Control | Expected | Observed | Pass |", "|---|---|---|---|---|"]
    L += [f"| {s} | {c} | {e} | {o} | {'yes' if p else '**NO**'} |" for s, c, e, o, p in RESULTS]
    print("\n".join(L))
    n_fail = sum(1 for r in RESULTS if not r[4])
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
    if a.out_md:
        with open(a.out_md, "w") as fh:
            fh.write("\n".join(L) + "\n")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
