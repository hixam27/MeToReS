#!/usr/bin/env python3
"""paper_numbers.py — print every number the manuscript needs from one MeToReS run.

Usage (from the MeToReS folder, metores environment active):
    python validation/paper_numbers.py results_forhot group > paper_numbers.txt
    python validation/paper_numbers.py results_forhot group ".snakemake/log/2026-11-*.snakemake.log"

The optional third argument (a glob of Snakemake log files) adds wall-clock
time per rule (Table S6). MeToReS rules have no benchmark directives, so peak
memory per rule is not recorded.

Read-only: it opens result files and prints values; it changes nothing.
"""
import glob
import json
import os
import sys

import pandas as pd

R = sys.argv[1] if len(sys.argv) > 1 else "results"
GC = sys.argv[2] if len(sys.argv) > 2 else "group"
LOGS = sys.argv[3] if len(sys.argv) > 3 else None


def section(t):
    print(f"\n=== {t} ===")


def first(pattern):
    hits = sorted(glob.glob(os.path.join(R, pattern)))
    return hits[0] if hits else None


def count_fasta(p):
    n = 0
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


# 1. Catalogue and assembly QC
section("1. Gene catalogue (assembly QC)")
f = first("qc/assembly_qc/assembly_qc_summary.json")
if f:
    d = json.load(open(f))
    def walk(x, pre=""):
        if isinstance(x, dict):
            for k, v in x.items():
                walk(v, f"{pre}{k}.")
        elif not isinstance(x, list):
            if any(s in pre.lower() for s in ("total_genes", "n50", "total_bp", "gc", "mean", "median")):
                print(f"  {pre[:-1]} = {x}")
    walk(d)
cat = first("gene_catalog/all_longest_orfs_cds_rmdup_id.fasta")
merged = first("gene_catalog/all_orfs_merged.fasta")
if cat:
    print(f"  catalogue genes (count) = {count_fasta(cat):,}")
if merged:
    print(f"  ORFs before clustering (>= min length) = {count_fasta(merged):,}")

# 2. Salmon mapping and library types
section("2. Salmon mapping to the catalogue")
rows = []
for m in sorted(glob.glob(os.path.join(R, "quant", "*_quant", "aux_info", "meta_info.json"))):
    s = os.path.basename(os.path.dirname(os.path.dirname(m))).replace("_quant", "")
    j = json.load(open(m))
    rows.append((s, j.get("percent_mapped"), j.get("num_processed"), j.get("num_mapped"),
                 (j.get("library_types") or [None])[0]))
if rows:
    df = pd.DataFrame(rows, columns=["sample", "pct_mapped", "reads_processed", "reads_mapped", "libtype"])
    print(df.to_string(index=False))
    print(f"  mean % mapped = {df.pct_mapped.mean():.1f}  (min {df.pct_mapped.min():.1f}, max {df.pct_mapped.max():.1f})")
    print(f"  total read pairs processed = {df.reads_processed.sum():,}")
for p in sorted(glob.glob(os.path.join(R, "safeguards", "libtype_probe", "*.txt"))):
    pass
probes = sorted(glob.glob(os.path.join(R, "safeguards", "libtype_probe", "*.txt")))
if probes:
    det = pd.Series([open(p).read().strip() for p in probes]).value_counts()
    print("  auto-detected library types: " + ", ".join(f"{k} {v}" for k, v in det.items()))

# 3. Gene filtering and DESeq2
section("3. Gene filtering and differential expression")
fm = first(f"deg/{GC}_deg_results/transcript_abundance_quantification_table_filter.csv")
if fm:
    n = sum(1 for _ in open(fm)) - 1
    print(f"  genes after the expression filter (filtered matrix) = {n:,}")
gs = first(f"deg/{GC}_deg_results/gene_set_summary.json")
if gs:
    g = json.load(open(gs))
    print(f"  gene set (deseq2.gene_set) = {g['gene_set']}; genes passed to DESeq2 = {g['n_genes_selected']:,}")
    print(f"  with protein-level evidence (eggNOG-mapper hit) = {g['n_with_protein_evidence']:,} "
          f"({g['pct_genes_with_protein_evidence']}% of filtered genes, "
          f"{g['pct_counts_with_protein_evidence']}% of their counts; per sample "
          f"{g['pct_counts_min']}-{g['pct_counts_max']}%)")
s = first(f"deg/{GC}_deg_results/all_contrasts_summary.csv")
if s:
    t = pd.read_csv(s)
    print(t.to_string(index=False))
    if {"Upregulated", "Downregulated"} <= set(t.columns):
        up, dn = t["Upregulated"].sum(), t["Downregulated"].sum()
        print(f"  DEGs summed over contrasts = {up + dn:,}  ({100 * dn / (up + dn):.1f}% down)")
tested = set()
for f2 in glob.glob(os.path.join(R, "deg", f"{GC}_deg_results", "*", "deseq2_results_all.csv")):
    tested |= set(pd.read_csv(f2, usecols=["GeneID"])["GeneID"].astype(str))
if tested:
    print(f"  genes tested (union over contrasts) = {len(tested):,}")
ids = set()
for f2 in glob.glob(os.path.join(R, "deg", f"{GC}_deg_results", "*", "differential_genes_id.txt")):
    ids |= {l.strip() for l in open(f2) if l.strip()}
if ids:
    print(f"  unique DEGs (union over contrasts) = {len(ids):,}")

# 4. Taxonomy
section("4. Kaiju taxonomy")
k = first("kaiju/kaiju_summary.json")
if k:
    kj = json.load(open(k))
    for key in ("classification_rate", "total_orfs", "classified", "unclassified"):
        if key in kj:
            print(f"  {key} = {kj[key]}")
ks = first("kaiju/kaiju_summary_phylum.tsv")
if ks:
    print("  top phyla (catalogue ORFs):")
    print(pd.read_csv(ks, sep="\t").head(8).to_string(index=False))

# 5. Functional annotation
section("5. eggNOG-mapper")
e = first(f"emapper_{GC}/all_genes_emapper.emapper.annotations")
if e:
    lines = [l for l in open(e) if not l.startswith("##")]
    hdr = lines[0].lstrip("#").rstrip("\n").split("\t")
    ann = pd.DataFrame([l.rstrip("\n").split("\t") for l in lines[1:]], columns=hdr)
    print(f"  genes with an eggNOG hit = {len(ann):,}")
    if fm:
        print(f"  = {100 * len(ann) / n:.1f}% of the {n:,} filtered genes")
    for col in ("COG_category", "GOs", "KEGG_ko"):
        if col in ann.columns:
            k_ = (ann[col].fillna("-") != "-").sum()
            print(f"  with {col}: {k_:,} ({100 * k_ / max(1, n if fm else len(ann)):.1f}% of filtered genes)")
    if tested and "KEGG_ko" in ann.columns:
        q = ann.columns[0]
        with_ko = set(ann.loc[ann["KEGG_ko"].fillna("-") != "-", q].astype(str))
        k_ = len(tested & with_ko)
        print(f"  tested genes with a KEGG orthologue = {k_:,} of {len(tested):,} ({100 * k_ / len(tested):.1f}%)")

# 6. Community analysis
section("6. ANOSIM and dbRDA")
a = first(f"downstream_{GC}/anosim_result.csv")
if a:
    print(pd.read_csv(a).to_string(index=False))
db = first(f"constrained_ordination_{GC}/dbrda_anova_results.csv")
if db:
    print(pd.read_csv(db).to_string(index=False))
vp = first(f"constrained_ordination_{GC}/variance_partitioning_fractions.csv")
if vp:
    print(pd.read_csv(vp).to_string(index=False))

# 7. Covariate integration
section("7. envfit")
ef = first(f"chemical_{GC}/envfit_results.csv")
if ef:
    t = pd.read_csv(ef)
    print(t.to_string(index=False))
    print(f"  variables tested = {len(t)}; p < 0.05: {(t.pvalue < 0.05).sum()}"
          + (f"; padj < 0.05: {(t.padj < 0.05).sum()}" if "padj" in t.columns else ""))
section("8. Gene-level associations")
ml = first(f"maaslin2_{GC}/maaslin2_significant.tsv")
if ml:
    t = pd.read_csv(ml, sep="\t")
    print(f"  MaAsLin2 significant associations = {len(t)}")
mla = first(f"maaslin2_{GC}/maaslin2_all_elements.tsv")
if mla:
    t = pd.read_csv(mla, sep=chr(9))
    print(f"  MaAsLin2 tests = {len(t):,}")
    if {"element", "qval", "qval_across_elements"} <= set(t.columns):
        print(f"  lowest q across variables = {t['qval_across_elements'].min():.3g}")
        for e, g in t.groupby("element"):
            print(f"    {e}: {g['feature'].nunique():,} genes tested; lowest q within variable "
                  f"{g['qval'].min():.3g}; q < 0.05 within variable: {(g['qval'] < 0.05).sum()}")
pp = first(f"bioindicators_{GC}/partial_correlations_significant.csv")
if pp:
    print(f"  partial correlation significant pairs = {len(pd.read_csv(pp))}")
ppa = first(f"bioindicators_{GC}/partial_correlations_all.csv")
if ppa:
    t = pd.read_csv(ppa)
    print(f"  partial correlation tests = {len(t):,}"
          + (f"; lowest adjusted P = {t['padj'].min():.3g}" if "padj" in t.columns else ""))
rc = first(f"chemical_{GC}/element_compound_correlations.csv")
if rc:
    t = pd.read_csv(rc)
    print(f"  resistance–variable pairs tested = {len(t)}"
          + (f"; padj < 0.05: {(t.padj < 0.05).sum()}" if "padj" in t.columns else ""))
fc = first(f"chemical_{GC}/full_correlation_matrix.csv")
if fc:
    t = pd.read_csv(fc)
    print(f"  full resistance x variable matrix rows = {len(t):,}"
          + (f"; padj < 0.05: {(t.padj < 0.05).sum()}" if "padj" in t.columns else ""))

# 9. Resistome
section("9. Resistome")
ra = first(f"resistance_{GC}/resistance_annotation.tsv")
if ra:
    t = pd.read_csv(ra, sep="\t", low_memory=False)
    for col in ("is_ARG", "is_MRG"):
        if col in t.columns:
            print(f"  {col} genes = {t[col].astype(str).str.lower().isin(['true', '1']).sum():,}")
aro = first(f"resistance_{GC}/arg_abundance_by_aro.csv")
if aro:
    print(f"  ARO terms with expression = {len(pd.read_csv(aro))}")
fam = first(f"resistance_{GC}/arg_abundance_by_family.csv")
if fam:
    print(f"  AMR gene families with expression = {len(pd.read_csv(fam))}")
cov = first(f"resistance_{GC}/resistome_coverage.csv")
if cov:
    print("  resistome coverage report:")
    print(pd.read_csv(cov).head(15).to_string(index=False))

# 10. Safeguards
section("10. Safeguards: key numbers")
pg = first("safeguards/prevalent_genes.txt")
if pg:
    n_pg = sum(1 for l in open(pg) if l.strip())
    n_sel = json.load(open(gs)).get("n_genes_selected") if gs else None
    print(f"  S1 genes passing the prevalence rule = {n_pg:,}"
          + (f" of the {n_sel:,} genes passed to DESeq2" if n_sel else ""))
lc = first("safeguards/library_complexity.json")
if lc:
    j = json.load(open(lc))
    print(f"  S2 low-complexity flags = {j.get('n_flagged')}; library-type mismatches = "
          f"{len(j.get('libtype_mismatch', []))}; contrasts mixing types = "
          f"{len(j.get('contrasts_mixing_libtypes', []))}")
sf = first("safeguards/maaslin2_singular_fits.json")
if sf:
    j = json.load(open(sf))
    print(f"  S3 singular fits = {j.get('pct_singular_all')}% of {j.get('n_tests'):,} tests; hits = {j.get('n_hits')}")
ao = first("safeguards/assembly_origin_enrichment.json")
if ao:
    j = json.load(open(ao))
    print(f"  S4 hit sets flagged = {j.get('n_sets_flagged')} of {j.get('n_sets_tested')}")
    for f_ in j.get("flagged", []):
        print("    " + ", ".join(f"{k}={v}" for k, v in f_.items() if not isinstance(v, (list, dict))))
section("10b. Safeguards summary (first 60 lines)")
sg = first("safeguards/SAFEGUARDS_SUMMARY.md")
if sg:
    print("".join(open(sg).readlines()[:60]))

# 11. Run time per rule (Table S6), from Snakemake logs
if LOGS:
    import re
    from datetime import datetime
    section("11. Wall-clock time per rule (Snakemake logs)")
    ts_re = re.compile(r"^\[(\w{3} \w{3} +\d+ \d\d:\d\d:\d\d \d{4})\]")
    jobs, rows = {}, []
    for lf in sorted(glob.glob(LOGS)):
        now, rule = None, None
        for line in open(lf, errors="replace"):
            m = ts_re.match(line)
            if m:
                now = datetime.strptime(" ".join(m.group(1).split()), "%a %b %d %H:%M:%S %Y")
                continue
            s = line.strip()
            m = re.match(r"^(?:local)?rule (\S+):$", s)
            if m:
                rule = m.group(1)
                continue
            m = re.match(r"^jobid: (\d+)$", s)
            if m and rule and now:
                jobs[(lf, m.group(1))] = (rule, now)
                rule = None
                continue
            m = re.match(r"^Finished job(?:id)?:? (\d+)", s)
            if m and now and (lf, m.group(1)) in jobs:
                r_, start = jobs.pop((lf, m.group(1)))
                rows.append((r_, start, now))
    if rows:
        d = pd.DataFrame(rows, columns=["rule", "start", "end"])
        d["hours"] = (d["end"] - d["start"]).dt.total_seconds() / 3600
        out = d.groupby("rule").agg(jobs=("hours", "size"), job_hours_sum=("hours", "sum"),
                                    job_hours_max=("hours", "max"), first_start=("start", "min"),
                                    last_end=("end", "max")).sort_values("first_start")
        out["span_hours"] = (out["last_end"] - out["first_start"]).dt.total_seconds() / 3600
        num = ["job_hours_sum", "job_hours_max", "span_hours"]
        out[num] = out[num].round(2)
        print(out.to_string())
        print(f"  whole run: {d['start'].min()} to {d['end'].max()} "
              f"({(d['end'].max() - d['start'].min()).total_seconds() / 3600:.1f} h, across "
              f"{len(glob.glob(LOGS))} log file(s))")
    else:
        print("  no finished jobs found in the log files")
