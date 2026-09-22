#!/usr/bin/env python3
"""generate_pipeline_summary.py — consolidates key pipeline metrics
(assembly, quantification, differential expression, taxonomy, and
environmental correlations) into one human-readable summary file. All
values are read from existing pipeline outputs.
"""
import argparse
import os
import sys
import json
import glob
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--results_dir", required=True, help="Pipeline's top-level results/ directory")
ap.add_argument("--group_col", required=True, help="GROUP_COL used for folder naming (e.g. 'country')")
ap.add_argument("--out_file", required=True)
a = ap.parse_args()

R = a.results_dir
GC = a.group_col


def safe_read(path, kind="csv", **kwargs):
    if not os.path.exists(path):
        return None
    try:
        if kind == "csv":
            return pd.read_csv(path, **kwargs)
        if kind == "tsv":
            return pd.read_csv(path, sep="\t", **kwargs)
        if kind == "json":
            with open(path) as fh:
                return json.load(fh)
    except Exception as e:
        print(f"WARNING: could not read {path}: {e}", file=sys.stderr)
        return None


lines = []
def add(s=""):
    lines.append(s)


add("# MeToReS Pipeline Summary")
add("")
add("Key metrics from this pipeline run. Each value links to its source file below.")
add("")
add("## Key Metrics")
add("")

# =============================================================================
# 1. Assembly / gene catalogue
# =============================================================================
catalog_fasta = os.path.join(R, "gene_catalog", "all_longest_orfs_cds_rmdup_id.fasta")
n_genes = None
if os.path.exists(catalog_fasta):
    with open(catalog_fasta) as fh:
        n_genes = sum(1 for line in fh if line.startswith(">"))

n50 = None
n50_source = None
n50_warning = None

# seqkit_stats.tsv is computed directly on the current catalog fasta, so
# its sequence count is guaranteed to match n_genes above. MetaQUAST can
# apply its own internal filtering, so its own "# contigs" total may cover
# only a subset of the real catalog — checked explicitly below rather than
# trusted silently, since that subset's N50 would then be a different,
seqkit_stats = os.path.join(R, "qc", "assembly_qc", "seqkit_stats.tsv")
sk = safe_read(seqkit_stats, kind="tsv")
if sk is not None:
    n50_cols = [c for c in sk.columns if str(c).strip().upper() == "N50"]
    if n50_cols:
        n50 = sk[n50_cols[0]].iloc[0]
        n50_source = seqkit_stats

metaquast_report = glob.glob(os.path.join(R, "qc", "assembly_qc", "metaquast", "**", "report.tsv"),
                              recursive=True)
if metaquast_report:
    mq = safe_read(metaquast_report[0], kind="tsv", header=None, index_col=0)
    if mq is not None:
        contig_rows = [i for i in mq.index if isinstance(i, str) and i.strip() == "# contigs"]
        if contig_rows and n_genes is not None:
            mq_n_contigs = int(mq.loc[contig_rows[0]].iloc[0])
            if mq_n_contigs != n_genes:
                n50_warning = (f"MetaQUAST's own sequence count ({mq_n_contigs:,}) does not "
                               f"match the real catalogue size ({n_genes:,}) — MetaQUAST is "
                               f"reporting on a subset, likely due to internal filtering. Its "
                               f"N50 is NOT used below for this reason.")
        if n50 is None:
            n50_rows = [i for i in mq.index if isinstance(i, str) and i.strip().upper() == "N50"]
            if n50_rows:
                n50 = mq.loc[n50_rows[0]].iloc[0]
                n50_source = metaquast_report[0]

add(f"Gene catalogue size (N)      : {n_genes:,}" if n_genes is not None else
    "Gene catalogue size (N)      : NOT FOUND")
add(f"  Source: {catalog_fasta}")
add(f"Catalogue N50 (X, bp)        : {n50}" if n50 is not None else
    "Catalogue N50 (X, bp)        : NOT FOUND (checked seqkit_stats.tsv and MetaQUAST report.tsv)")
if n50_source:
    add(f"  Source: {n50_source}")
if n50_warning:
    add(f"  NOTE: {n50_warning}")
add("")

# =============================================================================
# 2. Salmon mapping rate (average across samples)
# =============================================================================
meta_info_files = glob.glob(os.path.join(R, "quant", "*_quant", "aux_info", "meta_info.json"))
mapped_pcts = []
for f in meta_info_files:
    d = safe_read(f, kind="json")
    if d and "percent_mapped" in d:
        mapped_pcts.append(d["percent_mapped"])

if mapped_pcts:
    avg_mapped = sum(mapped_pcts) / len(mapped_pcts)
    add(f"Mean reads mapped to catalogue (Y%) : {avg_mapped:.1f}% "
        f"(n={len(mapped_pcts)} samples, range {min(mapped_pcts):.1f}-{max(mapped_pcts):.1f}%)")
else:
    add("Mean reads mapped to catalogue (Y%) : NOT FOUND")
add(f"  Source: results/quant/*_quant/aux_info/meta_info.json ({len(meta_info_files)} files found)")
add("")

# =============================================================================
# 3. DESeq2 — DEGs and contrast count
# =============================================================================
deg_summary_path = None
deg_candidates = glob.glob(os.path.join(R, "deg", "*_deg_results", "all_contrasts_summary.csv"))
if deg_candidates:
    deg_summary_path = deg_candidates[0]
    deg_df = safe_read(deg_summary_path)
else:
    deg_df = None

if deg_df is not None and "Total_DEGs" in deg_df.columns:
    total_degs = int(deg_df["Total_DEGs"].sum())
    n_contrasts = len(deg_df)
    add(f"Total DEGs (N) across contrasts (k) : {total_degs:,} DEGs across {n_contrasts} contrasts")
    add(f"  Source: {deg_summary_path}")
    add("  Per-contrast breakdown:")
    for _, row in deg_df.iterrows():
        add(f"    {row['Contrast']}: {row['Total_DEGs']} DEGs "
            f"({row['Upregulated']} up / {row['Downregulated']} down, "
            f"{row['Genes_tested']} genes tested)")
else:
    add("Total DEGs (N) across contrasts (k) : NOT FOUND")
add("")

# =============================================================================
# 4. Kaiju — classification rate and dominant phyla
# =============================================================================
kaiju_json = safe_read(os.path.join(R, "kaiju", "kaiju_summary.json"), kind="json")
if kaiju_json:
    add(f"ORFs classified (Z%)         : {kaiju_json.get('classification_rate', 'NOT FOUND')}%")
    add(f"  Source: results/kaiju/kaiju_summary.json")
    # Exclude "NA" (ORFs with a domain/kingdom-level hit but no resolved
    # phylum) — not a real phylum, shouldn't occupy a "dominant phyla" slot.
    # Filtered before truncating to top 5, so a genuine phylum fills that
    # slot instead of being silently pushed out.
    all_phyla = kaiju_json.get("top10_phylums", [])
    top_phyla = [p for p in all_phyla if p.get("phylum", "").upper() != "NA"][:5]
    if top_phyla:
        phyla_str = ", ".join(f"{p['phylum']} ({p['percentage']}%)" for p in top_phyla)
        add(f"Dominant phyla                : {phyla_str}")
else:
    add("ORFs classified (Z%)         : NOT FOUND")
add("")

# =============================================================================
# 5. envfit — significant elements
# =============================================================================
envfit_candidates = glob.glob(os.path.join(R, f"chemical_{GC}", "envfit_results.csv"))
if envfit_candidates:
    envfit_df = safe_read(envfit_candidates[0])
    if envfit_df is not None and "pvalue" in envfit_df.columns:
        n_total = len(envfit_df)
        sig = envfit_df[envfit_df["pvalue"] < 0.05]
        n_sig = len(sig)
        sig_elements = sig["Element"].tolist() if "Element" in sig.columns else []
        add(f"Significant elements in envfit (n of {n_total}) : {n_sig}")
        if sig_elements:
            add(f"  Elements: {', '.join(sig_elements)}")
        add(f"  Source: {envfit_candidates[0]}")
    else:
        add("Significant elements in envfit : NOT FOUND (envfit_results.csv missing 'pvalue' column)")
else:
    add("Significant elements in envfit : NOT FOUND (chemical module likely not enabled)")
add("")

# =============================================================================
# 6. Expression-concentration correlations — ALL THREE candidate sources,
# clearly labeled, since which one is the relevant "significant pairs"
# metric depends on the analysis design being used, not fixed in advance.
# =============================================================================
add("Significant expression-concentration relationships (three candidate")
add("sources reported below — pick the one matching your intended")
add("statistical design before citing a single number):")
add("")

ppcor_candidates = glob.glob(os.path.join(R, f"bioindicators_{GC}", "partial_correlations_significant.csv"))
if ppcor_candidates:
    ppcor_df = safe_read(ppcor_candidates[0])
    n = len(ppcor_df) if ppcor_df is not None else 0
    add(f"  (a) Partial correlation (ppcor), gene-element pairs, padj<0.05 : {n}")
    add(f"      Source: {ppcor_candidates[0]}")
else:
    add("  (a) Partial correlation (ppcor) : NOT FOUND (bioindicators module likely not enabled)")

maaslin_candidates = glob.glob(os.path.join(R, f"maaslin2_{GC}", "maaslin2_significant.tsv"))
if maaslin_candidates:
    maaslin_df = safe_read(maaslin_candidates[0], kind="tsv")
    n = len(maaslin_df) if maaslin_df is not None else 0
    add(f"  (b) MaAsLin2, gene-element associations, qval<0.05           : {n}")
    add(f"      Source: {maaslin_candidates[0]}")
else:
    add("  (b) MaAsLin2 : NOT FOUND (maaslin2 module likely not enabled)")

rescor_candidates = glob.glob(os.path.join(R, f"chemical_{GC}", "element_compound_correlations.csv"))
if rescor_candidates:
    rescor_df = safe_read(rescor_candidates[0])
    if rescor_df is not None and "padj" in rescor_df.columns:
        n_total = len(rescor_df)
        n_sig = int((rescor_df["padj"] < 0.05).sum())
        add(f"  (c) Resistance-chemical correlation, direct element-compound pairs, "
            f"padj<0.05 : {n_sig} (of {n_total} tested)")
    else:
        # Older output without a padj column — report what's available rather
        # than silently treat row count as a significance count.
        n_total = len(rescor_df) if rescor_df is not None else 0
        add(f"  (c) Resistance-chemical correlation : {n_total} pairs computed "
            f"(no padj column found — cannot report significant count)")
    add(f"      Source: {rescor_candidates[0]}")
else:
    add("  (c) Resistance-chemical correlation : NOT FOUND (resistance/chemical module likely not enabled)")
add("")

# =============================================================================
# Write out
# =============================================================================
os.makedirs(os.path.dirname(a.out_file), exist_ok=True)
with open(a.out_file, "w") as fh:
    fh.write("\n".join(lines) + "\n")

print(f"Pipeline summary written: {a.out_file}")
print(f"({len(lines)} lines)")
