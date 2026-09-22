#!/usr/bin/env python3
"""test_reference_bias_maaslin2.py — extends the same assembly-origin
enrichment test from test_reference_bias.py to MaAsLin2's significant
gene-element associations, per Lorenzo's §3 follow-up: confirming
`country` is fitted as (1 | country) rules out one failure mode, but
MaAsLin2 runs on the same per-gene, pooled-catalogue matrix as everything
else, so it can carry the same assembly-origin confound.

Genes are deduplicated before testing (a gene's assembly origin is a
property of the gene, not of which element it happened to be tested
against, and the same gene can appear in multiple element-association
rows). Split by the sign of `coef` as a directional check, mirroring the
up/down negative-control structure already used for the DEG test — if
enrichment is present regardless of direction, that's a different,
weaker signal than a direction-specific pattern.
"""
import argparse
import pandas as pd
from scipy.stats import fisher_exact

ap = argparse.ArgumentParser()
ap.add_argument("--catalog", required=True)
ap.add_argument("--samplesheet", required=True)
ap.add_argument("--maaslin2_sig", required=True, help="maaslin2_significant.tsv")
ap.add_argument("--ref_level", default="Greece_Natural")
a = ap.parse_args()

# ── Step 1: map each catalog gene to its origin sample (same logic as
# test_reference_bias.py) ────────────────────────────────────────────────
ss = pd.read_csv(a.samplesheet)
sra_to_sample = dict(zip(ss["sra"], ss["sample"]))
sra_ids_sorted = sorted(sra_to_sample.keys(), key=len, reverse=True)

gene_origin = {}
unmatched = 0
with open(a.catalog) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        gene_id = line[1:].strip().split()[0]
        matched = None
        for sra in sra_ids_sorted:
            if gene_id.startswith(sra):
                matched = sra_to_sample[sra]
                break
        if matched:
            gene_origin[gene_id] = matched
        else:
            unmatched += 1

print(f"Catalog genes with recovered origin sample: {len(gene_origin):,}")
if unmatched:
    print(f"WARNING: {unmatched:,} headers did not match any known SRA prefix.")

ref_samples = set(ss.loc[ss["deseq2_group"] == a.ref_level, "sample"])
origin_is_ref = {g: (s in ref_samples) for g, s in gene_origin.items()}
baseline_rate = sum(origin_is_ref.values()) / len(origin_is_ref)
print(f"Catalogue-wide baseline: {baseline_rate*100:.2f}% of genes "
      f"originally assembled from a {a.ref_level} sample\n")

# ── Step 2: load MaAsLin2 significant results ────────────────────────────
sig = pd.read_csv(a.maaslin2_sig, sep="\t")
required = {"feature", "coef"}
missing = required - set(sig.columns)
if missing:
    raise SystemExit(f"maaslin2_significant.tsv missing expected column(s): {missing}. "
                      f"Available: {list(sig.columns)}")

print(f"Total significant gene-element rows: {len(sig)}")

# MaAsLin2 (via R) sanitizes feature names used as data-frame identifiers,
# converting "-" to "." — e.g. "...CR215-S02-001T0001..." becomes
# "...CR215.S02.001T0001...". Reverse this before matching against the
# catalog's real headers, which still use the original "-".
sig["feature"] = sig["feature"].str.replace(".", "-", regex=False)

print(f"Unique genes (deduplicated for origin testing): {sig['feature'].nunique()}\n")


def test_set(gene_ids, label):
    gene_ids = set(gene_ids) & set(gene_origin.keys())
    if not gene_ids:
        print(f"{label}: no genes with recoverable origin — skipping")
        return
    n_ref = sum(origin_is_ref[g] for g in gene_ids)
    n_total = len(gene_ids)
    pct = 100 * n_ref / n_total

    rest = set(gene_origin.keys()) - gene_ids
    n_ref_rest = sum(origin_is_ref[g] for g in rest)
    table = [[n_ref, n_total - n_ref],
             [n_ref_rest, len(rest) - n_ref_rest]]

    _, p_enriched = fisher_exact(table, alternative="greater")
    _, p_depleted = fisher_exact(table, alternative="less")

    flag = ""
    if p_enriched < 0.05 and pct > baseline_rate * 100:
        flag = "  <-- enriched"
    elif p_depleted < 0.05 and pct < baseline_rate * 100:
        flag = "  <-- depleted"

    print(f"{label:35} n={n_total:5d}  %ref-origin={pct:6.2f}%  "
          f"baseline={baseline_rate*100:6.2f}%  "
          f"p(enriched)={p_enriched:.2e}  p(depleted)={p_depleted:.2e}{flag}")


print("=== All significant associations (deduplicated genes) ===")
test_set(sig["feature"].unique(), "All significant (153 rows)")

print("\n=== Split by direction of association (coef sign) ===")
pos_genes = sig.loc[sig["coef"] > 0, "feature"].unique()
neg_genes = sig.loc[sig["coef"] < 0, "feature"].unique()
test_set(pos_genes, "Positive association (coef>0)")
test_set(neg_genes, "Negative association (coef<0)")
