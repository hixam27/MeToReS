#!/usr/bin/env python3
"""test_random_effect_absorption.py — tests Lorenzo's proposed mechanism:
MaAsLin2's (1 | country) random effect absorbs genes detected in
essentially one country, regardless of WHICH country — predicting
depletion for every single-country-origin category, not specifically
Greece_Natural, and no directional mirror (since the mechanism is about
detection breadth, not concentration).

Two parts:
  (1) Origin composition of the significant hits against baseline, across
      all SEVEN deseq2_group origins (matches Lorenzo's specific request;
      note this is one level finer than the model's actual 6-level
      country random effect, since Greece splits into two deseq2_group
      levels but one country).
  (2) Detection breadth: how many of the SIX country groups (matching
      MaAsLin2's actual fitted random effect exactly) each gene is
      detected in (nonzero TPM in >=1 sample from that country),
      comparing significant hits against the full tested background.
"""
import argparse
import pandas as pd
from scipy.stats import chisquare, mannwhitneyu

ap = argparse.ArgumentParser()
ap.add_argument("--catalog", required=True)
ap.add_argument("--samplesheet", required=True)
ap.add_argument("--maaslin2_sig", required=True)
ap.add_argument("--filtered_matrix", required=True,
                 help="transcript_abundance_quantification_table_filter.csv, for detection breadth")
a = ap.parse_args()

# ── Step 1: gene -> origin sample -> origin deseq2_group mapping ────────
ss = pd.read_csv(a.samplesheet)
sra_to_sample = dict(zip(ss["sra"], ss["sample"]))
sra_ids_sorted = sorted(sra_to_sample.keys(), key=len, reverse=True)
sample_to_group = dict(zip(ss["sample"], ss["deseq2_group"]))
sample_to_country = dict(zip(ss["sample"], ss["country"]))

gene_origin_sample = {}
with open(a.catalog) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        gene_id = line[1:].strip().split()[0]
        for sra in sra_ids_sorted:
            if gene_id.startswith(sra):
                gene_origin_sample[gene_id] = sra_to_sample[sra]
                break

gene_origin_group = {g: sample_to_group[s] for g, s in gene_origin_sample.items()}

# Catalog-wide baseline distribution across all 7 deseq2_group origins
baseline_counts = pd.Series(gene_origin_group.values()).value_counts()
baseline_pct = baseline_counts / baseline_counts.sum()
print("=== Part 1: Origin composition, all 7 deseq2_group origins ===\n")
print("Catalogue-wide baseline:")
print((baseline_pct * 100).round(2).to_string())
print()

# ── Significant genes' origin composition ────────────────────────────
sig = pd.read_csv(a.maaslin2_sig, sep="\t")
sig["feature"] = sig["feature"].str.replace(".", "-", regex=False)
sig_genes = sig["feature"].unique()
sig_genes_with_origin = [g for g in sig_genes if g in gene_origin_group]
print(f"Significant genes with recoverable origin: {len(sig_genes_with_origin)} of {len(sig_genes)}\n")

sig_origin_counts = pd.Series([gene_origin_group[g] for g in sig_genes_with_origin]).value_counts()
sig_origin_counts = sig_origin_counts.reindex(baseline_counts.index, fill_value=0)

print(f"{'Origin':22} {'Observed':>10} {'Expected':>10} {'Obs %':>8} {'Baseline %':>12}")
for origin in baseline_counts.index:
    obs = sig_origin_counts[origin]
    exp = baseline_pct[origin] * len(sig_genes_with_origin)
    print(f"{origin:22} {obs:10d} {exp:10.2f} {100*obs/len(sig_genes_with_origin):7.2f}% {100*baseline_pct[origin]:11.2f}%")

expected = (baseline_pct * len(sig_genes_with_origin)).reindex(baseline_counts.index)
chi2, p_chi2 = chisquare(f_obs=sig_origin_counts.values, f_exp=expected.values)
print(f"\nChi-square goodness-of-fit (composition differs from baseline across all 7 origins): "
      f"chi2={chi2:.2f}, p={p_chi2:.2e}")
print("(Lorenzo's mechanism predicts EVERY single-origin category depleted, not just Greece_Natural — "
      "check whether all 7 rows above show obs < expected, not just Greece_Natural.)\n")

# ── Step 2: detection breadth across the 6 country groups ───────────────
print("=== Part 2: Detection breadth (number of countries with nonzero TPM) ===\n")

mat = pd.read_csv(a.filtered_matrix)
mat = mat.rename(columns={mat.columns[0]: "ORF_id"})
sample_cols = [c for c in mat.columns if c in sample_to_country]

country_to_samples = {}
for s in sample_cols:
    c = sample_to_country[s]
    country_to_samples.setdefault(c, []).append(s)

print(f"Countries: {list(country_to_samples.keys())}")

breadth = pd.Series(0, index=mat["ORF_id"])
for country, samples in country_to_samples.items():
    detected_in_country = (mat.set_index("ORF_id")[samples] > 0).any(axis=1)
    breadth += detected_in_country.astype(int)

breadth_sig = breadth.reindex(sig_genes_with_origin).dropna()
breadth_background = breadth  # all tested genes

print(f"\nDetection breadth (countries out of {len(country_to_samples)}), significant hits (n={len(breadth_sig)}):")
print(breadth_sig.describe().round(2).to_string())
print(f"\nDetection breadth, ALL tested genes background (n={len(breadth_background)}):")
print(breadth_background.describe().round(2).to_string())

stat, p_mw = mannwhitneyu(breadth_sig, breadth_background, alternative="greater")
print(f"\nMann-Whitney U (significant hits have BROADER detection than background): "
      f"U={stat:.1f}, p={p_mw:.2e}")
print("(Lorenzo's mechanism predicts significant hits are the broadly-detected subset — "
      "significant p here with hits skewed toward higher breadth supports it.)")
