#!/usr/bin/env python3
"""test_poland_origin_by_element.py — follow-up to the random-effect-
absorption test: the 91 MaAsLin2 significant genes are 8.6x
over-represented for Poland origin (46% vs 5.37% baseline). This checks
whether that skew concentrates in elements where Poland has an unusual
(extreme-ranked) concentration, the same style of check already run for
Greece_Natural, or whether it's uniform across elements regardless of
Poland's concentration rank for each one.
"""
import argparse
import pandas as pd
from scipy.stats import fisher_exact

ap = argparse.ArgumentParser()
ap.add_argument("--catalog", required=True)
ap.add_argument("--samplesheet", required=True)
ap.add_argument("--maaslin2_sig", required=True)
ap.add_argument("--chemical", required=True)
a = ap.parse_args()

# ── gene -> origin deseq2_group mapping ──────────────────────────────
ss = pd.read_csv(a.samplesheet)
sra_to_sample = dict(zip(ss["sra"], ss["sample"]))
sra_ids_sorted = sorted(sra_to_sample.keys(), key=len, reverse=True)
sample_to_group = dict(zip(ss["sample"], ss["deseq2_group"]))

gene_origin_group = {}
with open(a.catalog) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        gene_id = line[1:].strip().split()[0]
        for sra in sra_ids_sorted:
            if gene_id.startswith(sra):
                gene_origin_group[gene_id] = sample_to_group[sra_to_sample[sra]]
                break

baseline_poland_rate = sum(v == "Poland" for v in gene_origin_group.values()) / len(gene_origin_group)
print(f"Catalogue-wide baseline: {baseline_poland_rate*100:.2f}% Poland-origin\n")

# ── Poland's concentration rank per element ──────────────────────────
chem = pd.read_csv(a.chemical)
meta_cols = {"SampleID", "Country", "Season", "Replicate"}
elem_cols = [c for c in chem.columns if c not in meta_cols]

print("=== Poland's concentration rank per element ===")
print(f"{'Element':10} {'Poland mean':>14} {'Rank':>8} {'of':>4} {'Extreme?':>10}")
poland_extreme_elements = set()
poland_rank = {}
for el in elem_cols:
    sub = chem[["SampleID", el]].copy()
    sub[el] = pd.to_numeric(sub[el], errors="coerce")
    sub = sub.merge(ss[["sample", "deseq2_group"]], left_on="SampleID", right_on="sample")
    group_means = sub.groupby("deseq2_group")[el].mean().dropna().sort_values()
    if "Poland" not in group_means.index:
        continue
    rank = list(group_means.index).index("Poland") + 1
    n_groups = len(group_means)
    is_extreme = rank <= 2 or rank >= n_groups - 1
    poland_rank[el] = (rank, n_groups)
    if is_extreme:
        poland_extreme_elements.add(el)
    print(f"{el:10} {group_means['Poland']:14.3f} {rank:8d} {n_groups:4d} {'EXTREME' if is_extreme else '':>10}")

print(f"\nElements where Poland is extreme (rank 1-2 or bottom-2, i.e. highest or lowest): "
      f"{sorted(poland_extreme_elements)}\n")

# ── Per-element breakdown: Poland-origin share of significant hits ──────
sig = pd.read_csv(a.maaslin2_sig, sep="\t")
sig["feature"] = sig["feature"].str.replace(".", "-", regex=False)

print("=== Poland-origin share of significant hits, per element ===")
print(f"{'Element':8} {'n_sig':>6} {'n_Poland':>9} {'%Poland':>9} {'baseline':>9} {'Fisher p':>10} {'Poland extreme here?':>22}")
for el, grp in sig.groupby("element"):
    genes = [g for g in grp["feature"].unique() if g in gene_origin_group]
    if not genes:
        continue
    n_poland = sum(gene_origin_group[g] == "Poland" for g in genes)
    n_total = len(genes)
    pct = 100 * n_poland / n_total

    rest = [g for g in gene_origin_group if g not in genes]
    n_poland_rest = sum(gene_origin_group[g] == "Poland" for g in rest)
    table = [[n_poland, n_total - n_poland], [n_poland_rest, len(rest) - n_poland_rest]]
    _, p = fisher_exact(table, alternative="greater")

    extreme = "yes" if el in poland_extreme_elements else "no"
    flag = "  <-- enriched" if p < 0.05 else ""
    print(f"{el:8} {n_total:6d} {n_poland:9d} {pct:8.2f}% {baseline_poland_rate*100:8.2f}% {p:10.2e} {extreme:>22}{flag}")
