#!/usr/bin/env python3
"""test_reference_bias_maaslin2_per_element.py — Lorenzo's follow-up:
tests the mirror-enrichment mechanism (positive-association genes
depleted, negative-association genes enriched for Greece_Natural origin)
PER ELEMENT rather than pooled, since the pooled test can hide or dilute
a real per-element pattern, and since his mechanism explicitly assumes
"Greece_Natural is the low-concentration reference site" — which may not
hold identically for every element. That assumption is checked directly
here first, before trusting any per-element enrichment/depletion result
built on top of it.
"""
import argparse
import pandas as pd
from scipy.stats import fisher_exact

ap = argparse.ArgumentParser()
ap.add_argument("--catalog", required=True)
ap.add_argument("--samplesheet", required=True)
ap.add_argument("--maaslin2_sig", required=True)
ap.add_argument("--chemical", required=True, help="chemical_data_<fraction>.csv, to check Greece_Natural's rank per element")
ap.add_argument("--ref_level", default="Greece_Natural")
a = ap.parse_args()

# ── Step 1: gene -> origin sample mapping (same as before) ──────────────
ss = pd.read_csv(a.samplesheet)
sra_to_sample = dict(zip(ss["sra"], ss["sample"]))
sra_ids_sorted = sorted(sra_to_sample.keys(), key=len, reverse=True)

gene_origin = {}
with open(a.catalog) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        gene_id = line[1:].strip().split()[0]
        for sra in sra_ids_sorted:
            if gene_id.startswith(sra):
                gene_origin[gene_id] = sra_to_sample[sra]
                break

ref_samples = set(ss.loc[ss["deseq2_group"] == a.ref_level, "sample"])
origin_is_ref = {g: (s in ref_samples) for g, s in gene_origin.items()}
baseline_rate = sum(origin_is_ref.values()) / len(origin_is_ref)
print(f"Catalogue-wide baseline: {baseline_rate*100:.2f}% {a.ref_level}-origin\n")

# ── Step 2: is Greece_Natural actually the low-concentration group,
# per element? Check directly rather than assume. ───────────────────────
chem = pd.read_csv(a.chemical)
meta_cols = {"SampleID", "Country", "Season", "Replicate"}
elem_cols = [c for c in chem.columns if c not in meta_cols]

ss_ref_samples = ref_samples
chem_ref_mean = {}
chem_rank = {}
print("=== Is Greece_Natural the low-concentration group, per element? ===")
print(f"{'Element':10} {'Greece_Natural mean':>20} {'Rank (1=lowest)':>18} {'Groups compared':>18}")
for el in elem_cols:
    sub = chem[["SampleID", el]].copy()
    sub[el] = pd.to_numeric(sub[el], errors="coerce")
    sub = sub.merge(ss[["sample", "deseq2_group"]], left_on="SampleID", right_on="sample")
    group_means = sub.groupby("deseq2_group")[el].mean().dropna().sort_values()
    if a.ref_level not in group_means.index:
        continue
    rank = list(group_means.index).index(a.ref_level) + 1
    chem_ref_mean[el] = group_means[a.ref_level]
    chem_rank[el] = (rank, len(group_means))
    print(f"{el:10} {group_means[a.ref_level]:20.3f} {rank:10d} of {len(group_means):<6} {len(group_means):>10}")

low_conc_elements = {el for el, (rank, n) in chem_rank.items() if rank <= 2}
print(f"\nElements where Greece_Natural ranks in the bottom 2 (of however many groups): "
      f"{sorted(low_conc_elements)}")
print(f"Elements where it does NOT: {sorted(set(chem_rank) - low_conc_elements)}\n")

# ── Step 3: per-element positive vs negative test ────────────────────────
sig = pd.read_csv(a.maaslin2_sig, sep="\t")
sig["feature"] = sig["feature"].str.replace(".", "-", regex=False)

print("=== Per-element mirror test (positive vs negative association) ===")
print(f"{'Element':8} {'Dir':4} {'n':>4} {'%ref-origin':>12} {'baseline':>9} {'Fisher p (2-sided context)':>28}")
for el, grp in sig.groupby("element"):
    for direction, sub in [("POS", grp[grp["coef"] > 0]), ("NEG", grp[grp["coef"] < 0])]:
        genes = set(sub["feature"].unique()) & set(gene_origin.keys())
        if not genes:
            continue
        n_ref = sum(origin_is_ref[g] for g in genes)
        n_total = len(genes)
        pct = 100 * n_ref / n_total
        rest = set(gene_origin.keys()) - genes
        n_ref_rest = sum(origin_is_ref[g] for g in rest)
        table = [[n_ref, n_total - n_ref], [n_ref_rest, len(rest) - n_ref_rest]]
        _, p_more = fisher_exact(table, alternative="greater")
        _, p_less = fisher_exact(table, alternative="less")
        p_min = min(p_more, p_less)
        tag = "enriched" if p_more < p_less else "depleted"
        flag = f"  <-- {tag}" if p_min < 0.05 else ""
        low_conc_flag = " [Greece_Natural = low-conc]" if el in low_conc_elements else " [Greece_Natural NOT low-conc]"
        print(f"{el:8} {direction:4} {n_total:4d} {pct:11.2f}% {baseline_rate*100:8.2f}% "
              f"p={p_min:.2e}{flag}{low_conc_flag if direction=='POS' else ''}")
