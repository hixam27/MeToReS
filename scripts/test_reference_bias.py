#!/usr/bin/env python3
"""test_reference_bias.py — tests Lorenzo's "pooled-catalogue reference
bias" hypothesis: are down-regulated genes (higher in Greece_Natural)
disproportionately genes that were originally ASSEMBLED FROM a
Greece_Natural sample, rather than genuine differential expression?

Every gene header in the final catalog retains its original SRA-sample
prefix (added by individual_assembly, preserved through CD-HIT-EST
clustering), so origin-per-gene is directly recoverable without
recomputing anything.

Up-regulated gene sets serve as a built-in negative control: if the
Greece_Natural-origin enrichment is specific to the down-regulated
("higher in reference") direction, that's strong evidence for the
reference-bias mechanism, not a general artifact of gene-set size.
"""
import argparse
import os
import glob
import pandas as pd
from scipy.stats import fisher_exact

ap = argparse.ArgumentParser()
ap.add_argument("--catalog", required=True)
ap.add_argument("--samplesheet", required=True)
ap.add_argument("--deg_dir", required=True, help="DEG_RESULTS dir containing per-contrast subfolders")
ap.add_argument("--ref_level", default="Greece_Natural")
a = ap.parse_args()

# ── Step 1: map each catalog gene to its origin sample ──────────────────
ss = pd.read_csv(a.samplesheet)
sra_to_sample = dict(zip(ss["sra"], ss["sample"]))
# Sort SRA IDs longest-first so a prefix match can't accidentally match a
# shorter SRA ID that happens to be a substring of a longer one.
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
    print(f"WARNING: {unmatched:,} headers did not match any known SRA prefix "
          f"— check --catalog/--samplesheet are from the same run.")

ref_samples = set(ss.loc[ss["deseq2_group"] == a.ref_level, "sample"])
print(f"'{a.ref_level}' samples: {sorted(ref_samples)}")

origin_is_ref = {g: (s in ref_samples) for g, s in gene_origin.items()}
baseline_rate = sum(origin_is_ref.values()) / len(origin_is_ref)
print(f"Catalogue-wide baseline: {baseline_rate*100:.2f}% of genes "
      f"originally assembled from a {a.ref_level} sample\n")

# ── Step 2: per-contrast enrichment test ─────────────────────────────────
def load_ids(path):
    if not os.path.exists(path):
        return set()
    with open(path) as fh:
        return set(l.strip() for l in fh if l.strip())

contrast_dirs = sorted(glob.glob(os.path.join(a.deg_dir, f"*_vs_{a.ref_level}")))

print(f"{'Contrast':35} {'Set':6} {'n':>6} {'%ref-origin':>12} {'baseline':>10} {'Fisher p':>10}")
for cdir in contrast_dirs:
    contrast = os.path.basename(cdir)
    for direction, fname in [("DOWN", "downregulated_genes_id.txt"),
                              ("UP", "upregulated_genes_id.txt")]:
        ids = load_ids(os.path.join(cdir, fname))
        ids = ids & set(gene_origin.keys())
        if not ids:
            continue
        n_ref = sum(origin_is_ref[g] for g in ids)
        n_total = len(ids)
        pct = 100 * n_ref / n_total

        # 2x2 table: [in-set & ref-origin, in-set & not-ref-origin],
        #            [rest-of-catalog & ref-origin, rest & not-ref-origin]
        rest = set(gene_origin.keys()) - ids
        n_ref_rest = sum(origin_is_ref[g] for g in rest)
        table = [[n_ref, n_total - n_ref],
                 [n_ref_rest, len(rest) - n_ref_rest]]
        _, p = fisher_exact(table, alternative="greater")

        flag = "  <-- enriched" if p < 0.05 and pct > baseline_rate * 100 else ""
        print(f"{contrast:35} {direction:6} {n_total:6d} {pct:11.2f}% {baseline_rate*100:9.2f}% {p:10.2e}{flag}")
    print()
