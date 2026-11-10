#!/usr/bin/env python3
"""select_coding_genes.py — choose the genes that enter differential expression.

A de novo catalogue also holds open reading frames that Prodigal predicts on
non-coding RNA (rRNA, tRNA/tmRNA, other ncRNAs) and on library artefacts. In
total-RNA libraries these can carry most of the counts of the filtered genes.

  --gene_set coding  keep the filtered genes with protein-level evidence, i.e.
                     a hit in the eggNOG-mapper annotations (default)
  --gene_set all     copy the filtered matrix unchanged (behaviour before v1.1.0)

The selected matrix keeps the input's lines byte for byte (header plus the
selected rows), so the "all" mode output is identical to its input.

The JSON summary is written in both modes: genes and Salmon counts (NumReads)
in the filtered set and in its protein-evidence subset, overall and per sample.
"""
import argparse
import json
import os
import shutil
import sys

import pandas as pd

ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("--matrix", required=True,
                help="filtered matrix (transcript_abundance_quantification_table_filter.csv)")
ap.add_argument("--annotations", required=True, help="eggNOG-mapper .emapper.annotations")
ap.add_argument("--samplesheet", required=True, help="samplesheet with 'sample' and 'sra' columns")
ap.add_argument("--quant_dir", required=True, help="directory with <sra>_quant/quant.sf")
ap.add_argument("--gene_set", choices=["coding", "all"], default="coding")
ap.add_argument("--out_matrix", required=True)
ap.add_argument("--out_json", required=True)
a = ap.parse_args()


def annotated_ids(path):
    """Query IDs of the eggNOG-mapper annotation rows (comment lines start with '#')."""
    ids = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            ids.add(line.split("\t", 1)[0].strip())
    return ids


def row_id(line):
    return line.split(",", 1)[0].strip().strip('"')


annot = annotated_ids(a.annotations)

with open(a.matrix) as fh:
    header = fh.readline()
    rows = fh.readlines()
filtered = [row_id(r) for r in rows]
filtered_set = set(filtered)
coding_set = filtered_set & annot

if a.gene_set == "coding" and not coding_set:
    sys.exit("ERROR: no filtered gene has an eggNOG-mapper hit, so no gene would enter "
             "DESeq2. Check the eggNOG-mapper output, or set deseq2.gene_set: \"all\".")

os.makedirs(os.path.dirname(os.path.abspath(a.out_matrix)), exist_ok=True)
if a.gene_set == "all":
    shutil.copyfile(a.matrix, a.out_matrix)
    n_written = len(rows)
else:
    keep = [r for r, g in zip(rows, filtered) if g in coding_set]
    with open(a.out_matrix, "w") as fh:
        fh.write(header)
        fh.writelines(keep)
    n_written = len(keep)

# Salmon counts per sample, filtered set vs protein-evidence subset
ss = pd.read_csv(a.samplesheet)
per_sample = {}
tot_f = tot_c = 0.0
for _, r in ss.iterrows():
    qf = os.path.join(a.quant_dir, f"{r['sra']}_quant", "quant.sf")
    q = pd.read_csv(qf, sep="\t", usecols=["Name", "NumReads"])
    q = q[q["Name"].isin(filtered_set)]
    f_counts = float(q["NumReads"].sum())
    c_counts = float(q.loc[q["Name"].isin(coding_set), "NumReads"].sum())
    tot_f += f_counts
    tot_c += c_counts
    per_sample[str(r["sample"])] = {
        "counts_filtered": round(f_counts, 1),
        "counts_with_protein_evidence": round(c_counts, 1),
        "pct_counts_with_protein_evidence": round(100 * c_counts / f_counts, 1) if f_counts else None,
    }

pcts = [v["pct_counts_with_protein_evidence"] for v in per_sample.values()
        if v["pct_counts_with_protein_evidence"] is not None]
summary = {
    "gene_set": a.gene_set,
    "evidence": "eggNOG-mapper hit",
    "n_filtered": len(filtered_set),
    "n_with_protein_evidence": len(coding_set),
    "pct_genes_with_protein_evidence": round(100 * len(coding_set) / len(filtered_set), 1)
                                       if filtered_set else None,
    "n_genes_selected": n_written,
    "counts_filtered": round(tot_f, 1),
    "counts_with_protein_evidence": round(tot_c, 1),
    "pct_counts_with_protein_evidence": round(100 * tot_c / tot_f, 1) if tot_f else None,
    "pct_counts_min": min(pcts) if pcts else None,
    "pct_counts_max": max(pcts) if pcts else None,
    "per_sample": per_sample,
}
with open(a.out_json, "w") as fh:
    json.dump(summary, fh, indent=2)

print(f"Gene set '{a.gene_set}': {n_written:,} of {len(filtered_set):,} filtered genes enter DESeq2")
print(f"Filtered genes with protein-level evidence (eggNOG-mapper hit): {len(coding_set):,} "
      f"({summary['pct_genes_with_protein_evidence']}% of genes, "
      f"{summary['pct_counts_with_protein_evidence']}% of counts; per sample "
      f"{summary['pct_counts_min']}-{summary['pct_counts_max']}%)")
