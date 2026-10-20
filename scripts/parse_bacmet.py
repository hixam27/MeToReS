#!/usr/bin/env python3
"""parse_bacmet.py — turn DIAMOND-vs-BacMet hits into an ORF-keyed MRG table.

Takes the DIAMOND blastp outfmt-6 hits and the BacMet gene->compound mapping
file, keeps the single best hit per ORF, extracts the BacMet ID from the
subject id, and joins the gene name / compound so each ORF gets a metal/biocide
resistance annotation.
"""
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--hits", required=True, help="DIAMOND outfmt6 hits (bacmet_hits.tsv)")
ap.add_argument("--mapping", required=True, help="BacMet2_EXP.753.mapping.txt")
ap.add_argument("--out", required=True, help="output TSV (bacmet_annotations.tsv)")
a = ap.parse_args()

cols = ["qseqid", "sseqid", "pident", "length", "qcovhsp", "evalue", "bitscore"]
h = pd.read_csv(a.hits, sep="\t", names=cols)

if h.empty:
    pd.DataFrame(columns=["ORF_id", "BacMet_ID", "Gene_name", "Compound",
                          "pident", "qcovhsp", "evalue", "bitscore"]).to_csv(
        a.out, sep="\t", index=False)
    raise SystemExit("No BacMet hits passed thresholds.")

# Enforce single best hit per ORF.
h = h.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
# BacMet ID is the first '|'-delimited token of the subject id (e.g. BAC0001|...)
h["BacMet_ID"] = h["sseqid"].astype(str).str.split("|").str[0]

m = pd.read_csv(a.mapping, sep="\t")
# Column names in the mapping file vary by release — detect them defensively.
id_col   = next((c for c in m.columns if "bacmet" in c.lower() or c.lower() == "id"), m.columns[0])
gene_col = next((c for c in m.columns if "gene" in c.lower()), None)
comp_col = next((c for c in m.columns if "compound" in c.lower()), None)
keep = [id_col] + [c for c in (gene_col, comp_col) if c]
ren = {id_col: "BacMet_ID"}
if gene_col: ren[gene_col] = "Gene_name"
if comp_col: ren[comp_col] = "Compound"
m = m[keep].rename(columns=ren)

out = (h.merge(m, on="BacMet_ID", how="left")
         .rename(columns={"qseqid": "ORF_id"}))
final_cols = ["ORF_id", "BacMet_ID", "Gene_name", "Compound",
              "pident", "qcovhsp", "evalue", "bitscore"]
out = out[[c for c in final_cols if c in out.columns]]
out.to_csv(a.out, sep="\t", index=False)
print(f"BacMet MRG calls: {len(out)}")
