#!/usr/bin/env python3
"""parse_card.py — tidy RGI (CARD) output into an ORF-keyed annotation table.

RGI protein mode writes a TSV whose query header is in 'ORF_ID' and whose
confidence tier is in 'Cut_Off' (Perfect / Strict / Loose). We keep only the
tiers requested (default Perfect,Strict) and emit a compact table keyed by the
catalog ORF id, so it can be left-joined onto the gene catalog downstream.
"""
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--rgi", required=True, help="RGI .txt output (card_rgi.txt)")
ap.add_argument("--keep", default="Perfect,Strict",
                help="comma-separated Cut_Off tiers to keep")
ap.add_argument("--out", required=True, help="output TSV (card_annotations.tsv)")
a = ap.parse_args()

keep = {k.strip() for k in a.keep.split(",")}

df = pd.read_csv(a.rgi, sep="\t")

# RGI protein mode: query header is in 'ORF_ID'; confidence tier in 'Cut_Off'.
if "Cut_Off" in df.columns:
    df = df[df["Cut_Off"].isin(keep)].copy()

# The ORF id is the first whitespace-delimited token of the FASTA header.
df["ORF_id"] = df["ORF_ID"].astype(str).str.split().str[0]

rename = {
    "ORF_id": "ORF_id",
    "Best_Hit_ARO": "Best_Hit_ARO",
    "ARO": "ARO",
    "Drug Class": "Drug_Class",
    "Resistance Mechanism": "Resistance_Mechanism",
    "AMR Gene Family": "AMR_Gene_Family",
    "Cut_Off": "Cut_Off",
    "Best_Identities": "pident",
}
cols = [c for c in rename if c in df.columns]
out = df[cols].rename(columns=rename).drop_duplicates("ORF_id")
out.to_csv(a.out, sep="\t", index=False)
print(f"CARD ARG calls ({a.keep}): {len(out)}")
