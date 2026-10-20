#!/usr/bin/env python3
"""parse_card.py — tidy RGI (CARD) output into an ORF-keyed annotation table.

RGI protein mode writes a TSV whose query header is in 'ORF_ID', whose
confidence tier is in 'Cut_Off' (Perfect / Strict / Loose) and whose CARD
model is in 'Model_type'. Three filters are applied:

  --keep        confidence tiers to keep (default Perfect,Strict)
  --models      CARD model types counted as ARG calls (default: protein
                homolog model). Other model types — protein variant,
                overexpression and rRNA models — are only meaningful when
                specific mutations are present, which a homology screen of
                assembled transcripts does not establish. Their hits are
                written to a separate table (--out_other) instead of being
                mixed into the ARG calls.
  --min_pident  minimum percent identity to the CARD reference (default 0:
                rely on CARD's curated per-model Strict/Perfect bit-score
                cut-offs). Set e.g. 80 for a high-identity sensitivity check.

Every kept hit gets an 'identity_tier' (<40, 40-60, 60-80, >=80 %), and the
tier counts are printed: low-identity Strict hits are distant homologs
(often efflux pumps, transporters, regulators), and a report should say how
much of the resistome sits there.

The output is keyed by catalogue ORF id, for joining onto the gene catalogue.
"""
import argparse

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--rgi", required=True, help="RGI .txt output (card_rgi.txt)")
ap.add_argument("--keep", default="Perfect,Strict",
                help="comma-separated Cut_Off tiers to keep")
ap.add_argument("--models", default="protein homolog model",
                help="comma-separated Model_type values counted as ARG calls")
ap.add_argument("--min_pident", type=float, default=0.0,
                help="minimum Best_Identities (percent) to keep a hit")
ap.add_argument("--out", required=True, help="output TSV (card_annotations.tsv)")
ap.add_argument("--out_other", required=True,
                help="output TSV for hits from other model types")
a = ap.parse_args()

keep = {k.strip() for k in a.keep.split(",") if k.strip()}
models = {m.strip().lower() for m in a.models.split(",") if m.strip()}

df = pd.read_csv(a.rgi, sep="\t")
n_raw = len(df)

if "Cut_Off" in df.columns:
    df = df[df["Cut_Off"].isin(keep)].copy()
n_tier = len(df)

if "Best_Identities" in df.columns:
    df = df[pd.to_numeric(df["Best_Identities"], errors="coerce") >= a.min_pident].copy()
n_pid = len(df)

# The ORF id is the first whitespace-delimited token of the FASTA header.
df["ORF_id"] = df["ORF_ID"].astype(str).str.split().str[0]

rename = {
    "ORF_id": "ORF_id",
    "Best_Hit_ARO": "Best_Hit_ARO",
    "ARO": "ARO",
    "Drug Class": "Drug_Class",
    "Resistance Mechanism": "Resistance_Mechanism",
    "AMR Gene Family": "AMR_Gene_Family",
    "Model_type": "Model_type",
    "Cut_Off": "Cut_Off",
    "Best_Identities": "pident",
}
cols = [c for c in rename if c in df.columns]
tidy = df[cols].rename(columns=rename)

if "Model_type" in tidy.columns:
    is_kept_model = tidy["Model_type"].astype(str).str.strip().str.lower().isin(models)
else:
    is_kept_model = pd.Series(True, index=tidy.index)

if "pident" in tidy.columns:
    tidy["identity_tier"] = pd.cut(pd.to_numeric(tidy["pident"], errors="coerce"),
                                   [-0.1, 40, 60, 80, 100.1], right=False,
                                   labels=["<40", "40-60", "60-80", ">=80"]).astype(str)
out = tidy[is_kept_model].drop_duplicates("ORF_id")
other = tidy[~is_kept_model].drop_duplicates("ORF_id")
out.to_csv(a.out, sep="\t", index=False)
other.to_csv(a.out_other, sep="\t", index=False)

print(f"RGI hits: {n_raw} raw -> {n_tier} in tiers {sorted(keep)} -> "
      f"{n_pid} with identity >= {a.min_pident:g}%")
print(f"CARD ARG calls ({', '.join(sorted(models))}): {len(out)}")
if "identity_tier" in out.columns and len(out):
    tiers = out["identity_tier"].value_counts()
    print("  identity tiers: " + ", ".join(
        f"{t}%: {int(tiers.get(t, 0))} ({100 * tiers.get(t, 0) / len(out):.0f}%)"
        for t in ["<40", "40-60", "60-80", ">=80"]))
if len(other):
    print(f"Hits from other model types, reported separately: {len(other)} "
          f"({other['Model_type'].value_counts().to_dict()})")
