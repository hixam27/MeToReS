#!/usr/bin/env python3
"""merge_resistance.py — merge ARG + MRG annotations onto the catalog and the
abundance matrix, then aggregate expression per ARO and per compound.

Outputs:
  --out_annotation  : per-ORF resistance annotation (ARG/MRG flags + details)
  --out_matrix      : resistance subset of the abundance matrix + annotation
  --out_by_aro      : ARG expression summed per best-hit ARO term
  --out_by_family   : ARG expression summed per AMR gene family (coarser,
                      more robust: near-identical ARO terms are pooled)
  --out_by_compound : MRG expression summed per metal/biocide compound
"""
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
# Unfiltered TPM matrix (transcript_abundance_quantification_table.csv), not
# the expression-filtered one: many resistance genes are expressed at low
# levels and would otherwise be dropped before annotation. Column 0 = ORF id.
ap.add_argument("--matrix", required=True)
ap.add_argument("--card", required=True)
ap.add_argument("--bacmet", required=True)
ap.add_argument("--taxonomy")               # optional Kaiju gene taxonomy
ap.add_argument("--emapper")                # optional eggNOG (context only)
ap.add_argument("--out_annotation", required=True)
ap.add_argument("--out_matrix", required=True)
ap.add_argument("--out_by_aro", required=True)
ap.add_argument("--out_by_family", required=True)
ap.add_argument("--out_by_compound", required=True)
a = ap.parse_args()

mat = pd.read_csv(a.matrix)
mat = mat.rename(columns={mat.columns[0]: "ORF_id"})
sample_cols = [c for c in mat.columns if c != "ORF_id"]

card = pd.read_csv(a.card, sep="\t");   card["is_ARG"] = True
bac  = pd.read_csv(a.bacmet, sep="\t"); bac["is_MRG"]  = True

# Both tables carry a percent-identity column named 'pident'. Give each a
# source-specific name before the outer merge, so the output does not end up
# with ambiguous pident_x / pident_y columns.
card = card.rename(columns={"pident": "card_pident"})
bac  = bac.rename(columns={"pident": "bacmet_pident"})
shared = (set(card.columns) & set(bac.columns)) - {"ORF_id"}
if shared:
    raise SystemExit(f"Unexpected columns shared by CARD and BacMet tables: {sorted(shared)}")

ann = card.merge(bac, on="ORF_id", how="outer")
ann["is_ARG"] = ann["is_ARG"].fillna(False)
ann["is_MRG"] = ann["is_MRG"].fillna(False)

# Optional taxonomic context (who carries the gene)
if a.taxonomy:
    tax = pd.read_csv(a.taxonomy, sep="\t")
    # Detect the gene/ORF id column defensively — don't assume column 0,
    # and don't assume its dtype matches ann["ORF_id"] (string).
    id_col = next((c for c in tax.columns
                   if c.lower() in ("orf_id", "gene_id", "geneid", "id", "name")),
                  tax.columns[0])
    tax = tax.rename(columns={id_col: "ORF_id"})
    tax["ORF_id"] = tax["ORF_id"].astype(str)
    ann["ORF_id"] = ann["ORF_id"].astype(str)
    ann = ann.merge(tax, on="ORF_id", how="left")

# (eggNOG join intentionally left out of the keyed table: column indices are
#  version-dependent. If wanted, read a.emapper with comment="#", header=None,
#  rename col 0 -> ORF_id and pull the KEGG_ko / Description columns for your
#  emapper version, then merge on ORF_id here.)

ann.to_csv(a.out_annotation, sep="\t", index=False)

# Attach per-sample abundances to the resistance subset
amat = ann.merge(mat, on="ORF_id", how="left")
amat.to_csv(a.out_matrix, index=False)

# Aggregate ARG expression per ARO (antibiotic resistance ontology term)
arg = amat[amat["is_ARG"]]
aro_key = "Best_Hit_ARO" if "Best_Hit_ARO" in arg.columns else "ARO"
if aro_key in arg.columns and not arg.empty:
    arg.groupby(aro_key)[sample_cols].sum().reset_index().to_csv(a.out_by_aro, index=False)
else:
    pd.DataFrame(columns=[aro_key] + sample_cols).to_csv(a.out_by_aro, index=False)

# Aggregate ARG expression per AMR gene family (a gene can list several,
# separated by ';' -> counted under each)
fam = arg.copy()
if "AMR_Gene_Family" in fam.columns and not fam.empty:
    fam["AMR_Gene_Family"] = fam["AMR_Gene_Family"].fillna("unclassified").astype(str)
    fam = fam.assign(AMR_Gene_Family=fam["AMR_Gene_Family"].str.split(";")).explode("AMR_Gene_Family")
    fam["AMR_Gene_Family"] = fam["AMR_Gene_Family"].str.strip()
    fam.groupby("AMR_Gene_Family")[sample_cols].sum().reset_index().to_csv(a.out_by_family, index=False)
else:
    pd.DataFrame(columns=["AMR_Gene_Family"] + sample_cols).to_csv(a.out_by_family, index=False)

# Aggregate MRG expression per compound (Compound may list several -> explode)
mrg = amat[amat["is_MRG"]].copy()
if "Compound" in mrg.columns and not mrg.empty:
    mrg["Compound"] = mrg["Compound"].fillna("unclassified").astype(str)
    mrg = mrg.assign(Compound=mrg["Compound"].str.split(r"[;,]")).explode("Compound")
    mrg["Compound"] = mrg["Compound"].str.strip()
    mrg.groupby("Compound")[sample_cols].sum().reset_index().to_csv(a.out_by_compound, index=False)
else:
    pd.DataFrame(columns=["Compound"] + sample_cols).to_csv(a.out_by_compound, index=False)

print(f"Resistance annotation: {int(ann['is_ARG'].sum())} ARG, "
      f"{int(ann['is_MRG'].sum())} MRG ORFs")
