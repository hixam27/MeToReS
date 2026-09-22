#!/usr/bin/env python3
"""resistome_coverage_report.py — for each BacMet compound / CARD ARO family,
reports how many annotated genes exist versus how many actually have
detectable (nonzero) expression in at least one sample. Per Lorenzo's
follow-up: makes the Fe/Cu-style numbers auditable rather than a black box,
on top of (not instead of) the merge_resistance.py matrix-source fix.

Uses abundance_with_resistance.csv (annotation + per-sample abundance
already joined) so gene-level coverage is computed directly, independent
of the aggregate ARO/compound tables built downstream from it.
"""
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--abundance_with_resistance", required=True)
ap.add_argument("--out_csv", required=True)
a = ap.parse_args()

df = pd.read_csv(a.abundance_with_resistance, low_memory=False)

# Confirmed exhaustive list from the real abundance_with_resistance.csv
# header — anything NOT in this set is treated as a sample column. An
# earlier, guessed version of this list was incomplete (missing ARO,
# Cut_Off, pident_x/y, BacMet_ID, qcovhsp, evalue, bitscore, status),
# which silently miscounted those as sample columns. Several of them
# (pident_x, bitscore, evalue, qcovhsp) are always nonzero for any
# annotated row by construction, which made every category appear
# trivially "100% detected" regardless of real per-sample abundance —
# caught by the sample count itself coming out wrong (23, not 14).
meta_like = {"ORF_id", "Best_Hit_ARO", "ARO", "Drug_Class", "Resistance_Mechanism",
             "AMR_Gene_Family", "Cut_Off", "pident_x", "is_ARG", "BacMet_ID",
             "Gene_name", "Compound", "pident_y", "qcovhsp", "evalue", "bitscore",
             "is_MRG", "status", "superkingdom", "phylum", "class", "order",
             "family", "genus", "species"}
sample_cols = [c for c in df.columns if c not in meta_like]
print(f"Detected {len(sample_cols)} sample columns: {sample_cols}")


def coverage_report(subset_df, category_col, label, min_detected_samples=4):
    sub = subset_df.copy()
    if category_col not in sub.columns:
        print(f"'{category_col}' not found — skipping {label}")
        return pd.DataFrame()
    sub[category_col] = sub[category_col].fillna("unclassified").astype(str)
    sub = sub.assign(Category=sub[category_col].str.split(r"[;,]")).explode("Category")
    sub["Category"] = sub["Category"].str.strip()

    rows = []
    for cat, grp in sub.groupby("Category"):
        n_genes = len(grp)
        abund = grp[sample_cols].apply(pd.to_numeric, errors="coerce")
        # "Detected" (nonzero in >=1 sample) is close to trivially true in a
        # deep metatranscriptome — not the figure a reader needs to judge
        # whether a correlation on this category is well-supported. The
        # stricter figure below matches the actual --min_detected_samples
        # threshold resistance_chemical_correlation.R requires before a
        # category is even tested, i.e. per-gene adequacy, not a per-category
        # "any signal anywhere" statement.
        n_samples_detected_per_gene = (abund > 0).sum(axis=1)
        n_genes_detected = (n_samples_detected_per_gene > 0).sum()
        n_genes_min_detected = (n_samples_detected_per_gene >= min_detected_samples).sum()
        n_samples_with_signal = (abund.sum(axis=0) > 0).sum()
        rows.append({
            "Source": label,
            "Category": cat,
            "n_genes_annotated": n_genes,
            "n_genes_detected_any_sample": n_genes_detected,
            "pct_genes_detected_any_sample": round(100 * n_genes_detected / n_genes, 1) if n_genes else 0,
            f"n_genes_detected_geq{min_detected_samples}_samples": n_genes_min_detected,
            f"pct_genes_detected_geq{min_detected_samples}_samples": round(100 * n_genes_min_detected / n_genes, 1) if n_genes else 0,
            "n_samples_with_any_signal": n_samples_with_signal,
        })
    return pd.DataFrame(rows)


def bool_col(df, col):
    if col in df.columns:
        return df[col] == True
    return pd.Series(False, index=df.index)


mrg_report = coverage_report(df[bool_col(df, "is_MRG")], "Compound", "MRG")
arg_col = "Best_Hit_ARO" if "Best_Hit_ARO" in df.columns else "ARO"
arg_report = coverage_report(df[bool_col(df, "is_ARG")], arg_col, "ARG")

out = pd.concat([mrg_report, arg_report], ignore_index=True)
out = out.sort_values(["Source", "n_genes_annotated"], ascending=[True, False])
out.to_csv(a.out_csv, index=False)

print(f"\nSaved: {a.out_csv} ({len(out)} categories)")
print("\n=== MRG coverage (compounds) ===")
print(mrg_report.sort_values("n_genes_annotated", ascending=False).to_string(index=False))
