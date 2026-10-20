#!/usr/bin/env python3
"""resistome_coverage_report.py — how much of each annotated resistance
category is actually expressed.

For every BacMet compound (MRG), CARD best-hit ARO and AMR gene family
(ARG), reports:

  n_genes_annotated                     genes annotated to the category
  n_genes_detected_any_sample           genes with abundance > 0 in >= 1 sample
  n_genes_detected_geqN_samples         genes with abundance > 0 in >= N samples
  n_samples_with_any_signal             samples where the category total is > 0
  passes_correlation_floor              n_samples_with_any_signal >= N, i.e.
                                        whether resistance_chemical_correlation.R
                                        tests the category at all

N is --min_detected_samples, the same floor the correlation step applies.
"Detected in one sample" is nearly always true in a deep metatranscriptome;
the per-gene >= N figure tells a reader whether a correlation on that
category rests on genes expressed across samples or on a few sporadic ones.

Sample columns are taken from the samplesheet, so annotation columns are
never mistaken for samples.
"""
import argparse

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--abundance_with_resistance", required=True)
ap.add_argument("--samplesheet", required=True)
ap.add_argument("--min_detected_samples", type=int, default=4)
ap.add_argument("--out_csv", required=True)
a = ap.parse_args()

df = pd.read_csv(a.abundance_with_resistance, low_memory=False)
ss = pd.read_csv(a.samplesheet)
sample_cols = [s for s in ss["sample"] if s in df.columns]
if not sample_cols:
    raise SystemExit("No samplesheet sample appears as a column of "
                     f"{a.abundance_with_resistance}")
print(f"{len(sample_cols)} sample columns; detection floor = "
      f"{a.min_detected_samples} samples")
N = a.min_detected_samples


def coverage_report(subset_df, category_col, label):
    if category_col not in subset_df.columns or subset_df.empty:
        print(f"{label}: no '{category_col}' data — skipped")
        return pd.DataFrame()
    sub = subset_df.copy()
    sub[category_col] = sub[category_col].fillna("unclassified").astype(str)
    sep = r";" if label == "ARG_family" else r"[;,]"
    sub = sub.assign(Category=sub[category_col].str.split(sep)).explode("Category")
    sub["Category"] = sub["Category"].str.strip()

    rows = []
    for cat, grp in sub.groupby("Category"):
        abund = grp[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
        per_gene = (abund > 0).sum(axis=1)
        n_genes = len(grp)
        n_any = int((per_gene > 0).sum())
        n_floor = int((per_gene >= N).sum())
        n_samples_signal = int((abund.sum(axis=0) > 0).sum())
        rows.append({
            "Source": label,
            "Category": cat,
            "n_genes_annotated": n_genes,
            "n_genes_detected_any_sample": n_any,
            "pct_genes_detected_any_sample": round(100 * n_any / n_genes, 1) if n_genes else 0,
            f"n_genes_detected_geq{N}_samples": n_floor,
            f"pct_genes_detected_geq{N}_samples": round(100 * n_floor / n_genes, 1) if n_genes else 0,
            "n_samples_with_any_signal": n_samples_signal,
            "passes_correlation_floor": n_samples_signal >= N,
        })
    return pd.DataFrame(rows)


def flag(col):
    if col not in df.columns:
        return pd.Series(False, index=df.index)
    return df[col].astype(str).str.lower().isin(["true", "1"])


arg_col = "Best_Hit_ARO" if "Best_Hit_ARO" in df.columns else "ARO"
reports = [
    coverage_report(df[flag("is_MRG")], "Compound", "MRG"),
    coverage_report(df[flag("is_ARG")], arg_col, "ARG"),
    coverage_report(df[flag("is_ARG")], "AMR_Gene_Family", "ARG_family"),
]
out = pd.concat([r for r in reports if not r.empty], ignore_index=True)
if not out.empty:
    out = out.sort_values(["Source", "n_genes_annotated"], ascending=[True, False])
out.to_csv(a.out_csv, index=False)

print(f"Saved: {a.out_csv} ({len(out)} categories)")
for src, grp in out.groupby("Source") if not out.empty else []:
    print(f"  {src}: {len(grp)} categories, {int(grp['passes_correlation_floor'].sum())} "
          f"pass the {N}-sample floor; {int(grp['n_genes_annotated'].sum())} genes annotated, "
          f"{int(grp[f'n_genes_detected_geq{N}_samples'].sum())} detected in >= {N} samples")
