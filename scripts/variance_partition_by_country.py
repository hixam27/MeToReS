#!/usr/bin/env python3
"""variance_partition_by_country.py — for each measured element, decomposes
total variance into between-country and within-country components (simple
one-way ANOVA eta-squared), per Lorenzo's §3 follow-up.

If chemistry is overwhelmingly a between-country contrast, fitting country
as a random effect (as MaAsLin2 does) absorbs most of the real chemical
signal along with the assembly-origin confound — leaving little genuine
within-country variance for the model to detect anything from. This
script quantifies that split directly rather than leaving it assumed.
"""
import argparse
import pandas as pd
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--chemical", required=True, help="chemical_data_<fraction>.csv")
ap.add_argument("--group_col", default="Country")
a = ap.parse_args()

df = pd.read_csv(a.chemical)
if a.group_col not in df.columns:
    raise SystemExit(f"'{a.group_col}' not found. Available: {list(df.columns)}")

meta_cols = {"SampleID", "Country", "Season", "Replicate"}
elem_cols = [c for c in df.columns if c not in meta_cols]

print(f"{'Element':10} {'n':>4} {'n_missing':>10} {'SS_between':>12} {'SS_within':>12} "
      f"{'eta2 (between %)':>18}")

results = []
for elem in elem_cols:
    sub = df[[a.group_col, elem]].copy()
    sub[elem] = pd.to_numeric(sub[elem], errors="coerce")
    n_missing = sub[elem].isna().sum()
    sub = sub.dropna()
    if sub[a.group_col].nunique() < 2 or len(sub) < 3:
        print(f"{elem:10} {len(sub):4d} {n_missing:10d} {'--':>12} {'--':>12} "
              f"{'insufficient data':>18}")
        continue

    grand_mean = sub[elem].mean()
    ss_total = ((sub[elem] - grand_mean) ** 2).sum()

    group_means = sub.groupby(a.group_col)[elem].mean()
    group_sizes = sub.groupby(a.group_col)[elem].size()
    ss_between = sum(group_sizes[g] * (group_means[g] - grand_mean) ** 2 for g in group_means.index)
    ss_within = ss_total - ss_between

    eta2 = ss_between / ss_total if ss_total > 0 else float("nan")

    print(f"{elem:10} {len(sub):4d} {n_missing:10d} {ss_between:12.2f} {ss_within:12.2f} "
          f"{eta2*100:17.1f}%")
    results.append({"Element": elem, "n": len(sub), "eta2_between_country": eta2})

if results:
    res_df = pd.DataFrame(results).sort_values("eta2_between_country", ascending=False)
    mean_eta2 = res_df["eta2_between_country"].mean()
    print(f"\nMean between-country variance share across {len(res_df)} elements: "
          f"{mean_eta2*100:.1f}%")
    print("\nHighest between-country elements (chemistry ~ entirely a country contrast):")
    print(res_df.head(5).to_string(index=False))
    print("\nLowest between-country elements (more genuine within-country variation):")
    print(res_df.tail(5).to_string(index=False))
