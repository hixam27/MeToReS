#!/usr/bin/env python3
"""safeguard_prevalence.py — S1: prevalence of genes behind association hits.

A gene detected in only two or three samples can produce a strong
gene–covariate association that is really presence/absence at one site.
With few samples per site this is the commonest way a single site
masquerades as a concentration–response.

Two subcommands:

  genes     list the genes that pass the prevalence rule. Used as the input
            filter for MaAsLin2 and partial correlation when
            safeguards.prevalence.apply_filter is true (the default).

  annotate  add, to a table of significant hits, the number and fraction of
            samples in which each gene is detected, and flag hits that fail
            the rule. Report only.

Prevalence rule: a gene passes if it is detected (abundance > 0) in at least
--min_fraction of all samples, OR in at least --min_in_group samples of one
group (--group_col). The group arm keeps genuinely site-specific genes once
sites have real replicates; with fewer replicates per group than
--min_in_group it never applies, and only the overall fraction counts.
"""
import argparse
import json
import math
import os
import sys

import pandas as pd


def detection(matrix_path):
    mat = pd.read_csv(matrix_path, index_col=0)
    mat.index = mat.index.astype(str)
    num = mat.apply(pd.to_numeric, errors="coerce")
    return num > 0


def add_rule_args(p):
    p.add_argument("--min_fraction", type=float, default=0.5)
    p.add_argument("--samplesheet", default=None)
    p.add_argument("--group_col", default=None)
    p.add_argument("--min_in_group", type=int, default=3)


ap = argparse.ArgumentParser()
sub = ap.add_subparsers(dest="cmd", required=True)

g = sub.add_parser("genes")
g.add_argument("--matrix", required=True)
add_rule_args(g)
g.add_argument("--out", required=True)

n = sub.add_parser("annotate")
n.add_argument("--matrix", required=True)
n.add_argument("--hits", required=True, help="significant-hit table (CSV or TSV)")
n.add_argument("--gene_col", required=True, help="column holding the gene id")
n.add_argument("--label", required=True, help="method label for the summary, e.g. MaAsLin2")
add_rule_args(n)
n.add_argument("--filter_applied", action="store_true",
               help="record that the same rule filtered the genes before testing")
n.add_argument("--out_tsv", required=True)
n.add_argument("--out_json", required=True)
a = ap.parse_args()

det = detection(a.matrix)
n_det, n_samples = det.sum(axis=1), det.shape[1]
min_n = math.ceil(a.min_fraction * n_samples)
passes = n_det >= min_n
max_in_group = pd.Series(0, index=det.index)
group_rule = "not used"
if a.samplesheet and a.group_col:
    ss = pd.read_csv(a.samplesheet)
    if a.group_col not in ss.columns:
        sys.exit(f"--group_col '{a.group_col}' not in samplesheet")
    grp = dict(zip(ss["sample"].astype(str), ss[a.group_col].astype(str)))
    cols = [c for c in det.columns if str(c) in grp]
    by_group = det[cols].T.groupby([grp[str(c)] for c in cols]).sum().T
    max_in_group = by_group.max(axis=1)
    sizes = pd.Series([grp[str(c)] for c in cols]).value_counts()
    eligible = int((sizes >= a.min_in_group).sum())
    group_rule = (f">= {a.min_in_group} samples within one '{a.group_col}' group "
                  f"({eligible} of {len(sizes)} groups have that many samples)")
    passes = passes | (max_in_group >= a.min_in_group)
rule_txt = f">= {min_n}/{n_samples} samples overall, or {group_rule}"

if a.cmd == "genes":
    keep = passes[passes].index
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("\n".join(keep) + ("\n" if len(keep) else ""))
    print(f"S1 prevalence filter: {len(keep)} of {len(n_det)} genes kept ({rule_txt})")
    sys.exit(0)

# annotate
sep = "\t" if a.hits.endswith((".tsv", ".txt")) else ","
hits = pd.read_csv(a.hits, sep=sep)
if a.gene_col not in hits.columns:
    sys.exit(f"Column '{a.gene_col}' not in {a.hits}: {list(hits.columns)}")
ids = hits[a.gene_col].astype(str)
hits["n_samples_detected"] = ids.map(n_det).astype("Int64")
hits["fraction_detected"] = (hits["n_samples_detected"] / n_samples).round(3)
hits["max_detected_in_one_group"] = ids.map(max_in_group).astype("Int64")
hits["low_prevalence"] = (~ids.map(passes).fillna(False).astype(bool)).astype("boolean")
hits.loc[hits["n_samples_detected"].isna(), "low_prevalence"] = pd.NA

os.makedirs(os.path.dirname(os.path.abspath(a.out_tsv)), exist_ok=True)
hits.to_csv(a.out_tsv, sep="\t", index=False)

n_hits = len(hits)
n_low = int(hits["low_prevalence"].fillna(False).sum())
n_unmatched = int(hits["n_samples_detected"].isna().sum())
genes = hits.drop_duplicates(a.gene_col)
summary = {
    "safeguard": "S1_prevalence",
    "method": a.label,
    "min_fraction": a.min_fraction,
    "min_samples": min_n,
    "rule": rule_txt,
    "filter_applied": bool(a.filter_applied),
    "n_samples": n_samples,
    "n_hits": n_hits,
    "n_hits_low_prevalence": n_low,
    "n_genes": int(len(genes)),
    "n_genes_low_prevalence": int(genes["low_prevalence"].fillna(False).sum()),
    "n_hits_not_in_matrix": n_unmatched,
    "detected_samples_distribution": {int(k): int(v) for k, v in
        hits["n_samples_detected"].dropna().value_counts().sort_index().items()},
}
with open(a.out_json, "w") as fh:
    json.dump(summary, fh, indent=2)

pct = 100 * n_low / n_hits if n_hits else 0
print(f"S1 prevalence ({a.label}): {n_low} of {n_hits} hits ({pct:.0f}%) come from genes "
      f"failing the prevalence rule ({rule_txt})")
if n_hits and pct >= 25:
    print(f"SAFEGUARD_WARNING: S1 {a.label}: {pct:.0f}% of hits are low-prevalence genes; "
          f"associations may reflect presence/absence at one site")
