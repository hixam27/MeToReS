#!/usr/bin/env python3
"""safeguard_library_complexity.py — S2: flag low-complexity libraries.

A library with very few distinct transcripts (for example, one dominated by
a handful of highly abundant sequences) can still have a normal read count,
so depth-based QC does not catch it. Entering such a sample as a replicate
distorts normalisation and produces presence/absence "differential
expression" driven by that one sample.

For every sample this reports three measures, each relative to the run
median, and flags the sample if any falls below --warn_fraction:

  genes_detected  genes with abundance > 0 in the quantification matrix
  orfs_assembled  ORFs predicted from the sample's own assembly
  reads_mapped    reads mapped by Salmon (num_mapped in meta_info.json)

It also records the library type each sample was quantified with
(samplesheet column 'libtype', else --default_libtype). A DESeq2 contrast
whose samples mix library types is
flagged: library type then acts as a batch effect (unstranded libraries also
count antisense reads). It can be added to the design formula
(e.g. "~ libtype + group") once replication allows.

With --probe_dir, each declared library type is also compared with the type
Salmon auto-detects on a read subsample ({probe_dir}/{sra}.txt, written by
the salmon_libtype_probe rule). A mismatch is flagged: a stranded type
declared for an unstranded library discards about half of its fragments,
and the reverse loses strand information.

Report only: nothing is excluded. To exclude a flagged sample, remove its
row from the samplesheet and re-run.
"""
import argparse
import json
import os
import sys

import pandas as pd


def count_fasta_records(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def salmon_mapped(quant_dir):
    meta = os.path.join(quant_dir, "aux_info", "meta_info.json")
    if not os.path.exists(meta):
        return None
    with open(meta) as fh:
        m = json.load(fh)
    return m.get("num_mapped")


ap = argparse.ArgumentParser()
ap.add_argument("--samplesheet", required=True)
ap.add_argument("--matrix", required=True,
                help="abundance matrix; first column gene id, one column per sample")
ap.add_argument("--assembly_dir", required=True,
                help="directory holding {sra}/orfs.fasta")
ap.add_argument("--quant_dir", required=True,
                help="directory holding {sra}_quant/aux_info/meta_info.json")
ap.add_argument("--warn_fraction", type=float, default=0.25)
ap.add_argument("--default_libtype", default="ISR")
ap.add_argument("--probe_dir", default=None,
                help="directory holding {sra}.txt with Salmon's auto-detected library type")
ap.add_argument("--deseq_group_col", default=None)
ap.add_argument("--contrasts", default="", help="space-separated DESeq2 contrasts (A_vs_B)")
ap.add_argument("--out_tsv", required=True)
ap.add_argument("--out_json", required=True)
ap.add_argument("--out_mqc", required=True, help="MultiQC custom-content YAML")
a = ap.parse_args()

ss = pd.read_csv(a.samplesheet)
mat = pd.read_csv(a.matrix, index_col=0)
missing = [s for s in ss["sample"] if s not in mat.columns]
if missing:
    sys.exit(f"Samples missing from {a.matrix}: {missing}")

rows = []
for _, r in ss.iterrows():
    s, sra = r["sample"], r["sra"]
    orfs = os.path.join(a.assembly_dir, str(sra), "orfs.fasta")
    rows.append({
        "sample": s,
        "sra": sra,
        "genes_detected": int((pd.to_numeric(mat[s], errors="coerce") > 0).sum()),
        "orfs_assembled": count_fasta_records(orfs) if os.path.exists(orfs) else None,
        "reads_mapped": salmon_mapped(os.path.join(a.quant_dir, f"{sra}_quant")),
        "libtype": (str(r["libtype"]).strip() if "libtype" in ss.columns and pd.notna(r["libtype"])
                    and str(r["libtype"]).strip() else a.default_libtype),
    })
df = pd.DataFrame(rows)


def detected_libtype(sra):
    if not a.probe_dir:
        return None
    path = os.path.join(a.probe_dir, f"{sra}.txt")
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        v = fh.read().strip()
    return v or None


df["libtype_detected"] = [detected_libtype(sra) for sra in df["sra"]]
df["libtype_mismatch"] = [
    bool(det) and decl.upper() != "A" and det != decl
    for decl, det in zip(df["libtype"], df["libtype_detected"])
]

metrics = ["genes_detected", "orfs_assembled", "reads_mapped"]
flag_cols = []
for m in metrics:
    vals = pd.to_numeric(df[m], errors="coerce")
    med = vals.median()
    df[f"{m}_pct_of_median"] = (vals / med * 100).round(1) if med and med > 0 else None
    df[f"{m}_low"] = vals < a.warn_fraction * med if med and med > 0 else False
    flag_cols.append(f"{m}_low")

df["flagged"] = df[flag_cols].any(axis=1)
df["reason"] = df.apply(
    lambda r: "; ".join(f"{m} {r[f'{m}_pct_of_median']}% of median"
                        for m in metrics if r[f"{m}_low"]), axis=1)

os.makedirs(os.path.dirname(os.path.abspath(a.out_tsv)), exist_ok=True)
df.to_csv(a.out_tsv, sep="\t", index=False)

# Contrasts that mix library types
mixed = []
if a.deseq_group_col and a.deseq_group_col in ss.columns and a.contrasts.strip():
    g = dict(zip(ss["sample"], ss[a.deseq_group_col].astype(str)))
    lt = dict(zip(df["sample"], df["libtype"]))
    for c in a.contrasts.split():
        if "_vs_" not in c:
            continue
        groups = c.split("_vs_", 1)
        samples = [s for s in df["sample"] if g.get(s) in groups]
        per_group = {grp: sorted({lt[s] for s in samples if g[s] == grp}) for grp in groups}
        if len({lt[s] for s in samples}) > 1:
            mixed.append({"contrast": c, "libtypes": per_group})

flagged = df[df["flagged"]]
summary = {
    "safeguard": "S2_library_complexity",
    "warn_fraction_of_median": a.warn_fraction,
    "n_samples": int(len(df)),
    "n_flagged": int(len(flagged)),
    "flagged": [{"sample": r["sample"], "reason": r["reason"]} for _, r in flagged.iterrows()],
    "libtypes": {r["sample"]: r["libtype"] for _, r in df.iterrows()},
    "contrasts_mixing_libtypes": mixed,
    "libtype_checked": bool(a.probe_dir),
    "libtype_mismatch": [{"sample": r["sample"], "declared": r["libtype"],
                          "detected": r["libtype_detected"]}
                         for _, r in df[df["libtype_mismatch"]].iterrows()],
    "medians": {m: (float(pd.to_numeric(df[m], errors="coerce").median())
                    if pd.to_numeric(df[m], errors="coerce").notna().any() else None)
                for m in metrics},
}
with open(a.out_json, "w") as fh:
    json.dump(summary, fh, indent=2)

# MultiQC custom content: one table, one row per sample
mqc = {
    "id": "metores_library_complexity",
    "section_name": "Library complexity (safeguard S2)",
    "description": (f"Genes detected, ORFs assembled and reads mapped per sample, as % of the "
                    f"run median. Samples below {a.warn_fraction*100:.0f}% on any measure are "
                    f"flagged: they can act as single-sample drivers of differential results."),
    "plot_type": "table",
    "pconfig": {"id": "metores_library_complexity_table", "namespace": "MeToReS"},
    "data": {r["sample"]: {
        "Genes detected (% median)": r["genes_detected_pct_of_median"],
        "ORFs assembled (% median)": r["orfs_assembled_pct_of_median"],
        "Reads mapped (% median)": r["reads_mapped_pct_of_median"],
        "Library type": r["libtype"],
        **({"Detected type": r["libtype_detected"] or "",
            "Type mismatch": "YES" if r["libtype_mismatch"] else "no"} if a.probe_dir else {}),
        "Flagged": "YES" if r["flagged"] else "no",
    } for _, r in df.iterrows()},
}
try:
    import yaml
    with open(a.out_mqc, "w") as fh:
        yaml.safe_dump(mqc, fh, sort_keys=False, allow_unicode=True)
except ImportError:  # JSON is valid YAML
    with open(a.out_mqc, "w") as fh:
        json.dump(mqc, fh, indent=2)

print(f"S2 library complexity: {len(flagged)} of {len(df)} samples flagged "
      f"(threshold {a.warn_fraction*100:.0f}% of run median)")
for _, r in flagged.iterrows():
    print(f"SAFEGUARD_WARNING: S2 low-complexity library {r['sample']}: {r['reason']}")
if df["libtype"].nunique() > 1:
    print("Library types in this run: " + ", ".join(
        f"{t} ({n})" for t, n in df["libtype"].value_counts().items()))
for _, r in df[df["libtype_mismatch"]].iterrows():
    print(f"SAFEGUARD_WARNING: S2 library type of {r['sample']}: declared {r['libtype']}, "
          f"Salmon detects {r['libtype_detected']}. Set its 'libtype' in the samplesheet "
          f"if the detection is right")
for m in mixed:
    print(f"SAFEGUARD_WARNING: S2 contrast {m['contrast']} mixes library types "
          f"{m['libtypes']}: library type is confounded with the comparison")
