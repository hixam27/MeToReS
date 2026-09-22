#!/usr/bin/env python3
"""
merge_taxonomy.py — Merge Kaiju taxonomy into filtered abundance matrix

Joins kaiju_gene_taxonomy.tsv with the filtered TPM abundance matrix,
adding full lineage columns (superkingdom → species) per gene.

Inputs:
  - Filtered abundance matrix CSV  (GeneID + sample TPM columns)
  - Kaiju gene taxonomy TSV        (kaiju_gene_taxonomy.tsv)

Outputs:
  - transcript_abundance_quantification_table_filter_taxonomy.csv
      Full matrix with taxonomy columns appended
  - taxonomy_merge_summary.json
      Merge statistics (how many genes got annotations)
"""

import argparse
import datetime
import json
import os
import sys

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: pandas is required. Install with: pip install pandas")


# =============================================================================
# Argument parsing
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="Merge Kaiju per-gene taxonomy into filtered abundance matrix"
    )
    p.add_argument("-m", "--matrix",   required=True,
                   help="Filtered abundance matrix CSV (GeneID as first column)")
    p.add_argument("-t", "--taxonomy", required=True,
                   help="Kaiju per-gene taxonomy TSV (kaiju_gene_taxonomy.tsv)")
    p.add_argument("-o", "--output",   required=True,
                   help="Output directory")
    return p.parse_args()


# =============================================================================
# Main
# =============================================================================
def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)

    print("=" * 62)
    print("  Taxonomy Merge — Kaiju + Abundance Matrix")
    print("=" * 62)
    print(f"  Matrix   : {args.matrix}")
    print(f"  Taxonomy : {args.taxonomy}")
    print(f"  Output   : {args.output}")
    print("=" * 62)

    # ── Validate inputs ───────────────────────────────────────────────────────
    for f in [args.matrix, args.taxonomy]:
        if not os.path.exists(f):
            sys.exit(f"ERROR: File not found: {f}")

    # ── Load filtered abundance matrix ────────────────────────────────────────
    print("\n=== STEP 1: Loading abundance matrix ===")
    matrix = pd.read_csv(args.matrix, header=0)

    if "GeneID" not in matrix.columns:
        sys.exit("ERROR: Matrix must have a 'GeneID' column as first column.")

    print(f"  Genes   : {len(matrix):,}")
    print(f"  Samples : {matrix.shape[1] - 1}")
    print(f"  Columns : {list(matrix.columns[:4])} ...")

    # ── Load Kaiju taxonomy table ─────────────────────────────────────────────
    print("\n=== STEP 2: Loading Kaiju taxonomy table ===")
    taxonomy = pd.read_csv(args.taxonomy, sep="\t", header=0)

    required_tax_cols = [
        "gene_id", "status",
        "superkingdom", "phylum", "class",
        "order", "family", "genus", "species"
    ]
    missing = [c for c in required_tax_cols if c not in taxonomy.columns]
    if missing:
        sys.exit(f"ERROR: Taxonomy TSV missing columns: {missing}")

    print(f"  Total ORFs in taxonomy : {len(taxonomy):,}")
    print(f"  Classified (C)         : {(taxonomy['status'] == 'C').sum():,}")
    print(f"  Unclassified (U)       : {(taxonomy['status'] == 'U').sum():,}")

    # Keep only the lineage columns we want to append
    tax_cols = [
        "gene_id",
        "superkingdom", "phylum", "class",
        "order", "family", "genus", "species",
        "status"
    ]
    taxonomy_slim = taxonomy[tax_cols].copy()

    # Rename status to kaiju_status to avoid ambiguity
    taxonomy_slim = taxonomy_slim.rename(columns={"status": "kaiju_status"})

    # ── Merge ─────────────────────────────────────────────────────────────────
    print("\n=== STEP 3: Merging ===")
    
    # Coerce both join keys to string to avoid int64/object mismatch
    matrix["GeneID"]        = matrix["GeneID"].astype(str)
    taxonomy_slim["gene_id"] = taxonomy_slim["gene_id"].astype(str)

    merged = matrix.merge(
        taxonomy_slim,
        left_on  = "GeneID",
        right_on = "gene_id",
        how      = "left"          # keep ALL genes from matrix
    ).drop(columns=["gene_id"])    # drop redundant join key

    # Fill unmatched genes (not in Kaiju output) with 'unclassified'
    tax_annotation_cols = [
        "superkingdom", "phylum", "class",
        "order", "family", "genus", "species", "kaiju_status"
    ]
    for col in tax_annotation_cols:
        merged[col] = merged[col].fillna("unclassified")

    # ── Reorder columns: GeneID | taxonomy | samples ──────────────────────────
    sample_cols = [c for c in matrix.columns if c != "GeneID"]

    col_order = (
        ["GeneID"]
        + ["superkingdom", "phylum", "class",
           "order", "family", "genus", "species", "kaiju_status"]
        + sample_cols
    )
    merged = merged[col_order]

    print(f"  Genes in matrix          : {len(matrix):,}")
    print(f"  Genes after merge        : {len(merged):,}")

    annotated = (merged["kaiju_status"] == "C").sum()
    unannotated = (merged["kaiju_status"] != "C").sum()
    print(f"  Classified (C)           : {annotated:,} "
          f"({annotated/len(merged)*100:.1f}%)")
    print(f"  Unclassified / no match  : {unannotated:,} "
          f"({unannotated/len(merged)*100:.1f}%)")

    # ── Breakdown by superkingdom ─────────────────────────────────────────────
    print("\n  Superkingdom breakdown:")
    sk_counts = (
        merged.groupby("superkingdom")
        .size()
        .sort_values(ascending=False)
    )
    for sk, n in sk_counts.items():
        print(f"    {sk:<30} : {n:>8,}  ({n/len(merged)*100:.1f}%)")

    # ── Save merged matrix ────────────────────────────────────────────────────
    print("\n=== STEP 4: Saving outputs ===")

    out_matrix = os.path.join(
        args.output,
        "transcript_abundance_quantification_table_filter_taxonomy.csv"
    )
    merged.to_csv(out_matrix, index=False)
    print(f"  Saved annotated matrix : {out_matrix}")
    print(f"  Dimensions             : {merged.shape[0]:,} genes "
          f"x {merged.shape[1]} columns")

    # ── JSON summary ──────────────────────────────────────────────────────────
    summary = {
        "timestamp"          : datetime.datetime.now().isoformat(),
        "input_matrix"       : args.matrix,
        "input_taxonomy"     : args.taxonomy,
        "total_genes"        : len(merged),
        "classified"         : int(annotated),
        "unclassified"       : int(unannotated),
        "classified_pct"     : round(annotated / len(merged) * 100, 2),
        "superkingdom_counts": {
            str(k): int(v) for k, v in sk_counts.items()
        },
        "output_matrix"      : out_matrix,
    }

    json_out = os.path.join(args.output, "taxonomy_merge_summary.json")
    with open(json_out, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"  Saved JSON summary     : {json_out}")

    # ── Final summary ─────────────────────────────────────────────────────────
    print("\n" + "=" * 62)
    print("  Taxonomy Merge Complete!")
    print("=" * 62)
    print(f"  Total genes   : {len(merged):,}")
    print(f"  Classified    : {annotated:,} ({annotated/len(merged)*100:.1f}%)")
    print(f"  Output        : {out_matrix}")
    print("=" * 62)
    sys.exit(0)


if __name__ == "__main__":
    main()
