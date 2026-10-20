#!/usr/bin/env python3
"""
6.DEG_analysis.py — Build the gene-by-sample TPM abundance matrix

Reads the per-sample Salmon quant.sf files and merges their TPM columns into a
single matrix keyed by gene ID (the pooled non-redundant gene catalogue IDs).

SCOPE — read this before extending the script:
    This script ONLY builds the abundance matrix. It is invoked by the
    Snakefile's `filter_genes` rule as:

        python 6.DEG_analysis.py -i <quant_dir> -s <samplesheet> \
               -o <out_dir> --group_col <col> --matrix_only

    Everything downstream is handled by dedicated rules:
      - gene filtering            -> filter_genes.R   (rule filter_genes)
      - differential expression   -> DEG_analysis.r   (rule DEG_analysis)
      - per-contrast sequences    -> seqkit           (rule DEG_analysis)


Output:
    transcript_abundance_quantification_table.csv  (GeneID + one TPM column per sample)
"""

import argparse
import os
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build the TPM abundance matrix from Salmon quant files"
    )
    parser.add_argument("-i", "--input_dir", required=True,
                        help="Directory containing the per-sample <sra>_quant/quant.sf files")
    parser.add_argument("-s", "--samplesheet", required=True,
                        help="Sample sheet CSV with at least: sample, sra, <group_col>")
    parser.add_argument("-o", "--output", required=True,
                        help="Output directory for the abundance matrix")
    parser.add_argument("--group_col", default="group",
                        help="Grouping column name in the samplesheet [default: group]")
    parser.add_argument("--matrix_only", action="store_true",
                        help="Accepted for compatibility with the Snakefile call. "
                             "Matrix generation is now this script's only mode, so "
                             "this flag is a no-op.")
    return parser.parse_args()


def create_abundance_matrix(input_dir, samplesheet, output_dir, group_col="group"):
    """Merge the per-sample Salmon TPM columns into one gene-by-sample matrix."""
    print("=" * 60)
    print("Creating abundance matrix from quant files...")
    print("=" * 60)

    samples_df = pd.read_csv(samplesheet)
    print(f"Loaded {len(samples_df)} samples from {samplesheet}")
    if group_col not in samples_df.columns:
        raise ValueError(
            f"Column '{group_col}' not in samplesheet. Available: {list(samples_df.columns)}"
        )
    print(f"Groups: {samples_df[group_col].unique()}")

    abundance_dict = {}
    missing_files = []

    for _, row in samples_df.iterrows():
        sra = row["sra"]
        sample_name = row["sample"]
        quant_file = os.path.join(input_dir, f"{sra}_quant", "quant.sf")

        if os.path.exists(quant_file):
            df = pd.read_csv(quant_file, sep="\t")
            abundance_dict[sample_name] = df.set_index("Name")["TPM"]
            print(f"  OK  Loaded {sample_name} ({sra}): {len(df)} genes")
        else:
            missing_files.append(quant_file)
            print(f"  --  Warning: {quant_file} not found")

    if missing_files:
        print(f"\nWARNING: {len(missing_files)} quant file(s) missing!")
        for f in missing_files:
            print(f"  - {f}")

    if not abundance_dict:
        raise ValueError("No quant files found! Check the input directory.")

    abundance_matrix = pd.DataFrame(abundance_dict)
    abundance_matrix.index.name = "GeneID"
    abundance_matrix = abundance_matrix.reset_index()

    if "GeneID" not in abundance_matrix.columns:
        raise ValueError(f"GeneID column not created! Columns: {list(abundance_matrix.columns)}")

    matrix_file = os.path.join(output_dir, "transcript_abundance_quantification_table.csv")
    abundance_matrix.to_csv(matrix_file, index=False)

    print("\n" + "=" * 60)
    print("Abundance Matrix Summary:")
    print("=" * 60)
    print(f"Matrix file: {matrix_file}")
    print(f"Dimensions: {abundance_matrix.shape[0]} genes x {abundance_matrix.shape[1]-1} samples")
    print(f"Samples: {list(abundance_matrix.columns[1:])}")
    print(f"Total non-zero entries: {(abundance_matrix.iloc[:, 1:] > 0).sum().sum()}")
    print("=" * 60 + "\n")

    return matrix_file


def main():
    args = parse_args()

    os.makedirs(args.output, exist_ok=True)
    print(f"\nOutput directory: {args.output}\n")

    matrix_file = create_abundance_matrix(
        args.input_dir, args.samplesheet, args.output, args.group_col
    )

    print(f"Matrix generated: {matrix_file}")
    print("Next step in the pipeline: filter_genes.R (rule filter_genes).")


if __name__ == "__main__":
    main()
