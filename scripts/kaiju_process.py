#!/usr/bin/env python3
"""
kaiju_process.py — Kaiju Taxonomic Classification

Two classification modes:
  1. Gene catalog (ORFs) — protein-level, integrates with DEG analysis
  2. Per-sample reads    — community composition per sample (NEW)

Outputs (gene catalog mode):
  - kaiju_classifications.tsv
  - kaiju_taxonomy_full.tsv
  - kaiju_gene_taxonomy.tsv       ← merged into abundance matrix
  - {taxon}_ids.csv               ← contamination filter for downstream
  - {taxon}_ids_detailed.tsv
  - kaiju_summary_{rank}.tsv
  - krona_input.txt
  - krona_taxonomy.html
  - taxonomy_barplot_{rank}.pdf
  - kaiju_summary.json
  - kaiju_mqc.yaml

Outputs (per-sample mode):
  - {sra}_kaiju.out
  - {sra}_kaiju_summary.tsv       ← phylum-level per-sample table
"""

import argparse
import datetime
import json
import os
import subprocess
import sys

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: pandas is required. Install with: pip install pandas")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib not available — skipping bar chart")


# =============================================================================
# Argument parsing
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="Kaiju taxonomic classification — gene catalog and per-sample reads"
    )
    # ── Mode selection ─────────────────────────────────────────────────────
    mode = p.add_subparsers(dest="mode", help="Classification mode")
    mode.required = True

    # ── Gene catalog mode ──────────────────────────────────────────────────
    catalog = mode.add_parser(
        "catalog",
        help="Classify assembled ORFs (gene catalog)"
    )
    catalog.add_argument("-i",  "--input",        required=True,
                         help="Input ORF FASTA (all_longest_orfs_cds_rmdup_id.fasta)")
    catalog.add_argument("-o",  "--output",       required=True,
                         help="Output directory")
    catalog.add_argument("-db", "--kaiju_db",     required=True,
                         help="Kaiju .fmi database file")
    catalog.add_argument("-n",  "--nodes",        required=True,
                         help="NCBI nodes.dmp")
    catalog.add_argument("-nm", "--names",        required=True,
                         help="NCBI names.dmp")
    catalog.add_argument("-t",  "--threads",      type=int, default=8)
    catalog.add_argument("-f",  "--taxon_filter", default="Streptophyta",
                         help="Contamination taxon [default: Streptophyta]")
    catalog.add_argument("-r",  "--rank",         default="phylum",
                         choices=["phylum","class","order","family","genus","species"])
    catalog.add_argument("--top_n", type=int, default=20)

    # ── Per-sample mode ────────────────────────────────────────────────────
    sample = mode.add_parser(
        "sample",
        help="Classify per-sample cleaned reads"
    )
    sample.add_argument("-1",  "--fastq1",    required=True,
                        help="Forward cleaned FASTQ (rmrRNA)")
    sample.add_argument("-2",  "--fastq2",    required=True,
                        help="Reverse cleaned FASTQ (rmrRNA)")
    sample.add_argument("-s",  "--sample_id", required=True,
                        help="Sample/SRA identifier")
    sample.add_argument("-o",  "--output",    required=True,
                        help="Output directory")
    sample.add_argument("-db", "--kaiju_db",  required=True,
                        help="Kaiju .fmi database file")
    sample.add_argument("-n",  "--nodes",     required=True,
                        help="NCBI nodes.dmp")
    sample.add_argument("-nm", "--names",     required=True,
                        help="NCBI names.dmp")
    sample.add_argument("-t",  "--threads",   type=int, default=8)
    sample.add_argument("-r",  "--rank",      default="phylum",
                        choices=["phylum","class","order","family","genus","species"])

    return p.parse_args()


# =============================================================================
# Utilities
# =============================================================================
def run(cmd, step=""):
    tag = f"[{step}] " if step else ""
    print(f"\n{tag}[{datetime.datetime.now():%H:%M:%S}] Running:\n  {cmd}")
    rc = subprocess.call(cmd, shell=True)
    if rc != 0:
        raise subprocess.CalledProcessError(rc, cmd)


def log(msg):
    print(f"  {msg}")


# =============================================================================
# Shared helpers
# =============================================================================
def add_taxon_names(kaiju_out, nodes, names, output_dir, suffix=""):
    full_tax_file = os.path.join(
        output_dir,
        f"kaiju_taxonomy_full{suffix}.tsv"
    )
    run(
        f"kaiju-addTaxonNames "
        f"-t {nodes} -n {names} "
        f"-i {kaiju_out} -o {full_tax_file} "
        f"-r superkingdom,phylum,class,order,family,genus,species ",
        step="kaiju-addTaxonNames"
    )
    return full_tax_file


def parse_taxonomy_table(full_tax_file):
    log(f"Parsing: {full_tax_file}")
    rows = []
    with open(full_tax_file) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            status   = parts[0].strip()
            gene_id = parts[1].strip()
            taxon_id  = parts[2].strip() if len(parts) > 2 else ""
            score    = parts[3].strip() if len(parts) > 3 else ""
            lineage  = parts[7].strip() if len(parts) >= 8 else ""
            lin_parts = [x.strip() for x in lineage.split(";")]
            while len(lin_parts) < 7:
                lin_parts.append("")

            def clean(v):
                return v if v else "unclassified"

            rows.append({
                "gene_id"      : gene_id,
                "status"       : status,
                "taxon_id"     : taxon_id,
                "score"        : score,
                "superkingdom" : clean(lin_parts[0]),
                "phylum"       : clean(lin_parts[1]),
                "class"        : clean(lin_parts[2]),
                "order"        : clean(lin_parts[3]),
                "family"       : clean(lin_parts[4]),
                "genus"        : clean(lin_parts[5]),
                "species"      : clean(lin_parts[6]),
            })

    df = pd.DataFrame(rows)
    classified   = (df["status"] == "C").sum()
    unclassified = (df["status"] == "U").sum()
    log(f"Total entries       : {len(df):,}")
    log(f"  Classified (C)    : {classified:,}")
    log(f"  Unclassified (U)  : {unclassified:,}")
    log(f"  Classification %  : {classified/len(df)*100:.1f}%")
    return df


# =============================================================================
# CATALOG mode functions
# =============================================================================
def run_kaiju_catalog(input_fasta, kaiju_db, nodes, output_dir, threads):
    out_file = os.path.join(output_dir, "kaiju_classifications.tsv")
    run(
        f"kaiju -z {threads} -t {nodes} -f {kaiju_db} "
        f"-i {input_fasta} -o {out_file} -v -x",
        step="kaiju-catalog"
    )
    return out_file


def save_gene_taxonomy(df, output_dir):
    out_file = os.path.join(output_dir, "kaiju_gene_taxonomy.tsv")
    df[["gene_id","status","superkingdom","phylum","class",
        "order","family","genus","species"]].to_csv(
        out_file, sep="\t", index=False
    )
    log(f"Saved per-gene taxonomy : {out_file}")
    return out_file


def generate_contamination_filter(df, taxon_filter, output_dir):
    rank_cols = ["superkingdom","phylum","class",
                 "order","family","genus","species"]
    mask = df[rank_cols].apply(
        lambda col: col.str.contains(taxon_filter, case=False, na=False)
    ).any(axis=1)

    contam_df  = df[mask]
    contam_ids = contam_df["gene_id"].tolist()

    log(f"Taxon filter        : '{taxon_filter}'")
    log(f"Contamination ORFs  : {len(contam_ids):,} "
        f"({len(contam_ids)/len(df)*100:.2f}%)")

    filter_file = os.path.join(output_dir,
                                f"{taxon_filter.lower()}_ids.csv")
    pd.DataFrame({"gene_id": contam_ids}).to_csv(
        filter_file, index=False, header=False
    )
    log(f"Saved contamination IDs : {filter_file}")

    detail_file = os.path.join(output_dir,
                                f"{taxon_filter.lower()}_ids_detailed.tsv")
    contam_df[["gene_id","superkingdom","phylum","class",
               "order","family","genus","species"]].to_csv(
        detail_file, sep="\t", index=False
    )
    return filter_file, contam_ids


def generate_summary(df, rank, output_dir, suffix=""):
    classified = df[df["status"] == "C"]
    summary = (
        classified.groupby(rank).size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
    )
    summary["percentage"] = (
        summary["count"] / len(classified) * 100
    ).round(2)
    out_file = os.path.join(output_dir, f"kaiju_summary_{rank}{suffix}.tsv")
    summary.to_csv(out_file, sep="\t", index=False)
    log(f"Saved {rank} summary : {out_file}")
    log(f"  Unique {rank}s     : {len(summary):,}")
    return summary


def generate_krona(df, output_dir):
    rank_cols  = ["superkingdom","phylum","class",
                  "order","family","genus","species"]
    krona_in   = os.path.join(output_dir, "krona_input.txt")
    krona_html = os.path.join(output_dir, "krona_taxonomy.html")

    classified = df[df["status"] == "C"]
    grouped    = (
        classified.groupby(rank_cols).size()
        .reset_index(name="count")
    )
    with open(krona_in, "w") as fh:
        for _, row in grouped.iterrows():
            lineage = "\t".join(str(row[r]) for r in rank_cols)
            fh.write(f"{row['count']}\t{lineage}\n")
    log(f"Saved Krona input : {krona_in}")

    rc = subprocess.call(
        f"ktImportText {krona_in} -o {krona_html}",
        shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )
    if rc == 0:
        log(f"Saved Krona HTML  : {krona_html}")
    else:
        log("WARNING: ktImportText not found — Krona HTML skipped")
    return krona_in


def generate_barplot(summary, rank, top_n, taxon_filter, output_dir,
                     title_suffix=""):
    if not HAS_MATPLOTLIB:
        return
    plot_data = summary.head(top_n).copy()
    colors = []
    for taxon in plot_data[rank]:
        if taxon_filter.lower() in str(taxon).lower():
            colors.append("#E64B35")
        elif taxon == "unclassified":
            colors.append("#B0B0B0")
        else:
            colors.append("#4DBBD5")

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.barh(plot_data[rank][::-1], plot_data["percentage"][::-1],
            color=colors[::-1], edgecolor="white", linewidth=0.5)
    for i, (_, row) in enumerate(plot_data[::-1].iterrows()):
        ax.text(row["percentage"] + 0.2, i,
                f"{row['count']:,}", va="center", ha="left",
                fontsize=8, color="grey")
    ax.set_xlabel("% of classified ORFs", fontsize=12)
    ax.set_title(
        f"Taxonomic Classification — {rank.capitalize()} level{title_suffix}\n"
        f"(Top {top_n} | red = {taxon_filter})",
        fontsize=13, fontweight="bold"
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    out_file = os.path.join(output_dir, f"taxonomy_barplot_{rank}.pdf")
    plt.savefig(out_file, dpi=150, bbox_inches="tight")
    plt.close()
    log(f"Saved bar chart : {out_file}")


def write_multiqc_custom(df, summary, rank, taxon_filter,
                          contam_ids, output_dir):
    classified_n = int((df["status"] == "C").sum())
    total_n      = len(df)
    top5 = summary.head(5).to_dict("records")
    top5_str = "; ".join(
        f"{r[rank]} ({r['percentage']}%)" for r in top5
    )
    content = {
        "id"          : "kaiju_taxonomy",
        "section_name": "Kaiju Taxonomic Classification",
        "description" : f"Protein-level classification. Filter: {taxon_filter}.",
        "plot_type"   : "table",
        "pconfig"     : {"id": "kaiju_table", "title": "Kaiju Summary"},
        "data"        : {
            "ORF Classification": {
                "Total ORFs"              : total_n,
                "Classified"              : classified_n,
                "Unclassified"            : total_n - classified_n,
                "Classification rate (%)" : round(classified_n/total_n*100, 1),
                f"Unique {rank}s"         : len(summary),
                f"Top 5 {rank}s"          : top5_str,
                f"{taxon_filter} ORFs"    : len(contam_ids),
                f"{taxon_filter} (%)"     : round(len(contam_ids)/total_n*100, 2),
            }
        }
    }
    mqc_file = os.path.join(output_dir, "kaiju_mqc.yaml")
    try:
        import yaml
        with open(mqc_file, "w") as fh:
            yaml.dump(content, fh, default_flow_style=False)
    except ImportError:
        with open(mqc_file, "w") as fh:
            json.dump(content, fh, indent=2)
    log(f"Saved MultiQC content : {mqc_file}")


def write_json_summary(df, summary, rank, taxon_filter,
                        contam_ids, output_dir):
    classified_n = int((df["status"] == "C").sum())
    total_n      = len(df)
    out = {
        "timestamp"           : datetime.datetime.now().isoformat(),
        "total_orfs"          : total_n,
        "classified"          : classified_n,
        "unclassified"        : total_n - classified_n,
        "classification_rate" : round(classified_n / total_n * 100, 2),
        f"unique_{rank}s"     : len(summary),
        f"top10_{rank}s"      : summary.head(10).to_dict("records"),
        "contamination_filter": taxon_filter,
        "contamination_orfs"  : len(contam_ids),
        "contamination_pct"   : round(len(contam_ids) / total_n * 100, 2),
    }
    json_file = os.path.join(output_dir, "kaiju_summary.json")
    with open(json_file, "w") as fh:
        json.dump(out, fh, indent=2)
    log(f"Saved JSON summary : {json_file}")
    return json_file


# =============================================================================
# PER-SAMPLE mode functions
# =============================================================================
def run_kaiju_sample(fastq1, fastq2, sample_id,
                      kaiju_db, nodes, output_dir, threads):
    """
    Classify paired-end reads directly — gives community composition
    per sample rather than per gene.
    -i / -j for paired-end input.
    """
    out_file = os.path.join(output_dir, f"{sample_id}_kaiju.out")
    run(
        f"kaiju -z {threads} -t {nodes} -f {kaiju_db} "
        f"-i {fastq1} -j {fastq2} "
        f"-o {out_file} -v -x",
        step=f"kaiju-sample-{sample_id}"
    )
    return out_file


def generate_sample_summary(kaiju_out, nodes, names,
                              sample_id, rank, output_dir):
    """
    Use kaiju2table to produce a phylum-level count table per sample.
    Output compatible with krona and MultiQC custom content.
    """
    summary_file = os.path.join(output_dir, f"{sample_id}_kaiju_summary.tsv")

    run(
        f"kaiju2table "
        f"-t {nodes} -n {names} "
        f"-r {rank} "
        f"-o {summary_file} "
        f"{kaiju_out}",
        step=f"kaiju2table-{sample_id}"
    )
    log(f"Saved sample summary : {summary_file}")

    # Print classification stats
    total      = sum(1 for _ in open(kaiju_out))
    classified = sum(1 for l in open(kaiju_out) if l.startswith("C"))
    log(f"  Total reads      : {total:,}")
    log(f"  Classified       : {classified:,} ({classified/total*100:.1f}%)")

    return summary_file


# =============================================================================
# Main
# =============================================================================
def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)

    # ── Gene catalog mode ──────────────────────────────────────────────────
    if args.mode == "catalog":
        print("=" * 62)
        print("  Kaiju — Gene Catalog Classification")
        print("=" * 62)
        print(f"  Input ORFs   : {args.input}")
        print(f"  Output dir   : {args.output}")
        print(f"  Taxon filter : {args.taxon_filter}")
        print(f"  Summary rank : {args.rank}")
        print("=" * 62)

        for f in [args.input, args.kaiju_db, args.nodes, args.names]:
            if not os.path.exists(f):
                sys.exit(f"ERROR: File not found: {f}")

        print("\n=== STEP 1: Running Kaiju ===")
        kaiju_out = run_kaiju_catalog(
            args.input, args.kaiju_db,
            args.nodes, args.output, args.threads
        )

        print("\n=== STEP 2: Adding full lineage ===")
        full_tax = add_taxon_names(
            kaiju_out, args.nodes, args.names, args.output
        )

        print("\n=== STEP 3: Parsing taxonomy table ===")
        df = parse_taxonomy_table(full_tax)

        print("\n=== STEP 4: Per-gene taxonomy table ===")
        save_gene_taxonomy(df, args.output)

        print(f"\n=== STEP 5: Contamination filter ({args.taxon_filter}) ===")
        _, contam_ids = generate_contamination_filter(
            df, args.taxon_filter, args.output
        )

        print(f"\n=== STEP 6: Summary at {args.rank} level ===")
        summary = generate_summary(df, args.rank, args.output)

        print("\n=== STEP 7: Krona chart ===")
        generate_krona(df, args.output)

        print("\n=== STEP 8: Bar chart ===")
        generate_barplot(summary, args.rank, args.top_n,
                         args.taxon_filter, args.output)

        print("\n=== STEP 9: MultiQC content ===")
        write_multiqc_custom(df, summary, args.rank,
                              args.taxon_filter, contam_ids, args.output)

        print("\n=== STEP 10: JSON summary ===")
        write_json_summary(df, summary, args.rank,
                            args.taxon_filter, contam_ids, args.output)

        classified_n = (df["status"] == "C").sum()
        print("\n" + "=" * 62)
        print("  Kaiju Catalog Classification Complete!")
        print("=" * 62)
        print(f"  Total ORFs      : {len(df):,}")
        print(f"  Classified      : {classified_n:,} ({classified_n/len(df)*100:.1f}%)")
        print(f"  {args.taxon_filter:<16} : {len(contam_ids):,} ORFs flagged")
        print("=" * 62)

    # ── Per-sample mode ────────────────────────────────────────────────────
    elif args.mode == "sample":
        print("=" * 62)
        print(f"  Kaiju — Per-Sample Read Classification: {args.sample_id}")
        print("=" * 62)
        print(f"  FASTQ R1   : {args.fastq1}")
        print(f"  FASTQ R2   : {args.fastq2}")
        print(f"  Output dir : {args.output}")
        print(f"  Rank       : {args.rank}")
        print("=" * 62)

        for f in [args.fastq1, args.fastq2, args.kaiju_db,
                  args.nodes, args.names]:
            if not os.path.exists(f):
                sys.exit(f"ERROR: File not found: {f}")

        print(f"\n=== Running Kaiju on sample {args.sample_id} ===")
        kaiju_out = run_kaiju_sample(
            args.fastq1, args.fastq2, args.sample_id,
            args.kaiju_db, args.nodes, args.output, args.threads
        )

        print(f"\n=== Generating {args.rank}-level summary ===")
        generate_sample_summary(
            kaiju_out, args.nodes, args.names,
            args.sample_id, args.rank, args.output
        )

        print("\n" + "=" * 62)
        print(f"  Sample {args.sample_id} classification complete!")
        print("=" * 62)

    sys.exit(0)


if __name__ == "__main__":
    main()
