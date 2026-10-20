#!/usr/bin/env python3
"""
assembly_qc.py — Gene catalogue Quality Control
Combines MetaQUAST + seqkit stats into a unified QC report
with structured JSON summary and MultiQC-compatible output.

NOTE ON TERMINOLOGY: this script runs on the POOLED NON-REDUNDANT GENE
CATALOGUE (predicted CDS from Prodigal, de-replicated with CD-HIT-EST),
not on raw contigs and not on a co-assembly. The pipeline assembles each
sample individually and pools the ORFs. Metrics are therefore reported per
GENE (ORF), not per contig. Per-sample contig statistics are produced
separately by seqkit stats in the Snakefile rule.

Outputs:
  - metaquast/                : Full MetaQUAST report directory
  - seqkit_stats.tsv          : Per-contig basic statistics
  - length_dist.pdf           : ORF length distribution plot
  - gc_dist.pdf               : GC content distribution plot
  - assembly_qc_summary.json  : Structured pass/warn/fail summary
  - assembly_qc_report.txt    : Human-readable report
  - assembly_qc_mqc.yaml      : MultiQC custom content section
"""

import argparse
import json
import os
import subprocess
import sys
import datetime

# ── optional plotting deps ────────────────────────────────────────────────────
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib/numpy not available — skipping distribution plots")

# ── optional BioPython for GC ─────────────────────────────────────────────────
try:
    from Bio import SeqIO
    from Bio.SeqUtils import gc_fraction
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("WARNING: Biopython not available — using built-in GC calculation")


# =============================================================================
# Argument parsing
# =============================================================================
def parse_args():
    p = argparse.ArgumentParser(
        description="Assembly QC: MetaQUAST + seqkit stats + custom metrics"
    )
    p.add_argument("-i", "--input",      required=True,
                   help="Input FASTA: pooled non-redundant gene catalogue (CDS)")
    p.add_argument("-o", "--output",     required=True,
                   help="Output directory for all QC results")
    p.add_argument("-r", "--references", default="",
                   help="Comma-separated reference genome FASTA paths "
                        "(optional). Leave empty for reference-free mode.")
    p.add_argument("-t", "--threads",    type=int, default=8,
                   help="Number of threads [default: 8]")
    p.add_argument("--min_contig",       type=int, default=200,
                   help="Minimum sequence length passed to MetaQUAST's "
                        "--min-contig [default: 200]")
    p.add_argument("--warn_n50",         type=int, default=500,
                   help="N50 below this value triggers WARN [default: 500]")
    p.add_argument("--warn_total_bp",    type=int, default=1_000_000,
                   help="Total catalogue size (bp) below this triggers WARN [default: 1000000]")
    p.add_argument("--warn_contigs",     type=int, default=100,
                   help="Gene count below this triggers WARN [default: 100]. "
                        "(Flag name kept as --warn_contigs for config compatibility.)")
    return p.parse_args()


# =============================================================================
# Utility
# =============================================================================
def run(cmd, step=""):
    tag = f"[{step}] " if step else ""
    print(f"\n{tag}[{datetime.datetime.now():%H:%M:%S}] Running:\n  {cmd}")
    rc = subprocess.call(cmd, shell=True)
    if rc != 0:
        print(f"ERROR: {tag}command returned exit code {rc}")
        raise subprocess.CalledProcessError(rc, cmd)


def gc_content_manual(seq):
    seq = seq.upper()
    gc  = seq.count("G") + seq.count("C")
    tot = len(seq) - seq.count("N")
    return (gc / tot * 100) if tot > 0 else 0.0


# =============================================================================
# FASTA parsing & metrics
# =============================================================================
def parse_fasta_metrics(fasta_path):
    """
    Returns (metrics dict, lengths list, gc_vals list).
    Computes N50, N90, L50, L90 and basic summary stats.
    """
    lengths = []
    gc_vals = []

    print(f"  Parsing FASTA: {fasta_path}")

    if HAS_BIOPYTHON:
        for rec in SeqIO.parse(fasta_path, "fasta"):
            seq = str(rec.seq)
            lengths.append(len(seq))
            gc_vals.append(gc_fraction(rec.seq) * 100)
    else:
        current_seq = []
        with open(fasta_path) as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if current_seq:
                        seq = "".join(current_seq)
                        lengths.append(len(seq))
                        gc_vals.append(gc_content_manual(seq))
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_seq:
                seq = "".join(current_seq)
                lengths.append(len(seq))
                gc_vals.append(gc_content_manual(seq))

    if not lengths:
        raise ValueError(f"No sequences found in {fasta_path}")

    lengths_sorted = sorted(lengths, reverse=True)
    total_bp       = sum(lengths_sorted)
    total_genes    = len(lengths_sorted)

    def nx_lx(sorted_lens, x):
        """Return (Nx length, Lx count) for given percentage x."""
        target = total_bp * x / 100
        cumsum = 0
        for i, ln in enumerate(sorted_lens, 1):
            cumsum += ln
            if cumsum >= target:
                return ln, i
        return sorted_lens[-1], total_genes

    n50, l50 = nx_lx(lengths_sorted, 50)
    n90, l90 = nx_lx(lengths_sorted, 90)

    sorted_for_median = sorted(lengths)
    median_len = sorted_for_median[total_genes // 2]

    metrics = {
        "total_genes"    : total_genes,
        "total_bp"       : total_bp,
        "n50"            : n50,
        "n90"            : n90,
        "l50"            : l50,
        "l90"            : l90,
        "longest_gene"   : lengths_sorted[0],
        "shortest_gene"  : lengths_sorted[-1],
        "mean_length"    : round(total_bp / total_genes, 1),
        "median_length"  : median_len,
        "mean_gc"        : round(sum(gc_vals) / len(gc_vals), 2),
    }

    return metrics, lengths, gc_vals


# =============================================================================
# Quality flags
# =============================================================================
def evaluate_metrics(metrics, args):
    flags   = {}
    overall = "PASS"

    checks = [
        ("n50",         metrics["n50"],         args.warn_n50,      "N50 (bp)"),
        ("total_bp",    metrics["total_bp"],    args.warn_total_bp, "Total catalogue (bp)"),
        ("total_genes", metrics["total_genes"], args.warn_contigs,  "Total genes"),
    ]

    for key, value, threshold, label in checks:
        status = "WARN" if value < threshold else "PASS"
        if status == "WARN":
            overall = "WARN"
        flags[key] = {
            "status"   : status,
            "value"    : value,
            "threshold": threshold,
            "label"    : label,
        }

    flags["overall"] = overall
    return flags


# =============================================================================
# Distribution plots
# =============================================================================
def plot_length_distribution(lengths, output_path):
    if not HAS_MATPLOTLIB:
        return
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].hist(lengths, bins=50, color="#4DBBD5", edgecolor="white", linewidth=0.5)
    axes[0].set_xlabel("ORF length (bp)", fontsize=12)
    axes[0].set_ylabel("Count", fontsize=12)
    axes[0].set_title("ORF Length Distribution", fontsize=13, fontweight="bold")
    median_val = sorted(lengths)[len(lengths) // 2]
    axes[0].axvline(median_val, color="#E64B35", linestyle="--",
                    linewidth=1.5, label=f"Median: {median_val:,} bp")
    axes[0].legend(fontsize=10)

    log_lengths = [max(0, __import__("math").log10(l + 1)) for l in lengths]
    axes[1].hist(log_lengths, bins=50, color="#00A087", edgecolor="white", linewidth=0.5)
    axes[1].set_xlabel("log₁₀(ORF length + 1)", fontsize=12)
    axes[1].set_ylabel("Count", fontsize=12)
    axes[1].set_title("ORF Length Distribution (log scale)",
                       fontsize=13, fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


def plot_gc_distribution(gc_vals, output_path):
    if not HAS_MATPLOTLIB:
        return
    mean_gc = sum(gc_vals) / len(gc_vals)
    fig, ax = plt.subplots(figsize=(8, 5))

    ax.hist(gc_vals, bins=50, color="#F39B7F", edgecolor="white", linewidth=0.5)
    ax.axvline(mean_gc, color="#E64B35", linestyle="--",
               linewidth=1.5, label=f"Mean GC: {mean_gc:.1f}%")
    ax.set_xlabel("GC content (%)", fontsize=12)
    ax.set_ylabel("Count",          fontsize=12)
    ax.set_title("GC Content Distribution", fontsize=13, fontweight="bold")
    ax.legend(fontsize=10)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {output_path}")


# =============================================================================
# Reports
# =============================================================================
def write_text_report(metrics, flags, output_path,
                       seqkit_tsv, metaquast_dir, references_used):
    now   = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    icon  = {"PASS": "✓ PASS", "WARN": "⚠ WARN", "FAIL": "✗ FAIL"}

    lines = [
        "=" * 62,
        "  Gene Catalogue QC Report",
        f"  Generated: {now}",
        "=" * 62,
        "",
        "── Basic Metrics ─────────────────────────────────────────────",
        f"  Total genes      : {metrics['total_genes']:>12,}",
        f"  Total bases      : {metrics['total_bp']:>12,}  bp",
        f"  Longest gene     : {metrics['longest_gene']:>12,}  bp",
        f"  Shortest gene    : {metrics['shortest_gene']:>12,}  bp",
        f"  Mean length      : {metrics['mean_length']:>12,}  bp",
        f"  Median length    : {metrics['median_length']:>12,}  bp",
        f"  Mean GC content  : {metrics['mean_gc']:>12.2f}  %",
        "",
        "── Catalogue Length Statistics ───────────────────────────────",
        f"  N50              : {metrics['n50']:>12,}  bp",
        f"  N90              : {metrics['n90']:>12,}  bp",
        f"  L50              : {metrics['l50']:>12,}  genes",
        f"  L90              : {metrics['l90']:>12,}  genes",
        "",
        "── Quality Flags ─────────────────────────────────────────────",
    ]

    for key, info in flags.items():
        if key == "overall":
            continue
        lines.append(
            f"  {icon[info['status']]}  {info['label']:<28}"
            f"  value={info['value']:>12,}"
            f"  threshold={info['threshold']:>10,}"
        )

    lines += [
        "",
        f"  Overall status : {icon[flags['overall']]}",
        "",
        "── Tools ─────────────────────────────────────────────────────",
        f"  seqkit stats : {seqkit_tsv}",
        f"  MetaQUAST    : {metaquast_dir}",
        f"  References   : {references_used or 'None (reference-free mode)'}",
        "",
        "=" * 62,
    ]

    with open(output_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print("\n" + "\n".join(lines))


def write_multiqc_custom(metrics, flags, output_path):
    """
    MultiQC custom content YAML — picked up automatically
    when MultiQC searches the assembly_qc output directory.
    """
    content = {
        "id"          : "assembly_qc",
        "section_name": "Gene Catalogue QC (MetaQUAST + seqkit)",
        "description" : "Quality metrics for the pooled non-redundant gene "
                        "catalogue (predicted CDS, CD-HIT-EST de-replicated).",
        "plot_type"   : "table",
        "pconfig"     : {"id": "assembly_qc_table", "title": "Gene Catalogue Metrics"},
        "data"        : {
            "Gene catalogue": {
                "Total genes"     : metrics["total_genes"],
                "Total bases (bp)": metrics["total_bp"],
                "N50 (bp)"        : metrics["n50"],
                "N90 (bp)"        : metrics["n90"],
                "L50"             : metrics["l50"],
                "L90"             : metrics["l90"],
                "Longest (bp)"    : metrics["longest_gene"],
                "Mean length (bp)": metrics["mean_length"],
                "Mean GC (%)"     : metrics["mean_gc"],
                "Status"          : flags["overall"],
            }
        }
    }

    try:
        import yaml
        with open(output_path, "w") as fh:
            yaml.dump(content, fh, default_flow_style=False)
    except ImportError:
        with open(output_path, "w") as fh:
            json.dump(content, fh, indent=2)

    print(f"  Saved MultiQC custom content: {output_path}")


# =============================================================================
# Main
# =============================================================================
def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)

    print("=" * 62)
    print("  Gene Catalogue QC")
    print("=" * 62)
    print(f"  Input    : {args.input}")
    print(f"  Output   : {args.output}")
    print(f"  Refs     : {args.references or 'None (reference-free)'}")
    print(f"  Threads  : {args.threads}")
    print("=" * 62)

    if not os.path.exists(args.input):
        sys.exit(f"ERROR: Input FASTA not found: {args.input}")

    # ── STEP 1: seqkit stats ──────────────────────────────────────────────────
    print("\n=== STEP 1: seqkit stats ===")
    seqkit_tsv = os.path.join(args.output, "seqkit_stats.tsv")
    run(f"seqkit stats -a -T {args.input} > {seqkit_tsv}", step="seqkit")

    # ── STEP 2: Custom FASTA metrics ─────────────────────────────────────────
    print("\n=== STEP 2: Computing gene catalogue metrics ===")
    metrics, lengths, gc_vals = parse_fasta_metrics(args.input)

    for label, key in [("Total genes  ", "total_genes"),
                        ("Total bases  ", "total_bp"),
                        ("N50          ", "n50"),
                        ("N90          ", "n90"),
                        ("L50          ", "l50"),
                        ("L90          ", "l90"),
                        ("Mean GC %    ", "mean_gc")]:
        val = metrics[key]
        print(f"  {label}: {val:,}")

    # ── STEP 3: Distribution plots ────────────────────────────────────────────
    print("\n=== STEP 3: Distribution plots ===")
    plot_length_distribution(lengths,  os.path.join(args.output, "length_dist.pdf"))
    plot_gc_distribution(gc_vals,      os.path.join(args.output, "gc_dist.pdf"))

    # ── STEP 4: MetaQUAST ────────────────────────────────────────────────────
    print("\n=== STEP 4: MetaQUAST ===")
    metaquast_dir   = os.path.join(args.output, "metaquast")
    ref_flag        = ""
    references_used = ""

    if args.references.strip():
        ref_paths = [r.strip() for r in args.references.split(",") if r.strip()]
        existing  = [r for r in ref_paths if os.path.exists(r)]
        missing   = [r for r in ref_paths if not os.path.exists(r)]

        if missing:
            print(f"  WARNING: {len(missing)} reference(s) not found — skipping:")
            for m in missing:
                print(f"    - {m}")
        if existing:
            ref_flag        = "-r " + ",".join(existing)
            references_used = ",".join(existing)
            print(f"  Using {len(existing)} reference(s)")
        else:
            print("  No valid references — running reference-free")
    else:
        print("  Running reference-free mode")

    run((
        f"metaquast.py {args.input} "
        f"{ref_flag} "
        f"-o {metaquast_dir} "
        f"--threads {args.threads} "
        f"--min-contig {args.min_contig} "
        f"--plots-format pdf "
        f"--silent"
    ), step="MetaQUAST")

    # ── STEP 5: Quality flags ─────────────────────────────────────────────────
    print("\n=== STEP 5: Quality flags ===")
    flags = evaluate_metrics(metrics, args)
    for key, info in flags.items():
        if key == "overall":
            continue
        icon = "✓" if info["status"] == "PASS" else "⚠"
        print(f"  {icon} {info['label']}: {info['value']:,}  [{info['status']}]")
    print(f"\n  Overall: {flags['overall']}")

    # ── STEP 6: Write all outputs ─────────────────────────────────────────────
    print("\n=== STEP 6: Writing reports ===")

    json_path = os.path.join(args.output, "assembly_qc_summary.json")
    with open(json_path, "w") as fh:
        json.dump({
            "timestamp"    : datetime.datetime.now().isoformat(),
            "input_fasta"  : args.input,
            "metrics"      : metrics,
            "quality_flags": flags,
            "tools"        : {
                "seqkit_stats" : seqkit_tsv,
                "metaquast_dir": metaquast_dir,
                "references"   : references_used or "reference-free",
            }
        }, fh, indent=2)
    print(f"  Saved: {json_path}")

    write_text_report(
        metrics, flags,
        os.path.join(args.output, "assembly_qc_report.txt"),
        seqkit_tsv, metaquast_dir, references_used
    )

    write_multiqc_custom(
        metrics, flags,
        os.path.join(args.output, "assembly_qc_mqc.yaml")
    )

    # ── Final summary ─────────────────────────────────────────────────────────
    print("\n" + "=" * 62)
    print("Gene Catalogue QC Complete!")
    print(f"  Status        : {flags['overall']}")
    print(f"  N50           : {metrics['n50']:,} bp")
    print(f"  Total genes   : {metrics['total_genes']:,}")
    print(f"  Total bases   : {metrics['total_bp']:,} bp")
    print(f"  Output dir    : {args.output}")
    print("=" * 62)

    sys.exit(0)   # Always exit 0 — warn but never block the pipeline


if __name__ == "__main__":
    main()
