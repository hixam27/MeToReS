#!/usr/bin/env python3
"""test_gene_selection.py — checks for deseq2.gene_set (select_coding_genes.py)
and its reporting in safeguard S2 and the pipeline summary.

Small synthetic inputs where the answer is known by construction:
  coding  keeps exactly the filtered genes with an eggNOG-mapper hit, lines
          unchanged; counts shares overall and per sample are right
  all     output is byte-identical to the filtered matrix; same shares reported
  none    no annotated gene: "coding" stops with an error, "all" still runs
  S2      the per-sample share appears in library_complexity.json and in
          SAFEGUARDS_SUMMARY.md, and does not flag anything
  summary PIPELINE_SUMMARY.md reports the share; with "all" also the share of
          DEGs per contrast with protein-level evidence

Usage (from the repository root):  python tests/test_gene_selection.py
Exits with status 1 if any check fails. Needs pandas only.
"""
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")
FAILS = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail and not ok else ""))
    if not ok:
        FAILS.append(name)


def run(args, expect_ok=True):
    p = subprocess.run([sys.executable] + args, capture_output=True, text=True)
    if expect_ok and p.returncode != 0:
        sys.stderr.write(p.stdout + p.stderr)
        raise SystemExit(f"command failed: {' '.join(args)}")
    return p


def build(w, annotated):
    """Two samples, six filtered genes (g1-g6) plus one gene outside the filter (g7)."""
    os.makedirs(os.path.join(w, "deg"), exist_ok=True)
    with open(os.path.join(w, "samplesheet.csv"), "w") as fh:
        fh.write("sample,sra,group\nS1,R1,A\nS2,R2,B\n")
    rows = ['"GeneID","S1","S2"\n'] + [f'"g{i}",{i}.5,{2 * i}.25\n' for i in range(1, 7)]
    with open(os.path.join(w, "deg", "filter.csv"), "w") as fh:
        fh.writelines(rows)
    counts = {"R1": [10, 20, 30, 40, 50, 60, 1000], "R2": [5, 5, 5, 5, 5, 75, 1000]}
    for sra, c in counts.items():
        d = os.path.join(w, "quant", f"{sra}_quant")
        os.makedirs(os.path.join(d, "aux_info"), exist_ok=True)
        with open(os.path.join(d, "quant.sf"), "w") as fh:
            fh.write("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
            fh.writelines(f"g{i + 1}\t300\t250\t1.0\t{v}\n" for i, v in enumerate(c))
        with open(os.path.join(d, "aux_info", "meta_info.json"), "w") as fh:
            json.dump({"num_mapped": sum(c), "percent_mapped": 50.0}, fh)
    with open(os.path.join(w, "ann.emapper.annotations"), "w") as fh:
        fh.write("## emapper\n#query\tseed_ortholog\tKEGG_ko\n")
        fh.writelines(f"{g}\tseed\t-\n" for g in annotated)
        fh.write("## done\n")
    return rows


def select(w, mode, expect_ok=True):
    return run([os.path.join(SCRIPTS, "select_coding_genes.py"),
                "--matrix", os.path.join(w, "deg", "filter.csv"),
                "--annotations", os.path.join(w, "ann.emapper.annotations"),
                "--samplesheet", os.path.join(w, "samplesheet.csv"),
                "--quant_dir", os.path.join(w, "quant"), "--gene_set", mode,
                "--out_matrix", os.path.join(w, "deg", f"coding_{mode}.csv"),
                "--out_json", os.path.join(w, "deg", f"gene_set_{mode}.json")], expect_ok)


def main():
    with tempfile.TemporaryDirectory() as w:
        # annotated: g2, g4, g6 in the filter; g7 annotated but not filtered
        rows = build(w, ["g2", "g4", "g6", "g7"])
        select(w, "coding")
        out = open(os.path.join(w, "deg", "coding_coding.csv")).readlines()
        check("coding keeps the annotated filtered genes, lines unchanged",
              out == [rows[0], rows[2], rows[4], rows[6]], f"got {out}")
        j = json.load(open(os.path.join(w, "deg", "gene_set_coding.json")))
        # S1: (20+40+60)/210 = 57.1 %; S2: (5+5+75)/100 = 85 %; overall 205/310 = 66.1 %
        check("coding JSON: genes", (j["n_filtered"], j["n_with_protein_evidence"], j["n_genes_selected"]) == (6, 3, 3), str(j))
        check("coding JSON: per-sample shares",
              (j["per_sample"]["S1"]["pct_counts_with_protein_evidence"],
               j["per_sample"]["S2"]["pct_counts_with_protein_evidence"]) == (57.1, 85.0), str(j["per_sample"]))
        check("coding JSON: overall share and range",
              (j["pct_counts_with_protein_evidence"], j["pct_counts_min"], j["pct_counts_max"]) == (66.1, 57.1, 85.0), str(j))
        check("genes outside the filter are ignored", j["counts_filtered"] == 310.0, str(j["counts_filtered"]))

        select(w, "all")
        same = open(os.path.join(w, "deg", "filter.csv"), "rb").read() == \
            open(os.path.join(w, "deg", "coding_all.csv"), "rb").read()
        check("all: output byte-identical to the filtered matrix", same)
        ja = json.load(open(os.path.join(w, "deg", "gene_set_all.json")))
        check("all: same shares reported, every gene selected",
              ja["pct_counts_with_protein_evidence"] == 66.1 and ja["n_genes_selected"] == 6, str(ja))

        # S2 reports the per-sample share and flags nothing because of it
        for sra in ("R1", "R2"):
            os.makedirs(os.path.join(w, "assembly", sra), exist_ok=True)
            with open(os.path.join(w, "assembly", sra, "orfs.fasta"), "w") as fh:
                fh.write(">o1\nATG\n>o2\nATG\n")
        sg = os.path.join(w, "safeguards"); os.makedirs(sg, exist_ok=True)
        run([os.path.join(SCRIPTS, "safeguard_library_complexity.py"),
             "--samplesheet", os.path.join(w, "samplesheet.csv"),
             "--matrix", os.path.join(w, "deg", "coding_coding.csv"),
             "--gene_set_json", os.path.join(w, "deg", "gene_set_coding.json"),
             "--assembly_dir", os.path.join(w, "assembly"), "--quant_dir", os.path.join(w, "quant"),
             "--default_libtype", "IU", "--out_tsv", os.path.join(sg, "lc.tsv"),
             "--out_json", os.path.join(sg, "lc.json"), "--out_mqc", os.path.join(sg, "lc_mqc.yaml")])
        s2 = json.load(open(os.path.join(sg, "lc.json")))
        pe = s2.get("protein_evidence", {})
        check("S2 JSON carries the per-sample share",
              pe.get("per_sample") == {"S1": 57.1, "S2": 85.0} and pe.get("gene_set") == "coding", str(pe))
        check("S2 does not flag on the share", s2["n_flagged"] == 0, str(s2["flagged"]))
        run([os.path.join(SCRIPTS, "safeguards_summary.py"), "--jsons", os.path.join(sg, "lc.json"),
             "--out", os.path.join(sg, "SAFEGUARDS_SUMMARY.md")])
        md = open(os.path.join(sg, "SAFEGUARDS_SUMMARY.md")).read()
        check("SAFEGUARDS_SUMMARY.md shows the share",
              "protein-level evidence" in md and "| S1 | 57.1 |" in md and "66.1% of the filtered-gene counts" in md)

        # Pipeline summary, "all" mode: share of DEGs with protein-level evidence
        res = os.path.join(w, "results"); dr = os.path.join(res, "deg", "group_deg_results")
        os.makedirs(os.path.join(dr, "B_vs_A"), exist_ok=True)
        os.makedirs(os.path.join(res, "emapper_group"), exist_ok=True)
        os.replace(os.path.join(w, "ann.emapper.annotations"),
                   os.path.join(res, "emapper_group", "all_genes_emapper.emapper.annotations"))
        os.replace(os.path.join(w, "deg", "gene_set_all.json"), os.path.join(dr, "gene_set_summary.json"))
        with open(os.path.join(dr, "all_contrasts_summary.csv"), "w") as fh:
            fh.write("Contrast,Total_DEGs,Upregulated,Downregulated,Genes_tested\nB_vs_A,4,3,1,6\n")
        with open(os.path.join(dr, "B_vs_A", "differential_genes.csv"), "w") as fh:
            fh.write('"GeneID","baseMean"\n"g1",1\n"g2",1\n"g3",1\n"g4",1\n')
        run([os.path.join(SCRIPTS, "generate_pipeline_summary.py"), "--results_dir", res,
             "--group_col", "group", "--out_file", os.path.join(res, "PIPELINE_SUMMARY.md")])
        ps = open(os.path.join(res, "PIPELINE_SUMMARY.md")).read()
        check("summary reports the protein-evidence share",
              "with protein-level evidence: 3 (50.0% of genes, 66.1% of counts; per sample 57.1-85.0%)" in ps, ps[:600])
        check("summary, all mode: share of DEGs with protein-level evidence",
              "50.0% with protein-level evidence" in ps, ps)

        # No annotated filtered gene
        build(w, ["g7"])
        p = select(w, "coding", expect_ok=False)
        check("coding with no annotated gene stops with an error",
              p.returncode != 0 and "no filtered gene has an eggNOG-mapper hit" in p.stderr, p.stderr)
        p = select(w, "all", expect_ok=False)
        check("all with no annotated gene still runs", p.returncode == 0, p.stderr)

    print(f"\n{'All checks passed' if not FAILS else f'{len(FAILS)} check(s) failed'}")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
