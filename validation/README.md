# Validation scripts

Scripts behind the validation reported in the MeToReS paper. They read
results folders and published files; none of them changes pipeline outputs.
Run them from the MeToReS folder with the `metores` environment active.

| Script | What it does |
|---|---|
| `compare_runs.py` | Compares two results folders metric by metric and table by table. `--mode reproducibility` with `--tolerance exact` (level 1a: re-run from fixed Salmon output; tables from gene selection onward must be identical byte for byte), `1b` (re-run from the fixed catalogue; small Salmon differences within set limits) or `level2` (end to end on another machine; genes compared through their KEGG orthologues, criteria listed at the top of the script and fixed before the level-2 run). |
| `paper_numbers.py` | Prints every number the manuscript takes from one run; with Snakemake log files, also wall-clock time per rule. |
| `forhot_benchmark.R` | Compares the ForHot run with the published read-based analysis (Söllinger *et al.*, 2022) at the KO level: the comparison criteria, their sensitivity to normalisation and design, and the assembly-origin diagnostics. |
| `forhot_checks.R` | Checks on the ForHot run: sample mapping against the published per-sample KO counts, rRNA and taxonomy of the filtered genes, and the share of reads they hold. |

**Published ForHot files** (DataverseNO doi:10.18710/HZBZQE, CC0) used by the
two R scripts: `001_metadata.txt`, `002a_KO_tax_table_raw.txt`,
`005a/b/c_DGE_Analysis_*.txt` and `006a_proteinBS_oxidP_mapping_file.txt`.

**Examples**

```bash
# Level 1a, 1b and 2 (second run in results_forhot_1a, _1b or a copy from another machine)
python validation/compare_runs.py results_forhot results_forhot_1a \
    --mode reproducibility --tolerance exact --group_col group --out validation/forhot_level1a.md
python validation/compare_runs.py results_forhot results_forhot_1b \
    --mode reproducibility --tolerance 1b --group_col group \
    --samplesheet examples/forhot/samplesheet_forhot.csv --out validation/forhot_level1b.md
python validation/compare_runs.py results_forhot results_forhot_external \
    --mode reproducibility --tolerance level2 --group_col group \
    --samplesheet examples/forhot/samplesheet_forhot.csv --formula "~ gradient + group" \
    --out validation/forhot_level2.md

# Numbers for the paper
python validation/paper_numbers.py results_forhot group > validation/paper_numbers.txt

# Comparison with the published analysis
Rscript validation/forhot_benchmark.R --results results_forhot \
    --samplesheet examples/forhot/samplesheet_forhot.csv \
    --published /path/to/HZBZQE --out validation/forhot_benchmark
```

Outputs written into `validation/` are ignored by git; only the scripts and
this README are tracked.
