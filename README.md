<p align="center">
  <img src="logos/metores.png" alt="MeToReS logo" width="220"/>
</p>

# MeToReS

Metatranscriptomics Toward Ecosystem Responses

A Snakemake pipeline for soil/environmental metatranscriptomics: individual
per-sample assembly → pooled non-redundant gene catalog → taxonomic
classification → differential expression → antimicrobial resistome
annotation → environmental chemistry integration.

---

## Pipeline overview

```
Raw reads
   │
   ▼
QC (FastQC) → Trimming (Trimmomatic) → rRNA removal (SortMeRNA)
   │
   ▼
Individual assembly per sample (MEGAHIT + Prodigal)
   │
   ▼
Pooled non-redundant gene catalog (CD-HIT-EST)
   │
   ├──► Taxonomic classification (Kaiju, ProGenomes)
   ├──► Quantification (Salmon)
   ├──► Differential expression (DESeq2, multi-contrast)
   ├──► Functional annotation (eggNOG-mapper)
   │
   ├──► [optional] Constrained ordination (db-RDA, vegan)
   ├──► [optional] Resistance annotation (CARD/RGI + BacMet/DIAMOND)
   │        └─ [optional] Resistome summary plots (top ARG/MRG families)
   ├──► [optional] Environmental chemistry integration
   │        ├─ envfit ordination
   │        ├─ partial correlation (ppcor)
   │        └─ MaAsLin2
   ├──► [optional] Bioindicator Venn (ppcor ∩ MaAsLin2 significant genes)
   └──► [optional] Co-selection heatmap (top resistome families × significant elements)
   │
   ▼
MultiQC report + per-analysis result tables and figures
```

Each `[optional]` block is gated by an `enabled:` flag in `config/config.yaml`
and only runs if turned on.

---

## Repository structure

```
MeToReS/
├── Snakefile
├── docs/
│   ├── logo.png                  # pipeline logo (shown at top of this README)
│   └── mobiles_logo.png          # MOBILES project logo (see Funding section)
├── installation/
│   ├── environment.yml          # main conda environment (everything except RGI)
│   ├── environment-rgi.yml      # RGI, kept isolated (see note below)
│   └── setup_databases.sh       # downloads/builds the large reference databases
├── config/
│   ├── config.yaml              # pipeline configuration — edit this per run
│   └── cluster.yaml             # SLURM resource profile (optional)
├── samplesheets/
│   ├── samplesheet_mobiles.csv  # sample metadata — replace with your own
│   └── *.xlsx                   # environmental chemistry data (if using that module)
├── scripts/                     # all R/Python analysis scripts, called by the Snakefile
├── databases/                   # reference databases (see Databases section below)
├── fastq/                       # raw/downloaded reads (your raw sequences)
├── results/                     # pipeline outputs (the results)
└── logs/                        # run logs
```

---

## Installation

### 1. Conda environments

Two environments, not one — see the comment block at the top of
`installation/environment-rgi.yml` for why RGI is kept separate (its own
install docs recommend this; it pins older BLAST/Python versions that
conflict with the rest of the stack).

```bash
conda env create -f installation/environment.yml
conda env create -f installation/environment-rgi.yml
```

Activate the main environment for everything except the resistance-annotation
step:

```bash
conda activate metores
```

### 2. Databases

Some reference data is small and permissively licensed enough to ship
directly in this repo (already present under `databases/`):

- `databases/bacmet/` — BacMet2 Experimentally Confirmed database (CC license)
- `databases/cog_funclass.tab` — NCBI COG functional category table
- `databases/trimmomatic_adapters/` — standard Trimmomatic adapter sequences (these ship with the Trimmomatic package itself, bundled here for convenience)

Everything else is too large (or, in CARD's case, licensed in a way that
doesn't clearly permit redistribution) to bundle, and is fetched by:

```bash
conda activate metores
bash installation/setup_databases.sh
```

This builds/downloads:
- **Kaiju** ProGenomes index (`kaiju-makedb`)
- **CARD** (`card.json`, for RGI)
- **eggNOG-mapper** database (large — budget real disk space and time)
- **SortMeRNA** rRNA reference database

---

## Compute resources

Per-rule `threads:` and `resources: mem_mb=...` values in the Snakefile are
tuned for a specific hardware profile (a shared server with a fixed
per-user core allocation). If you're running on different hardware,
adjust these — in particular, the thread counts on the per-sample rules
(`QC_rmrRNA`, `individual_assembly`, `kaiju_per_sample`,
`gene_expression_quant`) are set to divide evenly into the available core
count, so several samples can run concurrently without leaving cores
idle. Check your own CPU allocation before assuming the defaults suit
your setup, and adjust the `threads:` values (and, if you increase
concurrency substantially, the `mem_mb` budget) accordingly.

---

## Configuration

1. Copy your sample metadata into `samplesheets/` (must have at minimum
   `sample`, `sra` columns, plus whatever grouping/covariate columns your
   analysis needs — see `config.yaml`'s `group_col`).
2. Edit `config/config.yaml`:
   - `group_col` / `ref_level` — your primary grouping variable and reference
     level, used for output folder naming and the community-level analyses
     (constrained ordination, chemical ordination, partial correlation,
     MaAsLin2)
   - `deseq2.group_col` / `deseq2.ref_level` — optional override, used only
     for DESeq2's contrasts. Leave unset to reuse `group_col`/`ref_level`.
     Useful when DESeq2 needs a finer-grained grouping than the rest of the
     pipeline (e.g. splitting one level of `group_col` into two DESeq2-level
     groups without affecting the other analyses)
   - `deseq2.ref_level_filter` — optional `"column:value"` filter to
     disambiguate the reference level when its name alone doesn't uniquely
     identify the intended reference samples
   - `deseq2.formula` — DESeq2 design formula (the grouping variable must be
     the **last** term — the pipeline uses it to determine the contrast factor)
   - Toggle `enabled:` under `constrained_ordination`, `resistance`,
     `chemical`, `bioindicators`, `maaslin2`, `resistome_summary_plots`,
     `bioindicator_venn`, and `coselection_heatmap` for whichever optional
     modules you want to run. `resistome_summary_plots` requires
     `resistance`; `bioindicator_venn` requires both `bioindicators` and
     `maaslin2`; `coselection_heatmap` requires both `resistance` and
     `chemical`
   - `qc_warnings` — optional thresholds (`max_duplication_pct`,
     `max_rrna_pct`, `min_salmon_mapped_pct`) for informational
     `QC_WARNING:` messages printed to the log during FastQC, rRNA removal,
     and quantification. Purely advisory — never affects pipeline outcome
   - `chemical.fractions` — maps each extraction fraction present in your
     chemistry spreadsheet to an output name (e.g.
     `"I:bioavailable,II:reducible,III:oxidizable,IV:residual,PT:pseudototal"`
     for a standard BCR sequential extraction). All fractions listed here
     are parsed and saved; `chemical.fraction` then selects which single
     one feeds the downstream analyses (envfit, correlations, MaAsLin2) for
     a given run
3. All paths in `config.yaml` use a `{base_dir}` placeholder that resolves
   automatically to wherever you cloned this repo — no manual path editing
   needed for a standard setup.

---

## Running

```bash
conda activate metores

# Always dry-run first
snakemake -n --cores <N>

# Real run
snakemake --cores <N> --resources mem_mb=<total_MB> --keep-going \
    2>&1 | tee -a logs/pipeline_$(date +%Y%m%d).log
```

For cluster/SLURM execution, see `config/cluster.yaml` and Snakemake's
`--cluster` / `--profile` options.

### Re-running a specific rule after editing a script

Snakemake tracks file timestamps, not script content, so editing a script
alone won't trigger a re-run. Force it explicitly:

```bash
snakemake --cores <N> --resources mem_mb=<total_MB> --keep-going \
    -R <rule_name> 2>&1 | tee -a logs/pipeline_rerun.log
```

---

## Output structure

All outputs land under `results/`, organized by pipeline stage:

```
results/
├── qc/                              # FastQC, assembly QC, MultiQC
├── reads/{sra}/                     # rmrRNA, Trimming
├── assembly/{sra}/                  # Individual Assembly
├── gene_catalog/                    # Pooled non-redundant catalog
├── kaiju/                           # Taxonomy
├── smr_index/                       # Shared SortMeRNA index
├── quant/                           # Salmon quantification
├── deg/<group_col>_deg_results/     # DESeq2 outputs, per contrast
├── emapper_<group_col>/             # eggNOG functional annotation
├── downstream_<group_col>/          # cross-contrast summaries, ordination, heatmaps
├── constrained_ordination_<group_col>/   # [if enabled]
├── resistance_<group_col>/               # [if enabled]
│   └── summary_plots/                    # [if resistome_summary_plots enabled]
├── chemical/                             # [if enabled]
├── chemical_<group_col>/                 # [if enabled]
│   └── coselection/                      # [if coselection_heatmap enabled]
├── bioindicators_<group_col>/            # [if enabled]
│   └── venn/                             # [if bioindicator_venn enabled]
├── maaslin2_<group_col>/                 # [if enabled]
└── PIPELINE_SUMMARY.md              # key metrics from across the whole run,
                                      # in one human-readable file
```

---

## Licensing notes on bundled/fetched data

- **CARD** is fetched, not bundled, because its terms restrict commercial
  reproduction (non-commercial/research use is unrestricted). See
  https://card.mcmaster.ca/about.
  
- **BacMet** is bundled directly, under its CC Attribution license.

- Code in this repository is licensed under [LICENSE](LICENSE) — see that
  file for terms. This does not extend to the licenses of the reference
  databases themselves, which remain governed by their original sources.

---

## Citation

If you use MeToReS, please cite [manuscript citation once available].

MeToReS relies on a number of external tools and databases (assembly,
taxonomy, quantification, statistics, and resistome/chemistry
annotation). If you publish results produced with this pipeline, please
also cite the underlying tools and databases you used — see
[`CITATIONS.md`](CITATIONS.md) for the full list.

---

## Funding

<p align="center">
  <img src="logos/mobiles_logo.png" alt="MOBILES project logo" width="200"/>
</p>

This work was supported by the European Union through project 101135402 -
MOBILES - Monitoring and Detection of Biotic and Abiotic Pollutants by
Electronic, Plants and Microorganisms Based Sensors
[HORIZON-CL6-2023-ZEROPOLLUTION-01]

- Website: https://www.mobiles-project.eu/
- LinkedIn: https://www.linkedin.com/company/mobiles-project/
- YouTube: https://www.youtube.com/@MOBILES-project
- X (Twitter): https://x.com/mobiles_project
