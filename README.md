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
   ├──► Taxonomic classification (Kaiju, proGenomes)
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
   ├──► [optional] Co-selection heatmap (top resistome families × chemical elements)
   └──► Safeguards S1–S4 (low-replication checks)
   │
   ▼
MultiQC report, PIPELINE_SUMMARY.md, SAFEGUARDS_SUMMARY.md + result tables and figures
```

Each `[optional]` block is gated by an `enabled:` flag in `config/config.yaml`
and only runs if turned on.

---

## Repository structure

```
MeToReS/
├── Snakefile
├── LICENSE                       # MIT licence
├── CITATIONS.md                  # tools and databases to cite
├── FUNDING.md                    # funding acknowledgement
├── Architecture_Diagram.png      # overview of the pipeline steps
├── logos/                        # MeToReS and MOBILES logos
├── installation/
│   ├── environment.yml           # main conda environment (everything except RGI)
│   ├── environment.lock.yml      # exact versions used for validation (Linux 64-bit)
│   ├── environment-rgi.yml       # RGI, kept in its own environment (see below)
│   ├── environment-rgi.lock.yml  # exact RGI versions used for validation
│   └── setup_databases.sh        # downloads/builds the large reference databases
├── config/
│   ├── config.yaml               # configuration — adapt per run
│   ├── cluster.yaml              # example SLURM resources (optional)
│   └── Config_Structure.png      # map of the configuration options
├── samplesheets/
│   ├── samplesheet_template.csv  # sample sheet template
│   └── chemistry_template.xlsx   # chemistry workbook template (optional module)
├── scripts/                      # R/Python analysis scripts called by the Snakefile
├── tests/
│   └── test_safeguards.py        # positive and negative controls for safeguards S1–S4
├── databases/                    # reference databases (see Databases below)
├── fastq/                        # input reads
├── results/                      # pipeline outputs
└── logs/                         # run logs
```

---

## Installation

### Requirements

- **Linux (64-bit).** The pipeline and its conda environments are built and
  validated on Linux.
- **Memory:** at least ~100 GB of RAM for the Kaiju step, whose proGenomes
  index (89 GB) is loaded entirely into memory; the other steps need less.
- **Disk:** about 140 GB for the reference databases (table below), plus space
  for the results, which take several times the size of the raw reads.

### 1. Conda environments

MeToReS uses two environments: the main one, and a separate one for RGI
(CARD), as recommended by RGI's developers — RGI pins older BLAST/Python
versions that conflict with the rest of the stack.

```bash
conda env create -f installation/environment.yml
conda env create -f installation/environment-rgi.yml
conda activate metores
```

To rebuild exactly the software versions used to validate MeToReS v1.0.0,
use the lock files instead (Linux 64-bit only):

```bash
conda env create -f installation/environment.lock.yml
conda env create -f installation/environment-rgi.lock.yml
```

The pipeline finds the RGI environment by its name (`metores-rgi`); you
don't need to activate it.

### 2. Databases

Three small reference sets are included in `databases/`:

- `bacmet/` — BacMet2 experimentally confirmed database (CC BY licence)
- `cog_funclass.tab` — NCBI COG functional categories
- `trimmomatic_adapters/` — standard Trimmomatic adapter sequences

The large ones are downloaded by:

```bash
conda activate metores
bash installation/setup_databases.sh
```

| Database | Used for | Size on disk |
|---|---|---|
| Kaiju proGenomes (pre-built index) | taxonomic classification | ~89 GB |
| eggNOG-mapper | functional annotation | ~48 GB |
| SortMeRNA | rRNA removal | ~3 GB |
| CARD (`card.json`) | antibiotic resistance genes (RGI) | ~52 MB |

CARD is downloaded in its latest release; note the version you used when
reporting results, since it changes over time.

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

## Safeguards

With few samples per group, differences between sites or single samples can
look like biological effects. MeToReS runs four checks by default, collects
them in `results/safeguards/SAFEGUARDS_SUMMARY.md`, and prints each finding to
the log as a `SAFEGUARD_WARNING` line. Except for the S1 filter, they report
only and change no result.

| Check | What it does |
|---|---|
| S1 Prevalence | Filters the genes tested by MaAsLin2 and partial correlation to those detected in ≥ 50 % of samples, or in ≥ 3 samples of one group (`safeguards.prevalence`); reports the detection of every hit |
| S2 Library complexity and type | Flags samples with far fewer genes detected, ORFs assembled or reads mapped than the run median; checks each declared Salmon library type against auto-detection on a read subsample; flags DESeq2 contrasts that mix library types (also shown in the MultiQC report) |
| S3 Singular fits | Flags MaAsLin2 hits from models whose random effect was estimated at zero, so that the grouping was effectively ignored |
| S4 Assembly origin | Flags DESeq2 and MaAsLin2 hit sets enriched for genes assembled from one sample or site, compared with tested genes of similar expression and detection |

Thresholds are set in the `safeguards:` block of `config/config.yaml`.
`safeguards.singular_fits.within_between: true` adds a refit separating the
between-group and within-group variation of each element.

`python tests/test_safeguards.py` runs each safeguard on small synthetic data
with a known answer and checks that it flags what it should and nothing else
(about one minute; no databases needed).

---

## Input data requirements

MeToReS takes paired-end metatranscriptomic reads. Two properties of the
libraries change what its results can show:

- **Stranded libraries of one type.** The default library type is stranded
  (`salmon.libtype: ISR`). S2 checks each sample's declared type against
  Salmon's auto-detection and flags DESeq2 contrasts that mix types. A mixed
  design is a batch effect the pipeline reports but cannot remove; add
  `libtype` to `deseq2.formula` when it is not confounded with the groups.
- **Unique dual indexes.** Reads that leak between samples sharing a
  sequencing run cannot be recognised in the data. Without unique dual indexes,
  S4 cannot tell whether a hit set dominated by one sample's assembly reflects
  community composition or cross-sample leakage.

---

## Configuration

### Input files

**1. Samplesheet** — one row per sample. Start from
`samplesheets/samplesheet_template.csv`:

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | Sample name, used in every output (and matched to `SampleID` in the chemistry workbook) |
| `sra` | yes | Either the stem of your FASTQ files in fastq/, or a public run accession (SRR/ERR/DRR…); if the files are not in fastq/, the reads are downloaded automatically |
| grouping / covariate columns | as needed | Any columns your analysis uses (e.g. `group`, `season`); you refer to them by name in `config.yaml` |
| `deseq2_group` | optional | A separate grouping used only for the DESeq2 contrasts (see `deseq2.group_col`) |
| `libtype` | optional | Salmon library type for this sample (e.g. `IU` for an unstranded library); empty cells use `salmon.libtype` |

**2. Reads** — paired-end FASTQ files in `fastq/`, named
`<sra>_1.fastq.gz` and `<sra>_2.fastq.gz`.

**3. Chemistry workbook** (optional module) — start from
`samplesheets/chemistry_template.xlsx`: one row per sample, a `SampleID`
column matching `sample`, and one column per variable and fraction named
`Variable_Fraction` (e.g. `Cu_I`, `Cu_PT`). Empty cells are missing values.
Other descriptive columns (site, season...) belong in the samplesheet.
The fraction suffixes are declared in `chemical.fractions` (e.g.
`"I:bioavailable,PT:pseudototal"`), and `chemical.fraction` selects the one
used by the downstream analyses.

### config/config.yaml

`config/config.yaml` is a commented template. As provided, it points to the
example samplesheet (`samplesheets/samplesheet_template.csv`), so once the
databases are installed a dry run (`snakemake -n`) works straight away. Before
your first real run, adapt at least:

- `samplesheet` — path to your samplesheet
- `group_col` / `ref_level` — your main grouping column and its reference level
- `deseq2.group_col` / `deseq2.ref_level` — the DESeq2 grouping and reference
  (leave unset to reuse the two above)
- `deseq2.formula` — the DESeq2 design. Always write the grouping variable as
  `group` (the pipeline renames the grouping column to `group` internally) and
  put it last; other terms use their samplesheet column names, e.g.
  `"~ season + group"` or `"~ libtype + group"`
- `chemical.xlsx_file`, `chemical.sheet`, `chemical.rename_map` — your
  workbook, its sheet name, and any sample-name prefixes to rename (set
  `rename_map: ""` if none)
- the covariates of the optional modules: `constrained_ordination.variables`,
  `bioindicators.covariates`, `maaslin2.random_effects`

All paths use a `{base_dir}` placeholder that resolves to the folder where
you cloned the repository.

Optional modules (`constrained_ordination`, `resistance`, `chemical`,
`bioindicators`, `maaslin2`, `resistome_summary_plots`, `bioindicator_venn`,
`coselection_heatmap`) are switched on or off with their `enabled:` flag. In the
template, `resistance` and `resistome_summary_plots` are on; the modules that
need your own variables are off until you configure them.
`resistome_summary_plots` requires `resistance`; `bioindicator_venn` requires
`bioindicators` and `maaslin2`; `coselection_heatmap` requires `resistance`
and `chemical`.

### Main settings and defaults

| Setting | Default | What it does |
|---|---|---|
| `salmon.libtype` | `ISR` | Library type for all samples (stranded dUTP); per-sample override in the samplesheet. Safeguard S2 checks each declared type against Salmon's auto-detection on the first `salmon.probe_reads` (500,000) read pairs (`salmon.check_libtype`) |
| `gene_catalog.min_orf_len` | `150` | Minimum ORF length (bp) kept in the gene catalogue |
| `deseq2.sf_type` | `poscounts` | DESeq2 size-factor method, suited to sparse metatranscriptomic counts |
| `deseq2.independent_filtering` | `true` | DESeq2 independent filtering. Set to `false` to test every gene passing `deseq2.min_count`, so calls do not depend on a data-driven cutoff (fewer DEGs). `all_contrasts_summary.csv` reports the cutoff used and the number of borderline calls per contrast |
| `pvalue` / `fold_change` | `0.05` / `2` | Thresholds for calling differentially expressed genes |
| `resistance.rgi_models` | protein homolog model | CARD model types counted as ARGs; other types are reported separately |
| `resistance.rgi_min_pident` | `0` | Extra identity floor for CARD hits on top of CARD's curated Strict/Perfect cut-offs (e.g. `80` for a stricter check) |
| `resistance.min_pid` | `80` | Identity floor for BacMet hits (metal/biocide resistance genes) |
| `maaslin2.min_samples` / `maaslin2.min_fraction` | `10` / `0.7` | MaAsLin2 tests an element only if it is measured in at least this many samples and this fraction of samples; skipped elements are listed in `maaslin2_elements_tested.tsv` |
| `safeguards.prevalence.apply_filter` | `true` | Test in MaAsLin2 and partial correlation only genes detected in ≥ 50 % of samples, or in ≥ 3 samples of one group |
| `coselection_heatmap.element_selection` | `all` | Show all elements; `padj` or `pvalue` show only elements significant in envfit |
| `qc_warnings` | — | Thresholds for informational `QC_WARNING` messages (duplication, rRNA, mapping rate); never change results |

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

`config/cluster.yaml` lists example SLURM resources. It uses the `--cluster`
option of Snakemake ≤ 7; with Snakemake ≥ 8 (the version pinned here),
cluster execution goes through an executor plugin — see the
[Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

### Re-running a specific rule after editing a script

Snakemake does not track the content of the R/Python scripts, so editing a
script alone won't trigger a re-run. Force it explicitly:

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
├── safeguards/                      # S1–S4 reports + SAFEGUARDS_SUMMARY.md
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

## Authors

- **Hicham Klaina**¹
- **Lorenzo Pin**²
- **Anna Marzec-Grządziel**³
- **Francisco Falcone**¹
- **Giuseppe Diego Puglia**²
- **Jose Javier Astrain**¹

¹ Public University of Navarre, Pamplona, Spain ·
² Institute for Agriculture and Forestry Systems in the Mediterranean (ISAFoM), National Research Council of Italy (CNR), Catania, Italy ·
³ Institute of Soil Science and Plant Cultivation – State Research Institute, Puławy, Poland

Questions and bug reports: please open an issue on GitHub.

---

## Citation

A manuscript describing MeToReS is in preparation. Until it is published,
please cite the software:

> Klaina H, Pin L, Marzec-Grządziel A, Falcone F, Puglia GD, Astrain JJ.
> *MeToReS: Metatranscriptomics Toward Ecosystem Responses* (version 1.0.0)
> [Software]. https://github.com/hixam27/MeToReS

GitHub's "Cite this repository" link (in the "About" box on the right of the
repository page) gives this citation in APA and BibTeX formats.

MeToReS relies on a number of external tools and databases. If you publish
results produced with this pipeline, please also cite the underlying tools
and databases you used — see [`CITATIONS.md`](CITATIONS.md).

---

## Funding

<p align="center">
  <img src="logos/eu_emblem.png" alt="Funded by the European Union" height="80" align="middle"/>
  &nbsp;&nbsp;&nbsp;&nbsp;
  <img src="logos/mobiles_logo.png" alt="MOBILES project logo" height="120" align="middle"/>
</p>

Funded by the European Union. Project 101135402 — MOBILES — Monitoring and
Detection of Biotic and Abiotic Pollutants by Electronic, Plants and
Microorganisms Based Sensors [HORIZON-CL6-2023-ZEROPOLLUTION-01], funded by
the EU Horizon Europe programme.

Views and opinions expressed are however those of the author(s) only and do
not necessarily reflect those of the European Union or the European Research
Executive Agency. Neither the European Union nor the European Research
Executive Agency can be held responsible for them.

- Website: https://www.mobiles-project.eu/
- LinkedIn: https://www.linkedin.com/company/mobiles-project/
- YouTube: https://www.youtube.com/@MOBILES-project
- X (Twitter): https://x.com/mobiles_project
