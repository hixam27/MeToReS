<p align="center">
  <img src="logos/metores.png" alt="MeToReS logo" width="220"/>
</p>

# Citing MeToReS's dependencies

If you publish results produced with MeToReS, please cite the tools and
databases that your specific run actually used, in addition to the MeToReS
manuscript itself. Which of these apply depends on which optional modules
you enabled (see `config/config.yaml`).

## Core pipeline (always used)

- **FastQC** — Andrews S. FastQC: a quality control tool for high throughput
  sequence data. Babraham Bioinformatics.
  https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **Trimmomatic** — Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible
  trimmer for Illumina sequence data. Bioinformatics. 2014;30(15):2114–2120.
  doi:10.1093/bioinformatics/btu170
- **SortMeRNA** — Kopylova E, Noé L, Touzet H. SortMeRNA: fast and accurate
  filtering of ribosomal RNAs in metatranscriptomic data. Bioinformatics.
  2012;28(24):3211–3217. doi:10.1093/bioinformatics/bts611
- **MEGAHIT** — Li D, Liu C-M, Luo R, Sadakane K, Lam T-W. MEGAHIT: an
  ultra-fast single-node solution for large and complex metagenomics
  assembly via succinct de Bruijn graph. Bioinformatics.
  2015;31(10):1674–1676. doi:10.1093/bioinformatics/btv033
- **Prodigal** — Hyatt D, Chen G-L, LoCascio PF, Land ML, Larimer FW,
  Hauser LJ. Prodigal: prokaryotic gene recognition and translation
  initiation site identification. BMC Bioinformatics. 2010;11:119.
  doi:10.1186/1471-2105-11-119
- **CD-HIT / CD-HIT-EST** — Fu L, Niu B, Zhu Z, Wu S, Li W. CD-HIT:
  accelerated for clustering the next-generation sequencing data.
  Bioinformatics. 2012;28(23):3150–3152. doi:10.1093/bioinformatics/bts565
- **seqkit** — Shen W, Le S, Li Y, Hu F. SeqKit: a cross-platform and
  ultrafast toolkit for FASTA/Q file manipulation. PLoS ONE.
  2016;11(10):e0163962. doi:10.1371/journal.pone.0163962
- **MetaQUAST** — Mikheenko A, Saveliev V, Gurevich A. MetaQUAST:
  evaluation of metagenome assemblies. Bioinformatics. 2016;32(7):1088–1090.
  doi:10.1093/bioinformatics/btv697
- **Kaiju** — Menzel P, Ng KL, Krogh A. Fast and sensitive taxonomic
  classification for metagenomics with Kaiju. Nat Commun. 2016;7:11257.
  doi:10.1038/ncomms11257
- **ProGenomes3** (Kaiju's reference database) — Fullam A, Letunic I,
  Schmidt TSB, et al. proGenomes3: approaching one million accurately and
  consistently annotated high-quality prokaryotic genomes. Nucleic Acids
  Res. 2023;51(D1):D760–D766. doi:10.1093/nar/gkac1078
- **Salmon** — Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon
  provides fast and bias-aware quantification of transcript expression.
  Nat Methods. 2017;14(4):417–419. doi:10.1038/nmeth.4197
- **tximport** — Soneson C, Love MI, Robinson MD. Differential analyses for
  RNA-seq: transcript-level estimates improve gene-level inferences.
  F1000Research. 2015;4:1521. doi:10.12688/f1000research.7563.2
- **DESeq2** — Love MI, Huber W, Anders S. Moderated estimation of fold
  change and dispersion for RNA-seq data with DESeq2. Genome Biol.
  2014;15:550. doi:10.1186/s13059-014-0550-8
- **eggNOG-mapper v2** — Cantalapiedra CP, Hernández-Plaza A, Letunic I,
  Bork P, Huerta-Cepas J. eggNOG-mapper v2: functional annotation,
  orthology assignments, and domain prediction at the metagenomic scale.
  Mol Biol Evol. 2021;msab293. doi:10.1093/molbev/msab293
- **eggNOG 5.0** (database) — Huerta-Cepas J, Szklarczyk D, Heller D, et al.
  eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated
  orthology resource based on 5090 organisms and 2502 viruses. Nucleic
  Acids Res. 2019;47(D1):D309–D314. doi:10.1093/nar/gky1085
- **DIAMOND** — Buchfink B, Reuter K, Drost H-G. Sensitive protein
  alignments at tree-of-life scale using DIAMOND. Nat Methods.
  2021;18:366–368. doi:10.1038/s41592-021-01101-x
- **MultiQC** — Ewels P, Magnusson M, Lundin S, Käller M. MultiQC:
  summarize analysis results for multiple tools and samples in a single
  report. Bioinformatics. 2016;32(19):3047–3048.
  doi:10.1093/bioinformatics/btw354
- **Snakemake** — Mölder F, Jablonski KP, Letcher B, et al. Sustainable
  data analysis with Snakemake. F1000Research. 2021;10:33.
  doi:10.12688/f1000research.29032.2

## `constrained_ordination` module

- **vegan** (R package) — Oksanen J, Simpson G, Blanchet F, et al. vegan:
  Community Ecology Package. R package version 2.7-3. CRAN.
  https://CRAN.R-project.org/package=vegan

## `resistance` module

- **CARD / RGI** — Alcock BP, Huynh W, Chalil R, et al. CARD 2023:
  expanded curation, support for machine learning, and resistome
  prediction at the Comprehensive Antibiotic Resistance Database. Nucleic
  Acids Res. 2023;51(D1):D690–D699. doi:10.1093/nar/gkac920
- **BacMet** — Pal C, Bengtsson-Palme J, Rensing C, Kristiansson E,
  Larsson DGJ. BacMet: antibacterial biocide and metal resistance genes
  database. Nucleic Acids Res. 2014;42(D1):D737–D743.
  doi:10.1093/nar/gkt1252

## `chemical` / `bioindicators` / `maaslin2` modules

- **ppcor** (R package) — Kim S. ppcor: an R package for a fast
  calculation to semi-partial correlation coefficients. Commun Stat Appl
  Methods. 2015;22(6):665–674. doi:10.5351/CSAM.2015.22.6.665
- **MaAsLin2** — Mallick H, Rahnavard A, McIver LJ, et al. Multivariable
  association discovery in population-scale meta-omics studies. PLoS
  Comput Biol. 2021;17(11):e1009442. doi:10.1371/journal.pcbi.1009442

## R packages used throughout (`downstream_analysis.R` and others)

- **GO.db**, **ComplexUpset**, **RColorBrewer**, **ape** — cite via
  `citation("packagename")` in R for the exact, version-matched wording.
