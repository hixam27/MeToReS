# =============================================================================
# MeToReS — Snakefile
# =============================================================================

import pandas as pd
import os
import json
import subprocess

configfile: "config/config.yaml"

# =============================================================================
# Base paths
# =============================================================================
BASE_DIR    = os.environ.get("METORES_HOME", workflow.basedir)
SCRIPTS_DIR = f"{BASE_DIR}/scripts"


def resolve_path(value):
    if isinstance(value, str):
        return value.format(base_dir=BASE_DIR)
    if isinstance(value, list):
        return [resolve_path(v) for v in value]
    if isinstance(value, dict):
        return {k: resolve_path(v) for k, v in value.items()}
    return value


config = resolve_path(config)


def get_conda_env_path(env_name):
    result = subprocess.run(
        ["conda", "info", "--envs", "--json"],
        capture_output=True, text=True, check=True
    )
    envs = json.loads(result.stdout)["envs"]
    for e in envs:
        if os.path.basename(e) == env_name:
            return e
    raise RuntimeError(
        f"Conda environment '{env_name}' not found. "
        f"Run: conda env create -f environment-rgi.yml"
    )

R_DESEQ2         = f"{SCRIPTS_DIR}/DEG_analysis.r"
R_FILTER         = f"{SCRIPTS_DIR}/filter_genes.R"
R_DOWNSTREAM     = f"{SCRIPTS_DIR}/downstream_analysis.R"
R_CONSTRAINED    = f"{SCRIPTS_DIR}/constrained_ordination.R"
PY_DEG           = f"{SCRIPTS_DIR}/DEG_analysis.py"
PY_ASSEMBLY_QC   = f"{SCRIPTS_DIR}/assembly_qc.py"
PY_KAIJU         = f"{SCRIPTS_DIR}/kaiju_process.py"
PY_MERGE_TAX     = f"{SCRIPTS_DIR}/merge_taxonomy.py"
PY_PIPELINE_SUMMARY = f"{SCRIPTS_DIR}/generate_pipeline_summary.py"

samples_df      = pd.read_csv(config["samplesheet"])
SAMPLES         = samples_df["sample"].tolist()
SRAS            = samples_df["sra"].tolist()

GROUP_COL = config.get("group_col", "group")
if GROUP_COL not in samples_df.columns:
    raise ValueError(
        f"ERROR: Column '{GROUP_COL}' not found in samplesheet.\n"
        f"Available columns: {samples_df.columns.tolist()}\n"
        f"Set 'group_col' in config.yaml to one of the above."
    )

# ── DESeq2-specific grouping (decoupled from GROUP_COL) ──────────────────
# GROUP_COL above drives folder naming (emapper_{GROUP_COL}, etc.) and the
# community-level analyses (constrained_ordination, chemical_ordination,
# partial_correlation, maaslin2) — none of which actually read GROUP_COL as
# a variable; they take their grouping/covariate columns as independent,
# explicitly-configured strings (e.g. constrained_ordination.variables,
# maaslin2.random_effects), so they are unaffected by anything below.
#
# DEG_analysis, filter_genes, and resistance_deg_report, however, need to
# build/report on CONTRASTS, and sometimes a finer-grained grouping is
# needed there than the one used for community-level analyses — e.g. one
# country needs to be split into two DESeq2-level groups without splitting
# it everywhere else. config["deseq2"]["group_col"]/["ref_level"] override
# GROUP_COL/ref_level for exactly that purpose; if unset, both fall back to
# the global values, so this has no effect unless explicitly configured.
DESEQ2_CONFIG   = config.get("deseq2", {})
DESEQ_GROUP_COL = DESEQ2_CONFIG.get("group_col") or GROUP_COL
DESEQ_REF_LEVEL = DESEQ2_CONFIG.get("ref_level") or config.get("ref_level")

if DESEQ_GROUP_COL not in samples_df.columns:
    raise ValueError(
        f"ERROR: Column '{DESEQ_GROUP_COL}' (deseq2.group_col) not found in "
        f"samplesheet.\nAvailable columns: {samples_df.columns.tolist()}\n"
        f"Set 'deseq2.group_col' in config.yaml to one of the above, or "
        f"remove it to fall back to 'group_col'."
    )

DESEQ_GROUPS = samples_df[DESEQ_GROUP_COL].unique().tolist()
if DESEQ_REF_LEVEL not in DESEQ_GROUPS:
    raise ValueError(
        f"ERROR: deseq2 ref_level '{DESEQ_REF_LEVEL}' not found among "
        f"'{DESEQ_GROUP_COL}' values: {DESEQ_GROUPS}"
    )

comparison_name = "_vs_".join(sorted(DESEQ_GROUPS))
NON_REF_LEVELS  = [lvl for lvl in DESEQ_GROUPS if lvl != DESEQ_REF_LEVEL]
CONTRASTS       = [f"{lvl}_vs_{DESEQ_REF_LEVEL}" for lvl in NON_REF_LEVELS]
wildcard_constraints:
    direction = "up|down"

output_dir       = config["output_dir"]
fastq_dir        = config["fastq_dir"]

qc_dir           = f"{output_dir}/qc"
assembly_qc_dir  = f"{qc_dir}/assembly_qc"
multiqc_dir      = f"{qc_dir}/multiqc"

reads_dir        = f"{output_dir}/reads"
assembly_dir     = f"{output_dir}/assembly"
catalog_dir      = f"{output_dir}/gene_catalog"
kaiju_dir        = f"{output_dir}/kaiju"
smr_index_dir    = f"{output_dir}/smr_index"

quant_dir        = f"{output_dir}/quant"
deg_dir          = f"{output_dir}/deg"
emapper_dir      = f"{output_dir}/emapper_{GROUP_COL}"
downstream_dir   = f"{output_dir}/downstream_{GROUP_COL}"
DEG_RESULTS      = f"{deg_dir}/{GROUP_COL}_deg_results"

CLEAN_FASTQ1 = f"{reads_dir}/{{sra}}/rmrRNA/{{sra}}_rmrRNA.fastq.1.gz"
CLEAN_FASTQ2 = f"{reads_dir}/{{sra}}/rmrRNA/{{sra}}_rmrRNA.fastq.2.gz"


constrained_dir = f"{output_dir}/constrained_ordination_{GROUP_COL}"
CO_CONFIG       = config.get("constrained_ordination", {})
CO_ENABLED      = CO_CONFIG.get("enabled", False)
CO_VARIABLES    = CO_CONFIG.get("variables", [])
CO_DISTANCE     = CO_CONFIG.get("distance", "bray")
CO_CONDITIONED  = CO_CONFIG.get("conditioned", [])
CO_PERMUTATIONS = CO_CONFIG.get("permutations", 999)

RES_CONFIG    = config.get("resistance", {})
RES_ENABLED   = RES_CONFIG.get("enabled", False)
resistance_dir = f"{output_dir}/resistance_{GROUP_COL}"

RGI_BIN     = None
RGI_ENV_BIN = None
if RES_ENABLED:
    RGI_ENV_PREFIX = get_conda_env_path(RES_CONFIG.get("rgi_env_name", "metores-rgi"))
    RGI_BIN     = f"{RGI_ENV_PREFIX}/bin/rgi"
    RGI_ENV_BIN = f"{RGI_ENV_PREFIX}/bin"

PY_PARSE_CARD   = f"{SCRIPTS_DIR}/parse_card.py"
PY_PARSE_BACMET = f"{SCRIPTS_DIR}/parse_bacmet.py"
PY_MERGE_RES    = f"{SCRIPTS_DIR}/merge_resistance.py"
R_RES_REPORT    = f"{SCRIPTS_DIR}/resistance_deg_report.R"

# ── Resistome summary plots (top ARG/MRG families, faceted by country) ──
# Only meaningful if resistance annotation is enabled — reuses its outputs.
RESPLOT_CONFIG  = config.get("resistome_summary_plots", {})
RESPLOT_ENABLED = RESPLOT_CONFIG.get("enabled", False) and RES_ENABLED
RESPLOT_TOP_N   = RESPLOT_CONFIG.get("top_n", 10)
resistome_plot_dir = f"{resistance_dir}/summary_plots"
R_RESISTOME_PLOTS  = f"{SCRIPTS_DIR}/resistome_summary_plots.R"

CHEM_CONFIG   = config.get("chemical", {})
CHEM_ENABLED  = CHEM_CONFIG.get("enabled", False)
CHEM_FRACTION = CHEM_CONFIG.get("fraction", "bioavailable")

# All fraction *output filenames* are derived from this string.
CHEM_FRACTIONS_RAW   = CHEM_CONFIG.get("fractions", "I:bioavailable,PT:pseudototal")
CHEM_FRACTION_NAMES  = [pair.split(":")[1].strip()
                        for pair in CHEM_FRACTIONS_RAW.split(",") if pair.strip()]

chemical_dir       = f"{output_dir}/chemical"
chemical_ord_dir   = f"{output_dir}/chemical_{GROUP_COL}"

PY_PARSE_CHEM   = f"{SCRIPTS_DIR}/parse_chemical_data.py"
R_CHEM_ORD      = f"{SCRIPTS_DIR}/chemical_ordination.R"
R_CHEM_RES_CORR = f"{SCRIPTS_DIR}/resistance_chemical_correlation.R"

# ── Co-selection heatmap (top ARG/MRG vs. significant envfit elements) ──
# Only meaningful if both resistance annotation and chemical integration
# are enabled — needs outputs from both.
COSEL_CONFIG   = config.get("coselection_heatmap", {})
COSEL_ENABLED  = COSEL_CONFIG.get("enabled", False) and RES_ENABLED and CHEM_ENABLED
COSEL_TOP_N    = COSEL_CONFIG.get("top_n", 20)
COSEL_ENVFIT_ALPHA = COSEL_CONFIG.get("envfit_alpha", 0.05)
COSEL_CORR_ALPHA   = COSEL_CONFIG.get("corr_alpha", 0.05)
coselection_dir = f"{chemical_ord_dir}/coselection"
R_COSELECTION  = f"{SCRIPTS_DIR}/coselection_heatmap.R"

KAIJU_TAXON = config["kaiju"]["taxon_filter"]
KAIJU_RANK  = config["kaiju"]["rank"]


BIO_CONFIG   = config.get("bioindicators", {})
BIO_ENABLED  = BIO_CONFIG.get("enabled", False)
BIO_FRACTION = BIO_CONFIG.get("fraction", "bioavailable")
bioindic_dir = f"{output_dir}/bioindicators_{GROUP_COL}"
R_PARTIAL_CORR = f"{SCRIPTS_DIR}/partial_correlation.R"

ML_CONFIG   = config.get("maaslin2", {})
ML_ENABLED  = ML_CONFIG.get("enabled", False)
ML_FRACTION = ML_CONFIG.get("fraction", "bioavailable")
maaslin_dir = f"{output_dir}/maaslin2_{GROUP_COL}"
R_MAASLIN2  = f"{SCRIPTS_DIR}/run_maaslin2.R"

# ── High-confidence bioindicators (MaAsLin2 ∩ ppcor per element) ────────
# Only meaningful, and only wired into rule all, if BOTH bioindicators
# (ppcor) and maaslin2 are enabled — the intersection needs both methods'
# significant-gene tables to exist.
BIOVENN_CONFIG   = config.get("bioindicator_venn", {})
BIOVENN_ENABLED  = BIOVENN_CONFIG.get("enabled", False) and BIO_ENABLED and ML_ENABLED
BIOVENN_ELEMENTS = BIOVENN_CONFIG.get("elements", "Cu Co Cd Ni Pb")
biovenn_dir      = f"{output_dir}/bioindicators_{GROUP_COL}/venn"
R_BIOVENN        = f"{SCRIPTS_DIR}/bioindicator_venn.R"

rule all:
    input:
        f"{multiqc_dir}/fastqc/multiqc_report.html",
        expand(f"{assembly_dir}/{{sra}}/final.contigs.fa", sra=SRAS),
        expand(f"{assembly_dir}/{{sra}}/orfs.fasta",       sra=SRAS),
        f"{catalog_dir}/all_longest_orfs_cds_rmdup_id.fasta",
        f"{assembly_qc_dir}/assembly_qc_summary.json",
        f"{assembly_qc_dir}/assembly_qc_report.txt",
        f"{assembly_qc_dir}/length_dist.pdf",
        f"{assembly_qc_dir}/gc_dist.pdf",
        f"{assembly_qc_dir}/assembly_qc_mqc.yaml",
        f"{assembly_qc_dir}/individual_assembly_stats.tsv",
        f"{kaiju_dir}/kaiju_gene_taxonomy.tsv",
        f"{kaiju_dir}/kaiju_summary.json",
        f"{kaiju_dir}/kaiju_mqc.yaml",
        f"{kaiju_dir}/taxonomy_barplot_{KAIJU_RANK}.pdf",
        expand(f"{kaiju_dir}/{{sra}}_kaiju_summary.tsv", sra=SRAS),
        f"{DEG_RESULTS}/transcript_abundance_quantification_table_filter_taxonomy.csv",
        f"{DEG_RESULTS}/samples_PCA.pdf",
        f"{DEG_RESULTS}/vst_normalised_counts.csv",
        f"{DEG_RESULTS}/deseq2_dispersion.pdf",
        f"{DEG_RESULTS}/all_contrasts_summary.csv",
        expand(f"{DEG_RESULTS}/{{contrast}}/differential_genes.csv",
               contrast=CONTRASTS),
        expand(f"{DEG_RESULTS}/{{contrast}}/differential_genes_volcano.pdf",
               contrast=CONTRASTS),
        expand(f"{DEG_RESULTS}/{{contrast}}/upregulated_genes_id.txt",
               contrast=CONTRASTS),
        expand(f"{DEG_RESULTS}/{{contrast}}/downregulated_genes_id.txt",
               contrast=CONTRASTS),
        f"{emapper_dir}/all_genes_emapper.emapper.annotations",
        f"{multiqc_dir}/final/multiqc_report.html",
        f"{downstream_dir}/gene_pcoa.pdf",
        f"{downstream_dir}/anosim_result.csv",
        f"{downstream_dir}/upset_all_DEGs.pdf",
        f"{downstream_dir}/upset_upregulated.pdf",
        f"{downstream_dir}/all_contrast_DEG_counts.pdf",
        f"{downstream_dir}/kaiju_community_barplot.pdf",
        expand(f"{downstream_dir}/per_contrast/{{contrast}}/cog.pdf",
               contrast=CONTRASTS),
        expand(f"{downstream_dir}/per_contrast/{{contrast}}/go_rich_bar.pdf",
               contrast=CONTRASTS),
        expand(f"{downstream_dir}/per_contrast/{{contrast}}/ko_rich_bar.pdf",
               contrast=CONTRASTS),
        expand(f"{downstream_dir}/per_contrast/{{contrast}}/DEG_up_heatmap.pdf",
               contrast=CONTRASTS),
        expand(f"{downstream_dir}/per_contrast/{{contrast}}/DEG_down_heatmap.pdf",
               contrast=CONTRASTS),
        *([f"{constrained_dir}/dbrda_anova_results.csv",
           f"{constrained_dir}/variance_explained.csv"] if CO_ENABLED else []),
        *([f"{resistance_dir}/resistance_annotation.tsv",
           f"{resistance_dir}/abundance_with_resistance.csv",
           f"{resistance_dir}/arg_abundance_by_aro.csv",
           f"{resistance_dir}/mrg_abundance_by_compound.csv",
           f"{resistance_dir}/resistant_DEGs_annotated.tsv",
           f"{resistance_dir}/arg_expression_heatmap.pdf",
           f"{resistance_dir}/mrg_expression_heatmap.pdf"] if RES_ENABLED else []),
        *([f"{resistome_plot_dir}/resistome_ARG_barplot.pdf",
           f"{resistome_plot_dir}/resistome_MRG_barplot.pdf"] if RESPLOT_ENABLED else []),
        *([f"{coselection_dir}/coselection_heatmap.pdf"] if COSEL_ENABLED else []),
        *([f"{chemical_dir}/chemical_data_{CHEM_FRACTION}.csv",
           f"{chemical_ord_dir}/envfit_results.csv",
           f"{chemical_ord_dir}/chemical_pca.pdf",
           f"{chemical_ord_dir}/full_correlation_matrix.csv",
           f"{chemical_ord_dir}/element_compound_correlations.csv"] if CHEM_ENABLED else []),
        *([f"{bioindic_dir}/partial_correlations_significant.csv"] if BIO_ENABLED else []),
        *([f"{maaslin_dir}/maaslin2_significant.tsv"] if ML_ENABLED else []),
        *([f"{biovenn_dir}/high_confidence_bioindicators_summary.csv"] if BIOVENN_ENABLED else []),
           f"{output_dir}/PIPELINE_SUMMARY.md",


rule prefetch_sra2fastq:
    output:
        fastq1 = f"{fastq_dir}/{{sra}}_1.fastq.gz",
        fastq2 = f"{fastq_dir}/{{sra}}_2.fastq.gz"
    params:
        sra_id    = "{sra}",
        fastq_dir = fastq_dir
    resources:
        mem_mb  = 8000,
        runtime = 180
    threads: 4
    shell:
        """
        set -e
        mkdir -p {params.fastq_dir}
        SRA={params.sra_id}
        SRA_LEN=${{#SRA}}

        if   [ $SRA_LEN -eq 9  ]; then ENA_URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${{SRA:0:6}}/${{SRA}}"
        elif [ $SRA_LEN -eq 10 ]; then DIR="0${{SRA: -2}}"  && ENA_URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${{SRA:0:6}}/${{DIR}}/${{SRA}}"
        elif [ $SRA_LEN -eq 11 ]; then DIR="${{SRA: -3}}"    && ENA_URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${{SRA:0:6}}/${{DIR}}/${{SRA}}"
        else                           DIR="00${{SRA: -1}}"  && ENA_URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${{SRA:0:6}}/${{DIR}}/${{SRA}}"
        fi

        echo "SRA: ${{SRA}} | ENA URL: ${{ENA_URL}}"

        if wget --tries=3 --retry-connrefused --waitretry=5 \
               -O {output.fastq1} "${{ENA_URL}}/${{SRA}}_1.fastq.gz" 2>/dev/null; then
            wget --tries=3 --retry-connrefused --waitretry=5 \
                 -O {output.fastq2} "${{ENA_URL}}/${{SRA}}_2.fastq.gz"
            echo "Downloaded ${{SRA}} from ENA."
        else
            echo "ENA failed — trying NCBI SRA Toolkit..."
            rm -f {output.fastq1} {output.fastq2}
            rm -rf {params.fastq_dir}/${{SRA}}
            rm -f {params.fastq_dir}/${{SRA}}_1.fastq \
                  {params.fastq_dir}/${{SRA}}_2.fastq

            mkdir -p ~/.ncbi
            echo '/repository/remote/main/SDL.2/resolver-cgi = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"' >> ~/.ncbi/user-settings.mkfg
            echo '/tls/allow-all-certs = "true"' >> ~/.ncbi/user-settings.mkfg

            prefetch ${{SRA}} -O {params.fastq_dir}/ --max-size 50G || \
            prefetch ${{SRA}} -O {params.fastq_dir}/ --max-size 50G -t http

            fasterq-dump {params.fastq_dir}/${{SRA}}/${{SRA}}.sra \
                -O {params.fastq_dir}/ -e {threads} --split-files

            gzip -f {params.fastq_dir}/${{SRA}}_1.fastq
            gzip -f {params.fastq_dir}/${{SRA}}_2.fastq
            rm -rf {params.fastq_dir}/${{SRA}}/
            echo "Downloaded ${{SRA}} from NCBI SRA."
        fi
        """


rule QC_test:
    input:
        fastq1 = f"{fastq_dir}/{{sra}}_1.fastq.gz",
        fastq2 = f"{fastq_dir}/{{sra}}_2.fastq.gz"
    output:
        qc_dir_out = directory(f"{qc_dir}/{{sra}}_qc")
    params:
        max_duplication_pct = config.get("qc_warnings", {}).get("max_duplication_pct", 75)
    resources:
        mem_mb  = 4000,
        runtime = 60
    threads: 2
    shell:
        """
        mkdir -p {output.qc_dir_out}
        fastqc -t {threads} -o {output.qc_dir_out} --nogroup \
            {input.fastq1} {input.fastq2}

        # ── QC warning: sequence duplication (informational only, does
        # not affect pipeline outcome — see project notes on why per-sample
        # rules can only use fixed thresholds, not cross-sample statistics) ──
        # FastQC's "Total Deduplicated Percentage" is the % of the library
        # remaining after removing duplicates — i.e. HIGH value = LOW
        # duplication. Duplication% = 100 - that value. Checked on both
        # read directions; whichever is worse gets reported.
        for r in 1 2; do
            ZIP=$(ls {output.qc_dir_out}/*_${{r}}_fastqc.zip 2>/dev/null | head -1)
            if [ -n "$ZIP" ]; then
                DEDUP=$(unzip -p "$ZIP" '*/fastqc_data.txt' 2>/dev/null \
                    | grep "Total Deduplicated Percentage" \
                    | awk -F'\t' '{{print $2}}')
                if [ -n "$DEDUP" ]; then
                    DUP=$(awk -v d="$DEDUP" 'BEGIN{{printf "%.1f", 100-d}}')
                    OVER=$(awk -v dup="$DUP" -v t="{params.max_duplication_pct}" \
                        'BEGIN{{print (dup>t)?1:0}}')
                    if [ "$OVER" -eq 1 ]; then
                        echo "QC_WARNING: {wildcards.sra} read${{r}} duplication ${{DUP}}% exceeds threshold {params.max_duplication_pct}% (FastQC) — may indicate low-input/degraded library"
                    fi
                fi
            fi
        done
        """


rule multiqc_fastqc:
    input:
        fastqc_dirs = expand(f"{qc_dir}/{{sra}}_qc", sra=SRAS)
    output:
        html = f"{multiqc_dir}/fastqc/multiqc_report.html",
        data = directory(f"{multiqc_dir}/fastqc/multiqc_report_data")
    params:
        outdir = f"{multiqc_dir}/fastqc",
        title  = config["multiqc"]["title_fastqc"]
    resources:
        mem_mb  = 8000,
        runtime = 30
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        multiqc {qc_dir} \
            --outdir {params.outdir} \
            --title "{params.title}" \
            --filename multiqc_report \
            --module fastqc \
            --force --no-ansi \
        || echo "WARNING: MultiQC Run 1 completed with warnings"
        """


rule sortmerna_index:
    output:
        done = touch(f"{smr_index_dir}/.index_done")
    params:
        smr_bin = config["sortmerna"]["bin"],
        smr_db  = config["sortmerna"]["db"],
        smr_idx = config["sortmerna"]["idx_dir"],
        workdir = f"{smr_index_dir}/wd"
    resources:
        mem_mb  = 16000,
        runtime = 240
    threads: 16
    shell:
        """
        set -e
        mkdir -p {params.smr_idx} {params.workdir}
        rm -rf {params.workdir}/kvdb

        {params.smr_bin} \
            --ref     {params.smr_db} \
            --idx-dir {params.smr_idx} \
            --workdir {params.workdir} \
            --task 5 \
            --threads {threads}

        echo "SortMeRNA index ready:"
        ls -la {params.smr_idx}
        """


rule QC_rmrRNA:
    input:
        fastq1    = f"{fastq_dir}/{{sra}}_1.fastq.gz",
        fastq2    = f"{fastq_dir}/{{sra}}_2.fastq.gz",
        qc_done   = f"{qc_dir}/{{sra}}_qc",
        smr_index = f"{smr_index_dir}/.index_done"
    output:
        qc_control   = directory(f"{reads_dir}/{{sra}}/trimmed"),
        rmrna_dir    = directory(f"{reads_dir}/{{sra}}/rmrRNA"),
        clean_fastq1 = CLEAN_FASTQ1,
        clean_fastq2 = CLEAN_FASTQ2,
        smr_log      = f"{reads_dir}/{{sra}}/rmrRNA/{{sra}}_sortmerna.log"
    params:
        adapter_path = config["adapter_path"],
        smr_bin = config["sortmerna"]["bin"],
        smr_db       = config["sortmerna"]["db"],
        smr_idx      = config["sortmerna"]["idx_dir"],
        phred        = config["phred"],
        quality      = config["quality_cutoff"],
        min_len      = config["min_len"],
        workdir      = lambda w: f"{reads_dir}/{w.sra}/smr_wd",
        max_rrna_pct = config.get("qc_warnings", {}).get("max_rrna_pct", 15)
    resources:
        mem_mb  = 64000,
        runtime = 1440
    threads: 14
    shell:
        """
        set -e
        mkdir -p {output.qc_control} {output.rmrna_dir}

        trimmomatic -Xmx8g PE -threads {threads} {params.phred} \
            {input.fastq1} {input.fastq2} \
            {output.qc_control}/{wildcards.sra}_clean_forward_paired.fastq.gz \
            {output.qc_control}/{wildcards.sra}_clean_forward_unpaired.fastq.gz \
            {output.qc_control}/{wildcards.sra}_clean_reverse_paired.fastq.gz \
            {output.qc_control}/{wildcards.sra}_clean_reverse_unpaired.fastq.gz \
            ILLUMINACLIP:{params.adapter_path}-PE-2.fa:2:30:10 \
            LEADING:5 TRAILING:5 \
            SLIDINGWINDOW:{params.quality} MINLEN:{params.min_len}

        rm -rf {params.workdir}
        mkdir -p {params.workdir}

        {params.smr_bin} \
            --ref     {params.smr_db} \
            --reads   {output.qc_control}/{wildcards.sra}_clean_forward_paired.fastq.gz \
            --reads   {output.qc_control}/{wildcards.sra}_clean_reverse_paired.fastq.gz \
            --workdir {params.workdir} \
            --idx-dir {params.smr_idx} \
            --index 0 \
            --aligned {params.workdir}/rRNA \
            --other   {params.workdir}/non_rRNA \
            --fastx --paired_in --out2 \
            --threads {threads}

        FWD=$(ls {params.workdir}/non_rRNA_fwd.* 2>/dev/null | head -1)
        REV=$(ls {params.workdir}/non_rRNA_rev.* 2>/dev/null | head -1)
        if [ -z "$FWD" ] || [ -z "$REV" ]; then
            echo "ERROR: SortMeRNA output not found in {params.workdir}" >&2
            ls -la {params.workdir} >&2
            exit 1
        fi
        case "$FWD" in *.gz) mv "$FWD" {output.clean_fastq1} ;;
                       *)    gzip -c "$FWD" > {output.clean_fastq1} ;; esac
        case "$REV" in *.gz) mv "$REV" {output.clean_fastq2} ;;
                       *)    gzip -c "$REV" > {output.clean_fastq2} ;; esac

        if [ -f {params.workdir}/rRNA.log ]; then
            cp {params.workdir}/rRNA.log {output.smr_log}
        else
            touch {output.smr_log}
        fi

        # ── QC warning: rRNA content (informational only, does not affect
        # pipeline outcome) ── extracts the rRNA% SortMeRNA itself reports
        # ("Total reads passing E-value threshold = N (X.XX)"), the single
        # strongest individual quality signal identified for this dataset —
        # confirmed to cleanly separate known low-quality samples (22-36%)
        # from typical ones (1-8.5%) earlier in this project.
        RRNA_PCT=$(grep -oP "passing E-value threshold\s*=\s*\d+\s*\(\K[0-9.]+" \
            {output.smr_log} 2>/dev/null | head -1)
        if [ -n "$RRNA_PCT" ]; then
            OVER=$(awk -v r="$RRNA_PCT" -v t="{params.max_rrna_pct}" \
                'BEGIN{{print (r>t)?1:0}}')
            if [ "$OVER" -eq 1 ]; then
                echo "QC_WARNING: {wildcards.sra} rRNA content ${{RRNA_PCT}}% exceeds threshold {params.max_rrna_pct}% (SortMeRNA) — may indicate low-input/degraded library"
            fi
        fi

        rm -rf {params.workdir}

        echo "rRNA removal complete for {wildcards.sra}"
        """


rule individual_assembly:
    input:
        fastq1 = CLEAN_FASTQ1,
        fastq2 = CLEAN_FASTQ2
    output:
        contigs = f"{assembly_dir}/{{sra}}/final.contigs.fa",
        orfs    = f"{assembly_dir}/{{sra}}/orfs.fasta",
        prots   = f"{assembly_dir}/{{sra}}/proteins.fasta"
    params:
        outdir = f"{assembly_dir}/{{sra}}"
    resources:
        mem_mb  = 100000,
        runtime = 720
    threads: 14
    shell:
        """
        set -e
        mkdir -p {params.outdir}
        rm -rf {params.outdir}/megahit_output

        echo "Running MEGAHIT assembly for sample {wildcards.sra}..."
        megahit \
            -1 {input.fastq1} \
            -2 {input.fastq2} \
            -o {params.outdir}/megahit_output \
            --min-contig-len 200 \
            --k-min 21 --k-max 141 --k-step 12 \
            -t {threads} --memory 100e9

        cp {params.outdir}/megahit_output/final.contigs.fa {output.contigs}

        echo "Predicting ORFs with Prodigal..."
        prodigal \
            -i {output.contigs} \
            -d {output.orfs} \
            -a {output.prots} \
            -p meta -q

        sed -i "s/^>/>{wildcards.sra}_/" {output.orfs}
        sed -i "s/^>/>{wildcards.sra}_/" {output.prots}

        echo "Assembly of {wildcards.sra} complete."
        echo "  Contigs : $(grep -c '>' {output.contigs})"
        echo "  ORFs    : $(grep -c '>' {output.orfs})"
        """


rule pool_gene_catalog:
    input:
        orfs = expand(f"{assembly_dir}/{{sra}}/orfs.fasta", sra=SRAS)
    output:
        merged    = f"{catalog_dir}/all_orfs_merged.fasta",
        catalog   = f"{catalog_dir}/all_longest_orfs_cds_rmdup_id.fasta",
        clusters  = f"{catalog_dir}/all_longest_orfs_cds_rmdup_id.fasta.clstr"
    params:
        outdir     = catalog_dir,
        num_samples = len(SRAS)
    resources:
        mem_mb  = 64000,
        runtime = 480
    threads: 16
    shell:
        """
        set -e
        mkdir -p {params.outdir}

        echo "Pooling ORFs from {params.num_samples} samples..."
        cat {input.orfs} > {output.merged}
        echo "  Total ORFs (with redundancy): $(grep -c '>' {output.merged})"

        echo "Removing redundancy with CD-HIT-EST (95% identity)..."
        cd-hit-est \
            -i {output.merged} \
            -o {output.catalog} \
            -c 0.95 -n 10 -M 32000 -T {threads}

        echo "Non-redundant gene catalog:"
        echo "  Unique genes: $(grep -c '>' {output.catalog})"
        """


rule assembly_qc:
    input:
        contigs            = f"{catalog_dir}/all_longest_orfs_cds_rmdup_id.fasta",
        individual_contigs = expand(
            f"{assembly_dir}/{{sra}}/final.contigs.fa", sra=SRAS
        )
    output:
        json_summary     = f"{assembly_qc_dir}/assembly_qc_summary.json",
        text_report      = f"{assembly_qc_dir}/assembly_qc_report.txt",
        mqc_yaml         = f"{assembly_qc_dir}/assembly_qc_mqc.yaml",
        seqkit_tsv       = f"{assembly_qc_dir}/seqkit_stats.tsv",
        length_plot      = f"{assembly_qc_dir}/length_dist.pdf",
        gc_plot          = f"{assembly_qc_dir}/gc_dist.pdf",
        metaquast_dir    = directory(f"{assembly_qc_dir}/metaquast"),
        individual_stats = f"{assembly_qc_dir}/individual_assembly_stats.tsv"
    params:
        output_dir   = assembly_qc_dir,
        assembly_dir = assembly_dir,
        sras         = " ".join(SRAS),
        references   = config["assembly_qc"]["references"],
        min_contig   = config["assembly_qc"]["min_contig"],
        warn_n50     = config["assembly_qc"]["warn_n50"],
        warn_total   = config["assembly_qc"]["warn_total_bp"],
        warn_contigs = config["assembly_qc"]["warn_contigs"],
        threads      = config["threads"]
    resources:
        mem_mb  = 32000,
        runtime = 240
    threads: 8
    shell:
        """
        set -e
        mkdir -p {params.output_dir}

        python {PY_ASSEMBLY_QC} \
            -i {input.contigs} \
            -o {params.output_dir} \
            -r "{params.references}" \
            -t {params.threads} \
            --min_contig    {params.min_contig} \
            --warn_n50      {params.warn_n50} \
            --warn_total_bp {params.warn_total} \
            --warn_contigs  {params.warn_contigs}

        echo "Generating per-sample assembly statistics..."
        seqkit stats -T -a \
            {params.assembly_dir}/*/final.contigs.fa \
            > {output.individual_stats}.raw

        awk -v assembly_dir="{params.assembly_dir}" 'BEGIN{{FS=OFS="\\t"}}
            NR==1 {{ $1="sample"; print; next }}
            {{
                n = split($1, parts, "/")
                $1 = parts[n-1]
                print
            }}' {output.individual_stats}.raw > {output.individual_stats}
        rm -f {output.individual_stats}.raw

        echo "Per-sample assembly stats:"
        cat {output.individual_stats}
        """


rule kaiju_catalog:
    input:
        orfs       = f"{catalog_dir}/all_longest_orfs_cds_rmdup_id.fasta",
        qc_summary = f"{assembly_qc_dir}/assembly_qc_summary.json"
    output:
        gene_tax   = f"{kaiju_dir}/kaiju_gene_taxonomy.tsv",
        contam_ids = f"{kaiju_dir}/{KAIJU_TAXON.lower()}_ids.csv",
        summary    = f"{kaiju_dir}/kaiju_summary_{KAIJU_RANK}.tsv",
        json_out   = f"{kaiju_dir}/kaiju_summary.json",
        mqc_yaml   = f"{kaiju_dir}/kaiju_mqc.yaml",
        krona_in   = f"{kaiju_dir}/krona_input.txt",
        barplot    = f"{kaiju_dir}/taxonomy_barplot_{KAIJU_RANK}.pdf"
    params:
        output_dir   = kaiju_dir,
        kaiju_db     = config["kaiju"]["db"],
        nodes        = config["kaiju"]["nodes"],
        names        = config["kaiju"]["names"],
        taxon_filter = KAIJU_TAXON,
        rank         = KAIJU_RANK,
        top_n        = config["kaiju"]["top_n"],
        threads      = config["threads"]
    resources:
        mem_mb  = 64000,
        runtime = 480
    threads: 24
    shell:
        """
        set -e
        mkdir -p {params.output_dir}
        python {PY_KAIJU} catalog \
            -i  {input.orfs} \
            -o  {params.output_dir} \
            -db {params.kaiju_db} \
            -n  {params.nodes} \
            -nm {params.names} \
            -t  {params.threads} \
            -f  {params.taxon_filter} \
            -r  {params.rank} \
            --top_n {params.top_n}
        """


rule kaiju_per_sample:
    input:
        fastq1 = CLEAN_FASTQ1,
        fastq2 = CLEAN_FASTQ2
    output:
        kaiju_out = f"{kaiju_dir}/{{sra}}_kaiju.out",
        summary   = f"{kaiju_dir}/{{sra}}_kaiju_summary.tsv"
    params:
        output_dir   = kaiju_dir,
        kaiju_db     = config["kaiju"]["db"],
        nodes        = config["kaiju"]["nodes"],
        names        = config["kaiju"]["names"],
        rank         = KAIJU_RANK,
        threads      = config["threads"]
    resources:
        mem_mb  = 64000,
        runtime = 720
    threads: 14
    shell:
        """
        set -e
        mkdir -p {params.output_dir}
        python {PY_KAIJU} sample \
            -1  {input.fastq1} \
            -2  {input.fastq2} \
            -s  {wildcards.sra} \
            -o  {params.output_dir} \
            -db {params.kaiju_db} \
            -n  {params.nodes} \
            -nm {params.names} \
            -t  {params.threads} \
            -r  {params.rank}
        """


rule transcript_index:
    input:
        contigs    = f"{catalog_dir}/all_longest_orfs_cds_rmdup_id.fasta",
        qc_summary = f"{assembly_qc_dir}/assembly_qc_summary.json",
        kaiju_done = f"{kaiju_dir}/kaiju_summary.json"
    output:
        index = directory(f"{quant_dir}/catalog_index")
    resources:
        mem_mb  = 64000,
        runtime = 480
    threads: 16
    shell:
        """
        echo "Building Salmon index on pooled gene catalog..."
        salmon index -t {input.contigs} -i {output.index} -k 31
        echo "Salmon index built!"
        """


rule gene_expression_quant:
    input:
        fastq1 = CLEAN_FASTQ1,
        fastq2 = CLEAN_FASTQ2,
        index  = f"{quant_dir}/catalog_index"
    output:
        quant      = directory(f"{quant_dir}/{{sra}}_quant"),
        quant_file = f"{quant_dir}/{{sra}}_quant/quant.sf"
    params:
        threads = config["threads"],
        min_salmon_mapped_pct = config.get("qc_warnings", {}).get("min_salmon_mapped_pct", 20)
    resources:
        mem_mb  = 64000,
        runtime = 720
    threads: 14
    shell:
        """
        echo "Quantifying {wildcards.sra}..."
        salmon quant \
            -i {input.index} -l A \
            -1 {input.fastq1} -2 {input.fastq2} \
            -o {output.quant} \
            --threads {params.threads} \
            --validateMappings \
            --rangeFactorizationBins 4 \
            --seqBias --gcBias
        echo "Quantification complete for {wildcards.sra}!"

        # ── QC warning: mapping rate (informational only, does not affect
        # pipeline outcome) ── this metric was individually weaker than
        # rRNA%/duplication% for the samples flagged earlier this project
        # (not all of them showed low mapping specifically), but a low
        # mapping rate is still a reasonable general red flag for other
        # issues (contamination, wrong reference, etc.), so it's still
        # worth surfacing on its own.
        MAPPED_PCT=$(grep -o '"percent_mapped"[: ]*[0-9.]*' \
            {output.quant}/aux_info/meta_info.json 2>/dev/null \
            | grep -o '[0-9.]*$')
        if [ -n "$MAPPED_PCT" ]; then
            UNDER=$(awk -v m="$MAPPED_PCT" -v t="{params.min_salmon_mapped_pct}" \
                'BEGIN{{print (m<t)?1:0}}')
            if [ "$UNDER" -eq 1 ]; then
                echo "QC_WARNING: {wildcards.sra} Salmon mapping rate ${{MAPPED_PCT}}% below threshold {params.min_salmon_mapped_pct}% — check for contamination or reference mismatch"
            fi
        fi
        """


rule filter_genes:
    input:
        quant_files = expand(
            f"{quant_dir}/{{sra}}_quant/quant.sf", sra=SRAS
        ),
        samplesheet = config["samplesheet"]
    output:
        raw_matrix      = f"{DEG_RESULTS}/transcript_abundance_quantification_table.csv",
        filtered_matrix = f"{DEG_RESULTS}/transcript_abundance_quantification_table_filter.csv"
    params:
        quant_dir  = quant_dir,
        output_dir = DEG_RESULTS,
        group_col = DESEQ_GROUP_COL
    resources:
        mem_mb  = 16000,
        runtime = 60
    threads: 2
    shell:
        """
        set -e
        mkdir -p {params.output_dir}

        echo "Step 1: Building TPM abundance matrix..."
        python {PY_DEG} \
            -i {params.quant_dir} \
            -s {input.samplesheet} \
            -o {params.output_dir} \
            --group_col {params.group_col} \
            --matrix_only

        echo "Step 2: Dual-mode filtering (TPM + counts)..."
        Rscript {R_FILTER} \
            -i {output.raw_matrix} \
            -q {params.quant_dir} \
            -s {input.samplesheet} \
            -o {output.filtered_matrix} \
            -g {params.group_col}

        echo "Genes after filtering: $(tail -n +2 {output.filtered_matrix} | wc -l)"
        """


rule merge_taxonomy:
    input:
        filtered_matrix = f"{DEG_RESULTS}/transcript_abundance_quantification_table_filter.csv",
        gene_taxonomy   = f"{kaiju_dir}/kaiju_gene_taxonomy.tsv"
    output:
        annotated_matrix = f"{DEG_RESULTS}/transcript_abundance_quantification_table_filter_taxonomy.csv",
        merge_summary    = f"{DEG_RESULTS}/taxonomy_merge_summary.json"
    params:
        output_dir = DEG_RESULTS
    resources:
        mem_mb  = 8000,
        runtime = 30
    threads: 1
    shell:
        """
        set -e
        echo "Merging Kaiju taxonomy into abundance matrix..."
        python {PY_MERGE_TAX} \
            -m {input.filtered_matrix} \
            -t {input.gene_taxonomy} \
            -o {params.output_dir}
        echo "Taxonomy merge complete!"
        """


rule DEG_analysis:
    input:
        filtered_matrix = f"{DEG_RESULTS}/transcript_abundance_quantification_table_filter.csv",
        quant_files     = expand(
            f"{quant_dir}/{{sra}}_quant/quant.sf", sra=SRAS
        ),
        samplesheet = config["samplesheet"],
        reference   = f"{catalog_dir}/all_longest_orfs_cds_rmdup_id.fasta"
    output:
        pca_plot      = f"{DEG_RESULTS}/samples_PCA.pdf",
        vst_matrix    = f"{DEG_RESULTS}/vst_normalised_counts.csv",
        dispersion    = f"{DEG_RESULTS}/deseq2_dispersion.pdf",
        dds_rdata     = f"{DEG_RESULTS}/dds_object.RData",
        summary       = f"{DEG_RESULTS}/all_contrasts_summary.csv",
        contrast_csvs = expand(
            f"{DEG_RESULTS}/{{contrast}}/differential_genes.csv",
            contrast=CONTRASTS),
        contrast_volc = expand(
            f"{DEG_RESULTS}/{{contrast}}/differential_genes_volcano.pdf",
            contrast=CONTRASTS),
        contrast_up   = expand(
            f"{DEG_RESULTS}/{{contrast}}/upregulated_genes_id.txt",
            contrast=CONTRASTS),
        contrast_down = expand(
            f"{DEG_RESULTS}/{{contrast}}/downregulated_genes_id.txt",
            contrast=CONTRASTS),
        contrast_up_fasta = expand(
            f"{DEG_RESULTS}/{{contrast}}/upregulated_genes.fasta",
            contrast=CONTRASTS),
        contrast_down_fasta = expand(
            f"{DEG_RESULTS}/{{contrast}}/downregulated_genes.fasta",
            contrast=CONTRASTS)
    params:
        pvalue           = config["pvalue"],
        fold_change      = config["fold_change"],
        formula          = config["deseq2"]["formula"],
        ref_level        = DESEQ_REF_LEVEL,
        ref_level_filter = DESEQ2_CONFIG.get("ref_level_filter", ""),
        min_count        = config["deseq2"]["min_count"],
        output_dir       = DEG_RESULTS,
        group_col        = DESEQ_GROUP_COL
    resources:
        mem_mb  = 32000,
        runtime = 360
    threads: 4
    shell:
        """
        set -e
        echo "Running DESeq2: {params.group_col} (ref: {params.ref_level})"
        echo "Formula: {params.formula}"

        Rscript {R_DESEQ2} \
            -q {quant_dir} \
            -i {input.filtered_matrix} \
            -s {input.samplesheet} \
            -o {params.output_dir} \
            --group_col   "{params.group_col}" \
            --formula     "{params.formula}" \
            -p {params.pvalue} \
            -f {params.fold_change} \
            --ref_level   "{params.ref_level}" \
            --ref_level_filter "{params.ref_level_filter}" \
            --min_count   {params.min_count}

        echo "DESeq2 multi-contrast analysis complete!"
        cat {params.output_dir}/all_contrasts_summary.csv

        echo "Extracting sequences per contrast..."
        for contrast_dir in {params.output_dir}/*/; do
            contrast=$(basename "$contrast_dir")
            echo "  Processing: $contrast"

            if [ -s "$contrast_dir/upregulated_genes_id.txt" ]; then
                seqkit grep -f "$contrast_dir/upregulated_genes_id.txt" \
                    {input.reference} \
                    -o "$contrast_dir/upregulated_genes_dna.fasta"
                seqkit translate "$contrast_dir/upregulated_genes_dna.fasta" \
                    > "$contrast_dir/upregulated_genes.fasta"
            else
                touch "$contrast_dir/upregulated_genes_dna.fasta"
                touch "$contrast_dir/upregulated_genes.fasta"
            fi

            if [ -s "$contrast_dir/downregulated_genes_id.txt" ]; then
                seqkit grep -f "$contrast_dir/downregulated_genes_id.txt" \
                    {input.reference} \
                    -o "$contrast_dir/downregulated_genes_dna.fasta"
                seqkit translate "$contrast_dir/downregulated_genes_dna.fasta" \
                    > "$contrast_dir/downregulated_genes.fasta"
            else
                touch "$contrast_dir/downregulated_genes_dna.fasta"
                touch "$contrast_dir/downregulated_genes.fasta"
            fi
        done

        echo "Sequence extraction complete!"
        """


rule emapper_all_genes:
    input:
        reference       = f"{catalog_dir}/all_longest_orfs_cds_rmdup_id.fasta",
        filtered_matrix = f"{DEG_RESULTS}/transcript_abundance_quantification_table_filter.csv"
    output:
        annotations = f"{emapper_dir}/all_genes_emapper.emapper.annotations"
    params:
        eggnog_db   = config["eggnog_db"],
        threads     = config["threads"],
        emapper_dir = emapper_dir
    resources:
        mem_mb  = 32000,
        runtime = 1440
    threads: 32
    shell:
        """
        set -e
        mkdir -p {params.emapper_dir}

        echo "Extracting filtered gene IDs..."
        tail -n +2 {input.filtered_matrix} | cut -d',' -f1 | tr -d '"' \
            > {params.emapper_dir}/filtered_gene_ids.txt

        echo "Extracting CDS sequences..."
        seqkit grep -f {params.emapper_dir}/filtered_gene_ids.txt \
            {input.reference} \
            -o {params.emapper_dir}/filtered_cds.fasta

        echo "Translating to protein..."
        seqkit translate {params.emapper_dir}/filtered_cds.fasta \
            > {params.emapper_dir}/filtered_proteins.fasta

        echo "Running eggNOG-mapper on ALL filtered genes..."
        emapper.py \
            -i {params.emapper_dir}/filtered_proteins.fasta \
            --output all_genes_emapper \
            --output_dir {params.emapper_dir} \
            --data_dir {params.eggnog_db} \
            -m diamond \
            --cpu {params.threads} \
            --override
        """


rule multiqc_final:
    input:
        fastqc_dirs  = expand(f"{qc_dir}/{{sra}}_qc",            sra=SRAS),
        quant_files  = expand(
            f"{quant_dir}/{{sra}}_quant/quant.sf",                sra=SRAS
        ),
        assembly_mqc = f"{assembly_qc_dir}/assembly_qc_mqc.yaml",
        kaiju_mqc    = f"{kaiju_dir}/kaiju_mqc.yaml",
        deg_done     = f"{DEG_RESULTS}/all_contrasts_summary.csv",
        emapper_all  = f"{emapper_dir}/all_genes_emapper.emapper.annotations",
    output:
        html = f"{multiqc_dir}/final/multiqc_report.html",
        data = directory(f"{multiqc_dir}/final/multiqc_report_data")
    params:
        outdir = f"{multiqc_dir}/final",
        title  = config["multiqc"]["title_final"],
        search = f"{qc_dir} {quant_dir}"
    resources:
        mem_mb  = 16000,
        runtime = 60
    threads: 1
    shell:
        """
        mkdir -p {params.outdir}
        multiqc {params.search} \
            --outdir {params.outdir} \
            --title "{params.title}" \
            --filename multiqc_report \
            --force --no-ansi \
        || echo "WARNING: MultiQC Run 2 completed with warnings"
        echo "Full pipeline MultiQC report: {output.html}"
        """


rule downstream_analysis:
    input:
        matrix      = f"{DEG_RESULTS}/transcript_abundance_quantification_table_filter.csv",
        emapper_all = f"{emapper_dir}/all_genes_emapper.emapper.annotations",
        kaiju_cats  = expand(
            f"{kaiju_dir}/{{sra}}_kaiju_summary.tsv", sra=SRAS
        ),
        multiqc     = f"{multiqc_dir}/final/multiqc_report.html"
    output:
        pcoa        = f"{downstream_dir}/gene_pcoa.pdf",
        anosim      = f"{downstream_dir}/anosim_result.csv",
        upset_all   = f"{downstream_dir}/upset_all_DEGs.pdf",
        upset_up    = f"{downstream_dir}/upset_upregulated.pdf",
        deg_counts  = f"{downstream_dir}/all_contrast_DEG_counts.pdf",
        kaiju_bar   = f"{downstream_dir}/kaiju_community_barplot.pdf",

        cog         = expand(
            f"{downstream_dir}/per_contrast/{{contrast}}/cog.pdf",
            contrast=CONTRASTS),
        heatmap_up   = expand(
            f"{downstream_dir}/per_contrast/{{contrast}}/DEG_up_heatmap.pdf",
            contrast=CONTRASTS),
        heatmap_down = expand(
            f"{downstream_dir}/per_contrast/{{contrast}}/DEG_down_heatmap.pdf",
            contrast=CONTRASTS),
        go_bar      = expand(
            f"{downstream_dir}/per_contrast/{{contrast}}/go_rich_bar.pdf",
            contrast=CONTRASTS),
        ko_bar      = expand(
            f"{downstream_dir}/per_contrast/{{contrast}}/ko_rich_bar.pdf",
            contrast=CONTRASTS)
    params:
        workdir      = lambda wildcards: os.path.abspath(downstream_dir),
        deg_dir      = lambda wildcards: os.path.abspath(DEG_RESULTS),
        kaiju_dir    = lambda wildcards: os.path.abspath(kaiju_dir),
        emapper_dir  = lambda wildcards: os.path.abspath(emapper_dir),
        contrasts    = " ".join(CONTRASTS),
        taxon_filter = KAIJU_TAXON,
        samplesheet  = lambda wildcards: os.path.abspath(config["samplesheet"]),
        heatmap_min  = config.get("downstream_heatmap_min_sum", 50),
        top_n_kaiju  = config.get("downstream_top_n_kaiju", 15),
        cog_funclass = f"{BASE_DIR}/databases/cog_funclass.tab",
        group_col    = GROUP_COL
    resources:
        mem_mb  = 32000,
        runtime = 120
    threads: 4
    shell:
        """
        set -e
        mkdir -p {params.workdir}

        python3 - <<PYEOF
import pandas as pd, sys
df = pd.read_csv("{params.samplesheet}")
gcol = "{params.group_col}"
if "sample" not in df.columns or gcol not in df.columns:
    sys.exit(f"ERROR: samplesheet must have 'sample' and '{{gcol}}' columns")
out = df.rename(columns={{gcol: "group"}})[["sample", "group"]]
out.to_csv("{params.workdir}/sample_group.csv", index=False)
print(f"sample_group.csv written: {{len(out)}} samples, "
      f"{{out['group'].nunique()}} groups (from column '{{gcol}}')")
PYEOF

        Rscript {R_DOWNSTREAM} \
            --workdir      {params.workdir} \
            --deg_dir      {params.deg_dir} \
            --kaiju_dir    {params.kaiju_dir} \
            --emapper_dir  {params.emapper_dir} \
            --contrasts    "{params.contrasts}" \
            --taxon_filter {params.taxon_filter} \
            --heatmap_min_sum {params.heatmap_min} \
            --top_n_kaiju  {params.top_n_kaiju} \
            --cog_funclass {params.cog_funclass} \
            --samplesheet  {params.samplesheet}

        echo "Downstream analysis complete: {params.workdir}"
        """


rule constrained_ordination:
    input:
        vst         = f"{DEG_RESULTS}/vst_normalised_counts.csv",
        samplesheet = config["samplesheet"]
    output:
        summary = f"{constrained_dir}/dbrda_anova_results.csv",
        varexp  = f"{constrained_dir}/variance_explained.csv"
    params:
        workdir     = lambda wildcards: os.path.abspath(constrained_dir),
        vst         = lambda wildcards: os.path.abspath(f"{DEG_RESULTS}/vst_normalised_counts.csv"),
        samplesheet = lambda wildcards: os.path.abspath(config["samplesheet"]),
        variables   = " ".join(CO_VARIABLES),
        distance    = CO_DISTANCE,
        conditioned = " ".join(CO_CONDITIONED),
        permutations = CO_PERMUTATIONS
    resources:
        mem_mb  = 16000,
        runtime = 60
    threads: 2
    shell:
        """
        set -e
        mkdir -p {params.workdir}
        Rscript {R_CONSTRAINED} \
            --vst          {params.vst} \
            --samplesheet  {params.samplesheet} \
            --workdir      {params.workdir} \
            --variables    "{params.variables}" \
            --distance     {params.distance} \
            --conditioned  "{params.conditioned}" \
            --permutations {params.permutations}
        echo "Constrained ordination complete: {params.workdir}"
        """

rule catalog_proteins:
    input:
        catalog = f"{catalog_dir}/all_longest_orfs_cds_rmdup_id.fasta"
    output:
        faa = f"{resistance_dir}/catalog_proteins.faa"
    resources:
        mem_mb = 8000, runtime = 60
    threads: 4
    shell:
        """
        set -e
        mkdir -p {resistance_dir}
        seqkit translate --transl-table 11 --trim {input.catalog} > {output.faa}
        echo "Translated proteins: $(grep -c '>' {output.faa})"
        """


rule rgi_card:
    input:
        faa = f"{resistance_dir}/catalog_proteins.faa"
    output:
        rgi_txt = f"{resistance_dir}/card_rgi.txt",
        parsed  = f"{resistance_dir}/card_annotations.tsv"
    params:
        card_json   = RES_CONFIG.get("card_json"),
        outbase     = "card_rgi",
        keep        = RES_CONFIG.get("rgi_keep", "Perfect,Strict"),
        rgi_bin     = RGI_BIN,
        rgi_env_bin = RGI_ENV_BIN
    resources:
        mem_mb = 32000, runtime = 240
    threads: 32
    shell:
        """
        set -e
        cd {resistance_dir}
        export PATH="{params.rgi_env_bin}:$PATH"
        {params.rgi_bin} load --card_json {params.card_json} --local
        {params.rgi_bin} main \
            --input_sequence catalog_proteins.faa \
            --output_file {params.outbase} \
            --input_type protein \
            --alignment_tool DIAMOND \
            -n {threads} \
            --local --clean
        python {PY_PARSE_CARD} \
            --rgi {params.outbase}.txt \
            --keep "{params.keep}" \
            --out card_annotations.tsv
        """


rule bacmet_search:
    input:
        faa = f"{resistance_dir}/catalog_proteins.faa"
    output:
        dmnd   = temp(f"{resistance_dir}/bacmet_exp.dmnd"),
        hits   = f"{resistance_dir}/bacmet_hits.tsv",
        parsed = f"{resistance_dir}/bacmet_annotations.tsv"
    params:
        bacmet_fasta = RES_CONFIG.get("bacmet_exp"),
        bacmet_map   = RES_CONFIG.get("bacmet_map"),
        pid          = RES_CONFIG.get("min_pid", 80),
        qcov         = RES_CONFIG.get("min_qcov", 80),
        evalue       = RES_CONFIG.get("evalue", "1e-10")
    resources:
        mem_mb = 32000, runtime = 120
    threads: 32
    shell:
        """
        set -e
        mkdir -p {resistance_dir}
        diamond makedb --in {params.bacmet_fasta} --db {output.dmnd}
        diamond blastp \
            --query {input.faa} \
            --db {output.dmnd} \
            --outfmt 6 qseqid sseqid pident length qcovhsp evalue bitscore \
            --id {params.pid} \
            --query-cover {params.qcov} \
            --evalue {params.evalue} \
            --max-target-seqs 1 \
            --threads {threads} \
            --out {output.hits}
        python {PY_PARSE_BACMET} \
            --hits {output.hits} \
            --mapping {params.bacmet_map} \
            --out {output.parsed}
        """


rule merge_resistance:
    input:
        card     = f"{resistance_dir}/card_annotations.tsv",
        bacmet   = f"{resistance_dir}/bacmet_annotations.tsv",
        matrix   = f"{DEG_RESULTS}/transcript_abundance_quantification_table.csv",
        taxonomy = f"{kaiju_dir}/kaiju_gene_taxonomy.tsv",
        emapper  = f"{emapper_dir}/all_genes_emapper.emapper.annotations"
    output:
        annotation  = f"{resistance_dir}/resistance_annotation.tsv",
        ann_matrix  = f"{resistance_dir}/abundance_with_resistance.csv",
        by_aro      = f"{resistance_dir}/arg_abundance_by_aro.csv",
        by_compound = f"{resistance_dir}/mrg_abundance_by_compound.csv"
    resources:
        mem_mb = 16000, runtime = 30
    threads: 2
    shell:
        """
        set -e
        python {PY_MERGE_RES} \
            --matrix    {input.matrix} \
            --card      {input.card} \
            --bacmet    {input.bacmet} \
            --taxonomy  {input.taxonomy} \
            --emapper   {input.emapper} \
            --out_annotation {output.annotation} \
            --out_matrix     {output.ann_matrix} \
            --out_by_aro     {output.by_aro} \
            --out_by_compound {output.by_compound}
        """


rule resistance_deg_report:
    input:
        annotation  = f"{resistance_dir}/resistance_annotation.tsv",
        by_aro      = f"{resistance_dir}/arg_abundance_by_aro.csv",
        by_compound = f"{resistance_dir}/mrg_abundance_by_compound.csv",
        deg         = expand(f"{DEG_RESULTS}/{{contrast}}/differential_genes.csv",
                             contrast=CONTRASTS),
        samplesheet = config["samplesheet"]
    output:
        table       = f"{resistance_dir}/resistant_DEGs_annotated.tsv",
        arg_heatmap = f"{resistance_dir}/arg_expression_heatmap.pdf",
        mrg_heatmap = f"{resistance_dir}/mrg_expression_heatmap.pdf"
    params:
        deg_dir   = lambda wildcards: os.path.abspath(DEG_RESULTS),
        contrasts = " ".join(CONTRASTS),
        group_col = DESEQ_GROUP_COL,
        top_n     = RES_CONFIG.get("heatmap_top_n", 50),
        padj      = config.get("pvalue", 0.05)
    resources:
        mem_mb = 8000, runtime = 30
    threads: 2
    shell:
        """
        set -e
        Rscript {R_RES_REPORT} \
            --annotation     {input.annotation} \
            --deg_dir        {params.deg_dir} \
            --contrasts      "{params.contrasts}" \
            --by_aro         {input.by_aro} \
            --by_compound    {input.by_compound} \
            --samplesheet    {input.samplesheet} \
            --group_col      {params.group_col} \
            --out_table      {output.table} \
            --out_arg_heatmap {output.arg_heatmap} \
            --out_mrg_heatmap {output.mrg_heatmap} \
            --top_n          {params.top_n} \
            --padj           {params.padj}
        """


# =============================================================================
# Rule — Resistome summary plots (top ARG/MRG families, faceted by country)
# =============================================================================
rule resistome_summary_plots:
    input:
        abundance_with_resistance = f"{resistance_dir}/abundance_with_resistance.csv",
        samplesheet = config["samplesheet"]
    output:
        arg_barplot = f"{resistome_plot_dir}/resistome_ARG_barplot.pdf",
        mrg_barplot = f"{resistome_plot_dir}/resistome_MRG_barplot.pdf"
    params:
        out_dir = lambda w: os.path.abspath(resistome_plot_dir),
        top_n   = RESPLOT_TOP_N
    resources:
        mem_mb = 8000, runtime = 30
    threads: 2
    shell:
        """
        set -e
        Rscript {R_RESISTOME_PLOTS} \
            --abundance_with_resistance {input.abundance_with_resistance} \
            --samplesheet {input.samplesheet} \
            --top_n       {params.top_n} \
            --out_dir     {params.out_dir}
        """


# =============================================================================
# Rule — Co-selection heatmap (top ARG/MRG vs. significant envfit elements)
# =============================================================================
rule coselection_heatmap:
    input:
        full_corr     = f"{chemical_ord_dir}/full_correlation_matrix.csv",
        envfit        = f"{chemical_ord_dir}/envfit_results.csv",
        arg_abundance = f"{resistance_dir}/arg_abundance_by_aro.csv",
        mrg_abundance = f"{resistance_dir}/mrg_abundance_by_compound.csv"
    output:
        heatmap = f"{coselection_dir}/coselection_heatmap.pdf",
        data    = f"{coselection_dir}/coselection_heatmap_data.csv"
    params:
        out_dir       = lambda w: os.path.abspath(coselection_dir),
        top_n         = COSEL_TOP_N,
        envfit_alpha  = COSEL_ENVFIT_ALPHA,
        corr_alpha    = COSEL_CORR_ALPHA
    resources:
        mem_mb = 8000, runtime = 30
    threads: 2
    shell:
        """
        set -e
        Rscript {R_COSELECTION} \
            --full_corr      {input.full_corr} \
            --envfit         {input.envfit} \
            --arg_abundance  {input.arg_abundance} \
            --mrg_abundance  {input.mrg_abundance} \
            --top_n          {params.top_n} \
            --envfit_alpha   {params.envfit_alpha} \
            --corr_alpha     {params.corr_alpha} \
            --out_dir        {params.out_dir}
        """

rule parse_chemical_data:
    input:
        xlsx        = CHEM_CONFIG.get("xlsx_file"),
        samplesheet = config["samplesheet"]
    output:
        fraction_files = expand(f"{chemical_dir}/chemical_data_{{fraction}}.csv",
                                fraction=CHEM_FRACTION_NAMES),
        combined = f"{chemical_dir}/chemical_data_combined.csv",
        report   = f"{chemical_dir}/sample_matching_report.txt"
    params:
        sheet     = CHEM_CONFIG.get("sheet", "14_DATASET"),
        rename    = CHEM_CONFIG.get("rename_map", "CY_:CYP_"),
        fractions = CHEM_FRACTIONS_RAW
    resources:
        mem_mb = 4000, runtime = 15
    threads: 1
    shell:
        """
        set -e
        mkdir -p {chemical_dir}
        python {PY_PARSE_CHEM} \
            --xlsx        {input.xlsx} \
            --sheet       "{params.sheet}" \
            --samplesheet {input.samplesheet} \
            --rename      "{params.rename}" \
            --fractions   "{params.fractions}" \
            --out_dir     {chemical_dir}
        """


rule chemical_ordination:
    input:
        vst      = f"{DEG_RESULTS}/vst_normalised_counts.csv",
        chemical = f"{chemical_dir}/chemical_data_{CHEM_FRACTION}.csv"
    output:
        envfit_csv = f"{chemical_ord_dir}/envfit_results.csv",
        biplot     = f"{chemical_ord_dir}/envfit_biplot.pdf",
        chem_pca   = f"{chemical_ord_dir}/chemical_pca.pdf",
        loadings   = f"{chemical_ord_dir}/chemical_pca_loadings.csv"
    params:
        workdir  = lambda wildcards: os.path.abspath(chemical_ord_dir),
        vst      = lambda wildcards: os.path.abspath(f"{DEG_RESULTS}/vst_normalised_counts.csv"),
        chemical = lambda wildcards: os.path.abspath(f"{chemical_dir}/chemical_data_{CHEM_FRACTION}.csv"),
        fraction = CHEM_FRACTION,
        perms    = CHEM_CONFIG.get("envfit_permutations", 999)
    resources:
        mem_mb = 8000, runtime = 30
    threads: 2
    shell:
        """
        set -e
        Rscript {R_CHEM_ORD} \
            --vst            {params.vst} \
            --chemical       {params.chemical} \
            --workdir        {params.workdir} \
            --fraction_label {params.fraction} \
            --permutations   {params.perms}
        """


rule resistance_chemical_correlation:
    input:
        by_aro      = f"{resistance_dir}/arg_abundance_by_aro.csv",
        by_compound = f"{resistance_dir}/mrg_abundance_by_compound.csv",
        chemical    = f"{chemical_dir}/chemical_data_{CHEM_FRACTION}.csv"
    output:
        direct = f"{chemical_ord_dir}/element_compound_correlations.csv",
        full   = f"{chemical_ord_dir}/full_correlation_matrix.csv",
        heat   = f"{chemical_ord_dir}/correlation_heatmap.pdf"
    params:
        out_dir  = lambda wildcards: os.path.abspath(chemical_ord_dir),
        chemical = lambda wildcards: os.path.abspath(f"{chemical_dir}/chemical_data_{CHEM_FRACTION}.csv"),
        fraction = CHEM_FRACTION,
        method   = CHEM_CONFIG.get("correlation_method", "spearman")
    resources:
        mem_mb = 8000, runtime = 30
    threads: 2
    shell:
        """
        set -e
        Rscript {R_CHEM_RES_CORR} \
            --by_aro         {input.by_aro} \
            --by_compound    {input.by_compound} \
            --chemical       {params.chemical} \
            --fraction_label {params.fraction} \
            --method         {params.method} \
            --out_dir        {params.out_dir}
        """


rule partial_correlation:
    input:
        vst         = f"{DEG_RESULTS}/vst_normalised_counts.csv",
        chemical    = f"{chemical_dir}/chemical_data_{BIO_FRACTION}.csv",
        samplesheet = config["samplesheet"]
    output:
        sig = f"{bioindic_dir}/partial_correlations_significant.csv",
        all = f"{bioindic_dir}/partial_correlations_all.csv"
    params:
        out_dir    = lambda w: os.path.abspath(bioindic_dir),
        vst        = lambda w: os.path.abspath(f"{DEG_RESULTS}/vst_normalised_counts.csv"),
        chemical   = lambda w: os.path.abspath(f"{chemical_dir}/chemical_data_{BIO_FRACTION}.csv"),
        covariates = BIO_CONFIG.get("covariates", "season"),
        transform  = BIO_CONFIG.get("transform", "log10_zscore"),
        top_genes  = BIO_CONFIG.get("top_var_genes", 2000),
        method     = BIO_CONFIG.get("method", "spearman"),
        r_thr      = BIO_CONFIG.get("r_threshold", 0.6),
        padj_thr   = BIO_CONFIG.get("padj_threshold", 0.05)
    resources:
        mem_mb = 16000, runtime = 60
    threads: 4
    shell:
        """
        set -e
        Rscript {R_PARTIAL_CORR} \
            --vst {params.vst} --chemical {params.chemical} \
            --samplesheet {input.samplesheet} --out_dir {params.out_dir} \
            --covariates "{params.covariates}" --transform {params.transform} \
            --top_var_genes {params.top_genes} --method {params.method} \
            --r_threshold {params.r_thr} --padj_threshold {params.padj_thr}
        """


rule maaslin2:
    input:
        vst         = f"{DEG_RESULTS}/vst_normalised_counts.csv",
        chemical    = f"{chemical_dir}/chemical_data_{ML_FRACTION}.csv",
        samplesheet = config["samplesheet"]
    output:
        sig = f"{maaslin_dir}/maaslin2_significant.tsv",
        all = f"{maaslin_dir}/maaslin2_all_elements.tsv"
    params:
        rbin      = ML_CONFIG.get("rscript_bin", "Rscript"),
        out_dir   = lambda w: os.path.abspath(maaslin_dir),
        vst       = lambda w: os.path.abspath(f"{DEG_RESULTS}/vst_normalised_counts.csv"),
        chemical  = lambda w: os.path.abspath(f"{chemical_dir}/chemical_data_{ML_FRACTION}.csv"),
        elements  = ML_CONFIG.get("elements", ""),
        random    = ML_CONFIG.get("random_effects", "country"),
        fixcov    = ML_CONFIG.get("fixed_covariates", ""),
        transform = ML_CONFIG.get("transform", "log10_zscore"),
        topgenes  = ML_CONFIG.get("top_var_genes", 2000),
        maxsig    = ML_CONFIG.get("max_significance", 0.05)
    resources:
        mem_mb = 32000, runtime = 180
    threads: 4
    shell:
        """
        set -e
        {params.rbin} {R_MAASLIN2} \
            --vst {params.vst} --chemical {params.chemical} \
            --samplesheet {input.samplesheet} --out_dir {params.out_dir} \
            --elements "{params.elements}" \
            --random_effects "{params.random}" \
            --fixed_covariates "{params.fixcov}" \
            --transform {params.transform} \
            --top_var_genes {params.topgenes} \
            --max_significance {params.maxsig} \
            --cores {threads}
        """


# =============================================================================
# Rule — High-confidence bioindicators (MaAsLin2 ∩ ppcor per element)
# =============================================================================
rule bioindicator_venn:
    input:
        maaslin2     = f"{maaslin_dir}/maaslin2_significant.tsv",
        partial_corr = f"{bioindic_dir}/partial_correlations_significant.csv"
    output:
        summary = f"{biovenn_dir}/high_confidence_bioindicators_summary.csv"
    params:
        out_dir  = lambda w: os.path.abspath(biovenn_dir),
        elements = BIOVENN_ELEMENTS
    resources:
        mem_mb = 4000, runtime = 15
    threads: 1
    shell:
        """
        set -e
        Rscript {R_BIOVENN} \
            --maaslin2     {input.maaslin2} \
            --partial_corr {input.partial_corr} \
            --elements     "{params.elements}" \
            --out_dir      {params.out_dir}
        """

# =============================================================================
# Rule — Pipeline summary
# =============================================================================
# Consolidates key pipeline metrics (gene catalogue size, N50, mapping
# rate, DEG counts, Kaiju classification, dominant phyla, significant
# envfit elements, significant expression-concentration correlations)
# into one human-readable file. Every number is read from an existing
# pipeline output — this rule recomputes nothing, it only aggregates and presents.
rule pipeline_summary:
    input:
        f"{assembly_qc_dir}/assembly_qc_summary.json",
        f"{kaiju_dir}/kaiju_summary.json",
        f"{DEG_RESULTS}/all_contrasts_summary.csv",
        f"{multiqc_dir}/final/multiqc_report.html",
        *([f"{chemical_ord_dir}/envfit_results.csv"] if CHEM_ENABLED else []),
        *([f"{bioindic_dir}/partial_correlations_significant.csv"] if BIO_ENABLED else []),
        *([f"{maaslin_dir}/maaslin2_significant.tsv"] if ML_ENABLED else []),
        *([f"{chemical_ord_dir}/element_compound_correlations.csv"]
          if (RES_ENABLED and CHEM_ENABLED) else [])
    output:
        summary = f"{output_dir}/PIPELINE_SUMMARY.md"
    params:
        results_dir = lambda w: os.path.abspath(output_dir),
        group_col   = GROUP_COL
    resources:
        mem_mb = 4000, runtime = 15
    threads: 1
    shell:
        """
        set -e
        python {PY_PIPELINE_SUMMARY} \
            --results_dir {params.results_dir} \
            --group_col   {params.group_col} \
            --out_file    {output.summary}
        """
