#!/usr/bin/env python3
"""compare_runs.py — compare two MeToReS results folders, metric by metric.

Two uses:

  regression     (default) "before" is a run of the previous code, "after" a run
                 of the new code. Every metric that changed is listed with the
                 fixes expected to change it; a change that no fix explains is
                 marked UNEXPECTED. Outputs that only the new code produces are
                 checked for presence.

  reproducibility  two runs of the same code. With --tolerance exact (level
                 1a: re-run from fixed quant.sf files) any difference is
                 UNEXPECTED, except assembly-derived numbers (multi-threaded
                 MEGAHIT), and every output table from gene selection onward
                 must be identical byte for byte (absolute paths of the two
                 results folders are ignored). With --tolerance 1b (re-run from the fixed
                 catalogue, Salmon included) small differences from Salmon's
                 multi-threaded EM are accepted within the limits in TOL_1B.
                 Both also run table-level checks: per-contrast DESeq2 fold
                 change correlation and DEG-set overlap, MaAsLin2 coefficient
                 correlation and significant set, envfit r² correlation and
                 significant set.

                 With --tolerance level2 (end to end, assembly included) gene
                 IDs differ between runs, so genes are compared through their
                 annotations. Criteria, fixed before the level-2 run:
                   - gene catalogue size within ± 0.5 %;
                   - Salmon mapping rate within ± 1 percentage point per sample;
                   - KO-level fold changes: Salmon counts of the tested genes
                     summed per KO, DESeq2 with the run's design (--formula)
                     and poscounts size factors, fold changes shrunk with
                     ashr, per contrast; KOs with ≥ --ko_min_count normalised
                     counts in both runs: Pearson r ≥ 0.98 and Spearman
                     ρ ≥ 0.95 (ρ guards against a few presence/absence KOs
                     carrying r). The log2 ratio of KO group means is
                     reported without a criterion;
                   - DEG KO sets under a fixed testing rule (BH over every gene
                     with a p-value, i.e. no independent filtering; same padj
                     and fold-change cutoffs as the pipeline), as (KO,
                     direction) pairs:
                       * no robust flips, every contrast: a pair called in one
                         run has adjusted P ≤ 0.10 in the other (best gene of
                         that KO, same direction);
                       * Jaccard ≥ 0.90, only where the larger set has ≥ 50
                         pairs (on smaller sets one borderline KO decides);
                   - top 20 KOs, AMR gene families and phyla (expression-
                     weighted, catalogue Kaiju taxonomy) by abundance: ≥ 90 %
                     of each run's top 20 within the other run's top 25, so
                     near-ties around rank 20 do not fail the check (AMR
                     families only when the resistome module ran);
                   - safeguards: S2 flags identical (low-complexity samples,
                     library-type mismatches, contrasts mixing types); S4
                     status (flagged or not) the same in all DESeq2 hit sets
                     but at most one in six, rounded down (10 of 12 with six
                     contrasts; 2 of 2 with one contrast).
                 Everything else, DEG counts and borderline calls included, is
                 reported without a criterion.

Usage:
  python validation/compare_runs.py results_before_fixes results \\
      --group_col country --out validation/regression_report.md \\
      [--size_factors --samplesheet samplesheets/samplesheet.csv]
  python validation/compare_runs.py results_run1 results_run2 \\
      --mode reproducibility [--tolerance 1b] --out validation/reproducibility_report.md
  python validation/compare_runs.py results_run1 results_run2 \
      --mode reproducibility --tolerance level2 \
      --samplesheet samplesheets/samplesheet.csv --deg_group_col deseq2_group \
      --out validation/reproducibility_level2.md

ForHot demonstration (v1.1.0), group column "group", design ~ gradient + group:
  python validation/compare_runs.py results_forhot results_forhot_1a \
      --mode reproducibility --tolerance exact --group_col group \
      --out validation/forhot_level1a.md
  python validation/compare_runs.py results_forhot results_forhot_1b \
      --mode reproducibility --tolerance 1b --group_col group \
      --samplesheet examples/forhot/samplesheet_forhot.csv --out validation/forhot_level1b.md
  python validation/compare_runs.py results_forhot results_forhot_external \
      --mode reproducibility --tolerance level2 --group_col group \
      --samplesheet examples/forhot/samplesheet_forhot.csv --formula "~ gradient + group" \
      --out validation/forhot_level2.md
"""
import argparse
import glob
import hashlib
import json
import math
import os
import shutil
import subprocess
import tempfile

import numpy as np
import pandas as pd

# Which fixes can legitimately change each metric (regression mode).
EXPECTED = {
    "catalogue":   ["F3 min ORF length", "assembly re-run"],
    "orfs":        ["assembly re-run (MEGAHIT is not fully deterministic)"],
    "mapping":     ["B1 Salmon -l ISR", "F3 catalogue"],
    "kaiju":       ["F3 catalogue"],
    "deg":         ["F6 poscounts", "B1 Salmon", "F3 catalogue"],
    "dbrda":       ["F1 adjusted R² / aliased full model", "F9 VST blind", "B1", "F3", "F6"],
    "envfit":      ["envfit per-element fit", "F9 VST blind", "B1", "F3", "F6"],
    "assoc":       ["S1 prevalence filter", "F9 VST blind", "F6 poscounts", "B1", "F3"],
    "resistance":  ["F2 CARD model type / identity floor", "F3 catalogue"],
    "rescorr":     ["F2", "F3", "F8 shared detection floor"],
}
NONDETERMINISTIC = {"catalogue", "orfs"}   # reproducibility mode


# Level-1b limits (re-run from a fixed catalogue). Metrics not listed must be
# identical: the catalogue, its Kaiju classification and the CARD/BacMet calls
# do not depend on Salmon.
def _abs(lim):
    return lambda b, a: abs(a - b) <= lim, f"± {lim}"


def _rel(frac):
    return (lambda b, a: abs(a - b) <= max(1, frac * abs(b))), f"± {frac:.0%}"


def _side(alpha):
    return (lambda b, a: (b < alpha) == (a < alpha)), f"same side of {alpha}"


TOL_1B = [  # (metric-key prefix, key suffix or None, rule)
    ("map:", None, _abs(0.1)),                 # percentage points
    ("deg:pct_down", None, _abs(1.0)),
    ("deg:", None, _rel(0.01)),                # DEG / up / down / genes tested
    ("dbrda:", ":pct", _abs(0.5)),
    ("dbrda:", ":adj", _abs(0.005)),
    ("dbrda:", ":p", _side(0.05)),
    ("envfit", None, _abs(1)),
    ("ml", None, _abs(1)),
    ("pp", None, _abs(1)),
    ("res:aro", None, _abs(1)),
    ("rescorr:", None, _abs(1)),
    ("geneset:", ":pct", _abs(0.5)),           # share of counts with protein evidence
    ("geneset:", None, _rel(0.01)),            # filtered / annotated gene numbers
]

# Level-2 limits (end to end). Only the metrics below have a criterion; every
# other metric is reported without one (see REPORT_ONLY_L2).
def _relpct(frac):
    return (lambda b, a: abs(a - b) <= frac * abs(b)), f"± {frac:.1%}"


def _equal():
    return (lambda b, a: b == a), "identical"


TOL_L2 = [
    ("catalogue", None, _relpct(0.005)),
    ("map:", None, _abs(1.0)),                 # percentage points
]
L2_EXACT = ("lib:",)                           # library type must not change
L2_KO_FC_R, L2_KO_FC_RHO = 0.98, 0.95
L2_KO_JACCARD, L2_KO_JACCARD_MIN_SET = 0.90, 50
L2_FLIP_PADJ = 0.10                            # "robust flip": called in one run, padj > this in the other
L2_TOP_N, L2_TOP_WINDOW, L2_TOP_OVERLAP = 20, 25, 0.90
L2_S4_MAX_DISAGREE_FRAC = 1 / 6              # S4 sets allowed to differ, rounded down

# Reported without a criterion at 1b and level 2: these describe how far
# calls can move, not whether they reproduce.
REPORT_ONLY = ("degrep:",)

# Table-level limits: (min correlation, min DEG-set Jaccard, max symmetric
# difference of small significant sets)
TABLE_LIMITS = {"exact": (1 - 1e-12, 1.0, 0), "1b": (0.999, 0.99, 1)}


def tolerance_rule(key, level="1b"):
    rules = TOL_L2 if level == "level2" else TOL_1B
    if level == "level2" and key.startswith(L2_EXACT):
        return (lambda b, a: False), "identical"
    for prefix, suffix, rule in rules:
        if key.startswith(prefix) and (suffix is None or key.endswith(suffix)):
            return rule
    return None


def read_json(p):
    try:
        with open(p) as fh:
            return json.load(fh)
    except Exception:
        return None


def read_csv(p, **kw):
    try:
        return pd.read_csv(p, **kw)
    except Exception:
        return None


def count_fasta(p):
    if not os.path.exists(p):
        return None
    n = 0
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def one(globpat):
    hits = sorted(glob.glob(globpat))
    return hits[0] if hits else None


def size_factors(R):
    """Effective DESeq2 size factor per sample: median over genes with a
    non-zero count of raw count / normalised count. Raw counts are Salmon's
    NumReads (what tximport hands to DESeq2); normalised counts are DESeq2's."""
    nf = one(os.path.join(R, "deg", "*_deg_results", "deseq2_normalized_counts.csv"))
    if not nf:
        return {}
    norm = pd.read_csv(nf, index_col=0)
    norm.index = norm.index.astype(str)
    ss_like = {os.path.basename(os.path.dirname(f)).replace("_quant", ""): f
               for f in glob.glob(os.path.join(R, "quant", "*_quant", "quant.sf"))}
    ssheet = SAMPLESHEET
    name_to_sra = {}
    if ssheet:
        s = pd.read_csv(ssheet)
        name_to_sra = dict(zip(s["sample"], s["sra"].astype(str)))
    out = {}
    for sample in norm.columns:
        qf = ss_like.get(name_to_sra.get(sample, sample))
        if not qf:
            continue
        q = pd.read_csv(qf, sep="\t", usecols=["Name", "NumReads"], index_col=0)["NumReads"]
        raw = q.reindex(norm.index).fillna(0).round()
        ok = (raw > 0) & (norm[sample] > 0)
        if ok.any():
            out[sample] = float((raw[ok] / norm.loc[ok, sample]).median())
    return out


def collect(R, gc):
    """Return {metric_id: (category, label, value)} for one results folder."""
    m = {}

    def put(key, cat, label, val):
        m[key] = (cat, label, val)

    put("catalogue", "catalogue", "Gene catalogue size",
        count_fasta(os.path.join(R, "gene_catalog", "all_longest_orfs_cds_rmdup_id.fasta")))
    for f in sorted(glob.glob(os.path.join(R, "assembly", "*", "orfs.fasta"))):
        s = os.path.basename(os.path.dirname(f))
        put(f"orfs:{s}", "orfs", f"ORFs assembled · {s}", count_fasta(f))

    for f in sorted(glob.glob(os.path.join(R, "quant", "*_quant", "aux_info", "meta_info.json"))):
        s = os.path.basename(os.path.dirname(os.path.dirname(f))).replace("_quant", "")
        d = read_json(f) or {}
        put(f"map:{s}", "mapping", f"Salmon % mapped · {s}", d.get("percent_mapped"))
        put(f"lib:{s}", "mapping", f"Salmon library type · {s}",
            d.get("library_types", [None])[0] if isinstance(d.get("library_types"), list) else d.get("library_types"))

    gs = read_json(one(os.path.join(R, "deg", "*_deg_results", "gene_set_summary.json")) or "") or {}
    if gs:
        put("geneset:n_filtered", "deg", "Filtered genes", gs.get("n_filtered"))
        put("geneset:n_protein", "deg", "Filtered genes with protein-level evidence",
            gs.get("n_with_protein_evidence"))
        put("geneset:counts:pct", "deg", "% of filtered-gene counts with protein-level evidence",
            gs.get("pct_counts_with_protein_evidence"))

    k = read_json(os.path.join(R, "kaiju", "kaiju_summary.json")) or {}
    put("kaiju", "kaiju", "Kaiju: % catalogue ORFs classified", k.get("classification_rate"))

    deg = read_csv(one(os.path.join(R, "deg", "*_deg_results", "all_contrasts_summary.csv")) or "")
    if deg is not None and "Contrast" in deg.columns:
        for _, r in deg.iterrows():
            c = r["Contrast"]
            for col, lab in (("Total_DEGs", "DEGs"), ("Upregulated", "up"),
                             ("Downregulated", "down"), ("Genes_tested", "genes tested")):
                if col in deg.columns:
                    put(f"deg:{c}:{col}", "deg", f"DESeq2 {lab} · {c}", r[col])
            for col, lab in (("Borderline", "borderline calls"),
                             ("Filter_baseMean", "independent-filtering cutoff"),
                             ("PresenceAbsence", "% presence/absence DEGs")):
                if col in deg.columns:
                    put(f"degrep:{c}:{col}", "deg", f"DESeq2 {lab} · {c}", r[col])

    if deg is not None and {"Upregulated", "Downregulated"} <= set(deg.columns):
        up, down = deg["Upregulated"].sum(), deg["Downregulated"].sum()
        if up + down:
            put("deg:pct_down", "deg", "DESeq2 % of DEGs down-regulated (all contrasts)",
                round(100 * down / (up + down), 1))
    if SIZE_FACTORS:
        sf = size_factors(R)
        if sf:
            put("deg:sf_min", "deg", "DESeq2 effective size factor, min", round(min(sf.values()), 3))
            put("deg:sf_max", "deg", "DESeq2 effective size factor, max", round(max(sf.values()), 3))
            put("deg:sf_ratio", "deg", "DESeq2 effective size factor, max/min",
                round(max(sf.values()) / min(sf.values()), 1))

    db = read_csv(os.path.join(R, f"constrained_ordination_{gc}", "dbrda_anova_results.csv"))
    if db is not None and "Model" in db.columns:
        for _, r in db.iterrows():
            put(f"dbrda:{r['Model']}:pct", "dbrda", f"dbRDA % constrained · {r['Model']}", r.get("Pct_constrained"))
            put(f"dbrda:{r['Model']}:p", "dbrda", f"dbRDA p · {r['Model']}", r.get("P_value"))
            if "Adj_R2" in db.columns:
                put(f"dbrda:{r['Model']}:adj", "dbrda", f"dbRDA adjusted R² · {r['Model']}", r.get("Adj_R2"))

    ef = read_csv(os.path.join(R, f"chemical_{gc}", "envfit_results.csv"))
    if ef is not None and "pvalue" in ef.columns:
        put("envfit", "envfit", "envfit: elements with p < 0.05", int((ef["pvalue"] < 0.05).sum()))

    ml = read_csv(os.path.join(R, f"maaslin2_{gc}", "maaslin2_significant.tsv"), sep="\t")
    if ml is not None:
        put("ml", "assoc", "MaAsLin2: significant associations", len(ml))
        if "feature" in ml.columns:
            put("ml_genes", "assoc", "MaAsLin2: genes with ≥1 association", ml["feature"].nunique())
    pp = read_csv(os.path.join(R, f"bioindicators_{gc}", "partial_correlations_significant.csv"))
    if pp is not None:
        put("pp", "assoc", "Partial correlation: significant pairs", len(pp))

    ann = read_csv(os.path.join(R, f"resistance_{gc}", "resistance_annotation.tsv"), sep="\t", low_memory=False)
    if ann is not None:
        for col, lab in (("is_ARG", "ARG ORFs (CARD)"), ("is_MRG", "MRG ORFs (BacMet)")):
            if col in ann.columns:
                put(f"res:{col}", "resistance", lab,
                    int(ann[col].astype(str).str.lower().isin(["true", "1"]).sum()))
    aro = read_csv(os.path.join(R, f"resistance_{gc}", "arg_abundance_by_aro.csv"))
    if aro is not None:
        put("res:aro", "resistance", "ARO terms with expression", len(aro))
    rc = read_csv(os.path.join(R, f"chemical_{gc}", "element_compound_correlations.csv"))
    if rc is not None:
        put("rescorr:n", "rescorr", "Resistance–chemistry pairs tested", len(rc))
        if "padj" in rc.columns:
            put("rescorr:sig", "rescorr", "Resistance–chemistry pairs, padj < 0.05", int((rc["padj"] < 0.05).sum()))
    return m


NEW_OUTPUTS = [
    ("safeguards/SAFEGUARDS_SUMMARY.md", "S1–S4 safeguard summary"),
    ("safeguards/library_complexity.tsv", "S2 library complexity"),
    ("safeguards/assembly_origin_enrichment.tsv", "S4 assembly origin"),
    ("resistance_{gc}/resistome_coverage.csv", "F8 resistome coverage"),
    ("resistance_{gc}/arg_abundance_by_family.csv", "F7 ARG by gene family"),
    ("resistance_{gc}/card_annotations_other_models.tsv", "F2 other CARD model types"),
    ("qc/multiqc/multiqc_config.yaml", "B4 MultiQC config file"),
    ("qc/multiqc/final/multiqc_report_data/multiqc_sortmerna.txt", "B4 SortMeRNA in MultiQC"),
    ("qc/multiqc/final/multiqc_report_data/multiqc_trimmomatic.txt", "F4 Trimmomatic in MultiQC"),
    ("kaiju/contamination_ids.csv", "F5 contamination list (fixed name)"),
]


def same(a, b):
    if a is None or b is None:
        return a is None and b is None
    try:
        fa, fb = float(a), float(b)
        if math.isnan(fa) or math.isnan(fb):
            return math.isnan(fa) and math.isnan(fb)
        return abs(fa - fb) <= 1e-9 * max(1.0, abs(fa), abs(fb))
    except (TypeError, ValueError):
        return str(a) == str(b)


def fmt(v):
    if v is None:
        return "—"
    if isinstance(v, float):
        return f"{v:,.4g}"
    if isinstance(v, int) or (isinstance(v, str) is False and float(v).is_integer()):
        try:
            return f"{int(v):,}"
        except (TypeError, ValueError):
            pass
    return str(v)


def delta(a, b):
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        return ""
    d = fb - fa
    pct = f" ({100 * d / fa:+.1f}%)" if fa else ""
    return f"{d:+,.4g}{pct}"


def _corr(x, y):
    x, y = pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")
    ok = x.notna() & y.notna()
    if ok.sum() < 3:
        return None
    if x[ok].equals(y[ok]):
        return 1.0
    return float(x[ok].corr(y[ok]))


def _spearman(x, y):
    x, y = pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")
    ok = x.notna() & y.notna()
    if ok.sum() < 3:
        return None
    return float(x[ok].rank().corr(y[ok].rank()))


# Folders written from gene selection onward (what a level-1a re-run rewrites).
STAGE_DIRS = ("deg/{gc}_deg_results", "downstream_{gc}", "constrained_ordination_{gc}", "chemical",
              "chemical_{gc}", "maaslin2_{gc}", "bioindicators_{gc}", "resistance_{gc}", "safeguards")
TABLE_EXT = (".csv", ".tsv", ".txt", ".json")


def _digest(path, roots):
    with open(path, "rb") as fh:
        data = fh.read()
    if path.endswith(".json"):                     # summaries may record the run folder
        for r in roots:
            data = data.replace(r.encode(), b"<RESULTS>")
    return hashlib.md5(data).hexdigest()


def identity_check(Rb, Ra, gc):
    """Level 1a: tables written from gene selection onward, byte for byte."""
    files = []
    for d in STAGE_DIRS:
        base = os.path.join(Rb, d.format(gc=gc))
        for root, _, names in os.walk(base):
            files += [os.path.relpath(os.path.join(root, n), Rb) for n in names if n.endswith(TABLE_EXT)]
    roots = [os.path.abspath(Rb), os.path.abspath(Ra)]
    same_, differ, missing = 0, [], []
    for rel in sorted(files):
        fa = os.path.join(Ra, rel)
        if not os.path.exists(fa):
            missing.append(rel)
        elif _digest(os.path.join(Rb, rel), roots) == _digest(fa, roots):
            same_ += 1
        else:
            differ.append(rel)
    detail = f"{same_} of {len(files)} identical"
    if differ:
        detail += f"; differ: {', '.join(differ[:5])}" + (" …" if len(differ) > 5 else "")
    if missing:
        detail += f"; missing in second run: {', '.join(missing[:5])}" + (" …" if len(missing) > 5 else "")
    return ("Output tables identical byte for byte (gene selection onward)", detail,
            "all identical", same_ == len(files) and len(files) > 0)


def table_checks(Rb, Ra, gc, level):
    """Return [(check, value, criterion, passed)] comparing tables of two runs."""
    rmin, jmin, dmax = TABLE_LIMITS[level]
    out = [identity_check(Rb, Ra, gc)] if level == "exact" else []

    def add_corr(label, r):
        if r is not None:
            out.append((label, f"r = {r:.6f}", f"r ≥ {rmin:g}" if level != "exact" else "r = 1",
                        r >= rmin))

    def add_set(label, sb, sa):
        d = len(sb ^ sa)
        out.append((label, f"{len(sb)} vs {len(sa)}, {d} differ", f"≤ {dmax} differ", d <= dmax))

    for fb in sorted(glob.glob(os.path.join(Rb, "deg", "*_deg_results", "*", "deseq2_results_all.csv"))):
        rel = os.path.relpath(fb, Rb)
        c = os.path.basename(os.path.dirname(fb))
        tb, ta = read_csv(fb), read_csv(os.path.join(Ra, rel))
        if tb is None or ta is None or "GeneID" not in tb.columns:
            continue
        m = tb.set_index("GeneID")[["log2FoldChange"]].join(
            ta.set_index("GeneID")[["log2FoldChange"]], lsuffix="_b", rsuffix="_a", how="inner")
        add_corr(f"DESeq2 log2FC correlation · {c}", _corr(m["log2FoldChange_b"], m["log2FoldChange_a"]))
        db_ = read_csv(os.path.join(os.path.dirname(fb), "differential_genes.csv"))
        da_ = read_csv(os.path.join(Ra, os.path.dirname(rel), "differential_genes.csv"))
        if db_ is not None and da_ is not None:
            sb, sa = set(db_["GeneID"]), set(da_["GeneID"])
            j = len(sb & sa) / len(sb | sa) if (sb | sa) else 1.0
            out.append((f"DEG set overlap (Jaccard) · {c}", f"{j:.4f} ({len(sb)} vs {len(sa)})",
                        f"≥ {jmin:g}", j >= jmin))
        if level == "1b" and {"pvalue", "log2FoldChange"} <= set(tb.columns) and "pvalue" in ta.columns:
            def fixed(t):
                t = t[t["pvalue"].notna()]
                q = bh(t["pvalue"].to_numpy())
                return set(t["GeneID"].astype(str)[(q < FIXED_ALPHA) & (t["log2FoldChange"].abs() > FIXED_LFC)])
            fb_, fa_ = fixed(tb), fixed(ta)
            jf = len(fb_ & fa_) / len(fb_ | fa_) if (fb_ | fa_) else 1.0
            out.append((f"DEG set overlap, fixed rule (Jaccard) · {c}",
                        f"{jf:.4f} ({len(fb_)} vs {len(fa_)}; {len(fb_ ^ fa_)} calls differ, "
                        f"{100 * len(fb_ ^ fa_) / max(1, len(fb_ | fa_)):.1f}%)", "report only", None))

    mb = read_csv(os.path.join(Rb, f"maaslin2_{gc}", "maaslin2_all_elements.tsv"), sep="\t")
    ma = read_csv(os.path.join(Ra, f"maaslin2_{gc}", "maaslin2_all_elements.tsv"), sep="\t")
    if mb is not None and ma is not None and {"feature", "metadata", "coef"} <= set(mb.columns):
        k = ["feature", "metadata"]
        m = mb.set_index(k)[["coef"]].join(ma.set_index(k)[["coef"]], lsuffix="_b", rsuffix="_a", how="inner")
        add_corr("MaAsLin2 coefficient correlation", _corr(m["coef_b"], m["coef_a"]))
    sb_ = read_csv(os.path.join(Rb, f"maaslin2_{gc}", "maaslin2_significant.tsv"), sep="\t")
    sa_ = read_csv(os.path.join(Ra, f"maaslin2_{gc}", "maaslin2_significant.tsv"), sep="\t")
    if sb_ is not None and sa_ is not None and {"feature", "metadata"} <= set(sb_.columns):
        add_set("MaAsLin2 significant set",
                set(zip(sb_["feature"], sb_["metadata"])), set(zip(sa_["feature"], sa_["metadata"])))

    eb = read_csv(os.path.join(Rb, f"chemical_{gc}", "envfit_results.csv"))
    ea = read_csv(os.path.join(Ra, f"chemical_{gc}", "envfit_results.csv"))
    if eb is not None and ea is not None and {"Element", "r2", "pvalue"} <= set(eb.columns):
        m = eb.set_index("Element")[["r2"]].join(ea.set_index("Element")[["r2"]],
                                                 lsuffix="_b", rsuffix="_a", how="inner")
        add_corr("envfit r² correlation", _corr(m["r2_b"], m["r2_a"]))
        add_set("envfit elements with p < 0.05",
                set(eb.loc[eb["pvalue"] < 0.05, "Element"]), set(ea.loc[ea["pvalue"] < 0.05, "Element"]))
    return out


# ── Level 2: annotation-level comparison ──────────────────────────────────
def bh(p):
    """Benjamini–Hochberg adjusted p-values (as R's p.adjust 'BH')."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)[::-1]
    q = np.minimum.accumulate(p[o] * n / np.arange(n, 0, -1))
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def _read_emapper_ko(path, genes):
    """gene -> list of KOs (KEGG_ko column, 'ko:' prefix removed), for `genes`."""
    if not path or not os.path.exists(path):
        return None
    skip = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#query"):
                break
            skip += 1
    parts = []
    for ch in pd.read_csv(path, sep="\t", skiprows=skip, usecols=["#query", "KEGG_ko"],
                          dtype=str, chunksize=1_000_000):
        ch = ch[ch["#query"].isin(genes) & ch["KEGG_ko"].notna() & (ch["KEGG_ko"] != "-")]
        parts.append(ch)
    if not parts:
        return pd.DataFrame(columns=["gene", "KO"])
    d = pd.concat(parts).rename(columns={"#query": "gene", "KEGG_ko": "KO"})
    d["KO"] = d["KO"].str.replace("ko:", "", regex=False).str.split(",")
    return d.explode("KO").drop_duplicates()


def _read_gene_phylum(path, genes):
    if not path or not os.path.exists(path):
        return None
    parts = []
    for ch in pd.read_csv(path, sep="\t", usecols=["gene_id", "phylum"], dtype=str,
                          chunksize=2_000_000):
        parts.append(ch[ch["gene_id"].isin(genes)])
    d = pd.concat(parts) if parts else pd.DataFrame(columns=["gene_id", "phylum"])
    d = d[d["phylum"].notna() & ~d["phylum"].str.upper().isin(["NA", "NAN", "UNCLASSIFIED", ""])]
    return d.set_index("gene_id")["phylum"]


R_KO_DESEQ = r"""
suppressPackageStartupMessages(library(DESeq2))
a <- commandArgs(trailingOnly = TRUE)
cts <- read.csv(a[1], row.names = 1, check.names = FALSE)
cd  <- read.csv(a[2], row.names = 1, check.names = FALSE, colClasses = "character")
design <- as.formula(a[3]); gcol <- a[4]; contrasts <- strsplit(a[5], ",", fixed = TRUE)[[1]]
cd <- cd[colnames(cts), , drop = FALSE]
for (v in all.vars(design)) cd[[v]] <- factor(cd[[v]])
cts <- round(as.matrix(cts)); cts <- cts[rowSums(cts) > 0, , drop = FALSE]
dds <- DESeqDataSetFromMatrix(cts, cd, design)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, quiet = TRUE)
out <- do.call(rbind, lapply(contrasts, function(cn) {
  lv <- strsplit(cn, "_vs_", fixed = TRUE)[[1]]
  r <- results(dds, contrast = c(gcol, lv[1], lv[2]))
  s <- lfcShrink(dds, contrast = c(gcol, lv[1], lv[2]), type = "ashr", res = r, quiet = TRUE)
  data.frame(contrast = cn, KO = rownames(s), baseMean = s$baseMean,
             lfc = s$log2FoldChange, lfc_mle = r$log2FoldChange)
}))
write.csv(out, a[6], row.names = FALSE)
"""


def ko_shrunken_fc(R, kc_raw, samplesheet, formula, group_col, contrasts):
    """Shrunken KO-level log2FC per contrast: DESeq2 (poscounts) + ashr on the
    KO-summed Salmon counts of one run. Returns {contrast: Series} or an error string."""
    rs = shutil.which("Rscript")
    if rs is None:
        return "Rscript not found"
    s = pd.read_csv(samplesheet).astype(str).set_index("sample")
    with tempfile.TemporaryDirectory() as tmp:
        cf, mf, rf, of = (os.path.join(tmp, n) for n in ("counts.csv", "coldata.csv", "ko.R", "out.csv"))
        kc_raw.to_csv(cf)
        s.reindex(kc_raw.columns).to_csv(mf)
        with open(rf, "w") as fh:
            fh.write(R_KO_DESEQ)
        p = subprocess.run([rs, rf, cf, mf, formula, group_col, ",".join(contrasts), of],
                           capture_output=True, text=True)
        if p.returncode != 0 or not os.path.exists(of):
            return "DESeq2 on KO counts failed: " + (p.stderr.strip().splitlines() or ["?"])[-1]
        res = pd.read_csv(of)
    return {c: g.set_index("KO")["lfc"] for c, g in res.groupby("contrast")}


def level2_tables(R, gc, groups, alpha, lfc, samplesheet=None, formula=None, group_col=None):
    """Annotation-level tables for one run."""
    out = {}
    nf = one(os.path.join(R, "deg", "*_deg_results", "deseq2_normalized_counts.csv"))
    if not nf:
        return out
    norm = pd.read_csv(nf, index_col=0)
    norm.index = norm.index.astype(str)
    genes = set(norm.index)
    deg_dir = os.path.dirname(nf)
    ko = _read_emapper_ko(one(os.path.join(R, f"emapper_{gc}", "*.emapper.annotations")), genes)
    if ko is not None and len(ko):
        kc = norm.loc[norm.index.intersection(ko["gene"].unique())]
        kc = kc.reindex(ko["gene"].values).set_axis(ko["KO"].values).groupby(level=0).sum()
        out["ko_counts"] = kc
        # Raw Salmon counts of the tested genes, summed per KO (for DESeq2 on KOs)
        if samplesheet:
            s = pd.read_csv(samplesheet)
            to_sra = dict(zip(s["sample"].astype(str), s["sra"].astype(str)))
            raw = {}
            for sample in norm.columns:
                qf = os.path.join(R, "quant", f"{to_sra.get(sample, sample)}_quant", "quant.sf")
                if os.path.exists(qf):
                    q = pd.read_csv(qf, sep="\t", usecols=["Name", "NumReads"], index_col=0)["NumReads"]
                    raw[sample] = q.reindex(ko["gene"].values).fillna(0).to_numpy()
            if len(raw) == len(norm.columns):
                out["ko_raw"] = pd.DataFrame(raw, index=ko["KO"].values).groupby(level=0).sum()
        out["ko_map"] = ko.groupby("gene")["KO"].apply(list)
        # KO-level fold change per contrast and fixed-rule DEG KO sets
        out["ko_fc"], out["ko_deg"], out["ko_padj"] = {}, {}, {}
        for f in sorted(glob.glob(os.path.join(deg_dir, "*", "deseq2_results_all.csv"))):
            c = os.path.basename(os.path.dirname(f))
            if "_vs_" not in c:
                continue
            lvl, ref = c.split("_vs_", 1)
            sa = [s for s in kc.columns if groups.get(s) == lvl]
            sb = [s for s in kc.columns if groups.get(s) == ref]
            if sa and sb:
                out["ko_fc"][c] = np.log2((kc[sa].mean(axis=1) + 1) / (kc[sb].mean(axis=1) + 1))
            res = pd.read_csv(f, usecols=lambda x: x in ("GeneID", "log2FoldChange", "pvalue"))
            res = res[res["pvalue"].notna()].copy()
            if res.empty:
                continue
            res["padj_fixed"] = bh(res["pvalue"].to_numpy())
            best = {}
            for g, l2, q in zip(res["GeneID"].astype(str), res["log2FoldChange"], res["padj_fixed"]):
                d = "up" if l2 > 0 else "down"
                for k in out["ko_map"].get(g, []):
                    if q < best.get((k, d), 1.0):
                        best[(k, d)] = q
            out["ko_padj"][c] = best
            sig = res[(res["padj_fixed"] < alpha) & (res["log2FoldChange"].abs() > lfc)]
            pairs = set()
            for g, l2 in zip(sig["GeneID"].astype(str), sig["log2FoldChange"]):
                for k in out["ko_map"].get(g, []):
                    pairs.add((k, "up" if l2 > 0 else "down"))
            out["ko_deg"][c] = pairs
    fam = read_csv(os.path.join(R, f"resistance_{gc}", "arg_abundance_by_family.csv"))
    if fam is not None and "AMR_Gene_Family" in fam.columns:
        fam = fam[fam["AMR_Gene_Family"].astype(str).str.lower() != "unclassified"]
        out["family"] = fam.set_index("AMR_Gene_Family").select_dtypes("number").sum(axis=1)
    ph = _read_gene_phylum(os.path.join(R, "kaiju", "kaiju_gene_taxonomy.tsv"), genes)
    if ph is not None and len(ph):
        tot = norm.loc[norm.index.intersection(ph.index)].sum(axis=1)
        out["phylum"] = tot.groupby(ph.reindex(tot.index).values).sum()
    if "ko_counts" in out:
        out["ko_total"] = out["ko_counts"].sum(axis=1)
    if "ko_raw" in out and formula and group_col and out.get("ko_fc"):
        out["ko_shrunk"] = ko_shrunken_fc(R, out["ko_raw"], samplesheet, formula, group_col,
                                          sorted(out["ko_fc"]))
    return out


def safeguard_checks(Rb, Ra):
    rows = []
    jb = read_json(os.path.join(Rb, "safeguards", "library_complexity.json"))
    ja = read_json(os.path.join(Ra, "safeguards", "library_complexity.json"))
    if jb is None or ja is None:
        rows.append(("S2 flags identical", "library_complexity.json missing in one run", "identical", False))
    else:
        def s2(j):
            return {"low-complexity": {f["sample"] for f in j.get("flagged", [])},
                    "type mismatch": {m["sample"] for m in j.get("libtype_mismatch", [])},
                    "mixed contrasts": {m["contrast"] for m in j.get("contrasts_mixing_libtypes", [])}}
        b, a = s2(jb), s2(ja)
        diff = {k: b[k] ^ a[k] for k in b if b[k] ^ a[k]}
        rows.append(("S2 flags identical",
                     "; ".join(f"{k}: {len(b[k])}" for k in b) if not diff else
                     "differ: " + "; ".join(f"{k} {sorted(v)}" for k, v in diff.items()),
                     "identical", not diff))
    tb = read_csv(os.path.join(Rb, "safeguards", "assembly_origin_enrichment.tsv"), sep="\t")
    ta = read_csv(os.path.join(Ra, "safeguards", "assembly_origin_enrichment.tsv"), sep="\t")
    if tb is None or ta is None or "enriched" not in tb.columns:
        rows.append(("S4 status per DESeq2 hit set", "assembly_origin_enrichment.tsv missing in one run",
                     "≤ 1 in 6 sets differ", False))
    else:
        def s4(t):
            t = t[t["analysis"] == "DESeq2"]
            return t.groupby("set")["enriched"].apply(lambda x: bool(x.astype(str).str.lower().eq("true").any()))
        b, a = s4(tb), s4(ta)
        sets = sorted(set(b.index) | set(a.index))
        dis = [st for st in sets if b.get(st) != a.get(st)]
        allowed = int(len(sets) * L2_S4_MAX_DISAGREE_FRAC)
        rows.append(("S4 status per DESeq2 hit set",
                     f"{len(sets) - len(dis)}/{len(sets)} sets agree"
                     + (f" (differ: {', '.join(dis)})" if dis else ""),
                     f"≤ {allowed} of {len(sets)} sets differ", len(dis) <= allowed))
    return rows


def level2_checks(Rb, Ra, gc, samplesheet, group_col, alpha, lfc, ko_min, formula=None):
    rows = []
    groups = {}
    if samplesheet:
        s = pd.read_csv(samplesheet)
        if group_col in s.columns:
            groups = dict(zip(s["sample"].astype(str), s[group_col].astype(str)))
    if not groups:
        rows.append(("KO-level checks", f"no '{group_col}' column in --samplesheet",
                     "needs --samplesheet / --deg_group_col", False))
    formula = formula or f"~ {group_col}"
    tb = level2_tables(Rb, gc, groups, alpha, lfc, samplesheet, formula, group_col)
    ta = level2_tables(Ra, gc, groups, alpha, lfc, samplesheet, formula, group_col)

    def ko_common(fb, fa):
        keep_b = tb["ko_total"].reindex(fb.index).fillna(0) >= ko_min
        keep_a = ta["ko_total"].reindex(fa.index).fillna(0) >= ko_min
        return fb.index[keep_b].intersection(fa.index[keep_a])

    crit = f"r ≥ {L2_KO_FC_R} and ρ ≥ {L2_KO_FC_RHO}"
    for c in sorted(set(tb.get("ko_fc", {})) | set(ta.get("ko_fc", {}))):
        label = f"KO log2FC correlation, shrunken (DESeq2 + ashr on KO counts) · {c}"
        sb_, sa_ = tb.get("ko_shrunk"), ta.get("ko_shrunk")
        if isinstance(sb_, str) or isinstance(sa_, str) or sb_ is None or sa_ is None:
            why = sb_ if isinstance(sb_, str) else sa_ if isinstance(sa_, str) else \
                "needs --samplesheet with an 'sra' column and the quant.sf files of both runs"
            rows.append((label, why, crit, False))
        elif c not in sb_ or c not in sa_:
            rows.append((label, "missing in one run", crit, False))
        else:
            fb, fa = sb_[c], sa_[c]
            common = ko_common(fb, fa)
            r, rho = _corr(fb.reindex(common), fa.reindex(common)), _spearman(fb.reindex(common), fa.reindex(common))
            rows.append((label, (f"r = {r:.4f}, ρ = {rho:.4f} ({len(common):,} KOs)"
                                 if r is not None and rho is not None else "too few KOs"), crit,
                         r is not None and rho is not None and r >= L2_KO_FC_R and rho >= L2_KO_FC_RHO))
        fb, fa = tb.get("ko_fc", {}).get(c), ta.get("ko_fc", {}).get(c)
        if fb is not None and fa is not None:
            common = ko_common(fb, fa)
            r, rho = _corr(fb.reindex(common), fa.reindex(common)), _spearman(fb.reindex(common), fa.reindex(common))
            if r is not None and rho is not None:
                rows.append((f"KO log2 ratio of group means, correlation · {c}",
                             f"r = {r:.4f}, ρ = {rho:.4f} ({len(common):,} KOs)", "report only", None))
    for c in sorted(set(tb.get("ko_deg", {})) | set(ta.get("ko_deg", {}))):
        sb, sa = tb.get("ko_deg", {}).get(c, set()), ta.get("ko_deg", {}).get(c, set())
        qb, qa = tb.get("ko_padj", {}).get(c, {}), ta.get("ko_padj", {}).get(c, {})
        flips = ([p for p in sb if qa.get(p, 1.0) > L2_FLIP_PADJ]
                 + [p for p in sa if qb.get(p, 1.0) > L2_FLIP_PADJ])
        rows.append((f"DEG KO calls, robust flips · {c}",
                     f"{len(flips)} of {len(sb | sa)} KO×direction calls"
                     + (f" (e.g. {', '.join(f'{k} {d}' for k, d in sorted(flips)[:3])})" if flips else ""),
                     f"0 calls with adjusted P > {L2_FLIP_PADJ} in the other run", not flips))
        j = len(sb & sa) / len(sb | sa) if (sb | sa) else 1.0
        big = max(len(sb), len(sa)) >= L2_KO_JACCARD_MIN_SET
        rows.append((f"DEG KO set overlap, fixed rule (Jaccard) · {c}",
                     f"{j:.4f} ({len(sb)} vs {len(sa)} KO×direction)",
                     f"≥ {L2_KO_JACCARD}" if big else
                     f"report only (< {L2_KO_JACCARD_MIN_SET} pairs)",
                     (j >= L2_KO_JACCARD) if big else None))
    for key, lab in (("ko_total", "KOs"), ("family", "AMR gene families"), ("phylum", "phyla")):
        vb, va = tb.get(key), ta.get(key)
        if key == "family" and vb is None and va is None:
            rows.append((f"Top {L2_TOP_N} {lab} by abundance", "resistome module not run", "—", None))
            continue
        if vb is None or va is None or not len(vb) or not len(va):
            rows.append((f"Top {L2_TOP_N} {lab} by abundance", "missing in one run",
                         f"≥ {L2_TOP_OVERLAP:.0%} within the other's top {L2_TOP_WINDOW}", False))
            continue
        k = min(L2_TOP_N, len(vb), len(va))
        rb, ra = vb.sort_values(ascending=False).index, va.sort_values(ascending=False).index
        topb, topa = set(rb[:k]), set(ra[:k])
        winb, wina = set(rb[:L2_TOP_WINDOW]), set(ra[:L2_TOP_WINDOW])
        ov = min(len(topb & wina), len(topa & winb)) / k
        rows.append((f"Top {L2_TOP_N} {lab} by abundance",
                     f"{len(topb & topa)}/{k} shared in the top {k}; "
                     f"{min(len(topb & wina), len(topa & winb))}/{k} within the other's top {L2_TOP_WINDOW}",
                     f"≥ {L2_TOP_OVERLAP:.0%} within the other's top {L2_TOP_WINDOW}",
                     ov >= L2_TOP_OVERLAP))

    # Safeguard reproduction
    rows += safeguard_checks(Rb, Ra)

    # Reported without a criterion
    eb = read_csv(os.path.join(Rb, f"chemical_{gc}", "envfit_results.csv"))
    ea = read_csv(os.path.join(Ra, f"chemical_{gc}", "envfit_results.csv"))
    if eb is not None and ea is not None and {"Element", "r2", "pvalue"} <= set(eb.columns):
        m = eb.set_index("Element")[["r2"]].join(ea.set_index("Element")[["r2"]],
                                                 lsuffix="_b", rsuffix="_a", how="inner")
        r = _corr(m["r2_b"], m["r2_a"])
        if r is not None:
            rows.append(("envfit r² correlation", f"r = {r:.4f}", "report only", None))
        pb = set(eb.loc[eb["pvalue"] < 0.05, "Element"])
        pa = set(ea.loc[ea["pvalue"] < 0.05, "Element"])
        rows.append(("envfit elements with p < 0.05", f"{len(pb)} vs {len(pa)}, {len(pb ^ pa)} differ",
                     "report only", None))
    return rows


ap = argparse.ArgumentParser()
ap.add_argument("before")
ap.add_argument("after")
ap.add_argument("--group_col", default="country")
ap.add_argument("--mode", choices=["regression", "reproducibility"], default="regression")
ap.add_argument("--tolerance", choices=["exact", "1b", "level2"], default="exact",
                help="reproducibility mode: exact (level 1a), the level-1b limits, "
                     "or the level-2 (end-to-end) criteria")
ap.add_argument("--deg_group_col", default=None,
                help="level2: samplesheet column with the DESeq2 groups (default: --group_col)")
ap.add_argument("--alpha", type=float, default=0.05,
                help="level2: adjusted-P cutoff of the fixed testing rule (pipeline 'pvalue')")
ap.add_argument("--fold_change", type=float, default=2.0,
                help="level2: raw fold-change cutoff of the fixed testing rule (pipeline 'fold_change')")
ap.add_argument("--formula", default=None,
                help="level2: DESeq2 design of the runs (deseq2.formula), e.g. '~ gradient + group' "
                     "[default: '~ <deg_group_col>']")
ap.add_argument("--ko_min_count", type=float, default=10,
                help="level2: KOs with fewer total normalised counts are left out of the "
                     "fold-change correlation")
ap.add_argument("--out", default="comparison_report.md")
ap.add_argument("--size_factors", action="store_true",
                help="also compute effective DESeq2 size factors (reads every quant.sf; slow)")
ap.add_argument("--samplesheet", default=None,
                help="samplesheet mapping sample names to quant folders (needed with --size_factors)")
a = ap.parse_args()
SIZE_FACTORS = a.size_factors
SAMPLESHEET = a.samplesheet
FIXED_ALPHA, FIXED_LFC = a.alpha, math.log2(a.fold_change)   # fixed testing rule (1b report, level 2)

B, A = collect(a.before, a.group_col), collect(a.after, a.group_col)
keys = list(dict.fromkeys(list(B) + list(A)))

rows, n_changed, n_unexpected = [], 0, 0
for k in keys:
    cat, label, vb = B.get(k, (A[k][0], A[k][1], None))
    va = A.get(k, (cat, label, None))[2]
    if k not in B:
        status = "new metric"
    elif k not in A:
        status = "missing in after"
        n_unexpected += 1
    elif same(vb, va):
        status = "same"
    elif a.mode == "reproducibility" and a.tolerance != "exact" and k.startswith(REPORT_ONLY):
        n_changed += 1
        status = "report only"
    elif a.mode == "reproducibility" and a.tolerance == "level2" and tolerance_rule(k, "level2") is None:
        n_changed += 1
        status = "report only"
    else:
        n_changed += 1
        if a.mode == "regression":
            status = "expected: " + ", ".join(EXPECTED.get(cat, [])) if EXPECTED.get(cat) else "UNEXPECTED"
        elif cat in NONDETERMINISTIC and a.tolerance == "exact":
            status = "expected: non-deterministic assembly"
        else:
            rule = tolerance_rule(k, a.tolerance) if a.tolerance in ("1b", "level2") else None
            try:
                ok = rule is not None and rule[0](float(vb), float(va))
            except (TypeError, ValueError):
                ok = False
            status = f"within tolerance ({rule[1]})" if ok else "UNEXPECTED"
        if status == "UNEXPECTED":
            n_unexpected += 1
    rows.append((label, fmt(vb), fmt(va), "" if status == "same" else delta(vb, va), status))

L = [f"# MeToReS run comparison ({a.mode})", "",
     f"- before: `{os.path.abspath(a.before)}`",
     f"- after:  `{os.path.abspath(a.after)}`", "",
     f"**{len(keys)} metrics compared · {n_changed} changed · "
     f"{n_unexpected} unexpected.**", "",
     "| Metric | Before | After | Change | Status |", "|---|---|---|---|---|"]
L += [f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]} |" for r in rows]

if a.mode == "reproducibility":
    if a.tolerance == "level2":
        checks = level2_checks(a.before, a.after, a.group_col, a.samplesheet,
                               a.deg_group_col or a.group_col, a.alpha,
                               math.log2(a.fold_change), a.ko_min_count, a.formula)
    else:
        checks = table_checks(a.before, a.after, a.group_col, a.tolerance)
    n_fail = sum(1 for c in checks if c[3] is False)
    n_unexpected += n_fail
    L += ["", f"## Table-level checks ({a.tolerance})", "",
          "| Check | Result | Criterion | Pass |", "|---|---|---|---|"]
    L += [f"| {c[0]} | {c[1]} | {c[2]} | "
          f"{'—' if c[3] is None else ('yes' if c[3] else '**NO**')} |" for c in checks]
    L[5] = (f"**{len(keys)} metrics compared · {n_changed} changed · "
            f"{len(checks)} table checks, {n_fail} failed · {n_unexpected} unexpected.**")

if a.mode == "regression":
    L += ["", "## Outputs added by the new code", "", "| Output | Present in after |", "|---|---|"]
    for path, lab in NEW_OUTPUTS:
        p = os.path.join(a.after, path.format(gc=a.group_col))
        ok = os.path.exists(p)
        if not ok:
            n_unexpected += 1
        L.append(f"| {lab} (`{path.format(gc=a.group_col)}`) | {'yes' if ok else '**NO**'} |")
    L[5] = (f"**{len(keys)} metrics compared · {n_changed} changed · "
            f"{n_unexpected} unexpected (incl. missing new outputs).**")

os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
with open(a.out, "w") as fh:
    fh.write("\n".join(L) + "\n")
print(f"{len(keys)} metrics · {n_changed} changed · {n_unexpected} unexpected -> {a.out}")
