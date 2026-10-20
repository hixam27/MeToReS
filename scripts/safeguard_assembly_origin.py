#!/usr/bin/env python3
"""safeguard_assembly_origin.py — S4: assembly-origin enrichment of hit sets.

Every gene in the pooled catalogue keeps the identifier of the sample it was
assembled from ({sra}_k141_...). A gene assembled from one sample's reads
tends to be quantified best in that same sample. With few samples per
group, this can create differential or association signals that follow
*where a gene was assembled* rather than biology.

For each hit set — the up- and down-regulated genes of every DESeq2
contrast, and the positive and negative MaAsLin2 associations (pooled and
per element) — this compares the origin of the hit genes with the origin of
all genes tested in that analysis (Fisher's exact test, per origin sample
and per origin group). P-values are BH-corrected across every set and
origin tested in the run, so that many small hit sets do not produce
chance flags. An origin that is strongly over-represented among the hits is
flagged.

A flag does not prove an artefact: a gene genuinely specific to one site
is also assembled mainly from that site. It means the hits cannot be
separated from assembly origin, and should be read as site-structured.

Matched null. Hit genes are usually more highly expressed and detected in
more samples than the average tested gene, and such genes may assemble
preferentially from deep, complex samples. Comparing hits with all tested
genes would then flag an origin for reasons of expression level alone. With
--matrix, tested genes are stratified by detection count (number of samples
with abundance > 0) and by decile of mean log abundance, and each origin is
tested with a Mantel–Haenszel test across strata: hits are compared only
with tested genes of the same detection and expression level. The common
(Mantel–Haenszel) odds ratio is the effect size; a set is flagged when
q < alpha and that odds ratio is >= --min_odds_ratio. The unmatched Fisher
test is kept in the table for comparison.

For DESeq2 sets the column 'origin_in_contrast' says whether the
over-represented origin belongs to one of the two groups being compared.
If it does, the self-assembly effect above can produce the signal directly.
If it does not (genes assembled from a third site), the hits are instead
dominated by lineages best assembled from that site: a community-composition
signal, worth checking taxonomically, but not the self-assembly effect.

Report only: no result is changed.
"""
import argparse
import json
import os

import math

import pandas as pd


def fisher_exact(table):
    """Two-sided Fisher's exact test for a 2x2 table [[a, b], [c, d]].
    Returns (sample odds ratio, p). Same definition as scipy.stats.fisher_exact,
    implemented here to avoid a SciPy dependency."""
    (a, b), (c, d) = table
    r1, c1, n = a + b, a + c, a + b + c + d
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)

    def logp(x):
        return (math.lgamma(r1 + 1) + math.lgamma(n - r1 + 1) + math.lgamma(c1 + 1)
                + math.lgamma(n - c1 + 1) - math.lgamma(n + 1) - math.lgamma(x + 1)
                - math.lgamma(r1 - x + 1) - math.lgamma(c1 - x + 1)
                - math.lgamma(n - r1 - c1 + x + 1))
    lp_obs = logp(a)
    tol = 1e-7
    p = sum(math.exp(logp(x)) for x in range(lo, hi + 1) if logp(x) <= lp_obs + tol)
    odds = (a * d / (b * c)) if b * c > 0 else (math.inf if a * d > 0 else math.nan)
    return odds, min(1.0, p)


def mantel_haenszel(strata):
    """strata: list of 2x2 tables [[a, b], [c, d]] (rows: hit / not hit,
    columns: from origin / not). Returns (common odds ratio, expected a,
    two-sided p of the CMH test with continuity correction)."""
    num = den = sum_a = sum_e = sum_v = 0.0
    for (a_, b_), (c_, d_) in strata:
        n = a_ + b_ + c_ + d_
        if n < 2:
            continue
        r1, r2, c1, c2 = a_ + b_, c_ + d_, a_ + c_, b_ + d_
        num += a_ * d_ / n
        den += b_ * c_ / n
        sum_a += a_
        sum_e += r1 * c1 / n
        sum_v += r1 * r2 * c1 * c2 / (n * n * (n - 1))
    odds = num / den if den > 0 else (math.inf if num > 0 else math.nan)
    if sum_v <= 0:
        return odds, sum_e, 1.0
    chi = max(0.0, abs(sum_a - sum_e) - 0.5) ** 2 / sum_v
    return odds, sum_e, math.erfc(math.sqrt(chi / 2))


def bh(p):
    p = pd.Series(p, dtype=float)
    n = p.notna().sum()
    if n == 0:
        return p
    order = p.dropna().sort_values()
    ranks = pd.Series(range(1, len(order) + 1), index=order.index)
    q = (order * n / ranks).iloc[::-1].cummin().iloc[::-1].clip(upper=1)
    return q.reindex(p.index)


ap = argparse.ArgumentParser()
ap.add_argument("--samplesheet", required=True)
ap.add_argument("--deg_dir", help="DESeq2 results directory (one sub-directory per contrast)")
ap.add_argument("--contrasts", default="", help="space-separated contrast names")
ap.add_argument("--deg_group_col", default=None, help="samplesheet column used for DESeq2 groups")
ap.add_argument("--maaslin_all", default=None)
ap.add_argument("--maaslin_group_col", default=None,
                help="samplesheet column summarising origin for MaAsLin2 sets")
ap.add_argument("--alpha", type=float, default=0.05)
ap.add_argument("--min_odds_ratio", type=float, default=2.0)
ap.add_argument("--min_set_size", type=int, default=5)
ap.add_argument("--matrix", default=None,
                help="gene x sample abundance matrix; enables the matched null")
ap.add_argument("--n_expr_bins", type=int, default=10)
ap.add_argument("--out_tsv", required=True)
ap.add_argument("--out_json", required=True)
a = ap.parse_args()

ss = pd.read_csv(a.samplesheet)
sras = sorted(ss["sra"].astype(str), key=len, reverse=True)
sra_to_sample = dict(zip(ss["sra"].astype(str), ss["sample"]))


stratum_of = None
if a.matrix:
    m = pd.read_csv(a.matrix, index_col=0)
    m.index = m.index.astype(str)
    m = m.apply(pd.to_numeric, errors="coerce").fillna(0)
    det = (m > 0).sum(axis=1)
    lmean = (m.clip(lower=0) + 1).apply(lambda col: col.map(math.log)).mean(axis=1)
    ebin = pd.qcut(lmean.rank(method="first"), a.n_expr_bins, labels=False)
    stratum_of = (det.astype(str) + ":" + ebin.astype(str)).to_dict()
    print(f"S4: matched null on {len(stratum_of)} genes "
          f"(detection count x {a.n_expr_bins} expression bins)")


def origin_sample(gene):
    for s in sras:
        if gene.startswith(s + "_"):
            return sra_to_sample[s]
    return "unknown"


def test_set(set_name, analysis, hits, background, group_col, contrast_groups=None):
    """hits, background: iterables of gene ids (hits must be a subset)."""
    bg = pd.Series(sorted(set(map(str, background))))
    hit_set = set(map(str, hits)) & set(bg)
    if len(hit_set) < a.min_set_size:
        return []
    o_sample = bg.map(origin_sample)
    levels = [("sample", o_sample)]
    if group_col and group_col in ss.columns:
        g_map = dict(zip(ss["sample"], ss[group_col]))
        levels.append((group_col, o_sample.map(lambda s: g_map.get(s, "unknown"))))
    is_hit = bg.isin(hit_set)
    stratum = bg.map(stratum_of).fillna("unmatched") if stratum_of else None
    g_of = dict(zip(ss["sample"], ss[group_col])) if group_col and group_col in ss.columns else {}
    rows = []
    for level, origin in levels:
        set_rows = []
        for o in sorted(origin.unique()):
            in_o = origin == o
            n11 = int((is_hit & in_o).sum()); n12 = int((is_hit & ~in_o).sum())
            n21 = int((~is_hit & in_o).sum()); n22 = int((~is_hit & ~in_o).sum())
            odds_u, p_u = fisher_exact([[n11, n12], [n21, n22]])
            exp_u = is_hit.sum() * in_o.mean()
            if stratum is not None:
                t = pd.crosstab([stratum, is_hit], in_o).reindex(columns=[True, False], fill_value=0)
                tables = []
                for st in stratum.unique():
                    def cell(h, io):
                        try:
                            return int(t.loc[(st, h), io])
                        except KeyError:
                            return 0
                    tables.append([[cell(True, True), cell(True, False)],
                                   [cell(False, True), cell(False, False)]])
                odds, exp_m, p = mantel_haenszel(tables)
            else:
                odds, exp_m, p = odds_u, exp_u, p_u
            set_rows.append({
                "analysis": analysis, "set": set_name, "level": level, "origin": o,
                "n_hits": int(is_hit.sum()), "hits_from_origin": n11,
                "pct_hits_from_origin": round(100 * n11 / is_hit.sum(), 1),
                "pct_background_from_origin": round(100 * in_o.mean(), 2),
                "expected_hits": round(exp_m, 1),
                "odds_ratio": odds, "p": p,
                "expected_hits_unmatched": round(exp_u, 1),
                "odds_ratio_unmatched": odds_u, "p_unmatched": p_u,
                "null": "matched" if stratum is not None else "all tested genes",
                "origin_in_contrast": (None if contrast_groups is None else
                                       bool((g_of.get(o, o) if level == "sample" else o)
                                            in contrast_groups))})
        rows += set_rows
    return rows


rows = []
# ── DESeq2 contrasts ────────────────────────────────────────────────────
for c in a.contrasts.split():
    d = os.path.join(a.deg_dir, c)
    res_all = os.path.join(d, "deseq2_results_all.csv")
    degs = os.path.join(d, "differential_genes.csv")
    if not (os.path.exists(res_all) and os.path.exists(degs)):
        print(f"S4: skipping {c} (results not found)")
        continue
    bg = pd.read_csv(res_all, usecols=[0]).iloc[:, 0]
    dg = pd.read_csv(degs)
    gid = dg.columns[0]
    for direction in ("Up", "Down"):
        hits = dg.loc[dg["change"].astype(str).str.lower() == direction.lower(), gid]
        groups = set(c.split("_vs_", 1)) if "_vs_" in c else None
        rows += test_set(f"{c}:{direction.lower()}", "DESeq2", hits, bg, a.deg_group_col, groups)

# ── MaAsLin2 ────────────────────────────────────────────────────────────
if a.maaslin_all and os.path.exists(a.maaslin_all):
    ml = pd.read_csv(a.maaslin_all, sep="\t")
    qcol = "qval_across_elements" if "qval_across_elements" in ml.columns else "qval"
    bg = ml["feature"].astype(str).unique()
    sig = ml[ml[qcol] < a.alpha]
    for sign, lab in ((1, "positive"), (-1, "negative")):
        s = sig[(sig["coef"] > 0) if sign > 0 else (sig["coef"] < 0)]
        rows += test_set(f"all_elements:{lab}", "MaAsLin2", s["feature"], bg, a.maaslin_group_col)
        for el, se in s.groupby("element"):
            rows += test_set(f"{el}:{lab}", "MaAsLin2", se["feature"], bg, a.maaslin_group_col)

cols = ["analysis", "set", "level", "origin", "n_hits", "hits_from_origin",
        "pct_hits_from_origin", "pct_background_from_origin", "expected_hits",
        "odds_ratio", "p", "q", "enriched", "origin_in_contrast",
        "expected_hits_unmatched", "odds_ratio_unmatched", "p_unmatched", "null"]
out = pd.DataFrame(rows, columns=[c for c in cols if c not in ("q", "enriched")])
# One multiple-testing family per run: BH across all sets and origins.
out["q"] = bh(out["p"]).values if len(out) else []
out["enriched"] = ((out["q"] < a.alpha) & (out["odds_ratio"] >= a.min_odds_ratio)
                   & (out["hits_from_origin"] >= 3)) if len(out) else []
out = out[cols]
os.makedirs(os.path.dirname(os.path.abspath(a.out_tsv)), exist_ok=True)
out.to_csv(a.out_tsv, sep="\t", index=False)

flags = out[out["enriched"]]
sets = out[["analysis", "set", "n_hits"]].drop_duplicates()
summary = {
    "safeguard": "S4_assembly_origin",
    "alpha": a.alpha, "min_odds_ratio": a.min_odds_ratio, "min_set_size": a.min_set_size,
    "null": "matched (detection count x expression decile, Mantel-Haenszel)" if stratum_of
            else "all tested genes (Fisher)",
    "n_sets_tested": int(len(sets)),
    "n_sets_flagged": int(flags[["analysis", "set"]].drop_duplicates().shape[0]),
    "flagged": [
        {"analysis": r["analysis"], "set": r["set"], "level": r["level"], "origin": r["origin"],
         "hits_from_origin": int(r["hits_from_origin"]), "n_hits": int(r["n_hits"]),
         "pct_hits": float(r["pct_hits_from_origin"]),
         "pct_background": float(r["pct_background_from_origin"]),
         "expected_hits": float(r["expected_hits"]),
         "odds_ratio": (None if r["odds_ratio"] == float("inf") else round(float(r["odds_ratio"]), 2)),
         "q": float(r["q"]),
         "origin_in_contrast": (None if pd.isna(r["origin_in_contrast"]) else bool(r["origin_in_contrast"]))}
        for _, r in flags.sort_values("q").iterrows()],
}
with open(a.out_json, "w") as fh:
    json.dump(summary, fh, indent=2)

print(f"S4 assembly origin: {summary['n_sets_flagged']} of {summary['n_sets_tested']} hit sets "
      f"have an over-represented origin (q < {a.alpha}, OR >= {a.min_odds_ratio}; "
      f"null: {summary['null']})")
for (an, st), grp in flags.groupby(["analysis", "set"], sort=False):
    top = grp.sort_values("q").iloc[0]
    print(f"SAFEGUARD_WARNING: S4 {an} {st}: {top['hits_from_origin']}/{top['n_hits']} hits "
          f"({top['pct_hits_from_origin']}%) assembled from {top['level']} '{top['origin']}' "
          f"vs {top['expected_hits']} expected (OR = {top['odds_ratio']:.1f}, q = {top['q']:.1e})"
          + ("" if pd.isna(top["origin_in_contrast"]) else
             ("; origin is one of the compared groups" if top["origin_in_contrast"]
              else "; origin is outside the compared groups (composition, not self-assembly)")))
