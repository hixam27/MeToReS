#!/usr/bin/env python3
"""safeguards_summary.py — collect the safeguard reports (S1–S4) into one
human-readable file, results/safeguards/SAFEGUARDS_SUMMARY.md.

Each safeguard writes a small JSON summary; this script only reads and
presents them. Nothing is recomputed.
"""
import argparse
import json
import os

ap = argparse.ArgumentParser()
ap.add_argument("--jsons", nargs="+", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()

reports = []
for p in a.jsons:
    with open(p) as fh:
        reports.append(json.load(fh))
by = {}
for r in reports:
    by.setdefault(r["safeguard"], []).append(r)

L = ["# MeToReS — safeguard report", "",
     "These checks flag patterns that, at low replication, can make site or",
     "sample differences look like biological effects. They change no result.",
     "A flag means *read the affected results with care*, not *the result is wrong*.", ""]
n_flags = 0

# S2
for r in by.get("S2_library_complexity", []):
    L += ["## S2 · Low-complexity libraries", ""]
    thr = r["warn_fraction_of_median"] * 100
    if r["n_flagged"]:
        n_flags += 1
        L += [f"**{r['n_flagged']} of {r['n_samples']} samples flagged** "
              f"(below {thr:.0f}% of the run median on at least one measure):", ""]
        L += [f"- `{f['sample']}`: {f['reason']}" for f in r["flagged"]]
        L += ["", "A flagged sample can drive differential results on its own. To exclude it,",
              "remove its row from the samplesheet and re-run.", ""]
    else:
        L += [f"No sample below {thr:.0f}% of the run median.", ""]
    mism = r.get("libtype_mismatch", [])
    if mism:
        n_flags += 1
        L += [f"**{len(mism)} sample(s) declared with a library type that Salmon does not detect**",
              "on a read subsample. A stranded type declared for an unstranded library discards",
              "about half of its fragments; the reverse loses strand information. Check the",
              "library preparation and set `libtype` in the samplesheet:", ""]
        L += [f"- `{m['sample']}`: declared {m['declared']}, detected {m['detected']}" for m in mism]
        L += [""]
    elif r.get("libtype_checked"):
        L += ["Declared library types match Salmon's auto-detection for every sample.", ""]
    lts = r.get("libtypes", {})
    if len(set(lts.values())) > 1:
        L += ["**Library types differ between samples:** " + "; ".join(
            f"{t}: {', '.join(sm for sm, v in lts.items() if v == t)}"
            for t in sorted(set(lts.values()))), ""]
    mixed = r.get("contrasts_mixing_libtypes", [])
    if mixed:
        n_flags += 1
        L += [f"**{len(mixed)} DESeq2 contrast(s) mix library types**, so library type is",
              "confounded with the comparison (unstranded libraries also count antisense",
              "reads). Add it to the design (`deseq2.formula: \"~ libtype + group\"`) once",
              "replication allows:", ""]
        L += [f"- {m['contrast']}: " + "; ".join(f"{g_}: {', '.join(t)}" for g_, t in m["libtypes"].items())
              for m in mixed]
        L += [""]

# S1
s1 = by.get("S1_prevalence", [])
if s1:
    L += ["## S1 · Prevalence of genes behind association hits", ""]
    for r in s1:
        pct = 100 * r["n_hits_low_prevalence"] / r["n_hits"] if r["n_hits"] else 0
        if r["n_hits"] and pct >= 25:
            n_flags += 1
        rule = r.get("rule", f">= {r['min_samples']}/{r['n_samples']} samples")
        L += [f"- **{r['method']}**: {r['n_hits_low_prevalence']} of {r['n_hits']} hits "
              f"({pct:.0f}%) come from genes failing the prevalence rule ({rule})."]
    filt = bool(s1[0].get("filter_applied"))
    L += ["", ("The prevalence filter is on: only genes passing the rule were tested."
               if filt else
               "The prevalence filter is off (`safeguards.prevalence.apply_filter: false`): "
               "low-prevalence genes were tested."),
          "Low-prevalence hits are often presence/absence at a single site. Per-hit",
          "detail, including detection counts: `*_hits_prevalence.tsv`.", ""]

# S3
for r in by.get("S3_singular_fits", []):
    L += ["## S3 · Singular random-effect fits (MaAsLin2)", ""]
    if not r.get("applicable", True):
        L += [f"Not applicable: {r.get('note', '')}.", ""]
        continue
    if r.get("warning"):
        n_flags += 1
    L += [f"- Singular fits ({r['random_effects']} variance estimated at 0): "
          f"**{r['pct_singular_all']}%** of all tests, **{r['pct_singular_hits']}%** of hits "
          f"({r['n_hits_singular']}/{r['n_hits']}); ratio {r['ratio_hits_vs_all']}."]
    if not r.get("refit_reliable", True):
        L += [f"- Note: the refit did not reproduce MaAsLin2's p-values "
              f"(rho = {r['refit_spearman_rho']}); treat these figures with caution."]
    if r.get("warning"):
        L += ["", "Singular fits are much more common among hits than among all tests: these",
              "models ignored the grouping, so the hits may be group differences rather",
              "than covariate effects."]
    wb = r.get("within_between")
    if wb:
        L += [f"- Within/between split by `{wb['grouping']}`: {wb['n_hits_between_p05']} of "
              f"{wb['n_hits']} hits have a between-group term at p < 0.05; "
              f"{wb['n_hits_within_q05']} have a within-group term at q < 0.05."]
    L += [""]

# S4
for r in by.get("S4_assembly_origin", []):
    L += ["## S4 · Assembly-origin enrichment of hit sets", ""]
    if r["n_sets_flagged"]:
        n_flags += 1
        L += [f"**{r['n_sets_flagged']} of {r['n_sets_tested']} hit sets** have an origin "
              f"over-represented among their genes (q < {r['alpha']}, OR ≥ {r['min_odds_ratio']}).",
              f"Null: {r.get('null', 'all tested genes')}. Strongest origin per set:", "",
              "| Analysis | Set | Origin | Hits from origin | Expected | Odds ratio | % of tested genes | q |",
              "|---|---|---|---|---|---|---|---|"]
        seen = set()
        for f in r["flagged"]:
            key = (f["analysis"], f["set"])
            if key in seen:
                continue
            seen.add(key)
            note = ""
            if f.get("origin_in_contrast") is True:
                note = " (compared group)"
            elif f.get("origin_in_contrast") is False:
                note = " (third site)"
            orr = f.get("odds_ratio")
            L += [f"| {f['analysis']} | {f['set']} | {f['origin']}{note} | "
                  f"{f['hits_from_origin']}/{f['n_hits']} ({f['pct_hits']}%) | "
                  f"{f.get('expected_hits', '—')} | {'∞' if orr is None else orr} | "
                  f"{f['pct_background']} | {f['q']:.1e} |"]
        L += ["", "Hits concentrated on genes assembled from one sample or site cannot be",
              "separated from that site. *Compared group*: the self-assembly effect can",
              "produce the signal directly. *Third site*: lineages or genes best assembled",
              "from that site (community composition, or a deep assembly contributing",
              "widely expressed genes). *Expected* is the count under the null used.",
              "Full table: `assembly_origin_enrichment.tsv`.", ""]
    else:
        L += [f"No over-represented origin in {r['n_sets_tested']} hit sets.", ""]

L[6:6] = [f"**Safeguards raising a flag: {n_flags}.**", ""]
os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
with open(a.out, "w") as fh:
    fh.write("\n".join(L) + "\n")
print(f"Safeguard summary written: {a.out} ({n_flags} safeguard(s) raised a flag)")
