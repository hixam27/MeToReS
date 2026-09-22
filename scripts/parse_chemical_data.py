#!/usr/bin/env python3
"""parse_chemical_data.py — load the chemical dataset, harmonize sample names against the
RNA samplesheet, and emit clean per-fraction CSVs keyed by the RNA-side sample name.

Fractions are detected generically by column suffix (e.g. "Cu_I", "Cu_II",
"Cu_PT") rather than hardcoded to a fixed pair, so any number of fractions
present in the sheet is handled without code changes.

Handles known naming mismatches between the chemistry team's labels and the
RNA samplesheet (e.g. CY_* -> CYP_*) via a configurable rename map, and
reports any samples that still don't match after renaming so they can be
caught before being silently dropped.
"""
import argparse
import sys
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--xlsx", required=True, help="Chemical_data_*.xlsx")
ap.add_argument("--sheet", default="14_DATASET",
                help="Sheet name containing the RNA-subset samples")
ap.add_argument("--samplesheet", required=True, help="RNA samplesheet.csv")
ap.add_argument("--rename", default="CY_:CYP_",
                help="Comma-separated old:new prefix pairs applied to SampleID, "
                     "e.g. 'CY_:CYP_,OLD_:NEW_'")
ap.add_argument("--fractions", default="I:bioavailable,PT:pseudototal",
                help="Comma-separated suffix:output_name pairs identifying "
                     "each fraction block by its column suffix, e.g. "
                     "'I:bioavailable,II:reducible,III:oxidizable,"
                     "IV:residual,PT:pseudototal'")
ap.add_argument("--out_dir", required=True)
ap.add_argument("--header_source_sheet", default="FULL_DATASET",
                help="Sheet to borrow column headers from if --sheet has no header row")
a = ap.parse_args()

# ── Load ──────────────────────────────────────────────────────────────
raw_peek = pd.read_excel(a.xlsx, sheet_name=a.sheet, header=None, nrows=1)
first_row = raw_peek.iloc[0].tolist()

def looks_like_header(row):
    return all(isinstance(v, str) for v in row if pd.notna(v)) and \
           any(str(v).strip().lower() == "sampleid" for v in row if pd.notna(v))

if looks_like_header(first_row):
    chem = pd.read_excel(a.xlsx, sheet_name=a.sheet)
else:
    print(f"NOTE: '{a.sheet}' has no header row; borrowing column names "
          f"from '{a.header_source_sheet}'.", file=sys.stderr)
    header_df = pd.read_excel(a.xlsx, sheet_name=a.header_source_sheet, nrows=0)
    expected_cols = list(header_df.columns)
    chem = pd.read_excel(a.xlsx, sheet_name=a.sheet, header=None)
    if chem.shape[1] != len(expected_cols):
        sys.exit(f"ERROR: '{a.sheet}' has {chem.shape[1]} columns but "
                 f"'{a.header_source_sheet}' header has {len(expected_cols)}. "
                 f"Cannot safely align columns — check the file structure.")
    chem.columns = expected_cols

if "SampleID" not in chem.columns:
    sys.exit(f"ERROR: '{a.sheet}' has no SampleID column. "
              f"Found: {list(chem.columns)}")

ss = pd.read_csv(a.samplesheet)
if "sample" not in ss.columns:
    sys.exit("ERROR: samplesheet must have a 'sample' column")
rna_samples = set(ss["sample"].astype(str))

# ── Apply rename map (prefix substitutions) ─────────────────────────
rename_pairs = []
for pair in a.rename.split(","):
    pair = pair.strip()
    if not pair:
        continue
    old, new = pair.split(":")
    rename_pairs.append((old, new))

def apply_renames(sid):
    for old, new in rename_pairs:
        if sid.startswith(old):
            return new + sid[len(old):]
    return sid

chem["SampleID_original"] = chem["SampleID"].astype(str)
chem["SampleID"] = chem["SampleID_original"].apply(apply_renames)

# ── Match report ─────────────────────────────────────────────────────
chem_samples = set(chem["SampleID"])
matched = chem_samples & rna_samples
unmatched_chem = chem_samples - rna_samples
unmatched_rna = rna_samples - chem_samples

report_lines = []
report_lines.append(f"Chemical samples (sheet '{a.sheet}'): {len(chem_samples)}")
report_lines.append(f"RNA samplesheet samples: {len(rna_samples)}")
report_lines.append(f"Matched: {len(matched)}")
report_lines.append("")
if rename_pairs:
    report_lines.append("Rename rules applied:")
    for old, new in rename_pairs:
        n_affected = (chem["SampleID_original"].str.startswith(old)).sum()
        if n_affected:
            report_lines.append(f"  {old}* -> {new}*  ({n_affected} samples)")
    report_lines.append("")
if unmatched_chem:
    report_lines.append(f"Chemical samples with NO match in RNA samplesheet ({len(unmatched_chem)}):")
    for s in sorted(unmatched_chem):
        report_lines.append(f"  {s}")
    report_lines.append("")
if unmatched_rna:
    report_lines.append(f"RNA samples with NO match in chemical data ({len(unmatched_rna)}):")
    for s in sorted(unmatched_rna):
        report_lines.append(f"  {s}")
    report_lines.append("")

report_text = "\n".join(report_lines)
print(report_text)

import os
os.makedirs(a.out_dir, exist_ok=True)
with open(os.path.join(a.out_dir, "sample_matching_report.txt"), "w") as f:
    f.write(report_text)

if unmatched_rna:
    print("\nWARNING: some RNA samples have no chemical data and will be "
          "excluded from chemical integration analyses.", file=sys.stderr)

chem_matched = chem[chem["SampleID"].isin(rna_samples)].copy()
if chem_matched.empty:
    sys.exit("ERROR: no samples matched between chemical data and RNA "
              "samplesheet after renaming. Check --rename and the sheet name.")

meta_cols = ["SampleID", "Country", "Season", "Replicate"]
meta_cols = [c for c in meta_cols if c in chem_matched.columns]

# ── Parse --fractions into (suffix, output_name) pairs, longest suffix
fraction_pairs = []
for pair in a.fractions.split(","):
    pair = pair.strip()
    if not pair:
        continue
    suffix, out_name = pair.split(":")
    fraction_pairs.append((suffix.strip(), out_name.strip()))
fraction_pairs.sort(key=lambda p: len(p[0]), reverse=True)

# ── Split into one CSV per fraction, generic over however many are
# configured. A column is claimed by the FIRST (longest) matching suffix
claimed = set()
saved = []
for suffix, out_name in fraction_pairs:
    tag = f"_{suffix}"
    cols = [c for c in chem_matched.columns
            if c.endswith(tag) and c not in claimed]
    claimed.update(cols)
    if not cols:
        print(f"WARNING: no columns found with suffix '{tag}' "
              f"(fraction '{out_name}') — skipping.", file=sys.stderr)
        continue

    frac_df = chem_matched[meta_cols + cols].copy()
    frac_df.columns = meta_cols + [c[:-len(tag)] for c in cols]

    out_file = os.path.join(a.out_dir, f"chemical_data_{out_name}.csv")
    frac_df.to_csv(out_file, index=False)
    saved.append((out_name, len(cols)))

chem_matched.to_csv(os.path.join(a.out_dir, "chemical_data_combined.csv"), index=False)

print(f"\nSaved {len(chem_matched)} matched samples:")
for out_name, n_elements in saved:
    print(f"  chemical_data_{out_name}.csv  ({n_elements} elements)")
print(f"  chemical_data_combined.csv      (all columns)")
