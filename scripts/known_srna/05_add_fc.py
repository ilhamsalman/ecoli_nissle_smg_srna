#!/usr/bin/env python3
import pandas as pd
import sys
import math

if len(sys.argv) != 4:
    sys.exit("Usage: 05_add_fc.py <known_hits.csv> <de_table.(csv|tsv)> <out.csv>")

known_csv, de_file, out_csv = sys.argv[1:]

# Read known hits
known = pd.read_csv(known_csv)

# Extract ANNOgesic Name from qseqid (format: sRNA_XXXXX::CP007799.1:... )
if "qseqid" not in known.columns:
    sys.exit("ERROR: known hits file must contain 'qseqid' column")

known["srna_id"] = known["qseqid"].astype(str).str.split("::").str[0]

# Read DE table (auto sep)
sep = "\t" if de_file.endswith(".tsv") or de_file.endswith(".txt") else ","
de = pd.read_csv(de_file, sep=sep)

# Identify ID column in DE table
candidate_id_cols = ["srna_id", "Name", "gene", "id", "feature", "rowname"]
id_col = None
for c in candidate_id_cols:
    if c in de.columns:
        id_col = c
        break
if id_col is None:
    sys.exit("ERROR: cannot find an ID column in DE table. Expected one of: " + ", ".join(candidate_id_cols))

# Identify log2FC column
candidate_l2fc = ["log2FoldChange", "log2FC", "log2_fc", "LFC"]
l2fc_col = None
for c in candidate_l2fc:
    if c in de.columns:
        l2fc_col = c
        break
if l2fc_col is None:
    sys.exit("ERROR: cannot find log2FC column in DE table. Expected one of: " + ", ".join(candidate_l2fc))

# Optional columns
padj_col = "padj" if "padj" in de.columns else None
pval_col = "pvalue" if "pvalue" in de.columns else None

# Normalize DE IDs to match srna_id
de = de.copy()
de["srna_id"] = de[id_col].astype(str)

# Keep only relevant columns
keep_cols = ["srna_id", l2fc_col]
if padj_col: keep_cols.append(padj_col)
if pval_col: keep_cols.append(pval_col)
de_small = de[keep_cols].drop_duplicates("srna_id")

# Merge
out = known.merge(de_small, on="srna_id", how="left")

# Add FC = 2^(log2FC)
def to_fc(x):
    try:
        if pd.isna(x):
            return pd.NA
        return float(2 ** float(x))
    except Exception:
        return pd.NA

out["FC"] = out[l2fc_col].apply(to_fc)

# Rename log2FC column to a standard name in output
if l2fc_col != "log2FC":
    out = out.rename(columns={l2fc_col: "log2FC"})

out.to_csv(out_csv, index=False)
print(f"Wrote: {out_csv}")
print("Rows:", len(out), "with FC annotated:", out["FC"].notna().sum())
