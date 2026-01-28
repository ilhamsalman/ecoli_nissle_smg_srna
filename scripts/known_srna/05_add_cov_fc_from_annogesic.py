#!/usr/bin/env python3
import pandas as pd
import sys
import numpy as np

if len(sys.argv) != 5:
    sys.exit(
        "Usage: 05_add_cov_fc_from_annogesic.py <known_hits.csv> <annogesic_srna_table.tsv> "
        "<phase:exp|sta> <out.csv>\n"
        "Example: 05_add_cov_fc_from_annogesic.py known.csv srna_exponential.csv exp out.csv"
    )

known_csv, annogesic_tsv, phase, out_csv = sys.argv[1:]
phase = phase.lower().strip()
if phase not in {"exp", "sta"}:
    sys.exit("ERROR: phase must be 'exp' or 'sta'")

PSEUDOCOUNT = 1.0  # conservative, avoids divide-by-zero

# ----------------------------
# Load known hits
# ----------------------------
known = pd.read_csv(known_csv)
if "qseqid" not in known.columns:
    sys.exit("ERROR: known hits file must contain 'qseqid' column")

known["srna_id"] = known["qseqid"].astype(str).str.split("::").str[0]

# ----------------------------
# Load ANNOgesic sRNA table
# ----------------------------
srna = pd.read_csv(annogesic_tsv, sep="\t")

if "Name" not in srna.columns:
    sys.exit("ERROR: ANNOgesic table must contain 'Name' column")

# Coverage columns depend on phase naming
# exp example:
# Avg_coverage:MG_exp_1_plus/MG_exp_1_minus
# Avg_coverage:MG_exp_2_plus/MG_exp_2_minus
# Avg_coverage:NG_exp_1_plus/NG_exp_1_minus
# Avg_coverage:NG_exp_2_plus/NG_exp_2_minus
#
# sta often uses MG_sta_*, NG_sta_* (or MG_stationary_*)
# We'll search by patterns.

def find_cov_cols(df, mg_or_ng, phase):
    # Try common patterns first
    patterns = [
        f"Avg_coverage:{mg_or_ng}_{phase}_1_plus/{mg_or_ng}_{phase}_1_minus",
        f"Avg_coverage:{mg_or_ng}_{phase}_2_plus/{mg_or_ng}_{phase}_2_minus",
    ]
    # If not found, fall back to fuzzy match
    found = []
    for p in patterns:
        if p in df.columns:
            found.append(p)

    if len(found) == 2:
        return found

    # Fuzzy fallback: find columns containing mg_or_ng and phase and "Avg_coverage"
    cand = [c for c in df.columns
            if c.startswith("Avg_coverage:") and (mg_or_ng in c) and (phase in c)]
    # keep only replicate-like columns (contain _1_ or _2_)
    cand = [c for c in cand if ("_1_" in c or "_2_" in c)]
    cand = sorted(cand)

    # choose first two if >2
    if len(cand) >= 2:
        return cand[:2]

    return []

mg_cols = find_cov_cols(srna, "MG", phase)
ng_cols = find_cov_cols(srna, "NG", phase)

if len(mg_cols) != 2 or len(ng_cols) != 2:
    msg = (
        "ERROR: Could not identify 2 MG and 2 NG coverage columns.\n"
        f"Phase={phase}\n"
        f"MG candidates found: {mg_cols}\n"
        f"NG candidates found: {ng_cols}\n"
        "Tip: run: python3 - <<'PY'\n"
        "import pandas as pd\n"
        "df=pd.read_csv('FILE', sep='\\t')\n"
        "print([c for c in df.columns if c.startswith('Avg_coverage:')])\n"
        "PY\n"
    )
    sys.exit(msg)

# Build small coverage table keyed by sRNA Name
cov = srna[["Name"] + mg_cols + ng_cols].copy()
cov = cov.rename(columns={"Name": "srna_id"})

# Ensure numeric
for c in mg_cols + ng_cols:
    cov[c] = pd.to_numeric(cov[c], errors="coerce")

cov["MG_mean_cov"] = cov[mg_cols].mean(axis=1, skipna=True)
cov["NG_mean_cov"] = cov[ng_cols].mean(axis=1, skipna=True)

# Coverage-based FC
cov["FC_cov_MG_over_NG"] = (cov["MG_mean_cov"] + PSEUDOCOUNT) / (cov["NG_mean_cov"] + PSEUDOCOUNT)
cov["log2FC_cov_MG_over_NG"] = np.log2(cov["FC_cov_MG_over_NG"])

# Merge
out = known.merge(cov[[
    "srna_id",
    mg_cols[0], mg_cols[1],
    ng_cols[0], ng_cols[1],
    "MG_mean_cov", "NG_mean_cov",
    "FC_cov_MG_over_NG", "log2FC_cov_MG_over_NG"
]], on="srna_id", how="left")

out.to_csv(out_csv, index=False)
print("Phase:", phase)
print("MG cols:", mg_cols)
print("NG cols:", ng_cols)
print(f"Wrote: {out_csv}")
print("Rows:", len(out))
print("Rows with coverage annotated:", out["MG_mean_cov"].notna().sum())
