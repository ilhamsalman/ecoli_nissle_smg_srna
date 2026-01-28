import pandas as pd

exp = pd.read_csv("srna_exponential.csv", sep="\t")
sta = pd.read_csv("srna_stationary.csv", sep="\t")
shared = pd.read_csv("shared_srna_minrecip_0.8.csv")

# columns to carry over if present
extra_cols = [
    "sRNA_types",
    "Overlap_CDS_forward", "Overlap_nts_forward",
    "Overlap_CDS_reverse", "Overlap_nts_reverse",
    "nr_hit_top3|ID|e-value|score",
    "sRNA_hit|e-value|score",
]

def pick_existing(df, cols):
    return [c for c in cols if c in df.columns]

exp_cols = ["Name"] + pick_existing(exp, extra_cols)
sta_cols = ["Name"] + pick_existing(sta, extra_cols)

exp_anno = exp[exp_cols].copy().rename(columns={"Name": "exponential_Name"})
sta_anno = sta[sta_cols].copy().rename(columns={"Name": "stationary_Name"})

# rename to keep provenance clear
exp_anno = exp_anno.rename(columns={c: f"exponential_{c}" for c in exp_anno.columns if c != "exponential_Name"})
sta_anno = sta_anno.rename(columns={c: f"stationary_{c}" for c in sta_anno.columns if c != "stationary_Name"})

out = shared.merge(sta_anno, on="stationary_Name", how="left").merge(exp_anno, on="exponential_Name", how="left")
out.to_csv("shared_srna_minrecip_0.8_annotated_plus_locusinfo.csv", index=False)

print("Wrote: shared_srna_minrecip_0.8_annotated_plus_locusinfo.csv")
