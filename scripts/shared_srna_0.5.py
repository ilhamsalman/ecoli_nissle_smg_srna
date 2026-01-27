# shared_srna.py
import pandas as pd

MIN_RECIP = 0.5

exp = pd.read_csv("srna_exponential.csv", sep="\t")
sta = pd.read_csv("srna_stationary.csv", sep="\t")

def overlap(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)

matches = []

for (genome, strand), sta_g in sta.groupby(["Genome", "Strand"]):
    exp_g = exp[(exp["Genome"] == genome) & (exp["Strand"] == strand)]
    if exp_g.empty:
        continue

    for _, s in sta_g.iterrows():
        s_start, s_end = int(s["Start"]), int(s["End"])
        s_len = s_end - s_start + 1

        for _, e in exp_g.iterrows():
            e_start, e_end = int(e["Start"]), int(e["End"])
            e_len = e_end - e_start + 1

            ov = overlap(s_start, s_end, e_start, e_end)
            if ov == 0:
                continue

            recip_sta = ov / s_len
            recip_exp = ov / e_len

            if recip_sta >= MIN_RECIP and recip_exp >= MIN_RECIP:
                matches.append({
                    "Genome": genome,
                    "Strand": strand,

                    "stationary_Name": s["Name"],
                    "stationary_Start": s_start,
                    "stationary_End": s_end,
                    "stationary_Length": s_len,

                    "exponential_Name": e["Name"],
                    "exponential_Start": e_start,
                    "exponential_End": e_end,
                    "exponential_Length": e_len,

                    "overlap_nt": ov,
                    "recip_stationary": recip_sta,
                    "recip_exponential": recip_exp,
                })

out = pd.DataFrame(matches)
out.to_csv("shared_srna_minrecip_0.5_with_lengths.csv", index=False)

print("Shared matches (min_recip=0.5):", out.shape[0])
print("Wrote: shared_srna_minrecip_0.5_with_lengths.csv")
