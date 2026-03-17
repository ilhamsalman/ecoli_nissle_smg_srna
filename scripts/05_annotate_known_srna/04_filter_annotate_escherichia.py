#!/usr/bin/env python3
import pandas as pd
import sys

if len(sys.argv) != 4:
    sys.exit("Usage: 04_filter_annotate_escherichia.py <blast.tsv> <BSRD.fa> <out.csv>")

blast_tsv, bsrd_fa, out_csv = sys.argv[1:]

MIN_IDENT = 80.0
MIN_ALN_LEN = 50
MAX_EVALUE = 1e-5

cols = ["qseqid","sseqid","pident","length","evalue","bitscore"]

# Parse BSRD FASTA headers: >ID|Organism|Name
records = []
with open(bsrd_fa, "r", encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith(">"):
            h = line[1:].strip()
            parts = h.split("|", 2)
            records.append({
                "bsrd_id": parts[0],
                "organism": parts[1] if len(parts) > 1 else "",
                "known_sRNA_name": parts[2] if len(parts) > 2 else ""
            })
lookup = pd.DataFrame(records).drop_duplicates("bsrd_id")

df = pd.read_csv(blast_tsv, sep="\t", names=cols)

df = df[(df.pident >= MIN_IDENT) &
        (df.length >= MIN_ALN_LEN) &
        (df.evalue <= MAX_EVALUE)].copy()

df["bsrd_id"] = df["sseqid"].astype(str).str.split("|").str[0]
df = df.merge(lookup, on="bsrd_id", how="left")

df = df[df["organism"].astype(str).str.contains("Escherichia", na=False)].copy()

df = df.sort_values(["qseqid","bitscore","evalue"], ascending=[True, False, True])
best = df.groupby("qseqid", as_index=False).first()

best.to_csv(out_csv, index=False)
print(f"Thresholds: identity>={MIN_IDENT}% aln_len>={MIN_ALN_LEN} evalue<={MAX_EVALUE}")
print(f"Wrote: {out_csv}")
print("Known sRNAs (unique query IDs):", best["qseqid"].nunique())
