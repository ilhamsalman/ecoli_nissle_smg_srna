#!/usr/bin/env python3
import pandas as pd
import sys

if len(sys.argv) != 3:
    sys.exit("Usage: 01_make_bed_from_annogesic.py <srna_table.tsv> <out.bed>")

infile, outbed = sys.argv[1], sys.argv[2]
df = pd.read_csv(infile, sep="\t")

required = {"Genome","Start","End","Name","Strand"}
missing = required - set(df.columns)
if missing:
    sys.exit(f"ERROR: missing columns in {infile}: {', '.join(sorted(missing))}")

bed = pd.DataFrame({
    "chrom": df["Genome"].astype(str),
    "start": df["Start"].astype(int) - 1,
    "end": df["End"].astype(int),
    "name": df["Name"].astype(str),
    "score": 0,
    "strand": df["Strand"].astype(str),
})
bed.to_csv(outbed, sep="\t", header=False, index=False)
print(f"Wrote: {outbed} (rows={len(bed)})")
