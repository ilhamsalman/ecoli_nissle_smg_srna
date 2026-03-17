#!/usr/bin/env python3
import argparse
import pandas as pd

def read_saf(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    need = {"GeneID","Chr","Start","End","Strand"}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"[ERR] SAF missing columns {missing}: {path}")
    df["Start"] = pd.to_numeric(df["Start"], errors="raise").astype(int)
    df["End"]   = pd.to_numeric(df["End"], errors="raise").astype(int)
    df["Strand"]= df["Strand"].astype(str).str.strip()
    return df.rename(columns={"GeneID":"Geneid"})

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sig_csv", required=True, help="voom_*_significant.csv")
    ap.add_argument("--saf", required=True, help="matching SAF used for counts")
    ap.add_argument("--out_bed", required=True)
    args = ap.parse_args()

    sig = pd.read_csv(args.sig_csv, dtype=str)
    for col in ["Geneid","sRNA_types"]:
        if col not in sig.columns:
            raise SystemExit(f"[ERR] Missing column '{col}' in {args.sig_csv}. Found: {list(sig.columns)}")

    # keep ONLY antisense by annotation field
    sig = sig[sig["sRNA_types"].astype(str).str.strip().str.lower() == "antisense"].copy()
    if sig.empty:
        raise SystemExit(f"[ERR] No rows with sRNA_types == Antisense in {args.sig_csv}")

    saf = read_saf(args.saf)

    # attach coordinates
    m = sig.merge(saf[["Geneid","Chr","Start","End","Strand"]], on="Geneid", how="left")
    if m["Chr"].isna().any():
        missing = m[m["Chr"].isna()]["Geneid"].unique().tolist()[:20]
        raise SystemExit(f"[ERR] Some significant antisense sRNAs lack SAF coords (Geneid mismatch). Examples: {missing}")

    # BED is 0-based, half-open: [start-1, end)
    bed = pd.DataFrame({
        "chrom":  m["Chr"],
        "start":  (m["Start"].astype(int) - 1).clip(lower=0),
        "end":    m["End"].astype(int),
        "name":   m["Geneid"],
        # score: keep logFC if available else 0
        "score":  pd.to_numeric(m["logFC"], errors="coerce").fillna(0).round(6),
        "strand": m["Strand"]
    })

    bed.to_csv(args.out_bed, sep="\t", header=False, index=False)
    print(f"[OK] wrote {len(bed)} BED6 rows -> {args.out_bed}")

if __name__ == "__main__":
    main()
