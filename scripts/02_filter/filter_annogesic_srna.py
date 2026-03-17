#!/usr/bin/env python3
import argparse
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_tsv", required=True, help="ANNOgesic sRNA table (tab-delimited)")
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--min_len", type=int, default=50)
    ap.add_argument("--max_len", type=int, default=500)
    ap.add_argument("--cov_quantile", type=float, default=0.25,
                    help="Quantile of Best_avg_coverage used as phase-specific min coverage")
    ap.add_argument("--min_samples", type=int, default=2,
                    help="Require >=min_samples with Avg_coverage >= min_cov")
    args = ap.parse_args()

    df = pd.read_csv(args.in_tsv, sep="\t")
    cov_cols = [c for c in df.columns if c.startswith("Avg_coverage:")]

    # Make sure coverage is numeric (sometimes comes in as strings)
    df["Best_avg_coverage"] = pd.to_numeric(df["Best_avg_coverage"], errors="coerce")
    for c in cov_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Length filter
    mask_len = df["sRNA_length"].between(args.min_len, args.max_len)

    # Phase-specific coverage threshold
    min_cov = float(df["Best_avg_coverage"].quantile(args.cov_quantile))

    # Reproducibility across samples
    mask_cov = (df[cov_cols] >= min_cov).sum(axis=1) >= args.min_samples

    out = df[mask_len & mask_cov].copy()

    # Save + print a tiny log to stdout
    print(f"[filter] input={len(df)} kept={len(out)} "
          f"min_len={args.min_len} max_len={args.max_len} "
          f"min_cov(q={args.cov_quantile})={min_cov:.4f} min_samples={args.min_samples}")

    out.to_csv(args.out_tsv, sep="\t", index=False)

if __name__ == "__main__":
    main()
