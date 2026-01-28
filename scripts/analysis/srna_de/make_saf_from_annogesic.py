#!/usr/bin/env python3
import argparse
import pandas as pd

def read_table(path: str) -> pd.DataFrame:
    # try tab first, then delimiter auto-detect
    try:
        return pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return pd.read_csv(path, sep=None, engine="python", dtype=str)

def pick_col(df, candidates):
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    raise ValueError(f"Could not find any of {candidates} in columns: {list(df.columns)}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True)
    ap.add_argument("--outfile", required=True)
    ap.add_argument("--phase_label", required=True, help="exp or sta (prepended to GeneID)")
    ap.add_argument("--chr_value", default=None, help="Force Chr/contig name if table lacks it")
    args = ap.parse_args()

    df = read_table(args.infile)

    name_col   = pick_col(df, ["Name", "ID"])
    start_col  = pick_col(df, ["Start"])
    end_col    = pick_col(df, ["End"])
    strand_col = pick_col(df, ["Strand"])

    chr_col = None
    if args.chr_value is None:
        for cand in ["Genome", "Chr", "Chromosome", "Contig", "seqname"]:
            if cand.lower() in [c.lower() for c in df.columns]:
                chr_col = pick_col(df, [cand])
                break
    if chr_col is None and args.chr_value is None:
        raise ValueError("No chromosome/contig column found. Provide --chr_value <contig_name>.")

    saf = pd.DataFrame()
    saf["GeneID"] = df[name_col].astype(str).str.strip()

    # Make IDs unique if duplicates exist
    coord_suffix = (
        df[start_col].astype(str).str.strip() + "_" +
        df[end_col].astype(str).str.strip() + "_" +
        df[strand_col].astype(str).str.strip()
    )
    saf["GeneID"] = args.phase_label + "__" + saf["GeneID"]
    dup = saf["GeneID"].duplicated(keep=False)
    saf.loc[dup, "GeneID"] = saf.loc[dup, "GeneID"] + "__" + coord_suffix[dup]

    saf["Chr"] = args.chr_value if args.chr_value is not None else df[chr_col].astype(str).str.strip()
    saf["Start"] = pd.to_numeric(df[start_col], errors="raise").astype(int)
    saf["End"]   = pd.to_numeric(df[end_col], errors="raise").astype(int)
    saf["Strand"]= df[strand_col].astype(str).str.strip()

    # Basic sanity filters
    saf = saf[saf["Strand"].isin(["+", "-"])].copy()
    saf = saf[(saf["Start"] > 0) & (saf["End"] >= saf["Start"])].copy()

    saf.to_csv(args.outfile, sep="\t", index=False)
    print(f"Wrote {len(saf)} rows to {args.outfile}")

if __name__ == "__main__":
    main()
