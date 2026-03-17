#!/usr/bin/env python3

import argparse
import pandas as pd
from collections import defaultdict, deque


def read_saf(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype={"GeneID": str, "Chr": str, "Strand": str})
    required = {"GeneID", "Chr", "Start", "End", "Strand"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {missing}. Found: {list(df.columns)}")

    df["Start"] = pd.to_numeric(df["Start"], errors="raise").astype(int)
    df["End"] = pd.to_numeric(df["End"], errors="raise").astype(int)
    df["Strand"] = df["Strand"].astype(str).str.strip()

    df = df[df["Strand"].isin(["+", "-"])].copy()
    df = df[(df["Start"] > 0) & (df["End"] >= df["Start"])].copy()
    df["len"] = df["End"] - df["Start"] + 1
    return df


def overlap_len(a_start, a_end, b_start, b_end) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def build_overlap_pairs(exp: pd.DataFrame, sta: pd.DataFrame, min_recip: float) -> pd.DataFrame:
    pairs = []
    sta_groups = sta.groupby(["Chr", "Strand"], sort=False)

    for (chr_, strand), sta_g in sta_groups:
        exp_g = exp[(exp["Chr"] == chr_) & (exp["Strand"] == strand)]
        if exp_g.empty:
            continue

        exp_g = exp_g.sort_values("Start").reset_index(drop=True)

        for _, s in sta_g.iterrows():
            s_start, s_end, s_len = int(s["Start"]), int(s["End"]), int(s["len"])
            cand = exp_g[(exp_g["Start"] <= s_end) & (exp_g["End"] >= s_start)]
            if cand.empty:
                continue

            for _, e in cand.iterrows():
                e_start, e_end, e_len = int(e["Start"]), int(e["End"]), int(e["len"])
                ov = overlap_len(s_start, s_end, e_start, e_end)
                if ov <= 0:
                    continue

                recip_sta = ov / s_len
                recip_exp = ov / e_len
                if recip_sta >= min_recip and recip_exp >= min_recip:
                    pairs.append({
                        "Chr": chr_,
                        "Strand": strand,
                        "sta_GeneID": s["GeneID"],
                        "sta_Start": s_start,
                        "sta_End": s_end,
                        "exp_GeneID": e["GeneID"],
                        "exp_Start": e_start,
                        "exp_End": e_end,
                    })

    return pd.DataFrame(pairs)


def connected_components(pairs: pd.DataFrame):
    adj = defaultdict(set)

    for _, r in pairs.iterrows():
        exp_node = "EXP::" + str(r["exp_GeneID"])
        sta_node = "STA::" + str(r["sta_GeneID"])
        adj[exp_node].add(sta_node)
        adj[sta_node].add(exp_node)

    visited = set()
    comps = []

    for node in list(adj.keys()):
        if node in visited:
            continue

        q = deque([node])
        visited.add(node)
        exp_nodes, sta_nodes = set(), set()

        while q:
            cur = q.popleft()
            if cur.startswith("EXP::"):
                exp_nodes.add(cur.replace("EXP::", "", 1))
            else:
                sta_nodes.add(cur.replace("STA::", "", 1))

            for nb in adj[cur]:
                if nb not in visited:
                    visited.add(nb)
                    q.append(nb)

        comps.append({"exp_nodes": exp_nodes, "sta_nodes": sta_nodes})

    return comps


def summarize_coords(ids, df, prefix):
    if not ids:
        return {
            f"{prefix}_n_features": 0,
            f"{prefix}_start_min": pd.NA,
            f"{prefix}_end_max": pd.NA,
            f"{prefix}_GeneIDs": "",
        }
    sub = df[df["GeneID"].isin(ids)]
    return {
        f"{prefix}_n_features": int(sub["GeneID"].nunique()),
        f"{prefix}_start_min": int(sub["Start"].min()),
        f"{prefix}_end_max": int(sub["End"].max()),
        f"{prefix}_GeneIDs": ";".join(sorted(set(sub["GeneID"]))),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp_saf", required=True)
    ap.add_argument("--sta_saf", required=True)
    ap.add_argument("--min_recip", type=float, default=0.5)
    ap.add_argument("--out_prefix", required=True)
    args = ap.parse_args()

    exp = read_saf(args.exp_saf)
    sta = read_saf(args.sta_saf)

    pairs = build_overlap_pairs(exp, sta, args.min_recip)
    out_loci = f"{args.out_prefix}_loci_minrecip_{args.min_recip}.csv"

    if pairs.empty:
        pd.DataFrame([{
            "min_recip": args.min_recip,
            "n_loci": 0,
            "note": "No reciprocal-overlap pairs found; no shared loci."
        }]).to_csv(out_loci, index=False)
        print("No reciprocal-overlap pairs found.")
        print(f"Wrote: {out_loci}")
        return

    comps = connected_components(pairs)

    locus_rows = []
    for i, c in enumerate(comps, start=1):
        locus_id = f"LOCUS_{i:05d}"

        mask = pairs["exp_GeneID"].isin(c["exp_nodes"]) | pairs["sta_GeneID"].isin(c["sta_nodes"])
        chr_strand = pairs.loc[mask, ["Chr", "Strand"]].drop_duplicates()
        chr_ = chr_strand.iloc[0]["Chr"]
        strand = chr_strand.iloc[0]["Strand"]

        exp_sum = summarize_coords(c["exp_nodes"], exp, "exp")
        sta_sum = summarize_coords(c["sta_nodes"], sta, "sta")

        locus_start = min(x for x in [exp_sum["exp_start_min"], sta_sum["sta_start_min"]] if pd.notna(x))
        locus_end   = max(x for x in [exp_sum["exp_end_max"], sta_sum["sta_end_max"]] if pd.notna(x))

        locus_rows.append({
            "locus_id": locus_id,
            "Chr": chr_,
            "Strand": strand,
            "locus_start": int(locus_start),
            "locus_end": int(locus_end),
            "min_recip": args.min_recip,
            "shared_locus": True,
            **exp_sum,
            **sta_sum,
        })

    loci = pd.DataFrame(locus_rows)
    loci.to_csv(out_loci, index=False)

    print("=== Shared sRNA by SAF reciprocal overlap (loci only) ===")
    print(f"min_recip = {args.min_recip}")
    print(f"EXP features total: {len(exp)}")
    print(f"STA features total: {len(sta)}")
    print(f"Shared LOCI (connected components): {len(loci)}")
    print(f"Wrote: {out_loci}")


if __name__ == "__main__":
    main()
