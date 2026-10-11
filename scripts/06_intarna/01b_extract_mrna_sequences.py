#!/usr/bin/env python3

import argparse
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio import SeqIO


def read_table(path, sheet=None):
    path = Path(path)

    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name=sheet)

    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)

    return pd.read_csv(path, sep="\t")


def clean_text(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def load_annotation(path):
    locus_map = defaultdict(list)

    for record in SeqIO.parse(path, "genbank"):
        for feature in record.features:

            if feature.type not in {"CDS", "gene", "pseudogene"}:
                continue

            tags = feature.qualifiers.get("locus_tag", [])
            if not tags:
                continue

            if feature.location.strand == 1:
                strand = "+"
            elif feature.location.strand == -1:
                strand = "-"
            else:
                continue

            locus = tags[0]

            locus_map[locus].append({
                "type": feature.type,
                "start0": int(feature.location.start),
                "end0": int(feature.location.end),
                "strand": strand,
            })

    return locus_map


def best_feature(items):
    priority = {
        "CDS": 0,
        "pseudogene": 1,
        "gene": 2,
    }

    return sorted(
        items,
        key=lambda x: priority.get(x["type"], 99)
    )[0]


def write_fasta(path, header, seq):
    with open(path, "w") as handle:
        handle.write(header + "\n")

        seq = str(seq)

        for i in range(0, len(seq), 70):
            handle.write(seq[i:i + 70] + "\n")


def main():

    parser = argparse.ArgumentParser(
        description="Extract phase-specific mRNA target sequences for IntaRNA."
    )

    parser.add_argument(
        "--table",
        required=True,
        help="Curated TSV/CSV/XLSX target table."
    )

    parser.add_argument(
        "--genome",
        required=True,
        help="Reference genome FASTA."
    )

    parser.add_argument(
        "--annotation",
        required=True,
        help="Reference GenBank annotation."
    )

    parser.add_argument(
        "--phase",
        required=True,
        help="Phase label, e.g. EXP or STA."
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for FASTA files."
    )

    parser.add_argument(
        "--sheet",
        default=None,
        help="Excel sheet name if input table is XLSX."
    )

    parser.add_argument(
        "--upstream",
        type=int,
        default=200,
        help="Upstream sequence length for annotation-derived targets."
    )

    parser.add_argument(
        "--downstream",
        type=int,
        default=500,
        help="Downstream sequence length for annotation-derived targets."
    )

    args = parser.parse_args()

    # ------------------------------------------------------------
    # Load curated table
    # ------------------------------------------------------------

    df = read_table(
        args.table,
        args.sheet
    )

    required_columns = {
        "Locus Tag",
        "Gene",
        "Fold Change",
        "Gene Function",
    }

    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        raise SystemExit(
            f"ERROR: missing required columns: "
            f"{sorted(missing_columns)}"
        )

    # ------------------------------------------------------------
    # Clean key fields
    # ------------------------------------------------------------

    df["Locus Tag"] = (
        df["Locus Tag"]
        .astype(str)
        .str.strip()
    )

    df["Gene"] = (
        df["Gene"]
        .astype(str)
        .str.strip()
    )

    # ------------------------------------------------------------
    # Validate locus tags
    # ------------------------------------------------------------

    duplicate_loci = df[
        df["Locus Tag"].duplicated(keep=False)
    ]

    if len(duplicate_loci):

        print("ERROR: duplicated locus tags:")

        print(
            duplicate_loci[
                [
                    "Locus Tag",
                    "Gene",
                    "Fold Change"
                ]
            ].to_string(index=False)
        )

        raise SystemExit(1)

    # ------------------------------------------------------------
    # Build safe gene-based FASTA names
    #
    # Unique gene:
    #   gadA.fasta
    #
    # Duplicated gene:
    #   priA__ECOLIN_26822.fasta
    # ------------------------------------------------------------

    gene_counts = (
        df["Gene"]
        .value_counts()
        .to_dict()
    )

    def make_target_id(gene, locus):

        gene = clean_text(gene)
        locus = clean_text(locus)

        if not gene:
            raise ValueError(
                f"Missing gene symbol for locus {locus}"
            )

        if gene_counts.get(gene, 0) > 1:
            return f"{gene}__{locus}"

        return gene

    # ------------------------------------------------------------
    # Load genome
    # ------------------------------------------------------------

    genome_record = SeqIO.read(
        args.genome,
        "fasta"
    )

    genome = genome_record.seq
    genome_length = len(genome)

    # ------------------------------------------------------------
    # Determine coordinate mode
    #
    # STA coordinate-ready table has:
    #   Start_0based
    #   End_0based
    #   Strand
    #
    # EXP falls back to GenBank lookup.
    # ------------------------------------------------------------

    coordinate_columns = {
        "Start_0based",
        "End_0based",
        "Strand",
    }

    use_table_coordinates = (
        coordinate_columns
        <= set(df.columns)
    )

    if use_table_coordinates:

        print(
            "Using validated target coordinates "
            "from curated table."
        )

        locus_map = None

    else:

        print(
            "Using GenBank annotation to generate "
            "target coordinates."
        )

        locus_map = load_annotation(
            args.annotation
        )

    # ------------------------------------------------------------
    # Output directory
    # ------------------------------------------------------------

    outdir = Path(args.outdir)

    outdir.mkdir(
        parents=True,
        exist_ok=True
    )

    # Remove old FASTA files before regeneration.
    for fasta in outdir.glob("*.fasta"):
        fasta.unlink()

    manifest = []
    unresolved = []

    # ------------------------------------------------------------
    # Extract sequences
    # ------------------------------------------------------------

    for _, row in df.iterrows():

        locus = clean_text(
            row["Locus Tag"]
        )

        gene = clean_text(
            row["Gene"]
        )

        fc = row["Fold Change"]

        function = clean_text(
            row["Gene Function"]
        )

        target_id = make_target_id(
            gene,
            locus
        )

        # --------------------------------------------------------
        # Use validated coordinates from curated table
        # --------------------------------------------------------

        if use_table_coordinates:

            try:
                start = int(
                    row["Start_0based"]
                )

                end = int(
                    row["End_0based"]
                )

            except Exception:

                unresolved.append({
                    "Locus Tag": locus,
                    "Gene": gene,
                    "Reason": "invalid_coordinate_value",
                })

                continue

            strand = clean_text(
                row["Strand"]
            )

            source = clean_text(
                row.get(
                    "Coordinate_Source",
                    "curated_table"
                )
            )

        # --------------------------------------------------------
        # Otherwise derive coordinates from current annotation
        # --------------------------------------------------------

        else:

            if locus not in locus_map:

                unresolved.append({
                    "Locus Tag": locus,
                    "Gene": gene,
                    "Reason":
                        "locus_not_found_in_annotation",
                })

                continue

            feature = best_feature(
                locus_map[locus]
            )

            strand = feature["strand"]

            if strand == "+":

                start = max(
                    0,
                    feature["start0"]
                    - args.upstream
                )

                end = min(
                    genome_length,
                    feature["start0"]
                    + args.downstream
                )

            else:

                start = max(
                    0,
                    feature["end0"]
                    - args.downstream
                )

                end = min(
                    genome_length,
                    feature["end0"]
                    + args.upstream
                )

            source = "current_annotation"

        # --------------------------------------------------------
        # Validate coordinates
        # --------------------------------------------------------

        if start < 0:

            unresolved.append({
                "Locus Tag": locus,
                "Gene": gene,
                "Reason": f"start_out_of_range:{start}",
            })

            continue

        if end > genome_length:

            unresolved.append({
                "Locus Tag": locus,
                "Gene": gene,
                "Reason": f"end_out_of_range:{end}",
            })

            continue

        if start >= end:

            unresolved.append({
                "Locus Tag": locus,
                "Gene": gene,
                "Reason":
                    f"invalid_interval:{start}-{end}",
            })

            continue

        if strand not in {"+", "-"}:

            unresolved.append({
                "Locus Tag": locus,
                "Gene": gene,
                "Reason":
                    f"invalid_strand:{strand}",
            })

            continue

        # --------------------------------------------------------
        # Extract sequence
        # --------------------------------------------------------

        seq = genome[start:end]

        if strand == "-":
            seq = seq.reverse_complement()

        # --------------------------------------------------------
        # Write gene-based FASTA
        # --------------------------------------------------------

        fasta_path = (
            outdir
            / f"{target_id}.fasta"
        )

        header = (
            f">{target_id} "
            f"locus_tag={locus} "
            f"gene={gene} "
            f"{genome_record.id}:"
            f"{start}-{end}({strand}) "
            f"logFC={fc} "
            f"phase={args.phase} "
            f"len={len(seq)}nt"
        )

        write_fasta(
            fasta_path,
            header,
            seq
        )

        manifest.append({
            "Target_ID": target_id,
            "Locus Tag": locus,
            "Gene": gene,
            "Fold Change": fc,
            "Gene Function": function,
            "Start_0based": start,
            "End_0based": end,
            "Strand": strand,
            "Target_Length": len(seq),
            "Coordinate_Source": source,
            "FASTA": str(fasta_path),
        })

    # ------------------------------------------------------------
    # Save reports
    # ------------------------------------------------------------

    manifest_df = pd.DataFrame(
        manifest
    )

    unresolved_df = pd.DataFrame(
        unresolved
    )

    manifest_path = (
        outdir
        / "mrna_fasta_manifest.tsv"
    )

    unresolved_path = (
        outdir
        / "unresolved_targets.tsv"
    )

    manifest_df.to_csv(
        manifest_path,
        sep="\t",
        index=False
    )

    unresolved_df.to_csv(
        unresolved_path,
        sep="\t",
        index=False
    )

    # ------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------

    print()
    print("=" * 70)
    print(
        f"{args.phase.upper()} mRNA FASTA EXTRACTION"
    )
    print("=" * 70)

    print(
        f"Curated table rows : {len(df)}"
    )

    print(
        f"FASTAs generated   : {len(manifest_df)}"
    )

    print(
        f"Unresolved targets : {len(unresolved_df)}"
    )

    if len(manifest_df):

        print()
        print("Coordinate sources:")

        print(
            manifest_df[
                "Coordinate_Source"
            ]
            .value_counts()
            .to_string()
        )

        print()
        print("Target length summary:")

        print(
            manifest_df[
                "Target_Length"
            ]
            .describe()
            .to_string()
        )

        duplicate_gene_targets = (
            manifest_df[
                manifest_df[
                    "Target_ID"
                ].str.contains(
                    "__ECOLIN_",
                    regex=False
                )
            ]
        )

        print()
        print(
            "Gene-name collisions resolved "
            f"with locus tag: "
            f"{len(duplicate_gene_targets)}"
        )

        if len(duplicate_gene_targets):

            print(
                duplicate_gene_targets[
                    [
                        "Target_ID",
                        "Locus Tag",
                        "Gene"
                    ]
                ]
                .to_string(index=False)
            )

    print()
    print(
        f"Manifest: {manifest_path}"
    )

    print(
        f"Unresolved report: "
        f"{unresolved_path}"
    )

    if len(unresolved_df):

        print()
        print("ERROR: unresolved targets remain:")

        print(
            unresolved_df.to_string(
                index=False
            )
        )

        raise SystemExit(1)

    print()
    print("DONE")


if __name__ == "__main__":
    main()
