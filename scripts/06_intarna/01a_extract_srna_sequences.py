#!/usr/bin/env python3
"""
Extract selected sRNA sequences for IntaRNA.

Inputs
------
Selected ID files:
    results/sRNA_detection/06_intarna/selected_exp_ids.txt
    results/sRNA_detection/06_intarna/selected_sta_ids.txt

BED files:
    results/sRNA_detection/05_known_srna/srna_exponential.bed
    results/sRNA_detection/05_known_srna/srna_stationary.bed

Genome:
    ref/EcN_CP007799.1.fa

Outputs
-------
    results/sRNA_detection/06_intarna/sequences/srnas/exp/
    results/sRNA_detection/06_intarna/sequences/srnas/sta/
"""

import os
import re
import sys
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    print("ERROR: Biopython is required.")
    print("Install with: conda install -c conda-forge biopython")
    sys.exit(1)


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

GENOME_FASTA = "ref/EcN_CP007799.1.fa"

BED_EXP = (
    "results/sRNA_detection/05_known_srna/"
    "srna_exponential.bed"
)

BED_STA = (
    "results/sRNA_detection/05_known_srna/"
    "srna_stationary.bed"
)

ID_EXP = Path(
    "results/sRNA_detection/06_intarna/"
    "selected_exp_ids.txt"
)

ID_STA = Path(
    "results/sRNA_detection/06_intarna/"
    "selected_sta_ids.txt"
)

OUT_ROOT = Path(
    "results/sRNA_detection/06_intarna/"
    "sequences/srnas"
)

OUT_EXP = OUT_ROOT / "exp"
OUT_STA = OUT_ROOT / "sta"


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def ensure_dirs():
    OUT_EXP.mkdir(
        parents=True,
        exist_ok=True
    )

    OUT_STA.mkdir(
        parents=True,
        exist_ok=True
    )


def load_genome(fasta_path):

    if not os.path.exists(
        fasta_path
    ):
        raise FileNotFoundError(
            f"Genome FASTA not found: {fasta_path}"
        )

    genome = SeqIO.to_dict(
        SeqIO.parse(
            fasta_path,
            "fasta"
        )
    )

    if not genome:
        raise RuntimeError(
            f"No records found in FASTA: {fasta_path}"
        )

    return genome


def read_ids(path):

    if not path.exists():
        raise FileNotFoundError(
            f"ID file not found: {path}"
        )

    ids = []

    for line in path.read_text().splitlines():

        line = line.strip()

        if not line:
            continue

        ids.append(line)

    if not ids:
        raise RuntimeError(
            f"No IDs found in: {path}"
        )

    return ids


def get_numeric_suffix(name):

    match = re.search(
        r"(\d{4,})",
        name
    )

    return (
        match.group(1)
        if match
        else None
    )


def extract_sequence(
    genome,
    chrom,
    start,
    end,
    strand
):

    if chrom not in genome:
        raise KeyError(
            f"Chromosome '{chrom}' "
            "not found in genome FASTA"
        )

    seq = genome[chrom].seq[
        start:end
    ]

    if strand == "-":
        seq = seq.reverse_complement()

    return seq


def write_fasta(
    seq,
    seq_id,
    out_path,
    desc=""
):

    with open(
        out_path,
        "w"
    ) as handle:

        handle.write(
            f">{seq_id} {desc}\n"
        )

        seq_text = str(seq)

        for i in range(
            0,
            len(seq_text),
            60
        ):
            handle.write(
                seq_text[i:i + 60]
                + "\n"
            )


# ------------------------------------------------------------
# BED matching
# ------------------------------------------------------------

def parse_bed(
    bed_path,
    selected_ids
):
    """
    Match selected IDs against BED entries.

    Matching order:
        1. exact ID
        2. numeric suffix
    """

    if not os.path.exists(
        bed_path
    ):
        raise FileNotFoundError(
            f"BED file not found: {bed_path}"
        )

    selected_set = set(
        selected_ids
    )

    suffix_to_id = {}

    for sid in selected_ids:

        suffix = get_numeric_suffix(
            sid
        )

        if suffix:
            suffix_to_id[
                suffix
            ] = sid

    coords = {}

    with open(
        bed_path
    ) as handle:

        for line in handle:

            line = line.strip()

            if (
                not line
                or line.startswith("#")
            ):
                continue

            parts = line.split("\t")

            if len(parts) < 6:
                continue

            (
                chrom,
                start,
                end,
                name,
                score,
                strand
            ) = parts[:6]

            name = name.strip()

            entry = {
                "chrom": chrom.strip(),
                "start": int(start),
                "end": int(end),
                "strand": strand.strip(),
            }

            # Exact ID match
            if name in selected_set:

                coords[
                    name
                ] = entry

                continue

            # Numeric suffix match
            suffix = get_numeric_suffix(
                name
            )

            if (
                suffix
                and suffix
                in suffix_to_id
            ):

                selected_id = (
                    suffix_to_id[
                        suffix
                    ]
                )

                if (
                    selected_id
                    not in coords
                ):
                    coords[
                        selected_id
                    ] = entry

    return coords


# ------------------------------------------------------------
# Extraction
# ------------------------------------------------------------

def extract_phase(
    phase,
    id_file,
    bed_file,
    out_dir,
    genome
):

    print(
        "=" * 70
    )

    print(
        f"{phase.upper()} sRNA EXTRACTION"
    )

    print(
        "=" * 70
    )

    selected_ids = read_ids(
        id_file
    )

    print(
        "Selected IDs:",
        len(selected_ids)
    )

    coords = parse_bed(
        bed_file,
        selected_ids
    )

    print(
        "IDs matched in BED:",
        len(coords)
    )

    missing = []
    generated = []

    for sid in selected_ids:

        if sid not in coords:

            missing.append(
                sid
            )

            continue

        info = coords[
            sid
        ]

        seq = extract_sequence(
            genome,
            info["chrom"],
            info["start"],
            info["end"],
            info["strand"]
        )

        safe_id = sid.replace(
            "__",
            "_"
        )

        outfile = (
            out_dir
            / f"{safe_id}.fasta"
        )

        description = (
            f"{info['chrom']}:"
            f"{info['start']}-"
            f"{info['end']}"
            f"({info['strand']}) "
            f"len={len(seq)}nt "
            f"phase={phase}"
        )

        write_fasta(
            seq,
            safe_id,
            outfile,
            description
        )

        generated.append(
            safe_id
        )

    print(
        "FASTAs generated:",
        len(generated)
    )

    print(
        "Missing:",
        len(missing)
    )

    if missing:

        print(
            "\nMissing IDs:"
        )

        for sid in missing:
            print(
                f"  {sid}"
            )

        raise SystemExit(
            "\nERROR: one or more selected "
            "sRNAs were not found in the BED file."
        )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    ensure_dirs()

    genome = load_genome(
        GENOME_FASTA
    )

    extract_phase(
        phase="EXP",
        id_file=ID_EXP,
        bed_file=BED_EXP,
        out_dir=OUT_EXP,
        genome=genome
    )

    extract_phase(
        phase="STA",
        id_file=ID_STA,
        bed_file=BED_STA,
        out_dir=OUT_STA,
        genome=genome
    )

    print(
        "\nDONE"
    )


if __name__ == "__main__":
    main()
