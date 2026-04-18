#!/usr/bin/env python3
"""
Extract intergenic sRNA candidate sequences for IntaRNA.

Inputs
------
Candidate directory:
    results/sRNA_detection/04_srna_de_filtered/

Expected candidate CSV files:
    voom_exp_s2_candidates.csv
    voom_sta_s2_candidates.csv

BED files:
    results/sRNA_detection/05_known_srna/srna_exponential.bed
    results/sRNA_detection/05_known_srna/srna_stationary.bed

Genome FASTA:
    ref/EcN_CP007799.1.fa

Outputs
-------
FASTA files:
    results/sRNA_detection/06_intarna/sequences/srnas/exp/
    results/sRNA_detection/06_intarna/sequences/srnas/sta/

ID lists:
    results/sRNA_detection/06_intarna/intergenic_exp_ids.txt
    results/sRNA_detection/06_intarna/intergenic_sta_ids.txt

Run from project root:
    python3 scripts/06_intarna/01_extract_sequences.py
"""

import csv
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

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------

GENOME_FASTA = "ref/EcN_CP007799.1.fa"

CANDIDATE_DIR = "results/sRNA_detection/04_srna_de_filtered"
EXP_CSV = f"{CANDIDATE_DIR}/voom_exp_s2_candidates.csv"
STA_CSV = f"{CANDIDATE_DIR}/voom_sta_s2_candidates.csv"

BED_EXP = "results/sRNA_detection/05_known_srna/srna_exponential.bed"
BED_STA = "results/sRNA_detection/05_known_srna/srna_stationary.bed"

OUT_ROOT = Path("results/sRNA_detection/06_intarna/sequences/srnas")
OUT_EXP = OUT_ROOT / "exp"
OUT_STA = OUT_ROOT / "sta"

ID_EXP = Path("results/sRNA_detection/06_intarna/intergenic_exp_ids.txt")
ID_STA = Path("results/sRNA_detection/06_intarna/intergenic_sta_ids.txt")


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def ensure_dirs():
    OUT_EXP.mkdir(parents=True, exist_ok=True)
    OUT_STA.mkdir(parents=True, exist_ok=True)
    ID_EXP.parent.mkdir(parents=True, exist_ok=True)


def load_genome(fasta_path):
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Genome FASTA not found: {fasta_path}")
    genome = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    if not genome:
        raise RuntimeError(f"No records found in FASTA: {fasta_path}")
    return genome


def resolve_default_chrom(genome):
    if "CP007799.1" in genome:
        return "CP007799.1"
    for key in genome:
        if "CP007799" in key or "EcN" in key:
            return key
    return next(iter(genome.keys()))


def extract_sequence(genome, chrom, start, end, strand):
    if chrom not in genome:
        raise KeyError(f"Chromosome '{chrom}' not found in genome FASTA")
    seq = genome[chrom].seq[start:end]
    if strand == "-":
        seq = seq.reverse_complement()
    return seq


def write_fasta(seq, seq_id, out_path, desc=""):
    with open(out_path, "w") as fh:
        fh.write(f">{seq_id} {desc}\n")
        seq_str = str(seq)
        for i in range(0, len(seq_str), 60):
            fh.write(seq_str[i:i+60] + "\n")


def get_numeric_suffix(name):
    m = re.search(r"(\d{4,})", name)
    return m.group(1) if m else None


# ---------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------

def parse_candidate_csv(csv_path):
    """
    Keep only Intergenic candidates.
    Requires at least:
        Geneid, sRNA_types
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Candidate CSV not found: {csv_path}")

    rows = []
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = set(reader.fieldnames or [])
        required = {"Geneid", "sRNA_types"}
        missing = required - fieldnames
        if missing:
            raise RuntimeError(
                f"{csv_path} missing required columns: {sorted(missing)}"
            )

        for row in reader:
            srna_type = row["sRNA_types"].strip().strip('"')
            if srna_type != "Intergenic":
                continue

            rows.append({
                "id": row["Geneid"].strip().strip('"'),
                "logFC": row.get("logFC", "").strip(),
                "AveExpr": row.get("AveExpr", "").strip(),
                "adjP": row.get("adj.P.Val", "").strip(),
                "direction": row.get("direction", "").strip().strip('"'),
                "phase": row.get("phase", "").strip().strip('"'),
            })

    return rows


def parse_srna_bed(bed_path, target_ids):
    """
    BED expected format:
        chrom  start  end  name  score  strand

    Matching:
        1. exact ID match
        2. numeric suffix match
    """
    if not os.path.exists(bed_path):
        raise FileNotFoundError(f"sRNA BED not found: {bed_path}")

    target_ids = [x.strip().strip('"') for x in target_ids]
    target_set = set(target_ids)

    suffix_to_target = {}
    for tid in target_ids:
        sfx = get_numeric_suffix(tid)
        if sfx:
            suffix_to_target[sfx] = tid

    coords = {}
    bed_count = 0

    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 6:
                continue

            chrom, start, end, name, score, strand = parts[:6]
            name = name.strip()
            bed_count += 1

            entry = {
                "chrom": chrom.strip(),
                "start": int(start),
                "end": int(end),
                "strand": strand.strip(),
            }

            if name in target_set:
                coords[name] = entry
                continue

            sfx = get_numeric_suffix(name)
            if sfx and sfx in suffix_to_target:
                original_id = suffix_to_target[sfx]
                if original_id not in coords:
                    coords[original_id] = entry

    print(f"  BED entries read: {bed_count}")
    print(f"  Candidate IDs matched: {len(coords)}")
    return coords


# ---------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------

def extract_phase(phase_name, csv_path, bed_path, out_dir, id_list_path, genome, default_chrom):
    print("=" * 80)
    print(f"PHASE: {phase_name}")
    print("=" * 80)
    print(f"CSV : {csv_path}")
    print(f"BED : {bed_path}")
    print(f"OUT : {out_dir}")

    candidates = parse_candidate_csv(csv_path)
    print(f"  Intergenic candidates in CSV: {len(candidates)}")

    if not candidates:
        id_list_path.write_text("")
        print("  No intergenic candidates found.")
        return

    coords = parse_srna_bed(bed_path, [c["id"] for c in candidates])

    extracted = []
    missing = []

    for c in candidates:
        sid = c["id"]

        if sid not in coords:
            missing.append(sid)
            continue

        info = coords[sid]
        chrom = info["chrom"] if info["chrom"] else default_chrom
        seq = extract_sequence(genome, chrom, info["start"], info["end"], info["strand"])

        safe_id = sid.replace("__", "_")
        out_file = out_dir / f"{safe_id}.fasta"

        desc_parts = [
            f"{chrom}:{info['start']}-{info['end']}({info['strand']})",
            f"len={len(seq)}nt",
            "type=Intergenic"
        ]
        if c["logFC"]:
            desc_parts.append(f"logFC={c['logFC']}")
        if c["AveExpr"]:
            desc_parts.append(f"AveExpr={c['AveExpr']}")
        if c["adjP"]:
            desc_parts.append(f"adjP={c['adjP']}")
        if c["direction"]:
            desc_parts.append(f"direction={c['direction']}")

        write_fasta(seq, safe_id, out_file, "  ".join(desc_parts))
        extracted.append(safe_id)

    id_list_path.write_text("\n".join(extracted) + ("\n" if extracted else ""))

    print(f"  FASTA written: {len(extracted)}")
    print(f"  ID list written: {id_list_path}")

    if missing:
        print(f"  WARNING: {len(missing)} candidate IDs not found in BED:")
        for mid in missing:
            print(f"    - {mid}")


def main():
    ensure_dirs()

    genome = load_genome(GENOME_FASTA)
    default_chrom = resolve_default_chrom(genome)

    print(f"Genome loaded: {GENOME_FASTA}")
    print(f"Default chromosome: {default_chrom}")

    extract_phase(
        phase_name="exponential",
        csv_path=EXP_CSV,
        bed_path=BED_EXP,
        out_dir=OUT_EXP,
        id_list_path=ID_EXP,
        genome=genome,
        default_chrom=default_chrom,
    )

    extract_phase(
        phase_name="stationary",
        csv_path=STA_CSV,
        bed_path=BED_STA,
        out_dir=OUT_STA,
        id_list_path=ID_STA,
        genome=genome,
        default_chrom=default_chrom,
    )

    print("=" * 80)
    print("DONE")
    print("=" * 80)
    print(f"FASTA root: {OUT_ROOT}")


if __name__ == "__main__":
    main()