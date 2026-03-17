#!/usr/bin/env python3
"""
Extract sRNA and mRNA sequences for IntaRNA interaction prediction.

Priorities:
  P1 — Novel leuC antisense sRNAs (exp_sRNA_00170, exp_sRNA_00172)
         sRNA coords : hardcoded
         mRNA targets: leuC (hardcoded)

  P2 — Novel intergenic sRNAs, exponential phase (exp__sRNA_02152)
         sRNA coords : looked up in srna_exponential.bed
         mRNA targets: ALL genes in DEG_exp_with_fc.bed

  P3 — Novel intergenic sRNAs, stationary phase (sta__sRNA_00458 etc.)
         sRNA coords : looked up in srna_stationary.bed
         mRNA targets: ALL genes in DEG_sta_with_fc.bed

Run from project root:
    conda activate intarna_env   (or ncrna_env)
    python3 scripts/07_intarna/01_extract_sequences.py
"""

import os
import re
import sys
import csv
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    print("ERROR: Biopython not installed. Run: conda install -c conda-forge biopython")
    sys.exit(1)

# =============================================================================
# Configuration — edit paths here if your layout differs
# =============================================================================

GENOME_FASTA = "ref/EcN_CP007799.1.fa"

# DE result CSVs — used to identify intergenic sRNA IDs
DE_EXP_CSV = "results/sRNA_detection/04_srna_de_filtered/voom_exp_s1_significant.csv"
DE_STA_CSV = "results/sRNA_detection/04_srna_de_filtered/voom_sta_s1_significant.csv"

# ANNOgesic sRNA BED files — source of intergenic sRNA genomic coordinates
# These are searched automatically if the primary path is not found (see find_bed_file())
BED_SRNA_EXP = "results/sRNA_detection/05_known_srna/srna_exponential.bed"
BED_SRNA_STA = "results/sRNA_detection/05_known_srna/srna_stationary.bed"

# DEG BED files — source of mRNA target coordinates
# Format: chrom  start  end  LOCUS|gene  logFC  strand
BED_DEG_EXP = "results/sRNA_detection/06_srna_deg/inputs/DEG_exp_with_fc.bed"
BED_DEG_STA = "results/sRNA_detection/06_srna_deg/inputs/DEG_sta_with_fc.bed"

# Output directories
SRNA_OUT = "results/sRNA_detection/07_intarna/sequences/srnas"
MRNA_OUT = "results/sRNA_detection/07_intarna/sequences/mrnas"

# mRNA extraction window: 200 nt upstream of CDS + first 500 nt of CDS
UTR5_EXTEND  = 200
CDS_WINDOW   = 500

# logFC placeholder values to skip (genes with no quantifiable FC)
FC_SKIP_VALUES = {999.0, -999.0}

# =============================================================================
# Priority 1: Novel leuC antisense sRNAs — hardcoded coordinates
# =============================================================================

PRIORITY1_SRNAS = {
    'exp_sRNA_00170': {'chr': 'CP007799.1', 'start': 356892, 'end': 357211, 'strand': '+'},
    'exp_sRNA_00172': {'chr': 'CP007799.1', 'start': 357086, 'end': 357292, 'strand': '+'},
}

PRIORITY1_MRNAS = {
    'leuC': {
        'chr': 'CP007799.1', 'start': 354739, 'end': 356095, 'strand': '-',
        'locus': 'ECOLIN_00365', 'note': 'leucine biosynthesis — P1 antisense target',
    },
}

# =============================================================================
# Genome utilities
# =============================================================================

def load_genome(fasta_file):
    print(f"Loading genome: {fasta_file}")
    if not os.path.exists(fasta_file):
        print(f"ERROR: not found: {fasta_file}")
        sys.exit(1)
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    print(f"  Chromosomes: {list(genome.keys())}")
    return genome


def resolve_chrom(genome):
    if 'CP007799.1' in genome:
        return 'CP007799.1'
    for k in genome:
        if 'CP007799' in k or 'EcN' in k:
            print(f"  Chromosome resolved to: {k}")
            return k
    print("WARNING: CP007799 not found — using first chromosome.")
    return list(genome.keys())[0]


def extract_sequence(genome, chrom, start, end, strand):
    if chrom not in genome:
        print(f"  WARNING: chromosome '{chrom}' not in genome.")
        return None
    seq = genome[chrom].seq[start:end]
    if strand == '-':
        seq = seq.reverse_complement()
    return seq


def save_fasta(sequence, name, output_file, description=""):
    with open(output_file, 'w') as f:
        f.write(f">{name} {description}\n")
        s = str(sequence)
        for i in range(0, len(s), 60):
            f.write(s[i:i+60] + "\n")


def extract_mrna_region(genome, chrom, start, end, strand):
    """5'UTR window + first CDS_WINDOW nt of CDS."""
    if strand == '+':
        ext_start = max(0, start - UTR5_EXTEND)
        ext_end   = min(end, start + CDS_WINDOW)
    else:
        ext_end   = end + UTR5_EXTEND
        ext_start = max(start, end - CDS_WINDOW)
    seq = extract_sequence(genome, chrom, ext_start, ext_end, strand)
    return seq, ext_start, ext_end


def find_bed_file(primary_path, phase_keyword):
    """
    Try to locate a sRNA BED file. Checks:
      1. The primary path as configured
      2. Common alternative locations under results/sRNA_detection/
      3. Any *.bed file whose name contains the phase keyword

    Returns the resolved path string, or None if not found.
    """
    # 1. Configured path
    if os.path.exists(primary_path):
        return primary_path

    print(f"  INFO: '{primary_path}' not found — searching for alternatives...")

    # 2. Common alternative locations
    search_roots = [
        "results/sRNA_detection",
        "results",
        ".",
    ]
    phase_variants = [phase_keyword, phase_keyword.replace('exponential', 'exp'),
                      phase_keyword.replace('stationary', 'sta')]

    for root in search_roots:
        if not os.path.isdir(root):
            continue
        for bed_file in Path(root).rglob("*.bed"):
            name_lower = bed_file.name.lower()
            if any(v in name_lower for v in phase_variants):
                # Skip DEG BED files
                if 'deg' in name_lower or 'with_fc' in name_lower:
                    continue
                print(f"  Found alternative: {bed_file}")
                return str(bed_file)

    print(f"  ERROR: Could not locate sRNA BED for phase '{phase_keyword}'.")
    print(f"  Expected path : {primary_path}")
    print(f"  Searched under: {search_roots}")
    print(f"  Fix: copy your ANNOgesic sRNA BED to the expected path, or update")
    print(f"       BED_SRNA_EXP / BED_SRNA_STA at the top of this script.")
    print()
    print(f"  To find your BED file, run:")
    print(f"    find results/ -name '*.bed' | grep -v DEG | grep -v with_fc")
    return None


def create_output_dirs():
    for d in [SRNA_OUT, MRNA_OUT]:
        Path(d).mkdir(parents=True, exist_ok=True)
        print(f"  {d}/")


# =============================================================================
# BED file parsers
# =============================================================================

def parse_srna_bed(bed_file, target_ids):
    """
    Parse an ANNOgesic sRNA BED file.
    BED columns: chrom  start  end  name  score  strand
    Returns {original_target_id: {start, end, strand}}.

    Uses fuzzy matching to handle ID format differences between the DE CSV
    (e.g. exp__sRNA_02152) and the BED file (e.g. sRNA_02152, exp_sRNA_02152).
    Matching is done on the numeric suffix (e.g. 02152).
    """
    coords = {}
    if not os.path.exists(bed_file):
        print(f"  WARNING: sRNA BED not found: {bed_file}")
        return coords

    # Build a lookup: numeric suffix -> original target_id
    def get_numeric_suffix(name):
        m = re.search(r'(\d{4,})', name)
        return m.group(1) if m else None

    suffix_to_target = {}
    for t in target_ids:
        t_clean = t.strip().strip('"')
        sfx = get_numeric_suffix(t_clean)
        if sfx:
            suffix_to_target[sfx] = t_clean

    # Also keep exact match set for direct lookup
    target_set = {t.strip().strip('"') for t in target_ids}

    # First pass: collect all BED entries
    all_bed_names = []
    bed_entries = {}
    with open(bed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 6:
                continue
            chrom, start, end, name, score, strand = parts[:6]
            name = name.strip()
            all_bed_names.append(name)
            bed_entries[name] = {
                'start': int(start), 'end': int(end), 'strand': strand.strip()
            }

    # Show sample BED IDs on first call (diagnostic)
    if all_bed_names:
        sample = all_bed_names[:5]
        print(f"  BED file: {bed_file}  ({len(all_bed_names)} entries)")
        print(f"  Sample BED IDs: {sample}")

    # Try exact match first, then suffix-based fuzzy match
    matched_direct = 0
    matched_fuzzy  = 0

    for name, entry in bed_entries.items():
        if name in target_set:
            coords[name] = entry
            matched_direct += 1
        else:
            sfx = get_numeric_suffix(name)
            if sfx and sfx in suffix_to_target:
                original_id = suffix_to_target[sfx]
                if original_id not in coords:   # don't overwrite exact match
                    coords[original_id] = entry
                    matched_fuzzy += 1

    if matched_direct or matched_fuzzy:
        print(f"  Matched: {matched_direct} exact + {matched_fuzzy} fuzzy = {matched_direct+matched_fuzzy} total")

    return coords


def parse_deg_bed(bed_file):
    """
    Parse a DEG BED file with format:
        chrom  start  end  LOCUS|gene  logFC  strand

    Returns list of dicts with keys:
        chrom, start, end, locus, gene, logFC, strand, safe_name

    Skips entries where logFC is a placeholder (999 / -999).
    """
    entries = []
    if not os.path.exists(bed_file):
        print(f"  WARNING: DEG BED not found: {bed_file}")
        return entries

    with open(bed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 6:
                continue

            chrom, start, end, name_field, score_field, strand = parts[:6]

            # Parse logFC — skip placeholders
            try:
                logFC = float(score_field)
            except ValueError:
                logFC = None
            if logFC in FC_SKIP_VALUES:
                continue

            # Parse LOCUS|gene name field
            if '|' in name_field:
                locus, gene = name_field.split('|', 1)
            else:
                locus = name_field
                gene  = name_field

            # Build a filesystem-safe name: prefer gene symbol, fall back to locus
            raw_name = gene.strip() if gene.strip() else locus.strip()
            # Sanitise: replace characters that break filenames
            safe_name = raw_name.replace('/', '_').replace(' ', '_').replace('|', '_')
            if not safe_name:
                safe_name = locus.strip().replace(' ', '_')

            entries.append({
                'chrom':     chrom.strip(),
                'start':     int(start),
                'end':       int(end),
                'locus':     locus.strip(),
                'gene':      gene.strip(),
                'safe_name': safe_name,
                'logFC':     logFC,
                'strand':    strand.strip(),
            })

    return entries


# =============================================================================
# DE CSV parser — intergenic sRNA IDs
# =============================================================================

def parse_de_csv_for_intergenic(csv_file):
    """Return list of dicts for all Intergenic sRNAs in a DE CSV."""
    results = []
    if not os.path.exists(csv_file):
        print(f"  WARNING: DE CSV not found: {csv_file}")
        return results
    with open(csv_file, newline='') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row['sRNA_types'].strip().strip('"') == 'Intergenic':
                results.append({
                    'id':        row['Geneid'].strip().strip('"'),
                    'logFC':     float(row['logFC']),
                    'adjP':      float(row['adj.P.Val']),
                    'direction': row['direction'].strip().strip('"'),
                    'phase':     row['phase'].strip().strip('"'),
                })
    return results


# =============================================================================
# Extraction helpers
# =============================================================================

def extract_intergenic_srnas(genome, chrom, phase, de_csv, srna_bed):
    """
    Extract FASTA for all significant intergenic sRNAs in one phase.
    Returns list of safe names that were successfully extracted.
    """
    print(f"  DE CSV  : {de_csv}")
    intergenic = parse_de_csv_for_intergenic(de_csv)

    if not intergenic:
        print(f"  No intergenic sRNAs found.")
        return []

    print(f"  Found {len(intergenic)} intergenic sRNA(s):")
    for s in intergenic:
        print(f"    {s['id']}   logFC={s['logFC']:+.3f}   FDR={s['adjP']:.2e}   {s['direction']}")

    # Auto-locate BED file if the configured path is missing
    resolved_bed = find_bed_file(srna_bed, phase)
    if resolved_bed is None:
        return []

    coords    = parse_srna_bed(resolved_bed, [s['id'] for s in intergenic])
    extracted = []
    missing   = []

    for srna in intergenic:
        sid = srna['id']
        if sid not in coords:
            print(f"  ⚠ Not in BED: {sid}")
            missing.append(sid)
            continue

        c         = coords[sid]
        seq       = extract_sequence(genome, chrom, c['start'], c['end'], c['strand'])
        if seq is None:
            missing.append(sid)
            continue

        safe_name = sid.replace('__', '_')    # exp__sRNA_02152 → exp_sRNA_02152
        out       = f"{SRNA_OUT}/{safe_name}.fasta"
        desc      = (f"{chrom}:{c['start']}-{c['end']}({c['strand']}) "
                     f"len={len(seq)}nt  phase={phase}  "
                     f"logFC={srna['logFC']:+.3f}  {srna['direction']}")
        save_fasta(seq, safe_name, out, desc)
        print(f"  ✓ {out}  ({len(seq)} nt)")
        extracted.append(safe_name)

    if missing:
        print(f"\n  ⚠ Could not extract {len(missing)} sRNA(s): {missing}")
        print( "    Check that the sRNA BED file contains these IDs.")

    return extracted


def extract_deg_mrnas(genome, chrom, phase, deg_bed, mrna_subdir=None):
    """
    Extract mRNA regions for all DEGs in a BED file.
    Saves to MRNA_OUT (or MRNA_OUT/mrna_subdir if provided).
    Returns list of safe_names extracted.
    """
    print(f"  DEG BED : {deg_bed}")
    entries = parse_deg_bed(deg_bed)

    if not entries:
        print(f"  No DEG entries found.")
        return []

    print(f"  Found {len(entries)} DEG mRNA target(s)")

    out_dir = Path(MRNA_OUT) if mrna_subdir is None else Path(MRNA_OUT) / mrna_subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    extracted  = []
    skipped    = []
    name_seen  = {}    # track duplicates: safe_name → count

    for gene in entries:
        # Handle duplicate gene names (e.g. kfiB appears twice)
        safe_name = gene['safe_name']
        if safe_name in name_seen:
            name_seen[safe_name] += 1
            safe_name = f"{safe_name}_{name_seen[safe_name]}"
        else:
            name_seen[gene['safe_name']] = 1

        seq, s, e = extract_mrna_region(
            genome, chrom,
            gene['start'], gene['end'], gene['strand']
        )
        if seq is None:
            skipped.append(safe_name)
            continue

        out  = str(out_dir / f"{safe_name}.fasta")
        fc_s = f"{gene['logFC']:+.2f}" if gene['logFC'] is not None else "NA"
        desc = (f"{chrom}:{s}-{e}({gene['strand']}) "
                f"{gene['locus']}|{gene['gene']}  "
                f"logFC={fc_s}  phase={phase}  len={len(seq)}nt")
        save_fasta(seq, safe_name, out, desc)
        extracted.append(safe_name)

    print(f"  ✓ Extracted {len(extracted)} mRNA targets  →  {out_dir}/")
    if skipped:
        print(f"  ⚠ Skipped {len(skipped)}: {skipped}")

    return extracted


# =============================================================================
# Main
# =============================================================================

def main():
    print("=" * 80)
    print("sRNA / mRNA Sequence Extraction for IntaRNA")
    print("=" * 80)
    print()

    print("Creating output directories...")
    create_output_dirs()
    print()

    genome = load_genome(GENOME_FASTA)
    chrom  = resolve_chrom(genome)
    print()

    # =========================================================================
    # Priority 1 — Novel leuC antisense sRNAs vs leuC mRNA
    # =========================================================================

    print("=" * 80)
    print("PRIORITY 1 — Novel leuC Antisense sRNAs")
    print("=" * 80)
    print()

    for name, info in PRIORITY1_SRNAS.items():
        print(f"  Extracting sRNA: {name}  [{info['start']}-{info['end']} ({info['strand']})]")
        seq = extract_sequence(genome, chrom, info['start'], info['end'], info['strand'])
        if seq:
            out  = f"{SRNA_OUT}/{name}.fasta"
            desc = (f"{chrom}:{info['start']}-{info['end']}({info['strand']}) "
                    f"len={len(seq)}nt  novel antisense leuC")
            save_fasta(seq, name, out, desc)
            print(f"  ✓ {out}  ({len(seq)} nt)")
    print()

    for name, info in PRIORITY1_MRNAS.items():
        print(f"  Extracting mRNA: {name}  [{info['start']}-{info['end']} ({info['strand']})]")
        seq, s, e = extract_mrna_region(genome, chrom, info['start'], info['end'], info['strand'])
        if seq:
            out  = f"{MRNA_OUT}/{name}.fasta"
            desc = (f"{chrom}:{s}-{e}({info['strand']}) "
                    f"{info['locus']}  {info['note']}  len={len(seq)}nt")
            save_fasta(seq, name, out, desc)
            print(f"  ✓ {out}  ({len(seq)} nt)")
    print()

    # =========================================================================
    # Priority 2 — Novel intergenic sRNAs, exponential phase
    #              Targets: all DEG exp mRNAs
    # =========================================================================

    print("=" * 80)
    print("PRIORITY 2 — Novel Intergenic sRNAs  (exponential phase)")
    print("=" * 80)
    print()

    print("  [sRNAs]")
    exp_intergenic = extract_intergenic_srnas(
        genome, chrom,
        phase    = 'exponential',
        de_csv   = DE_EXP_CSV,
        srna_bed = BED_SRNA_EXP,
    )
    print()

    print("  [mRNA targets — DEG exp]")
    exp_mrna_targets = extract_deg_mrnas(
        genome, chrom,
        phase    = 'exponential',
        deg_bed  = BED_DEG_EXP,
        mrna_subdir = 'deg_exp',
    )
    print()

    # =========================================================================
    # Priority 3 — Novel intergenic sRNAs, stationary phase
    #              Targets: all DEG sta mRNAs
    # =========================================================================

    print("=" * 80)
    print("PRIORITY 3 — Novel Intergenic sRNAs  (stationary phase)")
    print("=" * 80)
    print()

    print("  [sRNAs]")
    sta_intergenic = extract_intergenic_srnas(
        genome, chrom,
        phase    = 'stationary',
        de_csv   = DE_STA_CSV,
        srna_bed = BED_SRNA_STA,
    )
    print()

    print("  [mRNA targets — DEG sta]")
    sta_mrna_targets = extract_deg_mrnas(
        genome, chrom,
        phase    = 'stationary',
        deg_bed  = BED_DEG_STA,
        mrna_subdir = 'deg_sta',
    )
    print()

    # =========================================================================
    # Write intergenic ID lists (consumed by 02_run_intarna.sh)
    # =========================================================================

    id_dir = Path("results/sRNA_detection/07_intarna")
    id_dir.mkdir(parents=True, exist_ok=True)

    exp_id_file = id_dir / "intergenic_exp_ids.txt"
    sta_id_file = id_dir / "intergenic_sta_ids.txt"

    exp_id_file.write_text('\n'.join(exp_intergenic) + ('\n' if exp_intergenic else ''))
    sta_id_file.write_text('\n'.join(sta_intergenic) + ('\n' if sta_intergenic else ''))

    print(f"✓ Intergenic sRNA ID lists:")
    print(f"    {exp_id_file}  ({len(exp_intergenic)} sRNA(s))")
    print(f"    {sta_id_file}  ({len(sta_intergenic)} sRNA(s))")
    print()

    # =========================================================================
    # Summary
    # =========================================================================

    srna_files = sorted(Path(SRNA_OUT).glob("*.fasta"))
    all_intergenic = set(exp_intergenic + sta_intergenic)

    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    print(f"sRNA sequences:  {len(srna_files)}")
    for f in srna_files:
        tag = "  [intergenic]" if f.stem in all_intergenic else "  [antisense]"
        print(f"    {f.name}{tag}")

    print()
    mrna_root = Path(MRNA_OUT)
    total_mrna = sum(1 for _ in mrna_root.rglob("*.fasta"))
    print(f"mRNA sequences:  {total_mrna}")
    print(f"    {mrna_root}/leuC.fasta  (P1 antisense target)")
    if exp_mrna_targets:
        print(f"    {mrna_root}/deg_exp/   ({len(exp_mrna_targets)} DEG exp targets)")
    if sta_mrna_targets:
        print(f"    {mrna_root}/deg_sta/   ({len(sta_mrna_targets)} DEG sta targets)")

    print()
    print("Next step:  bash scripts/07_intarna/02_run_intarna.sh")
    print("=" * 80)


if __name__ == "__main__":
    main()
