#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# STEP 2 — Convert DEG tables (with coordinates) to BED6
# Input tables MUST already contain genomic coordinates + strand.
# ============================================================

# --- Adjust paths ---
EXP_TSV="inputs/DEG_exp_with_coordinates.tsv"
STA_TSV="inputs/DEG_sta_with_coordinates.tsv"

EXP_BED="inputs/DEG_exp_with_fc.bed"
STA_BED="inputs/DEG_sta_with_fc.bed"

# Expected columns in DEG_*_with_coordinates.tsv (header row):
# 1 LocusTag   2 Gene(optional)   3 log2FC(or FC score)   4 Function(optional)
# 5 chr        6 start   7 end   8 strand

# --- Exponential ---
awk 'BEGIN{FS="\t"; OFS="\t"}
NR==1{next}
$5!="" && $6!="" && $7!="" && $8!=""{
    chrom=$5; start=$6; end=$7; strand=$8;

    locus=$1;
    gene=$2; if(gene=="") gene="NA";
    name=locus"|"gene;

    fc=$3;
    if(fc=="INF") fc=999;
    else if(fc=="-INF") fc=-999;

    gsub(/\..*/, "", start);
    gsub(/\..*/, "", end);

    if(start>=0 && end>start)
        print chrom, start, end, name, fc, strand;
}' "$EXP_TSV" > "$EXP_BED"

# --- Stationary ---
awk 'BEGIN{FS="\t"; OFS="\t"}
NR==1{next}
$5!="" && $6!="" && $7!="" && $8!=""{
    chrom=$5; start=$6; end=$7; strand=$8;

    locus=$1;
    gene=$2; if(gene=="") gene="NA";
    name=locus"|"gene;

    fc=$3;
    if(fc=="INF") fc=999;
    else if(fc=="-INF") fc=-999;

    gsub(/\..*/, "", start);
    gsub(/\..*/, "", end);

    if(start>=0 && end>start)
        print chrom, start, end, name, fc, strand;
}' "$STA_TSV" > "$STA_BED"

echo "[OK] wrote:"
echo "  $EXP_BED  (lines: $(wc -l < "$EXP_BED"))"
echo "  $STA_BED  (lines: $(wc -l < "$STA_BED"))"
