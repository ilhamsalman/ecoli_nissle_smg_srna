#!/usr/bin/env bash
set -euo pipefail
# Usage: 02_extract_fasta.sh <genome.fa> <in.bed> <out.fa>
GENOME=$1
BED=$2
OUT=$3
bedtools getfasta -s -name -fi "$GENOME" -bed "$BED" -fo "$OUT"
echo "Wrote: $OUT"
