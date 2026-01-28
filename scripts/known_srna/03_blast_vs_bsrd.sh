#!/usr/bin/env bash
set -euo pipefail
# Usage: 03_blast_vs_bsrd.sh <query.fa> <bsrd_db_prefix> <out.tsv>
QUERY=$1
DB=$2
OUT=$3
blastn \
  -query "$QUERY" \
  -db "$DB" \
  -out "$OUT" \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" \
  -evalue 1e-5 \
  -max_target_seqs 10
echo "Wrote: $OUT"
