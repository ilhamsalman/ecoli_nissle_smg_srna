#!/bin/bash
set -euo pipefail

INPUT="meta/samples_nissle_GSE147272.tsv"
OUTDIR="fastq"

mkdir -p "$OUTDIR"

tail -n +2 "$INPUT" | while IFS=$'\t' read -r sample SRR phase gravity; do
    echo "=== Downloading $sample ($SRR) ==="
    prefetch "$SRR"
    fasterq-dump "$SRR" --split-files --threads 8 --outdir "$OUTDIR"
done
