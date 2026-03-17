#!/usr/bin/env bash
set -euo pipefail

INPUT="meta/samples_nissle_GSE147272.tsv"
OUTDIR="fastq"
THREADS=8

mkdir -p "$OUTDIR"

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: Missing required tool '$1'"; exit 127; }; }
need prefetch
need fasterq-dump

tail -n +2 "$INPUT" | while IFS=$'\t' read -r sample SRR phase gravity; do
  echo "=== Downloading $sample ($SRR) ==="

  prefetch "$SRR"
  fasterq-dump "$SRR" --split-files --threads "$THREADS" --outdir "$OUTDIR"

  # Rename SRR outputs to sample-based names (keeps pipeline consistent)
  if [[ -f "$OUTDIR/${SRR}_1.fastq" && -f "$OUTDIR/${SRR}_2.fastq" ]]; then
    mv -f "$OUTDIR/${SRR}_1.fastq" "$OUTDIR/${sample}_1.fastq"
    mv -f "$OUTDIR/${SRR}_2.fastq" "$OUTDIR/${sample}_2.fastq"
  else
    echo "ERROR: Expected $OUTDIR/${SRR}_1.fastq and ${SRR}_2.fastq not found"
    exit 1
  fi
done
