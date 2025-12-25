#!/usr/bin/env bash
set -euo pipefail

INPUT="fastq"
OUTPUT="trimmed"
THREADS=8

mkdir -p "$OUTPUT"

for r1 in ${INPUT}/*_1.fastq.gz; do
    base=$(basename "$r1" _1.fastq.gz)
    r2="${INPUT}/${base}_2.fastq.gz"

    echo "Trimming $base"

    fastp \
      -i "$r1" \
      -I "$r2" \
      -o "${OUTPUT}/${base}_1.trimmed.fastq.gz" \
      -O "${OUTPUT}/${base}_2.trimmed.fastq.gz" \
      --detect_adapter_for_pe \
      --cut_front \
      --cut_tail \
      --cut_window_size 4 \
      --cut_mean_quality 20 \
      --length_required 18 \
      --thread "$THREADS" \
      --html "${OUTPUT}/${base}_fastp.html" \
      --json "${OUTPUT}/${base}_fastp.json"
done
