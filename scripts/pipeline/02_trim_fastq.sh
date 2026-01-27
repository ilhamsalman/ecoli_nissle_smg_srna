#!/usr/bin/env bash
set -euo pipefail

THREADS=8
IN_DIR="fastq"
OUT_DIR="trimmed"
QC_DIR="qc/fastp"

mkdir -p "$OUT_DIR" "$QC_DIR"

shopt -s nullglob

# Find R1 files
r1_files=("$IN_DIR"/*_1.fastq)
if [ ${#r1_files[@]} -eq 0 ]; then
  echo "ERROR: No files matching $IN_DIR/*_1.fastq"
  echo "Check your filenames with: ls -lh $IN_DIR | head"
  exit 1
fi

for r1 in "${r1_files[@]}"; do
  base=$(basename "$r1" _1.fastq)
  r2="$IN_DIR/${base}_2.fastq"

  if [ ! -f "$r2" ]; then
    echo "WARNING: missing pair for $r1"
    echo "         expected: $r2"
    continue
  fi

  echo "=== Trimming ${base} ==="

  fastp \
    -i "$r1" -I "$r2" \
    -o "$OUT_DIR/${base}_1.trimmed.fastq.gz" \
    -O "$OUT_DIR/${base}_2.trimmed.fastq.gz" \
    -w "$THREADS" \
    --detect_adapter_for_pe \
    --html "$QC_DIR/${base}.html" \
    --json "$QC_DIR/${base}.json"
done

echo "Done. Outputs:"
echo "  $OUT_DIR/*.trimmed.fastq.gz"
echo "  $QC_DIR/*.html  $QC_DIR/*.json"
