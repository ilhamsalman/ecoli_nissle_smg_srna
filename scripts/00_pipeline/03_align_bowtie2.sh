#!/usr/bin/env bash
set -euo pipefail

THREADS=8
IDX="ref/EcN_CP007799.1"          # bowtie2 index prefix
TRIMDIR="trimmed"
OUTDIR="alignment/bam"
LOGDIR="alignment/logs"

mkdir -p "$OUTDIR" "$LOGDIR"

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: Missing required tool '$1'"; exit 127; }; }

need bowtie2
need samtools

# Check index exists
if ! ls "${IDX}".*.bt2 >/dev/null 2>&1 && ! ls "${IDX}".*.bt2l >/dev/null 2>&1; then
  echo "ERROR: Bowtie2 index not found for prefix: $IDX"
  echo "Expected files like: ${IDX}.1.bt2 (or .bt2l)"
  echo "If you only have a FASTA, build the index with:"
  echo "  bowtie2-build ref/EcN_CP007799.1.fasta $IDX"
  exit 1
fi

shopt -s nullglob
r1s=("$TRIMDIR"/*_1.trimmed.fastq.gz)
if [[ ${#r1s[@]} -eq 0 ]]; then
  echo "ERROR: No trimmed R1 files found in $TRIMDIR"
  echo "Run trimming first: scripts/1.pipeline/02_trim_fastq.sh"
  exit 1
fi

for r1 in "${r1s[@]}"; do
  base="$(basename "$r1" _1.trimmed.fastq.gz)"
  r2="$TRIMDIR/${base}_2.trimmed.fastq.gz"
  [[ -f "$r2" ]] || { echo "ERROR: missing mate for $r1 (expected $r2)"; exit 1; }

  outbam="$OUTDIR/${base}.sorted.bam"
  log="$LOGDIR/${base}.bowtie2.log"

  echo "=== Aligning ${base} ==="
  bowtie2 -x "$IDX" --very-sensitive -p "$THREADS" -1 "$r1" -2 "$r2" 2> "$log" \
  | samtools view -bS - \
  | samtools sort -@ "$THREADS" -o "$outbam" -

  samtools index "$outbam"
done

echo "Done. BAMs: $OUTDIR  Logs: $LOGDIR"
