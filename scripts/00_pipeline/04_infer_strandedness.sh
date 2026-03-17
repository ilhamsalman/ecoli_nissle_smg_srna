#!/usr/bin/env bash
set -euo pipefail

# Infer library strandedness using RSeQC infer_experiment.py
# Input: one representative BAM + annotation (GTF/BED)
# Output: qc/strandedness/infer_experiment.<bam>.txt

BAM="alignment/bam/SRR11357510.sorted.bam"
ANNOT="ref/EcN_ASM71459v1.gtf"
OUTDIR="qc/strandedness"

mkdir -p "$OUTDIR"

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: Missing required tool '$1'"; exit 127; }; }
need infer_experiment.py

if [[ ! -f "$BAM" ]]; then
  echo "ERROR: BAM not found: $BAM"
  echo "Tip: check available BAMs: ls -lh alignment/bam | head"
  exit 1
fi

if [[ ! -f "$ANNOT" ]]; then
  echo "ERROR: Annotation not found: $ANNOT"
  echo "Tip: check ref/: ls -lh ref | head"
  exit 1
fi

out="$OUTDIR/infer_experiment.$(basename "$BAM").txt"

echo "==> Running infer_experiment.py"
echo "    BAM:   $BAM"
echo "    ANNOT: $ANNOT"
infer_experiment.py -i "$BAM" -r "$ANNOT" | tee "$out"

echo "Saved: $out"
