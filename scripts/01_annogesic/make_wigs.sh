#!/usr/bin/env bash
set -euo pipefail

# Generate RF strand-specific WIGs from fragment BAMs for ANNOgesic
# Usage:
#   bash scripts/01_annogesic/make_wigs.sh annogesic_exp
#   bash scripts/01_annogesic/make_wigs.sh annogesic_sta

PROJ="${1:-}"
if [[ -z "$PROJ" ]]; then
  echo "ERROR: provide project dir (annogesic_exp or annogesic_sta)"
  exit 1
fi

BAMDIR="$PROJ/input/BAMs/BAMs_map_reference_genomes/fragment"
OUTDIR="$PROJ/input/wigs/fragment"
mkdir -p "$OUTDIR"

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: Missing required tool '$1'"; exit 127; }; }
need bedtools

shopt -s nullglob
bams=("$BAMDIR"/*.bam)
shopt -u nullglob
[[ ${#bams[@]} -gt 0 ]] || { echo "ERROR: No BAMs in $BAMDIR"; exit 1; }

for bam in "${bams[@]}"; do
  base="$(basename "$bam" .bam)"
  echo "==> $PROJ : $base"

  # RF: RNA plus strand <- genomic minus
  bedtools genomecov -bg --split -strand - -ibam "$bam" \
    | sort -k1,1 -k2,2n \
    | awk 'BEGIN{OFS="\t"} { if ($1!=c) {c=$1; print "variableStep chrom="c} print $2+1, $4 }' \
    > "$OUTDIR/${base}_plus.wig"

  # RF: RNA minus strand <- genomic plus
  bedtools genomecov -bg --split -strand + -ibam "$bam" \
    | sort -k1,1 -k2,2n \
    | awk 'BEGIN{OFS="\t"} { if ($1!=c) {c=$1; print "variableStep chrom="c} print $2+1, $4 }' \
    > "$OUTDIR/${base}_minus.wig"
done

echo "Done: $OUTDIR"
