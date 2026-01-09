#!/usr/bin/env bash
set -euo pipefail

BAMDIR="annogesic_exp/input/BAMs/BAMs_map_reference_genomes/fragment"
OUTDIR="annogesic_exp/input/wigs/fragment"

mkdir -p "$OUTDIR"

for bam in "$BAMDIR"/*.bam; do
  base=$(basename "$bam" .bam)
  echo "== Processing $base =="

  # RNA +  (alignment -)
  bedtools genomecov -bg -split -strand - -ibam "$bam" \
    | sort -k1,1 -k2,2n \
    | awk 'BEGIN{OFS="\t"} {if ($1!=c){c=$1; print "variableStep chrom="c} print $2+1, $4}' \
    > "$OUTDIR/${base}.plus.wig"

  # RNA -  (alignment +)
  bedtools genomecov -bg -split -strand + -ibam "$bam" \
    | sort -k1,1 -k2,2n \
    | awk 'BEGIN{OFS="\t"} {if ($1!=c){c=$1; print "variableStep chrom="c} print $2+1, $4}' \
    > "$OUTDIR/${base}.minus.wig"
done

echo "WIG generation complete."
