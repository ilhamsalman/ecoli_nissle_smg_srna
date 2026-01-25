#!/usr/bin/env bash
set -euo pipefail

# Step (ANNOgesic prep): Generate strand-specific WIG coverage tracks from BAM files
#
# Purpose:
#   ANNOgesic uses strand-specific coverage tracks. This script converts BAM alignments
#   into WIG files for each sample.
#
# Strandedness:
#   This script assumes an RF (reverse-forward) stranded RNA-seq library.
#   Therefore:
#     RNA plus-strand coverage  <- alignments on the minus genomic strand  (-strand -)
#     RNA minus-strand coverage <- alignments on the plus genomic strand   (-strand +)
#
# Input:
#   Sorted BAM files in: annogesic_exp/input/BAMs/BAMs_map_reference_genomes/fragment/*.bam
# Output:
#   WIG files in:        annogesic_exp/input/wigs/fragment/
#     <sample>.plus.wig  (RNA + strand coverage)
#     <sample>.minus.wig (RNA - strand coverage)
#


shopt -s nullglob
bams=("$BAMDIR"/*.bam)
shopt -u nullglob

if [[ ${#bams[@]} -eq 0 ]]; then
  echo "ERROR: No BAM files found in: $BAMDIR"
  exit 1
fi

for bam in "${bams[@]}"; do
  base="$(basename "$bam" .bam)"
  echo "==> Processing $base"

  # RNA plus strand (RF) corresponds to alignments on genomic minus strand
  bedtools genomecov -bg --split -strand - -ibam "$bam" \
    | sort -k1,1 -k2,2n \
    | awk 'BEGIN{OFS="\t"} { if ($1!=c) {c=$1; print "variableStep chrom="c} print $2+1, $4 }' \
    > "$OUTDIR/${base}.plus.wig"

  # RNA minus strand (RF) corresponds to alignments on genomic plus strand
  bedtools genomecov -bg --split -strand + -ibam "$bam" \
    | sort -k1,1 -k2,2n \
    | awk 'BEGIN{OFS="\t"} { if ($1!=c) {c=$1; print "variableStep chrom="c} print $2+1, $4 }' \
    > "$OUTDIR/${base}.minus.wig"
done

echo "WIG generation complete: $OUTDIR"

