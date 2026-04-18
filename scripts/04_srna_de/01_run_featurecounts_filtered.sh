#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-8}"

EXP_SAF="results/sRNA_detection/03_shared_filtered/srna_exp.filtered.saf"
STA_SAF="results/sRNA_detection/03_shared_filtered/srna_sta.filtered.saf"

OUTDIR="results/sRNA_detection/04_srna_de_filtered"
mkdir -p "$OUTDIR"

EXP_BAMS=(
  annogesic_exp/input/BAMs/BAMs_map_reference_genomes/fragment/MG_exp_1.bam
  annogesic_exp/input/BAMs/BAMs_map_reference_genomes/fragment/MG_exp_2.bam
  annogesic_exp/input/BAMs/BAMs_map_reference_genomes/fragment/NG_exp_1.bam
  annogesic_exp/input/BAMs/BAMs_map_reference_genomes/fragment/NG_exp_2.bam
)

STA_BAMS=(
  annogesic_sta/input/BAMs/BAMs_map_reference_genomes/fragment/MG_sta_1.bam
  annogesic_sta/input/BAMs/BAMs_map_reference_genomes/fragment/MG_sta_2.bam
  annogesic_sta/input/BAMs/BAMs_map_reference_genomes/fragment/NG_sta_1.bam
  annogesic_sta/input/BAMs/BAMs_map_reference_genomes/fragment/NG_sta_2.bam
)

run_fc () {
  local SAF="$1"
  local OUT="$2"
  shift 2
  local BAMS=("$@")

  featureCounts -T "$THREADS" -F SAF -a "$SAF" -o "$OUT" \
    -s 2 -p -B -C \
    "${BAMS[@]}"
}

echo "== featureCounts EXP (-s 2) =="
run_fc "$EXP_SAF" "$OUTDIR/counts_exp_s2.txt" "${EXP_BAMS[@]}"

echo "== featureCounts STA (-s 2) =="
run_fc "$STA_SAF" "$OUTDIR/counts_sta_s2.txt" "${STA_BAMS[@]}"

echo "Done. Summaries:"
grep -H "^Assigned" "$OUTDIR"/counts_*_s2.txt.summary || true