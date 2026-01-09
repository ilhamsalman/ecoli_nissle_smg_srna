#!/usr/bin/env bash
set -euo pipefail

# ANNOgesic exponential-phase analysis
# E. coli Nissle 1917 (GSE147272)
# Strand-specific RNA-seq (reverse, RF)

conda activate annogesic_env

# Transcript detection
annogesic transcript \
  -pj annogesic_exp \
  -g annogesic_exp/input/references/annotations/annotation.gff \
  -fl \
    MG_exp_1_plus.wig:frag:1:a:+ \
    MG_exp_1_minus.wig:frag:1:a:- \
    MG_exp_2_plus.wig:frag:1:b:+ \
    MG_exp_2_minus.wig:frag:1:b:- \
    NG_exp_1_plus.wig:frag:2:a:+ \
    NG_exp_1_minus.wig:frag:2:a:- \
    NG_exp_2_plus.wig:frag:2:b:+ \
    NG_exp_2_minus.wig:frag:2:b:- \
  -rf 1_2 2_2

# sRNA prediction
annogesic srna \
  -pj annogesic_exp \
  -a annogesic_exp/output/transcripts/gffs/annotation_transcript.gff \
  -g annogesic_exp/input/references/annotations/annotation.gff \
  -f annogesic_exp/input/references/fasta_files/genome.fasta \
  -d none \
  -fl \
    annogesic_exp/input/wigs/fragment/MG_exp_1_plus.wig:frag:1:a:+ \
    annogesic_exp/input/wigs/fragment/MG_exp_1_minus.wig:frag:1:a:- \
    annogesic_exp/input/wigs/fragment/MG_exp_2_plus.wig:frag:1:b:+ \
    annogesic_exp/input/wigs/fragment/MG_exp_2_minus.wig:frag:1:b:- \
    annogesic_exp/input/wigs/fragment/NG_exp_1_plus.wig:frag:2:a:+ \
    annogesic_exp/input/wigs/fragment/NG_exp_1_minus.wig:frag:2:a:- \
    annogesic_exp/input/wigs/fragment/NG_exp_2_plus.wig:frag:2:b:+ \
    annogesic_exp/input/wigs/fragment/NG_exp_2_minus.wig:frag:2:b:- \
  -rf 1_2 2_2

# Terminator prediction using sRNA GFF
annogesic terminator \
  -pj annogesic_exp \
  -a annogesic_exp/output/sRNAs/gffs/best_candidates/annotation_sRNA.gff \
  -g annogesic_exp/input/references/annotations/annotation.gff \
  -f annogesic_exp/input/references/fasta_files/genome.fasta \
  -fl \
    annogesic_exp/input/wigs/fragment/MG_exp_1_plus.wig:frag:1:a:+ \
    annogesic_exp/input/wigs/fragment/MG_exp_1_minus.wig:frag:1:a:- \
    annogesic_exp/input/wigs/fragment/MG_exp_2_plus.wig:frag:1:b:+ \
    annogesic_exp/input/wigs/fragment/MG_exp_2_minus.wig:frag:1:b:- \
    annogesic_exp/input/wigs/fragment/NG_exp_1_plus.wig:frag:2:a:+ \
    annogesic_exp/input/wigs/fragment/NG_exp_1_minus.wig:frag:2:a:- \
    annogesic_exp/input/wigs/fragment/NG_exp_2_plus.wig:frag:2:b:+ \
    annogesic_exp/input/wigs/fragment/NG_exp_2_minus.wig:frag:2:b:- \
  -rf 1_2 2_2

# Riboswitch / RNA thermometer prediction
annogesic riboswitch_thermometer \
  -pj annogesic_exp \
  -p thermometer \
  --annotation_files annogesic_exp/input/references/annotations/annotation.gff \
  --transcript_files annogesic_exp/output/transcripts/gffs/annotation_transcript.gff \
  --fasta_files annogesic_exp/input/references/fasta_files/genome.fasta \
  --rfam_path annogesic_exp/input/databases/Rfam.cm \
  --rna_thermometer_id_file annogesic_exp/input/RNA_thermometer_ID_file/Rfam_RNA_thermometer_ID.txt
