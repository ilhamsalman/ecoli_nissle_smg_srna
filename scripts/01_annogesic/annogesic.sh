#!/usr/bin/env bash
set -euo pipefail

# ANNOgesic workflow (exponential phase)
# NOTE: Run the same commands for stationary phase using annogesic_sta.
#   1) change PJ="annogesic_exp" -> PJ="annogesic_sta"
#   2) change all MG_exp_* -> MG_sta_* and NG_exp_* -> NG_sta_*
# Then rerun the script

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: Missing required tool '$1'"; exit 127; }; }
need annogesic

PJ="annogesic_exp"
GFF="$PJ/input/references/annotations/annotation.gff"
FASTA="$PJ/input/references/fasta_files/genome.fasta"
WIGDIR="$PJ/input/wigs/fragment"

# sanity checks (fail early, reproducible behavior)
[[ -f "$GFF" ]]   || { echo "ERROR: Missing $GFF"; exit 1; }
[[ -f "$FASTA" ]] || { echo "ERROR: Missing $FASTA"; exit 1; }
[[ -d "$WIGDIR" ]] || { echo "ERROR: Missing $WIGDIR"; exit 1; }

# Transcript detection
annogesic transcript \
  -pj "$PJ" \
  -g "$GFF" \
  -fl \
    "$WIGDIR/MG_exp_1_plus.wig:frag:1:a:+" \
    "$WIGDIR/MG_exp_1_minus.wig:frag:1:a:-" \
    "$WIGDIR/MG_exp_2_plus.wig:frag:1:b:+" \
    "$WIGDIR/MG_exp_2_minus.wig:frag:1:b:-" \
    "$WIGDIR/NG_exp_1_plus.wig:frag:2:a:+" \
    "$WIGDIR/NG_exp_1_minus.wig:frag:2:a:-" \
    "$WIGDIR/NG_exp_2_plus.wig:frag:2:b:+" \
    "$WIGDIR/NG_exp_2_minus.wig:frag:2:b:-" \
  -rf 1_2 2_2

# sRNA prediction
annogesic srna \
  -pj "$PJ" \
  -a "$PJ/output/transcripts/gffs/annotation_transcript.gff" \
  -g "$GFF" \
  -f "$FASTA" \
  -d none \
  -fl \
    "$WIGDIR/MG_exp_1_plus.wig:frag:1:a:+" \
    "$WIGDIR/MG_exp_1_minus.wig:frag:1:a:-" \
    "$WIGDIR/MG_exp_2_plus.wig:frag:1:b:+" \
    "$WIGDIR/MG_exp_2_minus.wig:frag:1:b:-" \
    "$WIGDIR/NG_exp_1_plus.wig:frag:2:a:+" \
    "$WIGDIR/NG_exp_1_minus.wig:frag:2:a:-" \
    "$WIGDIR/NG_exp_2_plus.wig:frag:2:b:+" \
    "$WIGDIR/NG_exp_2_minus.wig:frag:2:b:-" \
  -rf 1_2 2_2

RFAM_CM="$PJ/input/databases/Rfam.cm"
RIBO_IDS="$PJ/input/riboswitch_ID_file/Rfam_riboswitch_ID.txt"
THERMO_IDS="$PJ/input/RNA_thermometer_ID_file/Rfam_RNA_thermometer_ID.txt"

[[ -f "$RFAM_CM" ]]   || { echo "ERROR: Missing $RFAM_CM"; exit 1; }
[[ -f "$RIBO_IDS" ]]  || { echo "ERROR: Missing $RIBO_IDS"; exit 1; }
[[ -f "$THERMO_IDS" ]]|| { echo "ERROR: Missing $THERMO_IDS"; exit 1; }

annogesic riboswitch_thermometer \
  -pj "$PJ" \
  --annotation_files "$GFF" \
  --transcript_files "$PJ/output/transcripts/gffs/annotation_transcript.gff" \
  --fasta_files "$FASTA" \
  --rfam_path "$RFAM_CM" \
  --riboswitch_id_file "$RIBO_IDS" \
  --rna_thermometer_id_file "$THERMO_IDS"
