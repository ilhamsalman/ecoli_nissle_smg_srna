#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# STEP 3 — Antisense overlaps between significant antisense sRNAs and DEGs
# Uses: bedtools intersect -S (opposite strands)
# ============================================================

# --- Adjust paths ---
SRNA_EXP="sig_antisense_srna_exp.bed"
SRNA_STA="sig_antisense_srna_sta.bed"

DEG_EXP="inputs/DEG_exp_with_fc.bed"
DEG_STA="inputs/DEG_sta_with_fc.bed"

OUT_EXP="antisense_deg_exp_overlap.bed"
OUT_STA="antisense_deg_sta_overlap.bed"

bedtools intersect -a "$SRNA_EXP" -b "$DEG_EXP" -S -wa -wb > "$OUT_EXP"
bedtools intersect -a "$SRNA_STA" -b "$DEG_STA" -S -wa -wb > "$OUT_STA"

echo "[OK] overlaps:"
echo "  $OUT_EXP  (lines: $(wc -l < "$OUT_EXP"))"
echo "  $OUT_STA  (lines: $(wc -l < "$OUT_STA"))"

awk -F'\t' '($6==$12){c++} END{print "[CHECK] same-strand overlaps EXP:", c+0}' "$OUT_EXP"
awk -F'\t' '($6==$12){c++} END{print "[CHECK] same-strand overlaps STA:", c+0}' "$OUT_STA"
