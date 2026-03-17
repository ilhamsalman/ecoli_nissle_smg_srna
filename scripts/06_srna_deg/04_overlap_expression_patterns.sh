#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# STEP 4 — Summarize overlaps + expression direction (OPPOSITE/SAME/NEUTRAL)
# Input: bedtools output from Step 3
# Output: tidy TSV per condition
# ============================================================

EXP_IN="antisense_deg_exp_overlap.bed"
STA_IN="antisense_deg_sta_overlap.bed"

EXP_OUT="antisense_deg_exp_analysis.tsv"
STA_OUT="antisense_deg_sta_analysis.tsv"

analyze_one () {
  local IN="$1"
  local OUT="$2"

  awk 'BEGIN{
      FS="\t"; OFS="\t";
      print "antisense_sRNA","asRNA_logFC","asRNA_strand",
            "DEG_locus","DEG_gene","DEG_logFC","DEG_strand",
            "overlap_bp","expression_pattern"
  }
  {
      asrna=$4; asfc=$5; asstr=$6;
      as_start=$2; as_end=$3;

      deg_start=$8; deg_end=$9;
      deg_name=$10; degfc=$11; degstr=$12;

      ostart = (as_start > deg_start) ? as_start : deg_start;
      oend   = (as_end   < deg_end)   ? as_end   : deg_end;
      olen   = oend - ostart;
      if (olen < 0) olen = 0;

      split(deg_name, parts, /\|/);
      locus = parts[1];
      gene  = (length(parts) > 1 && parts[2] != "") ? parts[2] : "NA";

      if (asfc+0 == 0 || degfc+0 == 0) {
        pat="NEUTRAL";
      } else if ((asfc+0 > 0 && degfc+0 < 0) || (asfc+0 < 0 && degfc+0 > 0)) {
        pat="OPPOSITE";
      } else {
        pat="SAME";
      }

      print asrna, asfc, asstr, locus, gene, degfc, degstr, olen, pat;
  }' "$IN" > "$OUT"

  echo "[OK] wrote $OUT (lines: $(wc -l < "$OUT"))"
}

analyze_one "$EXP_IN" "$EXP_OUT"
analyze_one "$STA_IN" "$STA_OUT"
