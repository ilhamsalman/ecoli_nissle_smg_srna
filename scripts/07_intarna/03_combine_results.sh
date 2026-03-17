#!/bin/bash
#
# Combine all IntaRNA prediction CSVs into one ranked TSV.
# Run from project root after 02_run_intarna.sh completes.
#
# Usage:
#   bash scripts/07_intarna/03_combine_results.sh
#

set -e

OUT_DIR="results/sRNA_detection/07_intarna"
combined_csv="${OUT_DIR}/all_predictions_combined.csv"
combined_tsv="${OUT_DIR}/all_predictions_combined.tsv"

echo "================================================================================"
echo "Combining IntaRNA Results"
echo "================================================================================"
echo ""

# Write header
echo "priority;srna;target;id1;start1;end1;id2;start2;end2;subseqDP;hybridDP;E" \
    > "$combined_csv"

n_files=0

for subdir in \
    "${OUT_DIR}/priority1_leuC_antisense" \
    "${OUT_DIR}/priority2_intergenic_exp" \
    "${OUT_DIR}/priority3_intergenic_sta"; do

    [ -d "$subdir" ] || continue

    dir_name=$(basename "$subdir")
    case "$dir_name" in
        priority1_leuC_antisense) label="P1_leuC_antisense" ;;
        priority2_intergenic_exp) label="P2_intergenic_exp" ;;
        priority3_intergenic_sta) label="P3_intergenic_sta" ;;
        *)                        label="$dir_name" ;;
    esac

    n_dir=0
    for csv_file in "$subdir"/*.csv; do
        [ -f "$csv_file" ] || continue
        fname=$(basename "$csv_file" .csv)
        srna=$(echo   "$fname" | awk -F'_vs_' '{print $1}')
        target=$(echo "$fname" | awk -F'_vs_' '{print $2}')
        tail -n +2 "$csv_file" | \
            awk -v p="$label" -v s="$srna" -v t="$target" \
                'BEGIN{FS=OFS=";"} NF>0 {print p,s,t,$0}' \
            >> "$combined_csv"
        n_dir=$(( n_dir + 1 ))
    done

    echo "  ${label}: ${n_dir} CSV files"
    n_files=$(( n_files + n_dir ))
done

echo ""
echo "Total CSV files merged: ${n_files}"

# Tab-delimited version
sed 's/;/\t/g' "$combined_csv" > "$combined_tsv"

n_predictions=$(tail -n +2 "$combined_csv" | grep -c '[^[:space:]]' || true)
echo "Total predictions:      ${n_predictions}"
echo ""

# ============================================================================
# Summary — top 20 by ΔG
# ============================================================================

echo "================================================================================"
echo "Top 20 interactions by ΔG"
echo "================================================================================"
echo ""

printf "%-30s %-24s %-20s %15s\n" "sRNA" "Priority" "Target" "ΔG (kcal/mol)"
echo "--------------------------------------------------------------------------"

tail -n +2 "$combined_tsv" | \
    awk 'NF>0 {print $2, $1, $3, $NF}' | \
    sort -k4,4n | \
    head -20 | \
    while read -r srna priority target energy; do
        printf "%-30s %-24s %-20s %15s\n" "$srna" "$priority" "$target" "$energy"
    done

echo ""
echo "  (top 20 of ${n_predictions} — full list: all_predictions_combined.tsv)"
echo ""
echo "  ΔG < -15      Strong   → high priority for validation"
echo "  ΔG -10 to -15 Moderate"
echo "  ΔG > -10      Weak"
echo ""
echo "================================================================================"
echo "✓ Done: ${combined_tsv}"
echo "================================================================================"
