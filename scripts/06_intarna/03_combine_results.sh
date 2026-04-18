#!/bin/bash

set -euo pipefail

OUT_DIR="results/sRNA_detection/06_intarna"
EXP_DIR="${OUT_DIR}/intarna_exp"
STA_DIR="${OUT_DIR}/intarna_sta"

echo "================================================================================"
echo "Combining IntaRNA Results (phase-specific)"
echo "================================================================================"
echo ""

if [[ ! -d "$EXP_DIR" ]]; then
    echo "ERROR: Missing EXP directory: $EXP_DIR"
    exit 1
fi

if [[ ! -d "$STA_DIR" ]]; then
    echo "ERROR: Missing STA directory: $STA_DIR"
    exit 1
fi

combine_one_phase() {
    local phase="$1"
    local src_dir="$2"
    local combined_csv="$3"
    local combined_tsv="$4"
    local ranked_tsv="$5"

    local n_files=0
    local n_predictions=0

    echo "--------------------------------------------------------------------------------"
    echo "${phase}"
    echo "--------------------------------------------------------------------------------"

    echo "phase;srna;target;id1;start1;end1;id2;start2;end2;subseqDP;hybridDP;E" > "$combined_csv"

    for csv_file in "$src_dir"/*.csv; do
        [[ -f "$csv_file" ]] || continue

        local fname
        local srna
        local target

        fname=$(basename "$csv_file" .csv)
        srna=$(echo "$fname" | awk -F'_vs_' '{print $1}')
        target=$(echo "$fname" | awk -F'_vs_' '{print $2}')

        [[ -s "$csv_file" ]] || continue

        tail -n +2 "$csv_file" | \
            awk -v p="$phase" -v s="$srna" -v t="$target" \
            'BEGIN{FS=OFS=";"} NF>0 {print p,s,t,$0}' \
            >> "$combined_csv"

        n_files=$((n_files + 1))
    done

    sed 's/;/\t/g' "$combined_csv" > "$combined_tsv"

    n_predictions=$(tail -n +2 "$combined_csv" | grep -c '[^[:space:]]' || true)

    if [[ "$n_predictions" -gt 0 ]]; then
        {
            head -n 1 "$combined_tsv"
            tail -n +2 "$combined_tsv" | sort -k12,12n
        } > "$ranked_tsv"
    else
        cp "$combined_tsv" "$ranked_tsv"
    fi

    echo "  CSV files merged : ${n_files}"
    echo "  Predictions      : ${n_predictions}"
    echo "  Output files:"
    echo "    ${combined_csv}"
    echo "    ${combined_tsv}"
    echo "    ${ranked_tsv}"
    echo ""

    if [[ "$n_predictions" -gt 0 ]]; then
        echo "  Top 10 interactions by ΔG"
        printf "  %-24s %-24s %15s\n" "sRNA" "Target" "ΔG (kcal/mol)"
        echo "  ---------------------------------------------------------------------"

        awk -F'\t' 'NR>1 && NR<=11 {printf "  %-24s %-24s %15s\n", $2, $3, $12}' "$ranked_tsv"
        echo ""
    fi
}

combine_one_phase \
    "EXP" \
    "$EXP_DIR" \
    "${OUT_DIR}/exp_predictions_combined.csv" \
    "${OUT_DIR}/exp_predictions_combined.tsv" \
    "${OUT_DIR}/exp_predictions_ranked.tsv"

combine_one_phase \
    "STA" \
    "$STA_DIR" \
    "${OUT_DIR}/sta_predictions_combined.csv" \
    "${OUT_DIR}/sta_predictions_combined.tsv" \
    "${OUT_DIR}/sta_predictions_ranked.tsv"

echo "phase;srna;target;id1;start1;end1;id2;start2;end2;subseqDP;hybridDP;E" \
    > "${OUT_DIR}/all_predictions_combined.csv"

tail -n +2 "${OUT_DIR}/exp_predictions_combined.csv" >> "${OUT_DIR}/all_predictions_combined.csv"
tail -n +2 "${OUT_DIR}/sta_predictions_combined.csv" >> "${OUT_DIR}/all_predictions_combined.csv"

sed 's/;/\t/g' "${OUT_DIR}/all_predictions_combined.csv" > "${OUT_DIR}/all_predictions_combined.tsv"

echo "Optional merged file:"
echo "  ${OUT_DIR}/all_predictions_combined.csv"
echo "  ${OUT_DIR}/all_predictions_combined.tsv"
echo ""
echo "================================================================================"
echo "Done"
echo "================================================================================"