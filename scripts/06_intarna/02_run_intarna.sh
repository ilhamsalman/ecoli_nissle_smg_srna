#!/bin/bash
#
# Run IntaRNA for selected intergenic sRNAs against phase-matched DEG mRNAs
#
# Current layout:
#   results/sRNA_detection/06_intarna/sequences/srnas/exp/*.fasta
#   results/sRNA_detection/06_intarna/sequences/srnas/sta/*.fasta
#   results/sRNA_detection/06_intarna/sequences/mrnas/deg_exp/*.fasta
#   results/sRNA_detection/06_intarna/sequences/mrnas/deg_sta/*.fasta
#
# Selected ID files:
#   results/sRNA_detection/06_intarna/selected_exp_ids.txt
#   results/sRNA_detection/06_intarna/selected_sta_ids.txt
#
# Usage:
#   bash scripts/06_intarna/02_run_intarna.sh
#   bash scripts/06_intarna/02_run_intarna.sh --threads 4
#

set -euo pipefail

# ============================================================================
# Arguments
# ============================================================================

THREADS=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)
            if [[ $# -lt 2 ]]; then
                echo "ERROR: --threads requires a value"
                exit 1
            fi
            THREADS="$2"
            shift 2
            ;;
        --threads=*)
            THREADS="${1#*=}"
            shift
            ;;
        *)
            echo "ERROR: Unknown argument: $1"
            exit 1
            ;;
    esac
done

if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -lt 1 ]]; then
    echo "ERROR: threads must be a positive integer"
    exit 1
fi

# ============================================================================
# Paths
# ============================================================================

BASE_DIR="results/sRNA_detection/06_intarna"

SRNA_DIR="${BASE_DIR}/sequences/srnas"
MRNA_DIR="${BASE_DIR}/sequences/mrnas"

SRNA_EXP_DIR="${SRNA_DIR}/exp"
SRNA_STA_DIR="${SRNA_DIR}/sta"

MRNA_EXP_DIR="${MRNA_DIR}/deg_exp"
MRNA_STA_DIR="${MRNA_DIR}/deg_sta"

EXP_IDS="${BASE_DIR}/selected_exp_ids.txt"
STA_IDS="${BASE_DIR}/selected_sta_ids.txt"

OUT_EXP="${BASE_DIR}/intarna_exp"
OUT_STA="${BASE_DIR}/intarna_sta"

mkdir -p "${OUT_EXP}" "${OUT_STA}"

# ============================================================================
# Header
# ============================================================================

echo "================================================================================"
echo "IntaRNA sRNA-mRNA Interaction Predictions"
echo "================================================================================"
echo "Started : $(date)"
echo "Threads : ${THREADS}"
echo ""

# ============================================================================
# Checks
# ============================================================================

if ! command -v IntaRNA >/dev/null 2>&1; then
    echo "ERROR: IntaRNA not found."
    echo "Install with: conda install -c bioconda intarna"
    exit 1
fi

echo "IntaRNA : $(IntaRNA --version 2>&1 | head -1)"
echo ""

for d in "${SRNA_EXP_DIR}" "${SRNA_STA_DIR}" "${MRNA_EXP_DIR}" "${MRNA_STA_DIR}"; do
    if [[ ! -d "$d" ]]; then
        echo "ERROR: Missing directory: $d"
        exit 1
    fi
done

for f in "${EXP_IDS}" "${STA_IDS}"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Missing selected ID file: $f"
        exit 1
    fi
done

# ============================================================================
# Helpers
# ============================================================================

JOB_DONE=0
JOB_TOTAL=0

count_fasta() {
    local dir="$1"
    find "$dir" -maxdepth 1 -type f -name "*.fasta" | wc -l
}

count_ids() {
    local file="$1"
    grep -c '[^[:space:]]' "$file" 2>/dev/null || echo 0
}

run_intarna() {
    local srna_fasta="$1"
    local mrna_fasta="$2"
    local output_csv="$3"

    local srna_name
    local mrna_name
    local energy

    srna_name=$(basename "$srna_fasta" .fasta)
    mrna_name=$(basename "$mrna_fasta" .fasta)

    if [[ ! -f "$srna_fasta" ]]; then
        echo "  WARNING: Missing sRNA FASTA: $srna_fasta"
        return
    fi

    if [[ ! -f "$mrna_fasta" ]]; then
        echo "  WARNING: Missing mRNA FASTA: $mrna_fasta"
        return
    fi

    IntaRNA \
        -q "$srna_fasta" \
        -t "$mrna_fasta" \
        --outMode=C \
        --out="$output_csv" \
        --seedBP=7 \
        --threads="${THREADS}" \
        2>/dev/null

    JOB_DONE=$((JOB_DONE + 1))

    if [[ -s "$output_csv" ]]; then
        energy=$(tail -n 1 "$output_csv" | awk -F';' '{print $NF}')
        printf "  [%4d/%d]  %-22s -> %-20s  dG=%s\n" \
            "$JOB_DONE" "$JOB_TOTAL" "$srna_name" "$mrna_name" "$energy"
    else
        printf "  [%4d/%d]  %-22s -> %-20s  (no interaction)\n" \
            "$JOB_DONE" "$JOB_TOTAL" "$srna_name" "$mrna_name"
    fi
}

run_batch() {
    local id_file="$1"
    local srna_dir="$2"
    local mrna_dir="$3"
    local out_dir="$4"
    local label="$5"

    local n_targets
    local srna_id
    local srna_fasta
    local mrna_fasta

    n_targets=$(count_fasta "$mrna_dir")

    echo "================================================================================"
    echo "${label}"
    echo "================================================================================"
    echo "sRNA dir  : ${srna_dir}"
    echo "mRNA dir  : ${mrna_dir}"
    echo "Output dir: ${out_dir}"
    echo "Targets   : ${n_targets}"
    echo ""

    if [[ "$n_targets" -eq 0 ]]; then
        echo "ERROR: No mRNA FASTA files found in ${mrna_dir}"
        exit 1
    fi

    while IFS= read -r srna_id || [[ -n "$srna_id" ]]; do
        [[ -z "$srna_id" ]] && continue

        srna_id=$(echo "$srna_id" | tr -d '\r')
        srna_fasta="${srna_dir}/${srna_id}.fasta"

        if [[ ! -f "$srna_fasta" ]]; then
            echo "  WARNING: Selected ID not found as FASTA: ${srna_fasta}"
            continue
        fi

        echo "  -- ${srna_id} vs ${n_targets} targets --"

        for mrna_fasta in "${mrna_dir}"/*.fasta; do
            [[ -f "$mrna_fasta" ]] || continue

            run_intarna \
                "$srna_fasta" \
                "$mrna_fasta" \
                "${out_dir}/${srna_id}_vs_$(basename "$mrna_fasta" .fasta).csv"
        done

        echo ""
    done < "$id_file"
}

# ============================================================================
# Job counts
# ============================================================================

N_EXP_IDS=$(count_ids "${EXP_IDS}")
N_STA_IDS=$(count_ids "${STA_IDS}")

N_EXP_MRNA=$(count_fasta "${MRNA_EXP_DIR}")
N_STA_MRNA=$(count_fasta "${MRNA_STA_DIR}")

JOB_TOTAL=$(( (N_EXP_IDS * N_EXP_MRNA) + (N_STA_IDS * N_STA_MRNA) ))

echo "Jobs queued:"
echo "  EXP : ${N_EXP_IDS} sRNAs x ${N_EXP_MRNA} mRNAs = $((N_EXP_IDS * N_EXP_MRNA))"
echo "  STA : ${N_STA_IDS} sRNAs x ${N_STA_MRNA} mRNAs = $((N_STA_IDS * N_STA_MRNA))"
echo "  ---------------------------------------------"
echo "  TOTAL = ${JOB_TOTAL}"
echo ""

if [[ "$JOB_TOTAL" -eq 0 ]]; then
    echo "ERROR: Total jobs = 0"
    echo "Check selected ID files and FASTA directories."
    exit 1
fi

# ============================================================================
# Run phase-matched predictions
# ============================================================================

run_batch \
    "${EXP_IDS}" \
    "${SRNA_EXP_DIR}" \
    "${MRNA_EXP_DIR}" \
    "${OUT_EXP}" \
    "EXPONENTIAL PHASE — selected sRNAs vs DEG exp mRNAs"

run_batch \
    "${STA_IDS}" \
    "${SRNA_STA_DIR}" \
    "${MRNA_STA_DIR}" \
    "${OUT_STA}" \
    "STATIONARY PHASE — selected sRNAs vs DEG sta mRNAs"

# ============================================================================
# Footer
# ============================================================================

echo "================================================================================"
echo "All IntaRNA jobs complete — $(date)"
echo "================================================================================"
echo ""
echo "Next step:"
echo "  bash scripts/06_intarna/03_combine_results.sh"
echo ""