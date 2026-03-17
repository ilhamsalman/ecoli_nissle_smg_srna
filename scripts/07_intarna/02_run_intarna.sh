#!/bin/bash
#
# Run IntaRNA sRNA-mRNA interaction predictions
#
# Priorities:
#   P1  Novel leuC antisense sRNAs    (exp_sRNA_00170, 00172 vs leuC)         2 jobs
#   P2  Novel intergenic sRNAs, exp   (exp_sRNA_02152 vs 89 DEG exp mRNAs)   89 jobs
#   P3  Novel intergenic sRNAs, sta   (6 sRNAs vs 209 DEG sta mRNAs)       1254 jobs
#
# Prerequisites:
#   conda activate intarna_env
#   python3 scripts/07_intarna/01_extract_sequences.py   <- must run first
#
# Usage:
#   bash scripts/07_intarna/02_run_intarna.sh
#   bash scripts/07_intarna/02_run_intarna.sh --threads 4
#

set -e

# ============================================================================
# Arguments
# ============================================================================

THREADS=1
for arg in "$@"; do
    case $arg in
        --threads=*) THREADS="${arg#*=}" ;;
        --threads)   shift; THREADS="$1" ;;
    esac
done

echo "================================================================================"
echo "IntaRNA sRNA-mRNA Interaction Predictions"
echo "================================================================================"
echo ""

if ! command -v IntaRNA &> /dev/null; then
    echo "ERROR: IntaRNA not found. Install: conda install -c bioconda intarna"
    exit 1
fi
echo "IntaRNA : $(IntaRNA --version 2>&1 | head -1)"
echo "Threads : ${THREADS}"
echo "Started : $(date)"
echo ""

# ============================================================================
# Paths
# ============================================================================

SRNA_DIR="results/sRNA_detection/07_intarna/sequences/srnas"
MRNA_DIR="results/sRNA_detection/07_intarna/sequences/mrnas"
OUT_DIR="results/sRNA_detection/07_intarna"

MRNA_DEG_EXP="${MRNA_DIR}/deg_exp"
MRNA_DEG_STA="${MRNA_DIR}/deg_sta"

INTERGENIC_EXP_IDS="${OUT_DIR}/intergenic_exp_ids.txt"
INTERGENIC_STA_IDS="${OUT_DIR}/intergenic_sta_ids.txt"

mkdir -p "${OUT_DIR}/priority1_leuC_antisense"
mkdir -p "${OUT_DIR}/priority2_intergenic_exp"
mkdir -p "${OUT_DIR}/priority3_intergenic_sta"

# ============================================================================
# Job counter
# ============================================================================

JOB_DONE=0
JOB_TOTAL=0

count_jobs() {
    local mrna_dir=$1
    if [ -d "$mrna_dir" ]; then
        find "$mrna_dir" -name "*.fasta" | wc -l
    else
        echo 0
    fi
}

# ============================================================================
# Core functions
# ============================================================================

run_intarna() {
    local srna_fasta=$1
    local mrna_fasta=$2
    local output_csv=$3
    local srna_name mrna_name energy

    srna_name=$(basename "$srna_fasta" .fasta)
    mrna_name=$(basename "$mrna_fasta" .fasta)

    [ ! -f "$srna_fasta" ] && { echo "  WARNING: Missing sRNA: $srna_fasta"; return; }
    [ ! -f "$mrna_fasta" ] && { echo "  WARNING: Missing mRNA: $mrna_fasta"; return; }

    IntaRNA \
        -q "$srna_fasta" \
        -t "$mrna_fasta" \
        --outMode=C \
        --out="$output_csv" \
        --seedBP=7 \
        --threads="${THREADS}" \
        2>/dev/null

    JOB_DONE=$(( JOB_DONE + 1 ))

    if [ -f "$output_csv" ] && [ -s "$output_csv" ]; then
        energy=$(tail -n 1 "$output_csv" | awk -F';' '{print $NF}')
        printf "  [%4d/%d]  %-28s -> %-20s  dG = %s kcal/mol\n" \
            "$JOB_DONE" "$JOB_TOTAL" "$srna_name" "$mrna_name" "$energy"
    else
        printf "  [%4d/%d]  %-28s -> %-20s  (no interaction)\n" \
            "$JOB_DONE" "$JOB_TOTAL" "$srna_name" "$mrna_name"
    fi
}

run_vs_mrna_dir() {
    local srna_fasta=$1
    local mrna_dir=$2
    local out_subdir=$3
    local srna_name n

    srna_name=$(basename "$srna_fasta" .fasta)

    [ ! -f "$srna_fasta" ] && { echo "  WARNING: Missing sRNA: $srna_fasta"; return; }

    if [ ! -d "$mrna_dir" ]; then
        echo "  WARNING: mRNA dir not found: $mrna_dir"
        echo "    Run 01_extract_sequences.py first."
        return
    fi

    n=$(find "$mrna_dir" -name "*.fasta" | wc -l)
    if [ "$n" -eq 0 ]; then
        echo "  WARNING: No mRNA targets in: $mrna_dir"
        return
    fi

    echo ""
    echo "  -- ${srna_name}  vs  ${n} targets --"

    for mrna_fasta in "$mrna_dir"/*.fasta; do
        [ -f "$mrna_fasta" ] || continue
        run_intarna \
            "$srna_fasta" \
            "$mrna_fasta" \
            "${out_subdir}/$(basename "$srna_fasta" .fasta)_vs_$(basename "$mrna_fasta" .fasta).csv"
    done
}

# ============================================================================
# Pre-calculate total jobs
# ============================================================================

N_EXP_IDS=0
N_STA_IDS=0
[ -f "$INTERGENIC_EXP_IDS" ] && N_EXP_IDS=$(grep -c '[^[:space:]]' "$INTERGENIC_EXP_IDS" 2>/dev/null || echo 0)
[ -f "$INTERGENIC_STA_IDS" ] && N_STA_IDS=$(grep -c '[^[:space:]]' "$INTERGENIC_STA_IDS" 2>/dev/null || echo 0)

N_DEG_EXP=$(count_jobs "$MRNA_DEG_EXP")
N_DEG_STA=$(count_jobs "$MRNA_DEG_STA")

JOB_TOTAL=$(( 2 + (N_EXP_IDS * N_DEG_EXP) + (N_STA_IDS * N_DEG_STA) ))

echo "Jobs queued:"
echo "  P1  2 sRNAs x 1 mRNA             =    2 jobs"
echo "  P2  ${N_EXP_IDS} sRNA(s) x ${N_DEG_EXP} DEG exp mRNAs  =  $(( N_EXP_IDS * N_DEG_EXP )) jobs"
echo "  P3  ${N_STA_IDS} sRNAs x ${N_DEG_STA} DEG sta mRNAs  = $(( N_STA_IDS * N_DEG_STA )) jobs"
echo "  --------------------------------"
echo "  TOTAL = ${JOB_TOTAL} jobs"
echo ""
echo "Estimated runtime: ~$(( JOB_TOTAL / 60 ))-$(( JOB_TOTAL * 2 / 60 )) minutes"
echo ""

# ============================================================================
# Priority 1 — Novel leuC antisense sRNAs
# ============================================================================

echo "================================================================================"
echo "PRIORITY 1 — Novel leuC Antisense sRNAs  (2 jobs)"
echo "================================================================================"
echo ""

for srna in exp_sRNA_00170 exp_sRNA_00172; do
    run_intarna \
        "${SRNA_DIR}/${srna}.fasta" \
        "${MRNA_DIR}/leuC.fasta" \
        "${OUT_DIR}/priority1_leuC_antisense/${srna}_vs_leuC.csv"
done

echo ""
echo "P1 done — elapsed: $(date)"
echo ""

# ============================================================================
# Priority 2 — Novel intergenic sRNAs, exponential
# ============================================================================

echo "================================================================================"
echo "PRIORITY 2 — Intergenic sRNAs (exp)  vs  ${N_DEG_EXP} DEG exp mRNAs  ($(( N_EXP_IDS * N_DEG_EXP )) jobs)"
echo "================================================================================"

if [ ! -f "$INTERGENIC_EXP_IDS" ] || [ ! -s "$INTERGENIC_EXP_IDS" ]; then
    echo "  WARNING: No IDs at: $INTERGENIC_EXP_IDS"
    echo "    Run 01_extract_sequences.py first."
else
    while IFS= read -r srna_id || [ -n "$srna_id" ]; do
        [ -z "$srna_id" ] && continue
        run_vs_mrna_dir \
            "${SRNA_DIR}/${srna_id}.fasta" \
            "${MRNA_DEG_EXP}" \
            "${OUT_DIR}/priority2_intergenic_exp"
    done < "$INTERGENIC_EXP_IDS"
fi

echo ""
echo "P2 done — elapsed: $(date)"
echo ""

# ============================================================================
# Priority 3 — Novel intergenic sRNAs, stationary
# ============================================================================

echo "================================================================================"
echo "PRIORITY 3 — Intergenic sRNAs (sta)  vs  ${N_DEG_STA} DEG sta mRNAs  ($(( N_STA_IDS * N_DEG_STA )) jobs)"
echo "================================================================================"
echo "  Note: largest block — ~$(( N_STA_IDS * N_DEG_STA / 60 ))-$(( N_STA_IDS * N_DEG_STA * 2 / 60 )) min estimated"
echo ""

if [ ! -f "$INTERGENIC_STA_IDS" ] || [ ! -s "$INTERGENIC_STA_IDS" ]; then
    echo "  WARNING: No IDs at: $INTERGENIC_STA_IDS"
    echo "    Run 01_extract_sequences.py first."
else
    while IFS= read -r srna_id || [ -n "$srna_id" ]; do
        [ -z "$srna_id" ] && continue
        run_vs_mrna_dir \
            "${SRNA_DIR}/${srna_id}.fasta" \
            "${MRNA_DEG_STA}" \
            "${OUT_DIR}/priority3_intergenic_sta"
    done < "$INTERGENIC_STA_IDS"
fi

echo ""
echo "P3 done — elapsed: $(date)"
echo ""

echo "================================================================================"
echo "All IntaRNA jobs complete — $(date)"
echo "================================================================================"
echo ""
echo "Next step:"
echo "  bash scripts/07_intarna/03_combine_results.sh"
echo ""
