#!/bin/bash
# ============================================================
# Conservation Analysis — Updated to include novel intergenic sRNAs
# Covers: SgrS, rydB, rdlD (known) +
#         exp_sRNA_00170, exp_sRNA_00172 (leuC antisense) +
#         exp_sRNA_02152 (intergenic exp) +
#         sta_sRNA_00458, sta_sRNA_00470, sta_sRNA_00504,
#         sta_sRNA_00850, sta_sRNA_00890, sta_sRNA_00946 (intergenic sta)
# ============================================================

set -euo pipefail

echo "=========================================="
echo " Conservation Analysis (Full sRNA Set)"
echo "=========================================="

# ── Paths ──────────────────────────────────────────────────
PROJECT="results/sRNA_detection"
INTARNA_SRNAS="${PROJECT}/07_intarna/sequences/srnas"
KNOWN_SRNAS="${PROJECT}/05_known_srna"          # where SgrS/rydB/rdlD FASTAs live
CONS_DIR="${PROJECT}/08_conservation"
GENOME_DIR="${CONS_DIR}/reference_genomes"
BLAST_DIR="${CONS_DIR}/blast_results"
QUERY="${CONS_DIR}/all_srnas.fasta"

mkdir -p "${GENOME_DIR}" "${BLAST_DIR}" "${CONS_DIR}/sequences"

# ── Step 1: Build query FASTA (intergenic sRNAs only — known/antisense
#            already analysed and hard-coded below) ──────────────────────────
echo ""
echo "[1/4] Building query FASTA for new intergenic sRNAs..."

> "${QUERY}"

NEW_SRNAS=(
    "exp_sRNA_02152"
    "sta_sRNA_00458"
    "sta_sRNA_00470"
    "sta_sRNA_00504"
    "sta_sRNA_00850"
    "sta_sRNA_00890"
    "sta_sRNA_00946"
)

for name in "${NEW_SRNAS[@]}"; do
    fasta=$(find "${INTARNA_SRNAS}" -iname "${name}.fasta" 2>/dev/null | head -1)
    if [ -n "$fasta" ]; then
        awk -v n="${name}" '/^>/{print ">" n; next} {print}' "$fasta" >> "${QUERY}"
        echo "   Added: ${name}"
    else
        echo "   WARNING: ${name}.fasta not found in ${INTARNA_SRNAS} — skipping"
    fi
done

# Count total sequences
NSEQ=$(grep -c "^>" "${QUERY}" || true)
echo ""
echo "   Total sequences in query: ${NSEQ}"

if [ "$NSEQ" -eq 0 ]; then
    echo "ERROR: No sequences found. Check that 01_extract_sequences.py ran successfully."
    exit 1
fi

# ── Step 2: Check / download reference genomes ────────────
echo ""
echo "[2/4] Checking reference genome BLAST databases..."

# Accessions
declare -A GENOMES
GENOMES["K12"]="U00096.3"
GENOMES["O157"]="AE005174.2"
GENOMES["CFT073"]="AE014075.1"

for strain in K12 O157 CFT073; do
    DB="${GENOME_DIR}/${strain}"
    if [ -f "${DB}.nhr" ] || [ -f "${DB}.nsq" ] || [ -f "${DB}.nsi" ]; then
        echo "   ${strain}: BLAST DB already exists — skipping download"
    else
        ACC="${GENOMES[$strain]}"
        echo "   ${strain}: Downloading ${ACC} from NCBI..."
        GENOME_FA="${GENOME_DIR}/${strain}.fa"
        # Try efetch first, fall back to wget
        if command -v efetch &>/dev/null; then
            efetch -db nucleotide -id "${ACC}" -format fasta > "${GENOME_FA}"
        else
            wget -q "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${ACC}&rettype=fasta&retmode=text" \
                 -O "${GENOME_FA}"
        fi
        echo "   ${strain}: Building BLAST DB..."
        makeblastdb -in "${GENOME_FA}" -dbtype nucl -out "${DB}" -title "${strain}" -parse_seqids
        echo "   ${strain}: Done"
    fi
done

# ── Step 3: Run BLAST ──────────────────────────────────────
echo ""
echo "[3/4] Running BLAST searches..."

for strain in K12 O157 CFT073; do
    DB="${GENOME_DIR}/${strain}"
    OUT="${BLAST_DIR}/${strain}_blast.txt"
    echo "   BLASTn vs ${strain}..."
    blastn -query  "${QUERY}" \
           -db     "${DB}" \
           -out    "${OUT}" \
           -outfmt "6 qseqid sseqid pident length evalue bitscore" \
           -evalue  1e-5 \
           -word_size 7 \
           -perc_identity 70 \
           -max_target_seqs 1 2>/dev/null
    NHITS=$(wc -l < "${OUT}")
    echo "   ${strain}: ${NHITS} hits"
done

# ── Step 4: Build summary table ───────────────────────────
echo ""
echo "[4/4] Building summary tables..."

TXT="${CONS_DIR}/conservation_summary.txt"
CSV="${CONS_DIR}/conservation_summary.csv"

# ── Hard-coded previous results (from screenshot — do not re-BLAST these) ──
# Format: sRNA  K12   O157  CFT073  Status
declare -A PREV_K12 PREV_O157 PREV_CFT PREV_STATUS
PREV_K12["SgrS"]="88.4%";   PREV_O157["SgrS"]="99.9%";   PREV_CFT["SgrS"]="100.0%";   PREV_STATUS["SgrS"]="Highly conserved"
PREV_K12["rydB"]="94.6%";   PREV_O157["rydB"]="94.2%";   PREV_CFT["rydB"]="100.0%";   PREV_STATUS["rydB"]="Highly conserved"
PREV_K12["rdlD"]="97.0%";   PREV_O157["rdlD"]="97.3%";   PREV_CFT["rdlD"]="100.0%";   PREV_STATUS["rdlD"]="Highly conserved"
PREV_K12["exp_sRNA_00170"]="97.5%"; PREV_O157["exp_sRNA_00170"]="96.9%"; PREV_CFT["exp_sRNA_00170"]="100.0%"; PREV_STATUS["exp_sRNA_00170"]="Highly conserved"
PREV_K12["exp_sRNA_00172"]="98.5%"; PREV_O157["exp_sRNA_00172"]="97.1%"; PREV_CFT["exp_sRNA_00172"]="100.0%"; PREV_STATUS["exp_sRNA_00172"]="Highly conserved"

PREV_SRNAS=("SgrS" "rydB" "rdlD" "exp_sRNA_00170" "exp_sRNA_00172")

# ── New intergenic sRNAs to BLAST ──────────────────────────
NEW_SRNAS=(
    "exp_sRNA_02152"
    "sta_sRNA_00458"
    "sta_sRNA_00470"
    "sta_sRNA_00504"
    "sta_sRNA_00850"
    "sta_sRNA_00890"
    "sta_sRNA_00946"
)

# Header
cat > "${TXT}" << 'EOF'
================================================================
sRNA Conservation Across E. coli Strains
================================================================

  Group                    | sRNA             | K-12      | O157:H7   | CFT073    | Status
  -------------------------|------------------|-----------|-----------|-----------|------------------
EOF

echo "Group,sRNA,K12_pident,O157_pident,CFT073_pident,Status" > "${CSV}"

# ── Helper: write one row ──────────────────────────────────
write_row() {
    local grp="$1" srna="$2" k12="$3" o157="$4" cft="$5" status="$6"
    printf "  %-25s | %-16s | %-9s | %-9s | %-9s | %s\n" \
        "$grp" "$srna" "$k12" "$o157" "$cft" "$status" >> "${TXT}"
    # strip % for CSV
    local k12n="${k12//%/}"; local o157n="${o157//%/}"; local cftn="${cft//%/}"
    echo "${grp},${srna},${k12n},${o157n},${cftn},${status}" >> "${CSV}"
}

# ── Write previous results (hard-coded) ───────────────────
echo "  # Known + antisense sRNAs (previous analysis)" >> "${TXT}"
for srna in "${PREV_SRNAS[@]}"; do
    case "$srna" in
        SgrS|rydB|rdlD)          grp="Known (stationary)" ;;
        exp_sRNA_00170|exp_sRNA_00172) grp="Novel antisense (exp)" ;;
    esac
    write_row "$grp" "$srna" \
        "${PREV_K12[$srna]}" "${PREV_O157[$srna]}" "${PREV_CFT[$srna]}" \
        "${PREV_STATUS[$srna]}"
done

echo "" >> "${TXT}"
echo "  # Novel intergenic sRNAs (new analysis)" >> "${TXT}"

# ── Write new intergenic results from BLAST ────────────────
for srna in "${NEW_SRNAS[@]}"; do
    # Skip if not in query (sequence extraction may have failed)
    if ! grep -q "^>${srna}$" "${QUERY}" 2>/dev/null; then
        echo "   SKIPPING ${srna}: not found in query FASTA"
        continue
    fi

    k12=$(awk  -v s="${srna}" '$1==s {printf "%.1f%%", $3; exit}' "${BLAST_DIR}/K12_blast.txt")
    o157=$(awk -v s="${srna}" '$1==s {printf "%.1f%%", $3; exit}' "${BLAST_DIR}/O157_blast.txt")
    cft=$(awk  -v s="${srna}" '$1==s {printf "%.1f%%", $3; exit}' "${BLAST_DIR}/CFT073_blast.txt")

    k12="${k12:-N/A}"; o157="${o157:-N/A}"; cft="${cft:-N/A}"

    count=0
    [[ "$k12"  != "N/A" ]] && count=$((count+1))
    [[ "$o157" != "N/A" ]] && count=$((count+1))
    [[ "$cft"  != "N/A" ]] && count=$((count+1))

    if   [ $count -eq 3 ]; then status="Highly conserved"
    elif [ $count -ge 1 ]; then status="Partially conserved"
    else                        status="Nissle-specific"
    fi

    case "$srna" in
        exp_sRNA_02152)                    grp="Novel intergenic (exp)" ;;
        sta_sRNA_00458|sta_sRNA_00504|sta_sRNA_00946) grp="Novel intergenic (sta)-TAT" ;;
        sta_sRNA_00470)                    grp="Novel intergenic (sta)-stpA" ;;
        *)                                 grp="Novel intergenic (sta)" ;;
    esac

    write_row "$grp" "$srna" "$k12" "$o157" "$cft" "$status"
done

cat >> "${TXT}" << 'EOF'

================================================================
INTERPRETATION:
  Highly conserved    = Hit in all 3 strains (identity shown)
  Partially conserved = Hit in 1-2 strains
  Nissle-specific     = No hit in any comparison strain (evalue < 1e-5)

TAT  = Predicted TAT secretion pathway regulator
stpA = Predicted global nucleoid regulator target
================================================================
EOF

echo ""
cat "${TXT}"
echo ""
echo "CSV saved: ${CSV}"
echo ""
echo "=========================================="
echo " CONSERVATION ANALYSIS COMPLETE"
echo "=========================================="
