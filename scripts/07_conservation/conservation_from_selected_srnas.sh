#!/usr/bin/env bash
set -euo pipefail

# ============================================================

ROOT_DIR="$(pwd)"
PROJECT_REL="results/sRNA_detection"
SRNA_BASE="${ROOT_DIR}/${PROJECT_REL}/06_intarna/sequences/srnas"
EXP_DIR="${SRNA_BASE}/exp"
STA_DIR="${SRNA_BASE}/sta"
OUT_DIR="${ROOT_DIR}/${PROJECT_REL}/07_conservation"
GENOME_DIR="${OUT_DIR}/reference_genomes"
BLAST_DIR="${OUT_DIR}/blast_results"
SEQ_DIR="${OUT_DIR}/sequences"
QUERY_FASTA="${SEQ_DIR}/selected_srnas_all.fasta"
EXP_QUERY_FASTA="${SEQ_DIR}/selected_exp_srnas.fasta"
STA_QUERY_FASTA="${SEQ_DIR}/selected_sta_srnas.fasta"
SUMMARY_TXT="${OUT_DIR}/conservation_summary.txt"
SUMMARY_CSV="${OUT_DIR}/conservation_summary.csv"
PASSING_HITS_CSV="${OUT_DIR}/passing_hits_detailed.csv"

# Excel files with selected sRNAs.
EXP_XLSX="${1:-selected_exp_srnas.xlsx}"
STA_XLSX="${2:-selected_sta_srnas.xlsx}"

# Adjustable filters.
MIN_PIDENT="${MIN_PIDENT:-70}"
MIN_COVERAGE="${MIN_COVERAGE:-80}"
EVALUE="${EVALUE:-1e-5}"
MAX_TARGET_SEQS="${MAX_TARGET_SEQS:-50}"
WORD_SIZE="${WORD_SIZE:-7}"

mkdir -p "${OUT_DIR}" "${BLAST_DIR}" "${SEQ_DIR}"

############################################
# CHECK DEPENDENCIES
############################################
need_cmd() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "ERROR: Required command not found: $1" >&2
        exit 1
    }
}

need_cmd python3
need_cmd blastn

if ! python3 - <<'PY' >/dev/null 2>&1
import openpyxl
print('ok')
PY
then
    echo "ERROR: Python package 'openpyxl' is required." >&2
    echo "Install with: pip install openpyxl" >&2
    exit 1
fi

############################################
# RESOLVE INPUT FILES
############################################
resolve_file() {
    local f="$1"
    if [[ -f "$f" ]]; then
        readlink -f "$f"
        return 0
    fi
    if [[ -f "${ROOT_DIR}/$f" ]]; then
        readlink -f "${ROOT_DIR}/$f"
        return 0
    fi
    echo ""
}

EXP_XLSX_RESOLVED="$(resolve_file "${EXP_XLSX}")"
STA_XLSX_RESOLVED="$(resolve_file "${STA_XLSX}")"

if [[ -z "${EXP_XLSX_RESOLVED}" ]]; then
    echo "ERROR: EXP Excel file not found: ${EXP_XLSX}" >&2
    exit 1
fi
if [[ -z "${STA_XLSX_RESOLVED}" ]]; then
    echo "ERROR: STA Excel file not found: ${STA_XLSX}" >&2
    exit 1
fi
if [[ ! -d "${EXP_DIR}" ]]; then
    echo "ERROR: EXP FASTA directory not found: ${EXP_DIR}" >&2
    exit 1
fi
if [[ ! -d "${STA_DIR}" ]]; then
    echo "ERROR: STA FASTA directory not found: ${STA_DIR}" >&2
    exit 1
fi
if [[ ! -d "${GENOME_DIR}" ]]; then
    echo "ERROR: Reference genome directory not found: ${GENOME_DIR}" >&2
    exit 1
fi

echo "=========================================="
echo " Conservation analysis from selected sRNAs"
echo "=========================================="
echo "EXP Excel      : ${EXP_XLSX_RESOLVED}"
echo "STA Excel      : ${STA_XLSX_RESOLVED}"
echo "EXP FASTA dir  : ${EXP_DIR}"
echo "STA FASTA dir  : ${STA_DIR}"
echo "Genome DB dir  : ${GENOME_DIR}"
echo "Output         : ${OUT_DIR}"
echo "Min % identity : ${MIN_PIDENT}"
echo "Min coverage   : ${MIN_COVERAGE}"
echo "E-value        : ${EVALUE}"
echo ""

############################################
# STEP 1. EXTRACT SELECTED IDs FROM EXCEL
############################################
echo "[1/5] Extracting selected sRNA IDs from Excel..."

EXP_LIST="${OUT_DIR}/selected_exp_ids.txt"
STA_LIST="${OUT_DIR}/selected_sta_ids.txt"
ALL_LIST="${OUT_DIR}/selected_all_ids.tsv"

python3 - "${EXP_XLSX_RESOLVED}" "${STA_XLSX_RESOLVED}" "${EXP_LIST}" "${STA_LIST}" "${ALL_LIST}" <<'PY'
import sys
import re
from openpyxl import load_workbook

exp_xlsx, sta_xlsx, exp_out, sta_out, all_out = sys.argv[1:]
pattern = re.compile(r'^(exp|sta)_sRNA_\d+$', re.IGNORECASE)

def extract_ids(path):
    wb = load_workbook(path, data_only=True)
    ws = wb[wb.sheetnames[0]]
    ids = []
    seen = set()
    for row in ws.iter_rows(values_only=True):
        for cell in row:
            if cell is None:
                continue
            value = str(cell).strip()
            if pattern.match(value) and value not in seen:
                ids.append(value)
                seen.add(value)
    return ids

exp_ids = extract_ids(exp_xlsx)
sta_ids = extract_ids(sta_xlsx)

with open(exp_out, 'w') as f:
    for x in exp_ids:
        f.write(x + '\n')
with open(sta_out, 'w') as f:
    for x in sta_ids:
        f.write(x + '\n')
with open(all_out, 'w') as f:
    for x in exp_ids:
        f.write(f'exp\t{x}\n')
    for x in sta_ids:
        f.write(f'sta\t{x}\n')

print(f'EXP selected IDs: {len(exp_ids)}')
print(f'STA selected IDs: {len(sta_ids)}')
PY

############################################
# STEP 2. BUILD QUERY FASTA FILES
############################################
echo ""
echo "[2/5] Building query FASTA files from selected IDs..."

: > "${QUERY_FASTA}"
: > "${EXP_QUERY_FASTA}"
: > "${STA_QUERY_FASTA}"
MISSING_REPORT="${OUT_DIR}/missing_fastas.txt"
: > "${MISSING_REPORT}"

normalize_fasta_to_id() {
    local fasta="$1"
    local seq_id="$2"
    awk -v n="${seq_id}" '
        BEGIN{printed=0}
        /^>/ {
            if (!printed) {
                print ">" n
                printed=1
            }
            next
        }
        {
            gsub(/\r/, "")
            if (length($0) > 0) print
        }
    ' "${fasta}"
}

add_fasta() {
    local phase="$1"
    local id="$2"
    local base_dir="$3"
    local phase_out="$4"

    local fasta=""
    for ext in fasta fa fna fasta.txt; do
        if [[ -f "${base_dir}/${id}.${ext}" ]]; then
            fasta="${base_dir}/${id}.${ext}"
            break
        fi
    done

    if [[ -z "${fasta}" ]]; then
        fasta=$(find "${base_dir}" -maxdepth 1 -type f \( -iname "${id}.fasta" -o -iname "${id}.fa" -o -iname "${id}.fna" \) | head -n 1 || true)
    fi

    if [[ -n "${fasta}" && -f "${fasta}" ]]; then
        normalize_fasta_to_id "${fasta}" "${id}" >> "${phase_out}"
        normalize_fasta_to_id "${fasta}" "${id}" >> "${QUERY_FASTA}"
        echo "Added ${phase}: ${id}"
    else
        printf "%s\t%s\tMISSING\n" "${phase}" "${id}" >> "${MISSING_REPORT}"
        echo "WARNING: FASTA not found for ${phase} / ${id}" >&2
    fi
}

while IFS= read -r id; do
    [[ -n "$id" ]] && add_fasta "exp" "$id" "${EXP_DIR}" "${EXP_QUERY_FASTA}"
done < "${EXP_LIST}"

while IFS= read -r id; do
    [[ -n "$id" ]] && add_fasta "sta" "$id" "${STA_DIR}" "${STA_QUERY_FASTA}"
done < "${STA_LIST}"

TOTAL_SEQ=$(grep -c '^>' "${QUERY_FASTA}" || true)
EXP_SEQ=$(grep -c '^>' "${EXP_QUERY_FASTA}" || true)
STA_SEQ=$(grep -c '^>' "${STA_QUERY_FASTA}" || true)

echo ""
echo "EXP sequences added : ${EXP_SEQ}"
echo "STA sequences added : ${STA_SEQ}"
echo "Total sequences     : ${TOTAL_SEQ}"

if [[ "${TOTAL_SEQ}" -eq 0 ]]; then
    echo "ERROR: No FASTA sequences were collected." >&2
    echo "Check selected IDs and file names in ${EXP_DIR} and ${STA_DIR}." >&2
    exit 1
fi

############################################
# STEP 3. VALIDATE EXISTING BLAST DATABASES
############################################
echo ""
echo "[3/5] Checking existing BLAST databases..."

check_db_prefix() {
    local prefix="$1"
    if [[ -f "${prefix}.nhr" || -f "${prefix}.nin" || -f "${prefix}.nsq" ]]; then
        return 0
    fi
    return 1
}

for strain in K12 O157 CFT073; do
    DB_PREFIX="${GENOME_DIR}/${strain}"
    if ! check_db_prefix "${DB_PREFIX}"; then
        echo "ERROR: BLAST DB for ${strain} not found under ${GENOME_DIR}" >&2
        echo "Expected files such as: ${DB_PREFIX}.nhr / ${DB_PREFIX}.nin / ${DB_PREFIX}.nsq" >&2
        exit 1
    fi
    echo "${strain}: existing BLAST DB found"
done

############################################
# STEP 4. RUN BLAST
############################################
echo ""
echo "[4/5] Running BLASTn with multi-hit output..."

for strain in K12 O157 CFT073; do
    DB_PREFIX="${GENOME_DIR}/${strain}"
    OUTFILE="${BLAST_DIR}/${strain}_blast.tsv"
    echo "BLASTn vs ${strain}"
    blastn \
        -query "${QUERY_FASTA}" \
        -db "${DB_PREFIX}" \
        -out "${OUTFILE}" \
        -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore' \
        -evalue "${EVALUE}" \
        -word_size "${WORD_SIZE}" \
        -perc_identity "${MIN_PIDENT}" \
        -max_target_seqs "${MAX_TARGET_SEQS}"
    echo "Raw hits in ${strain}: $(wc -l < "${OUTFILE}")"
done

############################################
# STEP 5. BUILD FILTERED SUMMARY
############################################
echo ""
echo "[5/5] Writing filtered conservation summary..."

python3 - \
    "${ALL_LIST}" \
    "${BLAST_DIR}" \
    "${SUMMARY_TXT}" \
    "${SUMMARY_CSV}" \
    "${PASSING_HITS_CSV}" \
    "${MIN_PIDENT}" \
    "${MIN_COVERAGE}" <<'PY'
import csv
import sys
from pathlib import Path

all_list, blast_dir, summary_txt, summary_csv, passing_csv, min_pident, min_coverage = sys.argv[1:]
blast_dir = Path(blast_dir)
min_pident = float(min_pident)
min_coverage = float(min_coverage)

queries = []
with open(all_list) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        phase, srna = line.split('\t')
        queries.append((phase, srna))

strains = ["K12", "O157", "CFT073"]
all_hits = {s: {} for s in strains}

def strand_from_coords(sstart, send):
    return '+' if int(sstart) <= int(send) else '-'

def normalized_coords(sstart, send):
    a, b = int(sstart), int(send)
    return (a, b) if a <= b else (b, a)

def parse_blast(path):
    hits = {}
    if not path.exists():
        return hits
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 11:
                continue
            qseqid, sseqid, pident, length, qlen, qstart, qend, sstart, send, evalue, bitscore = parts[:11]
            pident = float(pident)
            length = int(length)
            qlen = int(qlen)
            qstart = int(qstart)
            qend = int(qend)
            sstart = int(sstart)
            send = int(send)
            bitscore = float(bitscore)
            coverage = 100.0 * length / qlen if qlen else 0.0
            hit = {
                'qseqid': qseqid,
                'sseqid': sseqid,
                'pident': pident,
                'length': length,
                'qlen': qlen,
                'qstart': qstart,
                'qend': qend,
                'coverage': coverage,
                'sstart': sstart,
                'send': send,
                'strand': strand_from_coords(sstart, send),
                'coord_min': normalized_coords(sstart, send)[0],
                'coord_max': normalized_coords(sstart, send)[1],
                'coords': f"{min(sstart, send)}-{max(sstart, send)}",
                'evalue': evalue,
                'bitscore': bitscore,
            }
            hits.setdefault(qseqid, []).append(hit)
    return hits

for strain in strains:
    all_hits[strain] = parse_blast(blast_dir / f"{strain}_blast.tsv")

def passing_hits_for_query(strain, srna):
    raw = all_hits[strain].get(srna, [])
    passing = [
        h for h in raw
        if h['pident'] >= min_pident and h['coverage'] >= min_coverage
    ]
    passing.sort(key=lambda h: (h['coverage'], h['bitscore'], h['pident'], -h['evalue'].count('e')), reverse=True)
    return passing

# Write all passing hits.
with open(passing_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow([
        'phase', 'sRNA', 'strain', 'rank', 'subject', 'pident', 'coverage', 'length', 'qlen',
        'qstart', 'qend', 'sstart', 'send', 'strand', 'coords', 'evalue', 'bitscore'
    ])
    for phase, srna in queries:
        for strain in strains:
            passing = passing_hits_for_query(strain, srna)
            for rank, h in enumerate(passing, start=1):
                writer.writerow([
                    phase, srna, strain, rank, h['sseqid'],
                    f"{h['pident']:.2f}", f"{h['coverage']:.2f}", h['length'], h['qlen'],
                    h['qstart'], h['qend'], h['sstart'], h['send'], h['strand'], h['coords'],
                    h['evalue'], f"{h['bitscore']:.2f}"
                ])

rows = []
for phase, srna in queries:
    row = {'phase': phase, 'sRNA': srna}
    valid_strains = 0
    all_unambiguous = True
    mean_ident_vals = []
    mean_cov_vals = []
    notes = []

    for strain in strains:
        passing = passing_hits_for_query(strain, srna)
        row[f'{strain}_nhits'] = len(passing)
        if passing:
            best = passing[0]
            row[f'{strain}_pident'] = f"{best['pident']:.1f}%"
            row[f'{strain}_coverage'] = f"{best['coverage']:.1f}%"
            row[f'{strain}_coords'] = best['coords']
            row[f'{strain}_strand'] = best['strand']
            row[f'{strain}_subject'] = best['sseqid']
            row[f'{strain}_evalue'] = best['evalue']
            valid_strains += 1
            mean_ident_vals.append(best['pident'])
            mean_cov_vals.append(best['coverage'])
            if len(passing) > 1:
                all_unambiguous = False
                notes.append(f'{strain}:multiple_passing_hits={len(passing)}')
        else:
            row[f'{strain}_pident'] = 'N/A'
            row[f'{strain}_coverage'] = 'N/A'
            row[f'{strain}_coords'] = 'N/A'
            row[f'{strain}_strand'] = 'N/A'
            row[f'{strain}_subject'] = 'N/A'
            row[f'{strain}_evalue'] = 'N/A'

    if valid_strains == 3 and all_unambiguous:
        status = 'Conserved (high-confidence)'
    elif valid_strains == 3:
        status = 'Conserved (ambiguous loci)'
    elif valid_strains >= 1:
        status = 'Partially conserved'
    else:
        status = 'No supported hit'

    row['status'] = status
    row['mean_identity'] = f"{(sum(mean_ident_vals)/len(mean_ident_vals)):.1f}%" if mean_ident_vals else 'N/A'
    row['mean_coverage'] = f"{(sum(mean_cov_vals)/len(mean_cov_vals)):.1f}%" if mean_cov_vals else 'N/A'
    row['notes'] = ';'.join(notes) if notes else '-'
    rows.append(row)

with open(summary_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow([
        'phase', 'sRNA',
        'K12_pident', 'K12_coverage', 'K12_coords', 'K12_strand', 'K12_subject', 'K12_evalue', 'K12_nhits',
        'O157_pident', 'O157_coverage', 'O157_coords', 'O157_strand', 'O157_subject', 'O157_evalue', 'O157_nhits',
        'CFT073_pident', 'CFT073_coverage', 'CFT073_coords', 'CFT073_strand', 'CFT073_subject', 'CFT073_evalue', 'CFT073_nhits',
        'mean_identity', 'mean_coverage', 'status', 'notes'
    ])
    for r in rows:
        writer.writerow([
            r['phase'], r['sRNA'],
            r['K12_pident'], r['K12_coverage'], r['K12_coords'], r['K12_strand'], r['K12_subject'], r['K12_evalue'], r['K12_nhits'],
            r['O157_pident'], r['O157_coverage'], r['O157_coords'], r['O157_strand'], r['O157_subject'], r['O157_evalue'], r['O157_nhits'],
            r['CFT073_pident'], r['CFT073_coverage'], r['CFT073_coords'], r['CFT073_strand'], r['CFT073_subject'], r['CFT073_evalue'], r['CFT073_nhits'],
            r['mean_identity'], r['mean_coverage'], r['status'], r['notes']
        ])

with open(summary_txt, 'w') as f:
    f.write("==========================================================================\n")
    f.write("Filtered conservation summary for selected exp/sta sRNAs\n")
    f.write("==========================================================================\n\n")
    f.write(f"Total selected sequences: {len(rows)}\n")
    f.write(f"Minimum percent identity: {min_pident:.1f}%\n")
    f.write(f"Minimum query coverage : {min_coverage:.1f}%\n\n")
    f.write("{:<6}  {:<18} {:<12} {:<12} {:<12} {:<10} {:<28}\n".format(
        'Phase', 'sRNA', 'K12', 'O157', 'CFT073', 'MeanCov', 'Status'
    ))
    f.write("-" * 116 + "\n")
    for r in rows:
        f.write("{:<6}  {:<18} {:<12} {:<12} {:<12} {:<10} {:<28}\n".format(
            r['phase'], r['sRNA'],
            r['K12_pident'], r['O157_pident'], r['CFT073_pident'],
            r['mean_coverage'], r['status']
        ))
    f.write("\nInterpretation:\n")
    f.write("  Conserved (high-confidence) = one passing hit per strain after identity/coverage filtering\n")
    f.write("  Conserved (ambiguous loci)  = passing hit in all strains, but >=2 passing loci in at least one strain\n")
    f.write("  Partially conserved         = passing hit in 1-2 strains only\n")
    f.write("  No supported hit            = no hit passes the identity+coverage filters\n\n")
    f.write("Notes:\n")
    f.write("  - Coordinates and strand are reported for manual locus inspection.\n")
    f.write("  - This script does not perform annotation-based synteny validation.\n")
    f.write("  - Use passing_hits_detailed.csv to inspect all passing loci per strain.\n")

print(f'Summary written       : {summary_txt}')
print(f'Summary CSV written   : {summary_csv}')
print(f'Detailed hits written : {passing_csv}')
PY

echo ""
echo "Done. Main outputs:"
echo "  ${SUMMARY_TXT}"
echo "  ${SUMMARY_CSV}"
echo "  ${PASSING_HITS_CSV}"
echo "  ${QUERY_FASTA}"
echo "  ${MISSING_REPORT}"
echo ""
echo "=========================================="
echo " CONSERVATION ANALYSIS COMPLETE"
echo "=========================================="
