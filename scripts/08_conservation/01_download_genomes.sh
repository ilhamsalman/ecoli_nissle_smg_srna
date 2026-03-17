#!/bin/bash
# Step 1: Download E. coli Reference Genomes

echo "=========================================="
echo "Downloading E. coli Reference Genomes"
echo "=========================================="

# Paths
CONS_DIR="results/sRNA_detection/08_conservation"
GENOME_DIR="${CONS_DIR}/reference_genomes"

# Create directories
mkdir -p "${GENOME_DIR}"
mkdir -p "${CONS_DIR}/blast_results"

echo ""
echo "Downloading genomes to: ${GENOME_DIR}"
echo ""

# Download K-12 MG1655
echo "1. Downloading E. coli K-12 MG1655..."
wget -q -O "${GENOME_DIR}/K12.fna.gz" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
gunzip -f "${GENOME_DIR}/K12.fna.gz"
echo "   ✓ Downloaded K-12"

# Download O157:H7
echo "2. Downloading E. coli O157:H7..."
wget -q -O "${GENOME_DIR}/O157.fna.gz" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2/GCF_000008865.2_ASM886v2_genomic.fna.gz"
gunzip -f "${GENOME_DIR}/O157.fna.gz"
echo "   ✓ Downloaded O157:H7"

# Download CFT073
echo "3. Downloading E. coli CFT073..."
wget -q -O "${GENOME_DIR}/CFT073.fna.gz" \
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/445/GCF_000007445.1_ASM744v1/GCF_000007445.1_ASM744v1_genomic.fna.gz"
gunzip -f "${GENOME_DIR}/CFT073.fna.gz"
echo "   ✓ Downloaded CFT073"

echo ""
echo "Creating BLAST databases..."
echo ""

# Create BLAST databases
makeblastdb -in "${GENOME_DIR}/K12.fna" -dbtype nucl -out "${GENOME_DIR}/K12" > /dev/null 2>&1
echo "   ✓ K12 database created"

makeblastdb -in "${GENOME_DIR}/O157.fna" -dbtype nucl -out "${GENOME_DIR}/O157" > /dev/null 2>&1
echo "   ✓ O157 database created"

makeblastdb -in "${GENOME_DIR}/CFT073.fna" -dbtype nucl -out "${GENOME_DIR}/CFT073" > /dev/null 2>&1
echo "   ✓ CFT073 database created"

echo ""
echo "=========================================="
echo "STEP 1 COMPLETE"
echo "=========================================="
echo ""
echo "Downloaded:"
ls -lh "${GENOME_DIR}/"*.fna
echo ""
echo "Next: bash scripts/08_conservation/02_blast_and_summarize.sh"
