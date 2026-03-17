# Conservation Analysis

Check if sRNAs are conserved across *E. coli* strains.

---

## Quick Run

```bash
cd ~/projects/ecoli_nissle_smg_srna

# Prepare sequences
mkdir -p results/sRNA_detection/08_conservation/sequences
cp results/sRNA_detection/07_intarna/sequences/srnas/*.fasta \
   results/sRNA_detection/08_conservation/sequences/
cat results/sRNA_detection/08_conservation/sequences/*.fasta > \
    results/sRNA_detection/08_conservation/all_srnas.fasta
sed -i 's/ .*//' results/sRNA_detection/08_conservation/all_srnas.fasta

# Run analysis
bash scripts/08_conservation/01_download_genomes.sh
bash scripts/08_conservation/02_blast_and_summarize.sh

# View results
cat results/sRNA_detection/08_conservation/conservation_summary.txt
```

---

## Output

**File:** `results/sRNA_detection/08_conservation/conservation_summary.txt`

Shows identity (%) of each sRNA against:
- K-12 MG1655 (model strain)
- O157:H7 (pathogenic)
- CFT073 (uropathogenic)

**Status:**
- Highly conserved = >90% in multiple strains
- Partially conserved = >90% in one strain  
- Not conserved = Nissle-specific


```
Conservation analysis by BLASTN against E. coli K-12, O157:H7, 
and CFT073 genomes (E-value ≤ 1e-5).
```

