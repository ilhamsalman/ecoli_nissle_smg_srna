Escherichia coli Nissle 1917 — sRNA responses under simulated microgravity

This repository contains scripts, metadata, and curated analysis outputs for the identification and characterization of small RNAs (sRNAs) in Escherichia coli Nissle 1917 under simulated microgravity (SMG) compared with normal gravity (NG) conditions.

The analysis focuses on:

genome-wide sRNA detection using ANNOgesic

condition-specific comparison across growth phase (exponential vs stationary)

prediction of regulatory RNA elements (riboswitches and RNA thermometers)


Dataset

Primary study

GEO: GSE147272
Transcriptional profiling of probiotic E. coli Nissle 1917 under simulated microgravity

RNA-seq runs used (GEO/SRA Run Selector):

Simulated microgravity (SMG)

Exponential phase: SRR11357510, SRR11357511

Stationary phase: SRR11357512, SRR11357513

Normal gravity (NG)

Exponential phase: SRR11357514, SRR11357515

Stationary phase: SRR11357516, SRR11357517

## Sample mapping

| sample_id | SRR         | growth_phase | gravity |
|-----------|-------------|--------------|---------|
| MG_exp_1  | SRR11357510 | exponential  | MG |
| MG_exp_2  | SRR11357511 | exponential  | MG |
| MG_sta_1  | SRR11357512 | stationary   | MG |
| MG_sta_2  | SRR11357513 | stationary   | MG |
| NG_exp_1  | SRR11357514 | exponential  | NG |
| NG_exp_2  | SRR11357515 | exponential  | NG |
| NG_sta_1  | SRR11357516 | stationary   | NG |
| NG_sta_2  | SRR11357517 | stationary   | NG |

## Workflow 

This repository documents the analysis steps performed on GSE147272 RNA-seq data for sRNA identification in *E. coli* Nissle 1917.

1. **FASTQ retrieval from SRA (SRR accessions)**
   - Script: `scripts/download_fastq.sh`


2. **Quality control on raw reads**
   - Tool(s): FastQC, MultiQC
 
3. **Read trimming**
   - Script: `scripts/trim_fastp.sh`
   - Tool: fastp


4. **Quality control on trimmed reads**
   - Tool(s): FastQC, MultiQC
   

5. **Reference preparation**
   - Reference: CP007799.1
     

6. **Read alignment**
   - Tool: Bowtie2
   
7. **BAM processing**
   - Tool: samtools
 

8. **sRNA annotation / downstream analysis**
   - Script: `scripts/annogesic_exp.sh`
   - Tool: ANNOgesic
   - Output: analysis outputs (see results/ or analysis/ )
