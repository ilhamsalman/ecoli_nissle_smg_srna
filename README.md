# *Escherichia coli* Nissle 1917 — RNA regulatory responses under simulated microgravity

This repository contains analysis scripts, sample metadata, and selected outputs for a computational study of small regulatory RNAs (sRNAs) and other RNA regulatory elements in *Escherichia coli* Nissle 1917 (EcN) under simulated microgravity (SMG) versus normal gravity (NG). Exponential- and stationary-phase samples are analysed separately.

The analysis includes genome-wide sRNA prediction using ANNOgesic, sRNA expression analysis using featureCounts and limma–voom, comparison with previously annotated sRNAs, IntaRNA prediction of candidate sRNA–mRNA interactions, cross-strain conservation analysis, and identification of candidate RNA thermometer/riboswitch-module elements and their phase-matched DEG associations.

**Interpretation:** The reported sRNA expression candidates meet a **nominal** P-value threshold; none remained significant after multiple-testing correction. The RNA–RNA interactions and cis-regulatory associations are **computational predictions**, not experimentally demonstrated regulatory interactions.

## Dataset and reference genome

- **GEO:** [GSE147272](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147272) — source RNA-seq dataset from Yim *et al.* (2020).
- **EcN reference genome:** GenBank accession **CP007799.1** (the alignment accession is also written **NZ_CP007799.1**).
- **Experimental groups:** SMG and NG, with two RNA-seq runs per condition in each growth phase.

### RNA-seq run mapping

| Sample ID | SRA run | Growth phase | Condition |
| --- | --- | --- | --- |
| MG_exp_1 | SRR11357510 | Exponential | SMG |
| MG_exp_2 | SRR11357511 | Exponential | SMG |
| MG_sta_1 | SRR11357512 | Stationary | SMG |
| MG_sta_2 | SRR11357513 | Stationary | SMG |
| NG_exp_1 | SRR11357514 | Exponential | NG |
| NG_exp_2 | SRR11357515 | Exponential | NG |
| NG_sta_1 | SRR11357516 | Stationary | NG |
| NG_sta_2 | SRR11357517 | Stationary | NG |

The `MG_` sample-ID prefix denotes the simulated-microgravity samples; `NG_` denotes normal gravity. The repository also contains `meta/samples_nissle_GSE147272.tsv`.

## Analysis overview

| Stage | Location | Purpose |
| --- | --- | --- |
| 00 | `scripts/00_pipeline/` | Download, trim, align, and assess strandedness of RNA-seq reads |
| 01 | `scripts/01_annogesic/` | Generate coverage tracks and run ANNOgesic |
| 02 | `scripts/02_filter/` | Filter predicted sRNA loci |
| 03 | `scripts/03_shared_srna/` | Build SAF files and compare sRNA loci across growth phases |
| 04 | `scripts/04_srna_de/` | Feature counting and phase-specific limma–voom analysis |
| 05 | `scripts/05_annotate_known_srna/` | Compare candidate sRNAs with BSRD |
| 06 | `scripts/06_intarna/` | Prepare selected sRNA/mRNA sequences, predict and filter interactions |
| 07 | `scripts/07_conservation/` | Search selected sRNAs for cross-strain conservation |
| 08 | `scripts/08_cis_deg/` | Compare ANNOgesic-associated CDS locus tags with phase-specific DEGs |
| 09 | `scripts/09_figures/` | Reproduce the repository's main plotted figures |

See **[WORKFLOW.md](WORKFLOW.md)** for script names, required inputs, execution notes, and output locations. Run commands from the repository root unless a script says otherwise. Some stages require locally prepared reads, reference sequences, annotations, or other dependencies not distributed here.

## Key exploratory results

- **Predicted/filtered sRNAs:** 3,752 exponential-phase and 519 stationary-phase loci after the initial sRNA filters.
- **Tested for differential expression:** 697 exponential-phase and 317 stationary-phase loci after expression filtering.
- **Nominal candidate expression changes:** 44 exponential-phase and 35 stationary-phase sRNA records (79 records representing **78 distinct genomic loci**). **None survived FDR correction.**
- **IntaRNA prioritisation:** Eight selected exponential-phase and three selected stationary-phase intergenic sRNAs were assessed against phase-specific curated DEG targets. The interaction results are predictions, not validated binding events.
- **ANNOgesic RNA thermometer-module predictions:** 43 exponential and 14 stationary records.
- **ANNOgesic riboswitch-module predictions:** 27 exponential and 11 stationary records.
- **Same-phase cis-element–DEG associations:** Two locus-tag-based associations were identified: an exponential-phase TPP-riboswitch-module prediction associated with *rcnA* and a stationary-phase Phe_leader prediction associated with *manA*. These are annotation-based associations; physical overlap, ligand-dependent regulation, and causality were not demonstrated.

The module counts above include the Rfam classes emitted by ANNOgesic; they should **not** all be interpreted as validated thermometers or riboswitches. In particular, **Phe_leader** is a cis-regulatory leader RNA rather than an established RNA thermometer, while **suhB** (sRNA) and **drz-agam-1** (ribozyme) are not established riboswitches.

## Where to find the reproducible outputs

| Output | Directory/file |
| --- | --- |
| Phase-specific sRNA differential-expression tables | `results/sRNA_detection/04_srna_de_filtered/` |
| Selected sRNA IDs and spreadsheets | `results/sRNA_detection/06_intarna/selected_exp_ids.txt`, `selected_sta_ids.txt`, `selected_exp_srnas.xlsx`, `selected_sta_srnas.xlsx` |
| Curated mRNA target tables | `results/sRNA_detection/06_intarna/curated_mrna_tables/` |
| IntaRNA per-target predictions | `results/sRNA_detection/06_intarna/intarna_exp/` and `intarna_sta/` |
| Combined IntaRNA predictions | `results/sRNA_detection/06_intarna/combined_results/` |
| Filtered/prioritised interactions | `results/sRNA_detection/06_intarna/filtered_interactions/` |
| ANNOgesic cis-element source tables | `results/riboswitch_thermometer/exp/` and `sta/` |
| Cis-element–DEG comparisons and mapping audit | `results/riboswitch_thermometer/deg_comparison/` |
| Publication-figure source plots | `results/figures/` |

For Figure 5, `results/figures/cis_elements/figure_5_plot_counts.csv` records the class counts used to generate the RNA thermometer-module and riboswitch-module plots.

## Reproducibility notes

- Differential-expression fold changes are **SMG relative to NG**, and sRNA findings at nominal P < 0.05 are exploratory.
- Sequence extraction and IntaRNA analyses are phase matched (exponential-to-exponential; stationary-to-stationary). Where a gene symbol is duplicated, target FASTA names can include its locus tag (for example, `kfiA__ECOLIN_25915.fasta`) to avoid ambiguity.
- `scripts/06_intarna/03_combine_results.sh` currently generates combined tables directly under `results/sRNA_detection/06_intarna/`. **The committed final combined tables are organised under `combined_results/`.** Move regenerated tables there or adjust the script's output directory before reproducing that final layout. The filtering script accepts either location but checks the root first; avoid leaving stale root-level copies that could override newer `combined_results/` files.
- `scripts/08_cis_deg/compare_cis_elements_to_degs.py` matches ANNOgesic `Associated_CDS` locus tags to curated DEG annotations. It does **not** test genomic coordinate overlap or prove regulation.
- Raw sequencing reads and some reference/annotation inputs must be obtained separately. Review the paths and requirements in **[WORKFLOW.md](WORKFLOW.md)** before running scripts.

## Citation and reuse

If using these analyses, cite the underlying GEO dataset and relevant computational tools. This repository provides computational hypotheses for follow-up experimental investigation rather than validated sRNA functions.
