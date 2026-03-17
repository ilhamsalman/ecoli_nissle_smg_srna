
---

## `scripts/01_annogesic/README.md`

```md
# 01_annogesic — ANNOgesic sRNA discovery and annotation

This folder documents the ANNOgesic workflow used to detect transcripts and predict sRNAs, riboswitches, and RNA thermometers.

ANNOgesic was run independently for:
- Exponential phase: `annogesic_exp/`
- Stationary phase: `annogesic_sta/`

## Expected ANNOgesic project layout

For each phase directory (`annogesic_exp` and `annogesic_sta`), the following inputs are required:

### 1) Fragment BAMs
Location:
- `input/BAMs/BAMs_map_reference_genomes/fragment/*.bam`

### 2) Strand-specific WIG coverage tracks (RF)
Location:
- `input/wigs/fragment/*_{plus,minus}.wig`

These WIG files are generated from BAMs in the preprocessing pipeline:
- `scripts/01_annogesic/make_wigs.sh`

### 3) Reference genome + annotation
Locations:
- `input/references/fasta_files/genome.fasta`
- `input/references/annotations/annotation.gff`

### 4) Rfam covariance model and ID lists (for riboswitch/thermometer)
Locations:
- `input/databases/Rfam.cm`
- `input/riboswitch_ID_file/Rfam_riboswitch_ID.txt`
- `input/RNA_thermometer_ID_file/Rfam_RNA_thermometer_ID.txt`

## Running ANNOgesic

The ANNOgesic runner script is written for the exponential phase by default (`PJ="annogesic_exp"`).
To run stationary phase, rerun after editing:
- `PJ="annogesic_exp"` → `PJ="annogesic_sta"`
- Update WIG filenames in the `-fl` lists:
  - `MG_exp_*` → `MG_sta_*`
  - `NG_exp_*` → `NG_sta_*`


