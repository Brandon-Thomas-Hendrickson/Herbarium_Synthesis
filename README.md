# Herbarium_Synthesis

A genomics pipeline for processing and analysing herbarium specimen sequencing data.
Raw paired-end Illumina reads are quality-filtered, decontaminated, aligned to a
reference genome, and then analysed for population structure and admixture.

## Pipeline overview

```
Raw FASTQ
   │
   ▼  01_build_database / build_database
Build contamination DB (NCBI bacteria + fungi, bowtie2 index)
   │
   ▼  02_trim_decontaminate_align / per-sample rules
fastp trim → bowtie2 decontam → bwa align → markdup → coverage downsample
   │
   ▼  03_create_beagle / create_beagle
ANGSD → Beagle genotype likelihoods (SNP_pval < 1e-6, MAF ≥ 0.05)
   │
   ▼  04_ngsadmix / ngsrelate + filter_related + get_fourfold + ngsadmix
ngsRelate filter (r > 0.25) → degenotate 4-fold sites → NGSadmix K1-9
   │
   ▼  05_pcangsd / pcangsd
PCAngsd PCA on 4-fold degenerate site Beagle
```

## Quick start

Choose the version that suits your workflow:

| Version | Use when... | Location |
|---|---|---|
| **Modular bash** | You want to run steps manually, inspect outputs between steps, or submit individual SLURM jobs | [`Pipeline/`](Pipeline/README.md) |
| **Snakemake** | You want automatic dependency tracking, parallel job submission, and reproducibility | [`Snakemake_Pipeline/`](Snakemake_Pipeline/README.md) |

Both versions use the same underlying tools and parameters. See the respective
README files for setup and usage instructions.

## Conda environment

```bash
conda env create -f environment.yaml
conda activate Herbarium_Pipeline
```

This installs: `samtools`, `bwa`, `bowtie2`, `fastp`, `fastqc`, `bamtools`,
`qualimap`, `pcangsd`, `degenotate`, and supporting libraries.

ANGSD, NGSadmix, ngsRelate, and evalAdmix must be compiled separately.
See [`Pipeline/README.md`](Pipeline/README.md) for instructions.

## Repository layout

```
Herbarium_Synthesis/
├── Pipeline/                  # Modular bash pipeline (numbered scripts + config)
├── Snakemake_Pipeline/        # Snakemake workflow for SLURM clusters
├── Scripts/                   # Original / reference scripts
├── environment.yaml           # Conda environment specification
└── README.md
```
