# Herbarium Population Genomics Pipeline

A modular bash pipeline for processing herbarium specimen sequencing data, from raw FASTQ reads through contamination removal, alignment, and population genomics analysis.

## Overview

```
Raw FASTQ Reads
      │
      ▼
01_build_database.sh   ← Build NCBI bacterial/fungal contamination index (run once)
      │
      ▼
02_trim_decontaminate_align.sh   ← QC → trim → decontaminate → align → markdup → downsample
      │
      ├──► 06_bam_qc_summary.sh          ← aggregate bamtools + qualimap CSVs
      │
      ▼
03_create_beagle.sh   ← ANGSD genotype likelihoods → Beagle format
      │
      ├──► 07_hetero_calc.sh             ← per-individual ANGSD doSaf + realSFS
      │
      ▼
04_ngsadmix.sh   ← relatedness filter → 4-fold sites → NGSadmix K1-9 (5 seeds each)
      │
      ▼
05_pcangsd.sh   ← PCAngsd PCA on 4-fold degenerate sites

  08_chloroplast_downsample.sh   ← fixed-fraction chloroplast downsampling [optional]
```

---

## Prerequisites

### 1. Conda environment

```bash
conda env create -f ../environment.yaml
conda activate Herbarium_Pipeline
```

The environment includes: `fastp`, `fastqc`, `bowtie2`, `bwa`, `samtools`, `bamtools`, `qualimap`, `pcangsd`, `degenotate`, and supporting libraries.

### 2. ANGSD and NGSadmix

ANGSD and NGSadmix must be compiled separately (not available on conda):

```bash
git clone https://github.com/ANGSD/angsd.git
cd angsd && make
```

The compiled `angsd` binary and `misc/NGSadmix` binary paths are set in `config.sh`.

### 3. Reference genome

Download the *Trifolium repens* assembly (or your target species):

```
GCA_030408175.1_UTM_Trep_v1.0_genomic.fna
```

Place it under the directory pointed to by `$GENOME` in `config.sh`.

---

## Configuration

**Edit `config.sh` before running any step.** All paths and parameters are centralised there.

Key variables to set:

| Variable | Description |
|---|---|
| `BASE` | Root working directory |
| `GENOME` | Path to reference genome FASTA |
| `ANNOTATION` | Path to reference genome GFF annotation |
| `ANGSD` | Path to compiled ANGSD binary |
| `NGSADMIX` | Path to compiled NGSadmix binary |
| `SLURM_ACCOUNT` | HPC allocation account name |
| `SLURM_PARTITION` | HPC partition (e.g. `workq`, `checkpt`) |

---

## Pipeline Steps

### Step 1 – Build Contamination Database

```bash
bash 01_build_database.sh
```

Downloads all **complete** bacterial and fungal genome assemblies from NCBI GenBank and builds a bowtie2 index. Run **once** before processing any samples.

**Outputs:**
- `$GENOMEDIR/GB_BAC_FUN.fasta` – combined reference FASTA
- `$GENOMEDIR/all_genera_index.*` – bowtie2 index

> This step requires several hours and significant disk space (~50–200 GB).

---

### Step 2 – Trim, Decontaminate, and Align

```bash
bash 02_trim_decontaminate_align.sh
```

Processes all samples in `$RAW` that match the `*_1.fq.gz` / `*_2.fq.gz` naming convention:

1. **FastQC** – quality assessment of raw reads
2. **fastp** – adapter trimming, deduplication, quality filtering (`--cut_right`)
3. **bowtie2** – decontamination against bacterial/fungal index
4. **bwa mem** – alignment to reference genome
5. **samtools fixmate + markdup** – duplicate marking (query-sort → fixmate → coord-sort → markdup)
6. **bamtools / qualimap** – post-alignment QC statistics
7. **Coverage downsampling** – samples with > 4X mean coverage are downsampled to ~3X; samples below 2X are kept as-is

**Outputs per sample (in `$ALIGNDIR`):**
- `<sample>_downsampled.bam` + `.bai` – final BAM for downstream analysis
- `<sample>_bamstats.txt` – bamtools read statistics
- `<sample>_qualimap/` – qualimap coverage report

---

### Step 3 – Create Beagle Genotype Likelihoods

```bash
sbatch 03_create_beagle.sh    # HPC (recommended)
# or
bash 03_create_beagle.sh      # local
```

Runs ANGSD on all downsampled BAMs to produce Beagle-format genotype likelihoods. Filters for:
- SNP p-value < 1e-6
- Minor allele frequency ≥ 0.05
- Mapping quality ≥ 30
- Base quality ≥ 20

**Output:** `$ANGSDDIR/Herbarium.beagle.gz`

---

### Step 4 – Admixture Analysis (NGSadmix)

```bash
sbatch 04_ngsadmix.sh    # HPC (recommended)
# or
bash 04_ngsadmix.sh      # local
```

Full pipeline:

1. **ngsRelate** – compute pairwise kinship coefficients
2. **Relatedness filtering** – individuals with *r* > 0.25 removed (one per pair)
3. **degenotate** – identify 4-fold degenerate (putatively neutral) sites
4. **Beagle subsetting** – retain only 4-fold sites for admixture
5. **NGSadmix** – run K = 1 to `$K_MAX` with 5 random seeds per K
6. **evalAdmix** – cross-validation residuals for each run
7. **Log-likelihood summary** – `K_likelihood.txt` and `clumppak_ready.txt`

**Key outputs (in `$ADMIXDIR`):**
- `ADMIX{K}.{seed}.Q` – ancestry proportion matrices
- `evalADMIX{K}.{seed}` – cross-validation output
- `K_likelihood.txt` – K / seed / log-likelihood table
- `clumppak_ready.txt` – two-column (K, logL) file for CLUMPPAK/POPHELPER

**Selecting best K:** Plot mean log-likelihood vs K; choose the K at the inflection point (elbow). Use evalAdmix residuals to assess fit.

---

### Step 5 – PCA (PCAngsd)

```bash
sbatch 05_pcangsd.sh    # HPC (recommended)
# or
bash 05_pcangsd.sh      # local
```

Runs PCAngsd on the 4-fold filtered Beagle file to produce a covariance matrix for PCA.

**Outputs (in `$PCADIR`):**
- `herbarium_pca.cov` – covariance matrix
- `herbarium_pca.txt` – convergence log

**Visualisation in R:**

```r
cov <- as.matrix(read.table("herbarium_pca.cov"))
eig <- eigen(cov)
plot(eig$vectors[,1], eig$vectors[,2],
     xlab = "PC1", ylab = "PC2",
     pch = 19, main = "Herbarium PCA")
```

---

### Step 6 – BAM QC Summary

```bash
bash 06_bam_qc_summary.sh
```

Aggregates the per-sample `bamtools stats` and qualimap outputs produced by step 02 into two summary CSV files. Re-runs the underlying tools automatically for any sample whose output file is missing.

**Outputs (in `$QCDIR`):**
- `bamstats_summary.csv` – 13 read-level statistics per sample (total reads, mapped, duplicates, strand distribution, etc.)
- `qualimap_summary.csv` – 12 coverage/quality statistics per sample (mean coverage, duplication rate, GC%, error rate, etc.)

---

### Step 7 – Per-Individual Heterozygosity

```bash
sbatch 07_hetero_calc.sh    # HPC (recommended)
# or
bash 07_hetero_calc.sh      # local
```

Runs ANGSD's `doSaf` model on each sample to compute site allele frequency likelihoods, then runs `realSFS` to estimate the maximum-likelihood 1-dimensional SFS. For a single diploid individual the SFS gives a direct estimate of per-site heterozygosity.

**Set `$ANCESTRAL_GENOME` in `config.sh`** before running — the ancestral reference is required for polarising the SFS.

**Outputs per sample (in `$HETDIR`):**
- `<sample>.saf.idx` – site allele frequency index
- `<sample>.ml` – ML SFS; three values: `homref  het  homalt`

To compute heterozygosity from the `.ml` file:

```bash
awk '{print $2 / ($1 + $2 + $3)}' sample.ml
```

---

### Step 8 – Chloroplast Fixed-Fraction Downsampling *(optional)*

```bash
bash 08_chloroplast_downsample.sh
```

Downsamples chloroplast BAM files to user-specified fractions so all samples reach comparable depth for haplotype-network or phylogenetic analyses.

**Create `$CHLORO_DOWNSAMPLE_LIST`** (set path in `config.sh`) as a two-column TSV:

```
sample1.bam   0.5
sample2.bam   0.3
sample3.bam   1.0
```

The fraction is the proportion of reads to keep (e.g. 0.5 = 50%). To calculate fractions, divide your target depth by each sample's observed mean depth.

**Outputs (in `$CHLORO_DOWNSAMPLED_DIR`):**
- `<sample>_downsampled.bam` + `.bai`

---

## Directory Structure After Full Run

```
$BASE/
├── Herbarium_Sequences/
│   ├── FastQC/              quality reports (raw reads)
│   ├── Trimmed/             fastp output
│   ├── Contaminated/        reads mapping to contaminant DB
│   ├── Decontaminated/      clean FASTQ
│   ├── Aligned/             BAM files (*_downsampled.bam)
│   │   └── chloroplast/     chloroplast BAMs + downsampled/ subdirectory
│   ├── ANGSDout/            Beagle + MAF files
│   ├── ADMIX/               NGSadmix Q-matrices and evalAdmix output
│   ├── PCA/                 PCAngsd covariance matrix
│   ├── Heterozygosity/      per-individual .saf.idx and .ml files
│   └── QC_Summary/          bamstats_summary.csv, qualimap_summary.csv
├── Genome/
│   ├── GB_BAC_FUN.fasta     contamination database
│   └── all_genera_index.*   bowtie2 index
└── degen/
    └── degeneracy-all-sites.bed   degenotate output
```

---

## Notes and Tips

- All scripts source `config.sh` at startup. If you need to run a subset of samples, update `$SAMPLE_LIST` and `$BAM_LIST` accordingly.
- The contamination database step (`01_build_database.sh`) only needs to be run once and can be shared across projects.
- SLURM `#SBATCH` directives in scripts 03–05 and 07 use shell variable expansion for the account and partition; these are read at submission time from `config.sh`.
- For very large datasets, consider splitting ANGSD across chromosomes and merging the Beagle files afterwards.
- Steps 06 and 08 run quickly and do not need SLURM; run them directly with `bash`.
- Step 07 requires `$ANCESTRAL_GENOME` — set this in `config.sh` before running.
- Step 08 is independent of the main nuclear genome pipeline and can be run at any time after chloroplast BAMs are available.
