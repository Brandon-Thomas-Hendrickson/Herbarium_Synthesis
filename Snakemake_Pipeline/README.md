# Herbarium Snakemake Pipeline

A Snakemake workflow version of the herbarium population genomics pipeline, designed for execution on SLURM HPC clusters.

## Overview

The Snakemake pipeline encodes the same logic as the modular bash pipeline in `../Pipeline/`, but with automatic dependency tracking, parallel job submission, and reproducibility through Snakemake's workflow engine.

```
build_database ──► index_genome
                         │
     ┌───────────────────┘
     ▼
fastqc ──► trim ──► decontaminate ──► align ──► namesort ──► fixmate
     ──► coordsort ──► markdup ──► index_bam ──► downsample
                                                      │
                                               create_beagle
                                                      │
                                    ┌─────────────────┴──────────────────┐
                                    ▼                                    ▼
                              ngsrelate                            get_fourfold
                                    │                                    │
                              filter_related ──────────────────► filter_fourfold_beagle
                                                                         │
                                                           ┌─────────────┴─────────────┐
                                                           ▼                           ▼
                                              ngsadmix (K×seed)                   pcangsd
                                                           │
                                              evaladmix (K×seed)
                                                           │
                                              summarize_likelihood
```

---

## Prerequisites

### 1. Conda environment

```bash
conda env create -f envs/herbarium.yaml
conda activate Herbarium_Pipeline
```

### 2. Snakemake executor plugin for SLURM

```bash
pip install snakemake-executor-plugin-slurm
```

### 3. ANGSD, NGSadmix, ngsRelate, evalAdmix

These must be compiled separately (not available on conda). Set their paths in `config.yaml`:

```bash
# ANGSD
git clone https://github.com/ANGSD/angsd.git
cd angsd && make

# ngsRelate
git clone https://github.com/samtools/htslib
git clone https://github.com/ANGSD/ngsRelate
cd htslib && make -j4
cd ../ngsRelate && make HTSSRC=../htslib/

# evalAdmix
git clone https://github.com/GenisGE/evalAdmix.git
cd evalAdmix && make
```

---

## Configuration

**Edit `config.yaml` before running.** All paths and parameters are defined there.

Key settings:

| Key | Description |
|---|---|
| `base` | Root working directory |
| `genome` | Path to reference genome FASTA |
| `annotation` | Path to reference GFF annotation (for degenotate) |
| `angsd` | Path to compiled ANGSD binary |
| `ngsadmix` | Path to compiled NGSadmix binary |
| `ngsrelate` | Path to compiled ngsRelate binary |
| `evaladmix` | Path to compiled evalAdmix binary |
| `raw_dir` | Directory containing raw FASTQ files (`<sample>_1.fq.gz`) |
| `slurm_account` | HPC allocation account |
| `slurm_partition` | SLURM partition name |
| `k_max` | Maximum K for NGSadmix (default: 9) |
| `seeds` | Random seeds for NGSadmix replicates |

---

## Running the Pipeline

### Dry run (check rule graph without executing)

```bash
snakemake --profile slurm --configfile config.yaml -n
```

### Submit to SLURM

```bash
snakemake \
    --profile slurm \
    --configfile config.yaml \
    --cores all \
    --jobs 100
```

### Run a specific step only

```bash
# Only process alignment steps
snakemake \
    --profile slurm \
    --configfile config.yaml \
    --until markdup \
    --cores all

# Only rerun admixture
snakemake \
    --profile slurm \
    --configfile config.yaml \
    --forcerun ngsadmix \
    --cores all
```

### Visualise the DAG

```bash
snakemake --configfile config.yaml --dag | dot -Tpng > dag.png
```

---

## SLURM Profile

The `slurm/config.yaml` profile sets default SLURM resources for all rules. Individual rules override these via their `resources:` blocks in the Snakefile.

To modify global SLURM defaults:

```yaml
# slurm/config.yaml
default-resources:
  slurm_account: your_account
  slurm_partition: your_partition
  runtime: 120        # minutes
  mem_mb: 16000
```

SLURM log files are written to `logs/slurm/`.

---

## Outputs

| File | Description |
|---|---|
| `{align_dir}/*_downsampled.bam` | Final analysis-ready BAMs |
| `{angsd_dir}/Herbarium.beagle.gz` | Raw Beagle genotype likelihoods |
| `{angsd_dir}/Herbarium_filtered.beagle.gz` | Relatedness-filtered Beagle |
| `{angsd_dir}/Herbarium_4fold.beagle.gz` | 4-fold degenerate site Beagle |
| `{admix_dir}/ADMIX{K}.{seed}.Q` | Ancestry proportion matrix per run |
| `{admix_dir}/evalADMIX{K}.{seed}.corres` | evalAdmix cross-validation residuals |
| `{admix_dir}/K_likelihood.txt` | K / seed / log-likelihood table |
| `{admix_dir}/clumppak_ready.txt` | Two-column (K, logL) for CLUMPPAK/POPHELPER |
| `{pca_dir}/herbarium_pca.cov` | PCAngsd covariance matrix for PCA |

---

## Differences from `../Pipeline/` (Modular Bash)

| Feature | Bash Pipeline | Snakemake Pipeline |
|---|---|---|
| Dependency tracking | Manual (run in order) | Automatic |
| Parallel jobs | Manual | Automatic via SLURM |
| Re-running failed steps | Manual | `--rerun-incomplete` |
| Reproducibility | Moderate | High (checksums, provenance) |
| Cluster submission | Script-level `sbatch` | Per-rule SLURM jobs |

---

## Troubleshooting

**Jobs not submitting:**
Check that `snakemake-executor-plugin-slurm` is installed and the SLURM account/partition in `config.yaml` are correct.

**"Missing files" errors after job completion:**
Increase `latency-wait` in `slurm/config.yaml` (network file systems can be slow).

**Rerunning from a checkpoint:**
Snakemake automatically skips completed outputs. To force a full rerun:
```bash
snakemake --forceall --profile slurm --configfile config.yaml
```
