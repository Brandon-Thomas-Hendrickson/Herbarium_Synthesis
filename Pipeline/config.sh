#!/bin/bash
# =============================================================================
# config.sh  –  Herbarium Pipeline Configuration
# Edit this file before running any pipeline step.
# =============================================================================

# ---------------------------------------------------------------------------
# Base directory (all other paths are derived from this unless overridden)
# ---------------------------------------------------------------------------
BASE=/work/calicraw

# ---------------------------------------------------------------------------
# Reference genome
# ---------------------------------------------------------------------------
GENOME=$BASE/Genome/GCA_030408175.1_UTM_Trep_v1.0_genomic.fna
GENOME_FAI=${GENOME}.fai
ANNOTATION=$BASE/Genome/GCA_030408175.1_UTM_Trep_v1.0_genomic.gff

# ---------------------------------------------------------------------------
# Contamination database (output of 01_build_database.sh)
# ---------------------------------------------------------------------------
GENOMEDIR=$BASE/Genome
CONTGENOME=$GENOMEDIR/all_genera_index   # bowtie2 index prefix

# ---------------------------------------------------------------------------
# Ancestral genome (used by hetero_calc for doSaf polarisation)
# ---------------------------------------------------------------------------
ANCESTRAL_GENOME=$BASE/Genome/T_ancestral.fa

# ---------------------------------------------------------------------------
# Software paths
# ---------------------------------------------------------------------------
ANGSD=$BASE/Software/angsd/angsd
REALSFS=$BASE/Software/angsd/misc/realSFS
NGSADMIX=$BASE/Software/angsd/misc/NGSadmix
NGSRELATE=$BASE/Software/ngsRelate/ngsRelate
EVALADMIX=$BASE/Software/evalAdmix/evalAdmix
DEGENOTATE=degenotate.py                # assumes on PATH via conda

# ---------------------------------------------------------------------------
# Input / output directories
# ---------------------------------------------------------------------------
RAW=$BASE/Herbarium_Sequences              # raw paired FASTQ (*_1.fq.gz, *_2.fq.gz)
FASTQCDIR=$BASE/Herbarium_Sequences/FastQC
TRIMDIR=$BASE/Herbarium_Sequences/Trimmed
CONTDIR=$BASE/Herbarium_Sequences/Contaminated
DECONDIR=$BASE/Herbarium_Sequences/Decontaminated
ALIGNDIR=$BASE/Herbarium_Sequences/Aligned
ANGSDDIR=$BASE/Herbarium_Sequences/ANGSDout
ADMIXDIR=$BASE/Herbarium_Sequences/ADMIX
PCADIR=$BASE/Herbarium_Sequences/PCA
DEGENDIR=$BASE/degen
HETDIR=$BASE/Herbarium_Sequences/Heterozygosity  # per-individual heterozygosity output
QCDIR=$BASE/Herbarium_Sequences/QC_Summary       # aggregated QC CSV files

# ---------------------------------------------------------------------------
# Chloroplast downsampling (step 08 – optional)
# ---------------------------------------------------------------------------
CHLORO_BAM_DIR=$BASE/Herbarium_Sequences/Aligned/chloroplast
CHLORO_DOWNSAMPLE_LIST=$BASE/downsampling_chloroplast_list.txt  # <sample.bam> <fraction>
CHLORO_DOWNSAMPLED_DIR=$BASE/Herbarium_Sequences/Aligned/chloroplast/downsampled

# ---------------------------------------------------------------------------
# Sample list (created by 02_trim_decontaminate_align.sh; re-used by later steps)
# ---------------------------------------------------------------------------
SAMPLE_LIST=$BASE/Herbarium_list.txt
BAM_LIST=$BASE/bam_list.txt

# ---------------------------------------------------------------------------
# Read suffix convention
# ---------------------------------------------------------------------------
S1=_1.fq.gz
S2=_2.fq.gz

# ---------------------------------------------------------------------------
# Thread counts
# ---------------------------------------------------------------------------
THREADS_ALIGN=16     # fastp, bowtie2, bwa mem
THREADS_SORT=48      # samtools sort / fixmate / markdup
THREADS_ADMIX=48     # NGSadmix, evalAdmix, ngsRelate
THREADS_PCA=24       # PCAngsd

# ---------------------------------------------------------------------------
# ANGSD SNP / quality filters
# ---------------------------------------------------------------------------
SNP_PVAL=1e-6
MIN_MAF=0.05
MIN_MQ=30
MIN_BQ=20

# ---------------------------------------------------------------------------
# NGSadmix parameters
# ---------------------------------------------------------------------------
K_MAX=9
SEEDS=(21 1995 7142023 3169147 1964)

# ---------------------------------------------------------------------------
# SLURM account / partition (used by scripts with #SBATCH directives)
# ---------------------------------------------------------------------------
SLURM_ACCOUNT=loni_trpopgen03
SLURM_PARTITION=checkpt
