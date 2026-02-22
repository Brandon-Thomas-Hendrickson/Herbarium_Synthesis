#!/bin/bash
# =============================================================================
# 03_create_beagle.sh
# Generate Beagle-format genotype likelihoods from aligned BAM files using
# ANGSD.
#
# Input:  $BAM_LIST  – one BAM file path per line (created by step 02)
# Output: $ANGSDDIR/Herbarium.beagle.gz  – used by steps 04 and 05
#
# SLURM directives below; adjust account/partition in config.sh.
# =============================================================================

#SBATCH --job-name create_beagle
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH -t 3-00:00:00
#SBATCH --account ${SLURM_ACCOUNT}
#SBATCH --partition ${SLURM_PARTITION}

set -euo pipefail
source "$(dirname "$0")/config.sh"

mkdir -p "$ANGSDDIR"

echo "[ANGSD] Generating Beagle genotype likelihoods..."

$ANGSD \
    -ref   "$GENOME" \
    -fai   "$GENOME_FAI" \
    -bam   "$BAM_LIST" \
    -SNP_pval    "$SNP_PVAL" \
    -doMajorMinor 4 \
    -minQ        "$MIN_BQ" \
    -minMapQ     "$MIN_MQ" \
    -doPost      1 \
    -GL          1 \
    -doGlf       2 \
    -minMaf      "$MIN_MAF" \
    -doCounts    1 \
    -doMaf       2 \
    -P           8 \
    -out "$ANGSDDIR/Herbarium"

echo "Done. Beagle file: $ANGSDDIR/Herbarium.beagle.gz"
