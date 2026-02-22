#!/bin/bash
# =============================================================================
# 07_hetero_calc.sh
# Per-individual heterozygosity estimation using ANGSD doSaf + realSFS.
#
# For each sample, ANGSD computes the site allele frequency spectrum (SAF)
# and realSFS estimates the maximum-likelihood 1D SFS.  The ratio of
# heterozygous sites to total sites is a per-individual inbreeding/
# heterozygosity proxy.
#
# Input:  $BAM_LIST – one downsampled BAM path per line (from step 02)
# Output: $HETDIR/<sample>.saf.idx   – SAF index
#         $HETDIR/<sample>.ml        – ML SFS (use for heterozygosity)
#
# SLURM directives below; adjust account/partition in config.sh.
# =============================================================================

#SBATCH --job-name hetero_calc
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH -t 3-00:00:00
#SBATCH --account ${SLURM_ACCOUNT}
#SBATCH --partition ${SLURM_PARTITION}

set -euo pipefail
source "$(dirname "$0")/config.sh"

mkdir -p "$HETDIR"

for BAM_FILE in $(cat "$BAM_LIST"); do
    SAMPLE=$(basename "$BAM_FILE" .bam)
    echo "[hetero] Processing $SAMPLE..."

    # Write a single-sample BAM list to a temp file with a unique name
    TMPLIST=$(mktemp "$HETDIR/${SAMPLE}_bam_XXXX.txt")
    echo "$BAM_FILE" > "$TMPLIST"

    # Compute site allele frequency likelihoods (SAF)
    $ANGSD \
        -bam     "$TMPLIST" \
        -ref     "$GENOME" \
        -fai     "$GENOME_FAI" \
        -anc     "$ANCESTRAL_GENOME" \
        -out     "$HETDIR/${SAMPLE}" \
        -doSaf   1 \
        -GL      1 \
        -doCounts 1 \
        -minQ    "$MIN_BQ" \
        -minMapQ "$MIN_MQ" \
        -C       50 \
        -P       8

    rm -f "$TMPLIST"

    # Estimate maximum-likelihood 1D SFS (= heterozygosity when sample n=1)
    $REALSFS \
        "$HETDIR/${SAMPLE}.saf.idx" \
        -maxiter 2000 \
        -tole    1e-8 \
        > "$HETDIR/${SAMPLE}.ml"

    echo "  Done: $HETDIR/${SAMPLE}.ml"
done

echo ""
echo "Heterozygosity estimates written to $HETDIR/*.ml"
echo "Each .ml file contains three space-separated values:"
echo "  <hom_ref_count>  <het_count>  <hom_alt_count>"
echo "  Heterozygosity = het_count / (hom_ref_count + het_count + hom_alt_count)"
