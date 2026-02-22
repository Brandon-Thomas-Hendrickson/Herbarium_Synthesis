#!/bin/bash
# =============================================================================
# 05_pcangsd.sh
# Principal Component Analysis using PCAngsd.
#
# Input:  $ANGSDDIR/Herbarium_4fold.beagle.gz  (from step 04)
# Output: $PCADIR/herbarium_pca.cov   – covariance matrix (use in R/Python)
#         $PCADIR/herbarium_pca.txt   – progress log
#
# SLURM directives below; adjust account/partition in config.sh.
# =============================================================================

#SBATCH --job-name pcangsd
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH -t 3-00:00:00
#SBATCH --account ${SLURM_ACCOUNT}
#SBATCH --partition ${SLURM_PARTITION}

set -euo pipefail
source "$(dirname "$0")/config.sh"

mkdir -p "$PCADIR"

BEAGLE_4FOLD="$ANGSDDIR/Herbarium_4fold.beagle.gz"

echo "[PCAngsd] Running PCA..."

pcangsd \
    --beagle  "$BEAGLE_4FOLD" \
    --threads "$THREADS_PCA" \
    --out     "$PCADIR/herbarium_pca" \
    --maf     "$MIN_MAF" \
    --iter    200 \
    > "$PCADIR/herbarium_pca.txt" 2>&1

echo "Done."
echo "  Covariance matrix : $PCADIR/herbarium_pca.cov"
echo "  Log               : $PCADIR/herbarium_pca.txt"
echo ""
echo "To visualise in R:"
echo "  cov <- read.table('$PCADIR/herbarium_pca.cov')"
echo "  eig <- eigen(cov)"
echo "  plot(eig\$vectors[,1], eig\$vectors[,2])"
