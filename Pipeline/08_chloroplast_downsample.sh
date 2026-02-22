#!/bin/bash
# =============================================================================
# 08_chloroplast_downsample.sh   [OPTIONAL]
# Fixed-fraction downsampling for chloroplast BAM files.
#
# Unlike the adaptive coverage downsampling in step 02 (which targets a
# coverage range), this script uses a user-supplied fraction per sample so
# that all chloroplast samples reach the same relative depth for phylogenetic
# or haplotype-network analyses.
#
# Input:  $CHLORO_DOWNSAMPLE_LIST – two-column TSV: <sample.bam>  <fraction>
#           fraction is a value in (0, 1]; e.g. 0.5 = 50% of reads kept
#         $CHLORO_BAM_DIR         – directory containing the chloroplast BAMs
# Output: $CHLORO_DOWNSAMPLED_DIR/<sample>_downsampled.bam
#
# To create downsampling_chloroplast_list.txt:
#   Compute mean coverage for each chloroplast BAM, decide a target depth,
#   then calculate fraction = target / observed_coverage.
# =============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

mkdir -p "$CHLORO_DOWNSAMPLED_DIR"

if [[ ! -f "$CHLORO_DOWNSAMPLE_LIST" ]]; then
    echo "ERROR: Downsampling list not found: $CHLORO_DOWNSAMPLE_LIST"
    echo "Create a two-column file with <sample.bam>  <fraction> per line."
    exit 1
fi

while IFS=$'\t ' read -r SAMPLE FRACTION; do
    # Skip blank lines and comment lines
    [[ -z "$SAMPLE" || "$SAMPLE" =~ ^# ]] && continue

    SAMPLE_BASENAME=$(basename "$SAMPLE" .bam)
    INPUT_BAM="$CHLORO_BAM_DIR/$SAMPLE"
    OUTPUT_BAM="$CHLORO_DOWNSAMPLED_DIR/${SAMPLE_BASENAME}_downsampled.bam"

    if [[ ! -f "$INPUT_BAM" ]]; then
        echo "  WARNING: $INPUT_BAM not found – skipping."
        continue
    fi

    echo "[downsample] $SAMPLE_BASENAME → fraction $FRACTION"
    samtools view -s "$FRACTION" -b "$INPUT_BAM" -o "$OUTPUT_BAM"
    samtools index "$OUTPUT_BAM"
    echo "  Saved: $OUTPUT_BAM"

done < "$CHLORO_DOWNSAMPLE_LIST"

echo "Done. Downsampled chloroplast BAMs in $CHLORO_DOWNSAMPLED_DIR"
