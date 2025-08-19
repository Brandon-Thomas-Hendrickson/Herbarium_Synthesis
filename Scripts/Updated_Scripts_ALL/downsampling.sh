#!/bin/bash

FILESDIR="/work/calicraw/Projects/HerbariumStructure/file_lists"
# File containing sample names and downsampling percentages
DOWNSAMPLE_FILE="downsampling_chloroplast_list.txt"
BAMDIR="/work/calicraw/Projects/HerbariumStructure/aligned/processed/herbarium/chloroplast"
# Directory for the downsampled BAM files
DOWNSAMPLED_DIR="/work/calicraw/Projects/HerbariumStructure/aligned/processed/herbarium/chloroplast/downsampled"

# Loop through each line in the file
while read -r SAMPLE PERCENT; do
    # Extract the sample name without the directory or extension
    BASENAME=$(basename "$SAMPLE" .bam)

    # Run /work/calicraw/Software/samtools-1.21/samtools to downsample the BAM file
    /work/calicraw/Software/samtools-1.21/samtools view -s $PERCENT -b $BAMDIR/$SAMPLE > "$DOWNSAMPLED_DIR/${BASENAME}_downsampled.bam"

    echo "Downsampled $SAMPLE to $PERCENT and saved as ${BASENAME}_downsampled.bam"
done < $FILESDIR/$DOWNSAMPLE_FILE
