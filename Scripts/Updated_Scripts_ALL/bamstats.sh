#!/bin/bash

source /home/calicraw/miniconda3/etc/profile.d/conda.sh
conda activate stats_conda

# Define the directory containing the BAM files and the output CSV file
OUTPUT_CSV="/work/calicraw/Projects/HerbariumStructure/aligned/bamstats_report/spain/SPAIN_BAM_STATS.csv"
BAM_DIR="/work/calicraw/Projects/HerbariumStructure/aligned/processed/spain"

# Initialize the CSV file with the header
echo "Sample,Total reads,Mapped reads,Forward strand,Reverse strand,Failed QC,Duplicates,Paired-end reads,Proper-pairs,Both pairs mapped,Read 1,Read 2,Singletons" > $OUTPUT_CSV

cd /work/calicraw/Projects/HerbariumStructure/aligned
ls $BAM_DIR/*.bam > spain_bam_list.txt

# Loop through each BAM file in the directory
for BAM_FILE in $(cat spain_bam_list.txt); do
    # Extract the sample name from the BAM file name
    SAMPLE_NAME=$(basename $BAM_FILE .bam)
    
    # Run bamtools stats and capture the output
    STATS_OUTPUT=$(bamtools stats -in $BAM_FILE)
    
    # Parse the output to extract the desired statistics
    TOTAL_READS=$(echo "$STATS_OUTPUT" | grep "Total reads" | awk '{print $3}')
    MAPPED_READS=$(echo "$STATS_OUTPUT" | grep "Mapped reads" | awk '{print $3}')
    FORWARD_STRAND=$(echo "$STATS_OUTPUT" | grep "Forward strand" | awk '{print $3}')
    REVERSE_STRAND=$(echo "$STATS_OUTPUT" | grep "Reverse strand" | awk '{print $3}')
    FAILED_QC=$(echo "$STATS_OUTPUT" | grep "Failed QC" | awk '{print $3}')
    DUPLICATES=$(echo "$STATS_OUTPUT" | grep "Duplicates" | awk '{print $2}')
    PAIRED_END_READS=$(echo "$STATS_OUTPUT" | grep "Paired-end reads" | awk '{print $3}')
    PROPER_PAIRS=$(echo "$STATS_OUTPUT" | grep "Proper-pairs" | awk '{print $2}')
    BOTH_PAIRS_MAPPED=$(echo "$STATS_OUTPUT" | grep "Both pairs mapped" | awk '{print $4}')
    READ_1=$(echo "$STATS_OUTPUT" | grep "Read 1" | awk '{print $3}')
    READ_2=$(echo "$STATS_OUTPUT" | grep "Read 2" | awk '{print $3}')
    SINGLETONS=$(echo "$STATS_OUTPUT" | grep "Singletons" | awk '{print $2}')
    
    # Write the results to the CSV file
    echo "$SAMPLE_NAME,$TOTAL_READS,$MAPPED_READS,$FORWARD_STRAND,$REVERSE_STRAND,$FAILED_QC,$DUPLICATES,$PAIRED_END_READS,$PROPER_PAIRS,$BOTH_PAIRS_MAPPED,$READ_1,$READ_2,$SINGLETONS" >> $OUTPUT_CSV
done
