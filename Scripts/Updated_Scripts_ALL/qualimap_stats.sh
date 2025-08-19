#!/bin/bash

#!/bin/bash

source /home/calicraw/miniconda3/etc/profile.d/conda.sh
conda activate stats_conda

# Define the directory containing the BAM files and the output CSV file
OUTPUT_CSV="/work/calicraw/Projects/HerbariumStructure/aligned/qualimap_reports/spain/SPAIN_QUALI_STATS.csv"
BAM_DIR="/work/calicraw/Projects/HerbariumStructure/aligned/processed/spain"
QUALIDIR="/work/calicraw/Projects/HerbariumStructure/aligned/qualimap_reports/spain/qualistats_reports"

# Initialize the CSV file with the header
echo "Sample,number of reads,number of mapped reads,number of supplementary alignments,number of secondary alignments,number of duplicated reads,duplication rate,mean mapping quality,GC percentage,general error rate,mean coverageData,std coverageData" > $OUTPUT_CSV

cd /work/calicraw/Projects/HerbariumStructure/aligned
ls $BAM_DIR/*.bam > spain_bam_list.txt

# Loop through each BAM file in the directory
for BAM_FILE in $(cat spain_bam_list.txt); do
    # Extract the sample name from the BAM file name
    SAMPLE_NAME=$(basename $BAM_FILE .bam)
    
    # Run qualimap bamqc and capture the output
    qualimap bamqc -bam $BAM_DIR/$BAM_FILE -outdir $QUALIDIR/$SAMPLE_NAME'_qualimap_output'
    
    # Rename the genome_results.txt file
    mv $QUALIDIR/$SAMPLE_NAME'_qualimap_output/genome_results.txt' $QUALIDIR/$SAMPLE_NAME'_qualimap_stats.txt'   
    
    # Parse the output to extract the desired statistics
    STATS_OUTPUT=$(echo $QUALIDIR/$SAMPLE_NAME'_qualimap_stats.txt')
    TOTAL_READS=$(cat "$STATS_OUTPUT" | grep "number of reads" | awk '{print $5}' | sed 's.,..g')
    MAPPED_READS=$(cat "$STATS_OUTPUT" | grep "number of mapped reads" | awk '{print $6}' | sed 's.,..g')
    SUPPLEMENTARY_ALIGNMENTS=$(cat "$STATS_OUTPUT" | grep "number of supplementary alignments" | awk '{print $6}' | sed 's.,..g')
    SECONDARY_ALIGNMENTS=$(cat "$STATS_OUTPUT" | grep "number of secondary alignments" | awk '{print $6}' | sed 's.,..g')
    DUPLICATED_READS=$(cat "$STATS_OUTPUT" | grep "number of duplicated reads" | awk '{print $7}' | sed 's.,..g')
    DUPLICATION_RATE=$(cat "$STATS_OUTPUT" | grep "duplication rate" | awk '{print $4}' | sed 's.,..g')
    MEAN_MAPPING_QUALITY=$(cat "$STATS_OUTPUT" | grep "mean mapping quality" | awk '{print $5}' | sed 's.,..g')
    GC_PERCENTAGE=$(cat "$STATS_OUTPUT" | grep "GC percentage" | awk '{print $4}' | sed 's.,..g')
    GENERAL_ERROR_RATE=$(cat "$STATS_OUTPUT" | grep "general error rate" | awk '{print $5}' | sed 's.,..g')
    MEAN_COVERAGE=$(cat "$STATS_OUTPUT" | grep "mean coverageData" | awk '{print $4}' | sed 's.,..g')
    STD_COVERAGE=$(cat "$STATS_OUTPUT" | grep "std coverageData" | awk '{print $4}' | sed 's.,..g')
    
    # Write the results to the CSV file
    echo "$SAMPLE_NAME,$TOTAL_READS,$MAPPED_READS,$SUPPLEMENTARY_ALIGNMENTS,$SECONDARY_ALIGNMENTS,$DUPLICATED_READS,$DUPLICATION_RATE,$MEAN_MAPPING_QUALITY,$GC_PERCENTAGE,$GENERAL_ERROR_RATE,$MEAN_COVERAGE,$STD_COVERAGE" >> $OUTPUT_CSV
done
