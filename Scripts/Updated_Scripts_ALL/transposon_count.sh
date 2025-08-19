#!/bin/sh

#SBATCH --account loni_trpopgen03
#SBATCH --partition workq
#SBATCH --job-name transposons
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 3-00:00:00

sed 's,/work/calicraw/Projects/HerbariumStructure/aligned/processed/herbarium/,,' /work/calicraw/Projects/HerbariumStructure/file_lists/herbarium/bam_list_filtered_markedadded.txt > /work/calicraw/Projects/HerbariumStructure/file_lists/herbarium/bam_list_filtered_markedadded_TEMP.txt

source /home/calicraw/miniconda3/etc/profile.d/conda.sh
conda activate alignment_conda

#Extract bam and TE annotations
for i in $(cat /work/calicraw/Projects/HerbariumStructure/file_lists/herbarium/bam_list_filtered_markedadded_TEMP.txt); do
    bedtools intersect -bed -abam /work/calicraw/Projects/HerbariumStructure/aligned/processed/herbarium/$i -b /work/calicraw/Genome/UTM_Trep_v1.0_haploid_reference_repeatMasker_CP.bed -wa -wb > /work/calicraw/Projects/HerbariumStructure/transposons/output/${i%_marked.bam}'.txt';
    awk '{print $14, $15, $17}' /work/calicraw/Projects/HerbariumStructure/transposons/output/${i%_marked.bam}'.txt' | sort | uniq -c | sort -k1,1nr -k2,2 -k3,3 > /work/calicraw/Projects/HerbariumStructure/transposons/output/${i%_marked.bam}'_COUNTS.txt';
    awk '{print $4}' /work/calicraw/Projects/HerbariumStructure/transposons/output/${i%_marked.bam}'_COUNTS.txt' | sort | uniq -c | sort -nr > /work/calicraw/Projects/HerbariumStructure/transposons/output/${i%_marked.bam}'_MOTIFS.txt';
    rm /work/calicraw/Projects/HerbariumStructure/transposons/output/${i%_marked.bam}'.txt';
done

