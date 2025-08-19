#!/bin/bash

#SBATCH --account loni_trpopgen03
#SBATCH --partition workq
#SBATCH --job-name Heterozygosity
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH -t 3-00:00:00

# Create a BAM list file
bam_list="/work/calicraw/Projects/HerbariumStructure/file_lists/herbarium_GLUE_spain/bam_list_HERBONLY.txt"
GENOME="/work/calicraw/Genome/GCA_030408175.1_UTM_Trep_v1.0_genomic.fna.fai"
GENOMEREF="/work/calicraw/Genome/GCA_030408175.1_UTM_Trep_v1.0_genomic.fna"
output_dir="/work/calicraw/Projects/HerbariumStructure/heterozygosity/output"

ANGSD="/work/calicraw/Software/angsd/angsd"

# Loop through each BAM file in the list
for bam_file in $(cat $bam_list); do
    sample_name=$(basename "$bam_file" .bam)
    echo $bam_file > temp.txt
    # Run ANGSD for each individual BAM file
    $ANGSD -bam temp.txt -ref $GENOMEREF -fai $GENOME -out "$output_dir/${sample_name}" -anc /work/calicraw/T_ANCESTRAL/T_ancestral.fa -doSaf 1 -GL 1 -doCounts 1 -minQ 20 -minmapq 30 -C 50 -P 8
    /work/calicraw/Software/angsd/misc/realSFS "$output_dir/${sample_name}.saf.idx" -maxiter 2000 -tole 1e-8 > "$output_dir/${sample_name}.ml"
done
