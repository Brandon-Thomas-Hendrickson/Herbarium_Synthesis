#!/bin/sh

#SBATCH --account loni_trpopgen03
#SBATCH --partition workq
#SBATCH --job-name G1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH -t 3-00:00:00

ANGSDDIR="/work/calicraw/Projects/HerbariumStructure/ANGSD/output/herbarium/beagles"
ANGSD="/work/calicraw/Software/angsd/angsd"
GENOME="/work/calicraw/Genome/GCA_030408175.1_UTM_Trep_v1.0_genomic.fna.fai"
BAMLIST="/work/calicraw/Projects/HerbariumStructure/file_lists/herbarium/G1_bam_list.txt"
BAMDIR="/work/calicraw/Projects/HerbariumStructure/aligned/processed/herbarium/nuclear"
GENOMEREF="/work/calicraw/Genome/GCA_030408175.1_UTM_Trep_v1.0_genomic.fna"

$ANGSD -ref $GENOMEREF -fai $GENOME -bam $BAMLIST -SNP_pval 1e-6 -doMajorMinor 4 -minQ 20 -minMapQ 30 -doPost 1 -GL 1 -doGlf 2 -minMaf 0.05 -doCounts 1 -doMaf 2 -P 8 -out $ANGSDDIR/G1_herb
