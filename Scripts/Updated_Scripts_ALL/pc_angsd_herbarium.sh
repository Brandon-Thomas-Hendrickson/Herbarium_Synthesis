#!/bin/bash

#SBATCH --account loni_trpopgen03
#SBATCH --partition checkpt
#SBATCH --job-name Herb_PCA
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH -t 3-00:00:00

pcangsd --beagle /work/calicraw/Projects/HerbariumStructure/ANGSD/output/herbarium/beagles/herb_4foldremoved_filtered.beagle.gz --threads 24 --out /work/calicraw/Projects/HerbariumStructure/pcangsd/herb_only_filtered_MAF05 --maf 0.05 --iter 200 > /work/calicraw/Projects/HerbariumStructure/pcangsd/pc_progMAF05.txt 2>&1

