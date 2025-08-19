#!/bin/bash

#SBATCH --account loni_trpopgen03
#SBATCH --partition checkpt
#SBATCH --job-name Herb_PCA
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH -t 3-00:00:00


ADMIXDIR="/work/calicraw/Projects/HerbariumStructure/ngs_admix"
ANGSDDIR="/work/calicraw/Projects/HerbariumStructure/ANGSD/output/herbarium/beagles"
#Run NGSadmix
for K in {1..9}; do
  for seed in 21 1995 7142023 3169147 1964; do
        /work/calicraw/Software/angsd/misc/NGSadmix -likes $ANGSDDIR/herb_4foldremoved_filtered.beagle.gz -K $K -P 48 -o $ADMIXDIR/ADMIX${K}.${seed} -minMaf 0.05 -seed $seed 
  done
done

#From the log output of NGSadmix, extract the K and the log likelihood
for K in {1..9}; do
  for seed in 21 1995 7142023 3169147 1964; do
    grep -A 1 "best like=" $ADMIXDIR/ADMIX${K}.${seed}.log | awk -v K=$K -v seed=$seed '{print K, seed, $2}'  
  done
done > $ADMIXDIR/K_likelihood.txt

sed -i 's/like=//g' $ADMIXDIR/K_likelihood.txt

awk '{print $1,$3}' $ADMIXDIR/K_likelihood.txt > $ADMIXDIR/clumppak_ready.txt
