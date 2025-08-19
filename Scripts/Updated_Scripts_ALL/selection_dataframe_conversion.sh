#!/bin/sh

#SBATCH --account loni_trpopgen03
#SBATCH --partition workq
#SBATCH --job-name herb_selection
#SBATCH --nodes=1
#SBATCH -t 3-00:00:00

#Sum the probabilities of major homozygote and heterozygote for each individual in the beagle file
BEAGLEDIR="/work/calicraw/Projects/HerbariumStructure/ANGSD/output/herbarium/beagles"
SELECTIONDIR="/work/calicraw/Projects/HerbariumStructure/ANGSD/output/herbarium/selection"
BEAGLE="combined_file.beagle.gz"

gunzip -c $BEAGLEDIR/$BEAGLE | awk '
NR == 1 {
    # Print the header
    printf "%s,%s,%s", $1, $2, $3;
    for (i = 4; i <= NF; i += 3) {
        ind = $i;  
        printf ",%s", ind;
    }
    print "";
}
NR > 1 {
    # Print the marker and alleles
    printf "%s,%s,%s", $1, $2, $3;
    for (i = 4; i <= NF; i += 3) {
        sum = $i + $(i + 1);
        printf ",%s", sum;
    }
    print "";
}
' > $SELECTIONDIR/Major_allele_probability.txt
