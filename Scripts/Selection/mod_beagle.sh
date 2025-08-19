#Sum the probabilities of major homozygote and heterozygote for each individual in the beagle file
gunzip -c $BEAGLE | awk '
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
' > Major_allele_probability.txt