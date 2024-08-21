ANGSD=/work/calicraw/Software/angsd/angsd

$ANGSD -fai ./Genome/GCA_030408175.1_UTM_Trep_v1.0_genomic.fna.fai -bam herblist_bams.txt -doMaf 1 -SNP_pval 1e-6 -doMajorMinor 4 -minQ 20 -minMapQ 30 -doPost 1 -GL 1 -doGlf 2 -minMaf 0.05 -doCounts 1 -doMaf 2 -P 8 -out ./Herbarium_Sequences/ANGSDout/Herbarium_Cult_Native
