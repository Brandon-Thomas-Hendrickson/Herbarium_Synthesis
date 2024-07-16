# Get 4-Four degenerate sites
degenotate.py -a $ANNOTATION -g $GENOME -o degen/ -d " "

#Exract only 4-fold degenerate sites
awk '{if($5==4) print $1"_"$2}' degen/degeneracy-all-sites.bed degen/degen_4fold.bed

#Extract four fold degenerate sites from the beagle file
zcat ANGSDout/Herbarium_Cult_Native.beagle.gz | awk 'NR==FNR {a[$1]; next} $1 in a' markers_fourfolddegen.txt - | gzip > ANGSDout/Herbarium_Cult_Native_FOURFOLDT2.beagle.gz

#Run NGSadmix
for K in {1..9}; do
  for seed in 21 1995 7142023; do
	/work/calicraw/Software/NGSadmix -likes $ANGSDDIR/Herbarium_Cult_Native_FOURFOLDT2.beagle.gz -K $K -P 48 -o $ADMIXDIR/ADMIX${K}.${seed} -minMaf 0.05 -seed $seed
  done
done
