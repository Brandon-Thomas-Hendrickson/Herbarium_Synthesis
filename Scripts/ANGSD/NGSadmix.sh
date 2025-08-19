#Run NGSrelate for filtering out clones or closely related individuals
cd /work/calicraw/Software/
git clone --recursive https://github.com/SAMtools/htslib
git clone https://github.com/ANGSD/ngsRelate
cd htslib/;make -j2;cd ../ngsRelate;make HTSSRC=../htslib/

cd /work/calicraw 

#Run NGSrelate
zcat ./Herbarium_Sequences/ANGSDout/Herbarium_Cult_Native.mafs.gz | cut -f5 |sed 1d >./Herbarium_Sequences/ANGSDout/freq
./ngsrelate -g ./Herbarium_Sequences/ANGSDout/Herbarium_Cult_Native.beagle.gz -n 100 -f ./Herbarium_Sequences/ANGSDout/freq -O newres -m 1 -p 48 -i 5000

python remove_related_individuals.py

# Filter the beagle file to exclude individuals to remove
zcat $ANGSDDIR/Herbarium_Cult_Native.beagle.gz | awk 'NR==FNR {a[$1]; next} !($1 in a)' individuals_to_remove.txt - | gzip > $ANGSDDIR/Herbarium_Cult_Native_Filtered.beagle.gz

# Get 4-Four degenerate sites
degenotate.py -a $ANNOTATION -g $GENOME -o degen/ -d " "

#Exract only 4-fold degenerate sites
awk '{if($5==4) print $1"_"$2}' degen/degeneracy-all-sites.bed degen/degen_4fold.bed

#Extract four fold degenerate sites from the beagle file
zcat $ANGSDDIR/Herbarium_Cult_Native_Filtered.beagle.gz | awk 'NR==FNR {a[$1]; next} $1 in a' degen/degen_4fold.bed - | gzip > $ANGSDDIR/Herbarium_Cult_Native_FOURFOLDT2.beagle.gz

#Run NGSadmix
for K in {1..9}; do
  for seed in 21 1995 7142023 3169147 1964; do
	/work/calicraw/Software/NGSadmix -likes $ANGSDDIR/Herbarium_Cult_Native_FOURFOLDT2.beagle.gz -K $K -P 48 -o /work/calicraw/Herbarium_Sequences/ADMIX/ADMIX${K}.${seed} -minMaf 0.05 -seed $seed
  done
done

#Run evaladmix on NGSadmix output for cross validation
cd /work/calicraw/Software/
git clone https://github.com/GenisGE/evalAdmix.git
cd evalAdmix
make

cd /work/calicraw/Herbarium_Sequences/ADMIX

for K in {1..9}; do
  for seed in 21 1995 7142023 3169147 1964; do
  /work/calicraw/Software/evalAdmix/evalAdmix -q ADMIX${K}.${seed}.Q -P 48 -o evalADMIX${K}.${seed}
  done
done

#From the log output of NGSadmix, extract the K and the log likelihood
for K in {1..9}; do
  for seed in 21 1995 7142023 3169147 1964; do
    grep -A 1 "best like=" ADMIX${K}.${seed}.log | awk -v K=$K -v seed=$seed '{print K, seed, $2}'
  done
done > K_likelihood.txt

sed -i 's/like=//g' K_likelihood.txt

awk '{print $1,$3}' K_likelihood.txt > clumppak_ready.txt