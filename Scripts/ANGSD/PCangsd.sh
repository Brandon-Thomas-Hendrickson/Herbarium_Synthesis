for K in {1..9}; do
  pcangsd --beagle /work/calicraw/Herbarium_Sequences/ANGSDout/Herbarium_Cult_Native.beagle.gz \
  --threads 20 --out /work/calicraw/Herbarium_Sequences/ADMIX/PCangsdHerbADMIX${K} \
  --admix --admix_K ${K} --maf 0.05
done



