BASE = /work/calicraw
FASTP = $BASE/Software/fastp
RAW = $BASE/Herbarium_Sequences
TRIMDIR = $BASE/Herbarium_Sequences/Trimmed
CONTGENOME = $BASE/Genome/all_genera_index
CONTDIR = $BASE/Herbarium_Sequences/Contaminated
DECONDIR = $BASE/Herbarium_Sequences/Decontaminated
GENOME = $BASE/Genome/GCA_030408175.1_UTM_Trep_v1.0_genomic.fna
ALIGNDIR = $BASE/Herbarium_Sequences/Aligned
S1 = '_1.fq.gz'
S2 = '_2.fq.gz'

#Index the Genome
bwa index $GENOME

#Create Sequence List
ls $RAW | grep $S1 | sed 's/_1.fq.gz\|_2.fq.gz//' | sort -u > $BASE/Herbarium_list.txt

#Run FastQC
for i in $(cat $BASE/Herbarium_list.txt);
do fastqc $RAW/$i$S1 -o $BASE/Herbarium_Sequences/FastQC;
fastqc $RAW/$i$S2 -o $BASE/Herbarium_Sequences/FastQC;
done

#Trim Herbarium Sequence Reads
for i in $(cat $BASE/Herbarium_list.txt);
do $FASTP -i $RAW/$i$S1 -I $RAW/$i$S2 -o $TRIMDIR/$i$S1 -O 
$TRIMDIR/$i$S2 --trim_poly_x --trim_poly_g --dedup -h 
$TRIMDIR/$i'.html' -g -w 16;
done

#Run Decontamination Protocol
for i in $(cat /work/calicraw/Herbarium_list.txt);
do bowtie2 -p 16 -x $CONTGENOME -1 $TRIMDIR/$i$S1 -2 $TRIMDIR/$i$S2 
--un-conc-gz /work/calicraw/Herbarium_Sequences/Decontaminated/$i'_CLEAN' > $i'FUN_BAC_REMOVED'.sam;
samtools view -hbS -F 4 $i'FUN_BAC_REMOVED.sam' > $CONTDIR/$i'_CONTAMINATED_SEQUENCES.bam';
rm $i'FUN_BAC_REMOVED.sam';
mv $DECONDIR/$i'_CLEAN.1' $DECONDIR/$i'_1.fq.gz';
mv $DECONDIR/$i'_CLEAN.2' $DECONDIR/$i'_2.fq.gz';
done 

#Run Alignment
for FILE in $(cat $BASE/Herbarium_list.txt);
do bwa mem -t 16 $GENOME $DECONDIR/$FILE$S1 $DECONDIR/$FILE$S2 > $ALIGNDIR/$FILE'.sam';
samtools view -bS $ALIGNDIR/$FILE'.sam' | samtools sort -o $ALIGNDIR/$FILE'.bam';
rm $ALIGNDIR/$FILE'.sam';
done

#Mark Duplicates
for i in $(cat /work/calicraw/Herbarium_Sequences/herbseq_list.txt); 
do samtools markdup $ALIGNDIR/$FILE'.bam' $ALIGNDIR/$FILE'marked_.bam'; 
done

#Index Bam Files
for i in $(cat /work/calicraw/Herbarium_Sequences/herbseq_list.txt); 
do samtools index -b $ALIGNDIR/$FILE'marked_.bam'; 
done