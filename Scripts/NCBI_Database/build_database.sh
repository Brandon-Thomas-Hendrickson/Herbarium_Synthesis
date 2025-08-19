#!/bin/bash
$GENOMEDIR = ./Genome
#Build Contamination Database
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt; 
mv assembly_summary.txt assembly_summary_bac.txt #Bacteria 
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/fungi/assembly_summary.txt; 
mv assembly_summary.txt assembly_summary_fun.txt #Fungi
awk -F '\t' '{if($12=="Complete Genome") print $8,$20}' assembly_summary_bac.txt > assembly_summary_complete_genomes_bac.txt
awk -F '\t' '{if($12=="Complete Genome") print $8,$20}' assembly_summary_fun.txt > assembly_summary_complete_genomes_fun.txt
#Combine Fungi and Bac Texts
cat assembly_summary_complete_genomes_fun.txt >> assembly_summary_complete_genomes_bac.txt
#Filter by Genera
sort -u -k1,1 -t' ' assembly_summary_complete_genomes_bac.txt > assembly_combined.txt
awk '{print gsub(/ /,"")}'  assembly_combined.txt > tab_per_line.txt
paste tab_per_line.txt assembly_combined.txt > ass_tab.txt
for i in `seq 2 5`; do awk '{if($1=$i) print $2,$3,$($i+2)}' ass_tab.txt >> gen_ftp.txt ; done
awk '{print $3}' gen_ftp.txt >> ftp_list.txt
mkdir BAC_FUN_GENOMES
for i in $(cat ftp_list.txt); do wget $i/${i#*/*/*/*/*/*/*/*/*/}'_genomic.fna.gz' -P BAC_FUN_GENOMES;done
cd BAC_FUN_GENOMES
gunzip ./*.gz
cat ./*.fna > $GENOMEDIR/GB_BAC_FUN.fasta
 
#Build Genome Index of Bacterial and Fungi Genome Database 
###Must use Bowtie2 version 2.2.4
mkdir ./Software/bowtie2
cd bowtie2
wget 
https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.4/bowtie2-2.2.4-linux-x86_64.zip/download 
--no-check-certificate
unzip download
./Software/bowtie2/bowtie2-2.2.4-linux-x86_64/bowtie2-2.2.4/bowtie2-build --packed $GENOMEDIR/GB_BAC_FUN.fasta $GENOMEDIR/all_genera_index


