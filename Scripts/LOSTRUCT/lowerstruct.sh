# trepens aus lowerstruct
cd ~/work/calicraw

# install PCAngsd v1.2 locally
module load python
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd/
pip install --user -r requirements.txt
python setup.py build_ext --inplace
pip3 install -e .

# make bam list
cat samples-meta-443.txt | awk -F "\t" '{print $1}' \
| while read samp; do ls bam/$samp.trepens2021.realigned.bam; done > bamlist-443.txt

# make window list
# in R
options(scipen=999)

# read in a fasta index and make a list of 10Mb windows across the genome
dat = read.table("~/work/calicraw/LOSTRUCT/T_repens_Genome.fa.fai")[, 1:2]

# only keep chromosomes
dat = dat[dat$V2 > 10000000, ]

# a vector of all contigs
ctgs = dat$V1

for (c in 1:length(ctgs)){
clen = dat[dat$V1 == ctgs[c], 2]
dat.c = data.frame(ctgs[c], seq(1, clen, 10000000), seq(1, clen, 10000000) + 9999999)
dat.c[nrow(dat.c), 3] = clen
if (c == 1){ dat.out = dat.c } else { dat.out = rbind(dat.out, dat.c) }
}

write.table(dat.out, file = "~/work/calicraw/windows-chrs-only-10Mb.list", col.names = F, row.names = F, quote = F)
# 113 windows
###

##########
# a-PCA-admix.sh
##########
#!/bin/bash
#SBATCH --job-name=a-PCA-admix
#SBATCH --account=def-rieseber
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=03:00:00
#SBATCH --output=a-PCA-admix-%a.out
#SBATCH --error=a-PCA-admix-%a.err
#SBATCH --array=1-113

cd ~/work/calicraw

N=$SLURM_ARRAY_TASK_ID
WIND=$(cat windows-chrs-only-10Mb.list | head -n $N | tail -n 1 | awk '{print $1 ":" $2 "-" $3}')

module load nixpkgs/16.09
module load intel/2018.3
module load angsd/0.929
module load plink/1.9b_5.2-x86_64

angsd \
-nThreads 8 \
-bam bamlist-443.txt \
-r $WIND \
-out angsd/PCA-admix-$N \
-GL 2 \
-doMajorMinor 1 \
-doCounts 1 \
-doGLF 2 \
-SNP_pval 1e-6 \
-doMaf 2 \
-doGeno -1 \
-doPost 1 \
-minMapQ 30 \
-minQ 20 \
-trim 5 \
-minMaf 0.05 \
-minInd 333 \
-geno_minDepth 2 \
-setMinDepthInd 2 \
-uniqueOnly 1

##########

# an R script to divide up 10Mb windows into 100kb windows for lostruct
# minimum number of sites per window = 100

##########
# lostruct-100kb-winds.R
##########
options(scipen=999)

library(data.table)

setwd("/work/calicraw/LOSTRUCT")

N = commandArgs(trailingOnly = TRUE)[1]

# read in beagle file
beagle = fread(paste0("angsd/PCA-admix-", N, ".beagle.gz"), header = T)

# read in mafs file
mafs = fread(paste0("angsd/PCA-admix-", N, ".mafs.gz"), header = T)

# get the chromosome name
chr = mafs$chromo[1]

# get 100kb windows for 10Mb window
wmin = floor(min(mafs$position) / 100000) * 100000
wmax = ceiling(max(mafs$position) / 100000) * 100000

winds = cbind(seq(wmin, wmax - 100000, 100000) + 1,
	seq(wmin + 100000, wmax, 100000))

# subset beagle for each 100kb window and write out
for (i in 1:nrow(winds)){
beagle100k = beagle[mafs$position >= winds[i, 1] & mafs$position <= winds[i, 2],]

if (nrow(beagle100k >= 100)){
fwrite(beagle100k, file = paste0("lostruct/", chr, "_", winds[i, 1] + 49999, ".beagle.gz"),
	col.names = F, row.names = F, quote = F, sep = "\t", compress = "gzip")
}

}

##########

# run R script to make 100kb window beagle files
module load r/4.3.1

for N in {1..113}; do echo $N; Rscript lostruct-100kb-winds.R $N; done

###

# run pcangsd on each 100kb window
ls lostruct/*.beagle.gz | while read beag
do python pcangsd/pcangsd/pcangsd.py -t 8 -b $beag -o ${beag/.beagle.gz/}
rm ${beag/.beagle.gz/.args}
done

# an R script to take a covariance matrix for a window and calculate a row for a lostruct input file

##########
# eigenstuff.R
##########
file = commandArgs(trailingOnly = TRUE)[1]

COV = read.table(file)
COV.pca = eigen(COV)
eigenstuff = c(sum(COV*COV), COV.pca$values[1], COV.pca$values[2], COV.pca$vectors[, 1], COV.pca$vectors[, 2])

write.table(eigenstuff,
file = paste0(file, ".es"),
col.names = F,
row.names = F,
quote = F)

##########

# run eigenstuff on each 100kb window
module load r/4.3.1

cd lostruct/

ls *.cov | while read mat
do Rscript ../eigenstuff.R $mat
done

# combine windows for each chromosome
# and output a list of window locations
cd ~/work/calicraw

ls lostruct/*.es | awk -F "/" '{print $2}' | awk -F "_" '{print $1}' | sort -u \
| while read CHR
do ls lostruct/${CHR}_*.cov.es | awk -F "/" '{print $2}' | awk -F "_" '{print $2}' \
| awk -F "." '{print $1}' | sort -n > lostruct2/$CHR.locs
ls lostruct/${CHR}_*.cov.es | awk -F "/" '{print $2}' | awk -F "_" '{print $2}' \
| awk -F "." '{print $1}' | sort -n \
| while read LOC
do echo $(cat lostruct/${CHR}_${LOC}.cov.es | tr -s '\n' ' ')
done > lostruct2/$CHR.lostruct
done


# in R
library(lostruct)
library(ggplot2)
library(cowplot)

setwd("~/work/calicraw/lostruct2/")

CHRS = read.table("~/work/calicraw/ref/chromo-length.txt", header = F)

#hbs = as.data.frame(read.table("~/work/calicraw/haploblocks-6-1-23.txt", header = F))

for (n in 1:length(CHRS$V1)){

# List all files with the .lostruct extension
files <- list.files(pattern = "\\.lostruct$")

# Read each file into a matrix and store in a list
eigenstuff_list <- lapply(files, function(file) {
  as.matrix(read.table(file, header = FALSE))
})

# If you need to combine these matrices into a single matrix, you can use do.call and rbind
eigenstuff <- do.call(rbind, eigenstuff_list)

# and now we're back to lostruct code
windist = pc_dist(eigenstuff, npc = 2)

# 5D MDS
fit5d = cmdscale(windist, eig = TRUE, k = 5)
mds.coords = fit5d$points

outdf = as.data.frame(mds.coords)
colnames(outdf) = c("MDS1", "MDS2", "MDS3", "MDS4", "MDS5")

# List all files with the .locs extension
locs_files <- list.files(pattern = "\\.locs$")

# Initialize an empty data frame
outdf <- data.frame(WIND = numeric(), CHR = character(), stringsAsFactors = FALSE)

# Read each file and append to the data frame
for (file in locs_files) {
  wind_data <- scan(file)
  chr_name <- sub("\\.locs$", "", file)
  temp_df <- data.frame(WIND = wind_data, CHR = chr_name, stringsAsFactors = FALSE)
  outdf <- rbind(outdf, temp_df)
}

# Print the resulting data frame
print(outdf)

outdf = outdf[, c("CHR", "WIND", "MDS1", "MDS2", "MDS3", "MDS4", "MDS5")]

locs = outdf$WIND

#hbc = hbs[hbs$V1 == CHRS$V1[n], ]

# for pair of MDS i and j
for (i in 1:(ncol(mds.coords) - 1)){
for (j in (i + 1):ncol(mds.coords)){

# identify 2D outliers for MDS pair
mds.corners = corners(mds.coords[, c(i, j)], prop = .05)

corners.mds.i.j = rep(NA, nrow(mds.coords))

# add corner annotations
for (c in 1:ncol(mds.corners)){
corners.mds.i.j[mds.corners[, c]] = c
}

# name the vector with MDS i and j
outdf = cbind(outdf, corners.mds.i.j)
colnames(outdf)[ncol(outdf)] = paste0("corners.mds.", i, ".", j)

mds.i.c.i.j.plot = ggplot() +
	geom_point(aes(x = locs, y = mds.coords[, i]), size = 2, shape = 19) +
	geom_point(aes(x = locs[mds.corners[, 1]], y = mds.coords[mds.corners[, 1], i]), size = 2, shape = 19, color = "#1B9E77") +
	geom_point(aes(x = locs[mds.corners[, 2]], y = mds.coords[mds.corners[, 2], i]), size = 2, shape = 19, color = "#D95F02") +
	geom_point(aes(x = locs[mds.corners[, 3]], y = mds.coords[mds.corners[, 3], i]), size = 2, shape = 19, color = "#7570B3") +
#	geom_segment(aes(x = hbc[,2], y = 1, xend = hbc[,3], yend = 1), colour = "#a6cee3", size = 4) +
	ylim(-1, 1) +
	xlim(0, CHRS$V2[n]) +
	ggtitle("") +
	labs(x = "",
		y = paste0("MDS", i)) +
	theme_classic() +
	theme(legend.position = "none",
		plot.title = element_text(size = 28),
		axis.text = element_text(size = 24),
		axis.title = element_text(size = 28)
	)

mds.j.c.i.j.plot = ggplot() +
	geom_point(aes(x = locs, y = mds.coords[, j]), size = 2, shape = 19) +
	geom_point(aes(x = locs[mds.corners[, 1]], y = mds.coords[mds.corners[, 1], j]), size = 2, shape = 19, color = "#1B9E77") +
	geom_point(aes(x = locs[mds.corners[, 2]], y = mds.coords[mds.corners[, 2], j]), size = 2, shape = 19, color = "#D95F02") +
	geom_point(aes(x = locs[mds.corners[, 3]], y = mds.coords[mds.corners[, 3], j]), size = 2, shape = 19, color = "#7570B3") +
#	geom_segment(aes(x = hbc[,2], y = 1, xend = hbc[,3], yend = 1), colour = "#a6cee3", size = 4) +
	ylim(-1, 1) +
	xlim(0, CHRS$V2[n]) +
	ggtitle("") +
	labs(x = "",
		y = paste0("MDS", j)) +
	theme_classic() +
	theme(legend.position = "none",
		plot.title = element_text(size = 28),
		axis.text = element_text(size = 24),
		axis.title = element_text(size = 28)
	)

ggsave(plot_grid(mds.i.c.i.j.plot, mds.j.c.i.j.plot, ncol = 1, align = "hv"),
	file = paste0("~/work/calicraw/figs/", CHRS$V1[n], "-", i, "-", j, "-mdsmans.png"),
	device = "png",
	width = 1920,
	height = 240 * 2,
	units = "px",
	dpi = 72)

}
}

write.table(outdf, file = paste0(CHRS$V1[n], ".mds"), col.names = T, row.names = F, quote = F)

}
