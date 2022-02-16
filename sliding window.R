#sliding window 
library(tidyverse)
library(PopGenome)
library(ggplot2)


x = 10519265
y = 10626343


GENOME.class <- readVCF('chr21.five.populations.vcf.gz',frompos = x, topos = y,tid = 'chr21',numcols = 1000, gffpath = 'gff')

#GENOME.class <- readData("t", gffpath = 'gff', format = "VCF")  #read VCF
############################
#TRYING TO ADD GTF FILE TO DO THAT PLOT LETS SEE

#gffpath = 'testgff'

#######################

#NEED TO FIND OUT HOW TO SPLIT THE REF GFF TO ONLY CHR21 AND FEED IN HERE

############

#NEED TO SET POPULATIONS

#set.populations

CEPH <- read.csv('CEPH.csv')
Fin <- read.csv('Finnish.csv')
Br <- read.csv('British.csv')
Ibe <- read.csv('Iberian.csv')
Tos <- read.csv('Toscani.csv')


British <- Br$British
Finnish <- Fin$Finnish
CEPH <- CEPH$CEPH
Iberian <- Ibe$Iberian
Toscani <- Tos$Toscani

GENOME.class <-set.populations(GENOME.class,new.populations = list(c(British),c(Finnish),c(CEPH),c(Iberian),c(Toscani)))
GENOME.class@populations

get.sum.data(object = GENOME.class)


#total SNPS in this = 107079

GENOME.class@n.biallelic.sites + GENOME.class@n.polyallelic.sites

#set populations
GENOME.class <-set.populations(GENOME.class,new.populations = list(c(British),c(AASW)))

GENOME.class@populations

#set up sliding window 

vcf_size <- 107078

window_size <- 10000
window_jump <- 2500


# use seq to find the start points of each window
window_start <- seq(from = 1, to = vcf_size, by = window_jump)
# add the size of the window to each start point
window_stop <- window_start + window_size

#sum of windows starting before the end of vcf_size ( = 0 )
sum(window_start > vcf_size)

#some window stop positions occur past final point ( = 1 )
sum(window_stop > vcf_size)


window_start <- window_start[which(window_stop < vcf_size)]
window_stop <- window_stop[which(window_stop < vcf_size)]

#see how far short the final window is ( = 77 )

vcf_size - window_stop[length(window_stop)]

#JUST BE AWARE THE ANALYSIS DOES NOT INCLUDE ALL VARIANTS

#for now save start/stop positions in data_frame

windows <- data.frame(start = window_start, stop = window_stop, 
                      mid = window_start + (window_stop-window_start)/2)

#WHY TYPE 2 NOT TYPE 1?? ------ 
vcf_windows <- sliding.window.transform(GENOME.class, width = 10000, jump = 2500, type = 2)

# SO IF JUMP = WIDTH THEN AN ERROR OCCURS, THIS IS DUE TO THE FIRST CODE? 

vcf_windows <- diversity.stats(vcf_windows, pi= TRUE)

#Note that here we use mode = "nucleotide" to specify we want it to be calculated sliding averages of nucleotides, 
#rather than using haplotype data, which is the alternative.

vcf_windows <- F_ST.stats(vcf_windows)

#extract nucleotide diversity and correcct for window size 
nd <- vcf_windows@nuc.diversity.within/10000

pops <- c('British','Finnish','CEPH','Iberian','Toscani')

colnames(nd) <- paste0(pops, "_pi")

#extract fst values

fst <- t(vcf_windows@nuc.F_ST.pairwise)

#we use t() to transpose FST matrix so that each column is a pairwise compairson and each row is an estimate for a genome window

#extract dxy 

dxy <- get.diversity(vcf_windows, between = T)[[2]]/10000


# get column names 
x <- colnames(fst)
# replace all occurrences of pop1 with house
x <- sub("pop1", "British", x)
x <- sub("pop2", "Finnish", x)
x <- sub("pop3", "CEPH", x)
x <- sub("pop4", "Iberian", x)
x <- sub("pop5", "Toscani", x)

# does the same thing as above but by indexing the pops vector
#x <- sub("pop1", pops[1], x)
# look at x to confirm the replacement has occurred
x <- sub("/", "_", x)
x
#TO CHANGE ALL POPULATIONS USE CODE BELOW

paste0(x, "_fst")
paste0(x, "_dxy")

colnames(fst) <- paste0(x, "_fst")
colnames(dxy) <- paste0(x, "_dxy")

#need to add 1 more row to each result 

#maybe not needed below 
#a<-matrix(nrow=1,ncol=1,0)
#b<- matrix(nrow=1,ncol=2,0)

dxy <- rbind(dxy,a)
fst <- rbind(fst,a)
nd <- rbind(nd,b)

vcf_data <- as_tibble(data.frame(windows, nd, fst, dxy))


vcf_data %>% select(contains("pi")) %>% summarise_all(mean)


pi_g <- vcf_data %>% select(contains("pi")) %>% gather(key = "species", value = "pi")

pdf('myplot.pdf')
myplot <- ggplot(pi_g, aes(population, pi)) + geom_boxplot() + theme_light() + xlab(NULL)
print(myplot)
dev.off()

pdf('test.pdf')
a <- ggplot(vcf_data, aes(mid/10, British_CEPH_fst)) + geom_line(colour = "red")
a <- a + xlab("Position (Mb)") + ylab(expression(italic(F)[ST]))
print(a)
dev.off()



hs <- vcf_data %>% select(mid, British_pi, Finnish_pi, British_Finnish_fst, British_Finnish_dxy)

hs_g <- gather(hs, -mid, key = "stat", value = "value")

pdf('multiplot.pdf')
a <- ggplot(hs_g, aes(mid/10^6, value, colour = stat)) + geom_line()
a <- a + facet_grid(stat~., scales = "free_y")
a <- a + xlab("Position (Mb)")
print(a)
dev.off()



#TRY TAJIMA

#set the synomymous and non synonymous

GENOME.class <- set.synnonsyn(GENOME.class, ref.chr='chr21-ref.fa')

#number of synonymous changes & non-synonymous (function disrupting) changes 

sum(GENOME.class@region.data@synonymous[[1]]==1, na.rm=TRUE)

#non-synonymous

sum(GENOME.class@region.data@synonymous[[1]]==0, na.rm=TRUE)

#positive selection in this population? 

#na.rm = TRUE because NaN are SNPS in non-coding regions - GENOME.class@region.data@synonymous to see what I mean

#Taj for nonsyn

genes <- splitting.data(GENOME.class, subsites='gene')

genes <- neutrality.stats(genes,subsites='nonsyn',FAST = TRUE)

nonsynTaj <- genes@Tajima.D

#repeat for syn SNP

genes <- splitting.data(GENOME.class, subsites='gene')

genes <- neutrality.stats(genes,subsites='syn',FAST = TRUE)

synTaj <- genes@Tajima.D

#PLOT PLS

#tajimas positive if high variation - e.g balancing selection whereby many alleles can be selected for 
#tajimas negative if low frequency high variation - e.g our lactase gene have few variation in the population


#d ratio? look at the omicron paper - 

#high frequency variation 


plot(nonsynTaj,synTaj,main="2L: Genes: Tajima' D")
points(nonsynTaj,pch=15,col='Orange')
points(synTaj,pch=15,col='Purple')
