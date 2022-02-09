library(PopGenome)

#TEST1 FOLDER CONTAINS UNCOMPRESSED 2 POPULATION VCF

GENOME.class <- readData("test1", gffpath = 'testgff', format = "VCF")  #read VCF
############################
#TRYING TO ADD GTF FILE TO DO THAT PLOT LETS SEE

#gffpath = 'testgff'

#######################

#NEED TO FIND OUT HOW TO SPLIT THE REF GFF TO ONLY CHR21 AND FEED IN HERE

############

#NEED TO SET POPULATIONS

#set.populations

Br <- read.csv('British.csv')
Aa <- read.csv('AASW.csv')

British <- Br$British
AASW <- Aa$AASW

GENOME.class <-set.populations(GENOME.class,new.populations = list(c(British),c(AASW)))
GENOME.class@populations

#Thoughts - load in full unfiltered VCF here, then do the population filtering in here opposed to BCFTools, then we can get the nucleotide diversity stats? 
#Trying readdata() to load full VCF, keep getting errors? Investigate
#Trying to use phased data, this can be fed in https://www.rdocumentation.org/packages/PopGenome/versions/2.7.5/topics/readVCF
#Can use unphased data apparently, trying to download on bash - check progress

#NEUTRALITY STATS

GENOME.class <- neutrality.stats(GENOME.class) #run stats on the VCF

get.neutrality(GENOME.class)[[1]] #pull out stats from matrix

#                                Tajima.D n.segregating.sites Rozas.R_2   Fu.Li.F   Fu.Li.D Fu.F_S Fay.Wu.H Zeng.E Strobeck.S
#chr21.recalibrated.samples.vcf 0.3370584               69527        NA 0.4104502 0.3066595     NA      NaN    NaN         NA

#DIVERSITY STATS

#run this first to get FST stats which diversity stats rely on? 
#make a F_ST for each pair, will need to run the below on a genome.class subset

GENOME.class <- F_ST.stats(GENOME.class)

get.diversity(GENOME.class)[1]

#[[1]]
#nuc.diversity.within hap.diversity.within       Pi hap.F_ST.vs.all nuc.F_ST.vs.all
#chr21.recalibrated.samples.vcf             16948.74                    1 16948.74               0       0.1175796


get.diversity(GENOME.class)[2]

#[[1]]
#nuc.diversity.within hap.diversity.within       Pi hap.F_ST.vs.all nuc.F_ST.vs.all
#chr21.recalibrated.samples.vcf             18526.77                    1 18526.77               0       0.1175796

#TRY PLOT

#set the synomymous and non synonymous

GENOME.class <- set.synnonsyn(GENOME.class,ref.chr='chr21-ref.fa')

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

#add legend


#add legend to colour non-syn & syn

#####

NUCLEOTIDEDIVERSITY.class <- diversity.stats(GENOME.class)

NUCLEOTIDEDIVERSITY.class@region.stats@nucleotide.diversity

#pop 1    pop 2
#pop 1 16948.74       NA
#pop 2 20101.25 18526.77
