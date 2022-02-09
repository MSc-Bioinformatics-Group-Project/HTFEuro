library(PopGenome)
snp <- readData("test1", gffpath = 'testgff', format = "VCF")

# You can access the different "slots" by using the "@" sign:
snp@n.sites

# Set populations

Br <- read.csv('British.csv')
Aa <- read.csv('AASW.csv')

British <- Br$British
AASW <- Aa$AASW

snp <-set.populations(GENOME.class,new.populations = list(c(British),c(AASW)))
snp@populations

# Diversities and FST (by scaffold)
snp <- F_ST.stats(snp) # this does the calculations and 
# adds the results to the appropriate slots

get.F_ST(snp) # each line is a scaffold
snp@nucleotide.F_ST

get.diversity(snp)
get.diversity(snp)[[1]] # pop1 (B)
get.diversity(snp)[[2]] # pop2 (b)
snp@nuc.diversity.within

#HAD TO CHANGE TYPE = 1 TO ASSESS SNPS NOT WHOLE GENOME

win_snp <- sliding.window.transform(snp, 
                                    width=10000, jump=10000, 
                                    type=1,
                                    whole.data=FALSE)

# Measurements per window
win_snp <- F_ST.stats(win_snp)

win_snp@nucleotide.F_ST
win_snp@nuc.diversity.within[1][1]

# A simple plot

#NEED TO WORK OUT EXACTLY WHAT THESE ARE EH

win_fst <- win_snp@nucleotide.F_ST[,1]
British_div  <- win_snp@nuc.diversity.within[,1]
AA_div  <- win_snp@nuc.diversity.within[,2] 

plot(1:length(win_fst), win_fst)

par(mfrow=c(2,1))
win_fst <- win_snp@nucleotide.F_ST[,1]
plot(1:length(AA_div), AA_div)

win_fst <- win_snp@nucleotide.F_ST[,1]
plot(1:length(British_div), British_div)
