#big boy no variables 
library(tidyverse)
library(PopGenome)
library(ggplot2)

#### constant population order - ('British','Finnish','CEPH','Iberian','Toscani')
#read data
#set.populations
getpositiondata <- function(x,y){
  GENOME.class <<- readVCF('chopped.vcf.gz',frompos = x, topos = y,tid = 'chr21',numcols = 1000)
  
  Br <<- read.csv('British.csv')
  Fin <<- read.csv('Finnish.csv')
  CEPH <<- read.csv('CEPH.csv')
  Ibe <<- read.csv('Iberian.csv')
  Tos <<- read.csv('Toscani.csv')
  
  British <<- Br$British
  Finnish <<- Fin$Finnish
  CEPH <<- CEPH$CEPH
  Iberian <<- Ibe$Iberian
  Toscani <<- Tos$Toscani
  
  GENOME.class <<- set.populations(GENOME.class,new.populations = list(c(British),c(Finnish),c(CEPH),c(Iberian),c(Toscani)))
  
  return(GENOME.class)
}

#get # of SNPs

SNPCount <- function(){
  x <- GENOME.class@n.biallelic.sites + GENOME.class@n.polyallelic.sites
  return(x)
}

#get summary data

sumdata <- function(){
  x <- get.sum.data(object = GENOME.class)
  return(x)
}


#set up sliding window 

SlidingWindow <- function(window_size,window_jump,vcf_size){
  
  window_start <<- seq(from = 1, to = vcf_size, by = window_jump)
  window_stop <<- window_start + window_size
  
  window_start <<- window_start[which(window_stop < vcf_size)]
  window_stop <<- window_stop[which(window_stop < vcf_size)]
  
  
  windows <<- data.frame(start = window_start, stop = window_stop, 
                         mid = window_start + (window_stop-window_start)/2)
  
  vcf_windows <<- sliding.window.transform(GENOME.class, width = window_size, jump = window_jump, type = 2)
  return(vcf_windows)
  return(windows)
}
############

#z specifies what population, 1,2,3,4,5 which pertains to the constant population order 

neutrality <- function(z){
  GENOME.class <- neutrality.stats(GENOME.class) #
  return(get.neutrality(GENOME.class)[[z]])
}

#analyse 

analysis <- function(window_size){
  
  vcf_windows <<- diversity.stats(vcf_windows, pi= TRUE)
  
  vcf_windows <<- F_ST.stats(vcf_windows)
  
  nd <<- vcf_windows@nuc.diversity.within/window_size
  
  pops <<- c('British','Finnish','CEPH','Iberian','Toscani')
  
  colnames(nd) <<- paste0(pops, "_pi")
  
  fst <<- t(vcf_windows@nuc.F_ST.pairwise)
  
  dxy <<- get.diversity(vcf_windows, between = T)[[2]]/window_size
}

pre_plot <- function(){
  x <- colnames(fst)
  # replace all occurrences of pop1 with house
  x <- sub("pop1", "British", x)
  x <- sub("pop2", "Finnish", x)
  x <- sub("pop3", "CEPH", x)
  x <- sub("pop4", "Iberian", x)
  x <- sub("pop5", "Toscani", x)
  x <- sub("/", "_", x)
  
  paste0(x, "_fst")
  paste0(x, "_dxy")
  
  colnames(fst) <- paste0(x, "_fst")
  colnames(dxy) <- paste0(x, "_dxy")
  
  vcf_data <<- as_tibble(data.frame(windows, nd, fst, dxy))
  return(vcf_data)
}


#specify what summary stats you want

##### pi

#these must be fed a list c(a,b,c) with populations of interest i.e 'CEPH'

populations <-c('British','Toscani','CEPH')

sum_pi <- function(pops){
  vcf_data %>% select(contains("pi")) %>% select(contains(pops)) %>% summarise_all(mean)
}


#### fst

sum_fst <- function(pops){
  vcf_data %>% select(contains("fst")) %>% select(contains(pops)) %>% summarise_all(mean)
}


#### dxy

sum_dxy <- function(pops){
  vcf_data %>% select(contains("dxy")) %>% select(contains(pops)) %>% summarise_all(mean)
}

#pi_plot

plot_pi <- function(pops){
  pi_g <- vcf_data %>% select(contains("pi")) %>% select(contains(populations)) %>% gather(key = "population", value = "pi")
  pdf('myplot.pdf')
  myplot <- ggplot(pi_g, aes(population, pi)) + geom_boxplot() + theme_light() + xlab(NULL)
  print(myplot)
  dev.off()
}

