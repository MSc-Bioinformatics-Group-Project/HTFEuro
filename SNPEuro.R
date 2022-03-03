#big boy no variables 
library(tidyverse)
library(PopGenome)
library(ggplot2)
library(bigmemory)

#### constant population order - ('British','Finnish','CEPH','Iberian','Toscani')
#read data
#set.populations


getpositiondata <- function(x,y){
  GENOME.class <<- readVCF('chr21.five.populations.vcf.gz',frompos = x, topos = y,tid = 'chr21',numcols = 1000)
  
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
  
  
  windows_df <<- data.frame(start = window_start, stop = window_stop, 
                         mid = window_start + (window_stop-window_start)/2)
  
  vcf_windows <<- sliding.window.transform(GENOME.class, width = window_size, jump = window_jump, type = 2)
  return(vcf_windows)
  return(windows_df)
}
############

#z specifies what population, 1,2,3,4,5 which pertains to the constant population order 

neutrality <- function(z){
  GENOME.class <- neutrality.stats(GENOME.class) #
  return(get.neutrality(GENOME.class)[[z]])
}

#analyse 

analysis <- function(window_size){
  
  #run all stat methods and store as variables to put into dataframe 
  
  vcf_windows <<- diversity.stats(vcf_windows, pi= TRUE)
  
  vcf_windows <<- F_ST.stats(vcf_windows)
  
  vcf_windows <<- neutrality.stats(vcf_windows)
  
  nd <<- vcf_windows@nuc.diversity.within/window_size
  
  pops <<- c('British','Finnish','CEPH','Iberian','Toscani')
  
  colnames(nd) <<- paste0(pops, "_pi")
  
  fst <<- t(vcf_windows@nuc.F_ST.pairwise)
  
  dxy <<- get.diversity(vcf_windows, between = T)[[2]]/window_size
  
  hpw <<-- t(vcf_windows@hap.diversity.between)
  
  tajd <<- vcf_windows@Tajima.D
}

pre_plot <- function(){
  
  #rename columns of each stat to be labelled with the populations for plots 
  
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
  paste0(x, "_hpw")
  paste0(x, "_tajd")
  
  colnames(fst) <- paste0(x, "_fst")
  colnames(dxy) <- paste0(x, "_dxy")
  colnames(hpw) <- paste0(x, "_hpw")
  colnames(tajd) <- c('British_tajd','Finnish_tajd','CEPH_tajd','Iberian_tajd','Toscani_tajd')
  vcf_data <<- as_tibble(data.frame(windows_df, nd, fst, dxy, hpw,tajd))
  
  return(vcf_data)
}



#specify what summary stats you want, return the mean across all windows 

#these must be fed a list c(a,b,c) with populations of interest i.e 'CEPH'


##### p

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

### hpw 

sum_hpw <- function(pops){
  vcf_data %>% select(contains("hpw")) %>% select(contains(pops)) %>% summarise_all(mean)
}


### tajd - extra because some populations return "nan" which ruins calculation

sum_tajd <- function(pops){
  tajd <- as.tibble(vcf_data)
  tajd <- tajd %>% mutate_all(~replace(., is.nan(.), 0))
  tajd %>% select(contains(pops)) %>% select(contains("tajd")) %>% summarise_all(mean)
}


#pi_plot

plot_pi <- function(pops){
  pi_g <- vcf_data %>% select(contains("pi")) %>% select(contains(pops)) %>% gather(key = "population", value = "pi")
  jpeg('static/piplot.jpg')
  myplot <- ggplot(pi_g, aes(population, pi)) + geom_boxplot() + theme_light() + xlab(NULL)
  print(myplot)
  dev.off()
}


#plot dxy

plot_dxy<- function(pops){
  dxy <- vcf_data %>% select(contains(pops)) %>% select(contains("dxy"))
  midcolumn <- vcf_data %>% select(mid)
  dxy <- cbind(midcolumn,dxy)
  
  
  dxy_g <- gather(dxy, -mid, key = "stat", value = "value")
  
  
  dxyplot <- ggplot(dxy_g, aes(mid/10^6, value, colour = stat)) + geom_line()
 
  #dxyplot <- dxyplot + facet_grid(stat~., scales = "free_y")
  dxyplot <- dxyplot + xlab("Position (Mb)")
  
  jpeg('static/dxyplot.jpg')
  print(dxyplot)
  dev.off()
}

#plot fst 
plot_fst<- function(pops){
  fst <- vcf_data %>% select(contains(pops)) %>% select(contains("fst"))
  midcolumn <- vcf_data %>% select(mid)
  fst <- cbind(midcolumn,fst)
  fst_g <- gather(fst, -mid, key = "stat", value = "value")
  fstplot <- ggplot(fst_g, aes(mid/10^6, value, colour = stat)) + geom_line()
  
  #fstplot <- fstplot + facet_grid(stat~., scales = "fixed")
  fstplot <- fstplot + xlab("Position (Mb)")
  jpeg('static/fstplot.jpg')
  print(fstplot)
  dev.off()
}

#plot hpw

plot_hpw<- function(pops){
  hpw <- vcf_data %>% select(contains(pops)) %>% select(contains("hpw"))#
  midcolumn <- vcf_data %>% select(mid)
  hpw <- cbind(midcolumn,hpw)
  hpw_g <- gather(hpw, -mid, key = "stat", value = "value")
  hpwplot <- ggplot(hpw_g, aes(mid/10^6, value, colour = stat)) + geom_line()
  #hpwplot <- hpwplot + facet_grid(stat~., scales = "free_y")
  hpwplot <- hpwplot + xlab("Position (Mb)")
  jpeg('static/hpwplot.jpg')
  print(hpwplot)
  dev.off()
}

#plot tajd

plot_tajd<- function(pops){
  tajd <- vcf_data %>% select(contains(pops)) %>% select(contains("tajd"))
  tajd <- as.tibble(tajd)
  
  #some positions for some pops have "nan", replace with 0 to plot 
  tajd <- tajd %>% mutate_all(~replace(., is.nan(.), 0))
  
  midcolumn <- vcf_data %>% select(mid)
  tajd <- cbind(midcolumn,tajd)
  tajd_g <- gather(tajd, -mid, key = "stat", value = "value")
  tajdplot <- ggplot(tajd_g, aes(mid/10^6, value, colour = stat)) + geom_line()
  #tajdplot <- tajdplot + facet_grid(stat~., scales = "free_y")
  tajdplot <- tajdplot + xlab("Position (Mb)")
  jpeg('static/tajdplot.jpg')
  print(tajdplot)
  dev.off()
}
