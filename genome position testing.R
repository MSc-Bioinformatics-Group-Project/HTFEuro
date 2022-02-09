#have a go at reading VCF start end pos
#installed WhopGenome on terminal R CMD INSTALL Desktop/[pkg.tar.gz]

library('WhopGenome')

GENOMEPOSITION.class <- readVCF('chr21positiontest.vcf.gz',frompos = 10519265 , topos = 10526343,tid = 'chr21',numcols = 1000)

vcf_handle   <- vcf_open("CH22AASWPhased.vcf.gz")
