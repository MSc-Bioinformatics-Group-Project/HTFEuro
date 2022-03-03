# SNPEuro

SNPEuro(Single Nucleotide Polymorphism Europe) is created as part of the Bioinformatics Group Software Development Project during the Queen Mary University of London Bioinformatics MSc.

User should be able to retrive single nucleotide polymorphism (SNP) information given either a genomic coordinate (chromosome, start and end), SNP name (rs value), or gene name (or any aliases associated to it).

SNPEUro will return the following information for each SNP: name (rs value), genomic position, genotype frequencies, and allele frequency. Frequencies should be provided for each population separately.

If multiple SNPs are returned, the user should be able to select the population(s) and summary statistics of interest, and the application will calculate them and plot their distribution in sliding-windows along the region. If multiple populations are selected, then population genetic variation (FST value) for each pair of populations should be reported.


# Installing SNPEuro

### Instructions to install SNPEuro ###

1. `Git clone git@github.com:MSc-Bioinformatics-Group-Project/SNPEuro.git`
