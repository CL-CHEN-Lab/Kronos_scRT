#!/usr/local/bin/Rscript

#install all the required R packages

install.packages(pkgs =c('Cairo','doSNOW','foreach','gplots','LaplacesDemon','MASS','matrixStats','optparse',
                   'RColorBrewer','scales','tidyverse','gridExtra','ggcorrplot','GGally','Rtsne','ade4','umap'), dependencies = T , repos = "http://cran.us.r-project.org")


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(pkgs = c("DNAcopy","Biostrings",'GenomicRanges','Rbowtie2','Rsamtools'))
