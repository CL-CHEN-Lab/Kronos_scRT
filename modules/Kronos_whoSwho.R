#!/usr/local/bin/Rscript --slave

# this script is meant to select the treshold to select cycling cells

suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)
options(warn=1) 
option_list = list(
    make_option(
        c("-F", "--file"),
        type = "character",
        default = NULL,
        help = "per cell stat file path",
        metavar = "character"
    ),
    make_option(
        c("-W", "--whoSwho"),
        type = "character",
        default = NULL,
        help = "who's who file path",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "./output",
        help = "Output directory [default= %default]",
        metavar = "character"
    )
)

opt = parse_args( OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))


#check inputs
if (!'file' %in% names(opt)) {
    stop("Per cell stat file must be provided. See script usage (--help)")
}
if (!'whoSwho' %in% names(opt)) {
    stop("who's who file must be provided. See script usage (--help)")
}
#create output directory
if (str_extract(opt$out,'.$')!='/'){
    opt$out=paste0(opt$out,'/')
}

system(paste0('mkdir -p ./', opt$out))

#load data
data<-inner_join(read_csv(opt$file,
               col_types = cols()),
                   read_tsv(opt$whoSwho,
                        col_types = cols()), by = "Cell")

# write data
data%>%
    mutate(is_high_dimapd=ifelse(Phase=='S',T,F))%>%
    dplyr::select(-Phase)%>%
    write.csv(paste0(opt$out,'phased_',basename(opt$file)))
