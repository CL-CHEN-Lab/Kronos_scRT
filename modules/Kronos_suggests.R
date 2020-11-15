#!/usr/local/bin/Rscript --slave
#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)
options(warn = 1, scipen = 999)

option_list = list(
    make_option(
        c("-R", "--single_cell_reads_info"),
        type = "character",
        default = NULL,
        help = "Single cell reads info file produced by Kronos fastqtoBAM or mean reads number followed by PE or SE ( for example: 100PE ). ",
        metavar = "character"
    ),make_option(
        c("-G", "--genome_size"),
        type = "numeric",
        default = 2900,
        help = "Genome size in MB [default= %default]",
        metavar = "numeric"
    ),
    make_option(
        c("-P", "--expected_ploidy"),
        type = "numeric",
        default = 2,
        help = "Expected ploidy [default= %default]",
        metavar = "numeric"
    )
)
#recover inputs
opt = parse_args(object = OptionParser(option_list = option_list))

#load libraries
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))

#check inputs
if(!'single_cell_reads_info' %in% names(opt)){
    
    stop('Single cell reads info file or manual info not provided. See script usage (--help)')
    
}else if(file.exists(opt$single_cell_reads_info)){
    info_file=read_tsv(opt$single_cell_reads_info, col_names = T,col_types = cols())
    
    PE=info_file%>%
        mutate(PE=UNPAIRED_READS_EXAMINED > UNPAIRED_READS_EXAMINED)%>%
        group_by(PE)%>%
        summarise(n=n())%>%
        filter(n==max(n))%>%pull(PE)
    
    median=info_file%>%
        mutate(reads=(UNPAIRED_READS_EXAMINED - UNPAIRED_READ_DUPLICATES)+(1+PE)*(READ_PAIRS_EXAMINED -  READ_PAIR_DUPLICATES ))%>%
        select(reads)%>%
        mutate(reads=reads/(opt$genome_size*opt$expected_ploidy))%>%
        summarise(median=median(reads))%>%
        pull(median)
                       
                     
}else{
    opt$single_cell_reads_info = str_to_upper(opt$single_cell_reads_info)
    if(str_detect(string =opt$single_cell_reads_info,pattern = '[PS]E')){
    PE=str_detect(string =opt$single_cell_reads_info,pattern = 'PE')
    median=as.numeric(str_extract(string =opt$single_cell_reads_info,pattern = '[0-9]{1,1000}'))
    }else{
        stop('Provided manual info lack of PE or SE information. See script usage (--help)')
    }
}

size_bin=as.integer((20*(56*(1+PE))/median))

if(size_bin >= 5 ){
    size_bin=trunc(size_bin/5)*5
}else if (size_bin >= 1 & size_bin <5 ) {
    size_bin=trunc(size_bin)
}else if(size_bin < 1){
    size_bin=1
}


resolution=size_bin*500/20
min_n_reads=round((20*(56*(1+PE)))/size_bin)
line='\n\n\n-----------------------------------\n\n\n'
cat(paste0(line,'Bins size should be around ',size_bin,'kb ',ifelse(PE,'PE','SE'),' mode. \nFinal resolution for RT can be increased to ',resolution,'kb.\nChange minimum number of reads per cell to ',min_n_reads,' reads per megabase per ploidy.\n',line))
print('done')
