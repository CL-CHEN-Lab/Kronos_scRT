#!/usr/local/bin/Rscript --slave


suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)
options(warn = 1)

option_list = list(
    make_option(
        c("-F", "--file"),
        type = "character",
        default = NULL,
        help = "Per cell stat file , if multiple files are provided they have to be separated by a comma",
        metavar = "character"
    ),
    make_option(
        c("-T", "--tracks"),
        type = "character",
        default = NULL,
        help = "Tracks file,  if multiple files are provided they have to be separated by a comma",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "output",
        help = "Output directory [default= %default]",
        metavar = "character"
    )
)

#recover inputs
opt = parse_args(object = OptionParser(option_list = option_list))

suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))

#create directory

if (str_extract(opt$out,'.$')!='/'){
    opt$out=paste0(opt$out,'/')
}

system(paste0('mkdir -p ', opt$out))

#load files
opt$file = str_split(opt$file, ',')[[1]]

opt$tracks = str_split(opt$tracks, ',')[[1]]
if ('file' %in% names(opt)) {
    for (file in opt$file) {
        read_csv(file, col_types = cols()) %>%
            select(cell_id,
                   normalized_dimapd,
                   mean_ploidy,
                   ploidy_confidence,
                   is_high_dimapd,
                   is_noisy,
                   effective_reads_per_1Mbp) %>%
            `colnames<-`(c(
                'Cell',
                'normalized_dimapd',
                'mean_ploidy',
                'ploidy_confidence',
                'is_high_dimapd',
                'is_noisy',
                'coverage_per_1Mbp'
            )) %>%
            write_csv(paste0(opt$out, '/Kronos_format_', basename(file)))
    }
}
if ('tracks' %in% names(opt)) {
    for (tracks in opt$tracks) {
        read_tsv(tracks, skip = 2, col_types = cols()) %>%
            select(id, `#chrom`, start, end, copy_number) %>%
            `colnames<-`(c('Cell', 'chr', 'start', 'end', 'copy_number')) %>%
	mutate(reads= '10X')%>%
            write_tsv(paste0(opt$out, '/Kronos_format_', basename(tracks)))
    }
}

if (!('tracks' %in% names(opt) & 'file' %in% names(opt) )){
    stop('No input')
}

print('done')
