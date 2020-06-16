#!/usr/local/bin/Rscript --slave
#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)
options(warn = 1,scipen = 999)

option_list = list(
    make_option(
        c("-F", "--file"),
        type = "character",
        default = NULL,
        help = "*calculated_Twhith* file(s) created by Kronos RT. If multiple files are provided they have to be separated by a comma.",
        metavar = "character"
    ),
    make_option(
        c("-A", "--annotation"),
        type = "character",
        default = NULL,
        help = "BedGraph file. No header.",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "output",
        help = "Output directory [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("-f", "--output_file_base_name"),
        type = "character",
        default = "out",
        help = "Base name for the output file [default= %default]",
        metavar = "character"
    )
)

#recover inputs
opt = parse_args(object = OptionParser(option_list = option_list))

suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))
suppressPackageStartupMessages(library(GenomicRanges, quietly = TRUE))
suppressPackageStartupMessages(library(ggpubr, quietly = TRUE))
suppressPackageStartupMessages(library(matrixStats, quietly = TRUE))

if (!'file' %in% names(opt)) {
    stop("Variability file must be provided. See script usage (--help)")
}

if (str_extract(opt$out, '.$') != '/') {
    opt$out = paste0(opt$out, '/')
}

system(paste0('mkdir -p ', opt$out))

opt$file = str_split(opt$file, ',')[[1]]

#load file
data <-
    foreach(file = opt$file,
            .combine = 'rbind',
            .packages = 'tidyverse') %do% {
                read_tsv(file, col_types = cols())
            }


if ('annotation' %in% names(opt)) {
    Annotation_file = read_tsv(opt$annotation,
                               col_names = c('chr', 'start', 'end', 'annotation'),col_types = cols()) %>%
        makeGRangesFromDataFrame(
            keep.extra.columns = T,
            seqnames.field = 'chr',
            end.field = 'end',
            start.field = 'start'
        )
    
    
    Bin = data %>%
        dplyr::select(chr, start, end) %>%
        unique() %>%
        makeGRangesFromDataFrame(
            keep.extra.columns = T,
            seqnames.field = 'chr',
            end.field = 'end',
            start.field = 'start'
        )
    
    #find overlaps
    hits = findOverlaps(Bin, Annotation_file)
    
    #indo about overlapping regins
    overlaps <-
        pintersect(Annotation_file[subjectHits(hits)], Bin[queryHits(hits)])
    
    #convert bins that have a notation and those that don't into df
    changed = as_tibble(Bin[queryHits(hits)])
    not_changed = as_tibble(Bin[-queryHits(hits)])
    
    #Based on the overla define the predominant notation of each bin.
    #A category to be chosen has to be predominant (at least 60% of the tatal overlaps in the bin)
    Annotation = bind_rows(
        changed %>%
            mutate(size = width(overlaps),
                   annotation = overlaps$annotation) %>%
            group_by(seqnames, start, end) %>%
            summarise(annotation = weightedMean(annotation,size)) ,
        not_changed%>%
            dplyr::select(seqnames,start,end)%>%
            mutate(annotation=0)
        
        
    )
    Annotation = Annotation %>%
        ungroup() %>%
        mutate(seqnames = as.character(seqnames))
    
    data = data %>%
        inner_join(Annotation, by = c("chr" = "seqnames", "start", "end"))
    

    plot = data %>%
        ggplot() +
        geom_point(aes(x = annotation, y = TW, color = basename), alpha = 0.2) +
        geom_smooth(aes(x = annotation, y = TW), color = 'black') +
        facet_wrap( ~ basename)
    
    suppressMessages(ggsave(
        plot,
        filename = paste0(
            opt$out,
            '/',
            opt$output_file_base_name,
            '_variability_plot_all_bins.pdf'
        )
    ))    
    }
    









print('done')
