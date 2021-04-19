#!/usr/local/bin/Rscript
#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)
options(warn = 1,scipen = 999)

option_list = list(
  make_option(
    c("-F", "--file"),
    type = "character",
    default = NULL,
    help = "Variability file produced by Kronos RT, if multiple files are provided they have to be separated by a comma",
    metavar = "character"
  ),
  make_option(
    c("-R", "--Annotation"),
    type = "character",
    default = NULL,
    help = "Genome annotation. chr<TAB>start<TAB>end<TAB>annotation. No header.",
    metavar = "character"
  ),
  make_option(
    c("-r", "--Annotation2"),
    type = "character",
    default = NULL,
    help = "Second genome annotation. chr<TAB>start<TAB>end<TAB>annotation. No header.",
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
  ),
  make_option(
    c("-m", "--min_overlap"),
    type = "numeric",
    default = 100,
    help = "Min overlap to apply the annotation in bp [default= %default]",
    metavar = "numeric"
  )
)

#recover inputs
opt = parse_args(object = OptionParser(option_list = option_list))

#load packages
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))
suppressPackageStartupMessages(library(doSNOW, quietly = TRUE))
suppressPackageStartupMessages(library(GenomicRanges, quietly = TRUE))

#set plotting theme
theme_set(theme_bw())

#check inputs
if (!'file' %in% names(opt)) {
  stop("Variability file must be provided. See script usage (--help)")
}
if (!'Annotation' %in% names(opt)) {
  stop("Annotation file must be provided. See script usage (--help)")
}
if (str_extract(opt$out, '.$') != '/') {
  opt$out = paste0(opt$out, '/')
}

system(paste0('mkdir -p ', opt$out))

opt$file = str_split(opt$file, ',')[[1]]

#load files
data <-
  foreach(file = opt$file,
          .combine = 'rbind',
          .packages = 'tidyverse') %do% {
            read_tsv(file, col_types = cols())
          }

Annotation_file = read_tsv(opt$Annotation,
                           col_names = c('chr', 'start', 'end', 'annotation'), col_types = cols())%>%
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
hits = findOverlaps(Bin, Annotation_file,minoverlap = opt$min_overlap)

#indo about overlapping regins
overlaps <-
  pintersect(Annotation_file[subjectHits(hits)], Bin[queryHits(hits)])

#Add annotation 1
add = as_tibble(Bin[queryHits(hits)])
not_add = as_tibble(Bin[-queryHits(hits)])

#Based on the overlap define the predominant notation of each bin.
#A category to be chosen has to be predominant (at least 60% of the tatal overlaps in the bin)
Annotation = bind_rows(
  add %>%
    mutate(size = width(overlaps),
           Cat1 = overlaps$annotation) %>%
    group_by(seqnames, start, end, Cat1) %>%
    summarise(n = sum(size)) %>%
    ungroup() %>%
    group_by(seqnames, start, end) %>%
    mutate(
      select = (n / sum(n) >= 0.6),
      Cat1 = ifelse(any(select), Cat1, '_Unknown_')
    ) %>%
    dplyr::select(seqnames, start, end, Cat1),
  not_add %>% mutate(Cat1 = '_Unknown_') %>%
    dplyr::select(seqnames, start, end, Cat1)
)
Annotation = Annotation %>%
  ungroup() %>%
  mutate(seqnames = as.character(seqnames))

data = data %>%
  inner_join(Annotation, by = c("chr" = "seqnames", "start", "end"))

if ('Annotation2' %in% names(opt)) {
  
  Annotation_file = read_tsv(opt$Annotation2,
                             col_names = c('chr', 'start', 'end', 'annotation'), col_types = cols()) %>%
    makeGRangesFromDataFrame(
      keep.extra.columns = T,
      seqnames.field = 'chr',
      end.field = 'end',
      start.field = 'start'
    )
  
  hits = findOverlaps(Bin, Annotation_file,minoverlap = opt$min_overlap)
  
  overlaps <-
    pintersect(Annotation_file[subjectHits(hits)], Bin[queryHits(hits)])
  
  add = as_tibble(Bin[queryHits(hits)])
  not_add = as_tibble(Bin[-queryHits(hits)])
  
  Annotation2 = bind_rows(
    add %>%
      mutate(
        size = width(overlaps),
        Cat2 = overlaps$annotation
      ) %>%
      group_by(seqnames, start, end, Cat2) %>%
      summarise(n = sum(size)) %>%
      ungroup() %>%
      group_by(seqnames, start, end) %>%
      mutate(
        select = (n / sum(n) >= 0.6),
        Cat2 = ifelse(any(select), Cat2, '_Unknown_')
      ) %>%
      dplyr::select(seqnames, start, end, Cat2),
    not_add %>% mutate(Cat2 = '_Unknown_') %>%
      dplyr::select(seqnames, start, end, Cat2)
  )
  Annotation2 = Annotation2 %>%
    ungroup() %>%
    mutate(seqnames = as.character(seqnames))
  
  data = data %>%
    inner_join(Annotation2, by = c("chr" = "seqnames", "start", "end"))
  
  data = data %>%
    rbind(data %>%
            mutate(Cat1 = '_ALL_'),
          data %>%
            mutate(Cat2 = '_ALL_'),
          data %>%
            mutate(Cat1 = '_ALL_',
                   Cat2 = '_ALL_'))
}


data%>%write_tsv(paste0(opt$out,
                 '/',
                 opt$output_file_base_name,
                 '_scRT_scRT_variability_with_annotation.tsv'))

print('done')

