#!/usr/local/bin/Rscript
#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)
options(warn = 1, scipen = 999)

option_list = list(
    make_option(
        c("-F", "--File"),
        type = "character",
        default = NULL,
        help = "Replication timing files separated by a comma. Format: chr <TAB> start <TAB> end <TAB> group",
        metavar = "character"
    ),
    make_option(
        c("-s", "--sort"),
        type = "character",
        default = NULL,
        help = "Group names orders",
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

#load module
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(GGally, quietly = TRUE))
suppressPackageStartupMessages(library(ggcorrplot, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))


if (!'File' %in% names(opt)) {
    stop("RT files were not provided. See script usage (--help)")
    
} else{
    opt$File = str_split(opt$File, ',')[[1]]
}
if ('sort' %in% names(opt)) {
    opt$sort = str_split(opt$sort, ',')[[1]]
    
}

#create directory
if (str_extract(opt$out, '.$') != '/') {
    opt$out = paste0(opt$out, '/')
}

system(paste0('mkdir -p ', opt$out))

scRT = foreach(
    i = 1:length(opt$File),
    .packages = 'tidyverse',
    .combine = 'rbind'
) %do% {
    tmp = read_tsv(opt$File[i], col_types = cols())
    if ('sort' %in% names(opt)) {
        tmp %>%
            mutate(group = factor(group, levels = opt$sort))
        
    } else{
        tmp
    }
}


#create directory
if (str_extract(opt$out, '.$') != '/') {
    opt$out = paste0(opt$out, '/')
}



scRT = scRT %>% spread(group, RT) %>%
    drop_na() %>%
    dplyr::select(-chr, -start, -end) 

plot = ggcorrplot(
    scRT%>%
        cor(),
    lab = T,
    lab_col = 'white',legend.title = 'Pearson\ncorrelation',
    colors = c('#21908CFF', '#F0F921FF', '#BB3754FF')
)

suppressMessages(ggsave(
    plot = plot,
    filename = paste0(
        opt$out,
        opt$output_file_base_name,
        'pearson_correlation',
        '.pdf'
    )
))
suppressMessages( ggsave(
    plot = ggpairs(scRT, aes(alpha = 0.3),lower = list(continuous = wrap("density", alpha = 0.5), combo = "box_no_facet"))+theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1)),
    filename =paste0(
        opt$out,
        '/',
        opt$output_file_base_name,
        'pair_scatter_plot_RTs.pdf'
    )
))

print('done')

