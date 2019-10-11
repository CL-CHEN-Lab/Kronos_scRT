#!/usr/local/bin/Rscript
if (!suppressPackageStartupMessages(require(optparse, quietly = TRUE))) {
    install.packages("optparse", quiet = T)
    suppressPackageStartupMessages(library(optparse))
}
options(stringsAsFactors = FALSE)
options(warn = 1)
option_list = list(
    make_option(
        c("-S", "--S50s"),
        type = "character",
        default = NULL,
        help = "RT files with same binning",
        metavar = "character"
    ),
    make_option(
        c("-R", "--referenceRT"),
        type = "character",
        default = NULL,
        help = "Reference RT min=Late, max=Early, only one reference is allowed",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "output",
        help = "output directory [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("-k", "--keepXY"),
        type = "logical",
        default = FALSE,
        action = "store_true",
        help = "keeps XY chr in the analysis",
        metavar = "logical"
    ),
    make_option(
        c("--Reference"),
        type = "character",
        action = "store",
        help = "base name to use as a reference, if not provided the first basename in the S50 file will be used or , if provided , the reference RT even if this option is selected",
        metavar = "character"
    ),
    make_option(
        c('-D', "--deltaRT_threshold"),
        type = "double",
        default = 0.1,
        action = "store",
        help = "detalRT threshold to define changes",
        metavar = "double"
    ),
    make_option(
        c('-n', "--n_regions"),
        type = "integer",
        action = "store",
        default = 10,
        help = "number of regions to plot",
        metavar = "integer"
    ),
    make_option(
        c("-r", "--region"),
        type = "character",
        default = NULL,
        help = "Region to plot  chr:start-end (multiple regins can be separated by a comma)",
        metavar = "character"
    ),
    make_option(
        c("-f", "--basename_filter"),
        type = "character",
        default = NULL,
        help = "filter out unwanted samples for RT files",
        metavar = "character"
    )
)

opt = parse_args(object = OptionParser(option_list = option_list))

if (!suppressPackageStartupMessages(require(BiocManager, quietly = TRUE))) {
    install.packages("BiocManager", quiet = T)
}
if (!suppressPackageStartupMessages(require(tidyverse, quietly = TRUE))) {
    install.packages("tidyverse", quiet = T)
    suppressPackageStartupMessages(library(tidyverse))
}

if (!suppressPackageStartupMessages(require(foreach, quietly = TRUE))) {
    install.packages("foreach", quiet = T)
    suppressPackageStartupMessages(library(foreach))
}
if (!suppressPackageStartupMessages(require(Cairo, quietly = TRUE))) {
    install.packages("Cairo", quiet = T)
    suppressPackageStartupMessages(library(Cairo))
}
if (!suppressPackageStartupMessages(require(scales))) {
    install.packages("scales", quiet = T)
    suppressPackageStartupMessages(library(scales))
}
if (!suppressPackageStartupMessages(require(GenomicRanges))) {
    BiocManager::install("GenomicRanges")
    suppressPackageStartupMessages(library(GenomicRanges))
}

#create output directory
if (str_extract(opt$out,'.$')!='/'){
    opt$out=paste0(opt$out,'/')
}

system(paste0('mkdir -p ./', opt$out))

if (!'S50s' %in% names(opt)) {
    stop("No S50 file has been provided. See script usage (--help)")
}
#load files
opt$S50s = str_split(opt$S50s, ',', simplify = F)[[1]]

data <-
    foreach(
        i = opt$S50s,
        .combine = 'rbind',
        .packages = 'tidyverse'
    ) %do% {
        read_tsv(i,
                 col_types = cols())
    }

if('basename_filter' %in% names(opt)){
    opt$basename_filter=str_split(opt$basename_filter, ',', simplify = F)[[1]]
    data=data%>%
        filter(!basename %in% opt$basename_filter)
}

if ('referenceRT' %in% names(opt)) {
    referenceRT = read_tsv(opt$referenceRT,
                           col_types = cols())
    names(referenceRT) = c(names(referenceRT)[-4], 'referenceRT')
    data = inner_join(data, referenceRT, by = c('chr', 'start', 'end'))
    if (length(data$chr) == 0) {
        stop('Reference RT and S50 files do not have the same binning')
    }
    opt$Reference = 'referenceRT'
} else{
    UB = unique(data$basename)
    if ('Reference' %in% names(opt)) {
        if (!opt$Reference %in% UB | !'Reference' %in% names(opt))
            opt$Reference = UB[1]
        warning(
            paste0(
                'Reference basename not found, ',
                opt$Reference,
                ' will be used as reference'
            )
        )
        
    }
    data = spread(data, basename, RT) %>%
        gather(key = 'basename', value = 'RT', UB[!UB %in% opt$Reference]) %>%
        filter(complete.cases(.))
}


regions = data %>%
    group_by(basename) %>%
    summarise(n = n())

percentages = foreach (i = seq(0, 1, 0.01),
                       .combine = 'rbind',
                       .packages = 'tidyverse') %do% {
                           type = function(RT.x, RT.y, i) {
                               RT.x = ((RT.x - min(RT.x)) / (max(RT.x) - min(RT.x))) - 0.5
                               RT.y = ((RT.y - min(RT.y)) / (max(RT.y) - min(RT.y))) - 0.5
                               delta = RT.x - RT.y
                               result = ifelse((RT.x / RT.y < 0 &
                                                    RT.y > 0 & abs(delta) > i),
                                               'EtoL',
                                               ifelse((RT.x / RT.y < 0 &
                                                           RT.y < 0 & abs(delta) > i),
                                                      'LtoE',
                                                      ifelse((RT.x / RT.y >= 0 &
                                                                  delta < -i),
                                                             'toLater',
                                                             ifelse((RT.x / RT.y >= 0 &
                                                                         delta > i), 'toEarlier', 'unchanged')
                                                      )
                                               ))
                               return(result)
                           }
                           
                           data %>%
                               mutate(type = type(RT, select(., matches(opt$Reference)), i)) %>%
                               group_by(basename, type) %>%
                               summarise(counts = n()) %>%
                               inner_join(regions, by = "basename") %>%
                               mutate(percent = counts / n) %>%
                               select(basename, type, percent) %>%
                               mutate(th = i)
                           
                       }

p = percentages %>%
    ungroup() %>%
    spread(type, percent, fill = 0) %>%
    gather(type, percent, -basename, -th) %>%
    filter(type != 'unchanged') %>%
    mutate(line = paste0('ΔRT = ', basename, ' - ', opt$Reference)) %>%
    ggplot(aes(x = th, y = percent, color = type)) +
    geom_line() +
    facet_grid( ~ line) +
    scale_fill_manual(values = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a')) +
    scale_color_manual(values = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a')) +
    ylab('% ofchanging bins') + xlab('ΔRT threshold')+
    scale_y_continuous(labels = scales::percent_format())

suppressMessages(
    ggsave(
    paste0(opt$out, '/changes_distribution.pdf'),
    plot = p,
    limitsize = FALSE,
    device = cairo_pdf
)
)


type = function(RT.x, RT.y, i) {
    RT.x = ((RT.x - min(RT.x)) / (max(RT.x) - min(RT.x))) - 0.5
    RT.y = ((RT.y - min(RT.y)) / (max(RT.y) - min(RT.y))) - 0.5
    delta = RT.x - RT.y
    result = ifelse((RT.x / RT.y < 0 &
                         RT.y > 0 & abs(delta) > i),
                    'EtoL',
                    ifelse((RT.x / RT.y < 0 &
                                RT.y < 0 & abs(delta) > i),
                           'LtoE',
                           ifelse((RT.x / RT.y >= 0 &
                                       delta < -i),
                                  'toLater',
                                  ifelse((RT.x / RT.y >= 0 &
                                              delta > i), 'toEarlier', 'unchanged')
                           )
                    ))
    return(result)
}

changing = data %>%
    mutate(
        type = type(RT, select(., matches(opt$Reference)), opt$deltaRT_threshold),
        ΔRT = opt$deltaRT_threshold
    ) %>%
    `colnames<-`(c(names(data), 'type', 'ΔRT')) %>%
    filter(type != 'unchanged')

write_tsv(changing, paste0(opt$out, '/changing_regions.tsv'))

changing_names = changing
changing = makeGRangesFromDataFrame(
    df = changing,
    keep.extra.columns = F ,
    seqnames.field = 'chr',
    start.field = 'start',
    end.field = 'end'
)
changing = changing %>%
    GenomicRanges::reduce(min.gapwidth=5000000) %>%
    data.frame() %>%
    select(seqnames, start, end) %>%
    `colnames<-`(c('chr', 'start', 'end'))%>%
    mutate(start=start,
           end=end)

if (!'region' %in% names(opt)) {
    for (i in sample(1:length(changing$chr), opt$n_regions)) {
        to_plot=changing[i,]
        if(to_plot$end-to_plot$start < 20000000){
            add=20000000-(to_plot$end-to_plot$start)
            to_plot=to_plot%>%
                mutate(start=start-add/2,
                       end=end+add/2)
        }
        p = data %>%
            mutate(end = end - 1) %>%
            filter(
                chr == to_plot$chr,
                start >= to_plot$start ,
                end <= to_plot$end 
            ) %>%
            gather('type', 'pos', start, end) %>%
            ggplot() +
            geom_hline(yintercept = 0.5) +
            geom_line(aes_string(
                x = 'pos',
                y = opt$Reference,
                color = shQuote(opt$Reference)
            )) +
            geom_line(aes(
                x = pos,
                y = RT,
                color = basename
            )) + facet_grid(basename ~ chr) +
            geom_rect(
                data = changing_names %>%
                    filter(
                        chr == to_plot$chr ,
                        start >= to_plot$start  ,
                        end <= to_plot$end 
                    ),
                aes(
                    xmin = start,
                    xmax = end,
                    fill = type
                ),
                alpha = 0.2,
                ymin = -Inf,
                ymax = Inf,
                inherit.aes = F
            ) +
            coord_cartesian(y = c(0, 1)) +
            scale_fill_manual(
                values = c(
                    'EtoL' = 'darkgreen',
                    'LtoE' = 'darkred',
                    'toEarlier' = 'red',
                    'toLater' = 'green'
                )
            )+
            ylab('RT')+xlab('')+ 
            scale_x_continuous(
                labels = trans_format(format = 'Mb',
                                      trans = function(x) paste((x/1000000),'Mb',sep='')))+ 
            theme(legend.title=element_blank(),axis.text.x = element_text(hjust = 1,angle = 45),legend.position = 'top')
        p
        suppressMessages(
        ggsave(p, filename = paste0(opt$out, '/changing_region_', paste(changing[i,], collapse = '_'), '.pdf'))
        )
    }
} else{
    opt$region = data.frame(coord = str_split(opt$region, pattern = ',')[[1]]) %>%
        separate(coord, c('chr', 'pos'), ':') %>%
        separate(pos, c('start', 'end'), '-')
    
    for (i in 1:length(opt$region$chr)) {
        p = data %>%
            mutate(start = start + 1) %>%
            filter(chr == opt$region$chr[i],
                   start >= opt$region$start[i],
                   end <= opt$region$end[i]) %>%
            gather('type', 'pos', start, end) %>%
            ggplot() +
            geom_hline(yintercept = 0.5) +
            geom_line(aes_string(
                x = 'pos',
                y = opt$Reference,
                color = shQuote(opt$Reference)
            )) +
            geom_line(aes(
                x = pos,
                y = RT,
                color = basename
            )) + facet_grid(basename ~ chr) +
            geom_rect(
                data = changing_names %>%
                    filter(
                        chr == opt$region$chr[i] ,
                        start >= opt$region$start[i] ,
                        end <= opt$region$end[i]
                    ),
                aes(
                    xmin = start,
                    xmax = end,
                    fill = type
                ),
                alpha = 0.2,
                ymin = -Inf,
                ymax = Inf,
                inherit.aes = F
            ) +
            coord_cartesian(y = c(0, 1)) +
            scale_fill_manual(
                values = c(
                    'EtoL' = 'darkgreen',
                    'LtoE' = 'darkred',
                    'toEarlier' = 'red',
                    'toLater' = 'green'
                )
            )+
            ylab('RT')+xlab('')+ 
            scale_x_continuous(
                labels = trans_format(format = 'Mb',
                                      trans = function(x) paste((x/1000000),'Mb',sep='')))+ 
            theme(legend.title=element_blank(),axis.text.x = element_text(hjust = 1,angle = 45),legend.position = 'top')
        
        suppressMessages(ggsave(p,
               filename = paste0(
                   opt$out,
                   '/changing_region_',
                   paste(opt$region[i,], collapse = '_'),
                   '.pdf'
               )))
        
    }
    
}