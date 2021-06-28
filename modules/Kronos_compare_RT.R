#!/usr/local/bin/Rscript
suppressPackageStartupMessages(library(optparse))

options(stringsAsFactors = FALSE,scipen = 999)
options(warn = 1)
option_list = list(
    make_option(
        c("-R", "--RTs"),
        type = "character",
        default = NULL,
        help = "RT files with same binning",
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
        c('-D', "--deltaRT_threshold"),
        type = "double",
        default = 0.3,
        action = "store",
        help = "DeltaRT threshold to define changes [default= %default]",
        metavar = "double"
        ),
        make_option(
            c('-C', "--CrossingRT"),
            type = "logical",
            default = F,
            action = "store_false",
            help = "RT has to cross the 0.5 line to be considered as changing [default= %default]",
            metavar = "logical"
        ),
    make_option(
        c("-f", "--group_filter"),
        type = "character",
        default = NULL,
        help = "Filter out unwanted samples for RT files",
        metavar = "character"
    )
)

opt = parse_args(object = OptionParser(option_list = option_list))

#load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(Cairo))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GGally))

#set plotting theme
theme_set(theme_bw())

#create output directory
if (str_extract(opt$out,'.$')!='/'){
    opt$out=paste0(opt$out,'/')
}

system(paste0('mkdir -p ', opt$out))

if (!'RTs' %in% names(opt)) {
    stop("No RT file has been provided. See script usage (--help)")
}
#load files
opt$RTs = str_split(opt$RTs, ',', simplify = F)[[1]]

data <-
    foreach(
        i = opt$RTs,
        .combine = 'rbind',
        .packages = 'tidyverse'
    ) %do% {
        read_tsv(i,
                 col_types = cols())
    }

if('group_filter' %in% names(opt)){
    opt$group_filter=str_split(opt$group_filter, ',', simplify = F)[[1]]
    data=data%>%
        filter(!group %in% opt$group_filter)
}

data=data%>%
    spread(group,RT)%>%
    drop_na()

lines=names(data)[!names(data) %in% c('chr','start','end')]

distance=foreach(names1=1:(length(lines)-1),.combine = cbind)%:%
    foreach(names2=(names1+1):length(lines),.combine = cbind)%do%{
        
        data[lines[names1]]-data[lines[names2]]
        
    }
crossing=foreach(names1=1:(length(lines)-1),.combine = cbind)%:%
    foreach(names2=(names1+1):length(lines),.combine = cbind)%do%{
        
        sign((0.5-data[lines[names1]])/(0.5-data[lines[names2]]))
        
    }

keep=abs(distance)>=opt$deltaRT_threshold
crossing=crossing==-1
crossing=rowSums(crossing)>0
keep=rowSums(keep)>0
keep=  crossing & keep

clusters=kmeans(x =distance[keep,],centers =length(lines)^2 ,iter.max = 10000,nstart = )
data_cluster=data[keep,]
data_cluster$clusters=clusters$cluster

data_cluster$rowmean=rowMeans(distance[keep,])

data_cluster=data_cluster%>%
    arrange(clusters,rowmean)%>%
    group_by(clusters)%>%
    mutate(n=1:n())%>%
    gather(group,RT,names(data)[!names(data) %in% c('chr','start','end','clusters','n')])

p=data_cluster%>%
    ggplot(aes(x=group,y=n,hight=1,width=1,fill=RT))+
    geom_raster()+facet_grid(clusters ~ .,scales = 'free')+
    scale_fill_gradient2(low = '#005095',midpoint = 0.5,high = '#a7001b',mid = 'white')+
    theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title = element_blank(),panel.spacing = unit(0,'line'))


suppressMessages(
    ggsave(
        paste0(opt$out, '/Changing_bins_with_th_',opt$deltaRT_threshold,
               ifelse(opt$CrossingRT,'_and_crossing_0.5',''),
               '.pdf'),
        plot = p,
        limitsize = FALSE,
        device = cairo_pdf
    )
)
data_cluster%>%dplyr::select(-rowmean)%>% write_tsv(
paste0(opt$out, '/Changing_bins_with_th_',opt$deltaRT_threshold,
       ifelse(opt$CrossingRT,'_and_crossing_0.5',''),
       '.tsv'))

                          
print('done')
