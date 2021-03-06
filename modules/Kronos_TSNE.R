#!/usr/local/bin/Rscript
#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)
options(warn = 1, scipen = 999)

option_list = list(
    make_option(
        c("-C", "--scCNV"),
        type = "character",
        default = NULL,
        help = "*single_cells_CNV* file(s) created by Kronos RT. If multiple files are provided they have to be separated by a comma.",
        metavar = "character"
    ),
    make_option(
        c("--CNV_values"),
        type = "character",
        default = "B",
        help = "What type of date to plot for the sigle cell traks: ('B'=Binarized, 'CNV'=Copy number variation, 'log2'=log2(CNV_Cell/CNV_mean_G1/G2_cells)) [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("--per_Chr"),
        type = 'logical',
        default = F,
        action = 'store_TRUE',
        help = "Calculate TSNE on each chromosome",
        metavar = "logical"
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
        c("-c", "--cores"),
        type = "integer",
        default = 3,
        help = "Numbers of cores to use [default= %default]",
        metavar = "integer"
    )
)

#recover inputs
opt = parse_args(object = OptionParser(option_list = option_list))

#load packages
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(RColorBrewer, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))
suppressPackageStartupMessages(library(doSNOW, quietly = TRUE))
suppressPackageStartupMessages(library(Rtsne, quietly = TRUE))
suppressPackageStartupMessages(library(ade4, quietly = TRUE))
#output dir
system(paste0('mkdir -p ', opt$out))

#set plotting theme
theme_set(theme_bw())

#load files
if (!'scCNV' %in% names(opt)) {
        stop("scCNV file must be provided. See script usage (--help)")
}
    
opt$scCNV = str_split(opt$scCNV, ',')[[1]]
    



if('order'%in% names(opt)){
    opt$order = str_split(opt$order, ',')[[1]]
    
}

#load CNV files
scCNV=foreach(i=1:length(opt$scCNV),.packages = 'tidyverse',.combine = 'rbind')%do%{
    tmp=read_tsv(opt$scCNV[i],col_types = cols())
    if('order'%in% names(opt)){
        tmp%>%
            mutate(
                group=factor(group, levels=opt$order)
            ) 
        
    }else{
        tmp
    }
}

#create directory
if (str_extract(opt$out, '.$') != '/') {
    opt$out = paste0(opt$out, '/')
}


if(!opt$CNV_values %in% c('B','CNV','Log2','all')){
    opt$CNV_values='B'
    warning('unrecognized CNV value to plot. Binarized tracks restored')
}

#store chr info
chr=unique(scCNV$chr)


scCNV=scCNV%>%
    mutate(pos=paste0(chr,':',start,'-',end))%>%
    dplyr::select(
        'pos',
        'Cell',
        'PercentageReplication',
        'basename',
        'group',
        'data'=case_when(
            opt$CNV_values=='B' ~ 'Rep',
            opt$CNV_values=='CNV' ~ 'CN',
            opt$CNV_values=='Log2' ~ 'CN_bg'
        )
    )%>%
    mutate(data=as.numeric(data))%>%
    spread(pos,data)

#removue na colums
scCNV=scCNV[ , colSums(is.na(scCNV)) == 0]
mat=scCNV[,-c(1:4)]
scCNV=scCNV[,c(1:4)]

if(opt$CNV_values=='B') {
    if (opt$per_Chr) {
        # calculate distance per each chr
        cl = makeCluster(opt$cores)
        registerDoSNOW(cl)
        

        results = foreach(C = chr,.packages = 'ade4',.inorder = T) %dopar% {
            as.matrix(dist.binary(mat[, grepl(pattern = paste0(C, ':'), x = colnames(mat))], method = 2))
            
        }
            
        stopCluster(cl)
        
        names(results)=chr
        
    } else{
        results=list()
        results[['all Chr']]=as.matrix(dist.binary(mat, method = 2))
    }
    }else{
       
       if (opt$per_Chr) {
           # calculate distance per each chr
           cl = makeCluster(opt$cores)
           registerDoSNOW(cl)
           
           
           results = foreach(C = chr,.inorder = T) %dopar% {
               as.matrix(mat[, grepl(pattern = paste0(C, ':'), x = colnames(mat))])
           }
           
           stopCluster(cl)
           
           names(results)=chr
           
       } else{
           results=list()
           results[['all Chr']]=as.matrix(mat)
       }
        
   }

#free mem
    rm('mat')
    
    scCNV=foreach(C = names(results),.packages = c('Rtsne','tidyverse'),.combine = 'rbind') %do% {    
        Perplex = ceiling(ncol(results[[C]]) / 50)
    tsne <-
        Rtsne(
            X = results[[C]],
            dims = 2,
            perplexity = ifelse(Perplex < 10, 10, Perplex),
            check_duplicates = F,
            theta = 0.25,
            is_distance = opt$CNV_values=='B',
            verbose = F,
            max_iter = 5000,
            num_threads = opt$cores,
            partial_pca = T
        )

    scCNV%>%mutate(Chr=C,
                x=tsne$Y[,1],
                y=tsne$Y[,2])
    }
    
    write_tsv(scCNV,paste0(
        opt$out,
        '/',
        opt$output_file_base_name,
        '_tsne.txt'
    ))
X=foreach(C = names(results),.packages = c('Rtsne','tidyverse'),.combine = 'rbind') %do% {
    plot=scCNV%>%filter(Chr==C)%>%ggplot()+geom_point(aes(x,y,color=group,shape=group),alpha=0.5)+xlab('TSNE - 1')+ylab('TSNE - 2')+facet_wrap(~Chr)
    
    suppressMessages(ggsave(
        plot = plot,
        filename = paste0(
            opt$out,
            opt$output_file_base_name,
            '_',C,'_tsne_color_by_group.pdf'
        )
    ))
    
    plot=scCNV%>%filter(Chr==C)%>%ggplot()+geom_point(aes(x,y,color=basename,shape=group),alpha=0.5)+xlab('TSNE - 1')+ylab('TSNE - 2')+facet_wrap(~Chr)
    
    suppressMessages(ggsave(
        plot = plot,
        filename = paste0(
            opt$out,
            opt$output_file_base_name,
            '_',C,'_tsne_color_by_basename.pdf'
        )
    ))
    
    plot=scCNV%>%filter(Chr==C)%>%ggplot()+geom_point(aes(x,y,color=PercentageReplication,shape=group),alpha=0.5)+scale_color_gradient2(low = "#FFEA46FF", mid = "#7C7B78FF", high = "#00204DFF",lim=c(0,1),midpoint = 0.5)+xlab('TSNE - 1')+ylab('TSNE - 2')+facet_wrap(~Chr)
    suppressMessages(ggsave(
        plot = plot,
        filename = paste0(
            opt$out,
            opt$output_file_base_name,
            '_',C,'_tsne_color_by_rep_percentage.pdf'
        )
    ))
    C
}
print('done')
                 
