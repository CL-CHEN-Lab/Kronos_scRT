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
        help = "What type of date to plot for the single cell traks: ('B'=Binarized, 'CNV'=Copy number variation, 'log2'=log2(CNV_Cell/CNV_mean_G1/G2_cells)) [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("--per_Chr"),
        type = 'logical',
        default = F,
        action = 'store_TRUE',
        help = "Calculate TSNE/UMAP on each chromosome",
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
    ),
    make_option(
        c("-s", "--seed"),
        type = "integer",
        default = as.integer(Sys.Date()),
        help = "Set seed for reproducibility (optional).",
        metavar = "integer"
    ),
    make_option(
        c("-U", "--UMAP"),
        type = "logical",
        default = FALSE,
        action = "store_true",
        help = "Skip t-SNE, only plot UMAP.",
        metavar = "logical"
    ),
    make_option(
        c("-T", "--TSNE"),
        type = "logical",
        default = FALSE,
        action = "store_true",
        help = "Skip UMAP, only plot t-SNE.",
        metavar = "logical"
    ),
    make_option(
        c("--chr_prefix"),
        type = "character",
        action = 'store',
        help = "Chromosome prefix, if there is no prefix use none [default= %default]",
        default = "chr",
        metavar = "character"
    ),
    make_option(
        c("--chr_range"),
        type = "character",
        action = 'store',
        help = "Chromosomes to consider in the analysis (example 1:5,8,15:18,X) [default= %default]",
        default = "1:22",
        metavar = "character"
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
suppressPackageStartupMessages(library(umap, quietly = TRUE))
#output dir
system(paste0('mkdir -p ', opt$out))

#set plotting theme
theme_set(theme_bw())

#Set seed for reproducibility
set.seed(opt$seed)

#check files
if (!'scCNV' %in% names(opt)) {
    stop("scCNV file must be provided. See script usage (--help)")
}

opt$scCNV = str_split(opt$scCNV, ',')[[1]]

if('order'%in% names(opt)){
    opt$order = str_split(opt$order, ',')[[1]]
    
}

#load CNV files
scCNV=foreach(i=1:length(opt$scCNV),.packages = 'tidyverse',.combine = 'rbind')%do%{
    tmp=read_tsv(opt$scCNV[i],col_types = cols(chr='c'))
    if('order'%in% names(opt)){
        tmp%>%
            mutate(
                group=factor(group, levels=opt$order)
            ) 
        
    }else{
        tmp
    }
}

chr_list = paste0(ifelse(opt$chr_prefix=='none','',opt$chr_prefix), unlist(Convert_to_range(str_split(opt$chr_range,',')[[1]])))

#filter chr and convert into factor
scCNV=scCNV%>%
    filter(chr %in% chr_list)%>%
    mutate(chr = factor(x =  chr, levels = chr_list)) %>%
    drop_na()


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

#remove NA columns
scCNV=scCNV[,colSums(is.na(scCNV)) == 0]

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

scCNV_umap = scCNV
# TSNE
if (!opt$UMAP){
    scCNV=foreach(C = names(results),.packages = c('Rtsne','tidyverse'),.combine = 'rbind') %do% {    
        Perplex = ceiling(nrow(results[[C]]) / 50)
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
    X=foreach(C = names(results),.packages = c('tidyverse'),.combine = 'rbind') %do% {
        plot=scCNV%>%filter(Chr==C)%>%ggplot()+geom_point(aes(x,y,color=group,shape=group),alpha=0.4,size=2)+xlab('TSNE1')+ylab('TSNE2')+facet_wrap(~Chr)
        
        suppressMessages(ggsave(
            plot = plot, dpi = 300,
            filename = paste0(
                opt$out,
                opt$output_file_base_name,
                '_',C,'_tsne_color_by_group.pdf'
            )
        ))
        
        plot=scCNV%>%filter(Chr==C)%>%ggplot()+geom_point(aes(x,y,color=basename,shape=group),alpha=0.4,size=2)+xlab('TSNE1')+ylab('TSNE2')+facet_wrap(~Chr)
        
        suppressMessages(ggsave(
            plot = plot, dpi = 300,
            filename = paste0(
                opt$out,
                opt$output_file_base_name,
                '_',C,'_tsne_color_by_basename.pdf'
            )
        ))
        
        plot=scCNV%>%filter(Chr==C)%>%ggplot()+geom_point(aes(x,y,color=PercentageReplication,shape=group),alpha=0.4,size=2)+scale_color_gradient2(low = "#FFEA46FF", mid = "#7C7B78FF", high = "#00204DFF",lim=c(0,1),midpoint = 0.5)+xlab('TSNE1')+ylab('TSNE2')+facet_wrap(~Chr)
        suppressMessages(ggsave(
            plot = plot, dpi = 300,
            filename = paste0(
                opt$out,
                opt$output_file_base_name,
                '_',C,'_tsne_color_by_rep_percentage.pdf'
            )
        ))
        C
    }
}

# UMAP
if (!opt$TSNE){
    scCNV_umap=foreach(C = names(results),.packages = c('umap','tidyverse'),.combine = 'rbind') %do% {    
        if(opt$CNV_values=='B'){
            input_mat = 'dist'
        }else{
            input_mat = "data"
        }
        umap <- umap::umap(d = results[[C]], input = input_mat, random_state = opt$seed)
        scCNV_umap%>%mutate(Chr=C,
                            x=umap$layout[,1],
                            y=umap$layout[,2])
    }
    write_tsv(scCNV_umap,paste0(
        opt$out,
        '/',
        opt$output_file_base_name,
        '_umap.txt'
    ))
    X_umap=foreach(C = names(results),.packages = c('tidyverse'),.combine = 'rbind') %do% {
        plot=scCNV_umap%>%filter(Chr==C)%>%ggplot()+geom_point(aes(x,y,color=group,shape=group),alpha=0.4,size=2)+xlab('UMAP1')+ylab('UMAP2')+facet_wrap(~Chr)
        
        suppressMessages(ggsave(
            plot = plot, dpi = 300,
            filename = paste0(
                opt$out,
                opt$output_file_base_name,
                '_',C,'_umap_color_by_group.pdf'
            )
        ))
        
        plot=scCNV_umap%>%filter(Chr==C)%>%ggplot()+geom_point(aes(x,y,color=basename,shape=group),alpha=0.4,size=2)+xlab('UMAP1')+ylab('UMAP2')+facet_wrap(~Chr)
        
        suppressMessages(ggsave(
            plot = plot, dpi = 300,
            filename = paste0(
                opt$out,
                opt$output_file_base_name,
                '_',C,'_umap_color_by_basename.pdf'
            )
        ))
        
        plot=scCNV_umap%>%filter(Chr==C)%>%ggplot()+geom_point(aes(x,y,color=PercentageReplication,shape=group),alpha=0.4,size=2)+scale_color_gradient2(low = "#FFEA46FF", mid = "#7C7B78FF", high = "#00204DFF",lim=c(0,1),midpoint = 0.5)+xlab('UMAP1')+ylab('UMAP2')+facet_wrap(~Chr)
        suppressMessages(ggsave(
            plot = plot, dpi = 300,
            filename = paste0(
                opt$out,
                opt$output_file_base_name,
                '_',C,'_umap_color_by_rep_percentage.pdf'
            )
        ))
        C
    }
}

print('done')
