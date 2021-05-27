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

theme_set(theme_bw())


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
    colors =  c( '#BCAF6FFF', '#7C7B78FF','#00204DFF')
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
    plot = ggpairs(scRT,diag = list(continuous =function(data, mapping, ...){
        p <- ggplot(data,mapping)+
            geom_density(aes(y=..density../max(..density..)),
                         fill='black')+
            scale_x_continuous(breaks = c(0,0.5,1))+
            scale_y_continuous(breaks = c(0,0.5,1))
        return(p)
    }),
    upper = list(continuous =function(data, mapping, ...){
        
        data=tibble(
            xmin=-Inf,
            xmax=Inf,
            ymin=-Inf,
            ymax=Inf,
            Corr=cor(data[as_label(mapping$x)],data[as_label(mapping$y)])
        )
        
        p <- ggplot(data,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=Corr)) + 
            geom_rect()+
            annotate('text',0.5,0.5,label=paste("Corr:",round(data$Corr,3),sep = '\n'),color='white')+
            scale_fill_gradient2(low = '#BCAF6FFF',high = '#00204DFF',mid = '#7C7B78FF',midpoint = 0,limits=c(-1,1))+
            coord_cartesian(xlim = c(0,1),ylim = c(0,1))+
            scale_x_continuous(breaks = c(0,0.5,1))+
            scale_y_continuous(breaks = c(0,0.5,1))
        
        return(p)
    }),
    lower = list(continuous =function(data, mapping, ...){
        p <- ggplot(data = data, mapping = mapping) + 
            geom_hex(bins=50,aes(fill=..ndensity..))+
            scale_fill_gradientn('Density',colours =c("#FFEA46FF","#D3C164FF","#A69D75FF","#7C7B78FF","#575C6DFF","#233E6CFF","#00204DFF"))+
            coord_cartesian(xlim = c(0,1),ylim = c(0,1))+
            scale_x_continuous(breaks = c(0,0.5,1))+
            scale_y_continuous(breaks = c(0,0.5,1))+
            geom_abline(slope = 1,color='black',alpha=0.5)
        
        return(p)
    }),legend = c(2,1))+ theme(legend.position = "right",
                               axis.text.x = element_text(angle = 45,hjust = 1)),
    filename =paste0(
        opt$out,
        '/',
        opt$output_file_base_name,
        'pair_scatter_plot_RTs.pdf'
    )
))

print('done')

