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
        default = 0.1,
        action = "store",
        help = "DeltaRT threshold to define changes",
        metavar = "double"
    ),
    make_option(
        c('-n', "--n_regions"),
        type = "integer",
        action = "store",
        default = 10,
        help = "Number of regions to plot",
        metavar = "integer"
    ),
    make_option(
        c("-r", "--region"),
        type = "character",
        default = NULL,
        help = "Region to plot  chr:start-end (multiple regins can be separated by a comma), it supports units",
        metavar = "character"
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


percentages = function(data,mapping,...){
    x <- pull(data,var =  as_label(mapping$x))
    y <- pull(data,var =  as_label(mapping$y))
    tmp=foreach (i = seq(0, 1, 0.01),
                       .combine = 'rbind',
                       .packages = 'tidyverse') %do% {
                           type = function(RT.x, RT.y, i) {
                               RT.x = ((RT.x - min(RT.x)) / (max(RT.x) - min(RT.x))) - 0.5
                               RT.y = ((RT.y - min(RT.y)) / (max(RT.y) - min(RT.y))) - 0.5
                               delta = RT.x - RT.y
                               result = case_when(
                                   (RT.x / RT.y < 0 & RT.y > 0 & abs(delta) > i) ~  'EtoL',
                                   (RT.x / RT.y < 0 & RT.y < 0 & abs(delta) > i) ~  'LtoE',
                                   (RT.x / RT.y >= 0 & delta < -i) ~ 'toLater',
                                   (RT.x / RT.y >= 0 & delta > i) ~ 'toEarlier',
                                   T ~ 'unchanged'
                               )
                               return(result)
                           }
                           
                           tibble(X=x,Y=y,l=length(x)) %>%
                               mutate(type = type(X, Y, i)) %>%
                               group_by(type) %>%
                               summarise(counts = n(),
                                         l=unique(l)) %>%
                               mutate(percent = round(counts / l,2)) %>%
                               dplyr::select(type, percent) %>%
                               mutate(th = i)%>%
                               filter(type != 'unchanged')
                           
                       }
    
    mapping=ggplot2:::rename_aes(modifyList(mapping,aes(x = th,y = percent,color=type)))
    p=ggplot(data=tmp,mapping=mapping,... = ...)+geom_line()+
        scale_y_continuous(labels = scales::percent_format())+coord_cartesian(xlim = c(0,1),ylim = c(0,1))+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        scale_color_manual(values = c(
            'EtoL'='#005095',
            'LtoE'='#dfbd31',
            'toLater'='#83007e',
            'toEarlier'='#78bd3e'
        ))
    
 return(p)
}

n_samples=length(names(data[-c(1:3)]))

p=ggpairs(
    data %>% dplyr::select(-chr, -start, -end),
    upper =  list(continuous = percentages),
    lower =  list(continuous = percentages),
    title = 'ΔRT=RTx-RTy',
    legend = c(1, 2),
    xlab = 'ΔRT threshold',
    ylab = '% ofchanging bins',
    diag = list(
        continuous = function(data, mapping)
            ggplot() + theme_void() + annotate('text',x = 1,y=1,label=as_label(mapping$x))+
            annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha=0.2)
    ),columnLabels=NULL
) +
    ylab('% ofchanging bins') + xlab('ΔRT threshold')


suppressMessages(
    ggsave(
    paste0(opt$out, '/changes_distribution.pdf'),
    plot = p,
    limitsize = FALSE,
    device = cairo_pdf,
    height = n_samples*1.2 ,
    width = n_samples*1.2 
)
)

# identify changing regions
type = function(RT.x, RT.y, i) {
    RT.x = ((RT.x - min(RT.x)) / (max(RT.x) - min(RT.x))) - 0.5
    RT.y = ((RT.y - min(RT.y)) / (max(RT.y) - min(RT.y))) - 0.5
    delta = RT.x - RT.y
    result = case_when(
        (RT.x / RT.y < 0 & RT.y > 0 & abs(delta) > i) ~  'EtoL',
        (RT.x / RT.y < 0 & RT.y < 0 & abs(delta) > i) ~  'LtoE',
        (RT.x / RT.y >= 0 & delta < -i) ~ 'toLater',
        (RT.x / RT.y >= 0 & delta > i) ~ 'toEarlier',
        T ~ 'unchanged'
    )
    return(result)
}


changing = foreach(x = colnames(data)[-c(1:3)],
                   .combine = function(x, y)
                       inner_join(
                           x,
                           y,
                           by = c("chr", "start", "end",  "ΔRT")
                       )) %:%
    foreach(
        y = colnames(data)[-c(1:3)],
        .combine = function(x, y)
            inner_join(
                x,
                y,
                by = c("chr", "start", "end", "ΔRT")
            )
    ) %do% {

        if(x!=y){
            data  %>%
                mutate(
                    type=type(RT.x = data[[x]],RT.y = data[[y]], i=opt$deltaRT_threshold),
                    ΔRT = opt$deltaRT_threshold
                ) %>%
                dplyr::select(chr,start,end,ΔRT,type)%>%
                rename(type=paste0('ΔRT=',x,'-',y))
            
           
        }else{
            data%>%
                mutate( ΔRT = opt$deltaRT_threshold)%>%
                dplyr::select(chr, start, end, ΔRT)
        }
        
    }
    
write_tsv(changing, paste0(opt$out, '/changing_regions.tsv'))

changing_info=changing%>%
    dplyr::select(-chr,-start,-end,-ΔRT)

n=(length(names(data[-c(1:3)]))-1)*2
changing_info = as.matrix(changing_info!='unchanged')

changing_info=changing[rowSums(changing_info)==n,]


specific_changes=foreach(i=names(data)[-c(1:3)],.combine = 'rbind')%do%{
    x=changing_info%>%
        dplyr::select(contains(i))
    x= as.matrix(x!='unchanged')
    x=changing_info[rowSums(x)==n,]
    
    x%>%
        dplyr::select(chr,start,end)%>%
        mutate(group=i)
}

specific_changes = specific_changes %>% makeGRangesFromDataFrame(
    seqnames.field = 'chr',
    start.field = 'start',
    end.field = 'end',
    keep.extra.columns = T
)

specific_changes=GenomicRanges::split(specific_changes,specific_changes$group)%>%
    GenomicRanges::reduce(min.gapwidth=5000000) %>%
    as_tibble() %>%
    dplyr::select('chr'=seqnames, start, end, 'group'=group_name) 

write_tsv(specific_changes, paste0(opt$out, '/cell_type_specific_changing_regions.tsv'))


changing = 
    changing%>%
    gather(couple,change,-chr,-start,-end,-ΔRT)%>%
    filter(change!='unchanged')%>%
    makeGRangesFromDataFrame(
    seqnames.field = 'chr',
    start.field = 'start',
    end.field = 'end',
    keep.extra.columns=T
)

changing=GenomicRanges::split(changing,changing$couple)%>%
        GenomicRanges::reduce(min.gapwidth=5000000) %>%
    as_tibble() %>%
    dplyr::select('chr'=seqnames, start, end) %>% 
    mutate(chr=as.character(chr))

# regions plot to feed to ggpairs
plot_rt=function(data,mapping,...){
    mapping1 = ggplot2:::rename_aes(modifyList(mapping,aes(x = pos)))
    mapping2 = mapping1
    mapping2$y= mapping$x
    mapping3=aes(xmin=start,xmax=end,ymin=0,ymax=1,fill=fill)
    data1=data%>%
        gather('type', 'pos', start, end)
    data$fill=abs(data[as_label(mapping$x)]-data[as_label(mapping$y)])[,]
    data=data%>%
        mutate(fill=ifelse(fill > ΔRT, paste('|ΔRT| >',ΔRT), paste('|ΔRT| <',ΔRT)),
               fill=factor(fill, levels = c( paste('|ΔRT| <',opt$deltaRT_threshold),paste('|ΔRT| >',opt$deltaRT_threshold),'Column RT','Row RT')))
    
    scales=c(NA,'orange','blue','red')
    names(scales)=c(paste('|ΔRT| <',opt$deltaRT_threshold),paste('|ΔRT| >',opt$deltaRT_threshold),'Column RT', 'Row RT')    
    p= ggplot() +
        geom_hline(yintercept = 0.5) +
        geom_rect(data=data,mapping = mapping3, alpha=0.5)+
        geom_line(data=data1,mapping = mapping1, color='red' ) +
        geom_line(data=data1,mapping = mapping2, color='blue' )+ 
        scale_x_continuous(
            labels = trans_format(format = 'Mb',
                                  trans = function(x) paste((x/1000000),'Mb',sep='')))+ 
        theme(legend.title=element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))+
        scale_fill_manual(values = scales,drop=F)
    return(p)
}



if (!'region' %in% names(opt)) {
   
     #limits for the extention of the plotting 
    chr_limits=data%>%group_by(chr)%>%summarise(start_min=min(start),end_max=max(end))%>%ungroup()
    
    for (i in sample(1:length(changing$chr), opt$n_regions)) {
        
        to_plot =changing[i,]
        
        if(to_plot$end-to_plot$start != 50000000){
            add=50000000-(to_plot$end-to_plot$start)
            to_plot=to_plot%>%
                mutate(start=start-add/2,
                       end=end+add/2)%>%
                inner_join(chr_limits, by = "chr")%>%
                mutate(start=case_when(
                    start < start_min ~ start_min,
                    T~start
                ),
                end=case_when(
                    end > end_max ~ end_max,
                    T~end
                ))
        }
        
        # prepare name file 
        name_reg=min( str_length(str_extract(to_plot$start,'0{1,10}$')), str_length(str_extract(to_plot$end,'0{1,10}$')))
        name_reg= paste0(
            to_plot$chr,"_",
            case_when(
                is.na(name_reg) ~ paste0(to_plot$start,'bp_',to_plot$end,'bp'),
                name_reg < 3 ~ paste0(to_plot$start,'bp_',to_plot$end,'bp'),
                name_reg < 6 ~ paste0(to_plot$start/10^3,'Kb_',to_plot$end/10^3,'Kb'),
                name_reg >= 6 ~ paste0(to_plot$start/10^6,'Mb_',to_plot$end/10^6,'Mb')
            ) )
        
        p = data %>%
            mutate(end = end - 1,
                   ΔRT=opt$deltaRT_threshold) %>%
            filter(
                chr == to_plot$chr,
                start >= to_plot$start ,
                end <= to_plot$end 
            )
        
        p=p %>%
            ggpairs(upper = NULL,lower  = list( continuous=plot_rt),columns = names(p)[!names(p) %in% c('chr','start','end','ΔRT')],
                    title = name_reg,
                    legend = c(2,1),
                    xlab = '',
                    ylab = 'RT',
                    diag = list(
                        continuous = function(data, mapping)
                            ggplot() + theme_void()+
                            geom_polygon(data=tibble(x=c(-1,-1,1),y=c(1,-1,1)), aes(x = x, y = y),fill='#a7001b',alpha=0.5)+
                            geom_polygon(data=tibble(x=-c(-1,-1,1),y=-c(1,-1,1)), aes(x = x, y = y),fill='#005095',alpha=0.5) + 
                            annotate('text',x = 0,y=0,label=paste('bold(',as_label(mapping$x),')'),parse=T)
                    ),columnLabels=NULL )
        p
        suppressMessages(
        ggsave(p, filename = paste0(opt$out, '/changing_region_', name_reg, '.pdf'),
               device = cairo_pdf,
               height = n_samples*1.2 ,
               width = n_samples*1.2 )
        )
    }
} else{
    #load bed file if exist
    if(file.exists(opt$region)){
    opt$region =read_tsv(opt$region,col_names = c('chr','start','end') )%>%
        mutate(
            n_0_start=str_length(str_extract(start,'0{1,10}$')),
            n_0_end=str_length(str_extract(end,'0{1,10}$')),
            unit=min(factor(case_when(
                is.na(n_0_start) ~'bp',
                n_0_start < 3 ~ 'bp',
                n_0_start < 6 ~'Kb',
                n_0_start >= 6 ~ 'Mp'),levels = c('bp','Kb','Mp'), ordered=TRUE),
                factor(case_when(
                    is.na(n_0_end) ~'bp',
                    n_0_end < 3 ~ 'bp',
                    n_0_end < 6 ~'Kb',
                    n_0_end >= 6 ~ 'Mp'),levels = c('bp','Kb','Mp'), ordered=TRUE
                )),
        name_reg=paste(
            chr,
            case_when(
                unit == 'bp' ~ paste0(start,unit,'_',end,unit),
                unit == 'Kb' ~ paste0(start/10^3,unit,'_',end/10^3,unit),
                unit == 'Mb' ~ paste0(start/10^6,unit,'_',end/10^6,unit)
            ),sep = '_' ))%>%
        dplyr::select(-unit,-n_0_start,-n_0_end)
    }else{
    #reshape provided regions
    opt$region = tibble(coord = str_split(opt$region, pattern = ',')[[1]]) %>%
        mutate(name_reg=str_replace_all(coord,pattern = '[-:]',replacement = '_'))%>%
        separate(coord, c('chr', 'pos'), ':') %>%
        separate(pos, c('start', 'end'), '-')%>%
        mutate(
            start_unit=str_extract(start,pattern = '.{2}$'),
            start=as.numeric(str_remove(start, "[Bb][Pp]|[Kk][Bb]|[Mm][Bb]")) * case_when(
            grepl(x =start_unit,pattern =  '[Kk][Bb]') ~ 1000,
            grepl(x =start_unit, pattern = '[Mm][Bb]') ~ 1000000,
            grepl(x =start_unit, pattern = '[Bp][Pp]') ~ 1,
            grepl(x = start_unit,pattern =  '[0-9][0-9]') ~ 1
        ),            end_unit=str_extract(end,pattern = '.{2}$'),

               end=as.numeric(str_remove(end, "[Bb][Pp]|[Kk][Bb]|[Mm][Bb]")) * case_when(
                   grepl(x =end_unit,pattern =  '[Kk][Bb]') ~ 1000,
                   grepl(x =end_unit, pattern = '[Mm][Bb]') ~ 1000000,
                   grepl(x =end_unit, pattern = '[Bp][Pp]') ~ 1,
                   grepl(x = end_unit,pattern =  '[0-9][0-9]') ~ 1
               ))%>%
        dplyr::select(-start_unit,-end_unit)
    }
    
    for (i in 1:length(opt$region$chr)) {
        p = data %>%
            mutate(end = end - 1) %>%
            filter(chr == opt$region$chr[i],
                   start >= opt$region$start[i],
                   end <= opt$region$end[i]) 
        
        p=p %>%
            ggpairs(upper = NULL,lower  = list( continuous=plot_rt),columns = names(p)[!names(p) %in% c('chr','start','end','ΔRT')],
                    title = opt$region$name_reg[i],
                    xlab = '',
                    ylab = 'RT',
                    diag = list(
                        continuous = function(data, mapping)
                            ggplot() + theme_void()+
                            geom_polygon(data=tibble(x=c(-1,-1,1),y=c(1,-1,1)), aes(x = x, y = y),fill='#a7001b',alpha=0.5)+
                            geom_polygon(data=tibble(x=-c(-1,-1,1),y=-c(1,-1,1)), aes(x = x, y = y),fill='#005095',alpha=0.5) + 
                            annotate('text',x = 0,y=0,label=paste('bold(',as_label(mapping$x),')'),parse=T)
                    ),columnLabels=NULL )
        
        suppressMessages(ggsave(p,
               filename = paste0(
                   opt$out,
                   '/changing_region_',
                   opt$region$name_reg[i],
                   '.pdf'
               ),
               device = cairo_pdf,
               height = n_samples*1.2 ,
               width = n_samples*1.2 ))
        
    }
    
}

                          
print('done')
