#!/usr/local/bin/Rscript

# this script is meant to select the treshold to select cycling cells
try(library(tidyverse))
try(library(optparse))
try(library(doSNOW))
try(library(chunked))
try(library(gplots))
try(library(RColorBrewer))
try(library(foreach))

options(stringsAsFactors = FALSE)
options(warn=1) 
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
        c("-R", "--referenceRT"),
        type = "character",
        default = NULL,
        help = "Reference RT min=Late, max=Early, only one reference is allowed",
        metavar = "character"
    ),
    make_option(
        c("-C", "--chrSizes"),
        type = "character",
        default = NULL,
        help = "chromosome size file",
        metavar = "character"
    ),
    make_option(
        c("-r", "--region"),
        type = "character",
        default = NULL,
        help = "Region to plot  chr:start-end (multiple regins can be separated by a comma)",
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
        c("-b", "--base_name"),
        type = "character",
        default = "exp",
        help = "base name for files names [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("-S", "--threshold_Sphase"),
        type = "character",
        help = "Threshold to identify S-phase cells",
        metavar = "character"
    ),
    make_option(
        c("-G", "--threshold_G1G2phase"),
        type = "character",
        help = "Threshold to identify G1-phase cells. -S has to be selected and has to be bigger than -G",
        metavar = "character"
    ),
    make_option(
        c("-M", "--threshold_meanploidy"),
        type = "character",
        default = 1.5,
        help = "Threshold to discard cells with ploidy lower and higher than n times the median ploidy of the G1 phase. [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("-B", "--binsSize"),
        type = "integer",
        default = 500000,
        help = "RT resolution [default= %default] ",
        metavar = "integer"
    ),
    make_option(
        c("-k", "--keepXY"),
        type = "logical",
        default = TRUE,
        action = "store_false",
        help = "base name for files names",
        metavar = "logical"
    ),
    make_option(
        c("-c", "--cores"),
        type = "integer",
        default = 3,
        help = "numbers of parallel jobs to run [default= %default] ",
        metavar = "integer"
    ),
    make_option(
        c("-p", "--plot"),
        type = "logical",
        default = F,
        action = "store_true",
        help = "if selected prints some rondome regins, if -r is selected those regins are use to print RT [default= %default] ",
        metavar = "logical"
    ),
    make_option(
        c("--min"),
        type = "integer",
        default = 3,
        help = "minimum number of cells per bin [default= %default] ",
        metavar = "integer"
    ),
    make_option(
        c("--Var_against_reference"),
        type = "logical",
        default = F,action = "store_true",
        help = "Variability metrics are calculated usign reference RT in addiction to the calculated one [default= %default] ",
        metavar = "logical"
    )
    
)

opt = parse_args(object = OptionParser(option_list = option_list))

#check inputs

if (!'file' %in% names(opt)) {
    stop("Per cell stat file must be provided. See script usage (--help)")
}

if (!'tracks' %in% names(opt)) {
    stop("File containing cells CNV must be provided. See script usage (--help)")
}

if (!'chrSizes' %in% names(opt)) {
    stop("File containing Chromosomes sizes must be provided. See script usage (--help)")
}

if (opt$Var_against_reference) {
    if (!'referenceRT' %in% names(opt)) {
        warning("Reference genome not provided")
        opt$Var_against_reference=F
    }
}

#create directory
system(paste0('mkdir -p ./', opt$out))

#load files
opt$file=str_split(opt$file,',',simplify = F)[[1]]
opt$tracks=str_split(opt$tracks,',',simplify = F)[[1]]
opt$base_name=str_split(opt$base_name,',',simplify = F)[[1]]

if (length(opt$tracks) != length(opt$file)){
    stop("The number stat files does not match provided trakcs. See script usage (--help)")
}
if (length(opt$base_name) != length(opt$file)){
    index=rep(1:length(opt$file),length(opt$file),each=length(opt$base_name),length.out=length(opt$file))
    opt$base_name=paste(rep_len(opt$base_name,length(opt$file)),index,sep = '-')
    warning('basenames will be cyclicly recicled')
}

if ('threshold_Sphase' %in% names(opt)) {
    opt$threshold_Sphase=str_split(opt$threshold_Sphase,',',simplify = F)[[1]]
    if (length(opt$threshold_Sphase)<length(opt$tracks)){
        opt$threshold_Sphase=rep_len(opt$threshold_Sphase,length(opt$tracks))
        warning('Sphase thresholds will be cyclicly recicled for all the samples')
        
    }
    opt$threshold_Sphase=as.numeric(opt$threshold_Sphase)
}
if ('threshold_G1G2phase' %in% names(opt)) {
    opt$threshold_G1G2phase=str_split(opt$threshold_G1G2phase,',',simplify = F)[[1]]
    if (length(opt$threshold_G1G2phase)<length(opt$tracks)){
        opt$threshold_G1G2phase=rep_len(opt$threshold_G1G2phase,length(opt$tracks))
        warning('G1G2phase thresholds will be cyclicly recicled for all the samples')
        
    }
    opt$threshold_G1G2phase=as.numeric(opt$threshold_G1G2phase)
}

if ('threshold_meanploidy' %in% names(opt)) {
    opt$threshold_meanploidy=str_split(opt$threshold_meanploidy,',',simplify = F)[[1]]
    if (length(opt$threshold_meanploidy)<length(opt$tracks)){
        opt$threshold_meanploidy=rep_len(opt$threshold_meanploidy,length(opt$tracks))
        warning('meanploidy thresholds will be cyclicly recicled for all the samples')
    }
    opt$threshold_meanploidy=tibble(threshold_meanploidy=as.numeric(opt$threshold_meanploidy),
                                        basename=opt$base_name)
}

data <- foreach(i=1:length(opt$file),.combine = 'rbind',.packages = 'tidyverse',.verbose = T)%do%{
    file=read_csv(opt$file[i])%>%
        mutate(basename=opt$base_name[i])
    if ('threshold_Sphase' %in% names(opt)) {
        
    file=file%>%
        mutate(threshold_Sphase=opt$threshold_Sphase[i])
    }
    if ('threshold_G1G2phase' %in% names(opt)) {
        
        file=file%>%
            mutate(threshold_G1G2phase=opt$threshold_G1G2phase[i])
    }
    file
}

all_tracks <-foreach(i=1:length(opt$tracks),.combine = 'rbind',.packages = 'tidyverse')%do%{
    read_delim(opt$tracks[i], delim = '\t', skip = 2)%>%
        mutate(basename=opt$base_name[i])
} 

if(opt$keepXY) {
    Chr_Size <-
        read_delim(opt$chrSizes, delim = '\t', col_names = c('chr','size')) %>%
        filter(!chr %in% c('chrX', 'chrY'))
} else{
    Chr_Size <- read_delim(opt$chrSizes, delim = '\t', col_names = c('chr','size'))
}
if ('referenceRT' %in% names(opt)) {
    Reference_RT <- read_delim(opt$referenceRT, delim = '\t', col_names = c('chr','start','end','RT'))%>%
        mutate(RT=(RT-min(RT))/(max(RT)-min(RT)))
}
# calculate the new treshold


if ('threshold_Sphase' %in% names(opt)) {
    data = data %>%
        mutate(
            is_high_dimapd = ifelse(normalized_dimapd > threshold_Sphase, T, F),
            is_noisy = ifelse(is_high_dimapd, T, is_noisy)
        )
}
if ('threshold_G1G2phase' %in% names(opt)) {
    if (any(data$threshold_G1G2phase > data$threshold_Sphase) ){
        stop("G1/G2 phase threshold has to be smaller than S phase threshold. See script usage (--help)")
    }
    median_ploidy_G1_G2_cells = data%>%
        filter(is_noisy == F & normalized_dimapd < threshold_G1G2phase)%>%
        group_by(basename)%>%
        summarise(median_ploidy_G1_G2_cells=median(mean_ploidy))

    data = data %>%
        inner_join(opt$threshold_meanploidy, by='basename')%>%
        inner_join(median_ploidy_G1_G2_cells, by='basename')%>%
        filter(
            mean_ploidy > median_ploidy_G1_G2_cells / threshold_meanploidy ,
            mean_ploidy < median_ploidy_G1_G2_cells * threshold_meanploidy,
            !ploidy_confidence <= 2
        ) %>%
        mutate(Type = ifelse(
            as.logical(is_high_dimapd) == T &
                as.logical(is_noisy) == T,
            'S-phase',
            ifelse(
                as.logical(is_high_dimapd) == F &
                    as.logical(is_noisy) == F &
                    normalized_dimapd < threshold_G1G2phase,
                'G1/G2 cells',
                'unknown cells'
            )
        ))
} else{
    
    median_ploidy_G1_G2_cells = data%>%
        filter(is_noisy == F )%>%
        group_by(basename)%>%
        summarise(median_ploidy_G1_G2_cells=median(mean_ploidy))
    
    data  = data %>%
        inner_join(median_ploidy_G1_G2_cells)%>%
        filter(
            mean_ploidy > median_ploidy_G1_G2_cells / 1.4 ,
            mean_ploidy < median_ploidy_G1_G2_cells * 1.4,
            !ploidy_confidence <= 2
        ) %>%
        mutate(Type = ifelse(
            as.logical(is_high_dimapd) == T &
                as.logical(is_noisy) == T,
            'S-phase',
            ifelse(
                as.logical(is_high_dimapd) == F &
                    as.logical(is_noisy) == F,
                'G1/G2 cells',
                'unknown cells'
            )
        ))
    
}




p = data %>%
    ggplot(aes(mean_ploidy, normalized_dimapd, color = Type, shape=basename)) +
    geom_point(alpha = 0.3) +
    scale_color_manual(
        values = c(
            'G1/G2 cells' = 'darkred',
            'S-phase' = 'darkgreen',
            'unknown cells' = 'darkorange'
        )
    ) +
    theme(legend.position = 'top', legend.title = element_blank()) +
    geom_vline(aes(xintercept = median_ploidy_G1_G2_cells))+facet_wrap(~basename)

ggsave(p, filename = paste0(opt$out, '/', paste(opt$base_name,collapse = '_'), '_plot.pdf'))


# correct mean ploidy late S phase

data = data %>%
    group_by(basename)%>%
    mutate(
        to_multiply_to_ploidy = median_ploidy_G1_G2_cells / min(mean_ploidy),
        to_add_to_ploidy = max(mean_ploidy)-median_ploidy_G1_G2_cells,
        mean_ploidy_corrected = ifelse(
            as.logical(is_noisy) == T &
                mean_ploidy < median_ploidy_G1_G2_cells,
            to_multiply_to_ploidy * mean_ploidy + to_add_to_ploidy,
            mean_ploidy
        )
    )

p = data %>%
    ggplot(aes(mean_ploidy_corrected, normalized_dimapd, color = Type)) +
    geom_point(alpha = 0.3) +
    scale_color_manual(
        values = c(
            'G1/G2 cells' = 'darkred',
            'S-phase' = 'darkgreen',
            'unknown cells' = 'darkorange'
        )
    ) +
    theme(legend.position = 'top', legend.title = element_blank()) +
    geom_vline(aes(xintercept = median_ploidy_G1_G2_cells))+facet_wrap(~basename)


ggsave(p,filename = paste0(opt$out, '/', paste(opt$base_name,collapse = '_'), '_plot_sphase_corrected.pdf'))

#resolution
cl <- makeCluster(opt$cores)
registerDoSNOW(cl)
resolution = opt$binsSize

bins = foreach(chr = 1:length(Chr_Size$chr),
               .combine = 'rbind') %dopar% {
                   bins = seq(from = resolution,
                              to = Chr_Size$size[chr] ,
                              by =  resolution)
                   bins = data.frame(chr = as.character(Chr_Size$chr[chr]),
                                     start = bins - resolution,
                                     end = bins)
                   bins
               }

#select Sphase cells
selected_data = data %>%
    filter(Type == 'S-phase') %>%
    arrange(mean_ploidy_corrected) %>%
    mutate(index = 1:n()) %>%
    select(index, barcode, cell_id, mean_ploidy, mean_ploidy_corrected,basename,to_add_to_ploidy,to_multiply_to_ploidy)

# select G1/G2 cells
G1_G2_cells = data %>%
    filter(Type=='G1/G2 cells')


# select tracks of G1/G2 cells
G1_G2_cells_tracks = all_tracks %>%
    inner_join(G1_G2_cells, by = c('id' = 'cell_id','basename'))

#select Sphase traks
Sphase_tracks = all_tracks %>%
    inner_join(selected_data, by = c('id' = 'cell_id','basename')) %>%
    mutate(
        copy_number_corrected = ifelse(
            mean_ploidy == mean_ploidy_corrected,
            copy_number,
            to_multiply_to_ploidy * copy_number + to_add_to_ploidy
        )
    )

#calculate median CNV across a bin for S cells
signal_smoothed = foreach(
    chr = unique(bins$chr),
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach'),
    .verbose = T
) %dopar% {
    bins_in_chr = bins[bins$chr == chr,]
    bins_chr = foreach(
        bin = 1:length(bins_in_chr$chr),
        .combine = 'rbind',
        .packages = 'tidyverse'
    ) %do% {
        track_sign = Sphase_tracks %>%
            filter(
                `#chrom` == chr,
                (start >= bins_in_chr$start[bin] &
                     end <= bins_in_chr$end[bin]) |
                    (start <= bins_in_chr$start[bin] &
                         end >= bins_in_chr$start[bin]) |
                    (start <= bins_in_chr$end[bin] &
                         end >= bins_in_chr$end[bin])
            ) %>%
            group_by(index,basename) %>%
            summarise(
                chr = bins_in_chr$chr[bin],
                start = bins_in_chr$start[bin],
                end = bins_in_chr$end[bin],
                CN = median(copy_number_corrected)
            )
        
        track_sign
        
    }
    bins_chr%>%
        drop_na()
}

write_delim(
    x = signal_smoothed,
    path = paste0(opt$out, '/tmp.tsv'),
    delim = '\t',
    col_names = T
)
rm('signal_smoothed')

#calculate median CNV across a bin across all the G1/G2 cells
backgroud_smoothed = foreach(
    chr = unique(bins$chr),
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach'),
    .verbose = T
) %dopar% {
    bins_in_chr = bins[bins$chr == chr, ]
    bins_chr = foreach(
        bin = 1:length(bins_in_chr$chr),
        .combine = 'rbind',
        .packages = 'tidyverse'
    ) %do% {
        track_back = G1_G2_cells_tracks %>%
            filter(
                `#chrom` == chr,
                (start >= bins_in_chr$start[bin] &
                     end <= bins_in_chr$end[bin]) |
                    (start <= bins_in_chr$start[bin] &
                         end >= bins_in_chr$start[bin]) |
                    (start <= bins_in_chr$end[bin] &
                         end >= bins_in_chr$end[bin])
            ) %>%
            group_by(id,basename) %>%
            summarise(
                chr = bins_in_chr$chr[bin],
                start = bins_in_chr$start[bin],
                end = bins_in_chr$end[bin],
                background = median(copy_number)
            )
        track_back
    }
    bins_chr%>%
        drop_na()
}

backgroud_smoothed=backgroud_smoothed%>%
    group_by(chr,start,end,basename)%>%
    summarise(background=median(background))

write_delim(
    x = backgroud_smoothed,
    path = paste0(opt$out, '/tmp_bg.tsv'),
    delim = '\t',
    col_names = T
)

if ('referenceRT' %in% names(opt)) {
    # rebin reference RT
    Reference_RT = foreach(
        chr = unique(bins$chr),
        .combine = 'rbind',
        .packages = c('dplyr', 'foreach'),
        .verbose = T
    ) %dopar% {
        bins_in_chr = bins[bins$chr == chr, ]
        bins_chr = foreach(
            bin = 1:length(bins_in_chr$chr),
            .combine = 'rbind',
            .packages = 'dplyr'
        ) %do% {
            track_back = Reference_RT %>%
                filter(
                    chr == chr,
                    (start >= bins_in_chr$start[bin] &
                         end <= bins_in_chr$end[bin]) |
                        (start <= bins_in_chr$start[bin] &
                             end >= bins_in_chr$start[bin]) |
                        (start <= bins_in_chr$end[bin] &
                             end >= bins_in_chr$end[bin])
                ) %>%
                summarise(
                    chr = bins_in_chr$chr[bin],
                    start = bins_in_chr$start[bin],
                    end = bins_in_chr$end[bin],
                    RT = median(RT, na.rm = T)
                )
            track_back %>%
                filter(!is.na(RT))
        }
        bins_chr
    }
    #write output
    
    write_delim(
        x = Reference_RT,
        path = paste0(
            opt$out,
            '/',
            paste(opt$base_name,collapse = '_'),
            '_reference_replication_timing_',
            opt$binsSize,
            'bp.tsv'
        ),
        delim = '\t',
        col_names = T
    )
}
stopCluster(cl)



signal_smoothed = read_table_chunkwise(paste0(opt$out, '/tmp.tsv'),
                                       chunk_size = 100000,
                                       format = 'table') %>%
    inner_join(backgroud_smoothed, by = c("chr", "start", "end","basename")) %>%
    mutate(CN_bg = CN / background)

signal_smoothed = collect(signal_smoothed)

signal_smoothed=signal_smoothed%>%
    drop_na()%>%
    filter(is.finite(CN_bg))

# remouve control track
rm('backgroud_smoothed')


# exclude extreme values
quantile = quantile(signal_smoothed$CN_bg, c(0.025, 0.975))

#identify range within looking for a CNV treshold to define replicated and not replicated values
range = seq(quantile[[1]], quantile[[2]], (quantile[[2]] - quantile[[1]]) /100)

# identify threshold that minimazes the difference of the real data with a binary state (1 or 2)
selecte_th = foreach(i = range,
                     .combine = 'rbind',
                     .packages = 'tidyverse', .verbose = T) %do% {
                         summary = signal_smoothed %>%
                             mutate(Rep = ifelse(CN_bg >= i, 2, 1),
                                    Error = (Rep - CN_bg) ^ 2) %>%
                                 group_by(index,basename)  %>%
                             summarise(summary = sum(Error))

                             data.frame(
                                 th = i,
                                 basename=summary$basename,
                                 index = summary$index,
                                 sum_error = summary$summary
                             )
                        
                         
                     }

    selecte_th = selecte_th %>%
        group_by(index,basename) %>%
        filter(sum_error == min(sum_error)) %>%
        summarise(th = min(th))
    
    # mark replicated bins
    signal_smoothed = signal_smoothed %>%
        inner_join(selecte_th) %>%
        mutate(Rep = ifelse(CN_bg >= th, T, F))


#identify new distribution in the S phase based the ammount of replicated bins
new_index_list = signal_smoothed %>%
    group_by(index,basename) %>%
    summarise(perc_replication = mean(Rep)) %>%
    ungroup() %>%
    arrange(perc_replication) %>%
    group_by(basename)%>%
    mutate(newIndex = 1:n()) %>%
    select(index, newIndex,basename)

signal_smoothed = signal_smoothed %>%
    inner_join(new_index_list, by = c('index','basename'))

# bin cells in order to have eqully spaced bins with at least 3 cells in each bin
mean_cn = signal_smoothed %>%
    group_by(index,basename) %>%
    summarise(mean_CN = mean(Rep))
n = length(mean_cn$index)

while (T) {
    breaks = seq(min(mean_cn$mean_CN),
                 max(mean_cn$mean_CN),
                 (max(mean_cn$mean_CN) - min(mean_cn$mean_CN)) / n)
    
    counts =foreach(i=unique(mean_cn$basename),.combine = 'c')%do%{
        
        hist(x = mean_cn$mean_CN[mean_cn$basename==i],
                  breaks = breaks,
                  plot = F)$counts
        
    }
    
    if (!any(counts <= opt$min)) {
        break
    } else{
        n = round(n / 1.1)
    }
    if (n < 5 ){
        warning('it was not possible to equally distribute all the cells of all the samples in the same bins')
        break
    }
    
}

mean_cn$group = NA
group=0
for (i in breaks) {
    group=group+1
    if (!which(breaks == i) == length(breaks)) {
        mean_cn$group[mean_cn$mean_CN >= i &
                          mean_cn$mean_CN <= breaks[which(breaks == i) + 1]] = group
    }
}

signal_smoothed=signal_smoothed %>%
    inner_join(mean_cn, by = c('index','basename'))

write_delim(
    x = signal_smoothed,
    path = paste0(
        opt$out,
        '/',
        paste(opt$base_name,collapse = '_'),
        '_single_cells_CNV_',
        opt$binsSize,
        'bp.tsv'
    ),
    delim = '\t',
    col_names = T
)

plot=signal_smoothed%>%
    group_by(index,basename)%>%
    summarise(Rep_percentage=mean(Rep))%>%
    ggplot(aes(Rep_percentage,color=basename))+
    geom_density()+
    geom_vline(xintercept = breaks, col='grey')+
    scale_x_continuous(labels = scales::percent)

ggsave(
    plot,
    filename = paste0(
        opt$out,
        '/',
        paste(opt$base_name,collapse = '_'),
        'percentage_of_replicating_cells_and_binning.pdf'
    )
)


#calculate replication timing normalizing each bin by the number of cells in each bin and then calculating the average of the average
s50 = signal_smoothed %>%
    group_by(chr, start, end, group,basename) %>%
    summarise(Rep = mean(Rep)) %>%
    ungroup() %>%
    group_by(chr, start, end,basename) %>%
    summarise(RT = mean(Rep))

write_delim(
    x = s50,
    path = paste0(
        opt$out,
        '/',
        paste(opt$base_name,collapse = '_'),
        '_calculated_replication_timing_',
        opt$binsSize,
        'bp.tsv'
    ),
    delim = '\t',
    col_names = T
)




#plots
if (opt$plot) {
    system(paste0('mkdir ', opt$out, '/regions'))
    if (!'region' %in% names(opt)) {
        for (i in 1:length(Chr_Size$chr)) {
            region = round(runif(1, min = 1000000, max = Chr_Size$size[i] - 8000000),
                           0)
            Chr = Chr_Size$chr[i]
            Start = region
            End = region + 60000000
            
            trac_toplot = signal_smoothed %>%
                filter(
                    chr %in% Chr,
                    (start >= Start & end <= End) |
                        (start <= Start & end >= Start) |
                        (start <= End & end >= End)
                ) %>%
                mutate(
                    start = ifelse(start < Start, Start, start),
                    end = ifelse(end > End , End, end)
                )
            
            
            
            s50_toplot = s50 %>%
                ungroup() %>%
                filter(
                    chr %in% Chr,
                    (start >= Start & end <= End) |
                        (start <= Start & end >= Start) |
                        (start <= End & end >= End)
                ) %>%
                mutate(
                    RT = RT + 1,
                    start = ifelse(start < Start, Start, start),
                    end = ifelse(end > End , End, end)
                )
            
            plot = trac_toplot %>%
                ggplot() +
                geom_rect(
                    aes(
                        xmin = start,
                        xmax = end,
                        ymin = -newIndex,
                        ymax = -newIndex - 1,
                        fill = Rep + 1
                    )
                ) +geom_rect(
                    data = s50_toplot,
                    aes(
                        xmin = start,
                        xmax = end,
                        ymin = 11,
                        ymax = 20,
                        fill = RT
                    ),
                    inherit.aes = F
                ) +
                facet_grid( ~ chr, scale = 'free') +
                scale_fill_gradient(low = 'red',
                                    high = 'green',
                                    limits = c(1, 2)) +
                scale_x_continuous(
                    labels = function(x)
                        paste(x / 1000000, 'Mb', sep = ' ')
                ) +
                scale_y_discrete(
                    limits = c(15, 5, -seq(
                        10, max(trac_toplot$index), round((max(
                            trac_toplot$index
                        ) - 10) / 20)
                    )),
                    labels = c('RT', 'RT reference', seq(
                        10, max(trac_toplot$index), round((max(
                            trac_toplot$index
                        ) - 10) / 20)
                    ))
                ) +
                labs(y = 'S phase progression', fill = "CNV Binary\n") + theme(legend.position = 'top')+
                facet_wrap(~basename)
            
            if ('referenceRT' %in% names(opt)) {
                RT_toplot = Reference_RT %>%
                    filter(chr %in% Chr,
                           start >= Start,
                           end <= End) %>%
                    mutate(
                        RT = 1 + RT / max(RT, na.rm = T),
                        start = ifelse(start < Start, Start, start),
                        end = ifelse(end > End , End, end)
                    )
                plot = plot +
                    geom_rect(
                        data = RT_toplot,
                        aes(
                            xmin = start,
                            xmax = end,
                            ymin = 0,
                            ymax = 10,
                            fill = RT
                        ),
                        inherit.aes = F
                    )
            }
            
            
            ggsave(
                plot,
                filename = paste0(
                    opt$out,
                    '/regions/',
                    paste(opt$base_name,collapse = '_'),
                    '_plot_RT_',
                    Chr,
                    ':',
                    Start,
                    '-',
                    End,
                    '.pdf'
                )
            )
            
        }
        
    } else{
        #reshape regins
        opt$region = data.frame(coord = str_split(opt$region, pattern = ',')[[1]]) %>%
            separate(coord, c('chr', 'pos'), ':') %>%
            separate(pos, c('start', 'end'), '-')
        
        for (i in 1:length(opt$region$chr)) {
            Chr = opt$region$chr[i]
            Start = opt$region$start[i]
            End = opt$region$end[i]
            
            trac_toplot = signal_smoothed %>%
                filter(
                    chr %in% Chr,
                    (start >= Start & end <= End) |
                        (start <= Start & end >= Start) |
                        (start <= End & end >= End)
                ) %>%
                mutate(
                    start = ifelse(start < Start, Start, start),
                    end = ifelse(end > End , End, end)
                )
            
            
            s50_toplot = s50 %>%
                ungroup() %>%
                filter(
                    chr %in% Chr,
                    (start >= Start & end <= End) |
                        (start <= Start & end >= Start) |
                        (start <= End & end >= End)
                ) %>%
                mutate(
                    RT = RT + 1,
                    start = ifelse(start < Start, Start, start),
                    end = ifelse(end > End , End, end)
                )
            
            plot = trac_toplot %>%
                ggplot() +
                geom_rect(
                    aes(
                        xmin = start,
                        xmax = end,
                        ymin = -newIndex,
                        ymax = -newIndex - 1,
                        fill = Rep + 1
                    )
                )  +
                geom_rect(
                    data = s50_toplot,
                    aes(
                        xmin = start,
                        xmax = end,
                        ymin = 11,
                        ymax = 20,
                        fill = RT
                    ),
                    inherit.aes = F
                ) +
                facet_grid( ~ chr, scale = 'free') +
                scale_fill_gradient(low = 'red',
                                    high = 'green',
                                    limits = c(1, 2)) +
                scale_x_continuous(
                    labels = function(x)
                        paste(x / 1000000, 'Mb', sep = ' ')
                ) +
                scale_y_discrete(
                    limits = c(15, 5, -seq(
                        10, max(trac_toplot$index), round((max(
                            trac_toplot$index
                        ) - 10) / 20)
                    )),
                    labels = c('RT', 'RT reference', seq(
                        10, max(trac_toplot$index), round((max(
                            trac_toplot$index
                        ) - 10) / 20)
                    ))
                ) +
                labs(y = 'S phase progression', fill = "CNV Binary\n") + theme(legend.position = 'top')
            if ('referenceRT' %in% names(opt)) {
                RT_toplot = Reference_RT %>%
                    filter(chr %in% Chr,
                           start >= Start,
                           end <= End) %>%
                    mutate(
                        RT = 1 + RT / max(RT, na.rm = T),
                        start = ifelse(start < Start, Start, start),
                        end = ifelse(end > End , End, end)
                    )
                
                plot = plot +
                    geom_rect(
                        data = RT_toplot,
                        aes(
                            xmin = start,
                            xmax = end,
                            ymin = 0,
                            ymax = 10,
                            fill = RT
                        ),
                        inherit.aes = F
                    )
                
            }
            plot=plot+facet_wrap(~basename)
            ggsave(
                plot,
                filename = paste0(
                    opt$out,
                    '/regions/',
                    paste(opt$base_name,collapse = '_'),
                    '_plot_RT_',
                    Chr,
                    ':',
                    Start,
                    '-',
                    End,
                    '.pdf'
                )
            )
        }
    }
}


##### RT distributions
if ('referenceRT' %in% names(opt)) {
    
RTs = rbind(
    s50 %>%
        ungroup() %>%
        select(chr,start,end,RT,basename)%>%
        mutate(RT = RT + 1),
    
    Reference_RT %>%
        mutate(RT = 1 + RT / max(RT, na.rm = T),
               basename='Reference')
)
}else{
    RTs =  s50 %>%
            ungroup() %>%
            select(chr,start,end,RT,basename)%>%
                             mutate(RT = RT + 1)
}
pdf(
    paste0(
        opt$out,
        '/',paste(opt$base_name,collapse = '_')
        ,
        '_RT_distribution_plot_calculated_RT_vs_reference.pdf'
    )
)

RTs %>%
    ggplot(aes(RT, fill = basename)) + geom_density(alpha = 0.2, aes(y = ..count.. /
                                                                   sum(..count..)))

dev.off()
##### Correlation Calculated RT and reference

RT_type = unique(RTs$basename)

if (length(RT_type)!=1) {
    results = matrix(
        nrow = length(RT_type),
        ncol = length(RT_type),
        dimnames = list(RT_type, RT_type)
    )
    
    
    RTs = RTs %>%
        spread(key = basename, value = RT) %>%
        filter(complete.cases(.))
    
    for (i in RT_type) {
        for (h in RT_type) {
            x = RTs[, i] %>%
                pull()
            y = RTs[, h] %>%
                pull()
            results[i, h] = cor.test(x = x,
                                     y = y,
                                     method = 'pearson')$estimate
        }
    }
    
    
    pdf(paste0(
        opt$out,
        '/',
        paste(opt$base_name, collapse = '_'),
        '_correlation_plot_RTs.pdf'
    ))
    heatmap.2(
        results,
        cellnote = round(results, 2),
        notecex = 0.5,
        notecol = 'black',
        trace = "none",
        dendrogram = 'none',
        Colv = FALSE,
        Rowv = FALSE,
        cexRow = 1,
        cexCol = 1,
        lhei = c(1, 3),
        lwid = c(1, 3),
        breaks = seq(0, 1, length.out = 101),
        srtCol = 0,
        srtRow = 90,
        adjRow = 0.5,
        adjCol = 0.5,
        density.info = 'density',
        key.title = 'Pearson'
    )
    dev.off()
    
}
#matrix for the correlation
mat=signal_smoothed%>%
    mutate(group=str_pad(group,5,pad = '0'),
           index=str_pad(index,5,pad = '0'))%>%
    unite(index,c(basename,group,index),sep = ' _ ')%>%
    unite(pos,c(chr,start),sep = ':')%>%
    mutate(Rep=as.numeric(Rep))%>%
    select(pos,index,Rep)%>%
    spread(key = index,value = Rep)%>%
    column_to_rownames('pos')%>%
    filter(complete.cases(.))%>%
    as.matrix()
#correlation
results=cor(mat,mat)

groups=as.numeric(str_extract(colnames(mat),' [0-9]{5}'))
basenames=str_extract(colnames(mat),'[a-z|A-Z|0-9]{1,100} ')
basename_n=basenames
for (i in 1:length(unique(basename_n))){
    basename_n[basename_n==unique(basename_n)[i]]=i
}

write.table(results, file= paste0(
    opt$out,
    '/',
    paste(opt$base_name,collapse = '_'),
    '_correlation_per_cell.mx'
), row.names=TRUE, col.names=TRUE)

color=colorRampPalette(colors = c('blue','white','red'))
selcol <- colorRampPalette(brewer.pal(12,"Set3"))
selcol2 <- colorRampPalette(brewer.pal(9,"Set1"))
colors_groups=selcol(length(unique(groups)))
color_basebanes=selcol2(length(unique(basename_n)))
pdf(
    paste0(
        opt$out,
        '/',
        paste(opt$base_name,collapse = '_'),
        '_correlation_plot_per_cell.pdf'
    )
)

heatmap.2(
    results,
    trace = "none",
    dendrogram = 'none',
    Colv = F,
    Rowv = F,
    breaks = seq(-1, 1, length.out = 101),
    col=color(100),
    density.info =  'density',
    key.title = 'Pearson',
    RowSideColors=colors_groups[groups],
    ColSideColors=color_basebanes[as.numeric(basename_n)],
    labRow = FALSE, 
    labCol = FALSE
)
dev.off()

#calculate variace within groups and basename 

variance_pearson_per_group=tibble()
unique_groups=unique(groups)
unique_basename=unique(basenames)
for (h in unique_basename) {
    for (i in unique_groups) {
        w = which(groups == i & basenames == h)
        
        x = summary(results[w, w][upper.tri(results[w, w])])
        x["Variance"] = var(results[w, w][upper.tri(results[w, w])])
        nx = names(x)
        x = t(tibble(as.vector(x)))
        colnames(x) = nx
        rownames(x) = paste('Pearson Correlation in group', i, 'for', h)
        variance_pearson_per_group = rbind(variance_pearson_per_group, x)
    }
}
variance_pearson_per_group=variance_pearson_per_group%>%
    rownames_to_column()
write_delim(variance_pearson_per_group,paste0(
    opt$out,
    '/',
    paste(opt$base_name,collapse = '_'),
    '_variance_stat.txt'
), delim = '\t')

#remouve tmp file
system(paste0('rm ', opt$out, '/tmp.tsv'))
system(paste0('rm ', opt$out, '/tmp_bg.tsv'))

# filter out extreme bins
bins=signal_smoothed%>%
    group_by(basename)%>%
    filter(CN_bg < quantile(CN_bg,c(0.01))[[1]] |
               CN_bg > quantile(CN_bg,c(0.99))[[1]])%>%
    mutate(coor=paste0(chr,':',start,'-',end))%>%
    ungroup()%>%
    select(coor)%>%
    unique()

signal_smoothed=signal_smoothed%>%
    filter(!paste0(chr,':',start,'-',end) %in% bins$coor)


s50=s50%>%
    filter(!paste0(chr,':',start,'-',end) %in% bins$coor)%>%
    group_by(basename)%>%
    mutate(RT=(RT-min(RT))/(max(RT)-min(RT)))

#joing s50 with relative signals
signal_smoothed=signal_smoothed%>%
    inner_join(s50,by = c("basename", "chr", "start", "end"))%>%
    group_by(basename)%>%
    mutate(CN_bg=10*(CN_bg-min(CN_bg))/(max(CN_bg)-min(CN_bg)),
           RT=10*RT)


#test multiple time windows
rt=seq(0,10,0.05)

percentage=function(df,seq){
    x=df%>%
        mutate(rep=CN_bg <= seq)%>%
        group_by(chr,start,end,basename)%>%
        summarise(percentage=mean(rep),
                  RT=unique(RT))%>%
        mutate(time=RT-seq)
    return(x)
}

x=lapply(X =rt,FUN = percentage, df=signal_smoothed)

x=do.call('rbind',x)

x=x%>%
    filter(time > -10,
           time < 10)

# calculate distance in time per bin
s50_bin=x%>%
    mutate(dist=abs(0.5-percentage))%>%
    group_by(basename,chr,start,end)%>%
    filter(dist==min(dist))%>%
    mutate(time50=time)%>%
    select(-dist,-time,-percentage,-RT)

x=x%>%
    inner_join(s50_bin, by = c("chr", "start", "end", "basename"))%>%
    group_by(basename)%>%
    mutate(time=time-time50,
           Early=ifelse(RT > median(RT),'Early','Late'))

#calculate tresholds 25% 75% replication keeping in account early and late domains
t = foreach(
    basename = unique(x$basename),
    .combine = 'rbind',
    .packages = c('tidyverse', 'mgcv','foreach')
) %do% {
    temp=foreach(
        EL = unique(x$Early),
        .combine = 'rbind',
        .packages = c('tidyverse', 'mgcv','foreach')
    ) %do% {
        T25_75 = function(df, name,EL) {
            try(library(mgcv))
            
            model = gam(formula = percentage ~ s(time), data = df[, c('percentage', 'time')])
            min = min(df$time)
            max = max(df$time)
            data = predict.gam(model, newdata = data.frame(time = seq(min, max, 0.01)))
            result = data.frame(
                time = seq(min, max, 0.01),
                percentage = data,
                basename = name
            )
            t = inner_join(
                result %>%
                    mutate(distance = abs(percentage - 0.75)) %>%
                    group_by(basename) %>%
                    mutate(min = min(distance)) %>%
                    filter(distance == min) %>%
                    select(basename, time) %>%
                    `colnames<-`(c('basename', 't75')),
                result %>%
                    mutate(distance = abs(percentage - 0.25)) %>%
                    group_by(basename) %>%
                    mutate(min = min(distance)) %>%
                    filter(distance == min) %>%
                    select(basename, time) %>%
                    `colnames<-`(c('basename', 't25')),
                by = "basename")%>%
                mutate(Early=EL)
            
            return(t)
        }
        
        t = T25_75(df = x[x$basename == basename & x$Early==EL, ], basename,EL)
    }
    temp
}
#calculate tresholds 25% 75% replication without keeping in account early and late domains

p=ggplot(x,aes(y=percentage,x=time))+
    geom_point(alpha=0.2)+
    stat_smooth(method = 'gam' , formula =  y ~ s(x))+
    geom_hline(yintercept = c(0.75,0.25),color='yellow' )+
    geom_vline(data = t, aes(xintercept = t75 ),color='red',inherit.aes = T)+
    geom_vline(data = t, aes(xintercept = t25 ),color='red',inherit.aes = T)+
    geom_text(data = t,aes(label=paste('Twidth ~', round((t25-t75),1),'h'),
                           x= (t25 + (t75-t25)/2 ),y=-0, vjust=1.5), color='black',inherit.aes = T)+
    facet_grid(Early~basename)+
    scale_x_reverse()

ggsave(
    p,
    filename = paste0(
        opt$out,
        '/',
        paste(opt$base_name,collapse = '_'),
        '_variability_plot_Early_Late.pdf'
    )
)

t = foreach(
    basename = unique(x$basename),
    .combine = 'rbind',
    .packages = c('tidyverse', 'mgcv')
) %do% {
    T25_75 = function(df, name) {
        try(library(mgcv))
        
        model = gam(formula = percentage ~ s(time), data = df[, c('percentage', 'time')])
        min = min(df$time)
        max = max(df$time)
        data = predict.gam(model, newdata = data.frame(time = seq(min, max, 0.01)))
        result = data.frame(
            time = seq(min, max, 0.01),
            percentage = data,
            basename = name
        )
        t = inner_join(
            result %>%
                mutate(distance = abs(percentage - 0.75)) %>%
                group_by(basename) %>%
                mutate(min = min(distance)) %>%
                filter(distance == min) %>%
                select(basename, time) %>%
                `colnames<-`(c('basename', 't75')),
            result %>%
                mutate(distance = abs(percentage - 0.25)) %>%
                group_by(basename) %>%
                mutate(min = min(distance)) %>%
                filter(distance == min) %>%
                select(basename, time) %>%
                `colnames<-`(c('basename', 't25')),
            by = "basename")
        
        return(t)
    }
    
    t = T25_75(df = x[x$basename == basename, ], basename)
}



p=ggplot(x,aes(y=percentage,x=time))+
    geom_point(alpha=0.2)+
    stat_smooth(method = 'gam' , formula =  y ~ s(x))+
    geom_hline(yintercept = c(0.75,0.25),color='yellow' )+
    geom_vline(data = t, aes(xintercept = t75 ),color='red',inherit.aes = T)+
    geom_vline(data = t, aes(xintercept = t25 ),color='red',inherit.aes = T)+
    geom_text(data = t,aes(label=paste('Twidth ~', round((t25-t75),1),'h'),
                           x= (t25 + (t75-t25)/2 ),y=-0, vjust=1.5), color='black',inherit.aes = T)+
    facet_grid(~basename)+
    scale_x_reverse()

ggsave(
    p,
    filename = paste0(
        opt$out,
        '/',
        paste(opt$base_name,collapse = '_'),
        '_variability_plot.pdf'
    )
)    


if (opt$Var_against_reference) {
    
    Reference_RT=Reference_RT%>%
        filter(!paste0(chr,':',start,'-',end) %in% bins$coor)
   

#joing s50 with relative signals
signal_smoothed=signal_smoothed%>%
    select(-RT)%>%
    inner_join(Reference_RT,by = c("chr", "start", "end"))%>%
    group_by(basename)%>%
    mutate(CN_bg=10*(CN_bg-min(CN_bg))/(max(CN_bg)-min(CN_bg)),
           RT=10*RT)


#test multiple time windows
rt=seq(0,10,0.05)

percentage=function(df,seq){
    x=df%>%
        mutate(rep=CN_bg <= seq)%>%
        group_by(chr,start,end,basename)%>%
        summarise(percentage=mean(rep),
                  RT=unique(RT))%>%
        mutate(time=RT-seq)
    return(x)
}

x=lapply(X =rt,FUN = percentage, df=signal_smoothed)

x=do.call('rbind',x)

x=x%>%
    filter(time > -10,
           time < 10)

# calculate distance in time per bin
RT_ref_bins=x%>%
    mutate(dist=abs(0.5-percentage))%>%
    group_by(basename,chr,start,end)%>%
    filter(dist==min(dist))%>%
    mutate(time50=time)%>%
    select(-dist,-time,-percentage,-RT)

x=x%>%
    inner_join(RT_ref_bins, by = c("chr", "start", "end", "basename"))%>%
    group_by(basename)%>%
    mutate(time=time-time50,
           Early=ifelse(RT > median(RT),'Early','Late'))

#calculate tresholds 25% 75% replication keeping in account early and late domains
t = foreach(
    basename = unique(x$basename),
    .combine = 'rbind',
    .packages = c('tidyverse', 'mgcv','foreach')
) %do% {
    temp=foreach(
        EL = unique(x$Early),
        .combine = 'rbind',
        .packages = c('tidyverse', 'mgcv','foreach')
    ) %do% {
        T25_75 = function(df, name,EL) {
            try(library(mgcv))
            
            model = gam(formula = percentage ~ s(time), data = df[, c('percentage', 'time')])
            min = min(df$time)
            max = max(df$time)
            data = predict.gam(model, newdata = data.frame(time = seq(min, max, 0.01)))
            result = data.frame(
                time = seq(min, max, 0.01),
                percentage = data,
                basename = name
            )
            t = inner_join(
                result %>%
                    mutate(distance = abs(percentage - 0.75)) %>%
                    group_by(basename) %>%
                    mutate(min = min(distance)) %>%
                    filter(distance == min) %>%
                    select(basename, time) %>%
                    `colnames<-`(c('basename', 't75')),
                result %>%
                    mutate(distance = abs(percentage - 0.25)) %>%
                    group_by(basename) %>%
                    mutate(min = min(distance)) %>%
                    filter(distance == min) %>%
                    select(basename, time) %>%
                    `colnames<-`(c('basename', 't25')),
                by = "basename")%>%
                mutate(Early=EL)
            
            return(t)
        }
        
        t = T25_75(df = x[x$basename == basename & x$Early==EL, ], basename,EL)
    }
    temp
}
#calculate tresholds 25% 75% replication without keeping in account early and late domains

p=ggplot(x,aes(y=percentage,x=time))+
    geom_point(alpha=0.2)+
    stat_smooth(method = 'gam' , formula =  y ~ s(x))+
    geom_hline(yintercept = c(0.75,0.25),color='yellow' )+
    geom_vline(data = t, aes(xintercept = t75 ),color='red',inherit.aes = T)+
    geom_vline(data = t, aes(xintercept = t25 ),color='red',inherit.aes = T)+
    geom_text(data = t,aes(label=paste('Twidth ~', round((t25-t75),1),'h'),
                           x= (t25 + (t75-t25)/2 ),y=-0, vjust=1.5), color='black',inherit.aes = T)+
    facet_grid(Early~basename)+
    scale_x_reverse()

ggsave(
    p,
    filename = paste0(
        opt$out,
        '/',
        paste(opt$base_name,collapse = '_'),
        '_variability_plot_Early_Late_ref_RT.pdf'
    )
)

t = foreach(
    basename = unique(x$basename),
    .combine = 'rbind',
    .packages = c('tidyverse', 'mgcv')
) %do% {
    T25_75 = function(df, name) {
        try(library(mgcv))
        
        model = gam(formula = percentage ~ s(time), data = df[, c('percentage', 'time')])
        min = min(df$time)
        max = max(df$time)
        data = predict.gam(model, newdata = data.frame(time = seq(min, max, 0.01)))
        result = data.frame(
            time = seq(min, max, 0.01),
            percentage = data,
            basename = name
        )
        t = inner_join(
            result %>%
                mutate(distance = abs(percentage - 0.75)) %>%
                group_by(basename) %>%
                mutate(min = min(distance)) %>%
                filter(distance == min) %>%
                select(basename, time) %>%
                `colnames<-`(c('basename', 't75')),
            result %>%
                mutate(distance = abs(percentage - 0.25)) %>%
                group_by(basename) %>%
                mutate(min = min(distance)) %>%
                filter(distance == min) %>%
                select(basename, time) %>%
                `colnames<-`(c('basename', 't25')),
            by = "basename")
        
        return(t)
    }
    
    t = T25_75(df = x[x$basename == basename, ], basename)
}



p=ggplot(x,aes(y=percentage,x=time))+
    geom_point(alpha=0.2)+
    stat_smooth(method = 'gam' , formula =  y ~ s(x))+
    geom_hline(yintercept = c(0.75,0.25),color='yellow' )+
    geom_vline(data = t, aes(xintercept = t75 ),color='red',inherit.aes = T)+
    geom_vline(data = t, aes(xintercept = t25 ),color='red',inherit.aes = T)+
    geom_text(data = t,aes(label=paste('Twidth ~', round((t25-t75),1),'h'),
                           x= (t25 + (t75-t25)/2 ),y=-0, vjust=1.5), color='black',inherit.aes = T)+
    facet_grid(~basename)+
    scale_x_reverse()

ggsave(
    p,
    filename = paste0(
        opt$out,
        '/',
        paste(opt$base_name,collapse = '_'),
        '_variability_plot_ref_RT.pdf'
    )
)
}

    
  
 