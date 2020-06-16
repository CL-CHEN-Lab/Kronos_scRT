#!/usr/local/bin/Rscript --slave
#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)
options(warn = 1, scipen = 999)

option_list = list(
    make_option(
        c("-L", "--List"),
        type = "character",
        default = NULL,
        help = "A Tab separated file containing in each colum scRT_Tracks,scTW_Tracks and scCNV files paths. Alternative to -R,-T and -C options.",
        metavar = "character"
    ), make_option(
        c("-R", "--scRT_Tracks"),
        type = "character",
        default = NULL,
        help = "*calculated_replication_timing* file(s) created by Kronos RT. If multiple files are provided they have to be separated by a comma.  Alternative to -L option.",
        metavar = "character"
    ),make_option(
        c("-T", "--scTW_Tracks"),
        type = "character",
        default = NULL,
        help = "*calculated_Twhith* file(s) created by Kronos RT. If multiple files are provided they have to be separated by a comma. Alternative to -L option.",
        metavar = "character"
    ),
    make_option(
        c("-C", "--scCNV"),
        type = "character",
        default = NULL,
        help = "*single_cells_CNV* file(s) created by Kronos RT. If multiple files are provided they have to be separated by a comma.  Alternative to -L option..",
        metavar = "character"
    ),
    make_option(
        c("-s", "--order"),
        type = "character",
        default = NULL,
        help = "basenames separated by a comma in the desired order for plotting.",
        metavar = "character"
    ),
    make_option(
        c("--CNV_values"),
        type = "character",
        default = "B",
        help = "What type of date to plot for the sigle cell traks: ('B'=Binarized, 'CNV'=Copy number variation or 'log2'=log2(CNV_Cell/CNV_mean_G1/G2_cells) 'all'= one file per option) [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("-r", "--region"),
        type = "character",
        default = NULL,
        help = "Region to plot  chr:start-end (multiple regins can be separated by a comma) or provided as a bed file",
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
suppressPackageStartupMessages(library(RColorBrewer, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))
suppressPackageStartupMessages(library(gridExtra, quietly = TRUE))

if('List' %in% names(opt)) {
    opt$List = tryCatch(
        expr = read_tsv(opt$List, col_names = F, col_types = cols()),
        error = function(e) {
            stop('Settings file does not exitst. See script usage (--help)')
        }
    )
    
    opt$scRT_Tracks = opt$List$X1
    opt$scTW_Tracks = opt$List$X2
    opt$scCNV = opt$List$X3
    
} else{
    if (!'scRT_Tracks' %in% names(opt)) {
        stop("scRT_Tracks file must be provided. See script usage (--help)")
    }
    if (!'scTW_Tracks' %in% names(opt)) {
        stop("scTW_Tracks file must be provided. See script usage (--help)")
    }
    if (!'scCNV' %in% names(opt)) {
        stop("scCNV file must be provided. See script usage (--help)")
    }
    
    opt$scRT_Tracks = str_split(opt$scRT_Tracks, ',')[[1]]
    opt$scTW_Tracks = str_split(opt$scTW_Tracks, ',')[[1]]
    opt$scCNV = str_split(opt$scCNV, ',')[[1]]
    
}

if (!'region' %in% names(opt)) {
    stop("No reagion of interest has been selected. See script usage (--help)")
    
}

if('order'%in% names(opt)){
    opt$order = str_split(opt$order, ',')[[1]]
    
}


scRT_TW=foreach(i=1:length(opt$scRT_Tracks),.packages = 'tidyverse',.combine = 'rbind')%do%{
    tmp=inner_join(read_tsv(opt$scRT_Tracks[i],col_types = cols()),
    read_tsv(opt$scTW_Tracks[i],col_types = cols()), by = c("chr", "start", "end", "basename"))
    if('order'%in% names(opt)){
        tmp%>%
            mutate(
                basename=factor(basename, levels=opt$order)
                ) 
        
    }else{
        tmp
    }
}

scCNV=foreach(i=1:length(opt$scCNV),.packages = 'tidyverse',.combine = 'rbind')%do%{
    tmp=read_tsv(opt$scCNV[i],col_types = cols())
    if('order'%in% names(opt)){
        tmp%>%
            mutate(
                basename=factor(basename, levels=opt$order)
            ) 
        
    }else{
        tmp
    }
}

#create directory
if (str_extract(opt$out, '.$') != '/') {
    opt$out = paste0(opt$out, '/')
}
system(paste0('mkdir -p ', opt$out, '/regions'))
    
if (file.exists(opt$region)) {
        #load bed file if exist
        opt$region = read_tsv(opt$region, col_names = c('chr', 'start', 'end'),col_types = cols()) %>%
            mutate(
                n_0_start = str_length(str_extract(start, '0{1,10}$')),
                n_0_end = str_length(str_extract(end, '0{1,10}$')),
                unit = min(
                    factor(
                        case_when(
                            is.na(n_0_start) ~ 'bp',
                            n_0_start < 3 ~ 'bp',
                            n_0_start < 6 ~ 'Kb',
                            n_0_start >= 6 ~ 'Mp'
                        ),
                        levels = c('bp', 'Kb', 'Mp'),
                        ordered = TRUE
                    ),
                    factor(
                        case_when(
                            is.na(n_0_end) ~ 'bp',
                            n_0_end < 3 ~ 'bp',
                            n_0_end < 6 ~ 'Kb',
                            n_0_end >= 6 ~ 'Mp'
                        ),
                        levels = c('bp', 'Kb', 'Mp'),
                        ordered = TRUE
                    )
                ),
                name_reg = paste(
                    chr,
                    case_when(
                        unit == 'bp' ~ paste0(start, unit, '_', end, unit),
                        unit == 'Kb' ~ paste0(start / 10 ^ 3, unit, '_', end /
                                                  10 ^ 3, unit),
                        unit == 'Mb' ~ paste0(start / 10 ^ 6, unit, '_', end /
                                                  10 ^ 6, unit)
                    ),
                    sep = '_'
                )
            ) %>%
            dplyr::select(-unit, -n_0_start, -n_0_end)
    } else{
        #reshape regions
        opt$region = tibble(coord = str_split(opt$region, pattern = ',')[[1]]) %>%
            mutate(name_reg = str_replace_all(coord,
                                              pattern = '[-:]',
                                              replacement = '_')) %>%
            separate(coord, c('chr', 'pos'), ':') %>%
            separate(pos, c('start', 'end'), '-') %>%
            mutate(
                start_unit = str_extract(start, pattern = '.{2}$'),
                start = as.numeric(str_remove(
                    start, "[Bb][Pp]|[Kk][Bb]|[Mm][Bb]"
                )) * case_when(
                    grepl(x = start_unit, pattern =  '[Kk][Bb]') ~ 1000,
                    grepl(x = start_unit, pattern = '[Mm][Bb]') ~ 1000000,
                    grepl(x = start_unit, pattern = '[Bp][Pp]') ~ 1,
                    grepl(x = start_unit, pattern =  '[0-9][0-9]') ~ 1
                ),
                end_unit = str_extract(end, pattern = '.{2}$'),
                
                end = as.numeric(str_remove(
                    end, "[Bb][Pp]|[Kk][Bb]|[Mm][Bb]"
                )) * case_when(
                    grepl(x = end_unit, pattern =  '[Kk][Bb]') ~ 1000,
                    grepl(x = end_unit, pattern = '[Mm][Bb]') ~ 1000000,
                    grepl(x = end_unit, pattern = '[Bp][Pp]') ~ 1,
                    grepl(x = end_unit, pattern =  '[0-9][0-9]') ~ 1
                )
            ) %>%
            dplyr::select(-start_unit, -end_unit)
    }
    
if(!opt$CNV_values %in% c('B','CNV','Log2','all')){
    opt$CNV_values='B'
    warning('unrecognized CNV value to plot. Binarized tracks restored')
}
    max_TW=max(scRT_TW$TW)
for (i in 1:length(opt$region$chr)) {
    Chr = opt$region$chr[i]
    Start = opt$region$start[i]
    End = opt$region$end[i]
    name_reg = opt$region$name_reg[i]
    
    scRT_TW_toplot = scRT_TW %>%
        filter(
            chr %in% Chr,
            (start >= Start & end <= End) |
                (start <= Start & end >= Start) |
                (start <= End & end >= End)
        ) %>%
        mutate(start = ifelse(start < Start, Start, start),
               end = ifelse(end > End , End, end))
    
    scCNV_toplot = scCNV %>%
        filter(
            chr %in% Chr,
            (start >= Start & end <= End) |
                (start <= Start & end >= Start) |
                (start <= End & end >= End)
        ) %>%
        mutate(start = ifelse(start < Start, Start, start),
               end = ifelse(end > End , End, end))
    
    if (length(scRT_TW_toplot$chr) != 0) {
        Maxi = max(scCNV_toplot$newIndex)
        n=str_count(Maxi)
        if (opt$CNV_values == 'B') {
            p1 = ggplot() + geom_line(
                data = scRT_TW_toplot %>%
                    mutate(
                        mid = (start + end) / 2,
                        RT = RT * Maxi / 6 + Maxi /
                            20
                    ) %>%
                    gather(what, pos, start, end) %>%
                    arrange(mid, pos),
                aes(pos, RT, color = TW)
            ) +
                geom_rect(
                    data = scCNV_toplot %>%
                        mutate(Rep = ifelse(
                            Rep, "Replicated", "Unreplicated"
                        )),
                    
                    aes(
                        xmin = start,
                        xmax = end,
                        ymin = -newIndex,
                        ymax = -newIndex - 1,
                        fill = Rep
                    )
                ) +
                annotate(
                    "rect",
                    xmin = -Inf,
                    xmax = Inf,
                    ymin = Maxi / 100,
                    ymax = Maxi / 20 - Maxi / 100,
                    fill = 'white'
                ) +
                facet_grid( ~ basename) +
                scale_fill_manual(values = c(
                    "Replicated" = 'green',
                    "Unreplicated" = 'red'
                )) +
                scale_color_gradient(low = '#236AB9',
                                     high = '#FC7307',
                                     limits = c(0,max_TW)) +
                scale_y_continuous(
                    breaks = c(Maxi / 6 + Maxi / 20, Maxi / 3 + Maxi / 20, Maxi / 20),
                    labels = c('Early - 1', 'Mind - 0.5', 'Late - 0'),
                    name = 'RT',
                    sec.axis = sec_axis(
                        ~ .,
                        breaks = c(1.5, -seq(1, round(
                            Maxi / 10 ^ (n - 2)
                        ), 1) * 10 ^ (n - 2))-0.5,
                        labels = as.character(c(1,
                                                seq(
                                                    1, round(Maxi / 10 ^ (n - 2)), 1
                                                ) * 10 ^ (n - 2))),
                        name = 'Single Cells tracks ordered by S-phase progression'
                    )
                ) +
                scale_x_continuous(
                    labels = function(x)
                        paste(x / 1000000, 'Mb', sep = ' ')
                ) + theme(
                    legend.position = 'bottom',
                    axis.text.x = element_text(
                        angle = 45,
                        hjust = 1,
                        vjust = 1
                    ),
                    axis.title.y.right = element_text(hjust = 0.8),
                    axis.title.y.left  = element_text(hjust = 0.92)
                    
                ) + xlab(Chr)+labs(fill='')
            ggsave(
                plot = p1,
                filename = paste0(
                    opt$out,
                    '/regions/',
                    opt$output_file_base_name,
                    '_scPlot_binarized',
                    name_reg,
                    '.pdf'
                ))
        } else if (opt$CNV_values == 'CNV') {
            x = median(scCNV_toplot$background)
            p1 = ggplot() + geom_line(
                data = scRT_TW_toplot %>%
                    mutate(
                        mid = (start + end) / 2,
                        RT = RT * Maxi / 6 +
                            Maxi / 20
                    ) %>%
                    gather(what, pos, start, end) %>%
                    arrange(mid, pos),
                aes(pos, RT, color = TW)
            ) +
                geom_rect(
                    data = scCNV_toplot,
                    
                    aes(
                        xmin = start,
                        xmax = end,
                        ymin = -newIndex,
                        ymax = -newIndex - 1,
                        fill = CN
                    )
                ) +
                annotate(
                    "rect",
                    xmin = -Inf,
                    xmax = Inf,
                    ymin = Maxi / 100,
                    ymax = Maxi / 20 - Maxi / 100,
                    fill = 'white'
                ) +
                facet_grid( ~ basename) +
                scale_color_gradient(low = '#236AB9',
                                     high = '#FC7307',
                                     limits = c(0,max_TW)) +
                scale_y_continuous(
                    breaks = c(Maxi / 6 + Maxi / 20, Maxi / 3 + Maxi / 20, Maxi / 20),
                    labels = c('Early - 1', 'Mind - 0.5', 'Late - 0'),
                    name = 'RT',
                    sec.axis = sec_axis(
                        ~ .,
                        breaks = c(1.5, -seq(1, round(
                            Maxi / 10 ^ (n - 2)
                        ), 1) * 10 ^ (n - 2))-0.5,
                        labels = as.character(c(1,
                                                seq(
                                                    1, round(Maxi / 10 ^ (n - 2)), 1
                                                ) * 10 ^ (n - 2))),
                        name = 'Single Cells tracks ordered by S-phase progression'
                    )
                ) +
                scale_x_continuous(
                    labels = function(x)
                        paste(x / 1000000, 'Mb', sep = ' ')
                ) + theme(
                    legend.position = 'bottom',
                    axis.text.x = element_text(
                        angle = 45,
                        hjust = 1,
                        vjust = 1
                    ),
                    axis.title.y.right = element_text(hjust = 0.8),
                    axis.title.y.left  = element_text(hjust = 0.92)
                    
                ) + xlab(Chr) +
                
                scale_fill_gradient(low = 'blue',
                                    high = 'orange',
                                    limits = c(0, 3 * x)) +
                labs(fill= 'CNV')
            ggsave(
                plot = p1,
                filename = paste0(
                    opt$out,
                    '/regions/',
                    opt$output_file_base_name,
                    '_scPlot_CNV',
                    name_reg,
                    '.pdf'
                ))
            
        } else if (opt$CNV_values == 'log2') {
            p1 = ggplot() + geom_line(
                data = scRT_TW_toplot %>%
                    mutate(
                        mid = (start + end) / 2,
                        RT = RT * Maxi / 6 + Maxi /
                            20
                    ) %>%
                    gather(what, pos, start, end) %>%
                    arrange(mid, pos),
                aes(pos, RT, color = TW)
            ) +
                geom_rect(
                    data = scCNV_toplot ,
                    
                    aes(
                        xmin = start,
                        xmax = end,
                        ymin = -newIndex,
                        ymax = -newIndex - 1,
                        fill = CN_bg
                    )
                ) +
                annotate(
                    "rect",
                    xmin = -Inf,
                    xmax = Inf,
                    ymin = Maxi / 100,
                    ymax = Maxi / 20 - Maxi / 100,
                    fill = 'white'
                ) +
                facet_grid( ~ basename) +
                scale_color_gradient(low = '#236AB9',
                                     high = '#FC7307',
                                     limits = c(0,max_TW)) +
                scale_y_continuous(
                    breaks = c(Maxi / 6 + Maxi / 20, Maxi / 3 + Maxi / 20, Maxi / 20),
                    labels = c('Early - 1', 'Mind - 0.5', 'Late - 0'),
                    name = 'RT',
                    sec.axis = sec_axis(
                        ~ .,
                        breaks = c(1.5, -seq(1, round(
                            Maxi / 10 ^ (n - 2)
                        ), 1) * 10 ^ (n - 2))-0.5,
                        labels = as.character(c(1,
                                                seq(
                                                    1, round(Maxi / 10 ^ (n - 2)), 1
                                                ) * 10 ^ (n - 2))),
                        name = 'Single Cells tracks ordered by S-phase progression'
                    )
                ) +
                scale_x_continuous(
                    labels = function(x)
                        paste(x / 1000000, 'Mb', sep = ' ')
                ) + theme(
                    legend.position = 'bottom',
                    axis.text.x = element_text(
                        angle = 45,
                        hjust = 1,
                        vjust = 1
                    ),
                    axis.title.y.right = element_text(hjust = 0.8),
                    axis.title.y.left  = element_text(hjust = 0.92)
                    
                ) + xlab(Chr) +
                scale_fill_gradient(low = 'purple',
                                    high = 'yellow',
                                    limits = c(-1.5, 1.5),breaks=c(-1.5,0, 1.5)
                                    )+labs(fill = expression(
                                        over(
                                            S[CNV[i]], 
                                            bar(
                                                G1 / 2[CNV]
                                            )
                                        )
                                    ))
            ggsave(
                plot = p1,
                filename = paste0(
                    opt$out,
                    '/regions/',
                    opt$output_file_base_name,
                    '_scPlot_log2',
                    name_reg,
                    '.pdf'
                ))
            
        }else{
            p1 = ggplot() + geom_line(
                data = scRT_TW_toplot %>%
                    mutate(
                        mid = (start + end) / 2,
                        RT = RT * Maxi / 6 + Maxi /
                            20
                    ) %>%
                    gather(what, pos, start, end) %>%
                    arrange(mid, pos),
                aes(pos, RT, color = TW)
            ) +
                geom_rect(
                    data = scCNV_toplot %>%
                        mutate(Rep = ifelse(
                            Rep, "Replicated", "Unreplicated"
                        )),
                    
                    aes(
                        xmin = start,
                        xmax = end,
                        ymin = -newIndex,
                        ymax = -newIndex - 1,
                        fill = Rep
                    )
                ) +
                annotate(
                    "rect",
                    xmin = -Inf,
                    xmax = Inf,
                    ymin = Maxi / 100,
                    ymax = Maxi / 20 - Maxi / 100,
                    fill = 'white'
                ) +
                facet_grid( ~ basename) +
                scale_fill_manual(values = c(
                    "Replicated" = 'green',
                    "Unreplicated" = 'red'
                )) +
                scale_color_gradient(low = '#236AB9',
                                     high = '#FC7307',
                                     limits = c(0,max_TW)) +
                scale_y_continuous(
                    breaks = c(Maxi / 6 + Maxi / 20, Maxi / 3 + Maxi / 20, Maxi / 20),
                    labels = c('Early - 1', 'Mind - 0.5', 'Late - 0'),
                    name = 'RT',
                    sec.axis = sec_axis(
                        ~ .,
                        breaks = c(1.5, -seq(1, round(
                            Maxi / 10 ^ (n - 2)
                        ), 1) * 10 ^ (n - 2))-0.5,
                        labels = as.character(c(1,
                                                seq(
                                                    1, round(Maxi / 10 ^ (n - 2)), 1
                                                ) * 10 ^ (n - 2))),
                        name = 'Single Cells tracks ordered by S-phase progression'
                    )
                ) +
                scale_x_continuous(
                    labels = function(x)
                        paste(x / 1000000, 'Mb', sep = ' ')
                ) + theme(
                    legend.position = 'bottom',
                    axis.text.x = element_text(
                        angle = 45,
                        hjust = 1,
                        vjust = 1
                    ),
                    axis.title.y.right = element_text(hjust = 0.8),
                    axis.title.y.left  = element_text(hjust = 0.92)
                    
                ) + xlab(Chr)+labs(fill='')
            
            x = median(scCNV_toplot$background)
            p2 = ggplot() + geom_line(
                data = scRT_TW_toplot %>%
                    mutate(
                        mid = (start + end) / 2,
                        RT = RT * Maxi / 6 +
                            Maxi / 20
                    ) %>%
                    gather(what, pos, start, end) %>%
                    arrange(mid, pos),
                aes(pos, RT, color = TW)
            ) +
                geom_rect(
                    data = scCNV_toplot,
                    
                    aes(
                        xmin = start,
                        xmax = end,
                        ymin = -newIndex,
                        ymax = -newIndex - 1,
                        fill = CN
                    )
                ) +
                annotate(
                    "rect",
                    xmin = -Inf,
                    xmax = Inf,
                    ymin = Maxi / 100,
                    ymax = Maxi / 20 - Maxi / 100,
                    fill = 'white'
                ) +
                facet_grid( ~ basename) +
                scale_color_gradient(low = '#236AB9',
                                     high = '#FC7307',
                                     limits = c(0,max_TW)) +
                scale_y_continuous(
                    breaks = c(Maxi / 6 + Maxi / 20, Maxi / 3 + Maxi / 20, Maxi / 20),
                    labels = c('Early - 1', 'Mind - 0.5', 'Late - 0'),
                    name = 'RT',
                    sec.axis = sec_axis(
                        ~ .,
                        breaks = c(1.5, -seq(1, round(
                            Maxi / 10 ^ (n - 2)
                        ), 1) * 10 ^ (n - 2))-0.5,
                        labels = as.character(c(1,
                                                seq(
                                                    1, round(Maxi / 10 ^ (n - 2)), 1
                                                ) * 10 ^ (n - 2))),
                        name = 'Single Cells tracks ordered by S-phase progression'
                    )
                ) +
                scale_x_continuous(
                    labels = function(x)
                        paste(x / 1000000, 'Mb', sep = ' ')
                ) + theme(
                    legend.position = 'bottom',
                    axis.text.x = element_text(
                        angle = 45,
                        hjust = 1,
                        vjust = 1
                    ),
                    axis.title.y.right = element_text(hjust = 0.8),
                    axis.title.y.left  = element_text(hjust = 0.92)
                    
                ) + xlab(Chr) +
                
                scale_fill_gradient(low = 'blue',
                                    high = 'orange',
                                    limits = c(0, 3 * x)) +
                labs(fill= 'CNV')
            
            p3 = ggplot() + geom_line(
                data = scRT_TW_toplot %>%
                    mutate(
                        mid = (start + end) / 2,
                        RT = RT * Maxi / 6 + Maxi /
                            20
                    ) %>%
                    gather(what, pos, start, end) %>%
                    arrange(mid, pos),
                aes(pos, RT, color = TW)
            ) +
                geom_rect(
                    data = scCNV_toplot ,
                    
                    aes(
                        xmin = start,
                        xmax = end,
                        ymin = -newIndex,
                        ymax = -newIndex - 1,
                        fill = CN_bg
                    )
                ) +
                annotate(
                    "rect",
                    xmin = -Inf,
                    xmax = Inf,
                    ymin = Maxi / 100,
                    ymax = Maxi / 20 - Maxi / 100,
                    fill = 'white'
                ) +
                facet_grid( ~ basename) +
                scale_color_gradient(low = '#236AB9',
                                     high = '#FC7307',
                                     limits = c(0,max_TW)) +
                scale_y_continuous(
                    breaks = c(Maxi / 6 + Maxi / 20, Maxi / 3 + Maxi / 20, Maxi / 20),
                    labels = c('Early - 1', 'Mind - 0.5', 'Late - 0'),
                    name = 'RT',
                    sec.axis = sec_axis(
                        ~ .,
                        breaks = c(1.5, -seq(1, round(
                            Maxi / 10 ^ (n - 2)
                        ), 1) * 10 ^ (n - 2))-0.5,
                        labels = as.character(c(1,
                                                seq(
                                                    1, round(Maxi / 10 ^ (n - 2)), 1
                                                ) * 10 ^ (n - 2))),
                        name = 'Single Cells tracks ordered by S-phase progression'
                    )
                ) +
                scale_x_continuous(
                    labels = function(x)
                        paste(x / 1000000, 'Mb', sep = ' ')
                ) + theme(
                    legend.position = 'bottom',
                    axis.text.x = element_text(
                        angle = 45,
                        hjust = 1,
                        vjust = 1
                    ),
                    axis.title.y.right = element_text(hjust = 0.8),
                    axis.title.y.left  = element_text(hjust = 0.92)
                    
                ) + xlab(Chr) +
                scale_fill_gradient(low = 'purple',
                                    high = 'yellow',
                                    limits = c(-1.5, 1.5),breaks=c(-1.5,0, 1.5)
                                    )+labs(fill = expression(
                                        over(
                                            S[CNV[i]], 
                                            bar(
                                                G1 / 2[CNV]
                                            )
                                        )
                                    ))
            
            ggsave(
                plot = p1,
                filename = paste0(
                    opt$out,
                    '/regions/',
                    opt$output_file_base_name,
                    '_scPlot_binarized',
                    name_reg,
                    '.pdf'
                )
            )
            ggsave(
                plot = p2,
                filename = paste0(
                    opt$out,
                    '/regions/',
                    opt$output_file_base_name,
                    '_scPlot_CNV',
                    name_reg,
                    '.pdf'
                )
            )
            ggsave(
                plot = p3,
                filename = paste0(
                    opt$out,
                    '/regions/',
                    opt$output_file_base_name,
                    '_scPlot_log2',
                    name_reg,
                    '.pdf'
                )
            )
        }
        
        
    
    }
}
    
    print('done')
    
                  