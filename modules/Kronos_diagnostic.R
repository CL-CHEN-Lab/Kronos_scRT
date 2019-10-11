#!/usr/local/bin/Rscript --slave

# this script is meant to select the treshold to select cycling cells

if (!suppressPackageStartupMessages(require(optparse, quietly = TRUE))) {
    install.packages("optparse", quiet = T)
    suppressPackageStartupMessages(library(optparse, quietly = TRUE))
}

options(stringsAsFactors = FALSE)
options(warn=1) 
option_list = list(
    make_option(
        c("-f", "--file"),
        type = "character",
        default = NULL,
        help = "dataset file name",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "./output",
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
        type = "double",
        help = "Threshold to identify S-phase cells",
        metavar = "double"
    ),
    make_option(
        c("-G", "--threshold_G1G2phase"),
        type = "double",
        help = "Threshold to identify G1-phase cells. -S has to be selected and has to be bigger than -G",
        metavar = "double"
    )
)

opt = parse_args( OptionParser(option_list=option_list))

if (!suppressPackageStartupMessages(require(tidyverse, quietly = TRUE))) {
    install.packages("tidyverse", quiet = T)
    suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
}

#create output directory
if (str_extract(opt$out,'.$')!='/'){
    opt$out=paste0(opt$out,'/')
}

system(paste0('mkdir -p ./', opt$out))

#load data
data<-read_csv(opt$file,
               col_types = cols())

if (!'threshold_Sphase' %in% names(opt)){

    p = data %>%
        mutate(Type = ifelse(
            as.logical(is_high_dimapd) == T &
                as.logical(is_noisy) == T,
            'S-phase',
            ifelse(
                as.logical(is_high_dimapd) == F &
                    as.logical(is_noisy) == T,
                'unknown cells',
                'G1/G2 cells'
            )
        )) %>%
        ggplot(aes(mean_ploidy, normalized_dimapd, color = Type)) +
        geom_point(alpha = 0.3) +
        scale_color_manual(
            values = c(
                'G1/G2 cells' = 'darkred',
                'S-phase' = 'darkgreen',
                'unknown cells' = 'darkorange'
            )
        ) +
        theme(legend.position = 'top', legend.title = element_blank())
    
    system(paste('mkdir -p', opt$out))
    suppressMessages( ggsave(p, filename = paste0(opt$out, '/', opt$base_name, '_plot.pdf')))

}else if('threshold_Sphase' %in% names(opt) & 'threshold_G1G2phase' %in% names(opt) ){
    data = data %>%
        mutate(
            is_high_dimapd = ifelse(normalized_dimapd > opt$threshold_Sphase, T, F),
            is_noisy = ifelse(is_high_dimapd, T, is_noisy),
            Type = ifelse(
                as.logical(is_high_dimapd) == T & as.logical(is_noisy) == T,
                'S-phase',
                ifelse(
                    as.logical(is_high_dimapd) == F &
                        as.logical(is_noisy) == F &
                        normalized_dimapd < opt$threshold_G1G2phase,
                    'G1/G2 cells',
                    'unknown cells'
                    
                )
            )
        )
    
    median_ploidy_not_noisy = median(data$mean_ploidy[data$is_noisy == F &
                                                          data$normalized_dimapd < opt$threshold_G1G2phase   ])
    data = data %>%
        filter(
            mean_ploidy > median_ploidy_not_noisy / 1.3 ,
            mean_ploidy < median_ploidy_not_noisy * 1.8,
            !ploidy_confidence <= 2
        )
    median_ploidy_not_noisy = median(data$mean_ploidy[data$is_noisy == F&
                                                          data$normalized_dimapd < opt$threshold_G1G2phase])
    p = data %>%
        ggplot(aes(mean_ploidy, normalized_dimapd, color = Type)) +
        geom_point(alpha = 0.3) +
        scale_color_manual(
            values = c(
                'G1/G2 cells' = 'darkred',
                'S-phase' = 'darkgreen',
                'unknown cells' = 'darkorange'
            )
        ) +
        theme(legend.position = 'top', legend.title = element_blank()) +
        geom_vline(xintercept = median_ploidy_not_noisy)

    suppressMessages(ggsave(p,
           filename = paste0(opt$out, '/', opt$base_name, '_plot_th_', opt$threshold_Sphase,'-',opt$threshold_G1G2phase, '.pdf')))
}else{

    data = data %>%
        mutate(
            is_high_dimapd = ifelse(normalized_dimapd > opt$threshold_Sphase, T, F),
            is_noisy = ifelse(is_high_dimapd, T, is_noisy),
            Type = ifelse(
                as.logical(is_high_dimapd) == T & as.logical(is_noisy) == T,
                'S-phase',
                ifelse(
                    as.logical(is_high_dimapd) == F &
                        as.logical(is_noisy) == F,
                    'G1/G2 cells',
                    'unknown cells'
                )
            )
        )
    
    median_ploidy_not_noisy = median(data$mean_ploidy[data$is_noisy == F])
    data = data %>%
        filter(
            mean_ploidy > median_ploidy_not_noisy / 1.5 ,
            mean_ploidy < median_ploidy_not_noisy * 1.5,
            !ploidy_confidence <= 2
        )
    median_ploidy_not_noisy = median(data$mean_ploidy[data$is_noisy == F])
    p = data %>%
        ggplot(aes(mean_ploidy, normalized_dimapd, color = Type)) +
        geom_point(alpha = 0.3) +
        scale_color_manual(
            values = c(
                'G1/G2 cells' = 'darkred',
                'S-phase' = 'darkgreen',
                'unknown cells' = 'darkorange'
            )
        ) +
        theme(legend.position = 'top', legend.title = element_blank()) +
        geom_vline(xintercept = median_ploidy_not_noisy)

    suppressMessages( ggsave(p,
           filename = paste0(opt$out, '/', opt$base_name, '_plot_th_', opt$threshold_Sphase, '.pdf')))
}
