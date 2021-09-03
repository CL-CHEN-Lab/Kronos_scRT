#!/usr/local/bin/Rscript

# this script is meant to select the treshold to select cycling cells

suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)
options(warn=1) 
option_list = list(
    make_option(
        c("-F", "--file"),
        type = "character",
        default = NULL,
        help = "Dataset file name",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "./output",
        help = "Output directory [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("-b", "--base_name"),
        type = "character",
        default = "exp",
        help = "Base name for files names [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("-C", "--correct"),
        type = "logical",
        action = "store_true",
        default = F,
        help = "If True diagnostic corrects the S-phase progression and returns a setting file [default= %default]",
        metavar = "logical"
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
    ),
    make_option(
        c("-f", "--Sphase_first_part"),
        type = "double",
        help = "Correction parameter for the first part of the S-phase [0.95,1]",
        metavar = "double"
    ),
    make_option(
        c("-s", "--Sphase_second_part"),
        type = "double",
        help = "Correction parameter for the second part of the S-phase [0.5,0.55]",
        metavar = "double"
    ),
    make_option(
        c("-c", "--cores"),
        type = "integer",
        default = 3,
        help = "Numbers of parallel jobs to run [default= %default] ",
        metavar = "integer"
    ),
    make_option(
        c("-m", "--min_n_reads"),
        type = "double",
        default = 160,
        action = 'store',
        help = "Min n of reads per million per haploid genome to keep a cell in the analysis [default= %default]",
        metavar = "double"
    )
)

opt = parse_args( OptionParser(option_list=option_list))

#load packages
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(LaplacesDemon, quietly = TRUE))
suppressPackageStartupMessages(library(doSNOW, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))

#set plotting theme
theme_set(theme_bw())

#check inputs
if (!'file' %in% names(opt)) {
    stop("Per cell stat file must be provided. See script usage (--help)")
}

#create output directory
if (str_extract(opt$out,'.$')!='/'){
    opt$out=paste0(opt$out,'/')
}

system(paste0('mkdir -p ', opt$out))

#load data
data<-read_csv(opt$file,
               col_types = cols())
if(opt$correct==F){
    if (!'threshold_Sphase' %in% names(opt)){
        
        data=data %>%
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
            )) 
        }else{
            if(!'threshold_G1G2phase' %in% names(opt) ){
                opt$threshold_G1G2phase = opt$threshold_Sphase
            }
            
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
        }
    
    median_ploidy_not_noisy = median(data%>%filter(is_noisy == F)%>%pull(mean_ploidy))
    
    data=data%>%
        mutate(Type= case_when(
            coverage_per_1Mbp < opt$min_n_reads*median_ploidy_not_noisy ~ 'Low Coverage',
            ploidy_confidence < 2 & ploidy_confidence!=-100 ~ 'Low Ployidy confidence',
            mean_ploidy < median_ploidy_not_noisy / 1.5 ~ 'Too low ploidy compared to G1/G2',
            mean_ploidy > median_ploidy_not_noisy * 2 ~ 'Too high ploidy compared to G1/G2',
            T ~ Type
        ))
    
    p =data%>%
        ggplot(aes(mean_ploidy, normalized_dimapd, color = Type)) +
        geom_point(alpha = 0.3) +
        scale_color_manual(
            values = c(
                'Low Coverage' = "#ff7949",
                'Low Ployidy confidence' = "#70001e",
                'Too low ploidy compared to G1/G2' = "#01e7ab",
                'Too high ploidy compared to G1/G2' ="#a7001b",
                'G1/G2 cells' = "#005095",
                'S-phase' = "#78bd3e",
                'unknown cells' = "#dfbd31"  
            )
        ) +
        theme(legend.position = 'top', legend.title = element_blank())+
        xlab('Ploidy') + ylab('Variability')
    
    system(paste('mkdir -p ', opt$out))
    suppressMessages( ggsave(p, filename = paste0(opt$out, '/', opt$base_name, '_no_correction_plot.pdf')))

    print('done')
    quit()
}

if (!'threshold_Sphase' %in% names(opt)){

    data=data %>%
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
        )) 
    
    median_ploidy_not_noisy = median(data%>%
                                         filter(is_noisy == F)%>%pull(mean_ploidy))

    data=data%>%
        filter(coverage_per_1Mbp >= opt$min_n_reads*median_ploidy_not_noisy,
               ploidy_confidence > 2 | ploidy_confidence==-100,
                mean_ploidy > median_ploidy_not_noisy / 1.5 ,
                mean_ploidy < median_ploidy_not_noisy * 2
                   )
     p =data%>%
        ggplot(aes(mean_ploidy, normalized_dimapd, color = Type)) +
        geom_point(alpha = 0.3) +
        scale_color_manual(
            values = c(
                'G1/G2 cells' = "#005095",
                'S-phase' = "#78bd3e",
                'unknown cells' = "#dfbd31"  
            )
        ) +
        theme(legend.position = 'top', legend.title = element_blank())+
         xlab('Ploidy') + ylab('Variability')
    
    system(paste('mkdir -p ', opt$out))
    suppressMessages( ggsave(p, filename = paste0(opt$out, '/', opt$base_name, '_plot.pdf')))
    
    opt$threshold_G1G2phase=NA
    opt$threshold_Sphase=NA
    
}else{
    
    if(!'threshold_G1G2phase' %in% names(opt) ){
        opt$threshold_G1G2phase = opt$threshold_Sphase
    }
    
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
    
    median_ploidy_not_noisy = median(data%>%filter(is_noisy == F)%>%pull(mean_ploidy))
    data = data %>%
        filter(coverage_per_1Mbp >= opt$min_n_reads*median_ploidy_not_noisy,
            mean_ploidy > median_ploidy_not_noisy / 1.5 ,
            mean_ploidy < median_ploidy_not_noisy * 2,
            !ploidy_confidence <= 2 | ploidy_confidence==-100
        )
    
    median_ploidy_not_noisy = median(data$mean_ploidy[data$is_noisy == F &
                                                          data$normalized_dimapd < opt$threshold_G1G2phase])
    p = data %>%
        ggplot(aes(mean_ploidy, normalized_dimapd, color = Type)) +
        geom_point(alpha = 0.3) +
        scale_color_manual(
            values = c(
                'G1/G2 cells' = "#005095",
                'S-phase' = "#78bd3e",
                'unknown cells' = "#dfbd31"  
            )
        ) +
        theme(legend.position = 'top', legend.title = element_blank()) +
        geom_vline(xintercept = median_ploidy_not_noisy)+
        xlab('Ploidy') + ylab('Variability')

    suppressMessages(ggsave(p,
           filename = paste0(opt$out, '/', opt$base_name, '_plot_th_', opt$threshold_Sphase,'-',opt$threshold_G1G2phase, '.pdf')))
}

if('Sphase_first_part' %in% names(opt) & 'Sphase_second_part' %in% names(opt)){
    distributions=tibble(A=opt$Sphase_first_part,
                         B=opt$Sphase_second_part)
}else{
    if(xor('Sphase_first_part' %in% names(opt) , 'Sphase_second_part' %in% names(opt))){
        warning('One of the correction factors has not been provided, the program will automatically correct the Sphase progression')
    }
# correct mean ploidy 
cl=makeCluster(opt$cores)
registerDoSNOW(cl)

# test multiple parameters to correct S phase
distributions = foreach(
    a = seq(0.95, 1, by = 0.001),
    .combine = 'rbind',
    .packages = c('tidyverse', 'LaplacesDemon', 'foreach')
) %dopar% {
    dist = foreach(
        b = seq(0.5, 0.55, by = 0.001),
        .combine = 'rbind',
        .packages = c('tidyverse', 'LaplacesDemon')
    ) %do% {
        #adjust S-phase based on a and b
        x = data %>%
            filter(Type=='S-phase')%>%
            mutate(
                corrected_mean_ploidy = ifelse(
                    mean_ploidy >= median_ploidy_not_noisy,
                    mean_ploidy / a,
                    mean_ploidy / b
                )
            )%>%
            dplyr::select(corrected_mean_ploidy)%>%pull()

            # are the data unimodal?
            if (is.unimodal(x)) {
                # d is the distance betwen theoretical center of the Sphase (G1 median ploidy *1.5)
                #and the average ploidy of the corrected S-phase
                # d is devided by sd(x) in order to select parametes that keep the distribution as wide as possible
                tibble(
                    A = a,
                    B = b,
                    d = 1/sd(x),
                    unimodal = T
                )
                
            }else{
                tibble(
                    A = a,
                    B = b,
                    d = sd(x),
                    unimodal = F
                )
            }

    }
    dist
}

stopCluster(cl)

#select the minimum value of d
if (nrow(distributions)==0){
   stop('S phase correction parameters could not be established, please provide manual ones') 
}else{
distributions=distributions%>%
    filter(unimodal==ifelse(any(unimodal==T),T,F))%>%
    filter(d==min(d))
if (nrow(distributions) > 1){
    stop('S phase correction parameters could not be established, please provide manual ones.') 
}
}
}

data = data %>%
    mutate(
        mean_ploidy_corrected = ifelse(
            Type=='S-phase' &
                mean_ploidy < median_ploidy_not_noisy,
            mean_ploidy / distributions$B,
            ifelse(
                Type=='S-phase' &
                    mean_ploidy > median_ploidy_not_noisy,
                mean_ploidy /  distributions$A,
                mean_ploidy
            )),
            Type = ifelse(
                Type=='S-phase' &
                    mean_ploidy < median_ploidy_not_noisy,
                'S-phase second part',
                ifelse(
                    Type=='S-phase' &
                        mean_ploidy > median_ploidy_not_noisy,
                    'S-phase first part',
                    Type
                )
        )
    )
    

p = data %>%
    ggplot(aes(mean_ploidy_corrected, normalized_dimapd, color = Type)) +
    geom_point(alpha = 0.3) +
    scale_color_manual(
        values = c(
            'G1/G2 cells' = '#005095',
            'S-phase first part' = '#78bd3e',
            'S-phase second part'='#83007e',
            'unknown cells' = '#dfbd31'
        )
    ) +
    theme(legend.position = 'top', legend.title = element_blank()) +
    geom_vline(xintercept = median_ploidy_not_noisy)+
    xlab('Ploidy') + ylab('Variability')+
    geom_density(data=data %>%filter(Type %in% c( 'S-phase second part','S-phase first part')),
                     aes(x=mean_ploidy_corrected,y=(..density..)), color="black")+
    scale_y_continuous(sec.axis = sec_axis(trans = ~.,name ='S-phase density distribution'))+
    theme(legend.position = 'top', legend.title = element_blank()) 

suppressMessages( ggsave(p,
                         filename = paste0(opt$out, '/', opt$base_name, '_plot_th_', opt$threshold_Sphase,'-',opt$threshold_G1G2phase, '_Sphase_corrected_',distributions$A,'-',distributions$B,'.pdf')))


tibble(
    threshold_Sphase= opt$threshold_Sphase,
    threshold_G1G2phase = opt$threshold_G1G2phase,
    Sphase_first_part=distributions$A,
    Sphase_second_part=distributions$B,
    RPM_TH=round(opt$min_n_reads*median_ploidy_not_noisy)
)%>%
    write_tsv(paste0(opt$out,opt$base_name, '_settings.txt'))

print('done')
