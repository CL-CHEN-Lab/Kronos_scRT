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
        help = "Per cell stat file path",
        metavar = "character"
    ),
    make_option(
        c("-W", "--whoSwho"),
        type = "character",
        default = NULL,
        help = "Who's who file path ( tsv file with header: Cell \t Phase)",
        metavar = "character"
    ),
    make_option(
        c("-o", "--out"),
        type = "character",
        default = "./output",
        help = "Output directory [default= %default]",
        metavar = "character"
    )
)

opt = parse_args( OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))


#check inputs
if (!'file' %in% names(opt)) {
    stop("Per cell stat file must be provided. See script usage (--help)")
}
if (!'whoSwho' %in% names(opt)) {
    stop("who's who file must be provided. See script usage (--help)")
}
#does exist/ right format function
does_exist_right_format = Vectorize(function(File, delim = '\t', columns_to_check,message='does not have the proper format') {
    #checks if the file exist
    if (!file.exists(File)) {
        return(paste(File, 'does not exist'))
    } else{
        # if columns_to_check is numeric check the number of colums
        if(is.numeric(columns_to_check)){
            if (ncol(tryCatch(
                expr =  read_delim(
                    File,
                    col_types = cols(),
                    n_max = 0,
                    delim = delim
                ),
                error = function(x)
                    tibble()
            ))!=columns_to_check) {
                return(paste(File, message,'\n'))
            } else{
                return('')
            }
            
            # if columns_to_check is not numeric check columns names    
        }else{
            #checks if it has the right format
            if (!all(columns_to_check %in% colnames(tryCatch(
                expr =  read_delim(
                    File,
                    col_types = cols(),
                    n_max = 0,
                    delim = delim
                ),
                error = function(x)
                    tibble()
            )))) {
                return(paste(File, message,'\n'))
            } else{
                return('')
            }
        }
    }
    
}, vectorize.args = 'File')

#check per cell files 
results=paste(does_exist_right_format(File=opt$file,delim = ',',columns_to_check=c('Cell',
                                                                                   'normalized_dimapd',
                                                                                   'mean_ploidy',
                                                                                   'ploidy_confidence',
                                                                                   'is_high_dimapd',
                                                                                   'is_noisy',
                                                                                   'coverage_per_1Mbp'),
                                      message = ',provided as a per cell file, does not have the right format'),collapse = '')
if(results!='') {
    stop(results)
}

results=paste(does_exist_right_format(File=opt$file,delim = ',',columns_to_check=c('Cell',
                                                                                   'Phase'),
                                      message = ',provided as a staging file, does not have the right format'),collapse = '')
if(results!='') {
    stop(results)
}


#create output directory
if (str_extract(opt$out,'.$')!='/'){
    opt$out=paste0(opt$out,'/')
}

system(paste0('mkdir -p ', opt$out))

#load data
data<-inner_join(read_csv(opt$file,
               col_types = cols()),
                   read_tsv(opt$whoSwho,
                        col_types = cols()), by = "Cell")

# write data
data%>%
    mutate(Phase=str_to_upper(Phase))%>%
    mutate(is_high_dimapd=ifelse(Phase=='S',T,F),
           is_noisy=ifelse(Phase=='S',T,F))%>%
    dplyr::select(-Phase)%>%
    write_csv(paste0(opt$out,'phased_',basename(opt$file)))

print('done')
