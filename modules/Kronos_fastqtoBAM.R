#!/usr/local/bin/Rscript --slave
#parse input

suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)

option_list = list(
    make_option(
        c("-l", "--fastq_list"),
        type = "character",
        default = NULL,
        help = "A table formatted in the following way: sam_file_basename\tFastq_1\tFastq_2(optional for PE sequencing). Compressed files are not allowed. Alternative to -1/-2/-b",
        metavar = "character"
    ),
    make_option(
        c("-O", "--one"),
        type = "character",
        default = NULL,
        help = "Fastq files, for multiple files they have to be separated by a comma. Compressed files are not allowed.  Alternative to -l",
        metavar = "character"
    ),
    make_option(
        c("-T", "--two"),
        type = "character",
        default = NULL,
        help = "Fastq files (for paired ends), for multiple files they have to be separated by a comma. Compressed files are not allowed.  Alternative to -l",
        metavar = "character"
    ),
    make_option(
        c("-b", "--sam_file_basename"),
        type = "character",
        default = NULL,
        help = "Sam file name, for multiple files they have to be separated by a comma.  Alternative to -l",
        metavar = "character"
    ),
    make_option(
        c("-i", "--index"),
        type = "character",
        action = 'store',
        help = "Bowtie 2 index",
        metavar = "character"
    ),
    make_option(
        c("-c", "--cores"),
        type = "integer",
        default = 1,
        action = 'store',
        help = "Number of cores to use. [default= %default]",
        metavar = "integer"
    ),
    make_option(
        c("-o", "--output_dir"),
        type = "character",
        default = 'output/',
        action = 'store',
        help = "Output folder. [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("--path_to_trim_galore"),
        type = "character",
        default = 'trim_galore',
        action = 'store',
        help = "Path to trim_galore",
        metavar = "character"
    ),
    make_option(
        c("--trim_galore_extra_option"),
        type = "character",
        default = '',
        action = 'store',
        help = "extra options for trim_galore",
        metavar = "character"
    ),
    make_option(
        c("--path_to_cutadapt"),
        type = "character",
        default = 'cutadapt',
        action = 'store',
        help = "Path to cutadapt",
        metavar = "character"
    ),
    make_option(
        c("--path_to_java"),
        type = "character",
        default = 'java',
        action = 'store',
        help = "Path to java",
        metavar = "character"
    ),
    make_option(
        c("--path_to_picard"),
        type = "character",
        default = 'picard',
        action = 'store',
        help = "Path to picard",
        metavar = "character"
    ),
    make_option(
        c("--keep_intermediate_files"),
        type = "logical",
        default = FALSE,
        action = "store_true",
        help = "keep trimmed fastq and sorted bam before depduplication",
        metavar = "logical"
    )
)

opt = parse_args(OptionParser(option_list = option_list),
                 convert_hyphens_to_underscores = T)

#load needed packages
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(Rbowtie2, quietly = TRUE))
suppressPackageStartupMessages(library(Rsamtools, quietly = TRUE))
suppressPackageStartupMessages(library(doSNOW, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))

options(scipen = 9999)


# check inputs
if (str_extract(opt$output_dir, '.$') != '/') {
    opt$output_dir = paste0(opt$output_dir, '/')
}


#findpaths
opt$path_to_trim_galore = Sys.which(opt$path_to_trim_galore)

if (opt$path_to_trim_galore == '') {
    stop('trim_galore was not found, please provide path')
}

opt$path_to_cutadapt = Sys.which(opt$path_to_cutadapt)

if (opt$path_to_cutadapt == '') {
    stop('cutadapt was not found, please provide path')
}

opt$path_to_java = Sys.which(opt$path_to_java)

if (opt$path_to_java == '') {
    stop('java was not found, please provide path.')
}

if (opt$path_to_picard == 'picard') {
    opt$path_to_picard = Sys.which(opt$path_to_picard)
}
if (opt$path_to_picard == '') {
    stop('picard was not found, please provide path.')
}


No_trim_galore_report=str_detect(string = opt$trim_galore_extra_option,pattern = '--no_report_file')

#create output directories
if(!No_trim_galore_report){
    system(
        paste0(
            'mkdir -p ',
            opt$output_dir,
            ' ',
            opt$output_dir,
            'trimmed ',
            opt$output_dir,
            'sorted_bam ',
            opt$output_dir,
            'delDup ',
            opt$output_dir,
            'MarkDuplicates_logs ',
            opt$output_dir,
            'Bowtie2_logs ',
            opt$output_dir,
            'TrimGalore_logs '
        )
    )    
}else{
    system(
        paste0(
            'mkdir -p ',
            opt$output_dir,
            ' ',
            opt$output_dir,
            'trimmed ',
            opt$output_dir,
            'sorted_bam ',
            opt$output_dir,
            'delDup ',
            opt$output_dir,
            'MarkDuplicates_logs ',
            opt$output_dir,
            'Bowtie2_logs '
        )
    )
}


#load check inputs
if ('fastq_list' %in% names(opt)) {
    fastqs = read_tsv(opt$fastq_list, col_names = F,col_types = cols())
    
    if (length(fastqs) == 2) {
        opt$sam_file_basename = fastqs$X1
        opt$one = fastqs$X2
        SE = T
    } else if (length(fastqs) == 3) {
        opt$two = fastqs$X3
        SE = F
        if (length(opt$one) != length(opt$two)) {
            stop('not all the fastq are paired')
        }
    } else{
        stop("Check fastq files table format. See script usage (--help)")
    }
    
} else{
    if ('one' %in% names(opt)) {
        opt$one = str_split(opt$one, ',')[[1]]
    } else{
        stop("No fastq files input. See script usage (--help)")
    }
    if ('two' %in% names(opt)) {
        SE = F
        opt$two = str_split(opt$two, ',')[[1]]
        
        if (length(opt$one) != length(opt$two)) {
            stop('not all the fastq files are paired')
        }
    } else{
        SE = T
    }
    
    if ('sam_file_basename' %in% names(opt)) {
        opt$sam_file_basename = str_split(opt$sam_file_basename, ',')[[1]]
    }

}

if (!"index" %in% names(opt)) {
    stop("bowtie2 index not provided. See script usage (--help)")
}

cl <-
    makeCluster(ifelse(opt$cores > length(opt$one), length(opt$one), opt$cores))
registerDoSNOW(cl)

tmp=foreach(i = 1:length(opt$one),.combine = 'rbind',.packages = c('tidyverse','Rbowtie2','Rsamtools')) %dopar% {
    if (SE) {
        #cut adapters
       system(
            paste0(
                opt$path_to_trim_galore,
                ' ',
                opt$one[i],
                ' --output_dir ',
                opt$output_dir,
                'trimmed/',
                ' --path_to_cutadapt ',
                opt$path_to_cutadapt,
                ' ',
                opt$trim_galore_extra_option
            )
        )
        
        #move trim galore report if it exists 
        if (!No_trim_galore_report) {
            system(
                paste0(
                    'mv ',
                    opt$output_dir,
                    'trimmed/',
                    basename(opt$one[i]),
                    '*.txt ',
                    opt$output_dir,
                    'TrimGalore_logs/'
                )
            )
        }
        file_basename = str_remove(basename(opt$one[i]), pattern = '.fastq$|.fq$')
        input_bowtie = list.files(paste0(opt$output_dir, 'trimmed/'), pattern = file_basename)
        
        if ('sam_file_basename' %in% names(opt)) {
            file_basename = opt$sam_file_basename[i]
        }
        
        #align with bowtie2
        suppressMessages(
            bowtie2(
                bt2Index = opt$index,
                samOutput = paste0(opt$output_dir, 'sorted_bam/', file_basename, '.sam'),
                seq1 = paste0(opt$output_dir, 'trimmed/', input_bowtie),
                ... = paste0('-k 1 --met 60 --met-file ',opt$output_dir,
                             'Bowtie2_logs/',file_basename,'.txt'),
                overwrite = TRUE
            )
        )
        
        #convert to BAM, sort and index 
        asBam(paste0(opt$output_dir, 'sorted_bam/', file_basename, '.sam'))
        system(paste0(
            'rm ',
            opt$output_dir,
            'sorted_bam/',
            file_basename,
            '.sam'
        ))
        
        #remove duplicates
        system(
            paste0(
                opt$path_to_java,
                ' -jar ',
                opt$path_to_picard,
                ' MarkDuplicates I=',
                opt$output_dir,
                'sorted_bam/',
                file_basename,
                '.bam O=',
                opt$output_dir,
                'delDup/',
                file_basename,
                "_delDupl.bam M=",
                opt$output_dir,
                'MarkDuplicates_logs/',
                file_basename,
                "_delDupl.log REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"
            )
        )
        
        
        if(!opt$keep_intermediate_files){
            system(
                paste('rm', 
                      paste0(opt$output_dir,'sorted_bam/',file_basename, '* '),
                      paste0(opt$output_dir, 'trimmed/', input_bowtie)
                ))
        }
        
    } else {
    
        #cut adapters
        system(
            paste0(
                opt$path_to_trim_galore,
                ' --paired ',
                opt$one[i],
                ' ',
                opt$two[i],
                ' --output_dir ',
                opt$output_dir,
                'trimmed/',
                ' --path_to_cutadapt ',
                opt$path_to_cutadapt,
                ' ',
                opt$trim_galore_extra_option, ' '
                
            )
        )
        

        #move trim_galore_reports if it exists into final folder        
        if (!No_trim_galore_report) {
            system(
                paste0(
                    'mv ',
                    opt$output_dir,
                    'trimmed/',
                    basename(opt$one[i]),
                    '*.txt ',
                    opt$output_dir,
                    'TrimGalore_logs/'
                )
            )
            system(
                paste0(
                    'mv ',
                    opt$output_dir,
                    'trimmed/',
                    basename(opt$two[i]),
                    '*.txt ',
                    opt$output_dir,
                    'TrimGalore_logs/'
                )
            )
        }
        

        file_basename_one = str_remove(basename(opt$one[i]),  pattern = '.fastq$|.fq$')
        file_basename_two = str_remove(basename(opt$two[i]),  pattern = '.fastq$|.fq$')
        input_bowtie_1 = list.files(paste0(opt$output_dir, 'trimmed/'), pattern = file_basename_one)
        input_bowtie_2 = list.files(paste0(opt$output_dir, 'trimmed/'), pattern = file_basename_two)
        
        if ('sam_file_basename' %in% names(opt)) {
            file_basename = opt$sam_file_basename[i]
        } else{
            file_basename = file_basename_one
        }
        
        #align with bowtie2
        suppressMessages(
            bowtie2(
                bt2Index = opt$index,
                samOutput = paste0(opt$output_dir, 'sorted_bam/', file_basename, '.sam'),
                seq1 = paste0(opt$output_dir, 'trimmed/', input_bowtie_1),
                seq2 = paste0(opt$output_dir, 'trimmed/', input_bowtie_2),
                ... = paste0('-k 1 --met 60 --met-file ',opt$output_dir,
                             'Bowtie2_logs/',file_basename,'.txt'),
                overwrite = TRUE
            )
        )
        
        #convert to BAM, sort and index 
        asBam(paste0(opt$output_dir, 'sorted_bam/', file_basename, '.sam'))
        system(paste0(
            'rm ',
            opt$output_dir,
            'sorted_bam/',
            file_basename,
            '.sam'
        ))
        
        #remove duplicates
        system(
            paste0(
                opt$path_to_java,
                ' -jar ',
                opt$path_to_picard,
                ' MarkDuplicates I=',
                opt$output_dir,
                'sorted_bam/',
                file_basename,
                '.bam O=',
                opt$output_dir,
                'delDup/',
                file_basename,
                "_delDupl.bam M=",
                opt$output_dir,
                'MarkDuplicates_logs/',
                file_basename,
                "_delDupl.log REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"
            )
        )
        if(!opt$keep_intermediate_files){
            system(
                paste('rm', 
                     paste0(opt$output_dir,'sorted_bam/',file_basename, '* '),
                paste0(opt$output_dir, 'trimmed/', input_bowtie_1),
                paste0(opt$output_dir, 'trimmed/', input_bowtie_2)
            ))
        }
    }
    
    
    read_tsv(paste0(opt$output_dir,
             'MarkDuplicates_logs/',file_basename,'_delDupl.log'), skip = 6)%>%
        mutate(LIBRARY=file_basename)
    
    
    
}

if(!opt$keep_intermediate_files){
    system(
        paste('rm -r', 
              paste0(opt$output_dir,'sorted_bam/'),
              paste0(opt$output_dir, 'trimmed/')
        ))
}

stopCluster(cl)

tmp%>%write_tsv(paste(opt$output_dir,'Single_cells_reads_info.tsv'), append = T) 
print('done')
