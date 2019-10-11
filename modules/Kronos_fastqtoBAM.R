#!/usr/local/bin/Rscript --slave
#parse input

if (!suppressPackageStartupMessages(require(optparse, quietly = TRUE))) {
    install.packages("optparse", quiet = T)
    suppressPackageStartupMessages(library(optparse, quietly = TRUE))
}

options(stringsAsFactors = FALSE)

option_list = list(
    make_option(
        c("-1", "--one"),
        type = "character",
        default = NULL,
        help = "fastq files",
        metavar = "character"
    ),make_option(
        c("-2", "--two"),
        type = "character",
        default = NULL,
        help = "fastq files (for paired ends)",
        metavar = "character"
    ),make_option(
        c("-b", "--sam_file_basename"),
        type = "character",
        default = NULL,
        help = "sam file name",
        metavar = "character"
    ),
    make_option(
        c("-i", "--index"),
        type = "character",
        action = 'store',
        help = "bowtie 2 index",
        metavar = "character"
    ),
    make_option(
        c("-c", "--cores"),
        type = "integer",
        default = 3,
        action = 'store',
        help = "Number of cores to use. [default= %default]",
        metavar = "integer"
    ),
    make_option(
        c("-o", "--output_dir"),
        type = "character",
        default = 'output/',
        action = 'store',
        help = "output folder. [default= %default]",
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
    )
    
)

opt = parse_args(OptionParser(option_list = option_list),convert_hyphens_to_underscores = T)

#load needed packages
if (!suppressPackageStartupMessages(require(tidyverse, quietly = TRUE))) {
    install.packages("tidyverse", quiet = T)
    suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
}

if (!suppressPackageStartupMessages(require(DescTools, quietly = TRUE))) {
    install.packages("DescTools", quiet = T)
    suppressPackageStartupMessages(library(DescTools, quietly = TRUE))
}

if(!suppressPackageStartupMessages(require(Rbowtie2, quietly = TRUE))){
    BiocManager::install('Rbowtie2',quiet = T)
    suppressPackageStartupMessages( library(Rbowtie2, quietly = TRUE))
}

if(!suppressPackageStartupMessages(require(Rsamtools, quietly = TRUE))){
    install.packages("Rsamtools",quiet = T)
    suppressPackageStartupMessages( library(Rsamtools))
}
options(scipen = 9999)


# check imputs 
if(str_extract(opt$output_dir,'.$')!='/'){
    opt$output_dir=paste0(opt$output_dir,'/')
}
if(str_extract(opt$path_to_fastq,'.$')!='/'){
    opt$path_to_fastq=paste0(opt$path_to_fastq,'/')
}

#findpaths
opt$path_to_trim_galore=Sys.which(opt$path_to_trim_galore)

if(opt$path_to_trim_galore==''){
    stop('trim_galore was not found, please provide path')
}

opt$path_to_cutadapt=Sys.which(opt$path_to_cutadapt)

if(opt$path_to_cutadapt==''){
    stop('cutadapt was not found, please provide path.')
}

opt$path_to_java=Sys.which(opt$path_to_java)

if(opt$path_to_java==''){
    stop('java was not found, please provide path.')
}


opt$path_to_picard=Sys.which(opt$path_to_picard)

if(opt$path_to_picard==''){
    stop('picard was not found, please provide path.')
}

#create output directories
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
        'delDup/logs '
    )
)

if ('one' %in% names(opt) &
    !'two' %in% names(opt)) {
    if (!"index" %in% names(opt)) {
        stop("bowtie2 index not provided. See script usage (--help)")
    }
    
    #cut adapters
    system(
        paste0(
            opt$path_to_trim_galore,
            ' ',
            opt$path_to_fastq,
            opt$one,
            ' --output_dir ',
            opt$output_dir,
            'trimmed/',
            ' --path_to_cutadapt ',
            opt$path_to_cutadapt,
            ' --no_report_file'
            
        )
    )
    
    file_basename = str_remove(basename(opt$one), pattern = '.fastq|.fastq.gz|.fq|.fq.gz')
    input_bowtie = list.files(paste0(opt$output_dir, 'trimmed/'), pattern = file_basename)
    
    if ('sam_file_basename' %in% names(opt)) {
        file_basename = opt$sam_file_basename
    }
    
    #align with bowtie2
    suppressMessages(
        bowtie2(
            bt2Index = opt$index,
            samOutput = paste0(opt$output_dir, 'sorted_bam/', file_basename, '.sam'),
            seq1 = paste0(opt$output_dir, 'trimmed/', input_bowtie),
            ... = paste0('-k 1 '),
            overwrite = TRUE
        )
    )
    
    asBam(paste0(opt$output_dir, 'sorted_bam/', file_basename, '.sam'))
    system(paste0('rm ', opt$output_dir, 'sorted_bam/', file_basename, '.sam'))
    
    
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
            'delDup/logs/',
            file_basename,
            "_delDupl.log REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"
        )
    )
    
} else if (c('one', 'two') %in% names(opt)) {
    if (!"index" %in% names(opt)) {
        stop("bowtie2 index not provided. See script usage (--help)")
    }
    
    #cut adapters
    system(
        paste0(
            opt$path_to_trim_galore,
            ' ',
            opt$path_to_fastq,
            list$X1[i],
            ' --output_dir ',
            opt$output_dir,
            'trimmed/',
            ' --path_to_cutadapt ',
            opt$path_to_cutadapt,
            ' --no_report_file'
            
        )
    )
    
    
    file_basename_one = str_remove(basename(opt$one), pattern = '.fastq|.fastq.gz|.fq|.fq.gz')
    file_basename_two = str_remove(basename(opt$two), pattern = '.fastq|.fastq.gz|.fq|.fq.gz')
    input_bowtie_1 = list.files(paste0(opt$output_dir, 'trimmed/'), pattern = file_basename_one)
    input_bowtie_2 = list.files(paste0(opt$output_dir, 'trimmed/'), pattern = file_basename_two)
    
    if ('sam_file_basename' %in% names(opt)) {
        file_basename = opt$sam_file_basename
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
            ... = paste0('-k 1 '),
            overwrite = TRUE
        )
    )
    
    asBam(paste0(opt$output_dir, 'sorted_bam/', file_basename, '.sam'))
    system(paste0('rm ', opt$output_dir, 'sorted_bam/', file_basename, '.sam'))
    
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
            'delDup/logs/',
            file_basename,
            "_delDupl.log REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"
        )
    )
    
} else{
    stop("Check list of fastq files format. See script usage (--help)")
}

