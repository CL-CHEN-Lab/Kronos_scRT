#!/usr/local/bin/Rscript
#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)

option_list = list(
    make_option(
        c("-R", "--RefGenome"),
        type = "character",
        default = NULL,
        help = "Fasta file of genome of interst",
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
        c("-s", "--reads_size"),
        type = "integer",
        default = 40,
        action = 'store',
        help = "Length of the simulated reads. [default= %default bp]",
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
        c("-i", "--index"),
        type = "character",
        action = 'store',
        help = "Bowtie 2 index",
        metavar = "character"
    ),
    make_option(
        c( "--paired_ends"),
        type = "logical",
        action = 'store_true',
        help = "Generates paired ends reads [default: %default]",
        metavar = "logical",
        default = F
    ),
    make_option(
        c( "--insert_size"),
        type = "integer",
        action = 'store',
        help = "Insert size if paired end option is used. [default: %default]",
        metavar = "integer",
        default = '200'
    ),
    make_option(
        c("--bin_size"),
        type = "character",
        default = '20Kb',
        action = 'store',
        help = "Bins size. [default= %default ]",
        metavar = "character"
    ),
    make_option(
        c("-d","--dir_indexed_bam"),
        type = "character",
        action = 'store',
        help = "If provided parameters will be automatically estimated from the data.",
        metavar = "character"
    ),
    make_option(
      c("-u","--upper_mappability_th"),
      type = "double",
      default = 1.5,
      action = 'store',
      help = "Maximum mappability for a bin to be considered in the analisys  [default= %default]",
      metavar = "double"
    ),
    make_option(
      c("-l","--lower_mappability_th"),
      type = "double",
      action = 'store',
      default = 0.8,
      help = "Minimum mappability for a bin to be considered in the analisys  [default= %default]",
      metavar = "double"
    ),
    make_option(
      c("-B","--black_list"),
      type = "character",
      action = 'store',
      help = "Regions to ignore",
      metavar = "character"
    ),
    make_option(
        c("-x","--coverage"),
        type = "character",
        action = 'store',
        help = "Coverage for simulated genome. [default= %default]",
        default = '1x',
        metavar = "character"
    ),
    make_option(
        c("-e","--errorRate"),
        type = "character",
        action = 'store',
        help = "Simulated sequencing error rate (%) [default= %default]",
        default = "0.1%",
        metavar = "character"
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
      default = "1:22,X,Y",
      metavar = "character"
    )
)

opt = parse_args(OptionParser(option_list = option_list),convert_hyphens_to_underscores = T)

#load needed packages

suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))
suppressPackageStartupMessages( library(doSNOW, quietly = TRUE))
suppressPackageStartupMessages( library(Biostrings, quietly = TRUE))
suppressPackageStartupMessages( library(Rbowtie2, quietly = TRUE))
suppressPackageStartupMessages( library(Rsamtools, quietly = TRUE))
suppressPackageStartupMessages( library(GenomicRanges, quietly = TRUE))

options(scipen = 9999)

# create output directory
if (str_extract(opt$output_dir,'.$')!='/'){
    opt$output_dir=paste0(opt$output_dir,'/')
}

system(paste0('mkdir -p ', opt$output_dir))

# check inputs 
if(!"RefGenome" %in% names(opt)){
    stop("Fasta file not provided. See script usage (--help)")
}

if(!"index" %in% names(opt)){
    stop("Bowtie2 indexed genome not provided. See script usage (--help)")
}

opt$coverage=tryCatch(expr = as.numeric(str_remove(opt$coverage,'[Xx]')),
                      error = function(x) tryCatch(expr = as.numeric(opt$coverage),
                      error= stop('Wrong coverage format')))

# convert binsize to numeric

extract_unit=str_extract(opt$bin_size,pattern = '.{2}$')

if(grepl(x = extract_unit,pattern =  '[0-9][0-9]')){
  n_of_0=str_length(str_extract(opt$bin_size,'0{1,10}$'))
  BS=case_when(
    is.na(n_of_0) ~ paste0(opt$bin_size,'bp'),
    n_of_0 < 3 ~ paste0(opt$bin_size,'bp'),
    n_of_0 < 6 ~ paste0(str_remove(opt$bin_size,'0{3}$'),'Kb'),
    n_of_0 >= 6 ~paste0(str_remove(opt$bin_size,'0{6}$'),'Mp'))
}else{
  BS=opt$bin_size
}

opt$bin_size = as.numeric(str_remove(opt$bin_size, "[Bb][Pp]|[Kk][Bb]|[Mm][Bb]")) * case_when(
  grepl(x =extract_unit,pattern =  '[Kk][Bb]') ~ 1000,
  grepl(x =extract_unit, pattern = '[Mm][Bb]') ~ 1000000,
  grepl(x =extract_unit, pattern = '[Bp][Pp]') ~ 1,
  grepl(x = extract_unit,pattern =  '[0-9][0-9]') ~ 1
)

if(is.na(opt$bin_size)){
  stop('binsize have an incorrect format')
}

#loading reference fa
reference=readDNAStringSet(opt$RefGenome)

#exstimate paramenters
cl=makeCluster(opt$cores)
registerDoSNOW(cl)

if ('dir_indexed_bam' %in% names(opt)){
    if (str_extract(opt$dir_indexed_bam,'.$')!='/'){
        opt$dir_indexed_bam=paste0(opt$dir_indexed_bam,'/')
    }
    #sample 30 files (if available) to exstimate parameters
    list=list.files(opt$dir_indexed_bam,pattern = 'bam$')
    
    if(length(list) > 30 ){
        list=sample(list,30) 
    }
   
    parameters=foreach(i=list,.combine = 'rbind',.packages = 'Rsamtools')%dopar%{
        sapply(scanBam(paste0(opt$dir_indexed_bam,i),param=ScanBamParam(what=c('isize','qwidth')))[[1]],
                function(x) median(abs(x),na.rm = T))
    }
    
    parameters=as_tibble(parameters)%>%
        summarise(qwidth=round(median(qwidth)),
                  isize=round(median(isize)))
    
    if(parameters$isize!=0){
        opt$paired_ends=T
        opt$insert_size=parameters$isize
        opt$reads_size=parameters$qwidth
    }else{
        opt$paired_ends=F
        opt$reads_size=parameters$qwidth
    }
  
}

# select chrs of interest
# convert string into range
Convert_to_range = Vectorize(function(x){
  if (str_detect(x, ':')) {
    x = str_split(x, ':')[[1]]
    return(as.numeric(x[1]):as.numeric(x[2]))
  } else{
    return(x)
  }
})

#select chrs
chr_list = paste0(ifelse(opt$chr_prefix=='none','',opt$chr_prefix), unlist(Convert_to_range(str_split(opt$chr_range,',')[[1]])))
chr_list = chr_list[chr_list %in% unique(names(reference))]

# Identify regins in which the sequence is known
find_known_sequences = function(x) {
    library(tidyverse)
    y = str_locate_all(x, '.[TACG]+')
    return(tibble(start = y[[1]][, 1], end = y[[1]][, 2]))
}
#recover reads and mutate them
recover_and_mutate = function(simulated_reads, Chr_reference,errorRate=opt$errorRate) {
    library(tidyverse)
    # simulate mutations 0.1 % rate
    errorRate=as.numeric(str_remove(string = errorRate ,pattern ='%'))/100
    mutate = sample(1:str_length(Chr_reference),
                    errorRate*str_length(Chr_reference))

    #convert characters into int
    REF_mutated=utf8ToInt(as.character(Chr_reference[[1]]))

    # applies mutation
    mutate_sequence=Vectorize(function(Base){
        
        bases=c(65,67,84,71)
        bases=bases[bases != Base]
        
        return(sample(bases,1))
        
    },vectorize.args = 'Base')
    
    #mutation
    REF_mutated[mutate]= mutate_sequence(REF_mutated[mutate])
    
    #convert back to string 
    REF_mutated=intToUtf8(REF_mutated)
    
    #recover reads sequences
    simulated_reads = tibble(
        reads = str_sub(
            REF_mutated,
            start = simulated_reads$start,
            end = simulated_reads$end - 1
        ),
        order = simulated_reads$order
    )%>%
        mutate(reads=str_remove(reads,'N{1,1000}$'),
               reads=str_remove(reads,'^N{1,1000}'))
    
    return(simulated_reads)
}
#reshape reads for the fastq file
reshape_and_save = function(simulated_reads, file,Chr,x) {
    simulated_reads %>%
        arrange(order)%>%
        `colnames<-`(c('2','order') )%>%
        mutate(
            n=str_count(`2`),
            `1` = paste0('@read', 1:n(),Chr,'round',x),
            `3` = '+',
            `4` = unlist(lapply(n, function(x) paste0(rep('D', x), collapse = '')))
        ) %>%
        select(-n)%>%
        gather('pos', 'towrite', -order) %>%
        arrange(order, pos) %>%
        select(towrite) %>%
        write_delim(path = file,
                    col_names = F,append = T)
    return(0)
}

#calculate reverse complement for PE
rev_com = function(x) {
    library(tidyverse)
    dict = list(
        A = 'T',
        `T` = 'A',
        C = 'G',
        G = 'C',
        N = 'N'
    )
    x = str_extract_all(x, '.')
    x = unlist(lapply(x , function(x)
        paste0(
            sapply(rev(x), function(x)
                dict[[x]], simplify = T),
            collapse = ''
        )))
    return(x)
}


genome.Chromsizes = foreach(
    Chr = chr_list,
    .combine = 'rbind',
    .packages = c('Biostrings', 'foreach', 'tidyverse')
)%dopar%{
    #genome size
    genome.Chromsizes = tibble(chr = Chr,
                               size = width(reference[Chr]))
    
    position = find_known_sequences(reference[Chr])
    
    #look for seeds
    if (opt$paired_ends) {
        #initialize simulated reads
        size =  opt$insert_size
        
    }else{
        
        size = opt$reads_size 
        
        }
    
    
    
    tmp=foreach(
        variability=round(seq(-size,size,by = 2*size/(opt$coverage+1)))[(1:opt$coverage)+1],
    .packages = c('Biostrings', 'foreach', 'tidyverse')
) %do% {
    
    #look for seeds
    if (opt$paired_ends) {
        
        #initialize simulated reads
        position = position %>%
            filter(end - start > size)
        
        simulated_reads_1 = foreach(i = 1:length(position$start),
                                    .combine = 'c') %do% {
                                        seq(position$start[i]+variability, position$end[i], by = size)
                                    }
        simulated_reads_1 = tibble(start = unlist(simulated_reads_1),
                                   end = start + opt$reads_size) %>%
            mutate(order =  row_number())
        simulated_reads_2 = simulated_reads_1 %>%
            mutate(end = start + opt$insert_size,
                   start= end - opt$reads_size)
        
    } else{
        
        position = position %>%
            filter(end - start > size)
        #initialize simulated reads
        simulated_reads = foreach(i = 1:length(position$start)) %do% {
            seq(position$start[i]+variability, position$end[i], by = opt$reads_size)
        }
        simulated_reads = tibble(start = unlist(simulated_reads),
                                 end = start + opt$reads_size) %>%
            mutate(order =  row_number())
    }

    if (opt$paired_ends) {
        #recover strings reads and mutate some of them
        simulated_reads_1 = recover_and_mutate(simulated_reads = simulated_reads_1,
                                               Chr_reference = reference[Chr])
        simulated_reads_2 = recover_and_mutate(simulated_reads = simulated_reads_2,
                                               Chr_reference = reference[Chr])

        
        simulated_reads_2 = simulated_reads_2 %>%
            mutate(reads = rev_com(reads))
    } else{
        #recover strings reads and mutate some of them
        simulated_reads = recover_and_mutate(simulated_reads = simulated_reads,
                                             Chr_reference = reference[Chr])
    }
    
    if (opt$paired_ends) {
        #save fastq files
        reshape_and_save(
            simulated_reads_1,
            file = paste0(
                opt$output_dir,
                basename(opt$index),'_',
                Chr,
                '_simulated_reads_1.fq'
            ),
            Chr=Chr,
            x = variability
        )
        reshape_and_save(
            simulated_reads_2,
            file = paste0(
                opt$output_dir,
                basename(opt$index),'_',
                Chr,
                '_simulated_reads_2.fq'
            ),
            Chr=Chr,
            x = variability
        )
        
    } else{
        #save fastq file
        reshape_and_save(
            simulated_reads,
            file = paste0(
                opt$output_dir,
                basename(opt$index),'_',
                Chr,
                '_simulated_reads.fq'
            ),
            Chr=Chr,
            x = variability
        )    #remove from memory
        rm('simulated_reads')
    }
    variability
}
genome.Chromsizes
    }

stopCluster(cl)

if(opt$paired_ends){
    #merge files
    system(paste0('cat ',
        opt$output_dir,basename(opt$index),
        '_*_simulated_reads_2.fq > ', opt$output_dir,
        basename(opt$index),
        '_simulated_reads_2.fq; cat ',
        opt$output_dir,basename(opt$index),
        '_*_simulated_reads_1.fq > ', opt$output_dir,
        basename(opt$index),
        '_simulated_reads_1.fq'
    ))
    system(paste0('rm ',opt$output_dir, basename(opt$index),'_*_simulated_reads_1.fq'))
    system(paste0('rm ',opt$output_dir, basename(opt$index),'_*_simulated_reads_2.fq'))
    #align with bowtie2
    bowtie2(
        bt2Index = opt$index,
        samOutput = paste0(opt$output_dir, basename(opt$index), '_simulated_reads.sam'),
        seq1 = paste0(opt$output_dir, basename(opt$index), '_simulated_reads_1.fq'),
        seq2 = paste0(opt$output_dir, basename(opt$index), '_simulated_reads_2.fq'),
        ... = paste0('--phred33 --ignore-quals -p ', opt$cores),
        overwrite=TRUE
    )
    
    # remove from hd simulated_reads.fq
    system(paste0('rm ',opt$output_dir, basename(opt$index),'_simulated_reads_1.fq'))
    system(paste0('rm ',opt$output_dir, basename(opt$index),'_simulated_reads_2.fq'))
    
}else{
    #merge files
    system(paste0('cat ',
                  opt$output_dir,basename(opt$index),
                  '_*_simulated_reads.fq > ', opt$output_dir,
                  basename(opt$index),
                  '_simulated_reads.fq'
    ))
    
    system(paste0('rm ',opt$output_dir, basename(opt$index),'_*_simulated_reads.fq'))
    #align with bowtie2
    suppressMessages(bowtie2(
        bt2Index = opt$index,
        samOutput = paste0(opt$output_dir, basename(opt$index), '_simulated_reads.sam'),
        seq1 = paste0(opt$output_dir, basename(opt$index), '_simulated_reads.fq'),
        ... = paste0('--phred33 --ignore-quals -p ', opt$cores),
        overwrite=TRUE
    ))
    
    # remove from hd simulated_reads.fq
    system(paste0('rm ',opt$output_dir, basename(opt$index),'_simulated_reads.fq'))
}

# calculate bins
dir.bam<-asBam( paste0(opt$output_dir, basename(opt$index), '_simulated_reads.sam'))
system(paste0('rm ',opt$output_dir, basename(opt$index),'_simulated_reads.sam'))

if(opt$paired_ends){
    param1 <- ScanBamParam(what=c('rname','pos','isize','mapq'),
                           flag=scanBamFlag(hasUnmappedMate = T,isUnmappedQuery = F))
    param2 <- ScanBamParam(what=c('rname','pos','isize','mapq','mrnm'),
                           flag=scanBamFlag(isPaired = T,isUnmappedQuery = F))
    bins = rbind(
        as.data.frame(scanBam(dir.bam, param = param1)) %>%
            filter(mapq >= 30) %>%
            select('rname', 'pos') %>%
            `colnames<-`(c('chr', 'pos')) %>%
            mutate(read = 1),
        as.data.frame(scanBam(dir.bam, param = param2)) %>%
            filter(mapq >= 30) %>%
            mutate(read = ifelse(rname == mrnm &
                                     abs(isize) < opt$bin_size, 0.5, 1)) %>%
            select('rname', 'pos', 'read') %>%
            `colnames<-`(c('chr', 'pos', 'read'))
    )%>%
        drop_na()
    #parameter used to estiamte mappability th
    theoretical_reads = opt$bin_size/(opt$insert_size)

}else{
    param <- ScanBamParam(what=c('rname','pos','mapq'),
                          flag=scanBamFlag(isUnmappedQuery = F))
    
    bins =as.data.frame(scanBam(dir.bam,param=param))%>%
        filter(mapq >= 30)%>%
        mutate(read=1)%>%
        select('rname', 'pos', 'read')%>%
        `colnames<-`(c('chr', 'pos','read'))
    
    #parameter used to estiamte mappability th
    theoretical_reads = opt$bin_size/opt$reads_size
    
}

bins = foreach (Chr = genome.Chromsizes$chr,
                .combine = 'rbind',
                .packages = 'tidyverse') %do% {
                    size = genome.Chromsizes$size[genome.Chromsizes$chr == Chr]
                    bins_chr = tibble(chr = Chr,
                                      start = seq(0, size, by = opt$bin_size)) %>%
                        mutate(end = lead(start, n = 1, default =  size))
                    ## calculate reads per bin Selected reads
                    reads_proper <- bins$pos[bins$chr == Chr &
                                                 bins$read == 0.5 ]
                    reads_notproper <- bins$pos[bins$chr == Chr &
                                                    bins$read == 1]
                    if (length(reads_proper)!=0) {
                        reads_proper[reads_proper <= 0] <- 1
                        reads_proper <-
                            hist(reads_proper,
                                 breaks =  c(1, bins_chr$end),
                                 plot = F)
                        reads_proper <-
                            reads_proper$counts / 2
                    }
                    if (length(reads_notproper)!=0) {
                        reads_notproper[reads_notproper <= 0] <- 1
                        reads_notproper <-
                            hist(reads_notproper,
                                 breaks =  c(1, bins_chr$end),
                                 plot = F)
                        reads_notproper <-
                            reads_notproper$counts
                    }
                    if (length(reads_proper)!=0 &
                        length(reads_notproper)!=0) {
                        reads = reads_notproper + reads_proper
                    } else if (length(reads_proper)!=0 &
                               length(reads_notproper)==0) {
                        reads = reads_proper
                    } else if (length(reads_proper)==0 &
                               length(reads_notproper)!=0) {
                        reads = reads_notproper
                    } else{
                        reads = 0
                    }
                    ## Concatenate
                    bins_chr %>%
                        mutate(reads = reads)
                }

bins=bins%>%
    mutate(mappability=reads/(opt$coverage*theoretical_reads),
           mappability_th=ifelse(
               mappability >= opt$lower_mappability_th &
                   mappability <= opt$upper_mappability_th ,T,F
           ))%>%
    group_by(chr)%>%
    select(chr,start,end,mappability,mappability_th)

#delete file
system(paste0('rm ',opt$output_dir, basename(opt$index), '_simulated_reads.bam*'))

#calculate gc % peer bin
cl=makeCluster(opt$cores)
registerDoSNOW(cl)

bins=foreach(i=unique(bins$chr),.combine = 'rbind', .packages =c('Biostrings','tidyverse') )%dopar%{
    bins%>%
        filter(chr==i)%>%
        mutate(
            seq=str_sub(string =  reference[names(reference)==i],start=start+1, end=end),
            gc_frequency=str_count(seq,'G|C')/str_length(seq)
        )%>%
        select(-seq)
}

stopCluster(cl)

bins=bins %>%
    mutate(type=ifelse(opt$paired_ends,'PE','SE'))

if('black_list' %in% names(opt)){
  bl=read_tsv(opt$black_list,col_names = c('chr','start','end'),col_types = cols(chr='c'))%>%
    makeGRangesFromDataFrame()
  
  tbins=bins%>%makeGRangesFromDataFrame()
  
  hits=findOverlaps(query = tbins, subject = bl)
  overlaps=pintersect(bl[subjectHits(hits)],tbins[queryHits(hits)])
  overlaps=overlaps[width(overlaps) > opt$reads_size]
  bins=overlaps%>%as_tibble()%>%rename(seqnames='chr')%>%
    dplyr::select(chr,start,end,hit)%>%
    right_join(bins, by = c("chr", "start", "end"))
  
  bins=bins%>%
    mutate(mappability_th=ifelse(!is.na(hit),F,mappability_th))%>%
    dplyr::select(-hit)%>%
    arrange(chr,start)
  
  }


#write bisns with info
write_tsv(bins, paste0(opt$output_dir, basename(opt$index),'_',
                       BS, '_bins_',
                       opt$coverage,'X_coverage_',
                       opt$reads_size,'bp_reads_',
                       ifelse(opt$paired_ends,paste0('PE_',opt$insert_size,'bp_InsertSize')
                              ,'SE_'),
                       ifelse('black_list' %in% names(opt),'blacklisted_',''),
                              paste0('error_rate_',opt$errorRate,'_min_mappability_',opt$lower_mappability_th,
                                     '_max_mappability_',opt$upper_mappability_th),'.tsv'))

print('done')
