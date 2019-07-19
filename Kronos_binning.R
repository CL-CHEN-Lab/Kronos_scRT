#!/usr/local/bin/Rscript
#parse input
if(!require(optparse, quietly = TRUE)){
    install.packages("optparse",quiet = T)
    library(optparse, quietly = TRUE)
}

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
        help = "lengh of the simulated reads. [default= %default bp]",
        metavar = "integer"
    ),
    make_option(
        c("-n", "--reads_number"),
        type = "integer",
        default = 2000000,
        action = 'store',
        help = "number of the simulated reads per bin. [default= %default bp]",
        metavar = "integer"
    ),
    make_option(
        c("-b", "--reads_per_bin"),
        type = "integer",
        default = 200,
        action = 'store',
        help = "number of the simulated reads if variable bins size is used. [default= %default bp]",
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
        c("-i", "--index"),
        type = "character",
        action = 'store',
        help = "bowtie 2 index",
        metavar = "character"
    ),
    make_option(
        c( "--paired_ends"),
        type = "logical",
        action = 'store_true',
        help = "generates paired ends reads [default: %default]",
        metavar = "logical",
        default = F
    ),
    make_option(
        c( "--insert_size"),
        type = "character",
        action = 'store',
        help = "min:max insert size if paired end option is used. [default: %default]",
        metavar = "interval",
        default = '100:300'
    ),
    make_option(
        c( "--variable_bins"),
        type = "logical",
        action = 'store_true',
        help = "generates variable bins depending on mapped simulated reads [default: %default]",
        metavar = "logical",
        default = F
    ),
    make_option(
        c("--bin_size"),
        type = "integer",
        default = 20000,
        action = 'store',
        help = "Bins size if --variable_bins is not used. [default= %default bp]",
        metavar = "integer"
    )
)

opt = parse_args(OptionParser(option_list = option_list),convert_hyphens_to_underscores = T)
system(paste0('mkdir -p ', opt$output_dir))

#load needed packages
if (!require(BiocManager, quietly = TRUE)){
    install.packages("BiocManager",quiet = T)
}

if(!require(tidyverse, quietly = TRUE)){
    install.packages("tidyverse",quiet = T)
    library(tidyverse, quietly = TRUE )
}

if(!require(DescTools, quietly = TRUE)){
    install.packages("DescTools",quiet = T)
    library(DescTools, quietly = TRUE )
}

if(!require(foreach, quietly = TRUE)){
    install.packages("foreach",quiet = T)
    library(foreach, quietly = TRUE)
}

if(!require(doSNOW, quietly = TRUE)){
    install.packages("doSNOW",quiet = T)
    library(doSNOW, quietly = TRUE)
}

if(!require(Biostrings, quietly = TRUE)){
    BiocManager::install('Biostrings')
    library(Biostrings, quietly = TRUE)
}

if(!require(Rbowtie2, quietly = TRUE)){
    BiocManager::install('Rbowtie2')
    library(Rbowtie2, quietly = TRUE)
}

if(!require(chunked, quietly = TRUE)){
    install.packages("chunked",quiet = T)
    library(chunked)
}

if(!require(Rsamtools, quietly = TRUE)){
    install.packages("Rsamtools",quiet = T)
    library(Rsamtools)
}

options(scipen = 9999)
# check imputs 
if(!"RefGenome" %in% names(opt)){
    stop("Fastq file not provided. See script usage (--help)")
}

if(!"index" %in% names(opt)){
    stop("Bowtie2 indexed genome not provided. See script usage (--help)")
}

#loading reference fa
reference=readDNAStringSet(opt$RefGenome)

#chr info
chr=names(reference)
size=width(reference)
#weird r repsonse need +1 -1
cum_max=cumsum(size+1-1)

genome.Chromsizes=tibble(chr=chr, 
                         size=size,
                         cumulative_max=cum_max)%>%
    mutate(cumulative_min=lag(cumulative_max,1,default = 0))

cl=makeCluster(opt$cores)
# Identify reagins in wich the sequence is known
position = parLapply(reference, function(x) {
    library(tidyverse) 
    y=str_locate_all(x, '.[TACG]+')
    return(tibble(start=y[[1]][,1],end=y[[1]][,2]))
}, cl = cl)

# add colum with chr name
for(name in names(position)){
    position[names(position)==name][[1]]$chr=name
}

# convert position on chr in order to have non repetitive positions in absence of chr
position = bind_rows(position) %>%
    inner_join(genome.Chromsizes, by = 'chr') %>%
    mutate(start = start + cumulative_min,
           end = end + cumulative_min - opt$reads_size) %>%
    select(start, end)

#identify max position
Max = max(position$end)

reads_number=opt$reads_number

#look for seeds 
if(opt$paired_ends){
    #initialize simulated reads
    simulated_reads_1 = NULL
    simulated_reads_2 = NULL
    opt$insert_size=as.numeric(str_split(opt$insert_size,':')[[1]])
    while(reads_number != 0) {
        simulated_reads_temp_1 = sample(1:Max, reads_number)
        simulated_reads_temp_2 = simulated_reads_temp_1 + sample(opt$insert_size[1]:opt$insert_size[2], reads_number,replace = T)
        #check that has not been selected identical
        to_keep=!simulated_reads_temp_1 %in% simulated_reads_1 | !simulated_reads_temp_2 %in% simulated_reads_2
        simulated_reads_temp_1 = simulated_reads_temp_1[to_keep]
        simulated_reads_temp_2 = simulated_reads_temp_2[to_keep]
        
        simulated_reads_temp_1 = lapply(simulated_reads_temp_1, function(x, ranges, size) {
            #check values in inside ranges
            if (any(x - ranges$start >= 0 & x - ranges$end <= 0)) {
                return(x)
                
            } else{
                return(NA)
            }
            
        }, ranges = position, size = opt$reads_size)
        
        simulated_reads_temp_2 = lapply(simulated_reads_temp_2, function(x, ranges, size, already_selected) {
            #check values in inside ranges
            if (any(x - ranges$start >= 0 & x - ranges$end <= 0)) {
                return(x)
                
            } else{
                return(NA)
            }
            
        }, ranges = position, size = opt$reads_size)
        
        simulated_reads_temp_1 = unlist(simulated_reads_temp_1)
        simulated_reads_temp_2 = unlist(simulated_reads_temp_2)
        to_keep=!(is.na(simulated_reads_temp_1) | is.na(simulated_reads_temp_2))
        simulated_reads_temp_1 = simulated_reads_temp_1[to_keep]
        simulated_reads_temp_2 = simulated_reads_temp_2[to_keep]
        
        reads_number = reads_number - length(simulated_reads_temp_1)
        simulated_reads_1 = c(simulated_reads_1, simulated_reads_temp_1)
        simulated_reads_2 = c(simulated_reads_2, simulated_reads_temp_2)
    }
}else{
    #initialize simulated reads
    simulated_reads = NULL
while(reads_number != 0) {
    simulated_reads_temp = sample(1:Max, reads_number)
    #check that has not been selected identical
    simulated_reads_temp = simulated_reads_temp[!simulated_reads_temp %in% simulated_reads]
    simulated_reads_temp = lapply(simulated_reads_temp, function(x, ranges, size, already_selected) {
        #check values in inside ranges
        if (any(x - ranges$start >= 0 & x - ranges$end <= 0)) {
            return(x)
            
        } else{
            return(NA)
        }
        
    }, ranges = position, already_selected = simulated_reads, size = opt$reads_size)
    simulated_reads_temp = unlist(simulated_reads_temp)
    simulated_reads_temp = simulated_reads_temp[!is.na(simulated_reads_temp)]
    reads_number = reads_number - length(simulated_reads_temp)
    simulated_reads = c(simulated_reads, simulated_reads_temp)

}
}

#function to reallocate reads on chromosomes
reallocate_positions_on_chr=function(simulated_reads,genome.Chromsizes) {
simulated_reads = tibble(pos = simulated_reads) %>%
    mutate(chr = 'x',
           min = 'y')


for (i in 1:length(genome.Chromsizes$chr)) {
    simulated_reads = simulated_reads %>%
        mutate(
            chr = ifelse(
                pos > genome.Chromsizes$cumulative_min[i] &
                    pos <= genome.Chromsizes$cumulative_max[i] ,
                genome.Chromsizes$chr[i],
                chr
            ),
            min = ifelse(
                pos > genome.Chromsizes$cumulative_min[i] &
                    pos <= genome.Chromsizes$cumulative_max[i] ,
                pos - genome.Chromsizes$cumulative_min[i],
                min
            )
        )
}

simulated_reads = simulated_reads %>%
    mutate(min = as.numeric(min),
           start = min,
           end = min + opt$reads_size) %>%
    select(chr, start, end)

return(simulated_reads)
}

if(opt$paired_ends){
    simulated_reads_1=reallocate_positions_on_chr(simulated_reads = simulated_reads_1,genome.Chromsizes = genome.Chromsizes )
    simulated_reads_2=reallocate_positions_on_chr(simulated_reads = simulated_reads_2,genome.Chromsizes = genome.Chromsizes )
    simulated_reads_1$order=1:length(simulated_reads_1$chr)
    simulated_reads_2$order=1:length(simulated_reads_2$chr)
    
}else{
    simulated_reads=reallocate_positions_on_chr(simulated_reads = simulated_reads,genome.Chromsizes = genome.Chromsizes )
    simulated_reads$order=1:length(simulated_reads$chr)
    
}

#recover reads and mutate them
recover_and_mutate=function(simulated_reads,opt,reference){
    library(tidyverse,foreach,doSNOW)
# recover sequences 
cl=makeCluster(opt$cores)
registerDoSNOW(cl)
simulated_reads=foreach(Chr=genome.Chromsizes$chr,.combine = 'rbind',.packages = 'tidyverse')%dopar%{
    x=tibble(reads=str_sub(reference[names(reference)==Chr],start = simulated_reads$start[simulated_reads$chr==Chr],
            end=simulated_reads$end[simulated_reads$chr==Chr]-1),order=simulated_reads$order[simulated_reads$chr==Chr])
}
stopCluster(cl)

# simulate mutations 0.1 % rate
mutate=sample(1:length(simulated_reads$reads),0.001*length(simulated_reads$reads))
to_mutate=simulated_reads[mutate,]
simulated_reads=simulated_reads[-mutate,]

to_mutate = to_mutate %>%
    group_by(reads) %>%
    mutate(
        len = str_length(reads),
        mutate = sample(1:40, 1),
        before = str_sub(reads, 0, mutate - 1),
        after = str_sub(reads, mutate + 1, len),
        mutate_b = str_sub(reads, mutate, mutate),
        mutated_base = ifelse(
            mutate_b == 'A',
            sample(c('T', 'C', 'G'), 1),
            ifelse(
                mutate_b == 'C',
                sample(c('A', 'T', 'G'), 1),
                ifelse(
                    mutate_b == 'G',
                    sample(c('A', 'T', 'C'), 1),
                    ifelse(mutate_b == 'T' ,
                           sample(c('A', 'C', 'G'), 1),
                           sample(c('A', 'C', 'G', 'T'), 1))
                )
            )
        ),
        new_seq = paste0(before, mutated_base, after),
        check = str_length(new_seq)
    ) %>%
    ungroup() %>%
    select(new_seq, order) %>%
    `colnames<-`(names(simulated_reads))
simulated_reads=rbind(simulated_reads,to_mutate)
return(simulated_reads)
}

if(opt$paired_ends){
    #recover strings reads and mutate some of them 
    simulated_reads_1=recover_and_mutate(simulated_reads = simulated_reads_1,opt = opt,reference=reference)
    simulated_reads_2=recover_and_mutate(simulated_reads = simulated_reads_2,opt = opt,reference=reference)
    #calculate reverse complement for PE
    cl=makeCluster(opt$cores)
    simulated_reads_2$reads=parSapply(simulated_reads_2$reads, function(x) as.character( Biostrings::reverseComplement(Biostrings::DNAString(x))),cl=cl)
}else{
    #recover strings reads and mutate some of them 
    simulated_reads=recover_and_mutate(simulated_reads = simulated_reads,opt = opt,reference=reference)
}

#reshape reads for the fastq file
reshape_and_save =function(simulated_reads,file){
    simulated_reads %>%
        arrange(order)%>%
        select(-order)%>%
        `colnames<-`('2') %>%
        mutate(
            index = 1:n(),
            `1` = paste0('@read', 1:n()),
            `3` = '+',
            `4` = paste0(rep('D', opt$reads_size), collapse = '')
        ) %>%
        gather('pos', 'towrite', -index) %>%
        arrange(index, pos) %>%
        select(towrite)%>%
        write_delim(path = file,col_names = F) 
} 

if(opt$paired_ends){
    #save fastq files
    reshape_and_save(simulated_reads_1,file = paste0(opt$output_dir, basename(opt$index), '_simulated_reads_1.fq'))
    reshape_and_save(simulated_reads_2,file = paste0(opt$output_dir, basename(opt$index), '_simulated_reads_2.fq'))
    
    #remove from memory
    rm('simulated_reads_1')
    rm('simulated_reads_2')
    
    #align with bowtie2
    bowtie2(
        bt2Index = opt$index,
        samOutput = paste0(opt$output_dir, basename(opt$index), '_simulated_reads.sam'),
        seq1 = paste0(opt$output_dir, basename(opt$index), '_simulated_reads_1.fq'),
        seq2 = paste0(opt$output_dir, basename(opt$index), '_simulated_reads_2.fq'),
        ... = paste0('-k 1 --phred33 --ignore-quals -p ', opt$cores),
        overwrite=TRUE
    )
    
    # remove from hd simulated_reads.fq
    system(paste0('rm ',opt$output_dir, basename(opt$index),'_simulated_reads_1.fq'))
    system(paste0('rm ',opt$output_dir, basename(opt$index),'_simulated_reads_2.fq'))
    
}else{
    #save fastq file
    reshape_and_save(simulated_reads,file = paste0(opt$output_dir, basename(opt$index), '_simulated_reads.fq'))    #remove from memory
    rm('simulated_reads')
    
    #align with bowtie2
    bowtie2(
        bt2Index = opt$index,
        samOutput = paste0(opt$output_dir, basename(opt$index), '_simulated_reads.sam'),
        seq1 = paste0(opt$output_dir, basename(opt$index), '_simulated_reads.fq'),
        ... = paste0('-k 1 --phred33 --ignore-quals -p ', opt$cores),
        overwrite=TRUE
    )
    
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
)
}else{
    param <- ScanBamParam(what=c('rname','pos','mapq'),
                          flag=scanBamFlag(isUnmappedQuery = F))
    
    bins =as.data.frame(scanBam(dir.bam,param=param))%>%
        filter(mapq >= 30)%>%
        mutate(read=1)%>%
        select('rname', 'pos', 'read')%>%
        `colnames<-`(c('chr', 'pos','read'))
}

if(!opt$variable_bins){
    
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
        mutate(mappability=reads/max(reads),
               mappability_th=ifelse(
                   reads >= 0.8*(opt$reads_number/(sum(position$end-position$start)/opt$bin_size)),T,F
               ))%>%
        select(chr,start,end,mappability,mappability_th)
    
}else{

bins=bins%>%
    arrange(pos)%>%
    group_by(chr)%>%
    mutate(group=ceiling(cumsum(read)/opt$reads_number))%>%
    ungroup()%>%
    group_by(chr,group)%>%
    summarise(start=min(pos),
              end=max(pos)+opt$reads_size)%>%
    ungroup()%>%
    mutate(start_plus_1=lead(start),
           end=ifelse(is.na(start_plus_1),end,
                      ifelse(start_plus_1==end,end,
                             start_plus_1)))%>%
    select(chr,start,end)%>%
    inner_join(genome.Chromsizes)%>%
    group_by(chr)%>%
    mutate(start=ifelse(start==min(start),0,start),
           end=ifelse(start==max(start),size,end)
           )%>%
    ungroup()%>%
    mutate(mappability=1,
           mappability_th=T)%>%
    select(chr,start,end,mappability,mappability_th)
}
#delete file
system(paste0('rm ',opt$output_dir, basename(opt$index), '_simulated_reads.bam*'))

#calculate gc % peer bin
cl=makeCluster(opt$cores)
registerDoSNOW(cl)

bins=foreach(i=1:length(bins$chr),.combine = 'rbind', .packages =c('Biostrings','tidyverse') )%dopar%{
    seq=paste(subseq(x = reference[names(reference)==bins$chr[i]],start=bins$start[i]+1, end=bins$end[i]))
    cbind(bins[i,],tibble(gc_frequency=str_count(seq,'G|C')/str_length(seq)))
}

stopCluster(cl)

#free memory
rm('reference')
genome_size=sum(genome.Chromsizes$size)

cl=makeCluster(opt$cores)
registerDoSNOW(cl)

#calculate correction for dimap
simulated_data=foreach(i=c(1,seq(10,310,20)),.combine = 'rbind',.packages = c('tidyverse','foreach'), .verbose = T)%dopar%{
    reads=foreach(h=1:1000,.combine = 'rbind',.packages = 'tidyverse')%do%{
       bins%>%
            filter(mappability_th)%>%
            mutate(reads=rpois(n(),i))%>%
            group_by(chr)%>%
            mutate(reads_1 = lag(reads, 1))%>%
            drop_na() %>%
            ungroup()%>%
            mutate(mapd = (reads - reads_1)/mean(reads))%>%
            summarise(mapd = median(abs(mapd) - median(mapd)),
                      inv_cov=1/(sum(reads))
                      )%>%
            mutate(Cell=h,
                   mean_reads=i)
            
    }
    reads
}

stopCluster(cl)

LM_mapd_coverage = lm(formula = mapd ~ inv_cov, data = simulated_data)

# add values for normalization
bins = bins %>%
    mutate(m = LM_mapd_coverage$coefficients[['inv_cov']],
           k = LM_mapd_coverage$coefficients[['(Intercept)']])

#write bisns with info
write_tsv(bins, paste0(opt$output_dir, basename(opt$index), '_bins_',ifelse(opt$paired_ends,'PE','SE'),ifelse(opt$variable_bins,'_variable_bins','_fixed_bins'),'.tsv'))


