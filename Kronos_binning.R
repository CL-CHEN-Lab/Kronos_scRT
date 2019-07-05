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
        help = "Chr Size file",
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
        help = "number of the simulated reads. [default= %default bp]",
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

#initialize simulated reads
simulated_reads = NULL

#look for seeds 
while(opt$reads_number != 0) {
    simulated_reads_temp = sample(1:Max, opt$reads_number)
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
    opt$reads_number = opt$reads_number - length(simulated_reads_temp)
    simulated_reads = c(simulated_reads, simulated_reads_temp)
}

#reallocate positions on chr 
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

# recover sequences 
cl=makeCluster(opt$cores)
registerDoSNOW(cl)
simulated_reads=foreach(Chr=genome.Chromsizes$chr,.combine = 'rbind',.packages = 'tidyverse')%dopar%{
    x=tibble(reads=str_sub(reference[names(reference)==Chr],start = simulated_reads$start[simulated_reads$chr==Chr],
            end=simulated_reads$end[simulated_reads$chr==Chr]-1))
}
stopCluster(cl)

# simulate mutations 0.1 % rate
mutate=sample(1:length(simulated_reads$reads),0.001*length(simulated_reads$reads))
to_mutate=simulated_reads[mutate,]
simulated_reads=simulated_reads[-mutate,]

change_base=function(N){
    if (N=='A') {
        return(sample(c('T','C','G'),1))
    }else if (N=='C') {
        return(sample(c('A','T','G'),1))
    }else if(N=='G') {
        return(sample(c('A','T','C'),1))
    }else if(N=='T') {
        return(sample(c('A','C','G'),1))
    }else{
        return(sample(c('A','C','G','T'),1))
    }
}

to_mutate=to_mutate%>%
    group_by(reads)%>%
    mutate(len=str_length(reads),
           mutate=sample(1:40,1),
           before=str_sub(reads,0,mutate-1),
           after=str_sub(reads,mutate+1,len),
           mutate_b=str_sub(reads,mutate,mutate),
           mutated_base=change_base(mutate_b),
           new_seq=paste0(before,mutated_base,after),
           check=str_length(new_seq))%>%
    ungroup()%>%
    select(new_seq)%>%
    `colnames<-`(names(simulated_reads))
simulated_reads=rbind(simulated_reads,to_mutate)
rm('to_mutate')

#reshape reads for the fastq file
simulated_reads = simulated_reads %>%
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
    write_delim(path = paste0(opt$output_dir,'simulated_reads.fq'),col_names = F)

#remove from memory
rm('simulated_reads')

#align with bowtie2
opt$index='~/alliners_compilers/mm10_bowtie2/mm10'
bowtie2(
    bt2Index = opt$index,
    samOutput = paste0(opt$output_dir, 'simulated_reads.sam'),
    seq1 = paste0(opt$output_dir, 'simulated_reads.fq'),
    ... = paste0('-k 1 --phred33  -p ', opt$cores),
    overwrite=TRUE
)
# remove from hd simulated_reads.fq
system(paste0('rm ',opt$output_dir,'simulated_reads.fq'))

# calculate bins

bins = read_tsv(
    paste0(opt$output_dir, 'simulated_reads.sam'),
    col_names = F,
    comment = '@'
) %>%
    select(X3, X4) %>%
    `colnames<-`(c('chr', 'pos')) %>%
    mutate(read = 1)%>%
    arrange(pos)%>%
    group_by(chr)%>%
    mutate(group=ceiling(cumsum(read)/otp$reads_number))%>%
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
    select(chr,start,end)

#calculate gc % peer bin
cl=makeCluster(opt$cores)
registerDoSNOW(cl)

bins=foreach(i=1:length(bins$chr),.combine = 'rbind', .packages =c('Biostrings','tidyverse') )%dopar%{
    seq=paste(subseq(x = reference[names(reference)==bins$chr[i]],start=bins$start[i]+1, end=bins$end[i]))
    cbind(bins[i,],tibble(counts_percent=str_count(seq,'G|C')/str_length(seq)))
}

stopCluster(cl)

#free memory
rm('reference')

write_tsv(bins, paste0(opt$output_dir, 'bins_with_gc_percentage.tsv'))
system(paste0('rm ',opt$output_dir, 'simulated_reads.sam'))


