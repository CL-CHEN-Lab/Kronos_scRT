#!/usr/local/bin/Rscript --slave
#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)

option_list = list(
    make_option(
        c("-D", "--directory"),
        type = "character",
        default = NULL,
        help = "Single cell Bamfiles directory",
        metavar = "character"
    ),
    make_option(
        c("-B", "--bins"),
        type = "character",
        default = NULL,
        help = "File with bins produced by Kronos binning",
        metavar = "character"
    ),
    make_option(
        c("-X", "--keep_X"),
        type = "logical",
        action = 'store_TRUE',
        default = F,
        help = "Keep X chromosomes. [default= %default]",
        metavar = "logical"
    ),
    make_option(
        c("-Y", "--keep_Y"),
        type = "logical",
        default = F,
        action = 'store_TRUE',
        help = "Keep Y chromosome. [default= %default]",
        metavar = "logical"
    ),
    make_option(
        c("-m", "--min_n_reads"),
        type = "double",
        default = 200000,
        action = 'store',
        help = "Min n of reads to keep a cell in the analysis [default= %default]",
        metavar = "double"
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
        help = "Output folder. [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("-e", "--ExpName"),
        type = "character",
        default = 'Exp',
        action = 'store',
        help = "Experiment name. [default= %default]",
        metavar = "character"
    )
)

opt = parse_args(OptionParser(option_list = option_list),
                 convert_hyphens_to_underscores = T)

#load needed packages
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))
suppressPackageStartupMessages(library(doSNOW, quietly = TRUE))
suppressPackageStartupMessages(library(Rsamtools, quietly = TRUE))
suppressPackageStartupMessages(library(DescTools, quietly = TRUE))
suppressPackageStartupMessages(library(DNAcopy, quietly = TRUE))
suppressPackageStartupMessages(library(gplots, quietly = TRUE))
suppressPackageStartupMessages(library(MASS, quietly = TRUE))

if (str_extract(opt$output_dir, '.$') != '/') {
    opt$output_dir = paste0(opt$output_dir, '/')
}

system(paste0('mkdir -p ', opt$output_dir))

#check inputs
if (!'bins' %in% names(opt)) {
    stop("Bins with gc percentage not provided. See script usage (--help)")
}

if (!"directory" %in% names(opt)) {
    stop("Directory to bam files not provided. See script usage (--help)")
}

# load bins and gc percentage
bins = read_tsv(opt$bins, col_types = cols())

if (str_extract(opt$directory, '.$') != '/') {
    opt$directory = paste0(opt$directory, '/')
}


#find bam files
files = list.files(paste0(opt$directory))
files = files[files %like% '%.bam']

#chr info
genome.Chromsizes <-  bins %>%
    group_by(chr) %>%
    summarise(size = max(end))


if (opt$keep_X & opt$keep_Y) {
    chromosome = genome.Chromsizes$chr[genome.Chromsizes$chr %like% 'chr[0-9XY]{1,2}']
    genome_size = sum(genome.Chromsizes$size[genome.Chromsizes$chr  %like% 'chr[0-9]{1,2}'])
    genome_size = genome_size + genome.Chromsizes$size[genome.Chromsizes$chr  %like% 'chrX']
    genome_size = genome_size + genome.Chromsizes$size[genome.Chromsizes$chr  %like% 'chrY']
    bins = bins[bins$chr %in% chromosome, ]
} else if (!opt$keep_X & opt$keep_Y) {
    chromosome = genome.Chromsizes$chr[genome.Chromsizes$chr %like% 'chr[0-9Y]{1,2}']
    genome_size = sum(genome.Chromsizes$size[genome.Chromsizes$chr  %like% 'chr[0-9]{1,2}'])
    genome_size = genome_size + genome.Chromsizes$size[genome.Chromsizes$chr  %like% 'chrY']
    bins = bins[bins$chr %in% chromosome, ]
} else if (opt$keep_X & !opt$keep_Y) {
    chromosome = genome.Chromsizes$chr[genome.Chromsizes$chr %like% 'chr[0-9X]{1,2}']
    genome_size = sum(genome.Chromsizes$size[genome.Chromsizes$chr  %like% 'chr[0-9]{1,2}'])
    genome_size = genome_size + genome.Chromsizes$size[genome.Chromsizes$chr  %like% 'chrX']
    bins = bins[bins$chr %in% chromosome, ]
} else if (!opt$keep_X & !opt$keep_Y) {
    chromosome = genome.Chromsizes$chr[genome.Chromsizes$chr %like% 'chr[0-9]{1,2}']
    genome_size = sum(genome.Chromsizes$size[genome.Chromsizes$chr  %like% 'chr[0-9]{1,2}'])
    bins = bins[bins$chr %in% chromosome, ]
}

# SE or PE ?
type = bins %>%
    dplyr::select(type) %>%
    unique() %>%
    pull()

cl = makeCluster(opt$cores)
registerDoSNOW(cl)

bins = bins %>%
    group_by(chr) %>%
    mutate(bin = 1:n())

# calculating profile
files = foreach (file = files,
                 .packages = c('Rsamtools', 'tidyverse', 'foreach')) %dopar% {
                     if (type == 'PE') {
                         bins_median_size = median(bins$end - bins$start)
                         
                         # Singe reads
                         param1 <- ScanBamParam(
                             what = c('rname', 'pos', 'mapq'),
                             flag = scanBamFlag(
                                 isPaired = F,
                                 isUnmappedQuery = F,
                                 isDuplicate = F
                             ),
                             mapqFilter = 30
                         )
                         
                         #First in a pair proper paired reads
                         param2 <-
                             ScanBamParam(
                                 what = c('qname', 'rname', 'pos', 'mapq'),
                                 flag = scanBamFlag(
                                     isPaired = T,
                                     isUnmappedQuery = F,
                                     isFirstMateRead = T,
                                     isDuplicate = F
                                 ),
                                 mapqFilter = 30
                             )
                         
                         #Second in a pair proper paired reads
                         param3 <-
                             ScanBamParam(
                                 what = c('qname', 'rname', 'pos', 'mapq'),
                                 flag = scanBamFlag(
                                     isPaired = T,
                                     isUnmappedQuery = F,
                                     isSecondMateRead = T,
                                     isDuplicate = F
                                 ),
                                 mapqFilter = 30
                             )
                         # load first in a pair
                         FP = as.data.frame(scanBam(paste0(opt$directory, file), param = param2)) %>%
                             drop_na() %>%
                             group_by(rname) %>%
                             mutate(bin = ceiling(pos / bins_median_size))
                         
                         # load second in a pair
                         SP = as.data.frame(scanBam(paste0(opt$directory, file), param = param3)) %>%
                             drop_na() %>%
                             group_by(rname) %>%
                             mutate(bin = ceiling(pos / bins_median_size))
                         
                         # select identifiers reads
                         qname_fp = FP %>% ungroup() %>% dplyr::select(qname, rname, bin) %>% `colnames<-`(c('qname', 'rname', 'mate_bin'))
                         qname_sp = SP %>% ungroup() %>% dplyr::select(qname, rname, bin) %>% `colnames<-`(c('qname', 'rname', 'mate_bin'))
                         
                         # if a read in a pair has it's paired in the same bin it counts as 1/2 a read eles one read.
                         FP = FP %>%
                             left_join(qname_sp, by = c("qname", 'rname')) %>%
                             mutate(read = ifelse(bin != mate_bin |
                                                      is.na(mate_bin), 1, 0.5)) %>%
                             dplyr::select('rname', 'pos', 'read', 'bin') %>%
                             `colnames<-`(c('chr', 'pos', 'read', 'bin')) %>%
                             ungroup()
                         
                         SP = SP %>%
                             left_join(qname_fp, by = c("qname", 'rname')) %>%
                             mutate(read = ifelse(bin != mate_bin |
                                                      is.na(mate_bin), 1, 0.5)) %>%
                             dplyr::select('rname', 'pos', 'read', 'bin') %>%
                             `colnames<-`(c('chr', 'pos', 'read', 'bin')) %>%
                             ungroup()
                         
                         #load non paired reads
                         SR = as.data.frame(scanBam(paste0(opt$directory, file), param = param1)) %>%
                             drop_na() %>%
                             mutate(read = 1) %>%
                             dplyr::select('rname', 'pos', 'read') %>%
                             `colnames<-`(c('chr', 'pos', 'read'))
                         
                         # include non paired reads if existing
                         # sum all the reads in a bin 
                         if (length(SR$chr) != 0) {
                             sam = rbind(SR = SR  %>%
                                             group_by(chr) %>%
                                             mutate(bin = ceiling(pos /
                                                                      bins_median_size)),
                                         FP ,
                                         SP) %>%
                                 group_by(chr, bin) %>%
                                 summarise(reads = sum(read))
                             
                         } else{
                             sam = rbind(FP ,
                                         SP) %>%
                                 group_by(chr, bin) %>%
                                 summarise(reads = sum(read))
                         }
                     } else if (type == 'SE') {
                         param <- ScanBamParam(
                             what = c('rname', 'pos', 'mapq'),
                             flag = scanBamFlag(isUnmappedQuery = F,
                                                isDuplicate = F),
                             mapqFilter = 30
                         )
                         
                         # calculate reads in each bin
                         sam = as.data.frame(scanBam(paste0(opt$directory, file), param =
                                                         param)) %>%
                             mutate(read = 1) %>%
                             dplyr::select('rname', 'pos', 'read') %>%
                             `colnames<-`(c('chr', 'pos', 'read')) %>%
                             group_by(chr) %>%
                             mutate(bin = ceiling(pos / bins_median_size)) %>%
                             group_by(chr, bin) %>%
                             summarise(reads = sum(read)) %>%
                             ungroup()
                     } else{
                         stop('Bins file does not contain a correct sequencing type (SE/PE)')
                     }
                     if (sum(sam$reads) >= opt$min_n_reads) {
                         # save one file per cell with a coverage track
                         bins %>%
                             left_join(sam %>% mutate(chr = as.character(chr)), by =
                                           c('chr', 'bin')) %>%
                             mutate(reads = ifelse(is.na(reads), 0, reads),
                                    Cell = file) %>%
                             dplyr::select(-bin) %>%
                             write_tsv(paste0(opt$output_dir, file, '.tmp'))
                         paste0(file, '.tmp')
                     }
                     
                 }
stopCluster(cl)

#merge tenp files, rm single ones and load data
system(paste0('cat ', opt$output_dir, files[[1]], ' > ', opt$output_dir, 'data.tmp'))
system(paste0(
    'for i in ',
    paste0(opt$output_dir, files[-1], collapse = ' '),
    "; do sed '1d' $i >> ",
    opt$output_dir,
    'data.tmp; done'
))
system(paste0('rm ', paste0(opt$output_dir, files, collapse = ' ')))

#load data
data = read_tsv(paste0(opt$output_dir, 'data.tmp'), col_types = cols())

#deletefile
system(paste0('rm ', opt$output_dir, 'data.tmp'))

# correct for mappability and normalize for GC content (norm reads= reads * median reads  / median reads per interval of GC (rounded at 2 digits))
gc_correction_value = data %>%
    filter(mappability_th) %>%
    group_by(chr, start, end, mappability) %>%
    summarise(reads = sum(reads),
              gc_frequency = unique(gc_frequency)) %>%
    mutate(
        gc_frequency = round(gc_frequency, 1),
        reads_mappability = reads / mappability
    ) %>%
    group_by(gc_frequency) %>%
    mutate(mGC = median(reads_mappability)) %>%
    ungroup() %>%
    mutate(
        m = median(reads_mappability),
        gc_corretion_values = ifelse(mGC == 0, NA, m / mGC)
    ) %>%
    dplyr::select(chr, start, end, gc_corretion_values)

data = left_join(data, gc_correction_value, by = c('chr', 'start', 'end')) %>%
    mutate(
        reads_mappability = reads / mappability,
        gc_corrected_reads = reads_mappability * gc_corretion_values
    )

data_500Kb = data %>%
    group_by(Cell) %>%
    mutate(cum = ceiling(cumsum(end - start) / 500000)) %>%
    group_by(Cell, chr, cum) %>%
    summarise(
        start = min(start),
        end = max(end),
        reads = sum(reads, na.rm = T),
        gc_corrected_reads = sum(gc_corrected_reads, na.rm = T)
    )

coverage = data_500Kb %>%
    group_by(Cell) %>%
    summarise(coverage = 1000000 * sum(reads) / genome_size)


#calculate normalized_MAPD and mapd
mapd = data_500Kb %>%
    group_by(Cell, chr) %>%
    mutate(read_1n = lag(gc_corrected_reads, 1)) %>%
    drop_na() %>%
    ungroup() %>%
    group_by(Cell) %>%
    mutate(
        mean_n = mean(gc_corrected_reads),
        normalized_mapd = (gc_corrected_reads - read_1n) / mean_n
    ) %>%
    summarise(normalized_mapd = median(abs(normalized_mapd) - median(normalized_mapd)))


#calculate DiMApd
mapd = inner_join(mapd, coverage, by = "Cell") %>%
    mutate(normalized_dimapd = normalized_mapd * sqrt(coverage),
           d = coverage - median(coverage))

LM_mapd_coverage = lm(formula = normalized_dimapd ~ d , data = mapd)
mapd = mapd %>%
    mutate(
        normalized_dimapd = 1 + normalized_dimapd - (
            d * LM_mapd_coverage$coefficients[[2]] + LM_mapd_coverage$coefficients[[1]]
        ),
        is_high_dimapd = F
    ) %>%
    dplyr::select(-d)

#fit dimapd to gaussian dist
mem = 0
while (T) {
    fit <-
        fitdistr(mapd$normalized_dimapd[!mapd$is_high_dimapd], 'normal')
    
    mapd = mapd %>%
        mutate(is_high_dimapd = ifelse(
            pnorm(
                normalized_dimapd,
                mean = fit$estimate[[1]],
                sd = fit$estimate[[2]],
                lower.tail = F
            ) <= 0.01,
            T,
            F
        ))
    if (mem == sum(mapd$is_high_dimapd)) {
        break
    } else{
        mem = sum(mapd$is_high_dimapd)
    }
    
}

#spread data for segmentation
data = data %>%
    filter(mappability_th) %>%
    group_by(Cell) %>%
    dplyr::select(chr, start, end, gc_corrected_reads, Cell) %>%
    ungroup() %>%
    spread(Cell, gc_corrected_reads) %>%
    drop_na()

cl = makeCluster(opt$cores)
registerDoSNOW(cl)

segment.smoothed.CNA.object = foreach(
    file = names(data)[!names(data) %in% c('chr', 'start', 'end')],
    .combine = 'rbind',
    .packages = c('DNAcopy', 'tidyverse')
) %dopar% {
    # create object
    mt <- as.matrix(data[, names(data) %in% file])
    CNA.object <-
        CNA(mt,
            data$chr,
            data$start,
            sampleid = file)
    
    #free memory
    rm('mt')
    
    # smooth data
    smoothed.CNA.object <- smooth.CNA(CNA.object)
    
    # free memory
    rm('CNA.object')
    
    # segment
    segment.smoothed.CNA.object <- segment(smoothed.CNA.object)
    
    #free memory
    rm('smoothed.CNA.object')
    
    as_tibble(segment.smoothed.CNA.object$output)
}

#free memory
rm('data')


IDs = unique(segment.smoothed.CNA.object$ID)

MMS = segment.smoothed.CNA.object %>%
    summarise(
        min = quantile(seg.mean, 0.05)[[1]],
        max = quantile(seg.mean, 0.95)[[1]],
        step = (max - min) / 1000
    )

bin_size = bins[1, ] %>%
    mutate(bs = end - start) %>%
    pull(bs)

# identify CN based on minimum of the target function
CNV_correction = foreach(
    id = IDs,
    .combine = 'rbind',
    .packages = c('foreach', 'tidyverse')
) %dopar% {
    s = segment.smoothed.CNA.object %>%
        filter(ID == id)
    
    weitghts = (s$loc.end - s$loc.start) / bin_size
    
    possible_factors = foreach(i = seq(MMS$min, MMS$max, MMS$step),
                               .combine = 'rbind') %do% {
                                   TargetF = sqrt(sum((weitghts * sinpi(
                                       s$seg.mean / i
                                   ) ^ 2)))
                                   mean_cn = weighted.mean(round(s$seg.mean / i), weitghts)
                                   Variability = 100 * sd(rep(s$seg.mean, weitghts)) / weighted.mean(s$seg.mean, weitghts)
                                   TargetF = tibble(
                                       possible_factors = TargetF,
                                       X = i,
                                       mean_cn = mean_cn,
                                       Variability = Variability
                                   )
                                   
                               }
    
    Var = unique(possible_factors$Variability)
    min = possible_factors$possible_factors[which(diff(sign(
        diff(possible_factors$possible_factors)
    )) == 2) + 1]
    possible_factors = possible_factors %>%
        filter(possible_factors %in% min,
               mean_cn < 8)
    
    if (Var < 5) {
        selected = possible_factors$X[possible_factors$mean_cn[which(abs(possible_factors$mean_cn -
                                                                             2) == min(abs(possible_factors$mean_cn - 2)))]]
        mean_cn = possible_factors$mean_cn[possible_factors$mean_cn[which(abs(possible_factors$mean_cn -
                                                                                  2) == min(abs(possible_factors$mean_cn - 2)))]]
        PloConf = -2
    } else{
        selected = min(possible_factors$possible_factors)
        PloConf = nth(possible_factors$possible_factors[base::order(possible_factors$possible_factors, decreasing = T)], n = -2) -
            selected
        mean_cn = possible_factors$mean_cn[possible_factors$possible_factors ==
                                               selected]
        selected = possible_factors$X[possible_factors$possible_factors ==
                                          selected]
        
    }
    
    tibble(
        ID = id,
        X = selected,
        ploidy_confidence = PloConf,
        mean_ploidy = mean_cn
    )
}

stopCluster(cl)

# called CNV
CNV = segment.smoothed.CNA.object %>%
    inner_join(CNV_correction, by = "ID") %>%
    mutate(CNV = round(seg.mean / X, 0)) %>%
    dplyr::select(-X, -seg.mean, -num.mark, -ploidy_confidence, -mean_ploidy) %>%
    `colnames<-`(c('Cell', 'chr', 'start', 'end', 'copy_number'))


CNV_correction = CNV_correction %>%
    dplyr::select(-X)

# calculate mean CNV
mapd = mapd %>%
    mutate(Cell = str_replace_all(Cell, '-', '.')) %>%
    inner_join(CNV_correction, by = c('Cell' = 'ID')) %>%
    mutate(
        is_noisy = ifelse(is_high_dimapd | ploidy_confidence < 2, T, F),
        coverage_per_1Mbp = coverage
    ) %>%
    dplyr::select(
        Cell,
        normalized_dimapd,
        mean_ploidy,
        ploidy_confidence,
        is_high_dimapd,
        is_noisy,
        coverage_per_1Mbp
    )
#write file
mapd %>%
    write_csv(paste0(opt$output_dir, opt$ExpName, '_per_Cell_summary_metrics.csv'),
              col_names = T)

CNV %>%
    write_tsv(paste0(opt$output_dir, opt$ExpName, '_cnv_calls.bed'),
              col_names = T)
