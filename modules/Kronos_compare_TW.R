#!/usr/local/bin/Rscript --slave
#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)
options(warn = 1,scipen = 999)

option_list = list(
    make_option(
        c("-F", "--file"),
        type = "character",
        default = NULL,
        help = "Variability file produced by Kronos RT, if multiple files are provided they have to be separated by a comma",
        metavar = "character"
    ),
    make_option(
        c("-R", "--Annotation"),
        type = "character",
        default = NULL,
        help = "Genome annotation. chr<TAB>start<TAB>end<TAB>annotation. No header.",
        metavar = "character"
    ),
    make_option(
        c("-b", "--both_annotations"),
        type = "logical",
        default = F,
        action = "store_true",
        help = "Plot Twidth divided by both annotations",
        metavar = "logical"
    ),
    make_option(
        c("-r", "--Annotation2"),
        type = "character",
        default = NULL,
        help = "Second genome annotation. chr<TAB>start<TAB>end<TAB>annotation. No header.",
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
    ),
    make_option(
        c("-m", "--min_overlap"),
        type = "numeric",
        default = 100,
        help = "min overlap to apply the annotation [default= %default]",
        metavar = "numeric"
    )
)

#recover inputs
opt = parse_args(object = OptionParser(option_list = option_list))

#load packages
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))
suppressPackageStartupMessages(library(GenomicRanges, quietly = TRUE))
suppressPackageStartupMessages(library(ggpubr, quietly = TRUE))

#set plotting theme
theme_set(theme_bw())

#check inputs
if (!'file' %in% names(opt)) {
    stop("Variability file must be provided. See script usage (--help)")
}
if (!'Annotation' %in% names(opt)) {
    stop("Annotation file must be provided. See script usage (--help)")
}
if (str_extract(opt$out, '.$') != '/') {
    opt$out = paste0(opt$out, '/')
}

system(paste0('mkdir -p ', opt$out))

opt$file = str_split(opt$file, ',')[[1]]

#load files
data <-
    foreach(file = opt$file,
            .combine = 'rbind',
            .packages = 'tidyverse') %do% {
                read_tsv(file, col_types = cols())
            }

Annotation1_file <-
    foreach(file = opt$Annotation,
            .combine = 'rbind',
            .packages = 'tidyverse') %do% {
                read_tsv(file, col_types = cols(),col_names = c( "chr", "start", "end",'Annotation'))
            }

Annotation_file = read_tsv(opt$Annotation,
                               col_names = c('chr', 'start', 'end', 'annotation'), col_types = cols())%>%
        makeGRangesFromDataFrame(
            keep.extra.columns = T,
            seqnames.field = 'chr',
            end.field = 'end',
            start.field = 'start'
        )
    
    
    Bin = data %>%
        dplyr::select(chr, start, end) %>%
        unique() %>%
        makeGRangesFromDataFrame(
            keep.extra.columns = T,
            seqnames.field = 'chr',
            end.field = 'end',
            start.field = 'start'
        )
    
    #find overlaps
    hits = findOverlaps(Bin, Annotation_file,minoverlap = opt$min_overlap)
    
    #indo about overlapping regins
    overlaps <-
        pintersect(Annotation_file[subjectHits(hits)], Bin[queryHits(hits)])
    
    #Add annotation 1
    add = as_tibble(Bin[queryHits(hits)])
    not_add = as_tibble(Bin[-queryHits(hits)])
    
    #Based on the overlap define the predominant notation of each bin.
    #A category to be chosen has to be predominant (at least 60% of the tatal overlaps in the bin)
    Annotation = bind_rows(
        add %>%
            mutate(size = width(overlaps),
                   Cat1 = overlaps$annotation) %>%
            group_by(seqnames, start, end, Cat1) %>%
            summarise(n = sum(size)) %>%
            ungroup() %>%
            group_by(seqnames, start, end) %>%
            mutate(
                select = (n / sum(n) >= 0.6),
                Cat1 = ifelse(any(select), Cat1, '_Unknown_')
            ) %>%
            dplyr::select(seqnames, start, end, Cat1),
        not_add %>% mutate(Cat1 = '_Unknown_') %>%
            dplyr::select(seqnames, start, end, Cat1)
    )
    Annotation = Annotation %>%
        ungroup() %>%
        mutate(seqnames = as.character(seqnames))
    
    data = data %>%
        inner_join(Annotation, by = c("chr" = "seqnames", "start", "end"))
    
    if (opt$both_annotations &'Annotation2' %in% names(opt)) {
            
            Annotation_file = read_tsv(opt$Annotation2,
                                       col_names = c('chr', 'start', 'end', 'annotation'), col_types = cols()) %>%
                makeGRangesFromDataFrame(
                    keep.extra.columns = T,
                    seqnames.field = 'chr',
                    end.field = 'end',
                    start.field = 'start'
                )
            
            hits = findOverlaps(Bin, Annotation_file,minoverlap = opt$min_overlap)
            
            overlaps <-
                pintersect(Annotation_file[subjectHits(hits)], Bin[queryHits(hits)])
            
            add = as_tibble(Bin[queryHits(hits)])
            not_add = as_tibble(Bin[-queryHits(hits)])
            
            Annotation = bind_rows(
                add %>%
                    mutate(
                        size = width(overlaps),
                        Cat2 = overlaps$annotation
                    ) %>%
                    group_by(seqnames, start, end, Cat2) %>%
                    summarise(n = sum(size)) %>%
                    ungroup() %>%
                    group_by(seqnames, start, end) %>%
                    mutate(
                        select = (n / sum(n) >= 0.6),
                        Cat2 = ifelse(any(select), Cat2, '_Unknown_')
                    ) %>%
                    dplyr::select(seqnames, start, end, Cat2),
                not_add %>% mutate(Cat2 = '_Unknown_') %>%
                    dplyr::select(seqnames, start, end, Cat2)
            )
            Annotation = Annotation %>%
                ungroup() %>%
                mutate(seqnames = as.character(seqnames))
            
            data = data %>%
                inner_join(Annotation, by = c("chr" = "seqnames", "start", "end"))
            
        
        data = data %>%
            rbind(data %>%
                      mutate(Cat1 = '_ALL_'),
                  data %>%
                      mutate(Cat2 = '_ALL_'),
                  data %>%
                      mutate(Cat1 = '_ALL_',
                             Cat2 = '_ALL_'))
        

        #calculate tresholds 25% 75% replication keeping in account early and late domains
        
        T25_75 = function(df, name, Cat1, Cat2) {
            if (length(df$group) != 0) {
                model = tryCatch( nls(percentage ~ SSlogis(time, Asym, xmid, scal),
                                      data = df[, c('percentage', 'time')]%>%
                                          add_row(percentage=1,time=-10)%>%
                                          add_row(percentage=0,time=10),
                                      control = nls.control(maxiter = 100),
                                      algorithm = 'port',
                                      lower = c(Asym=1,xmid=0,scal=-1.5)
                ),
                  #If the data cannot be fitted with a Gauss-Cat2ton algorithm, try the
                  #Golub and Pereyra algorithm for the solution of a nonlinear least squares
                  #problem which assumes a number of the parameters are linear.
                  #Also, add a higher tolerance (1e-04 Vs 1e-05).
                              error = function(e) nls(percentage ~ SSlogis(time, Asym, xmid, scal),
                              data = df[, c('percentage', 'time')], algorithm = 'plinear',
                              control = nls.control(tol = 1e-04, warnOnly = T) ) )
                min = min(df$time)
                max = max(df$time)
                data = predict(
                    model,
                    newdata = data.frame(time = seq(
                        min, max, 0.01
                    )),
                    type = "l"
                )
                result = data.frame(
                    time = seq(min, max, 0.01),
                    percentage = data,
                    group = name
                )
                t = result %>%
                    mutate(
                        distance75 = abs(percentage - 0.75),
                        distance25 = abs(percentage - 0.25)
                    ) %>%
                    mutate(
                        min75 = min(distance75),
                        min25 = min(distance25),
                        t75 = distance75 == min75,
                        t25 = distance25 == min25
                    ) %>%
                    dplyr::select(group,
                                  time,
                                  percentage,
                                  t75,
                                  t25)  %>%
                    mutate(Cat1 = Cat1,
                           Cat2 = Cat2)
                
                return(t)
            } else{
                return(tibble())
            }
        }
        
        x=data%>%
            group_by(group,time,Cat1,Cat2)%>%
            summarise(percentage=mean(percentage)) 
        
        fitted_data = foreach(
            group = unique(x$group),
            .combine = 'rbind',
            .packages = c('tidyverse', 'foreach'),
            .errorhandling = 'remove'
        ) %do% {
            temp2 = foreach(
                Cat1 = unique(x$Cat1),
                .combine = 'rbind',
                .packages = c('tidyverse', 'foreach'),
                .errorhandling = 'remove'
            ) %do% {
                temp = foreach(
                    Cat2 = unique(x$Cat2),
                    .combine = 'rbind',
                    .packages = c('tidyverse', 'foreach'),
                    .errorhandling = 'remove'
                ) %do% {
                    
                    t = T25_75(df = x[x$group == group &
                                          x$Cat1 == Cat1 &
                                          x$Cat2 == Cat2, ], group, Cat1, Cat2)
                }
                temp
            }
            temp2
        }
        
        t = fitted_data %>% filter(t75 | t25) %>%
            gather('t', 'value', t25, t75) %>%
            filter(value) %>%
            dplyr::select(-percentage, -value) %>%
            spread(t, time)%>%
            mutate(Twidth = abs(t75 - t25))
        
        t %>% write_tsv(paste0(opt$out,
                               '/',
                               opt$output_file_base_name,
                               '_Twidth_categories.tsv'))
        p = ggplot(t) +
            geom_col(aes(' ', Twidth, fill = group), position = 'dodge') +
            ylab('Twidth') + xlab('')+facet_grid(Cat1~Cat2)+
            theme(axis.text.x = element_text(angle = 45, hjust=1))
        
        suppressMessages(ggsave(
            p,
            filename = paste0(opt$out,
                              '/',
                              opt$output_file_base_name,
                              '_Twidths_2_categories.pdf')
        ))
        
    }else{
        data = data %>%
            rbind(data %>%
                      mutate(Cat1 = '_All_'))
        
        x=data%>%
            group_by(group,time,Cat1)%>%
            summarise(percentage=mean(percentage)) 
        #T25_75 function
        T25_75 = function(df, name, EL) {
            model = tryCatch(
                nls(percentage ~ SSlogis(time, Asym, xmid, scal),
                    data = df[, c('percentage', 'time')]%>%
                        add_row(percentage=1,time=-10)%>%
                        add_row(percentage=0,time=10),
                    control = nls.control(maxiter = 100),
                    algorithm = 'port',
                    lower = c(Asym=1,xmid=0,scal=-1.5)
                ),
                #If the data cannot be fitted with a Gauss-Newton algorithm, try the
                #Golub and Pereyra algorithm for the solution of a nonlinear least squares
                #problem which assumes a number of the parameters are linear.
                #Also, add a higher tolerance (1e-04 Vs 1e-05).
                error = function(e)
                    nls(
                        percentage ~ SSlogis(time, Asym, xmid, scal),
                        data = df[, c('percentage', 'time')]%>%
                            add_row(percentage=1,time=-10)%>%
                            add_row(percentage=0,time=10),
                        algorithm = 'plinear',
                        control = nls.control(maxiter = 100,tol = 1e-04, warnOnly = T)
                    )
            )
            min = min(df$time)
            max = max(df$time)
            data = predict(model,
                           newdata = data.frame(time = seq(min, max, 0.01)),
                           type = "l")
            result = data.frame(
                time = seq(min, max, 0.01),
                percentage = data,
                group = name
            )
            t = result %>%
                mutate(
                    distance75 = abs(percentage - 0.75),
                    distance25 = abs(percentage - 0.25)
                ) %>%
                mutate(
                    min75 = min(distance75),
                    min25 = min(distance25),
                    t75 = distance75 == min75,
                    t25 = distance25 == min25
                ) %>%
                dplyr::select(group, time, percentage, t75, t25)  %>%
                mutate(Cat1 = EL)
            
            return(t)
        }
        
        #calculate tresholds 25% 75% replication keeping in account early and late domains
        fitted_data = foreach(
            group = unique(x$group),
            .combine = 'rbind',
            .packages = c('tidyverse', 'foreach')
        ) %do% {
            temp = foreach(
                EL = unique(x$Cat1),
                .combine = 'rbind',
                .packages = c('tidyverse', 'foreach')
            ) %do% {
                t = T25_75(df = x[x$group == group &
                                      x$Cat1 == EL, ], group, EL)
            }
            temp
        }
        
        t = fitted_data %>% filter(t75 | t25) %>%
            gather('t', 'value', t25, t75) %>%
            filter(value) %>%
            dplyr::select(-percentage, -value) %>%
            spread(t, time) %>%
            mutate(Twidth = abs(t75 - t25))
        
        t %>% write_tsv(paste0(opt$out,
                               '/',
                               opt$output_file_base_name,
                               '_Twidth.tsv'))
        
        
        plot=ggplot(x) +
            geom_point(aes(time,percentage,color=group))+
            geom_line(data=fitted_data,aes(time,percentage),color='blue')+
            scale_x_reverse()+
            geom_vline(data=t,aes(xintercept=t25),color='red')+
            geom_vline(data=t,aes(xintercept=t75),color='red')+
            geom_text(data=t,aes(label=paste('TW\n',Twidth)),x=Inf,y=0.5, hjust=1)+
            facet_grid(group~Cat1)
        
        ncat=length(unique(x$Cat1))
        nbasen=length(unique(x$group))
        suppressMessages(ggsave(
            plot,
            filename = paste0(opt$out,
                              '/',
                              opt$output_file_base_name,
                              '_Twidths_extended.pdf'),width = 2.2*ncat,height = 4*nbasen
        ))
        
        
        p = ggplot(t) +
            geom_col(aes(Cat1, Twidth, fill = group), position = 'dodge') +
            ylab('Twidth') + xlab('')+
            theme(axis.text.x = element_text(angle = 45, hjust=1))
        
        suppressMessages(ggsave(
            p,
            filename = paste0(opt$out,
                              '/',
                              opt$output_file_base_name,
                              '_Twidths.pdf')
        ))
        
        t %>% write_tsv(paste0(opt$out,
                               '/',
                               opt$output_file_base_name,
                               '_Twidth.tsv'))
    }
    

print('done')
