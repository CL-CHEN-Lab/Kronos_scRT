#!/usr/local/bin/Rscript --slave
#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)
options(warn = 1,scipen = 999)

option_list = list(
    make_option(
        c("-C", "--CNV"),
        type = "character",
        default = NULL,
        help = "scCNV file produced by Kronos RT, if multiple files are provided they have to be separated by a comma",
        metavar = "character"
    ),make_option(
        c("-T", "--RT"),
        type = "character",
        default = NULL,
        help = "RT file produced by Kronos RT, if multiple files are provided they have to be separated by a comma",
        metavar = "character"
    ),
    make_option(
        c("-R", "--regions"),
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
        c("-r", "--regions2"),
        type = "character",
        default = NULL,
        help = "Second genome annotation. chr<TAB>start<TAB>end<TAB>annotation. No header. If b option is activated it substiutes the RT division.",
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

suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))
suppressPackageStartupMessages(library(GenomicRanges, quietly = TRUE))
suppressPackageStartupMessages(library(ggpubr, quietly = TRUE))

if (!'CNV' %in% names(opt)) {
    stop("scCNV file must be provided. See script usage (--help)")
}
if (!'RT' %in% names(opt)) {
    stop("RT file must be provided. See script usage (--help)")
}
if (str_extract(opt$out, '.$') != '/') {
    opt$out = paste0(opt$out, '/')
}

system(paste0('mkdir -p ', opt$out))

opt$CNV = str_split(opt$CNV, ',')[[1]]
opt$RT = str_split(opt$RT, ',')[[1]]

#load files
data <-
    foreach(file = opt$CNV,
            .combine = 'rbind',
            .packages = 'tidyverse') %do% {
                read_tsv(file, col_types = cols())
            }

RT_file <-
    foreach(file = opt$RT,
            .combine = 'rbind',
            .packages = 'tidyverse') %do% {
                read_tsv(file, col_types = cols())
            }

data = data %>%
    dplyr::select(basename, chr, start, end, Rep, newIndex,mean_CN) %>%
    inner_join(RT_file, by = c("basename", "chr", "start", "end"))  %>%
    mutate(
        RT = 10 * RT ,
        time = round(10 - RT - 10*mean_CN,1),
        Cat_RT = case_when(
            RT < 2 ~ '1 - Very Early',
            RT >= 2 & RT < 4  ~ '2 - Early',
            RT >= 4 & RT < 6 ~ '3 - Mid ',
            RT >= 6 & RT < 8 ~ '4 - Late',
            RT >= 8  ~ '5 - Very Late'
        ),
        Cat_RT = factor(
            Cat_RT,
            levels = c(
                '0 - All',
                '1 - Very Early',
                '2 - Early',
                '3 - Mid ',
                '4 - Late',
                '5 - Very Late'
            )
        )
        )


if ('regions' %in% names(opt)) {
    Annotation_file = read_tsv(opt$regions,
                               col_names = c('chr', 'start', 'end', 'annotation'), col_types = cols()) %>%
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
    
    #convert bins that have a notation and those that don't into df
    changed = as_tibble(Bin[queryHits(hits)])
    not_changed = as_tibble(Bin[-queryHits(hits)])
    
    #Based on the overla define the predominant notation of each bin.
    #A category to be chosen has to be predominant (at least 60% of the tatal overlaps in the bin)
    Annotation = bind_rows(
        changed %>%
            mutate(size = width(overlaps),
                   Cat_RT = overlaps$annotation) %>%
            group_by(seqnames, start, end, Cat_RT) %>%
            summarise(n = sum(size)) %>%
            ungroup() %>%
            group_by(seqnames, start, end) %>%
            mutate(
                select = (n / sum(n) >= 0.6),
                Cat_RT = ifelse(any(select), Cat_RT, '_Unknown_')
            ) %>%
            dplyr::select(seqnames, start, end, Cat_RT),
        not_changed %>% mutate(Cat_RT = '_Unknown_') %>%
            dplyr::select(seqnames, start, end, Cat_RT)
    )
    Annotation = Annotation %>%
        ungroup() %>%
        mutate(seqnames = as.character(seqnames))
    
    data = data %>%
        mutate(old_Cat_RT = Cat_RT) %>%
        dplyr::select(-Cat_RT) %>%
        inner_join(Annotation, by = c("chr" = "seqnames", "start", "end"))
    
    if (opt$both_annotations) {
        if ('regions2' %in% names(opt)) {
            
            Annotation_file = read_tsv(opt$regions2,
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
            
            changed = as_tibble(Bin[queryHits(hits)])
            not_changed = as_tibble(Bin[-queryHits(hits)])
            
            Annotation = bind_rows(
                changed %>%
                    mutate(
                        size = width(overlaps),
                        old_Cat_RT = overlaps$annotation
                    ) %>%
                    group_by(seqnames, start, end, old_Cat_RT) %>%
                    summarise(n = sum(size)) %>%
                    ungroup() %>%
                    group_by(seqnames, start, end) %>%
                    mutate(
                        select = (n / sum(n) >= 0.6),
                        old_Cat_RT = ifelse(any(select), old_Cat_RT, '_Unknown_')
                    ) %>%
                    dplyr::select(seqnames, start, end, old_Cat_RT),
                not_changed %>% mutate(old_Cat_RT = '_Unknown_') %>%
                    dplyr::select(seqnames, start, end, old_Cat_RT)
            )
            Annotation = Annotation %>%
                ungroup() %>%
                mutate(seqnames = as.character(seqnames))
            
            data = data %>%
                mutate(old_Cat_RT = Cat_RT) %>%
                dplyr::select(-old_Cat_RT) %>%
                inner_join(Annotation, by = c("chr" = "seqnames", "start", "end"))
            
        }
        
        data = data %>%
            rbind(data %>%
                      mutate(Cat_RT = '_ALL_'),
                  data %>%
                      mutate(old_Cat_RT = '_ALL_'),
                  data %>%
                      mutate(Cat_RT = '_ALL_',
                             old_Cat_RT = '_ALL_'))
        
        # New Category and RT
        #calculate tresholds 25% 75% replication keeping in account early and late domains
        
        T25_75 = function(df, name, RT, New) {
            if (length(df$basename) != 0) {
                model = tryCatch( nls(percentage ~ SSlogis(time, Asym, xmid, scal),
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
                    basename = name
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
                    dplyr::select(basename,
                                  time,
                                  percentage,
                                  t75,
                                  t25)  %>%
                    mutate(Cat_RT = New,
                           old_Cat_RT = RT)
                
                return(t)
            } else{
                return(tibble())
            }
        }
        
        x=data%>%
            group_by(basename,time,Cat_RT,old_Cat_RT)%>%
            summarise(percentage=mean(Rep))
        
        
        fitted_data = foreach(
            basename = unique(x$basename),
            .combine = 'rbind',
            .packages = c('tidyverse', 'foreach'),
            .errorhandling = 'remove'
        ) %do% {
            temp2 = foreach(
                New = unique(x$Cat_RT),
                .combine = 'rbind',
                .packages = c('tidyverse', 'foreach'),
                .errorhandling = 'remove'
            ) %do% {
                temp = foreach(
                    RT = unique(x$old_Cat_RT),
                    .combine = 'rbind',
                    .packages = c('tidyverse', 'foreach'),
                    .errorhandling = 'remove'
                ) %do% {
                    
                    t = T25_75(df = x[x$basename == basename &
                                          x$Cat_RT == New &
                                          x$old_Cat_RT == RT, ], basename, RT, New)
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
        
        t %>%rename(Cat_RT='Annotation1'  , old_Cat_RT='Annotation2')%>%
        write_tsv(paste0(opt$out,
                               '/',
                               opt$output_file_base_name,
                               '_Twidth_2_categories.tsv'))
        
        
        p = ggplot(t) +
            geom_col(aes(' ', Twidth, fill = basename), position = 'dodge') +
            ylab('Twidth') + xlab('')+facet_grid(Cat_RT~old_Cat_RT)
        
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
                      mutate(Cat_RT = '_All_'))
        
        x=data%>%
            group_by(basename,time,Cat_RT)%>%
            summarise(percentage=mean(Rep))
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
                basename = name
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
                dplyr::select(basename, time, percentage, t75, t25)  %>%
                mutate(Cat_RT = EL)
            
            return(t)
        }
        
        #calculate tresholds 25% 75% replication keeping in account early and late domains
        fitted_data = foreach(
            basename = unique(x$basename),
            .combine = 'rbind',
            .packages = c('tidyverse', 'foreach')
        ) %do% {
            temp = foreach(
                EL = unique(x$Cat_RT),
                .combine = 'rbind',
                .packages = c('tidyverse', 'foreach')
            ) %do% {
                t = T25_75(df = x[x$basename == basename &
                                      x$Cat_RT == EL, ], basename, EL)
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
            geom_point(aes(time,percentage,color=basename))+
            geom_line(data=fitted_data,aes(time,percentage),color='blue')+
            scale_x_reverse()+
            geom_vline(data=t,aes(xintercept=t25),color='red')+
            geom_vline(data=t,aes(xintercept=t75),color='red')+
            geom_text(data=t,aes(label=paste('TW\n',Twidth)),x=Inf,y=0.5, hjust=1)+
            facet_grid(basename~Cat_RT)
        
        ncat=length(unique(x$Cat_RT))
        nbasen=length(unique(x$basename))
        suppressMessages(ggsave(
            plot,
            filename = paste0(opt$out,
                              '/',
                              opt$output_file_base_name,
                              '_Twidths_extended.pdf'),width = 2.2*ncat,height = 8*nbasen
        ))
        
        
        p = ggplot(t) +
            geom_col(aes(Cat_RT, Twidth, fill = basename), position = 'dodge') +
            ylab('Twidth') + xlab('')
        
        suppressMessages(ggsave(
            p,
            filename = paste0(opt$out,
                              '/',
                              opt$output_file_base_name,
                              '_Twidths.pdf')
        ))
    }
    
}else{

data = data %>%
    rbind(data %>%
              mutate(Cat_RT = '0 - All'))

x=data%>%
    group_by(basename,time,Cat_RT)%>%
    summarise(percentage=mean(Rep))

x %>%
    write_tsv(paste0(
        opt$out,
        '/',
        opt$output_file_base_name,
        '_scRT_variability.tsv'
    ))

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
        basename = name
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
        dplyr::select(basename, time, percentage, t75, t25)  %>%
        mutate(Cat_RT = EL)
    
    return(t)
}

#calculate tresholds 25% 75% replication keeping in account early and late domains
fitted_data = foreach(
    basename = unique(x$basename),
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach'),
    .errorhandling = 'remove'
) %do% {
    temp = foreach(
        EL = unique(x$Cat_RT),
        .combine = 'rbind',
        .packages = c('tidyverse', 'foreach'),
        .errorhandling = 'remove'
    ) %do% {
        t = T25_75(df = x[x$basename == basename &
                              x$Cat_RT == EL, ], basename, EL)
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
    geom_point(aes(time,percentage,color=basename))+
    geom_line(data=fitted_data,aes(time,percentage),color='blue')+
    scale_x_reverse()+
    geom_vline(data=t,aes(xintercept=t25),color='red')+
    geom_vline(data=t,aes(xintercept=t75),color='red')+
    geom_text(data=t,aes(label=paste('TW\n',Twidth)),x=Inf,y=0.5, hjust=1)+
    facet_grid(basename~Cat_RT)

ncat=length(unique(x$Cat_RT))
nbasen=length(unique(x$basename))
suppressMessages(ggsave(
    plot,
    filename = paste0(opt$out,
                      '/',
                      opt$output_file_base_name,
                      '_Twidths_extended.pdf'),width = 2.2*ncat,height = 8*nbasen
))


p = ggplot(t) +
    geom_col(aes(Cat_RT, Twidth, fill = basename), position = 'dodge') +
    ylab('Twidth') + xlab('')

suppressMessages(ggsave(
    p,
    filename = paste0(opt$out,
                      '/',
                      opt$output_file_base_name,
                      '_Twidths.pdf')
))

}
print('done')
