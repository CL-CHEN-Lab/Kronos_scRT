#!/usr/local/bin/Rscript --slave
#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)
options(warn = 1)

option_list = list(
    make_option(
        c("-F", "--file"),
        type = "character",
        default = NULL,
        help = "Variability file produced by Kronos RT, if multiple files are provided they have to be separated by a comma",
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
        c("-c", "--cores"),
        type = "integer",
        default = 3,
        help = "Numbers of parallel jobs to run [default= %default] ",
        metavar = "integer"
    )
)

#recover inputs
opt = parse_args(object = OptionParser(option_list = option_list))

suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))
suppressPackageStartupMessages(library(GenomicRanges, quietly = TRUE))
suppressPackageStartupMessages(library(ggpubr, quietly = TRUE))

if (!'file' %in% names(opt)) {
    stop("Variability file must be provided. See script usage (--help)")
}

opt$file = str_split(opt$file, ',')[[1]]

#load file
data <-
    foreach(file = opt$file,
            .combine = 'rbind',
            .packages = 'tidyverse') %do% {
                read_tsv(file, col_types = cols()) %>%
                    filter(Cat_RT != 'All')
            }


if ('regions' %in% names(opt)) {
    Annotation_file = read_tsv(opt$regions,
                               col_names = c('chr', 'start', 'end', 'annotation')) %>%
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
    hits = findOverlaps(Bin, Annotation_file)
    
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
    
    
    if (opt$both_annotations) {
        if ('regions2' %in% opt) {
            
            Annotation_file = read_tsv(opt$regions2,
                                       col_names = c('chr', 'start', 'end', 'annotation')) %>%
                makeGRangesFromDataFrame(
                    keep.extra.columns = T,
                    seqnames.field = 'chr',
                    end.field = 'end',
                    start.field = 'start'
                )
            
            hits = findOverlaps(Bin, Annotation_file)
            
            overlaps <-
                pintersect(Annotation_file[subjectHits(hits)], Bin[queryHits(hits)])
            
            changed = as_tibble(Bin[queryHits(hits)])
            not_changed = as_tibble(Bin[-queryHits(hits)])
            
            Annotation = bind_rows(
                changed %>%
                    mutate(
                        size = width(overlaps),
                        Cat_RT = overlaps$annotation
                    ) %>%
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
            
        }
        
        data = data %>%
            rbind(data %>%
                      mutate(Cat_RT = 'ALL'),
                  data %>%
                      mutate(old_Cat_RT = 'ALL'))
        
        # New Category and RT
        #calculate tresholds 25% 75% replication keeping in account early and late domains
        fitted_data = foreach(
            basename = unique(data$basename),
            .combine = 'rbind',
            .packages = c('tidyverse', 'foreach')
        ) %do% {
            temp2 = foreach(
                New = unique(data$Cat_RT),
                .combine = 'rbind',
                .packages = c('tidyverse', 'foreach')
            ) %do% {
                temp = foreach(
                    RT = unique(data$old_Cat_RT),
                    .combine = 'rbind',
                    .packages = c('tidyverse', 'foreach')
                ) %do% {
                    T25_75 = function(df, name, RT, New) {
                        if (length(df$chr) != 0) {
                            model = nls(percentage ~ SSlogis(time, Asym, xmid, scal),
                                        data = df[, c('percentage', 'time')])
                            min = -10
                            max = 10
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
                    
                    t = T25_75(df = data[data$basename == basename &
                                             data$Cat_RT == New &
                                             data$old_Cat_RT == RT, ], basename, RT, New)
                }
                temp
            }
            temp2
        }
        
        t = fitted_data %>% filter(t75 | t25) %>%
            gather('t', 'value', t25, t75) %>%
            filter(value) %>%
            dplyr::select(-percentage, -value) %>%
            spread(t, time)
        
        counts = data %>% select(chr, start, Cat_RT, old_Cat_RT, basename) %>%
            unique() %>%
            group_by(Cat_RT, old_Cat_RT, basename) %>%
            summarise(n = n())
        
        p = ggplot() +
            stat_density_2d(
                data = data,
                geom = "polygon",
                aes(
                    alpha = ..level..,
                    fill = basename,
                    y = percentage,
                    x = time
                )
            ) +
            geom_line(
                data = fitted_data,
                aes(y = percentage, x = time),
                color = 'blue',
                inherit.aes = F
            ) +
            geom_hline(yintercept = c(0.75, 0.25), color = 'yellow') +
            geom_vline(data = t,
                       aes(xintercept = t25),
                       color = 'red') +
            geom_vline(data = t,
                       aes(xintercept = t75),
                       color = 'red') +
            geom_text(
                data = t,
                aes(
                    label = paste('Twidth ~', round((t25 - t75), 1), 'h'),
                    x = (t25 + (t75 - t25) / 2),
                    y = Inf,
                    vjust = 1
                ),
                color = 'black',
                inherit.aes = F
            ) +
            facet_grid(basename + old_Cat_RT ~ Cat_RT) +
            ggplot2::scale_x_reverse() +
            scale_y_continuous(labels = scales::percent_format()) +
            ylab('Percentage of cells') +
            geom_text(data = counts,
                      aes(
                          label = paste('N bins: ', n),
                          x = 0,
                          y = -Inf
                      ),
                      vjust = 0)
        
        suppressMessages(ggsave(
            p,
            filename = paste0(
                opt$out,
                '/',
                paste(unique(data$basename), collapse = '_'),
                '_variability_2_categories.pdf'
            )
        ))
        
        suppressMessages(ggsave(
            p,
            filename = paste0(
                opt$out,
                '/',
                paste(unique(data$basename), collapse = '_'),
                '_variability_2_categories.jpg'
            )
        ))
        
        
    }
    
}

data = data %>%
    rbind(data %>%
              mutate(Cat_RT = 'ALL'))

fitted_data = foreach(
    basename = unique(data$basename),
    .combine = 'rbind',
    .packages = c('tidyverse', 'foreach')
) %do% {
    temp = foreach(
        EL = unique(data$Cat_RT),
        .combine = 'rbind',
        .packages = c('tidyverse', 'foreach')
    ) %do% {
        T25_75 = function(df, name, EL) {
            model = nls(percentage ~ SSlogis(time, Asym, xmid, scal),
                        data = df[, c('percentage', 'time')])
            min = -10
            max = 10
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
        
        t = T25_75(df = data[data$basename == basename &
                                 data$Cat_RT == EL, ], basename, EL)
    }
    temp
}

t = fitted_data %>% filter(t75 | t25) %>%
    gather('t', 'value', t25, t75) %>%
    filter(value) %>%
    dplyr::select(-percentage, -value) %>%
    spread(t, time)

counts = data %>% select(chr, start, Cat_RT, basename) %>%
    unique() %>%
    group_by(Cat_RT, basename) %>%
    summarise(n = n())

p = ggplot() +
    stat_density_2d(
        data = data,
        geom = "polygon",
        aes(
            alpha = ..level..,
            fill = basename,
            y = percentage,
            x = time
        )
    ) +
    geom_line(
        data = fitted_data,
        aes(y = percentage, x = time),
        color = 'blue',
        inherit.aes = F
    ) +
    geom_hline(yintercept = c(0.75, 0.25), color = 'yellow') +
    geom_vline(data = t,
               aes(xintercept = t25),
               color = 'red') +
    geom_vline(data = t,
               aes(xintercept = t75),
               color = 'red') +
    geom_text(
        data = t,
        aes(
            label = paste('Twidth ~', round((t25 - t75), 1), 'h'),
            x = (t25 + (t75 - t25) / 2),
            y = Inf,
            vjust = 1
        ),
        color = 'black',
        inherit.aes = F
    ) +
    facet_grid(basename ~ Cat_RT) +
    ggplot2::scale_x_reverse() +
    scale_y_continuous(labels = scales::percent_format()) +
    ylab('Percentage of cells') +
    geom_text(data = counts,
              aes(
                  label = paste('N bins: ', n),
                  x = 0,
                  y = -Inf
              ),
              vjust = 0)

suppressMessages(ggsave(
    p,
    filename = paste0(
        opt$out,
        '/',
        paste(unique(data$basename), collapse = '_'),
        '_variability_1_category.pdf'
    )
))

suppressMessages(ggsave(
    p,
    filename = paste0(
        opt$out,
        '/',
        paste(unique(data$basename), collapse = '_'),
        '_variability_1_category.pdf'
    )
))
