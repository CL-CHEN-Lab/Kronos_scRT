#!/usr/local/bin/Rscript
#parse input
suppressPackageStartupMessages(library(optparse, quietly = TRUE))

options(stringsAsFactors = FALSE)
options(warn = 1,scipen = 999)

option_list = list(
    make_option(
        c("-F", "--file"),
        type = "character",
        default = NULL,
        help = "Variability file with groups produced by Kronos annotate, if multiple files are provided they have to be separated by a comma",
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
        c("-p", "--pval"),
        type = "logical",
        action = "store_true",
        default = F,
        help = "Bootstrap pValue for the difference in TW between groups (it works only with one annotation) [default= %default]",
        metavar = "logical"
    ),make_option(
        c("-a", "--padj_method"),
        type = "character",
        default = 'none',
        help = "holm, hochberg, hommel, bonferroni, BY (Benjamini & Yekutieli ),fdr (false discovery rate),none [default= %default]",
        metavar = "character"
    ),
    make_option(
        c("-A", "--Annotation_to_use_for_pval"),
        type = "numeric",
        default = 1,
        help = "Annotatation to use to calculate pvalues (1=Cat1,2=Cat2,3=Cat1_Cat2) [default= %default]",
        metavar = "numeric"
    ),
    make_option(
        c("-B", "--between_groups"),
        type = "logical",
        action = "store_true",
        default = F,
        help = "If selected pvalues are calculated between samples instead of within samples.",
        metavar = "logical"
    ),
    make_option(
        c("-G", "--pairs_to_test"),
        type = "character",
        help = "Pairs of groups for which to calculate the pvalue. Groups in a pair have to be separeted by a comma while pairs are separated by a semicolun eg. A,B;A,C",
        metavar = "character"
    ),
    make_option(
        c("-H", "--pval_alternative_hypotesis"),
        default = 'two.sided',
        help = "greater, lower, two.sided. This option is active only if the option G is porided [default= %default]",
        metavar = "character",
        type = "character"
        
    ),
    make_option(
        c("-i", "--number_of_iterations"),
        default = 10000,
        help = "Number of iterations to calculate boostrap [default= %default]",
        metavar = "numeric",
        type = "numeric"
    ),
    make_option(
        c("-c", "--cores"),
        default = 3,
        help = "Number of for bootstrapping [default= %default]",
        metavar = "numeric",
        type = "numeric"
    )
)

#recover inputs
opt = parse_args(object = OptionParser(option_list = option_list))

#load packages
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(foreach, quietly = TRUE))
suppressPackageStartupMessages(library(doSNOW, quietly = TRUE))

#set plotting theme
theme_set(theme_bw())

#check inputs
if (!'file' %in% names(opt)) {
    stop("Variability file must be provided. See script usage (--help)")
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


if ('Cat2' %in% names(data) & !opt$pval) {
    
    # add all bins together 
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
                                  start = c(Asym=1,xmid=0,scal=-0.5)
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
        summarise(percentage=mean(percentage))%>%
        ungroup() 
    
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
    # count number of unique bins in each category
    unique_bins=data%>%dplyr::select(group,chr,start,end,Cat1,Cat2)%>%unique()%>%group_by(group,Cat1,Cat2)%>%summarise(`N of bins`=n())
    
    #merge with t
    t=t%>%left_join(unique_bins, by = c("group", "Cat1", "Cat2"))
    
    # write file 
    t %>% write_tsv(paste0(opt$out,
                           '/',
                           opt$output_file_base_name,
                           '_Twidth_categories.tsv'))
    

    
    p = ggplot(t) +
        geom_col(aes(' ', Twidth, fill = group), position = 'dodge') +
        ylab('Twidth') + xlab('')+facet_grid(Cat1~Cat2)+
        theme(axis.text.x = element_text(angle = 45, hjust=1))+
        geom_text(aes(' ',Twidth/2, label=paste0('n bins:\n',`N of bins`)),angle=90, hjust=0.5,size=2, vjust=0.5)
    
    suppressMessages(ggsave(
        p,
        filename = paste0(opt$out,
                          '/',
                          opt$output_file_base_name,
                          '_Twidths_2_categories.pdf')
    ))
    
}else{
    #add all bins together 
    data = data %>%
        rbind(data %>%
                  mutate(Cat1 = '_All_'))
    
    #select annotation to use for plotting and pvalue
    if (opt$Annotation_to_use_for_pval==2){
        data=data%>%dplyr::mutate(Cat1=Cat2)
    }else if(opt$Annotation_to_use_for_pval==3){
        data=data%>%dplyr::mutate(Cat1=paste(Cat1,Cat2,sep = '_'))
    }
    
    x=data%>%
        group_by(group,time,Cat1)%>%
        summarise(percentage=mean(percentage,na.rm=T)) %>%
        ungroup()
    #T25_75 function
    T25_75 = function(df, name, EL) {
        model = tryCatch(
            nls(percentage ~ SSlogis(time, Asym, xmid, scal),
                data = df[, c('percentage', 'time')]%>%
                    add_row(percentage=1,time=-10)%>%
                    add_row(percentage=0,time=10),
                control = nls.control(maxiter = 100),
                algorithm = 'port',
                start = c(Asym=1,xmid=0,scal=-0.5)
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
            .packages = c('tidyverse', 'foreach'),
            .errorhandling = 'remove'
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
    # count number of unique bins in each category
    unique_bins=data%>%dplyr::select(group,chr,start,end,Cat1)%>%unique()%>%group_by(group,Cat1)%>%summarise(`N of bins`=n())
    
    #merge with t and asssign 
    t = t %>% left_join(unique_bins, by = c("group", "Cat1")) 
    
    # write fiel 
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
    ####pvalue
    if(opt$pval & ! opt$between_groups){
        
        #identily pairs to test
        if ('pairs_to_test' %in% names(opt)){
            
            # if provided by the user separate pairs using ; and groups using ,
            groups= tibble( pairs=str_split(opt$pairs_to_test,pattern = ';')[[1]])%>%
                separate(pairs,into =c('Group1','Group2'),sep = ',' )
            
            
        }else{
            
            # if not provided create all possible combinations
            groups=sort(unique(data$Cat1))
            # remove all since the other groups are subgroups of it
            groups=groups[groups!='_All_']
            
            groups=foreach(i=1:(length(groups)-1),.combine = 'rbind')%:%
                foreach(h=(i+1):length(groups),.combine = 'rbind')%do%{
                    
                    tibble(
                        Group1=groups[i],
                        Group2=groups[h]
                    )
                    
                }
            
        }
        
        # base of the input inverte the order of element or calculate absolute
        stat_type = function(opt, x, y) {
            # if pairs are not provided a two sided pval is provided
            if ('pairs_to_test' %in% names(opt)) {
                return(switch (
                    opt$pval_alternative_hypotesis,
                    'greater' = x - y,
                    'lower' = y - x,
                    'two.sided' = abs(y - x),
                    stop('wrong alternative hypotesys provided')
                ))
            } else{
                return(abs(y - x))
            }
        }
        pval = foreach(
            g = unique(data$group),
            .combine = 'rbind',
            .packages = c('tidyverse', 'foreach'),
            .export = c('stat_type', 'T25_75')
        ) %do% {
            pval_TW = foreach (
                i = 1:length(groups$Group1),
                .combine = 'rbind' ,
                .packages = c('tidyverse', 'foreach'),
                .export = c('stat_type', 'T25_75'),
                .errorhandling = "remove"
            ) %do% {
                #remove colum Cat1 from data and filter for line name
                data_g = data %>%
                    filter(group == g)
                
                #bins belonging to the first group
                Group1 = data_g %>%
                    filter(Cat1 == groups$Group1[i])
                #bins belonging to the second group
                Group2 = data_g %>%
                    filter(Cat1 == groups$Group2[i])
                
                #n of bins belonging to the first group
                Group1_n = Group1%>%dplyr::select(chr,start,end)%>%unique()%>%pull(chr)%>%length()
                #n of bins belonging to the second group
                Group2_n = Group2%>%dplyr::select(chr,start,end)%>%unique()%>%pull(chr)%>%length()
                
                # calculate percent of replicated at each time
                G1x = Group1 %>%
                    group_by(group, time, Cat1) %>%
                    summarise(percentage = mean(percentage)) %>%
                    ungroup()
                
                G2x = Group2 %>%
                    group_by(group, time, Cat1) %>%
                    summarise(percentage = mean(percentage)) %>%
                    ungroup()
                
                
                # calculate twidth group1
                T_G1x = T25_75(df = G1x, g, groups$Group1[i]) %>%
                    filter(t75 | t25) %>%
                    gather('t', 'value', t25, t75) %>%
                    filter(value) %>%
                    dplyr::select(-percentage,-value) %>%
                    spread(t, time) %>%
                    mutate(Twidth = abs(t75 - t25)) %>%
                    pull(Twidth)
                
                #calculate TW group2
                T_G2x = T25_75(df = G2x, g, groups$Group2[i]) %>%
                    filter(t75 | t25) %>%
                    gather('t', 'value', t25, t75) %>%
                    filter(value) %>%
                    dplyr::select(-percentage,-value) %>%
                    spread(t, time) %>%
                    mutate(Twidth = abs(t75 - t25)) %>%
                    pull(Twidth)
                
                Real_difference = stat_type(opt, T_G1x, T_G2x)
                
                
                cl = makeCluster(opt$cores)
                registerDoSNOW(cl)
                
                Boot_Strapped = foreach(
                    index = 1:opt$number_of_iterations,
                    .combine = '+',
                    .packages = c('tidyverse'),
                    .export = c('stat_type', 'T25_75'),
                    .errorhandling = "remove"
                ) %dopar% {
                    
                    # bins belonging to both groups
                    Groups = data_g %>% filter(Cat1 %in% c(groups$Group1[i], groups$Group2[i]))%>%
                        dplyr::select(group,chr,start,end)%>%
                        unique()
                    
                    #randomly select bins
                    Group1 = Groups[sample(
                        x = 1:length(Groups$chr),
                        size = Group1_n,
                        replace = T
                    ),] %>%
                        inner_join(data_g, by = c("group","chr", "start", "end"))%>%
                        mutate(Cat1 = groups$Group1[i]) 
                    
                    Group2 = Groups[sample(
                        x = 1:length(Groups$chr),
                        size = Group2_n,
                        replace = T
                    ),]  %>%
                        inner_join(data_g, by = c("group","chr", "start", "end"))%>%
                        mutate(Cat1 = groups$Group2[i])
                    
                    # calculate percent of replicated at each time
                    G1x = Group1 %>%
                        group_by(group, time, Cat1) %>%
                        summarise(percentage = mean(percentage)) %>%
                        ungroup()
                    
                    G2x = Group2 %>%
                        group_by(group, time, Cat1) %>%
                        summarise(percentage = mean(percentage)) %>%
                        ungroup()
                    
                    T_G1x = T25_75(df = G1x, g, groups$Group1[i]) %>%
                        filter(t75 | t25) %>%
                        gather('t', 'value', t25, t75) %>%
                        filter(value) %>%
                        dplyr::select(-percentage,-value) %>%
                        spread(t, time) %>%
                        mutate(Twidth = abs(t75 - t25)) %>%
                        pull(Twidth)
                    
                    #calculate TW group2
                    T_G2x = T25_75(df = G2x, g, groups$Group2[i]) %>%
                        filter(t75 | t25) %>%
                        gather('t', 'value', t25, t75) %>%
                        filter(value) %>%
                        dplyr::select(-percentage,-value) %>%
                        spread(t, time) %>%
                        mutate(Twidth = abs(t75 - t25)) %>%
                        pull(Twidth)
                    
                    c(stat_type(opt, T_G1x, T_G2x) >=
                          Real_difference,1)
                    
                    
                }
                
                stopCluster(cl)
                
                #return pvalues
                tibble(
                    group = g,
                    Cat1 = groups$Group1[i],
                    Cat2 = groups$Group2[i],
                    #calculate pvalue over effective iterations
                    pval = Boot_Strapped[1]/Boot_Strapped[2],
                    # actual iterations 
                    iterations=Boot_Strapped[2]
                )
            }
            pval_TW
        }
        
    }else if (opt$pval &  opt$between_groups){
        
        #identily pairs to test
        if ('pairs_to_test' %in% names(opt)){
            
            # if provided by the user separate pairs using ; and groups using ,
            groups= tibble( pairs=str_split(opt$pairs_to_test,pattern = ';')[[1]])%>%
                separate(pairs,into =c('Group1','Group2'),sep = ',' )
            
            
        }else{
            
            # if not provided create all possible combinations
            groups=unique(data$group)
            
            groups=foreach(i=1:(length(groups)-1),.combine = 'rbind')%:%
                foreach(h=(i+1):length(groups),.combine = 'rbind')%do%{
                    
                    tibble(
                        Group1=groups[i],
                        Group2=groups[h]
                    )
                    
                }
            
        }
        
        # base of the input inverte the order of element or calculate absolute
        stat_type = function(opt, x, y) {
            # if pairs are not provided a two sided pval is provided
            if ('pairs_to_test' %in% names(opt)) {
                return(switch (
                    opt$pval_alternative_hypotesis,
                    'greater' = x - y,
                    'lower' = y - x,
                    'two.sided' = abs(y - x),
                    stop('wrong alternative hypotesys provided')
                ))
            } else{
                return(abs(y - x))
            }
        }
        
        pval = foreach(
            g = unique(data$Cat1),
            .combine = 'rbind',
            .packages = c('tidyverse', 'foreach'),
            .export = c('stat_type', 'T25_75')
        ) %do% {
            pval_TW = foreach (
                i = 1:length(groups$Group1),
                .combine = 'rbind' ,
                .packages = c('tidyverse', 'foreach'),
                .export = c('stat_type', 'T25_75'),
                .errorhandling = "remove"
            ) %do% {
                #remove colum Cat1 from data and filter for line name
                data_g = data %>%
                    filter(Cat1 == g,
                           group %in% c(groups$Group1[i], groups$Group2[i]))
                
                Bins_for_shuffling=data_g %>% dplyr::select(group,chr,start,end)%>%unique()
                
                #n of bins belonging to the first group
                Group1_n = sum(Bins_for_shuffling$group == groups$Group1[i])
                #n of bins belonging to the second group
                Group2_n = sum(Bins_for_shuffling$group == groups$Group2[i])
                
                #bins belonging to the first group
                Group1 = data_g %>%
                    filter(group == groups$Group1[i])
                #bins belonging to the second group
                Group2 = data_g %>%
                    filter(group == groups$Group2[i])
                
                # calculate percent of replicated at each time
                G1x = Group1 %>%
                    group_by(group, time, Cat1) %>%
                    summarise(percentage = mean(percentage)) %>%
                    ungroup()
                
                G2x = Group2 %>%
                    group_by(group, time, Cat1) %>%
                    summarise(percentage = mean(percentage)) %>%
                    ungroup()
                
                
                # calculate twidth group1
                T_G1x = T25_75(df = G1x, g, groups$Group1[i]) %>%
                    filter(t75 | t25) %>%
                    gather('t', 'value', t25, t75) %>%
                    filter(value) %>%
                    dplyr::select(-percentage,-value) %>%
                    spread(t, time) %>%
                    mutate(Twidth = abs(t75 - t25)) %>%
                    pull(Twidth)
                
                #calculate TW group2
                T_G2x = T25_75(df = G2x, g, groups$Group2[i]) %>%
                    filter(t75 | t25) %>%
                    gather('t', 'value', t25, t75) %>%
                    filter(value) %>%
                    dplyr::select(-percentage,-value) %>%
                    spread(t, time) %>%
                    mutate(Twidth = abs(t75 - t25)) %>%
                    pull(Twidth)
                
                Real_difference = stat_type(opt, T_G1x, T_G2x)
                
                
                cl = makeCluster(opt$cores)
                registerDoSNOW(cl)
                
                Boot_Strapped = foreach(
                    index = 1:opt$number_of_iterations,
                    .combine = '+',
                    .packages = c('tidyverse'),
                    .export = c('stat_type', 'T25_75'),
                    .errorhandling = "remove"
                ) %dopar% {
                    
                    Group1 = Bins_for_shuffling[sample(
                        x = 1:length(Bins_for_shuffling$chr),
                        size = Group1_n,
                        replace = T
                    ),] %>%
                        inner_join(data_g, by = c('group',"chr", "start", "end"))%>%
                        mutate(group=groups$Group1[i])
                    Group2 = Bins_for_shuffling[sample(
                        x = 1:length(Bins_for_shuffling$chr),
                        size = Group2_n,
                        replace = T
                    ),] %>%
                        inner_join(data_g, by = c('group',"chr", "start", "end"))%>%
                        mutate(group=groups$Group2[i])
                    
                    # calculate percent of replicated at each time
                    G1x = Group1 %>%
                        group_by(group, time, Cat1) %>%
                        summarise(percentage = mean(percentage)) %>%
                        ungroup()
                    
                    G2x = Group2 %>%
                        group_by(group, time, Cat1) %>%
                        summarise(percentage = mean(percentage)) %>%
                        ungroup()
                    
                    T_G1x = T25_75(df = G1x, g, groups$Group1[i]) %>%
                        filter(t75 | t25) %>%
                        gather('t', 'value', t25, t75) %>%
                        filter(value) %>%
                        dplyr::select(-percentage,-value) %>%
                        spread(t, time) %>%
                        mutate(Twidth = abs(t75 - t25)) %>%
                        pull(Twidth)
                    
                    #calculate TW group2
                    T_G2x = T25_75(df = G2x, g, groups$Group2[i]) %>%
                        filter(t75 | t25) %>%
                        gather('t', 'value', t25, t75) %>%
                        filter(value) %>%
                        dplyr::select(-percentage,-value) %>%
                        spread(t, time) %>%
                        mutate(Twidth = abs(t75 - t25)) %>%
                        pull(Twidth)
                    
                    c( stat_type(opt, T_G1x, T_G2x) >=
                           Real_difference,1)
                    
                    
                }
                
                stopCluster(cl)
                
                #return pvalues
                tibble(
                    Cat1 = g,
                    Group1 = groups$Group1[i],
                    Group2 = groups$Group2[i],
                    #calculate pvalue over effective iterations
                    pval = Boot_Strapped[1]/Boot_Strapped[2],
                    # actual iterations 
                    iterations=Boot_Strapped[2]
                )
            }
            pval_TW
        }
        
    }
    

    #convert Cat1 and group into factor
    t=t%>% 
        mutate(Cat1 =factor(Cat1, levels = unique(Cat1)),
              group = factor(group, levels = unique(group)))
    
    
    
    if(opt$pval & ! opt$between_groups){
        
        pval=pval%>%
            mutate(a=as.numeric(factor(Cat1,levels(t$Cat1))),
                   b=as.numeric(factor(Cat2,levels(t$Cat1))),
                   position=(a+b)/2,
                   x=ifelse(a>b,b,a),
                   xend=ifelse(a<b,b,a),
                   adj_pval=p.adjust(pval,method = opt$padj_method),
                   pval = ifelse(
                       pval == 0,
                       paste0(
                           '<',
                           format(
                               1 / iterations,
                               digits = 2,
                               scientific = T
                           )
                       ),pval),
                   adj_pval_label = ifelse(
                       adj_pval == 0,
                       paste0(
                           '<',
                           format(
                               1 / iterations,
                               digits = 2,
                               scientific = T
                           )
                       ),
                       format(
                           adj_pval,
                           digits = 2,
                           scientific = T
                       )
                   ),
                   Statistically=ifelse(
                       adj_pval < 0.05 & iterations > 2000,
                       'significant',
                       'non-significant'),
                   adj_pval = ifelse(
                       adj_pval == 0,
                       paste0(
                           '<',
                           format(
                               1 / iterations,
                               digits = 2,
                               scientific = T
                           )
                       ),
                       adj_pval
                   )
            )%>%
            group_by(group)%>%
            mutate( y=seq(1.1*max(t$Twidth),(1+n()/10)*max(t$Twidth),length.out = n()))%>%
            ungroup()
        if(opt$padj_method=='none'){
            
            pval%>% dplyr::select(group,Cat1,Cat2,pval,succesful_iterations=iterations)%>%write_tsv(paste0(opt$out,
                                                                                                                    '/',
                                                                                                                    opt$output_file_base_name,
                                                                                                                    '_Twidth_pvalues_within_groups.tsv'))
            
        }else{
        pval%>% dplyr::select(group,Cat1,Cat2,pval,adj_pval,succesful_iterations=iterations)%>%write_tsv(paste0(opt$out,
                                                                                                                '/',
                                                                                                                opt$output_file_base_name,
                                                                                                                '_Twidth_pvalues_within_groups.tsv'))
        }
        p=ggplot(t) +
            geom_col(aes(Cat1, Twidth, fill = group), position = 'dodge') +
            ylab('Twidth') + xlab('')+
            theme(axis.text.x = element_text(angle = 45, hjust=1))+
            geom_text(aes(Cat1,Twidth/2, label=paste0('n bins: ',`N of bins`)),angle=90, hjust=0.5,size=2, vjust=0.5)+
            geom_text(data=pval, aes(x=position,y=y,label=adj_pval_label,color=Statistically),vjust=-0.25)+
            geom_segment(data=pval,aes(
                x=x,
                xend=xend,
                y=y,
                yend=y,color=Statistically))+
            facet_grid(~group)+
            scale_color_manual(values = c('significant'='red','non-significant'='black'))
        
    }else if(opt$pval & opt$between_groups){
        
        pval=pval%>%
            mutate(a=as.numeric(factor(Group1,levels(t$group))),
                   b=as.numeric(factor(Group2,levels(t$group))),
                   position=(a+b)/2,
                   x=ifelse(a>b,b,a),
                   xend=ifelse(a<b,b,a),
                   adj_pval=p.adjust(pval,method = opt$padj_method),
                   pval = ifelse(
                       pval == 0,
                       paste0(
                           '<',
                           format(
                               1 / iterations,
                               digits = 2,
                               scientific = T
                           )
                       ),pval),
                   adj_pval_label = ifelse(
                       adj_pval == 0,
                       paste0(
                           '<',
                           format(
                               1 / iterations,
                               digits = 2,
                               scientific = T
                           )
                       ),
                       format(
                           adj_pval,
                           digits = 2,
                           scientific = T
                       )
                   ),
                   Statistically=ifelse(
                       adj_pval < 0.05 & iterations > 2000,
                       'significant',
                       'non-significant'),
                   adj_pval = ifelse(
                       adj_pval == 0,
                       paste0(
                           '<',
                           format(
                               1 / iterations,
                               digits = 2,
                               scientific = T
                           )
                       ),
                       adj_pval
                   )
            )%>%
            group_by(Cat1)%>%
            mutate( y=seq(1.1*max(t$Twidth),(1+n()/10)*max(t$Twidth),length.out = n()))
        if(opt$padj_method=='none'){
            pval%>% dplyr::select(Group1,Group2,Cat1,pval,succesful_iterations=iterations)%>% write_tsv(paste0(opt$out,
                                                                                                                        '/',
                                                                                                                        opt$output_file_base_name,
                                                                                                                        '_Twidth_pvalues_between_groups.tsv'))
        }else{
        pval%>% dplyr::select(Group1,Group2,Cat1,pval,adj_pval,succesful_iterations=iterations)%>% write_tsv(paste0(opt$out,
                                                                                                                    '/',
                                                                                                                    opt$output_file_base_name,
                                                                                                                    '_Twidth_pvalues_between_groups.tsv'))
        }
        
        p= ggplot(t) +
            geom_col(aes(group, Twidth, fill = group), position = 'dodge') +
            ylab('Twidth') + xlab('')+
            theme(axis.text.x = element_text(angle = 45, hjust=1))+
            geom_text(aes(group,Twidth/2, label=paste0('n bins: ',`N of bins`)),angle=90, hjust=0.5,size=2, vjust=0.5)+
            facet_grid(~Cat1)+
            geom_text(data=pval, aes(x=position,y=y,label=adj_pval_label,color=Statistically),vjust=-0.25)+
            geom_segment(data=pval,aes(
                x=x,
                xend=xend,
                y=y,
                yend=y,color=Statistically))+
            scale_color_manual(values = c('significant'='red','non-significant'='black'))
        
    }else{
        p = ggplot(t) +
            geom_col(aes(Cat1, Twidth, fill = group), position = 'dodge') +
            ylab('Twidth') + xlab('')+
            theme(axis.text.x = element_text(angle = 45, hjust=1))+
            geom_text(aes(Cat1,Twidth/2, label=paste0('n bins: ',`N of bins`)),angle=90, hjust=0.5, vjust=0.5)+
            facet_grid(~group)
            
    }
    suppressMessages(ggsave(
        p,
        filename = paste0(opt$out,
                          '/',
                          opt$output_file_base_name,
                          '_Twidths.pdf'),width = 14,height = 8
    ))
    
    t %>% write_tsv(paste0(opt$out,
                           '/',
                           opt$output_file_base_name,
                           '_Twidth.tsv'))
}


print('done')
