#  Kronos scRT: A uniform framework for single-cell replication timing analysis

DNA replication is a fundamental process in all living organisms. If all of the origins of replication were activated at the same time in a single mammalian cell it would take about 30 minutes for complete genome replication to occur. However, in nature DNA replication actually takes several hours due to the need to coordinate with other processes such as chromatin 3D organization and transcription. The cell-type specific program that regulates the spatiotemporal progression of DNA replication is referred to as replication timing (RT). In recent years, advances in high-throughput single-cell (sc) sequencing techniques have made it possible to analyze the RT at the single-cell level enabling detection of cell-to-cell variability. This work aims to create a uniform computational framework to investigate scRT using large single-cell whole genome sequencing datasets based on single cell copy number variation (scCNV) detection. This pipeline can be used to analyze datasets from various experiments including classical scWGA (single-cell Whole Genome Amplification), 10x Genomics scCNV solution and scHiC (single cell High-throughput Chromosome conformation capture) from the unsorted cells or cells sorted to enrich for the S phase population. The framework described here allows to increase the number of cells used to analyze scRT by 10 fold (>1000 cells) compared to the current existing analysis(1,2). Potentially, this pipeline can also be combined with the analysis of single-cell RNA-seq, CpG methylation and chromatin accessibility to study the relation between expression, chromatin architecture and RT at the single-cell level.

### Case study
To test our software we generated sc-gDNA sequencing data from MCF7 (ER-positive breast cancer) cell line. In order to increase the number of cells in S-phase, cells were sorted before using the 10x genomics platform for single cell copy number variation. In this specific case, the reduction of the G1/G2 phase impaired the automatic identification of the S-phase. However, the user can set a threshold in order to manually select G1/G2 and S phase cells.

### Cell cycle staging
Kronos and Cell Ranger (10x Genomics) can calculate cell ploidy and variability inside a cell (fig 1A). These parameters can be used to identify the cell cycle stage of each cell. Both programs rely on the same function to calculate cell ploidy. This formula, introduces two biases. Firstly, it is not possible to distinguish between G1 and G2 cells that co-occupy the same area (Fig. 1 A red population). Secondly, the S phase is split in two: with the first part progresses normally, while the second part is approaching the G1/G2 population from the left side of the plot (Fig. 1 A green population). Kronos diagnostic calculates two parameters to correct S-phase populations. Preferentially, the program tries to reunite the S phase in a monomodal distribution in which the ploidy variability is maximized, when this is not possible, parameters are chosen in order to create a bimodal distribution with a minimized ploidy variability (Fig. 2 B). The User can as well manually set these parameters.

![](https://github.com/CL-CHEN-Lab/Kronos_scRT/blob/master/img/1.png)

### Reconstructing the Replication Timing Program
Once the copy number has been adjusted the median profile of the G1/G2 population can be used to normalise the profile of each cell in S-phase. Data from each cell are then binarised using a threshold that minimises the euclidean distance between the real data and their binary counterpart (an example in fig 2A). Cells that poorly correlated with the rest of the sample are eliminated (fig 2B, Pearson correlation before and after filtering) and the rest is used to calculate the pseudo bulk RT profile (fig 2C, In the upper part of the plot referenceRT (population RT data) and RT (pseudo bulk RT calculated from scRT data): red=Early, blue=Late; below, Replication Tracks for individual cells, order from early to late from top to bottom (in blue=non replicated and in red=replicated). The pseudo bulk RT and the population RT have a very high correlation (Pearson correlation R=0.93).

![](https://github.com/CL-CHEN-Lab/Kronos_scRT/blob/master/img/2.png)

### Studying the DNA replication program
Kronos offers a series of diagnostic plots. The main program delivers immediately the Twidth value (1) that describes the cell to cell variability (fig 3A). Kronos can compare as well multiple samples and identify RT changing regions (fig 3B and 3C).

![](https://github.com/CL-CHEN-Lab/Kronos_scRT/blob/master/img/3.png)


### References:
1 - Dileep, V. & Gilbert, D. M. Single-cell replication profiling to measure stochastic variation in mammalian replication timing. Nat. Commun. 9, 427 (2018).

2 - Takahashi, S. et al. Genome-wide stability of the DNA replication program in single mammalian cells. Nat. Genet. 51, 529–540 (2019).

### Kronos: pipeline
![](https://github.com/CL-CHEN-Lab/Kronos_scRT/blob/master/img/schema.jpg)

### Run script

First download the repository in the location of your choice, either with git clone https://github.com/CL-CHEN-Lab/Kronos_scRT.git or by clicking on 'Clone or Download' -> 'Download ZIP' and unzip.

Make sure to make the main script Kronos executable :

     chmod +X Kronos

Run the script
    
    ./Kronos <command> [options]
        
    commands:
    
    InstRpacks      Install all the required R packages
    fastqtoBAM      Trims and maps reads
    suggests        Suggests parameters for binning, cell filtering and RT resolution starting from fastqtoBAM report or manual input
    binning         Calculates mappability and gc content for bins to be used with Kronos CNV
    CNV             Calculates copy number variation 
    10xtoKronos     Converts 10X genomics files in a format that Kronos can use
    WhoIsWho        Manually assign cell cycle stage
    diagnostic      Plotting tools to identify appropriate thresholds for Kronos RT
    RT              Calculates scReplication profiles and scRT
    compare RT      Compares RT results from multiple experiments
    compare TW      Compares variability from multiple experiments and/or over multiple regions
    population RT   Calcualtes population RT starting from single cell BAM files and Kronos diagnostic outputs
--InstRpacks

    ./Kronos InstRpacks

-- fastqtoBAM module

    ./Kronos fastqtoBAM [options]
    
    Options:
    -l CHARACTER, --fastq_list=CHARACTER                A table formatted in the following way: sam_file_basename\tFastq_1\tFastq_2(optional for PE sequencing). Compressed files are not allowed. Alternative to --one/--two/-b
    --one=CHARACTER                                     Fastq files, for multiple files they have to be separated by a comma. Compressed files are not allowed.  Alternative to -l
    --two=CHARACTER                                     Fastq files (for paired ends), for multiple files they have to be separated by a comma. Compressed files are not allowed.  Alternative to -l
    -b CHARACTER, --sam_file_basename=CHARACTER         Sam file name, for multiple files they have to be separated by a comma.  Alternative to -l
    -i CHARACTER, --index=CHARACTER                     Bowtie 2 index
    -c INTEGER, --cores=INTEGER                         Number of cores to use. [default= 3]
    -o CHARACTER, --output_dir=CHARACTER                Output folder. [default= output/]
    --path_to_trim_galore=CHARACTER                     Path to trim_galore
    --path_to_cutadapt=CHARACTER                        Path to cutadapt
    --path_to_java=CHARACTER                            Path to java
    --path_to_picard=CHARACTER                          Path to picard
    -h, --help                                          Show this help message and exit

-- suggests module

    ./Kronos suggests [options]


    Options:
	-R CHARACTER, --single_cell_reads_info=CHARACTER    Single cell reads info file produced by Kronos fastqtoBAM or mean reads number followed by PE or SE ( for example: 100PE ). 
    -G NUMERIC, --genome_size=NUMERIC                   Genome size in MB [default= 2900]
    -P NUMERIC, --expected_ploidy=NUMERIC               Expected ploidy [default= 2]
    -h, --help                                          Show this help message and exit

-- binning module

    ./Kronos binning [options]

    Options:
    -R CHARACTER, --RefGenome=CHARACTER                 Fasta file of genome of interst
    -c INTEGER, --cores=INTEGER                         Number of cores to use. [default= 3]
    -s INTEGER, --reads_size=INTEGER                    Lengh of the simulated reads. [default= 40 bp]
    -o CHARACTER, --output_dir=CHARACTER                Output folder. [default= output/]
    -i CHARACTER, --index=CHARACTER                     Bowtie 2 index
    --paired_ends                                       Generates paired ends reads [default: FALSE]
    --insert_size=INTEGER                               Insert size if paired end option is used. [default: 200]
    --bin_size=INTEGER                                  Bins size. [default= 20000 bp]
    -d CHARACTER, --dir_indexed_bam=CHARACTER           If provided parameters will be automatically estimated from the data.
    -u DOUBLE, --upper_mappability_th=DOUBLE            Maximum mappability for a bin to be considered in the analisys  [default= 1.5]
    -l DOUBLE, --lower_mappability_th=DOUBLE            Minimum mappability for a bin to be considered in the analisys  [default= 0.8]
    -B CHARACTER, --black_list=CHARACTER                Regions to ignore
    -x CHARACTER, --coverage=CHARACTER                  Coverage for simulated genome. [default= 1x]
    -h, --help                                          Show this help message and exit
    
-- CNV module

    ./Kronos CNV [options]

    Options:
    -D CHARACTER, --directory=CHARACTER                 Single cell Bamfiles directory
    -B CHARACTER, --bins=CHARACTER                      File with bins produced by Kronos binning
    -X, --keep_X                                        Keep X chromosomes [default= FALSE]
    -Y, --keep_Y                                        Keep Y chromosome [default= FALSE]
    -n DOUBLE, --min_n_reads=DOUBLE                     Min n of reads to keep a cell in the analysis [default= 2e+05]
    -c INTEGER, --cores=INTEGER                         Number of cores to use. [default= 3]
    -o CHARACTER, --output_dir=CHARACTER                Output folder. [default= output/]
    -e CHARACTER, --ExpName=CHARACTER                   Experiment name. [default= Exp]
    -p NUMERIC, --ploidy=NUMERIC                        User extimated ploidy (optional)
    -m NUMERIC, --min_CNV_accepted=NUMERIC              Min mean CNV accepted as result. [default= 0]
    -M NUMERIC, --max_CNV_accepted=NUMERIC              Max mean CNV accepted as result. [default= 8]
    -h, --help                                          Show this help message and exit

-- 10xtoKronos module

    ./Kronos 10xtoKronos [options]

    Options:
    -F CHARACTER, --file=CHARACTER                      Per cell stat file , if multiple files are provided they have to be separated by a comma
    -T CHARACTER, --tracks=CHARACTER                    Tracks file,  if multiple files are provided they have to be separated by a comma
    -o CHARACTER, --out=CHARACTER                       Output directory [default= output]
    -h, --help                                          Show this help message and exit

-- WhoIsWho module

    ./Kronos WhoIsWho [options]

    Options:
    -F CHARACTER, --file=CHARACTER                      per cell stat file path
    -W CHARACTER, --whoSwho=CHARACTER                   who's who file path
    -o CHARACTER, --out=CHARACTER                       Output directory [default= ./output]
    -h, --help                                          Show this help message and exit

-- diagnostic module

    ./Kronos diagnostic [options]

    Options:
    -F CHARACTER, --file=CHARACTER                      Dataset file name
    -o CHARACTER, --out=CHARACTER                       Output directory [default= ./output]
    -b CHARACTER, --base_name=CHARACTER                 Base name for files names [default= exp]
    -S DOUBLE, --threshold_Sphase=DOUBLE                Threshold to identify S-phase cells
    -G DOUBLE, --threshold_G1G2phase=DOUBLE             Threshold to identify G1-phase cells. -S has to be selected and has to be bigger than -G
    -f DOUBLE, --Sphase_first_part=DOUBLE               Correction parameter for the first part of the S-phase [0.95,1]
    -s DOUBLE, --Sphase_second_part=DOUBLE              Correction parameter for the second part of the S-phase [0.5,0.55]
    -c INTEGER, --cores=INTEGER                         Numbers of parallel jobs to run [default= 3] 
    -m DOUBLE, --min_n_reads=DOUBLE                     Min n of reads per million per aploid genome to keep a cell in the analysis [default= 112]
    -h, --help                                          Show this help message and exit

-- RT module

    ./Kronos RT [options]

    Options:
    -K CHARACTER, --Kronos_conf_file=CHARACTER          Kronos setting file. If provided -F,-T,-S,-b and -g are ignored. Tab file containing: Per cell stat file /t tracks file /t settings file /t basename (optional) /t group (optional)
    -F CHARACTER, --file=CHARACTER                      Per cell stat file , if multiple files are provided they have to be separated by a comma
    -T CHARACTER, --tracks=CHARACTER                    Tracks file,  if multiple files are provided they have to be separated by a comma
    -R CHARACTER, --referenceRT=CHARACTER               Reference RT min=Late, max=Early, only one reference is allowed
    --ref_name=CHARACTER                                Name for the reference track [default= Reference]
    -C CHARACTER, --chrSizes=CHARACTER                  Chromosome size file
    -r CHARACTER, --region=CHARACTER                    Region to plot  chr:start-end (multiple regins can be separated by a comma) or a bed file can be provided
    -o CHARACTER, --out=CHARACTER                       Output directory [default= output]
    -b CHARACTER, --base_name=CHARACTER                 Base name for files names [default= exp]
    -f CHARACTER, --output_file_base_name=CHARACTER     Base name for the output file [default= out]
    -g CHARACTER, --groups=CHARACTER                    Grouping names of multiple basenames [default= base_name]
    -S CHARACTER, --settings_file=CHARACTER             File generated by Kronos diagnostic
    -B INTEGER, --binsSize=INTEGER                      RT resolution [default= 500000] 
    -k, --keepXY                                        keep XY chromosomes in the analysis
    -c INTEGER, --cores=INTEGER                         Numbers of parallel jobs to run [default= 3] 
    -p, --plot                                          If selected prints some randome regins, if -r is selected those regins are use to print RT [default= FALSE] 
    --Var_against_reference                             Variability metrics are calculated usign reference RT in addiction to the calculated one [default= FALSE] 
    --min_correlation=DOUBLE                            Minimum correlation value between one cell and its best correlating cell for this cell to not be discarded [default= 0.25] 
    -h, --help                                          Show this help message and exit

-- compare  RT module

    ./Kronos compare RT [options]

    Options:
    -S CHARACTER, --S50s=CHARACTER                      RT files with same binning
    -R CHARACTER, --referenceRT=CHARACTER               Reference RT min=Late, max=Early, only one reference is allowed
    -o CHARACTER, --out=CHARACTER                       Output directory [default= output]
    -k, --keepXY                                        keeps XY chr in the analysis
    --Reference=CHARACTER                               Base name to use as a reference, if not provided the first basename in the S50 file will be used or , if provided , the reference RT even if this option is selected
    -D DOUBLE, --deltaRT_threshold=DOUBLE               DeltaRT threshold to define changes
    -n INTEGER, --n_regions=INTEGER                     number of regions to plot
    -r CHARACTER, --region=CHARACTER                    Region to plot  chr:start-end (multiple regins can be separated by a comma)
    -f CHARACTER, --basename_filter=CHARACTER           Filter out unwanted samples for RT files
    -h, --help                                          Show this help message and exit

-- compare  TW  module

    ./Kronos compare TW [options]

    Options:
    -C CHARACTER, --CNV=CHARACTER                       scCNV file produced by Kronos RT, if multiple files are provided they have to be separated by a comma
    -T CHARACTER, --RT=CHARACTER                        RT file produced by Kronos RT, if multiple files are provided they have to be separated by a comma
    -R CHARACTER, --regions=CHARACTER                   Genome annotation. chr<TAB>start<TAB>end<TAB>annotation. No header.
    -b, --both_annotations                              Plot Twidth divided by both annotations
    -r CHARACTER, --regions2=CHARACTER                  Second genome annotation. chr<TAB>start<TAB>end<TAB>annotation. No header. If option b is activated it substiutes the RT division.
    -o CHARACTER, --out=CHARACTER                       Output directory [default= output]
    -f CHARACTER, --output_file_base_name=CHARACTER     Base name for the output file [default= out]
    -h, --help                                          Show this help message and exit
                                      Show this help message and exit

-- Kronos population RT module

    ./Kronos population RT [options]
    
    Options:
    -F CHARACTER, --file=CHARACTER                 Per cell stat file, to merge multiple runs separate directories with a comma
    -S CHARACTER, --settings_file=CHARACTER        File generated by Kronos diagnostic, to merge multiple runs separate directories with a comma
    -D CHARACTER, --directory=CHARACTER            Single cell Bamfiles directory, to merge multiple runs separate directories with a comma
    -o CHARACTER, --out=CHARACTER                  Output directory [default= output]
    -b CHARACTER, --base_name=CHARACTER            Base name for files names [default= exp]
    -X, --keep_X                                   Keep X chromosomes. [default= FALSE]
    -Y, --keep_Y                                   Keep Y chromosome. [default= FALSE]
    -c INTEGER, --cores=INTEGER                    Numbers of parallel jobs to run [default= 3] 
    -R INTEGER, --bin_size=INTEGER                 Bins size in bp,multiple bin size can be provided separated by a comma. [default= 50Kb]
    -C CHARACTER, --chrSizes=CHARACTER             Chromosome size file
    -B CHARACTER, --black_list=CHARACTER           Regions to ignore
    -h, --help                                     Show this help message and exit

-- Kronos scPlots module

    ./Kronos_scPlots.R [options]

    Options:
    -L CHARACTER, --List=CHARACTER                  A Tab separated file containing in each colum scRT_Tracks and scCNV files paths. Alternative to -R,-T and -C options.
    -R CHARACTER, --scRT_Tracks=CHARACTER           *calculated_replication_timing* file(s) created by Kronos RT. If multiple files are provided they have to be separated by a comma.  Alternative to -L option.
    -C CHARACTER, --scCNV=CHARACTER                 *single_cells_CNV* file(s) created by Kronos RT. If multiple files are provided they have to be separated by a comma.  Alternative to -L option..
    -s CHARACTER, --order=CHARACTER                 basenames separated by a comma in the desired order for plotting.
    --CNV_values=CHARACTER                          What type of date to plot for the sigle cell traks: ('B'=Binarized, 'CNV'=Copy number variation, 'log2'=log2(CNV_Cell/CNV_mean_G1/G2_cells) or 'all'= one file per option) [default= B]
    -r CHARACTER, --region=CHARACTER                Region to plot  chr:start-end (multiple regins can be separated by a comma) or provided as a bed file
    -o CHARACTER, --out=CHARACTER                   Output directory [default= output]
    -f CHARACTER, --output_file_base_name=CHARACTER Base name for the output file [default= out]
    -h, --help                                      Show this help message and exit

### Requirements
    Programs:
     - trim_galore
     - picard
    
    R Packages:
     - Biostrings
     - Cairo
     - DNAcopy
     - doSNOW
     - foreach
     - GenomicRanges
     - gridExtra
     - gplots
     - LaplacesDemon
     - MASS
     - matrixStats
     - optparse
     - Rbowtie2
     - RColorBrewer
     - Rsamtools
     - scales
     - tidyverse

### Session Info

    R version 3.5.2 (2018-12-20)
    Platform: x86_64-apple-darwin15.6.0 (64-bit)
    Running under: macOS Mojave 10.14.5

    Matrix products: default
    BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
    LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    attached base packages:
    [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] MASS_7.3-51.4        DNAcopy_1.56.0       scales_1.0.0         Cairo_1.5-10         ggpubr_0.2.1         magrittr_1.5        
    [7] LaplacesDemon_16.1.1 Rbowtie2_1.4.0       RColorBrewer_1.1-2   matrixStats_0.54.0   gplots_3.0.1.1       doSNOW_1.0.16       
    [13] snow_0.4-3           iterators_1.0.10     foreach_1.4.4        forcats_0.4.0        stringr_1.4.0        dplyr_0.8.3         
    [19] purrr_0.3.2          readr_1.3.1          tidyr_0.8.3          tibble_2.1.3         ggplot2_3.2.0        tidyverse_1.2.1     
    [25] Rsamtools_1.34.1     Biostrings_2.50.2    XVector_0.22.0       GenomicRanges_1.34.0 GenomeInfoDb_1.18.2  IRanges_2.16.0      
    [31] S4Vectors_0.20.1     BiocGenerics_0.28.0  optparse_1.6.2      

    loaded via a namespace (and not attached):
    [1] httr_1.4.0             jsonlite_1.6           modelr_0.1.4           gtools_3.8.1           assertthat_0.2.1       GenomeInfoDbData_1.2.0
    [7] cellranger_1.1.0       yaml_2.2.0             pillar_1.4.2           backports_1.1.4        lattice_0.20-38        glue_1.3.1            
    [13] ggsignif_0.5.0         rvest_0.3.4            colorspace_1.4-1       pkgconfig_2.0.2        broom_0.5.2            haven_2.1.1           
    [19] zlibbioc_1.28.0        gdata_2.18.0           getopt_1.20.3          BiocParallel_1.16.6    generics_0.0.2         withr_2.1.2           
    [25] lazyeval_0.2.2         cli_1.1.0              crayon_1.3.4           readxl_1.3.1           nlme_3.1-140           xml2_1.2.0            
    [31] tools_3.5.2            hms_0.5.0              munsell_0.5.0          compiler_3.5.2         caTools_1.17.1.2       rlang_0.4.0           
    [37] grid_3.5.2             RCurl_1.95-4.12        rstudioapi_0.10        bitops_1.0-6           gtable_0.3.0           codetools_0.2-16      
    [43] R6_2.4.0               lubridate_1.7.4        zeallot_0.1.0          KernSmooth_2.23-15     stringi_1.4.3          Rcpp_1.0.1            
    [49] vctrs_0.2.0            tidyselect_0.2.5   

### Authors

Please contact the authors for any further questions:

Stefano Gnan stefano.gnan@curie.fr

Chunlong Chen chunlong.chen@curie.fr


