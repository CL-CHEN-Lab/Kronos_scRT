#  Kronos: A uniform framework for single-cell replication timing analysis

DNA replication is a fundamental process in all living organisms. If all of the origins of replication were activated at the same time in a single mammalian cell it would take about 30 minutes for complete genome replication to occur. However, in nature DNA replication actually takes several hours due to the need to coordinate with other processes such as chromatin 3D organization and transcription. The cell-type specific program that regulates the spatiotemporal progression of DNA replication is referred to as replication timing (RT). In recent years, advances in high-throughput single-cell (sc) sequencing techniques have made it possible to analyze the RT at the single-cell level enabling detection of cell-to-cell variability. This work aims to create a uniform computational framework to investigate scRT using large single-cell whole genome sequencing datasets based on single cell copy number variation (scCNV) detection. This pipeline can be used to analyze datasets from various experiments including classical scWGA (single-cell Whole Genome Amplification), 10x Genomics scCNV solution and scHiC (single cell High-throughput Chromosome conformation capture) from the unsorted cells or cells sorted to enrich for the S phase population. The framework described here allows to increase the number of cells used to analyze scRT by 10 fold (>1000 cells) compared to the current existing analysis. Potentially, this pipeline can also be combined with the analysis of single-cell RNA-seq, CpG methylation and chromatin accessibility to study the relation between expression, chromatin architecture and RT at the single-cell level.

### Case study
To test our software we generated sc-gDNA sequencing data from MCF7 (ER-positive breast cancer) cell line. In order to increase the number of cells in S-phase, cells were sorted before using the 10x genomics platform for single cell copy number variation. In this specific case, the reduction of the G1/G2 phase impaired the automatic identification of the S-phase. However, the user can set a threshold in order to manually select G1/G2 and S phase cells.

### Cell cycle staging
Kronos and Cell Ranger (10x Genomics) can calculate cell ploidy and variability inside a cell (fig 1A). These parameters can be used to identify the cell cycle stage of each cell. Both programs rely on the same function to calculate cell ploidy. This formula, as shown by the simulated data in figure 1B, introduces a bias assigning the same copy number to G1 and G2 cells as well as dividing the S phase into two fractions (fig 1A and 1B). Fitting the simulation into a linear model it is possible to estimate parameters to remove this bias from S phase cells (fig 1B). As a result, after correction S phase cells are distributed accordingly to their mean ploidy (fig 1C).

![](https://github.com/CL-CHEN-Lab/Kronos_scRT/blob/master/img/1.png)

### Reconstructing the Replication Timing Program
Once the copy number has been adjusted the median profile of the G1/G2 population can be used to normalise the profile of each cell in S-phase. Data from each cell are then binarised using a threshold that minimises the euclidean distance between the real data and their binary counterpart. Cells that poorly correlated with the rest of the sample are eliminated (fig 2A, Pearson correlation before and after filtering) and the rest is used to calculate the pseudo bulk RT profile (fig 2B, In the upper part of the plot referenceRT (population RT data) and RT (pseudo bulk RT calculated from scRT data): red=Early, blue=Late; below, Replication Tracks for individual cells, order from early to late from top to bottom (in blue=non replicated and in red=replicated). The pseudo bulk RT and the population RT have a very high correlation (Pearson correlation R=0.9).

![](https://github.com/CL-CHEN-Lab/Kronos_scRT/blob/master/img/2.png)

### Studying the DNA replication program
Kronos offers a series of diagnostic plots. The main program delivers immediately the Twidth value (Dileep et al 2018) that describes the cell to cell variability (fig 3A). Kronos can compare as well multiple samples and identify RT changing regions (fig 3B and 3C).

![](https://github.com/CL-CHEN-Lab/Kronos_scRT/blob/master/img/3.png)


Reference:
Dileep, Vishnu, and David M. Gilbert, ‘Single-Cell Replication Profiling to Measure Stochastic Variation in Mammalian Replication Timing’, Nature Communications, 9 (2018), 427 

### Kronos: pipeline
![](https://github.com/CL-CHEN-Lab/Kronos_scRT/blob/master/img/schema.jpg)

### Run script

First download the repository in the location of your choice, either with git clone git@github.com:CL-CHEN-Lab/Kronos_scRT.git or by clicking on 'Clone or Download' -> 'Download ZIP' and unzip.

Make sure to make the main script Kronos executable :

        chmod +X Kronos

Run the script
    
    ./Kronos <command> [options]
        
    commands:

    fastqtoBAM      Trims and maps reads
    binning         Calculates mappability and gc content for bins to be used with Kronos CNV
    CNV             Calculates copy number variation 
    10xtoKronos     Converts 10X genomics files in a format that Kronos can use
    diagnostic      Plotting tools to identify appropriate tresholds for Kronos RT
    RT              Calculates scReplication profiles and scRT
    compare         Compares results from multiple experiments

-- fastqtoBAM module

    ./Kronos fastqtoBAM [options]
    
    Options:
    --one=CHARACTER                                     Fastq files
    --two=CHARACTER                                     Fastq files (for paired ends)
    -b CHARACTER, --sam_file_basename=CHARACTER         Sam file name
    -i CHARACTER, --index=CHARACTER                     Bowtie 2 index
    -c INTEGER, --cores=INTEGER                         Number of cores to use. [default= 3]
    -o CHARACTER, --output_dir=CHARACTER                Output folder. [default= output/]
    --path_to_trim_galore=CHARACTER                     Path to trim_galore
    --path_to_cutadapt=CHARACTER                        Path to cutadapt
    --path_to_java=CHARACTER                            Path to java
    --path_to_picard=CHARACTER                          Path to picard
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
    -d CHARACTER, --dir_indexed_bam=CHARACTER           If provided parameters will be automatically estimated form the data.
    -h, --help                                          Show this help message and exit
    
-- CNV module

    ./Kronos CNV [options]

    Options:
    -D CHARACTER, --directory=CHARACTER                 Single cell Bamfiles directory
    -B CHARACTER, --bins=CHARACTER                      File with bins produced by Kronos binning
    -X, --keep_X                                        Keep X chromosomes. If active n of X chromosomes has to be provided, by default it will assume 1 X chromosomes if Kee_Y is not selected or 2 if it is not. Correct number can be passed using --number_of_X. [default= FALSE]
    --number_of_X=INTEGER                               Number of X chromosomes 
    -Y, --keep_Y                                        Keep Y chromosome [default= FALSE]
    -m DOUBLE, --min_n_reads=DOUBLE                     Min n of reads to keep a cell in the analysis [default= 2e+05]
    -c INTEGER, --cores=INTEGER                         Number of cores to use. [default= 3]
    -o CHARACTER, --output_dir=CHARACTER                Output folder. [default= output/]
    -e CHARACTER, --ExpName=CHARACTER                   Experiment name. [default= Exp]
    -h, --help                                          Show this help message and exit

-- 10xtoKronos module

    ./Kronos 10xtoKronos [options]

    Options:
    -F CHARACTER, --file=CHARACTER                      Per cell stat file , if multiple files are provided they have to be separated by a comma
    -T CHARACTER, --tracks=CHARACTER                    Tracks file,  if multiple files are provided they have to be separated by a comma
    -o CHARACTER, --out=CHARACTER                       Output directory [default= output]
    -h, --help                                          Show this help message and exit

-- diagnostic module

    ./Kronos diagnostic [options]

    Options:
    -f CHARACTER, --file=CHARACTER                      Dataset file name
    -o CHARACTER, --out=CHARACTER                       Output directory [default= ./output]
    -b CHARACTER, --base_name=CHARACTER                 Base name for files names [default= exp]
    -S DOUBLE, --threshold_Sphase=DOUBLE                Threshold to identify S-phase cells
    -G DOUBLE, --threshold_G1G2phase=DOUBLE             Threshold to identify G1-phase cells. -S has to be selected and has to be bigger than -G
    -h, --help                                          Show this help message and exit

-- RT module

    ./Kronos RT [options]

    Options:
    -F CHARACTER, --file=CHARACTER                      Per cell stat file , if multiple files are provided they have to be separated by a comma
    -T CHARACTER, --tracks=CHARACTER                    Tracks file,  if multiple files are provided they have to be separated by a comma
    -R CHARACTER, --referenceRT=CHARACTER               Reference RT min=Late, max=Early, only one reference is allowed
    -C CHARACTER, --chrSizes=CHARACTER                  Chromosome size file
    -r CHARACTER, --region=CHARACTER                    Region to plot  chr:start-end (multiple regins can be separated by a comma)
    -o CHARACTER, --out=CHARACTER                       Output directory [default= output]
    -b CHARACTER, --base_name=CHARACTER                 Base name for files names [default= exp]
    -g CHARACTER, --groups=CHARACTER                    Grouping names of multiple basenames [default= base_name]
    -S DOUBLE, --threshold_Sphase=DOUBLE                Threshold to identify S-phase cells
    -G DOUBLE, --threshold_G1G2phase=DOUBLE             Threshold to identify G1-phase cells. -S has to be selected and has to be bigger than -G
    -B INTEGER, --binsSize=INTEGER                      RT resolution [default= 500000] 
    -k, --keepXY                                        Keep XY chromosomes in the analysis
    -c INTEGER, --cores=INTEGER                         Numbers of parallel jobs to run [default= 3] 
    -p, --plot                                          If selected prints some randome regins, if -r is selected those regins are use to print RT [default= FALSE] 
    --Var_against_reference                             Variability metrics are calculated usign reference RT in addiction to the calculated one [default= FALSE] 
    --min_correlation=DOUBLE                            Minimum correlation value between one cell and its best correlating cell for this cell to not be discarded [default= 0.25] 
    -h, --help                                          Show this help message and exit

-- compare module

    ./Kronos compare [options]

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

### Requirements
    Programs:
        cutadapt
        trim_galore
        java
        picard
    
    R Packages:
     - Required packages are installed automatically by the modules

### Authors

Please contact the authors for any further questions:

Stefano Gnan stefano.gnan@curie.fr

Chunlong Chen chunlong.chen@curie.fr


