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

    fastqtoBAM      Trims and maps reads onto there a reference genome
    binning         Calculates mappability and gc content for bins to be used with Kronos CNV
    CNV             Calculates copy number variation 
    10xtoKronos     Converts 10X genomics files in a format that Kronos can use
    diagnostic      Plotting tools to identify appropriate tresholds for Kronos RT
    RT              Calculates scReplication profiles and scRT
    compare         Compares results from multiple experiments

### Authors

Please contact the authors for any further questions:

Stefano Gnan stefano.gnan@curie.fr

Chunlong Chen chunlong.chen@curie.fr


