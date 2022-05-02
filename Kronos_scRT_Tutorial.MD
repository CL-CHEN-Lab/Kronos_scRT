# Kronos scRT Tutorial

This Tutorial will guide you through the main steps of [Kronos scRT](https://github.com/CL-CHEN-Lab/Kronos_scRT), a software for single-cell Replication Timing (scRT) analysis.

## Table of content
  * [Installation & Preparation](#Installation-and-preparation)
  * [MCF-7 Tutorial](#Tutorial-on-MCF-7-cells)
    * [Binning](#Bin-the-genome)
    * [Copy-Number Variation](#Copy-number-variation)
    * [Diagnostic](#Correction-parameters-and-unveiling-cell-phases-with-the-Diagnostic-module)
    * [Replication Timing](#Replication-Timing)
    * [Dimensionality Reduction](#Dimensionality-reduction)
  * [Runtime](#Runtime)
  * [Citation](#Citation)


## Installation and preparation
You can download Kronos scRT to your desired directory by using [git](https://git-scm.com/downloads).
```bash
git clone https://github.com/CL-CHEN-Lab/Kronos_scRT
```
Make Kronos scRT excecutable and export its path (this makes it easier to call Kronos scRT).
```bash
chmod +X ./Kronos_scRT/Kronos
PATH=$PATH:$(pwd)/Kronos_scRT
```
Explore the help option listing the various modules.
```bash
Kronos -h
```
Install the required R packages.
```bash
Kronos InstRpacks
```

## Tutorial on MCF-7 cells

### Getting started
#### Prerequisites
For this tutorial you will need to use hg38 as a reference genome. If you don't already have a hg38 template, you can download it from [Illumina's iGenome website](https://support.illumina.com/sequencing/sequencing_software/igenome.html) which includes its bowtie2 index. Otherwise, you can download or create the [bowtie2 index](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) yourself.
You will also need [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/), [Cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html), [Java](https://www.oracle.com/java/technologies/downloads/) and [Picard](https://broadinstitute.github.io/picard/). 
We strongly recommend the use of a [blacklist](https://doi.org/10.1038/s41598-019-45839-z) to mask problematic genomic regions (centromeres, telomeres, microsatellites, etc...).
#### Data for the tutorial
We will need to download the related data which are comprised of whole-genome sequencing of 36 [MCF-7](https://en.wikipedia.org/wiki/MCF-7) cells. The reads have already been aligned to the reference genome (hg38 in this case) resulting in 1 BAM file per cell.
We can now download the compressed data with [curl](https://curl.se/download.html), extract with [tar](https://linux.die.net/man/1/tar) and move to the downloaded directory.
```bash
curl -O https://xfer.curie.fr/get/nil/WujPODodyAG/Kronos_scRT_Tutorial.tar.gz
tar -xf Kronos_scRT_Tutorial.tar.gz
cd Kronos_scRT_Tutorial
```
We will also need a [blacklist]() and the [MCF-7 Reference RT file](https://www2.replicationdomain.com/database.php) (to compare the resulting pseudo-bulk RT issued from our single cells).
```bash
## Blacklist
curl -O https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
## Reference RT
# For MacOS (save to current directory)
open 'https://www2.replicationdomain.com/download_data.php?dataset=Int67688290&build=hg38'
# For Linux (save to current directory)
xdg-open 'https://www2.replicationdomain.com/download_data.php?dataset=Int67688290&build=hg38'
# Remove header
cat RT_MCF-7_Adenocarcinoma_Int67688290_hg38.bedgraph | tail -n +12 > RT_MCF-7_Adenocarcinoma_Int67688290_hg38_no_header.bedgraph
```


#### Declaring paths
The first step is to declare all the paths/values we want to use. You will have to update the following paths to their actual location on your side:
  * Programs:
    * TrimGalore
    * Cutadapt
    * Java
    * Picard
  * Reference Genome Files:
    * Fasta
    * Bowtie 2 index
    * Blacklist
  * Directory Containing the BAM files which will be used in the analysis
  * Reference RT file
  * Number of processing units (cores) to simultaneously processs the data

<i>Tip: If you already have these tools installed, you can use the "which" command to find their location (e.g.</i> ```which cutadapt```<i>)</i>.

```bash
TrimGalore=$(which trim_galore)
Cutadapt=$(which cutadapt)
Java=$(which java)
Picard=Path_to_Picard/picard.jar

hg38Ref=../Homo_sapiens_UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
hg38BT2=../Homo_sapiens_UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome
hg38Blist=./hg38-blacklist.v2.bed

BAMDir=./BAM_files/MCF7_BAM
RefRT=./RT_MCF-7_Adenocarcinoma_Int67688290_hg38_no_header.bedgraph

cores=4
```

### Bin the genome
Now we can begin by cutting up our reference genome into 20 kilobase (kb) bins by providing:
 * Fasta file of the reference genome (-R).
 * Indexed genome from Bowtie 2 (-i).
 * Blacklist bed file (-B).
 * Bin size (--bin_size).
 * Since we have paired-end reads we need to flag this (--paired_ends).
 * Output directory which will be created automatically (-o).
 * Number of cores (-c).
```bash
Kronos binning -R $hg38Ref -i $hg38BT2 -B $hg38Blist --bin_size 20000 --paired_ends -o ./Binning -c $cores
```
| Binning Output File | Further Usage |
| ---| --- |
| TSV file with mappability and GC content of each 20kb bin. | CNV -B |

### Copy-number variation
Once the genome has been cut up in fixed-sized bins, we can proceed with counting the number of reads in each bin of each cell to estimate the copy-number by providing:
 * Directory containing the BAM files (-D).
 * Binning file previously calculated (-B).
 * Output directory (-o).
 * Name of the experiment (-e).
 * Number of cores (-c).
```bash
Kronos CNV -D $BAMDir -B ./Binning/*.tsv -o ./CNV -e MCF7 -c $cores
```
| CNV Output File | Further Usage |
| ---| --- |
| BED file containing the copy-number for each 20kb bin containing reads in each cell. | RT -T |
| CSV file containing statistics for each cell. | Diagnostic -F; WhoIsWho -F; RT -F |

### Correction parameters and unveiling cell phases with the Diagnostic module
The next step consists of categorising our cells into 2 major categories: Replicating (S phase) and Non-replicating (G1/G2 Phase). To do so, we have to find the bin-to-bin variability and cell ploidy cut-off values. There are several ways to do so:

#### 1. Bin-to-bin variability cut-offs for S/G phase
Due to the fact that the calculated copy-numbers are relative, cells belonging to S phase are split into early and late and will display in the upper-left and upper-right parts of the diagnostic plot, in relation to the G1/G2 population. We will call diagnostic a first time to visualise our cells on this graph by providing:
 * CSV cell statistics file from CNV (-F).
 * Output directory (-o).
 * Basename of the output files (-b).
 * Number of cores (-c).
 * Reads per megabase per haploid genome threshold (-m).
```bash
Kronos diagnostic -F ./CNV/MCF7_per_Cell_summary_metrics.csv -o ./Diagnostic -b MCF7_inspection -c $cores -m 117
```
<img src="https://user-images.githubusercontent.com/49646747/159041873-dd661953-d724-43cf-9f45-5f25657b7185.png" width="500"/>

Due to the limited number of cells here, automatic phase detection can not be achieved. Thus, Kronos scRT cannot distinguish replicating from non-replicating cells and all cells are color-coded as G1/G2-phase cells. By visually inspecting this plot, you can set:
  - Bin-to-bin copy-number variability cut-off parameters that will separate S cells (-S) from the G1/G2 (-G) population (y-axis).
  - Correction parameters for the first (-f) and second (-s) part of S phase which will essentially shift the mean ploidy of those cells accordingly.

Here we will only define the S phase variability threshold and we will let Kronos select the best ploidy correction parameters for the S phase cells (By not providing the -f and -s arguments, Kronos scRT will do this automatically). We can see that the cells start to branch out to both the left (second half of S phase) and the right (first part of S phase) somewhere between 0.8 and 0.9. Thus, we can run Diagnostic again and select 0.85 as the S phase variability cut-off by setting the required arguments:
 * Variability cut-off parameter to define S phase cells (-S).
 * Request correction of S phase and the final output file with all the selected parameters (-C).
```bash
Kronos diagnostic -F ./CNV/MCF7_per_Cell_summary_metrics.csv -o ./Diagnostic -b MCF7_manual -c $cores -m 117 -S 0.85 -C
```
<img src="https://user-images.githubusercontent.com/49646747/159042349-2b374f5d-3d65-4f4a-a2d1-2be6739c1a6d.png" width="500"/>

#### 2. Providing cell phases
 1. If the metadata is available (e.g. from FACS sorting), you can use this information to attribute a phase to each cell (S or G) and then call the diagnostic module. To demonstrate this, we have provided a file containing the cell phases that can be used with the WhoIsWho module with the following arguments:
 * CSV cell statistics file from CNV (-F).
 * Cell name/type correspondence TSV file (-W).
 * Output directory (-o).

This time, instead of providing the the diagnostic module with the cell statistics file from CNV (option -F), we will provide the new file from WhoIsWho with the forced cell phase information. We will let Kronos scRT calculate the ploidy correction values.
```bash
Kronos WhoIsWho -F ./CNV/MCF7_per_Cell_summary_metrics.csv -W ./MCF-7_Cell_Phase_Correspondence.tsv -o ./WhoIsWho
Kronos diagnostic -F ./WhoIsWho/phased_MCF7_per_Cell_summary_metrics.csv -o ./Diagnostic -b MCF7_who -c $cores -m 117
```

#### 3. Complete automatic parameter detection
When there are enough cells (numbers may depend by cell type and experiment) you can let Kronos scRT automatically calculate all the cut-off and correction parameters without providing any extra arguments (-S -G -f -s). In this tutorial, we do not have enough cells to demonstrate this but it could have been achieved by running the following command.
```bash
#Just for reference. Do not run.
Kronos diagnostic -F ./CNV/MCF7_per_Cell_summary_metrics.csv -o ./Diagnostic -b MCF7_auto -c $cores -m 117 -C
```

| WhoIsWho Output File | Further Usage |
| ---| --- |
| CSV file containing corrected cell statistics from the cell phase data. | Diagnostic -F |

| Diagnostic Output File | Further Usage |
| ---| --- |
| TXT file with the parameters used for the copy-number adjustments and variability cut-offs. | RT -S |

### Replication Timing
Before calculating RT, we need a file with the chromosome lengths of the autosome.
```bash
cut -f1,2 ../Homo_sapiens_UCSC_hg38/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa.fai > ./hg38_autosomal_sizes.tsv
```

Now we can get to the part where we calculate Replication Timing by providing the RT module with:
 * CSV cell statistics file from CNV (-F).
 * Copy-number BED file from CNV (-T).
 * Ploidy correction parameter and variability threshold file from diagnostic (-S).
 * Chromosome size file (-C).
 * Bulk RT reference file (-R).
 * RT reference name (--ref_name).
 * Bin size for CNV and RT calculation (-B).
 * Number of categories (e.g. Early/Mid/Late S phase) for Twidth calculation (-N).
 * Output directory (-o).
 * Flag for plot creation (-p).
 * Genomic coordinates of region(s) to plot (-r).
 * Basename of the output files (-b).
 * Supplementary basename (useful when combining replicates) of the output files (-f).
 * Number of cores (-c).

```bash
Kronos RT -F ./CNV/MCF7_per_Cell_summary_metrics.csv -T ./CNV/MCF7_cnv_calls.bed -C ./hg38_autosomal_sizes.tsv -o ./RT -b MCF7_Tutorial -f 250kb -S ./Diagnostic/MCF7_manual_settings.txt -c $cores -p -r chr2:30Mb-60Mb -R $RefRT -B 250000 --ref_name MCF7_Reference -N 3
```
| RT Output File | Further Usage |
| ---| --- |
| TSV file with the T-width values for each category (All/Early/Mid/Late). | |
| TSV file with the calculated pseudo-bulk RT values from the single cells. | |
| TSV file with the reference RT values in 250kb bins. | |
| TSV file containing scRT variability per region. | |
| TSV file containing single cell copy-numbers and binarized RT for each cell in 250kb bins. | DRed |
| TXT file listing the index for each cell name (e.g. in the plots) after filtering out cells. | |
| TXT file listing the index for each cell name (e.g. in the plots) before filtering out cells. | |

<img src="https://user-images.githubusercontent.com/49646747/159042526-eb6c2e33-8926-4beb-811a-e4dea3cae67e.png" width="400"/>
<img src="https://user-images.githubusercontent.com/49646747/159042695-20755c2b-887b-42d3-af0e-a95ad2237298.png" width="400"/>
<img src="https://user-images.githubusercontent.com/49646747/159042675-f95fb444-07d0-486e-b3be-5e8f303150e7.png" width="400"/>

### Dimensionality reduction
We can now visualize cell-to-cell RT variability using the DRed module and by providing the following parameters:
 * CNV/RT file(s) created by Kronos RT (-C).
 * Output directory (-o).
 * Basename of the output files (-f).
 * Values to plot i.e. binarized RT values and not CNV values for example (--CNV_values)
 * Seed for reproducibility (-s).
 * Skip t-SNE, only plot UMAP (-U).
 * Number of cores (-c).

```bash
Kronos DRed -C ./RT/250kb_single_cells_CNV_250Kb.tsv -o ./DimRed -f MCF7_Tutorial --CNV_values B -s 2022 -U -c $cores
```

## Runtime
We tested the execution time of each part of this tutorial on a macOS 12.3.1 system with 32GB of ram and a 2,3 GHz Quad-Core Intel Core i7 processor.
| Kronos Module | Runtime |
| ---| --- |
| Binning | 45min |
| CNV | 15min |
| Diagnostic | <1min |
| WhoIsWho | <1min |
| RT | 1min |
| DRed | <1min |
| **Total** | 1h2min |


## Citation
If you use Kronos scRT, please cite the following paper:

Gnan, S., Josephides, J.M., Wu, X. et al. Kronos scRT: a uniform framework for single-cell replication timing analysis. Nat Commun 13, 2329 (2022). https://doi.org/10.1038/s41467-022-30043-x

