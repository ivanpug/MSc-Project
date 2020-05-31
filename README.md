# Presentation
The following is an except of my Bioinformatics MSc Final Project.

The programming language used is R / Bioconductor
Other complementary tools and BASH Scripts have been used for data cleaning and dataset preparation.
In order to replicate the anaysis, it is recommended using RStudio (and possibly Anaconda for installing the additional tools) in GNU/Linux Operating System.

The original RData files can be loaded inside an R session (after installing Bioconductor and the required packages).
The files contained inside the "analysis" folder are better described at paragraph **3.4. Bioconductor Analysis**.

RData files size exceed the maximum file size upload limit for GitHub.
The content of **analysis** folder can be downloaded from my Nextcloud server at this [LINK](https://cloud.catalisi.it/s/TJKnztE3nasS825)

The report presented here is a copy of the original "report.Rmd" hosted in this repository. Unfortunately GitHub
doesn't support anymore rendering of Rmd Files ([Announce in 2018](https://yihui.org/en/2018/10/rmd-github/)) so I apologize for any text formatting errors. For a correct rendering of this report it's preferable to clone the repository and opening the original Rmd file inside RStudio.

## Instructions
Before starting, install the latest release of R.
Then get the latest version of Bioconductor by starting R and entering the commands:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.11")
```
Finally, install che required Bioconductor Packages:
```
BiocManager::install(c("scater","scran","limma","ensembldb"))
```
Now you can load the RData files inside **analysis** folder.

---
---
# Laboratory Report - Investigating the role of Nanog during peri-implatation of mouse embryo using scRNA-seq expression analysis
Author: Salvatore Ivan Puglisi (B125042)

Date: 05/2018

---

# 1. Background

The pluripotent state of mouse embryonic stem cells (ESCs) is commonly associated with high levels of certain transcription factors such as Nanog, Sox2 and Oct4, which are important for establishing and maintaining the pluripotent state.  The influence of specific ratios of Oct4 and Nanog in individual ESCs on the pluripotency establishment, has already been observed in-vitro [@munoz_descalzo_silvia_correlations_2012], suggesting the possibility of using the distribution of this ratio as a quantifiers to distinguish between three subpopulations in a ESCs culture: pluripotent, lineage-primed and differentiating cells.
Nanog is a divergent homeodomain protein found in mammalian pluripotent cells [@chambers_functional_2003], its presence is considered and hallmark of pluripotent cells in vivo and in-vitro, and its loss an early marker of differentiation [@chambers_nanog_2007]. Nanog is a core element of the pluripotency gene regulatory network (GRN) [@boyer_core_2005] and previous studies [@ptashne_transcription_2001] suggests a combinatorial control of transcription factors  to taking place in proximity of their target genes. What is still unclear is the magnitude of influence upon the GRN targets delivered by the expression variation of a single factor. 
A culture of mouse ESCs consists in a mixture of cells with different levels of Nanog expression, because individual mouse ESCs express different levels of Nanog. Individual Nanog low cells, if grown in culture, can reproduce the whole population of cells expressing different levels of Nanog. 
ESCs expressing Nanog at lower levels are more likely to differentiate than cells expressing high Nanog levels and so the Nanog low state may be considered a marker of early stage in differentiation.  
ESCs knockout for Nanog gene can be produced and maintained in culture. By inducing the expression of a Nanog transgene is possible to revert back the Nanog expression of the WT ESCs.
The inability of Nanog -/-cells to complete transcription-factor-based reprogramming mirrors the phenotype observed in Nanog null embryos, providing a model to study the unique role of Nanog during the acquisition of pluripotency in early development.

In a previous study [@festuccia_esrrb_2012] the effects of altering the expression of Nanog has been analyzed, generating in-vitro ESCs knockout for Nanog. This ESCs Nanog -/- model could provide a way to study the unique role of Nanog during the acquisition of pluripotency in early development. More than 5,000 genes were confirmed to bind Nanog using ChIP-seq experiments, but only a small subset of 64 genes showed a > 1.5 fold change in increasing or decreasing expression, after reinduction of Nanog activity in the ESCs Nanog -/-. Among this 64 Nanog identified target genes, ESRRB shows the strongest transcriptional induction and has been proven to substitute for Nanog function in pluripotent cells. Furthermore, the findings that ESRRB is a direct target of Nanog, together with the notion that ESRRB can positively regulate Nanog, demonstrate the existence of a positive feedback mechanism [@oliveri_global_2008].
Similar observations suggest that the presence of Nanog is necessary but not sufficient to alter the transcriptional rates of it's target genes and confirm the needing of multiple additional pluripotency transcription factors to bind the same targets. Apparently some pluripotency factors like Oct4 are essential to maintain the pluripotent state (Oct4 depletion leads to differentiation [@hall_oct4_2009])  while fluctuations in Nanog (and consequently ESRRB) confer flexibility to the GRN, tuning the expression of downstream genes and leading to cell fate decisions [@chambers_nanog_2007].
The pre-implantation mouse embryo at day E3.5 consists of an inner cell mass (ICM) not possessing distinct lineage identities. Networks of genes including several known pluripotency markers are observed exclusively at this stage. Implantation occurs at approximately embryonic day E4.5 and marks several key changes in the embryo. The embryonic epiblast is formed, combined with two extra-embryonic layers: the trophectoderm and primitive endoderm (PrE). the pluripotent epiblast dynamically changes post-implantation, developing into a transcriptionally distinct entity primed for differentiation [@mohammed_single-cell_2017].
ESCs are derived from the embryos at stage E3.5 of differentiation and might be expected to show similar expression profile for Nanog. Embryonic cells at stage E4.5 still express Nanog but at a lower average level.
Nanog expression is indeed able to affects the ability of cells to differentiate in-vitro but is unknown how the same genes responding to Nanog observed in-vitro, correlate with Nanog in the embryo at stages E3.5 and E4.5 of differentiation.
In this study the following question will be addressed, using single cell RNA expression data. 

- Is it possible to clarify how the 64 genes responding to Nanog identified using the ESCs Nanog -/- model [@festuccia_esrrb_2012] behave in the embryo?  
- Are these genes unregulated or regulated in the same way observed in the ESCs in-vitro models? 
- Is it possible to measure how the Nanog expression is distributed across individual cells? Does it fluctuates naturally as expected in the real mouse embryo?

Performing a correlation analysis between the Nanog expression profile and the correspondent fluctuation of its target genes in the mouse embryo, will be possible to compare the behaviour of these gene between the ESCs,  pre-implantation (E3.5)  and post-implantation cells (E4.5).

# 2. Purpose

The aim of this study is to investigate if the Nanog target genes proven to bind and respond to Nanog in the ESCs in-vitro model, are correlated with Nanog expression fluctuation in single mESC cells profiled using single cell RNA-seq (scRNA-seq).  We would then like to know how these genes are correlated to Nanog expression  during two stages of the early mouse embryo development: the pre-implantation inner cell mass at E3.5 and the post-implantation epiblast at E4.5.

# 3. Analysis 

An outline of the analytic methods used during my intership is reported below.
More detailed R scripts and higher resolution output graphics can be found into the "analysis" folder, inside the specific subdirectories "all", "e35", "e45".
It is also possible to load the R Workspace Files and run the script inside the R markdown notebook changing some parameters.

## 3.1. Dataset Description

The data analyzed in this project are originated from 

> Mohammed, H., Hernando-Herraez, I., Savino, A., Scialdone, A., Macaulay, I., Mulas, C., Chandra, T., Voet, T., Dean, W., Nichols, J., et al. (2017). Single-Cell Landscape of Transcriptional Heterogeneity and Cell Fate Decisions during Mouse Early Gastrulation. Cell Rep. 20, 1215-1228.

The authors isolated Single cells from C57Bl/6Babr Mus musculus embryos at E3.5, E4.5, E5.5 and E6.5 stages, subjecting them to single-cell RNA-seq protocol, using the platform Illumina HiSeq 2500 (generating 100 bp paired-end reads).

The references for the online repositories are reported below:

BioProject ID: [PRJNA392258](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA392258) 

GEO ID: [GSE100597](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100597)

SRA ID: [SRP110669](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=SRP110669)


Only E3.5 & E4.5 are considered in this study, to explore lineage and regulatory processes involved in early peri-implantation mouse embryos.

Total SRA Experiments: 
721 samples (SRR5763563 - SRR5764283)

E3.5 subset: 
99 samples (SRR5763563 - SRR5763661)

E4.5 subset: 
105 samples (SRR5763662 - SRR5763766)

![](images/dataset_info.png)

*Fig. 1 - Example summary generated by FastQC for one sample.*



## 3.2. Quality control & preprocessing of reads

The analysis has been carried out via SSH on the following server 

> bioinfmsc3.mvm.ed.ac.uk"

The UNIX environment has been prepared, upgrading R to the release 3.5.0 as a requirement for installing the ["scater"] [@mccarthy_scater:_2017] and other [Bioconductor] [@huber_orchestrating_2015] packages [@lun_scran_2018], [@ritchie_limma_2015].

```
source("https://bioconductor.org/biocLite.R")
biocLite(c("scater","scran","limma"))
```

Other command-line tools required for the data preprocessing and the genomic allignment have been installed via [Bioconda](https://bioconda.github.io/) [@gruning_bioconda:_2018].

```
wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh
sh Anaconda3-5.2.0-Linux-x86_64.sh # added Anaconda/bin to PATH
conda install -c bioconda sra-tools fastqc multiqc trim-galore kallisto
```
[SRA Toolkit](https://www.ncbi.nlm.nih.gov/books/NBK158900/) is a collection of tools and libraries for using data in the INSDC Sequence Read Archives.
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) [@andrews_fastqc_2018] and [MultiQC](http://multiqc.info/) [@ewels_multiqc:_2016] are the quality control tools which have been used to assess the overall quality of the reads.
[Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) [@krueger_trim_2018] is a wrapper of cutadapt and Fastqc, used to filter-out the low quality reads and trimming the adapter contamination detected during the quality control.
[kallisto](https://pachterlab.github.io/kallisto/) [@bray_near-optimal_2016] is a pseudoalignment program for rapidly determining the compatibility of high-throughput sequencing reads with targets sequences, without the need for alignment. It has been used for quantifying abundances of transcripts in the Mus musculus transcriptome from the scRNA-Seq data.

The main "project" folder has been created on my home directory, with the following subfolders:

- data -> containing the fastq files downloaded and processed from the SRA repositories
- references -> containing the M. musculus transcriptome and annotation files dowloaded from ENSEMBL repositories
- indices -> containing the index file generated by kallisto, to be used for the preudoalignment
- results -> containing the abundances files generated by kallisto quant
- analysis -> containing the RData files and the output of the Bioconductor analysis

The experiment "SRP110669" has been downloaded from SRA repository inside the "data" folder and the sample files converted to compressed FastQ files by the fastq-dump utility (part of SRA Toolkit).

```
wget -r -N -nd ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP110/SRP110669/

# using the option --split-files forward and reverse reads in every SRA file are separated in 2 Fastq files

fastq-dump --split-files --gzip --sra-id *.sra

# Run FastQC for all samples

$ fastqc *.fastq.gz
```

Looking at the HTML reports generated by FastQC the quality appears to be on average over 25 phred score along the lenght of each read.

![](images/qc_score_pre.png)

*Fig. 2 - MultiQC Sequence Quality Histograms. The mean quality value across each base position in the read for all samples. *


Anyway it is possible to observe that 5' end of the reads in any sample is affected by base content similarity, due to use of random k-mers (which are not so random) or transposases (nextera transposases, causing tagmentation) in library preparation.

![](images/per_base_seq_content.png)

*Fig. 3 - Per Base Sequence Content Plot from FastQC for one sample, showing the relative amount of each base at each position.*


The first run of FASTQC detected adapter contamination by nextera transposase.

![](images/adapter_content_pre.png)
 
*Fig. 4 - FastQC Adapter Content plot, showing a cumulative percentage count of the proportion of library which has seen each of the adapter sequences at each position* 

 
The non-random per base sequence content caused by random k-mers cannot be corrected but should not affect the analysis.
The removal of the adapter contamination is possible.

Using Trim Galore! has been performed a quality trimming/filtering (phred score >20) and sequence trimming (for nextera trasposase).

First I prepared a tab-separated list of the fastq paired files

> fastq_gz_trim_list.txt

```
SRR5763563_1.fastq.gz	SRR5763563_2.fastq.gz
SRR5763564_1.fastq.gz	SRR5763564_2.fastq.gz
SRR5763565_1.fastq.gz	SRR5763565_2.fastq.gz
SRR5763566_1.fastq.gz	SRR5763566_2.fastq.gz
...............      ................
```
and script to create a parallel job list that runs trim-galore!

> Make_job_list_trim.sh

```
#!/bin/bash

[ $# -ne 2 ] && { echo -en "\n*** This script generates jobs for GNU parallel. *** \n\n Error Nothing to do, usage: < input tab delimited list > < output run list file >\n\n" ; exit 1; }
set -o pipefail

# Get command-line args
INPUT_LIST=$1
OUTPUT=$2

# Set counter
COUNT=1
END=$(wc -l $INPUT_LIST | awk '{print $1}')

echo " "
echo " * Input file is: $INPUT_LIST"
echo " * Number of runs: $END"
echo " * Output job list for GNU parallel saved to: $OUTPUT"
echo " "

# Main bit of command-line for job

CMD="trim_galore --nextera --paired"

# In alternative, to generate also BAM files
#CMD="kallisto quant -t 30 -b 30 --bias --pseudobam --genomebam --gtf ~/project/annotation/Mus_musculus.GRCm38.92.gtf --chromosomes ~/project/reference/chr_len.txt -i ~/project/indices/kallist$

# Main Loop
[ -e $OUTPUT ] && rm $OUTPUT
while [ $COUNT -le $END ];
do
    LINE=( $(awk "NR==$COUNT" $INPUT_LIST) )
    # Make file list
    echo "Working on $COUNT of $END, Files ${LINE[0]} ${LINE[1]}"
    echo "$CMD ${LINE[0]} ${LINE[1]} " >> $OUTPUT
    ((COUNT++))
done
```

```
# run script to create job list
sh Make_job_list_trim.sh fastq_gz_trim_list.txt trim_jobs_temp.txt

# converting carriage return into new line
tr '\r' '\n' < trim_jobs_temp.txt > trim_jobs.txt | rm trim_jobs_temp.txt

# run parallel on 8 threads (saving log)
parallel --progress --jobs 8 --joblog trim_joblog.txt < trim_jobs.txt

# rename generated trimmed files
rename "_val_1" "" *.fq.gz
rename "_val_2" "" *.fq.gz
PRE-TRIM vs POST-TRIM MultiQC Comparison
```

![](images/adapter_content_post.png)

*Fig.5 - FastQC Adapter Content plot after adapter trimming and quality filtering.*


Nextera transposase residual sequences have been removed  removed but a slight heterogeneity in sequence length distribution has been generated across the samples.

![](images/seq_len_dist_post.png)

*Fig. 6 - Sequence Length Distribution plot, showing the distribution of fragment sizes (read lengths) found.*



## 3.3. Reads alignment to transcriptome (Kallisto)


Befor to perform the pseudo-alignment against the mouse transcriptome, kallisto requires to create an index from a FASTA formatted file of target sequences

M. musculus transcriptome (cDNA) has been dowloaded from the Ensembl FTP repository

```
wget ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz /reference
```

Using this file a kallisto index has been generated with a default k-mer size of 31.

```
kallisto index  -i indices/kallisto/cdna.kidx reference/Mus_musculus.GRCm38.cdna.all.fa.gz 
```
```
[build] loading fasta file reference/transcripts/Mus_musculus.GRCm38.cdna.all.fa
[build] k-mer length: 31
[build] warning: clipped off poly-A tail (longer than 10) from 589 target sequences
[build] warning: replaced 3 non-ACGUT characters in the input sequence with pseudorandom nucleotides
[build] counting k-mers ... done.
[build] building target de Bruijn graph ...  done
[build] creating equivalence classes ...  done
[build] target de Bruijn graph has 691348 contigs and contains 97656136 k-mers
```

In order to produce the quantification, Kallisto has to be run for every paired-end sample, generating the "abundance" files to be imported in R and evaluate the differential expression.

kallisto can run on multiple instances using GNU parallel, in order to speed-up the process. 
To do so I prepared a tab-separated list of the samples

```
SRR5763563    SRR5763563_1.fq    SRR5763563_2.fq
SRR5763564    SRR5763564_1.fq    SRR5763564_2.fq
SRR5763565    SRR5763565_1.fq    SRR5763565_2.fq
SRR5763566    SRR5763566_1.fq    SRR5763566_2.fq
SRR5763567    SRR5763567_1.fq    SRR5763567_2.fq
.........	    ..........	        ...........

```

And a bash script to ingest this list and generate the jobs for parallel...


> Make_job_list_kallisto.sh

```
#!/bin/bash

[ $# -ne 2 ] && { echo -en "\n*** This script generates jobs for GNU parallel. *** \n\n Error Nothing to do, usage: < input tab delimited list > < output run list file >\n\n" ; exit 1; }
set -o pipefail

# Get command-line args
INPUT_LIST=$1
OUTPUT=$2

# Set counter
COUNT=1
END=$(wc -l $INPUT_LIST | awk '{print $1}')

echo " "
echo " * Input file is: $INPUT_LIST"
echo " * Number of runs: $END"
echo " * Output job list for GNU parallel saved to: $OUTPUT"
echo " "

# Main bit of command-line for job

CMD="kallisto quant -t 30 -b 30 --bias -i ~/project/indices/kallisto/cdna.kidx"

# In alternative, to generate also BAM files
#CMD="kallisto quant -t 30 -b 30 --bias --pseudobam --genomebam --gtf ~/project/annotation/Mus_musculus.GRCm38.92.gtf --chromosomes ~/project/reference/chr_len.txt -i ~/project/indices/kallisto/cdna.kidx"

# Main Loop
[ -e $OUTPUT ] && rm $OUTPUT
while [ $COUNT -le $END ];
do
    LINE=( $(awk "NR==$COUNT" $INPUT_LIST) )
    # Make file list
    echo "Working on $COUNT of $END Sample ID: ${LINE[0]}, Files ${LINE[1]} ${LINE[2]}"
    echo "$CMD -o ~/project/results/kallisto/${LINE[0]} ${LINE[1]} ${LINE[2]} " >> $OUTPUT
    ((COUNT++))
done
```

Other than index file and the FastQ files, additional options can be specified for Kallisto:
- --bias - learns parameters for a model of sequences specific bias and corrects the abundances accordlingly.
- --b - number of bootstrap for estimating the technical variance in samples (all bootstrap are compressed inside abundance.h5 output)
- --pseudobam - save pseudoalignments to transcriptome to BAM file
- --genomebam - Project pseudoalignments to genome sorted BAM file (used in combination with -g to specify the annotation GTF and -c to give a list of chromosome lenghts)

```
#Generating the job list
sh Make_job_list_kallisto.sh fastq_gz_kal_list.txt kallisto_jobs_temp.txt 

# converting carriage line to newline
tr '\r' '\n' < kallisto_jobs_temp.txt > kallisto_jobs.txt | rm kallisto_jobs_temp.txt

# running parallel for 8 threads (saving log)
nohup parallel --progress --jobs 8 --joblog kallisto_joblog.txt < kallisto_jobs.txt &
```

Outputs of kallisto are in folder "project/results/kallisto/".


## 3.4. Bioconductor analysis

The same analysis pipeline has been performed for the abundance files referring to the stages E3.5 and E4.5
following the guidlines published here:

>Lun, A.T.L., McCarthy, D.J., and Marioni, J.C. (2016). A step-by-step workflow for low-level analysis of single-cell RNA-seq data. F1000Research 5, 2122.

All the R files and outputs are insite the "analysis" folder, organized in subfolders depending on which samples the analysis has been carried out:
- all -> contains the R output of the whole dataset analysis
- e35 -> contains the R output of the stage E3.5 analysis
- e45 -> contains the R output of the stage E4.5 analysis
- report -> contains the R output and the image files of this report

The exemple code used for the stage E4.5 is reported below.

*****

Set working directory

```{R}
setwd("~/project/analysis/e45")
```

### 3.4.1. Data praparation

Loading required packages

```{R include=FALSE}
library(scater)
library(tximport)
library(scran)
library(limma)
```

Preparing E4.5 samples files

```{R}
metadata <- read.delim("samples_info.txt", check.names=FALSE, header=TRUE)
samples <- as.character(metadata$sample[metadata$stage == "E4.5"])
kal_dir <- "../../results/kallisto"
files <- file.path(kal_dir, samples, "abundance.h5")
```

Creating tx2gene to convert transcript_ids into gene_names in kallisto results

```{R include=FALSE}
library(EnsDb.Mmusculus.v79)
txdb <- EnsDb.Mmusculus.v79
keys <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys, "GENEID", "TXNAME")
```

Creating singleCellExperiment from Kallisto results using Tximport, collapsing the rows by genes (txOut = FALSE)

```{R echo=TRUE}
sce <- readTxResults(samples = samples, files = files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE, txOut = FALSE)
```

Query ensembldb using the Ensembl transcript ID (rownames) to retriver the annotation info to store into sce
By default BiomaRt uses the "mmusculus_gene_ensembl" dataset, retrieving "ensembl_transcript_id", "ensembl_gene_id", "mgi_symbol", "chromosome_name", "transcript_biotype", "transcript_start", "transcript_end"

```{R include=FALSE}
library(biomaRt)
sce <- getBMFeatureAnnos(sce, filters = "ensembl_gene_id")
```

Rename rownames with gene symbols

```{R}
rowData(sce)$ensembl_gene_id <- rownames(sce)
new.names <- rowData(sce)$mgi_symbol
missing.name <- is.na(new.names)
new.names[missing.name] <- rowData(sce)$ensembl_gene_id[missing.name]
dup.name <- new.names %in% new.names[duplicated(new.names)]
new.names[dup.name] <- paste0(new.names, "_", rowData(sce)$ensembl_gene_id)[dup.name]
rownames(sce) <- new.names

head(rownames(sce))
```

Define references for control_features to be used in calculateQCMetrics()

```{R}
mito <- which(rowData(sce)$chromosome_name=="MT")
```

Calculate cell metrics

```{r echo=TRUE}
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))

names(colData(sce))
```

### 3.4.2. Cell-based QC

First let's have an overview of the features counts across the cells

```{R echo=TRUE, fig.height=10, fig.width=10}
par(mfrow=c(2,2))
hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="",
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="",
     breaks=20, col="grey80", ylab="Number of cells")
plot(sce$total_features, sce$total_counts/1e6, xlab="Number of expressed genes", ylab="Library size (millions)")
plot(sce$total_features, sce$pct_counts_Mt, xlab="Number of expressed genes", ylab="Mitochondrial proportion (%)")
```
*Fig. 7 - Behaviour of each QC metric compared to the total number of expressed features. Each point represents a cell in the E4.5 dataset.*


Low-quality cells need to be removed to ensure that technical effects do not distort downstream analysis results. Two common measures of cell quality are the library size and the number of expressed features in each library. The library size is defined as the total sum of counts across all features. Cells with relatively small library sizes are considered to be of low quality as the RNA has not been efficiently captured (i.e., converted into cDNA and amplified) during library preparation. The number of expressed features in each cell is defined as the number of features with non-zero counts for that cell. Any cell with very few expressed genes is likely to be of poor quality as the diverse transcript population has not been successfully captured. 

### 3.4.3. Identifying outliers for each metrics

Outliers are defined based on the median absolute deviation (MADs) from the median value of each metric across all cells. We remove cells with log-library sizes that are more than 3 MADs below the median log-library size & cells where the log-transformed number of expressed genes is 3 MADs below the median value. Another measure of quality is the proportion of reads mapped to genes in the mitochondrial genome. High proportions are indicative of poor-quality cells ([@ilicic_classification_2016], [@islam_quantitative_2014]), possibly because of increased apoptosis and/or loss of cytoplasmic RNA from lysed cells.

```{r echo=TRUE}
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)

feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)

mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=3, type="higher")
```

Subsetting by column will retain only the high-quality cells that pass each filter described. 
Let's see the number of cells removed by each filter as well as the total number of retained cells.

```{R echo=TRUE}
keep <- !(libsize.drop | feature.drop | mito.drop)
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), ByMito=sum(mito.drop), Remaining=sum(keep))
```

Now subset the SingleCellExperiment object to retain only the putative high-quality cells. 

```{R}
sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
```

### 3.4.4. Classification of cell-cycle phase

Using a method reported in literature [@scialdone_computational_2015] is possible to classify cells into cell cycle phases based on the gene expression data. Using a training dataset, the sign of the difference in expression between two genes was computed for each pair of genes. Pairs with changes in the sign across cell cycle phases were chosen as markers. Cells in a test dataset can then be classified into the appropriate phase, based on whether the observed sign for each marker pair is consistent with one phase or another.
The cyclone package contains a pre-trained set of marker pairs for mouse data, which can be loaded in the the readRDS function. We use the Ensembl identifiers for each gene in our dataset to match up with the names in the pre-trained set of gene pairs.

```{R echo=TRUE}
library(scran)
set.seed(100)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ensembl_gene_id)

plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)
```
*Fig. 8 - Cell cycle phase scores from applying the pair-based classifier on the dataset. Each point represents a cell, plotted according to its scores for G1 and G2/M phases.*


How many cells have been predicted as dividing?

```{R echo=TRUE}
sce$phases <- assignments$phases
table(sce$phases)
```

Each cell is assigned a score for each phase, with a higher score corresponding to a higher probability that the cell is in that phase. We focus on the G1 and G2/M scores as these are the most informative for classification. Cells are classified as being in G1 phase if the G1 score is above 0.5 and greater than the G2/M score; in G2/M phase if the G2/M score is above 0.5 and greater than the G1 score; and in S phase if neither score is above 0.5. We can use the assigned phase as a blocking factor in downstream analyses. This protects against cell cycle effects without discarding information.

### 3.4.5. Examining gene-level expression metrics

Inspecting the most highly expressed genes.
This should generally be dominated by constitutively expressed transcripts, such as those for ribosomal or mitochondrial proteins. The presence of other classes of features may be cause for concern if they are not consistent with expected biology. For example, the absence of ribosomal proteins and/or the presence of their pseudogenes are indicative of suboptimal alignment.

```{R echo=TRUE, fig.height=10, fig.width=8}
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))

plotQC(sce, type = "highest-expression", n=50) + fontsize
```
*Fig. 9 - Percentage of total counts assigned to the top 50 most highly-abundant features in the dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell, while circles are coloured according to whether the feature is labelled as a control feature.*


### 3.4.6. Filtering out low-abundance genes

The average count for each gene, is computed across all cells in the dataset using the calcAverage() function, which also performs some adjustment for library size differences between cells. Typically can be observed a peak of moderately expressed genes following a plateau of lowly expressed genes.

```{R echo=TRUE}
ave.counts <- calcAverage(sce, use_size_factors=FALSE)

hist(log10(ave.counts), breaks=100, main="", col="grey80", xlab=expression(Log[10]~"average count"))
abline(v=log10(1), col="blue", lwd=2, lty=2)
```
*Fig. 10 - Histogram of log-average counts for all genes in the dataset. The filter threshold is represented by the blue line.*


All the genes with average counts less than 1 should be removed

```{R echo=TRUE}
demo.keep <- ave.counts >= 1
summary(demo.keep)
```

The number of TRUE values corresponds to the number of retained rows/genes after filtering.
Apply the threshold and create a filtered sce

```{R}
filtered.sce <- sce[demo.keep,]
```

### 3.4.7. Examine number of cells that express each gene

Genes expressed in very few cells are often uninteresting as they are driven by amplification artifacts.

```{R echo=TRUE}
num.cells <- nexprs(sce, byrow=TRUE)

smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells", xlab=expression(Log[10]~"average count"))
```
*Fig. 11 - The number of cells expressing each gene in the dataset, plotted against the log-average count. Intensity of colour corresponds to the number of genes at any given location.*


Genes that are not expressed in any cell are removed to reduce teh computational work in downstream steps

```{R echo=TRUE}
to.keep <- num.cells > 0
sce <- sce[to.keep,]
summary(to.keep)
```

### 3.4.8. Normalization of cell-specific biases

Any systematic difference in count size across the non-DE majority of genes between two cells is assumed to represent bias and is removed by scaling. More specifically, "size factors" are calculated that represent the extent to which counts should be scaled in each library.
Single-cell data can be problematic for these bulk data-based methods due to the dominance of low and zero counts. To overcome this, counts are pooled from many cells to increase the count size for accurate size factor estimation. Pool-based size factors are then "deconvolved" into cell-based factors for cell-specific normalization.

```{R echo=TRUE}
# check how many samples remained
dim(sce)
```

```{R echo=TRUE}
# pooling groups of cells to calculate size factor
sce <- computeSumFactors(sce, sizes=c(10, 20, 30, 40, 50, 60, 70, 80, 94))
summary(sizeFactors(sce))
```

```{R echo=TRUE}
plot(sizeFactors(sce), sce$total_counts/1e6, log="xy",
     ylab="Library size (millions)", xlab="Size factor")
```
*Fig. 12 - Size factors from deconvolution, plotted against library sizes for all cells in the dataset. Axes are shown on a log-scale.*


Applying the size factors to normalize gene expression

```{R}
sce <- normalize(sce)
```

### 3.4.9. Modelling the technical noise in gene expression

Testing for significantly positive biological components
(The biological component for each gene is defined as the difference between its total variance and the fitted value of the trend)
First, the biological and technical components of the gene-specific variance are computed...

```{R}
var.fit <- trendVar(sce, parametric=TRUE, loess.args=list(span=0.3), use.spikes=FALSE)
```

...then biological and technical variance are decomposed

```{R echo=TRUE}
var.out <- decomposeVar(sce, var.fit)

head(var.out)
```

We can visually inspect the trend

```{r echo=TRUE}
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression", main="E4.5")
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
```
*Fig. 13 - Variance of normalized log-expression values for each gene in the dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances.*


The wave-like shape is typical of the mean-variance trend for log-expression values. A linear increase in the variance is observed as the mean increases from zero, as larger variances are possible when the counts increase. At very high abundances, the effect of sampling noise decreases due to the law of large numbers, resulting in a decrease in the variance.

###  3.4.10. Visualizing data in low-dimensional space

Once the technical noise is modelled, principal components analysis can be used to remove random technical noise. Consider that each cell represents a point in the high-dimensional expression space, where the spread of points represents the total variance. PCA identifies axes in this space that capture as much of this variance as possible. Each axis is a principal component (PC), where any early PC will explain more of the variance than a later PC.

It is assumed that biological processes involving co-regulated groups of genes will account for the most variance in the data. If this is the case, this process should be represented by one or more of the earlier PCs. In contrast, random technical noise affects each gene independently and will be represented by later PCs. The denoisePCA() function removes later PCs until the total discarded variance is equal to the sum of technical components for all genes used in the PCA.

```{R echo=TRUE}
sce <- denoisePCA(sce, technical=var.out, assay.type="logcounts")
dim(reducedDim(sce, "PCA"))
```

The function returns a SingleCellExperiment object containing the PC scores for each cell in the reducedDims slot. The aim is to eliminate technical noise and enrich for biological signal in the retained PCs. This improves resolution of the underlying biology during downstream procedures such as clustering.

Now relationships between cells can be visualized by constructing pairwise PCA plots for the first three components

```{r}
plotReducedDim(sce, use_dimred="PCA", ncomponents=3) + fontsize
```
*Fig. 14 - Pairwise PCA plots of the first three PCs in the E4.5 dataset, constructed from normalized log-expression values of genes with positive biological components*


t-SNE (t-stochastic neighbour embedding) tends to work better than PCA for separating cells in more diverse populations. This is because the former can directly capture non-linear relationships in high-dimensional space, whereas the latter must represent them on linear axes.

use_dimred="PCA" can be set to perform the t-SNE on the low-rank approximation of the data, allowing the algorithm to take advantage of the previous denoising step.

```{r echo=TRUE, fig.height=4, fig.width=10}
set.seed(100)
out5 <- plotTSNE(sce, run_args=list(use_dimred="PCA", perplexity=5)) + fontsize + ggtitle("Perplexity = 5")

set.seed(100)
out10 <- plotTSNE(sce, run_args=list(use_dimred="PCA", perplexity=10)) + fontsize + ggtitle("Perplexity = 10")

set.seed(100)
out20 <- plotTSNE(sce, run_args=list(use_dimred="PCA", perplexity=20)) + fontsize + ggtitle("Perplexity = 20")

multiplot(out5, out10, out20, cols=3)
```
*Fig. 15 - t-SNE plots constructed from the denoised PCs in the E4.5 dataset, using a range of perplexity values.*


Scripts should set a seed to ensure that the chosen results are reproducible. It is also advisable to test different settings of the "perplexity" parameter as this will affect the distribution of points in the low-dimensional space.

Now run runTSNE() with a perplexity of 20 to store the t-SNE results inside our SingleCellExperiment object. This avoids repeating the calculations whenever we want to create a new plot with plotTSNE(), as the stored results will be used instead.

```{R echo=TRUE}
set.seed(100)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=20)
reducedDimNames(sce)
```


### 3.4.11. Clustering cells into putative subpopulations

The denoised log-expression values are used to cluster cells into putative subpopulations. Specifically, hierarchical clustering is performed on the Euclidean distances between cells, using Ward's criterion to minimize the total variance within each cluster. This yields a dendrogram that groups together cells with similar expression patterns across the chosen genes.

```{r}
pcs <- reducedDim(sce, "PCA")
my.dist <- dist(pcs)
my.tree <- hclust(my.dist, method="ward.D2")
```

Clusters are explicitly defined by applying a dynamic tree cut (Langfelder, Zhang, and Horvath 2008) to the dendrogram. This exploits the shape of the branches in the dendrogram to refine the cluster definitions.

```{r}
library(dynamicTreeCut)
my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), 
    minClusterSize=10, verbose=0))
```

Let's see the distribution of cells in each cluster with respect to known factors.

```{r echo=TRUE}
table(my.clusters)
```

Visualize the cluster assignments for all cells on the t-SNE plot

```{r echo=TRUE}
sce$cluster <- factor(my.clusters)

plotTSNE(sce, colour_by="cluster") + fontsize
```
*Fig. 16 - t-SNE plot of the denoised PCs of the E4.5 dataset*


The separatedness of the clusters is checked using the silhouette width. Cells with large positive silhouette widths are closer to other cells in the same cluster than to cells in different clusters. Conversely, cells with negative widths are closer to other clusters than to other cells in the cluster to which it was assigned.

```{R echo=TRUE}
library(cluster)
clust.col <- scater:::.get_palette("tableau10medium")
sil <- silhouette(my.clusters, dist = my.dist)
sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
sil.cols <- sil.cols[order(-sil[,1], sil[,3])]

plot(sil, main = paste(length(unique(my.clusters)), "clusters"), border=sil.cols, col=sil.cols, do.col.sort=FALSE)

```

*Fig. 17 - Barplot of silhouette widths for cells in each cluster.*


Clusters 2 and 3 have very positive widths.


### 3.4.12. Detecting marker genes between clusters

Once putative subpopulations are identified by clustering, marker genes for each cluster can be identified using the findMarkers function. This fits a linear model to the log-expression values for each gene using limma. The aim is to test for DE in each cluster compared to the others. The top DE genes are likely to be good candidate markers as they can effectively distinguish between cells in different clusters.

For each cluster, the DE results of the relevant comparisons are consolidated into a single output table. This allows a set of marker genes to be easily defined by taking the top DE genes from each pairwise comparison between clusters. 

For example, to construct a marker set for cluster 2 from the top 10 genes of each comparison, marker.set is filtered to retain rows with Top less than or equal to 10. Other statistics are also reported for each gene, including the adjusted p-values (see below) and the log-fold changes relative to every other cluster.

```{r echo=TRUE}
markers <- findMarkers(sce, my.clusters)
marker.set <- markers[["2"]]

head(marker.set, 10)
```

Save the list of candidate marker genes for further examination.

```{r}
write.table(marker.set, file="e45_markers_cl2.tsv", sep="\t",  quote=FALSE, row.names=TRUE)
```

The expression profiles of the top candidates can be visualized to verify that the DE signature is robust. 
The clusters assignment for every sample and some markers of pluripotency are reported on the top rows.

```{R echo=TRUE, include=FALSE, fig.height=8, fig.width=20}
top.markers <- rownames(marker.set)[marker.set$Top <= 10]

library(pheatmap)

plotHeatmap(sce, features=top.markers, columns=order(sce$cluster), 
    colour_columns_by=c("cluster","Pou5f1", "Nanog", "Esrrb","Otx2"),
    cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-5, 5),
    main = "E4.5 \n top.markers + pluriset",
    cellwidth=10, width=18, height=10
    )
```


### 3.4.13 Identifying correlated gene pairs within the HVGs using Spearman's rho

HVGs are defined as genes with biological components that are significantly greater than zero at a false discovery rate (FDR) of 5%. These genes are interesting as they drive differences in the expression profiles between cells, and should be prioritized for further investigation. In addition, we only consider a gene to be a HVG if it has a biological component greater than or equal to 0.5. For transformed expression values on the log2 scale, this means that the average difference in true expression between any two cells will be at least 2-fold. (This reasoning assumes that the true log-expression values are Normally distributed with variance of 0.5. The root-mean-square of the difference between two values is treated as the average log2-fold change between cells and is equal to unity.) We rank the results by the biological component to focus on genes with larger biological variability.

```{r echo=TRUE}
hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),]
nrow(hvg.out)
```

```{r}
write.table(file="e45_hvg.tsv", hvg.out, sep="\t", quote=FALSE, col.names=NA)
head(hvg.out)
```

The distribution of expression values for the top HVGs can be seen to ensure that the variance estimate is not being dominated by one or two outlier cells

```{r echo=TRUE}
plotExpression(sce, rownames(hvg.out)[1:10]) + fontsize
```
*Fig. 18 - Violin plots of normalized log-expression values for the top 10 genes with the largest biological components in the E4.5 dataset*


Now HVGs that are highly correlated with one another can be identified. This distinguishes between HVGs caused by random noise and those involved in driving systematic differences between subpopulations. Correlations between genes are quantified by computing Spearman's rho, which accommodates non-linear relationships in the expression values. Gene pairs with significantly large positive or negative values of rho are identified using the correlatePairs function. 

Calculating correlations for all possible gene pairs would require too much computational time and increase the severity of the multiple testing correction. It may also prioritize uninteresting genes that have strong correlations but low variance, e.g., tightly co-regulated house-keeping genes.

In this study the main interested is finding the genes cotrelated to Nanog. 
In the E4.5 dataset after gene level summarization (tximport), the technical component of the variance for Nanog expression (6.29) results very high compared to the biological component (-2.76) causing a very low significance (p-value = 0.99), so Nanog is not included among the HVGs (FRD = 1). Nanog has been explicitly included in a coercive way into the CorrelatePairs input list together with the HVGs, only for evaluation purpose on the output list of Nanog correlated genes. (The CorrelatePairs is repeated below at transcript level where the variant Nanog-203 has a good biological variance and is included among the HVGs).

```{R echo=TRUE, include=FALSE}
# computing Spearman's rank test
set.seed(100)
var.cor <- correlatePairs(sce, iters=1e8, subset.row=c(rownames(hvg.out), "Nanog"))

# exporting table of all correlated genes
write.table(file="e45_cor.tsv", var.cor, sep="\t", quote=FALSE, row.names=FALSE)

# extracting only pairs including Nanog
nanog.rows <- var.cor$gene1 == "Nanog" | var.cor$gene2 == "Nanog"
var.cor.nanog <- var.cor[nanog.rows,]
# exporting table of Nanog correlated genes
write.table(file="e45_cor_nanog.tsv", var.cor.nanog, sep="\t", quote=FALSE, row.names=FALSE)
```

The significance of each correlation is determined using a permutation test. For each pair of genes, the null hypothesis is that the expression profiles of two genes are independent. Shuffling the profiles and recalculating the correlation yields a null distribution that is used to obtain a p-value for each observed correlation value. Correction for multiple testing across many gene pairs is performed by controlling the FDR at 5%. 

```{r}
sig.cor <- var.cor$FDR <= 0.05
summary(sig.cor)
```

### 3.4.14. Using Nanog correlated HVGs for further data exploration

Have a look to the top Nanog vorrelated genes, by restricting the list of correlated pairs by significance (FDR <= 5%) and high strength of association (rho > [0.4]).

```{r}
sig.cor.nanog <- var.cor.nanog$FDR <= 0.05 & abs(var.cor.nanog$rho) >= 0.4

summary(sig.cor.nanog)
```

The expression profiles of these top 227 Nanog correlated HVGs can be visualized with a heatmap. 
All expression values are mean-centred for each gene to highlight the relative differences in expression between cells.

```{r echo=TRUE, include=FALSE, fig.height=60, fig.width=20}
# define the list of genes significantly correlated to Nanog
chosen.nanog <- unique(c(var.cor.nanog$gene1[sig.cor.nanog], var.cor.nanog$gene2[sig.cor.nanog]))
chosen.nanog <- chosen.nanog[chosen.nanog != "Nanog"]

# plot the expression for these genes
library(pheatmap)
plotHeatmap(sce, features=chosen.nanog, columns=order(sce$cluster), 
    colour_columns_by=c("Nanog","cluster"),
    cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-5, 5),
    main = "E4.5 \n Top Nanog Correlated Genes (forced)",
    #filename="e45_heatmap_sig_cor_nanog.pdf",
    cellwidth=10, width=18, height=80
    )
```


### 3.4.15. Finding Nanog correlated genes in the epiblast subset of cells

From the heatmak of the top marker genes identified by tSNE, We can observe that the clusters 2 represent a population of pluripotent cells (Oct4+) with epiblast features (Sox2+, Fgf4+) while the rest of cells could be attributed the the primitive endoderm (Sox17+, Dab2+). It is also known by literature [@aksoy_oct4_2013] that at this stage Oct4 switches partner from Sox2 to Sox17 contributing to the commitment of the PrE differentiation. This hipotesis can be verified by comparing visually the expression levels of the markers for these two complementary populations.

```{r fig.height=5, fig.width=10}
tsne.fgf4 <- plotTSNE(sce, colour_by="Fgf4", size_by="Pou5f1") + fontsize
tsne.sox2 <- plotTSNE(sce, colour_by="Sox2", size_by="Pou5f1") + fontsize

multiplot(tsne.fgf4, tsne.sox2, cols=2)
```
*Fig. 19 - t-SNE plot of the denoised PCs of the E4.5 dataset showing the expression distribution of epiblast markers across samples*


```{r fig.height=5, fig.width=10}
tsne.dab2 <- plotTSNE(sce, colour_by="Dab2", size_by="Pou5f1") + fontsize
tsne.sox17 <- plotTSNE(sce, colour_by="Sox17", size_by="Pou5f1") + fontsize

multiplot(tsne.dab2, tsne.sox17, cols=2)
```
*Fig. 20 - t-SNE plot of the denoised PCs of the E4.5 dataset showing the expression distribution of primitive endoderm markers across samples*


So the analysis can be concentrated on the subset of cells of the cluster 2, showing features of the epiblast.

```{r}
epi.sce <- sce[,sce$cluster == 2]
dim(epi.sce)
```

It is possible to define a set of epiblast (Epi) and primitive endoderm (PrE) markers chosen among the top markers assigned to cluster 2 by tSNE and comparing the distribution of these markers between the full dataset and the putative epiblast subset.

```{r fig.height=5, fig.width=10}
episet <- c("Fgf4", "Tdgf1", "Igfbp2", "Morc1", "Slc7a3", "Nanog", "Esrrb", "Sox2", "Sox17")

episet.pre.sub <- plotExpression(sce, episet, colour_by="Pou5f1") + fontsize + ggtitle("E4.5")
episet.post.sub <- plotExpression(epi.sce, episet, colour_by="Pou5f1") + fontsize + ggtitle("E4.5\nEPIBLAST")
multiplot(episet.pre.sub, episet.post.sub, cols = 2)
```
*Fig. 21 - Violin plots of normalized log-expression values for a set of epiblast markers across cells before and after subsetting of tSNE cluster 2*


The strongest epiblast marker identified, showing the best separation between the high Oct4 and low Oct4 populations is Tdgf1, which is coding for a protein involved in Nodal signaling and plays a role in the determination of the epiblastic cells that subsequently give rise to the mesoderm.

It is required to normalize again for cell-specific biases in the subset

```{R}
epi.sce <- computeSumFactors(epi.sce, sizes=c(5, 10, 15, 20, 25))
summary(sizeFactors(epi.sce))
```


Apply the size factors to normalize gene expression

```{R}
epi.sce <- normalize(epi.sce)
```

Identify correlated gene pairs within the epiblast subpopulation

```{R}
set.seed(100)
epi.var.cor <- correlatePairs(epi.sce, iters=1e8, subset.row=c(rownames(hvg.out), "Nanog"))
epi.nanog.rows <- epi.var.cor$gene1 == "Nanog" | epi.var.cor$gene2 == "Nanog"
epi.var.cor.nanog <- epi.var.cor[epi.nanog.rows,]
write.table(file="epi.e45_cor.tsv", epi.var.cor, sep="\t", quote=FALSE, row.names=FALSE)
write.table(file="epi.e45_cor_nanog.tsv", epi.var.cor.nanog, sep="\t", quote=FALSE, row.names=FALSE)
head(epi.var.cor.nanog)
```

```{r}
epi.sig.cor <- epi.var.cor$FDR <= 0.05
summary(epi.sig.cor)
```

```{r}
epi.sig.cor.nanog <- epi.var.cor.nanog$FDR <= 0.05
summary(epi.sig.cor.nanog)
```


### 3.4.16. Analysis at transcript level


Create singleCellExperiment from Kallisto abundance files using Tximport without collapsing the rows by gene (txOut = TRUE)

```{R}
sce2 <- readTxResults(samples = samples, files = files, type = "kallisto", txOut = TRUE)
```

Retrieve annotation info from ensembldb and store it into sce2

```{R include=FALSE}
library(biomaRt)
sce2 <- getBMFeatureAnnos(sce2, filters = "ensembl_transcript_id")
```

Rename rownames to include gene symbols

```{R}
rownames(sce2) <- paste0(rowData(sce2)$mgi_symbol, "_",rownames(sce2))
head(rownames(sce2))
```

Define references for control_features to be used in calculateQCMetrics()

```{R}
mito2 <- which(rowData(sce2)$chromosome_name=="MT")
```

Calculate cell metrics

```{r}
sce2 <- calculateQCMetrics(sce2, feature_controls=list(Mt=mito2))

names(colData(sce2))
```

Define the quality control metrics and identifying outliers for each metric

```{r}
libsize.drop2 <- isOutlier(sce2$total_counts, nmads=3, type="lower", log=TRUE)

feature.drop2 <- isOutlier(sce2$total_features, nmads=3, type="lower", log=TRUE)

mito.drop2 <- isOutlier(sce2$pct_counts_Mt, nmads=3, type="higher")

keep2 <- !(libsize.drop2 | feature.drop2 | mito.drop2)

# Subsetting by column will retain only the high-quality cells that pass each filter

sce2 <- sce2[,!(libsize.drop2 | feature.drop2 | mito.drop2)]
```

Remove genes that are not expressed in any cell to reduce computational work in downstream steps

```{R}
num.cells2 <- nexprs(sce2, byrow=TRUE)
to.keep2 <- num.cells2 > 0
sce2 <- sce2[to.keep2,]
```

Normalization of cell-specific biases

```{R}
dim(sce2)
```

```{R}
sce2 <- computeSumFactors(sce2, sizes=c(5, 10, 20, 40, 60, 80, 94))
sce2 <- normalize(sce2)
```

Model the technical noise in gene expression

```{R}
var.fit2 <- trendVar(sce2, parametric=TRUE, loess.args=list(span=0.3), use.spikes=FALSE)
var.out2 <- decomposeVar(sce2, var.fit2)
```

Visualize data in low-dimensional space

```{R}
sce2 <- denoisePCA(sce2, technical=var.out2, assay.type="logcounts")
```

```{r fig.height=3.5, fig.width=10}
set.seed(100)
out5.2 <- plotTSNE(sce2, run_args=list(use_dimred="PCA", perplexity=5)) + fontsize + ggtitle("Perplexity = 5")

set.seed(100)
out10.2 <- plotTSNE(sce2, run_args=list(use_dimred="PCA", perplexity=10)) + fontsize + ggtitle("Perplexity = 10")

set.seed(100)
out20.2 <- plotTSNE(sce2, run_args=list(use_dimred="PCA", perplexity=20)) + fontsize + ggtitle("Perplexity = 20")

multiplot(out5.2, out10.2, out20.2, cols=3)
```
*Fig. 22 - t-SNE plots constructed from the denoised PCs in the E4.5 transcripts level dataset, using a range of perplexity values*


```{R}
set.seed(100)
sce2 <- runTSNE(sce2, use_dimred="PCA", perplexity=10)
```

Clustering cells into putative subpopulations

```{r}
pcs2 <- reducedDim(sce2, "PCA")
my.dist2 <- dist(pcs2)
my.tree2 <- hclust(my.dist2, method="ward.D2")
```

```{r}
library(dynamicTreeCut)
my.clusters2 <- unname(cutreeDynamic(my.tree2, distM=as.matrix(my.dist2), 
    minClusterSize=10, verbose=0))
```

```{r}
sce2$cluster <- factor(my.clusters2)

plotTSNE(sce2, colour_by="cluster") + fontsize
```
*Fig. 23 - t-SNE plot of the denoised PCs of the E4.5 transcripts level dataset*


```{R}
library(cluster)
clust.col <- scater:::.get_palette("tableau10medium")
sil2 <- silhouette(my.clusters2, dist = my.dist)
sil.cols2 <- clust.col[ifelse(sil2[,3] > 0, sil2[,1], sil2[,2])]
sil.cols2 <- sil.cols2[order(-sil2[,1], sil2[,3])]

plot(sil2, main = paste(length(unique(my.clusters2)), "clusters"), border=sil.cols2, col=sil.cols2, do.col.sort=FALSE)
```
*Fig. 24 - Barplot of silhouette widths for cells in each cluster for E4.5 transcipts level dataset*


Detecting marker transcripts between clusters

```{r}
markers2 <- findMarkers(sce2, my.clusters2)
marker.set2 <- markers2[["2"]]
top.markers2 <- rownames(marker.set2)[marker.set2$Top <= 10]
```

Define row filters for transcripts of pluripotency genes

```{R}
nanog <- rownames(sce2)[grep("Nanog_", rownames(sce2))]
esrrb <- rownames(sce2)[grep("Esrrb_", rownames(sce2))]
oct4 <- rownames(sce2)[grep("Pou5f1_", rownames(sce2))]
sox2 <- rownames(sce2)[grep("Sox2_", rownames(sce2))]

pluriset <- c(oct4, sox2, nanog, esrrb)
```

Visualize the expression distribution of cluster 2 top markers together with the pluripotency gene transcripts

```{R eval=FALSE, fig.height=10, fig.width=20, include=FALSE}
library(pheatmap)
plotHeatmap(sce2, features=top.markers2, columns=order(sce2$cluster), 
    colour_columns_by=c(pluriset, "cluster"),
    cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-5, 5), 
    cellwidth=10, width=20, height=10,
    #filename = "e45.mRNA_top_markers_heatmap_clusters.pdf",
    main = "E4.5 mRNAs \n top.markers + pluriset"
    )
```

Identifying HVGs from the normalized log-expression

```{R eval=FALSE}
hvg.out2 <- var.out2[which(var.out2$FDR <= 0.05 & var.out2$bio >= 0.5),]
hvg.out2 <- hvg.out2[order(hvg.out2$bio, decreasing=TRUE),]
write.table(file="e45_mRNA_HVGs.tsv", hvg.out2, sep="\t", quote=FALSE, col.names=NA)
```

Check if some Nanog transcripts have been included among the HVGs

```{r}
rownames(hvg.out2)[grep("Nanog_", rownames(hvg.out2))]
```

Identify correlated transcripts using Spearman's ranking

```{R eval=FALSE}
set.seed(100)
var.cor2 <- correlatePairs(sce2, iters=1e8, subset.row=rownames(hvg.out2))
write.table(file="e45_mRNA_cor.tsv", var.cor2, sep="\t", quote=FALSE, row.names=FALSE)
```

Retrieve transcripts correlated with Nanog from correlate pairs

```{r}
nanog.rows2 <- var.cor2$gene1 == "Nanog_ENSMUST00000012540.4" | var.cor2$gene2 == "Nanog_ENSMUST00000012540.4"
var.cor.nanog2 <- var.cor2[nanog.rows2,]
write.table(file="e45_mRNA_cor_nanog.tsv", var.cor.nanog2, sep="\t", quote=FALSE, row.names=FALSE)
```

Find the top Nanog correlated pairs by significance (FDR <= 5%) and high strength of association (rho > [0.4]).

```{r eval=FALSE}
sig.cor.nanog2 <- var.cor.nanog2$FDR <= 0.05 & abs(var.cor.nanog2$rho) >= 0.4

summary(sig.cor.nanog2)
```

Visualize the expression distribution of the top Nanog correlated transcripts

```{r eval=FALSE, fig.height=50, fig.width=20, include=FALSE}
chosen.nanog2 <- unique(c(var.cor.nanog2$gene1[sig.cor.nanog2], var.cor.nanog2$gene2[sig.cor.nanog2]))
chosen.nanog2 <- chosen.nanog2[chosen.nanog2 != "Nanog_ENSMUST00000012540.4"]

library(pheatmap)
plotHeatmap(sce2, features=chosen.nanog2, columns=order(sce2$cluster), 
    colour_columns_by=c("Nanog_ENSMUST00000012540.4", "cluster"),
    cluster_cols=FALSE, center=TRUE, symmetric=TRUE, zlim=c(-5, 5),
    main = "E4.5 \n Top Nanog Correlated mRNAs",
    #filename="e45.mRNA_heatmap_sig_cor_nanog.pdf",
    cellwidth=10, width=20, height=50
    )
```

subset SingleCellExperiment to Cluster 2 (Epiblast)

```{r}
epi.sce2 <- sce2[,sce2$cluster == 2]

dim(epi.sce2)
```

```{R}
epi.sce2 <- computeSumFactors(epi.sce2, sizes=c(5, 10, 15, 20, 24))
summary(sizeFactors(epi.sce2))
```

Applying the size factors to normalize gene expression

```{R}
epi.sce2 <- normalize(epi.sce2)
```

Identify correlated HVGs transcripts in the epiblast using Spearman's ranking

```{R}
set.seed(100)
epi.var.cor2 <- correlatePairs(epi.sce2, iters=1e8, subset.row=rownames(hvg.out2))
write.table(file="epi.e45_mRNA_cor.tsv", epi.var.cor2, sep="\t", quote=FALSE, row.names=FALSE)
head(epi.var.cor2)
```

Retrieve epiblast transcripts correlated with Nanog from correlate pairs

```{r}
epi.nanog.rows2 <- epi.var.cor2$gene1 == "Nanog_ENSMUST00000012540.4" | epi.var.cor2$gene2 == "Nanog_ENSMUST00000012540.4"
epi.var.cor.nanog2 <- epi.var.cor2[epi.nanog.rows2,]
write.table(file="epi.e45_mRNA_cor_nanog.tsv", epi.var.cor.nanog2, sep="\t", quote=FALSE, row.names=FALSE)
head(epi.var.cor2)
```

Find the top Nanog correlated pairs by significance (FDR <= 20%)

```{r}
epi.sig.cor.nanog2 <- epi.var.cor.nanog2$FDR <= 0.2 

summary(epi.sig.cor.nanog2)
```

Compute the correlate pair within the Nanog and Esrrb variants

```{r}
set.seed(100)
pluri.cor <- correlatePairs(sce2, subset.row=c(nanog, esrrb))
sig.pluri.cor <- pluri.cor$FDR <= 0.05
write.table(file="e45_mRNA_Nanog-Esrrb_cor.tsv", pluri.cor[sig.pluri.cor,], sep="\t", quote=FALSE, row.names=FALSE)
```

# 4. Output and discussion

From the preliminary quality control on the sequencing data performed with "FastQC", no particular issues have been identified, and the overall quality of the reads appeared to be good, except an adapter contamination deriving from the Illumina mate pair libraries constructed using the Nextera protocol.
The raw scRNA-seq data have been processed with "Trim Galore" for quality filtering of the reads with a low-quality base calls and the adapter trimming.

Following the QC, in alternative to the classic alignment based tools (such as tophat, STAR, bowtie, HISAT) the allignment free trasciptome quantification has been preferred due to the purpose of this study (involving a differential expression analysis, without necessity to identify novel transcripts). The pseudoallignment tool used "kallisto", break up reads into k-mers before assigning them to transcripts.
This results in a substantial gain in speed compared to the alignment based workflows. The workflows also differ in how the expression abundance is estimated, enabling quantification on transcript level [@everaert_benchmarking_2017].

After the transcriptome index preparation and the abundance quantification, the kallisto output has been imported into the R environmente using "tximport", generating SingleCellExperiment (sce) files for the Bioconductor analysis pipeline.
The analysis output reported here is focused on the peri-implatation stages E3.5 and E4.5.
The full dataset of the original study has been subjected to pseudoalignment with the mouse transcriptome and a principal component analysis has been performed, showing that the dataset separates by developmental stages.


![](images/all_PCA_stages.png)

*Fig. 25 - PCA plot for 2 dimensions, after QC and normalization of cell-specific biases, of the full dataset in the original study of Mohammed et al. 2017. The cells has been colored by developmental stage.*

By colouring differentially the cells on the PCA plot, it is possible to inspect the expression profile associated with gene markers reported in literature [@tam_gene_2007] allowing the identification of different lineages.

![](images/all_PCA_pluriset.png)

*Fig. 26 - PCA plot of the full dataset colored by gene expression levels (logcounts) of selected marker genes - Nanog (ICM/epiblast), Gata6 (PrE/VE), Esrrb (naive pluripotency). Cells are also sized by expression levels of Oct4 (core pluripotency).*


Examining the figures 30 and 31, the population of cells at E3.5 appear homogeneous, without a distinct lineage indentities. But colouring some marker genes (Fig. 31), it is clearly visible a separation of two populations with different features, among the pluripotent cells (high Oct4) of the stage E4.5, characterized by opposite expression profiles of Nanog and Gata6, respectively associated to the Epiblast (Epi) and primitive endoderm (PrE) in previous studies [@tam_gene_2007]. Esrrb expression distribution is concentrated on the stages E3.5 and E4.5 and the expression levels seem to follow the newely forming epiblast cell population.

Single cell consensus clustering (SC3) has been used to explore the data, find a reasonable estimate of the number of clusters and calculate the biological features based on the identified cell clusters:

![](images/all_sc3_de_k4.png)

*Fig. 27 - Heatmap showing the differential expression between clusters coincident to developmental stages, calculated using non-parametric Kruskal-Wallis test.*


After visually examining the consensus matrices generated for a range of clusters between 2 and 6, the optimal separation of 4 clusters has been chosen for downstream analysis, because segregates better with the four developmental stages constituing the dataset.
SC3 has been used to create a list of all differentially expressed (DE) genes between the clusters, with adjusted p-values < 0.01 and plots gene expression profiles of the 50 genes with the lowest p-values.
Interestingly, in coincidence with the stage E4.5, a subpopulation with strong positive DE for Dab2 is visible.
The endocytic adaptor protein Dab2 mediates directional vesicular trafficking required for the genesis of an apical polarity, which is known to be crucial for the sorting and positioning of the PrE cells at the surface of the inner cell mass (ICM) during the embryo implantation [@moore_primitive_2014].

Based on the mean expression values of the genes, marker genes for each cluster/stage have been defined.

![](images/all_sc3_markers_k4.png)

*Fig. 28 - Heatmap showing the marker genes identified for each cluster. To find marker genes, for each gene a binary classifier is constructed based on the mean cluster expression values. The classifier prediction is then calculated using the gene expression ranks. The area under the receiver operating characteristic (ROC) curve is used to quantify the accuracy of the prediction. A p-value is assigned to each gene by using the Wilcoxon signed rank test. By default the genes with the area under the ROC curve (AUROC) > 0.85 and with the p-value < 0.01 are selected and the top 10 marker genes of each cluster are visualized in this heatmap.*

Notably, an high expression of developmental pluripotency-associated protein 3 (Dppa3) characterizes the cluster 3, correspondent to the preimplantation E3.5 stage, while the cluster 4, correspondent to the post-implantation stage E4.5 express high levels of laminin (Lama1 and Lamc1) required for the attachment, migration and organization of cells into tissues.

Following the gereral characterization of the dataset, already described by the authors [@mohammed_single-cell_2017], specific SingleCellExperiment files have been generated from the abundance output for the stage E3.5 and E4.5, and the same Bioconductor pipeline [@lun_simplesinglecell_2018] has been performed on both datasets, in order to detect highly variable genes, significantly correlated genes and subpopulation-specific marker genes characteristic of each stage.

The approach used for dimensionality reduction is the t-stochastic neighbour embedding (tSNE) method [@maaten_visualizing_2008], which showed to be more suitable for capturing the non-linear relationships in the high dimensional space and separating subpopulations by expressing features.
A hierarchical clustering on the Euclidean distances between cells was performed, generating a dendrogram that groups together cells with similar expression patterns. Then, clusters have been explicitly defined by applying a dynamic tree cut [@langfelder_wgcna:_2008] to the dendrogram.

![](images/e35_tSNE.png)

![](images/e45_tSNE.png)

*Fig. 29 - tSNE plot of the denoised PCA of the E3.5 dataset (perplexity = 20).*


The separatedness of the clusters has been verified using the silhouette width (Fig. 17, Fig. 29) for each cluster marker genes have been identified using limma [@ritchie_limma_2015].
Five clusters have been identified on E3.5 but the separatedness is much lower in comparison to the 4 clusters of E4.5 (in particular cluster 2 and 3), so they are not likely to reflect actual populations on E3.5 dataset. The cluster 2 of E4.5 on the other hand, groups a subpopulation apart from the rest of the cells with different features.

For both datasets the clusters have been tested. The top differentially expressed genes are likely to be good candidate markers as they can effectively distinguish between cells in different clusters. The top markers candidates for cluster 1 of E3.5 (which contains the majority of cells and the largest positive silhouette width) and cluster 2 of E4.5, have been chosen to be visualized on a heatmap, together with the expression profile of some pluripotency markers (Oct4, Nanog, Esrrb, Otx2).

![30a](images/e35_top_markers_heatmap_clusters.png)


***


![30b](images/e45_heatmap_top-markers_clusters.png)

*Fig. 30 - Heatmap of mean-centred and normalized log-expression values for the top set of markers for cluster 1 in the E3.5 dataset. Markers of pluripotency and lineage differentiation are also added on top rows.*


While no pattern is recognizable on the heatmap for E3.5, the expression profiles for the markers for cluster 2 of E4.5 dataset result strongly associated to a subpopulation of cells expressing typical features of the epiblast phenotype.

This is in accordance to the biological knowledge on the embryo development. As reported in the original study, at stage E3.5 transcriptional variability occurring in the inner cell mass (ICM) in the absence of cell-type substructure indicates the existence of uncoordinated transcriptional heterogeneity or transcriptional noise [@mohammed_single-cell_2017].
But this also confirms previous findings that genes specific of PrE (Gata6, Sox17) and epiblast (Nanog, Sox2) initially result co-expressed in cells within the E3.5 ICM before becoming restricted to specific cell types by E4.5 [@ohnishi_cell_2014]. 

![31a](images/e35_tSNE_markers.png)


***


![31b](images/e45_tSNE_markers.png)

*Fig. 31 - t-SNE plot comparison between E3.5 and E4.5 datasets showing the expression distribution of epiblast and primitive endoderm markers across the cells.*


It is also known by studies exploring the cooperation of transcription factors in PrE induction in ESCs, that Sox17/Oct4 partner to co-select specific target genes during commitment to primitive endoderm [@aksoy_oct4_2013].
A model has been proposed. In pluripotent cells, Oct4 and Sox2 expression levels being high, both factors cooperate and target speci???cally the canonical motif to regulate the expression of genes involved in self-renewal and pluripotency (e.g., Nanog). When these cells are subjected to an endodermal differentiation signal such as FGF4 within the ICM, Sox17 levels increase leading to a switch of Oct4 from an interaction with Sox2 to an interaction with Sox17, and thereby targets speci???c genes containing a compressed motif to trigger the endodermal expression program.
The expression profile of these genes across the stages E3.5 and E4.5 in this dataset seems to confirm this model (Fig. 31), supporting the hypotesis that the cluster 2 identified by tSNE at stage E4.5 correspond to the epiblast cell type.

It's also important to notice that Esrrb is included among the top markers of cluster 2 showing an expression profile very similar to Nanog and others pluripotency genes on the dataset. Among the Nanog target genes identified by using a Nanog knockout system on mouse ESCs, Esrrb showed the strongest transcriptional induction. Moreover the existance of partial funcional overlap between Nanog and Esrrb has been proposed [@festuccia_esrrb_2012].

On the next step the attention has been focused on the genes that are driving heterogeneity across the population of cells, by identifying highly variable genes (HVGs): those genes with the largest biological components of the variance in expression.
In order to distinguish between HVGs caused by random noise and those involved in driving systematic differences between subpopulations, gene pairs with significantly positive or negative correlation have been quantified by computing Spearman's rho(rho = 1 means a perfect positive correlation and the value rho = -1 means a perfect negataive correlation), which accommodates non-linear relationships in the expression values, therefore has been preferred to the classical Pearson's correlation coefficients.
On the E3.5 dataset only 12 out of 2601 significantly correlated genes pairs (FDR <= 5%) include Nanog, while 341 out of 119754 genes resulted significantly correlated to Nanog in the E4.5 dataset.

![32a](images/e35_heatmap_sig_cor_nanog.png)


***


![32b](images/e45.mRNA_heatmap_sig_cor_nanog.png)

*Fig. 32 - Heatmap of mean-centred and normalized log-expression values for the top HVGs correlated to Nanog in the E3.5 and E4.5 datasets *


Among the list of genes positively correlated to Nanog at stage E3.5 it's important to report Rad51 (rho = 0.42), a key component of the DNA double-strand breaks repair through homologous recombination and has been shown to increase the reprogramming efficiency during the generation of induced pluripotent stem cells by facilitating mesenchymal-to-epithelial transition during the early phase of the reprogramming process (@lee_rad51_2016), and Fgf4 (rho = 0.39), which contributes to the ICM lineage segregation [@ohnishi_cell_2014].
Notably, Sox17 and Gata6 results negatively correlated to Nanog (rho = -0.41 / -0.36), anticipating the emergence of these two TFs among the PrE precursors at stage E4.5. 

The same correletion pattern can in fact been found on the list of Nanog correlatet genes at stage E4.5, with Gata6 showing a stronger negative association (rho = -0.54) on conjunction with specific PrE phenotype markers like Dab2 (rho = -0.71) and Sox17 (rho = -0.44).
Observing the heatmap on figure 32b it is possible to spot a group of genes with high positive correlation to Nanog overlapping with the top marker genes used to define the cluster 2 on the previous step (Tdgf1, Morc1, Slc7a3, Igfbp2) attributable to the epiblast subpopulation.
The presence of Esrrb within the Nanog correlated genes (rho = 0.57) and its coexpression with Nanog in the epiblast subpopulation confirms the functional overlap observed in the ESCs model [@festuccia_esrrb_2012].

It's also curious to notice at this stage, a positive correlation between Nanog and Otx2 (rho = 0.49) showing a similar expression profile characterized by high levels in the epiblast subpopulation (see Fig. 30b). Nanog is a direct target of Otx2 and previous studies suggested how Otx2 regulation of Nanog contributes to ICM differentiation of the epiblast, observing that stage E4.5 is characterized by contemporary and frequently complementary expression of Nanog and Otx2, resembling the identity of ESCs cultured in LIF/FBS [@acampora_loss_2016], where through mutual antagonism, these TFs specify the heterogeneous identity of ESCs predisposing them for optimal response to naive or primed inducing factors [@acampora_functional_2017].

More than 300 genes resulted significantly correlated to Nanog at stage E4.5. In order to understand which of these genes is driving the specification of the epiblast, the downstream analyais has been focused only on the cluster 2, subsetting the dataset to include only cells showing elevated expression for Oct4 and others naive pluripotency markers. However the Searman's rank correlation performed on this subset produced correlated a list of gene pair with a low significance (FDR greater than 35%) which cannot lead to strong biological conclusions. 

At the end of the study a transcription level comparison has been berformed, without summarizing the expression at gene level, in order to get some insights into the distribution of expression profile among the trascripts of interesting genes.

![33a](images/e35_mRNA_heatmap_top_marker.png)


***


![33b](images/e45.mRNA_top_markers_heatmap_clusters.png)

*Fig. 33 - Heatmap of mean-centred and normalized log-expression values for the top set of markers transcripts for E3.5 and E4.5 dataset. Marker transcripts of pluripotency and lineage differentiation are also added on top rows.*

Surprisingly, many alternative spliced variants for the naive pluripotency genes can be observed with this approach, as well as their differential expression profile. 
Six protein-coding transcripts were detected for Oct4 but only one able of biological function (ENSMUST00000025271.16).

Nanog includes 3 protein-coding transcripts. One of them is results highly expressed at stage E3.5, the other two variants have been characterized to show attenuated capacities for self-renewal and pluripotency in ESCs [@das_alternative_2011]. However two Nanog isoforms (ENSMUST00000012540.4 and ENSMUST00000012540.4) curiously showed a complementary expression profile at stage E4.5 for this dataset (Fig. 33b) which is not considered during the gene level analysis.

The mouse Esrrb gene has six coding exons, with evidence for six alternatively spliced Esrrb mRNAs in the ENSEMBL EST databases, but one of them is not protein coding (ENSMUST00000136464.2). the abundance for all the six isoforms have been  quantified by kallisto and the most abundant transcripts on the E4.5 stage is Esrrb-203 (ENSMUST00000110204.8, total logcounts = 170.7892) followed by Esrrb-206 (ENSMUST00000167891.1, total logcounts = 104.4746). 

To understand the biological relationship among these alternative transcripts, the Spearman's rho has been computed to find correlate pairs on the E4.5 dataset. 

```{r echo = FALSE, results = 'asis'}
library(knitr)
kable(pluri.cor[sig.pluri.cor,], caption = "Nanog vs Esrrb mRNA correlation")
```
*Fig. 34 - Spearman's rank correlate pairs computed for Nanog and Esrrb variants within the E4.5 dataset.*



# 5. Conclusions

On multiple studies ChIP-seq has been used to map the locations of speci???c TFs considerred crucial for the transcriptional regulatory networks in embryonic stem cells. These factors are known to play different roles in ES-cell biology as components of the LIF
and BMP signaling pathways, self-renewal regulators, and key reprogramming factors and a combinatorial control of transcription factors has been observer [@chen_integration_2008, @kim_extended_2008]. Examining the binding pro???les, has been found that a subset of binding sites was bound by many of these TFs, and by clustered the peak sites was possible to define multiple transcription factor-binding loci (MTL).
In particular a Nanog-Oct4-Sox2-speci???c MTL has been characterized, and exhibits features of enhanceosomes by enhancing transcription from a distance.

Other studies analyzed the effect of acute depletion of Oct4 or Nanog, in order to understand the mechanisms underlying the transcriptional modulation of the naive pluripotency of ESCs [@hall_oct4_2009],[@festuccia_esrrb_2012].

Less studies however have been carried out on the embryo.

The list of Nanog correlated genes identified in this single-cell dataset for the stages E3.5 and E4.5 of the mouse embryo development, will be used for a comparison study to a set of 64 genes responding to Nanog identified using an ESCs Nanog -/- model [@festuccia_esrrb_2012].
The Spearman's rank correlation coefficient assigned to each gene will be used to discover the strength of a link between two sets of data.
As a control for this comparison the expression variation of Esrrb in correlation with Nanog will be used as reference, being Esrrb reported experimentally to have the strongest transcriptional induction as a result of Nanog binding [@festuccia_esrrb_2012].

If the in-vitro ESCs Nanog -/- model will be proven to be consistent with the embryo expression analysis reported here, this study will hopefully lead to a better understanding of the pluripotency gene regulatory network.


# 6. Bibliography


@article{das_alternative_2011,
	title = {Alternative Splicing Produces Nanog Protein Variants with Different Capacities for Self-renewal and Pluripotency in Embryonic Stem Cells},
	volume = {286},
	issn = {0021-9258, 1083-351X},
	url = {http://www.jbc.org/content/286/49/42690},
	doi = {10.1074/jbc.M111.290189},
	abstract = {Embryonic stem ({ES}) cells are distinguished by their ability to undergo unlimited self-renewal although retaining pluripotency, the capacity to specify cells of all germ layers. Alternative splicing contributes to these biological processes by vastly increasing the protein coding repertoire, enabling genes to code for novel variants that may confer different biological functions. The homeodomain transcription factor Nanog acts collaboratively with core factors Oct4 and Sox2 to govern the maintenance of pluripotency. We have discovered that Nanog is regulated by alternative splicing. Two novel exons and six subexons have been identified that extend the known Nanog gene structure and protein coding capacity. Alternative splicing results in two novel Nanog protein variants with attenuated capacities for self-renewal and pluripotency in {ES} cells. Our previous results have implicated the C-terminal domain, including the tryptophan-rich ({WR}) domain of Nanog, to be important for the function of Nanog (Wang, J., Levasseur, D. N., and Orkin, S. H. (2008) Proc. Natl. Acad. Sci. U.S.A. 105, 6326–6331). Using point mutation analyses, serine 2 (Ser-2) of Nanog has been identified as critical for {ES} cell self-renewal and for stabilizing a pluripotent gene signature. An inducible conditional knock-out was created to test the ability of new Nanog variants to genetically complement Nanog null {ES} cells. These results reveal for the first time an expanded Nanog protein coding capacity. We further reveal that a short region of the N-terminal domain and a single phosphorylatable Ser-2 is essential for the maintenance of self-renewal and pluripotency, demonstrating that this region of the protein is highly regulated.},
	pages = {42690--42703},
	number = {49},
	journaltitle = {Journal of Biological Chemistry},
	shortjournal = {J. Biol. Chem.},
	author = {Das, Satyabrata and Jena, Snehalata and Levasseur, Dana N.},
	urldate = {2018-09-06},
	date = {2011-12-09},
	langid = {english},
	pmid = {21969378}
}

@article{acampora_loss_2016,
	title = {Loss of the Otx2-Binding Site in the Nanog Promoter Affects the Integrity of Embryonic Stem Cell Subtypes and Specification of Inner Cell Mass-Derived Epiblast},
	volume = {15},
	issn = {2211-1247},
	url = {https://www.cell.com/cell-reports/abstract/S2211-1247(16)30624-6},
	doi = {10.1016/j.celrep.2016.05.041},
	abstract = {{\textless}h2{\textgreater}Summary{\textless}/h2{\textgreater}{\textless}p{\textgreater}Mouse embryonic stem cells ({ESCs}) and the inner cell mass ({ICM})-derived epiblast exhibit naive pluripotency. {ESC}-derived epiblast stem cells ({EpiSCs}) and the postimplantation epiblast exhibit primed pluripotency. Although core pluripotency factors are well-characterized, additional regulators, including Otx2, recently have been shown to function during the transition from naive to primed pluripotency. Here we uncover a role for Otx2 in the control of the naive pluripotent state. We analyzed Otx2-binding activity in {ESCs} and {EpiSCs} and identified Nanog, Oct4, and Sox2 as direct targets. To unravel the Otx2 transcriptional network, we targeted the strongest Otx2-binding site in the Nanog promoter, finding that this site modulates the size of specific {ESC}-subtype compartments in cultured cells and promotes Nanog expression in vivo, predisposing {ICM} differentiation to epiblast. Otx2-mediated Nanog regulation thus contributes to the integrity of the {ESC} state and cell lineage specification in preimplantation development.{\textless}/p{\textgreater}},
	pages = {2651--2664},
	number = {12},
	journaltitle = {Cell Reports},
	shortjournal = {Cell Reports},
	author = {Acampora, Dario and Omodei, Daniela and Petrosino, Giuseppe and Garofalo, Arcomaria and Savarese, Marco and Nigro, Vincenzo and Di Giovannantonio, Luca Giovanni and Mercadante, Vincenzo and Simeone, Antonio},
	urldate = {2018-09-06},
	date = {2016-06-21},
	pmid = {27292645, 27292645},
	file = {Acampora et al. - 2016 - Loss of the Otx2-Binding Site in the Nanog Promote.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\J3ENHN6Z\\Acampora et al. - 2016 - Loss of the Otx2-Binding Site in the Nanog Promote.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\TL3BUA49\\S2211-1247(16)30624-6.html:text/html}
}

@article{lee_rad51_2016,
	title = {Rad51 Regulates Reprogramming Efficiency through {DNA} Repair Pathway},
	volume = {20},
	issn = {2465-9525},
	url = {https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5027222/},
	doi = {10.12717/DR.2016.20.2.163},
	abstract = {Rad51 is a key component of homologous recombination ({HR}) to repair {DNA} double-strand breaks and it forms Rad51 recombinase filaments of broken single-stranded {DNA} to promote {HR}. In addition to its role in {DNA} repair and cell cycle progression, Rad51 contributes to the reprogramming process during the generation of induced pluripotent stem cells. In light of this, we performed reprogramming experiments to examine the effect of co-expression of Rad51 and four reprogramming factors, Oct4, Sox2, Klf4, and c-Myc, on the reprogramming efficiency. Co-expression of Rad51 significantly increased the numbers of alkaline phosphatase-positive colonies and embryonic stem cell-like colonies during the process of reprogramming. Co-expression {ofRad}51 significantly increased the expression of epithelial markers at an early stage of reprogramming compared with control cells. Phosphorylated histone H2AX (γH2AX), which initiates the {DNA} double-strand break repair system, was highly accumulated in reprogramming intermediates upon co-expression of Rad51. This study identified a novel role of Rad51 in enhancing the reprogramming efficiency, possibly by facilitating mesenchymal-to-epithelial transition and by regulating a {DNA} damage repair pathway during the early phase of the reprogramming process.},
	pages = {163--169},
	number = {2},
	journaltitle = {Development \& Reproduction},
	shortjournal = {Dev Reprod},
	author = {Lee, Jae-Young and Kim, Dae-Kwan and Ko, Jeong-Jae and Kim, Keun Pil and Park, Kyung-Soon},
	urldate = {2018-09-06},
	date = {2016-06},
	pmid = {27660832},
	pmcid = {PMC5027222},
	file = {Lee et al. - 2016 - Rad51 Regulates Reprogramming Efficiency through D.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\P8UD4RB6\\Lee et al. - 2016 - Rad51 Regulates Reprogramming Efficiency through D.pdf:application/pdf}
}

@article{acampora_functional_2017,
	title = {Functional Antagonism between {OTX}2 and {NANOG} Specifies a Spectrum of Heterogeneous Identities in Embryonic Stem Cells},
	volume = {9},
	issn = {2213-6711},
	url = {https://www.cell.com/stem-cell-reports/abstract/S2213-6711(17)30426-5},
	doi = {10.1016/j.stemcr.2017.09.019},
	abstract = {{\textless}h2{\textgreater}Summary{\textless}/h2{\textgreater}{\textless}p{\textgreater}Embryonic stem cells ({ESCs}) cultured in leukemia inhibitory factor ({LIF}) plus fetal bovine serum ({FBS}) exhibit heterogeneity in the expression of naive and primed transcription factors. This heterogeneity reflects the dynamic condition of {ESCs} and their versatility to promptly respond to signaling effectors promoting naive or primed pluripotency. Here, we report that {ESCs} lacking \textit{Nanog} or overexpressing {\textless}i{\textgreater}Otx2{\textless}/i{\textgreater} exhibit an early primed identity in {LIF} + {FBS} and fail to convert into 2i-induced naive state. Conversely, {\textless}i{\textgreater}Otx2{\textless}/i{\textgreater}-null {ESCs} possess naive identity features in {LIF} + {FBS} similar to {\textless}i{\textgreater}Nanog{\textless}/i{\textgreater}-overexpressing {ESCs} and convert poorly into {FGF}-induced early primed state. When both {\textless}i{\textgreater}Nanog{\textless}/i{\textgreater} and {\textless}i{\textgreater}Otx2{\textless}/i{\textgreater} are inactivated, {ESCs} cultured in {LIF} + {FBS} exhibit primed identity and weakened ability to convert into naive state. These data suggest that, through mutual antagonism, {NANOG} and {OTX}2 specify the heterogeneous identity of {ESCs} cultured in {LIF} + {FBS} and individually predispose them for optimal response to naive or primed inducing factors.{\textless}/p{\textgreater}},
	pages = {1642--1659},
	number = {5},
	journaltitle = {Stem Cell Reports},
	shortjournal = {Stem Cell Reports},
	author = {Acampora, Dario and Giovannantonio, Luca Giovanni Di and Garofalo, Arcomaria and Nigro, Vincenzo and Omodei, Daniela and Lombardi, Alessia and Zhang, Jingchao and Chambers, Ian and Simeone, Antonio},
	urldate = {2018-09-04},
	date = {2017-11-14},
	pmid = {29056334, 29056334},
	file = {Acampora et al. - 2017 - Functional Antagonism between OTX2 and NANOG Speci.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\U366JDWD\\Acampora et al. - 2017 - Functional Antagonism between OTX2 and NANOG Speci.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\I5KPVK3S\\S2213-6711(17)30426-5.html:text/html}
}

@article{langfelder_wgcna:_2008,
	title = {{WGCNA}: an R package for weighted correlation network analysis},
	volume = {9},
	issn = {1471-2105},
	url = {https://doi.org/10.1186/1471-2105-9-559},
	doi = {10.1186/1471-2105-9-559},
	shorttitle = {{WGCNA}},
	abstract = {Correlation networks are increasingly being used in bioinformatics applications. For example, weighted gene co-expression network analysis is a systems biology method for describing the correlation patterns among genes across microarray samples. Weighted correlation network analysis ({WGCNA}) can be used for finding clusters (modules) of highly correlated genes, for summarizing such clusters using the module eigengene or an intramodular hub gene, for relating modules to one another and to external sample traits (using eigengene network methodology), and for calculating module membership measures. Correlation networks facilitate network based gene screening methods that can be used to identify candidate biomarkers or therapeutic targets. These methods have been successfully applied in various biological contexts, e.g. cancer, mouse genetics, yeast genetics, and analysis of brain imaging data. While parts of the correlation network methodology have been described in separate publications, there is a need to provide a user-friendly, comprehensive, and consistent software implementation and an accompanying tutorial.},
	pages = {559},
	number = {1},
	journaltitle = {{BMC} Bioinformatics},
	shortjournal = {{BMC} Bioinformatics},
	author = {Langfelder, Peter and Horvath, Steve},
	urldate = {2018-09-01},
	date = {2008-12-29},
	file = {Full Text PDF:C\:\\Users\\Ivan\\Zotero\\storage\\UGZIX6D5\\Langfelder and Horvath - 2008 - WGCNA an R package for weighted correlation netwo.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\2NUB3VXP\\1471-2105-9-559.html:text/html}
}

@article{maaten_visualizing_2008,
	title = {Visualizing Data using t-{SNE}},
	volume = {9},
	issn = {{ISSN} 1533-7928},
	url = {http://www.jmlr.org/papers/v9/vandermaaten08a.html},
	pages = {2579--2605},
	issue = {Nov},
	journaltitle = {Journal of Machine Learning Research},
	author = {Maaten, Laurens van der and Hinton, Geoffrey},
	urldate = {2018-09-01},
	date = {2008},
	file = {Full Text PDF:C\:\\Users\\Ivan\\Zotero\\storage\\E2W3JTKA\\Maaten and Hinton - 2008 - Visualizing Data using t-SNE.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\UHJZH4KV\\vandermaaten08a.html:text/html}
}

@article{everaert_benchmarking_2017,
	title = {Benchmarking of {RNA}-sequencing analysis workflows using whole-transcriptome {RT}-{qPCR} expression data},
	volume = {7},
	rights = {2017 The Author(s)},
	issn = {2045-2322},
	url = {https://www.nature.com/articles/s41598-017-01617-3},
	doi = {10.1038/s41598-017-01617-3},
	abstract = {{RNA}-sequencing has become the gold standard for whole-transcriptome gene expression quantification. Multiple algorithms have been developed to derive gene counts from sequencing reads. While a number of benchmarking studies have been conducted, the question remains how individual methods perform at accurately quantifying gene expression levels from {RNA}-sequencing reads. We performed an independent benchmarking study using {RNA}-sequencing data from the well established {MAQCA} and {MAQCB} reference samples. {RNA}-sequencing reads were processed using five workflows (Tophat-{HTSeq}, Tophat-Cufflinks, {STAR}-{HTSeq}, Kallisto and Salmon) and resulting gene expression measurements were compared to expression data generated by wet-lab validated {qPCR} assays for all protein coding genes. All methods showed high gene expression correlations with {qPCR} data. When comparing gene expression fold changes between {MAQCA} and {MAQCB} samples, about 85\% of the genes showed consistent results between {RNA}-sequencing and {qPCR} data. Of note, each method revealed a small but specific gene set with inconsistent expression measurements. A significant proportion of these method-specific inconsistent genes were reproducibly identified in independent datasets. These genes were typically smaller, had fewer exons, and were lower expressed compared to genes with consistent expression measurements. We propose that careful validation is warranted when evaluating {RNA}-seq based expression profiles for this specific gene set.},
	pages = {1559},
	number = {1},
	journaltitle = {Scientific Reports},
	author = {Everaert, Celine and Luypaert, Manuel and Maag, Jesper L. V. and Cheng, Quek Xiu and Dinger, Marcel E. and Hellemans, Jan and Mestdagh, Pieter},
	urldate = {2018-09-01},
	date = {2017-05-08},
	langid = {english},
	file = {Everaert et al. - 2017 - Benchmarking of RNA-sequencing analysis workflows .pdf:C\:\\Users\\Ivan\\Zotero\\storage\\8ASXXYGG\\Everaert et al. - 2017 - Benchmarking of RNA-sequencing analysis workflows .pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\W4EC5ASA\\s41598-017-01617-3.html:text/html}
}

@article{ohnishi_cell_2014,
	title = {Cell-to-cell expression variability followed by signal reinforcement progressively segregates early mouse lineages},
	volume = {16},
	rights = {2013 Nature Publishing Group},
	issn = {1476-4679},
	url = {http://www.nature.com/articles/ncb2881},
	doi = {10.1038/ncb2881},
	abstract = {It is now recognized that extensive expression heterogeneities among cells precede the emergence of lineages in the early mammalian embryo. To establish a map of pluripotent epiblast ({EPI}) versus primitive endoderm ({PrE}) lineage segregation within the inner cell mass ({ICM}) of the mouse blastocyst, we characterized the gene expression profiles of individual {ICM} cells. Clustering analysis of the transcriptomes of 66 cells demonstrated that initially they are non-distinguishable. Early in the segregation, lineage-specific marker expression exhibited no apparent correlation, and a hierarchical relationship was established only in the late blastocyst. Fgf4 exhibited a bimodal expression at the earliest stage analysed, and in its absence, the differentiation of {PrE} and {EPI} was halted, indicating that Fgf4 drives, and is required for, {ICM} lineage segregation. These data lead us to propose a model where stochastic cell-to-cell expression heterogeneity followed by signal reinforcement underlies {ICM} lineage segregation by antagonistically separating equivalent cells.},
	pages = {27--37},
	number = {1},
	journaltitle = {Nature Cell Biology},
	author = {Ohnishi, Yusuke and Huber, Wolfgang and Tsumura, Akiko and Kang, Minjung and Xenopoulos, Panagiotis and Kurimoto, Kazuki and Oleś, Andrzej K. and Araúzo-Bravo, Marcos J. and Saitou, Mitinori and Hadjantonakis, Anna-Katerina and Hiiragi, Takashi},
	urldate = {2018-08-30},
	date = {2014-01},
	langid = {english},
	file = {Ohnishi et al. - 2014 - Cell-to-cell expression variability followed by si.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\UPKLKC6V\\Ohnishi et al. - 2014 - Cell-to-cell expression variability followed by si.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\FXUGARC6\\ncb2881.html:text/html}
}

@online{r_core_team_r:_2018,
	title = {R: A Language and Environment for Statistical Computing},
	url = {http://www.R-project.org},
	author = {{R Core Team}},
	date = {2018-07-02}
}

@article{huber_orchestrating_2015,
	title = {Orchestrating high-throughput genomic analysis with Bioconductor},
	volume = {12},
	rights = {2015 Nature Publishing Group},
	issn = {1548-7105},
	url = {https://www.nature.com/articles/nmeth.3252},
	doi = {10.1038/nmeth.3252},
	abstract = {Bioconductor is an open-source, open-development software project for the analysis and comprehension of high-throughput data in genomics and molecular biology. The project aims to enable interdisciplinary research, collaboration and rapid development of scientific software. Based on the statistical programming language R, Bioconductor comprises 934 interoperable packages contributed by a large, diverse community of scientists. Packages cover a range of bioinformatic and statistical applications. They undergo formal initial review and continuous automated testing. We present an overview for prospective users and contributors.},
	pages = {115--121},
	number = {2},
	journaltitle = {Nature Methods},
	author = {Huber, Wolfgang and Carey, Vincent J. and Gentleman, Robert and Anders, Simon and Carlson, Marc and Carvalho, Benilton S. and Bravo, Hector Corrada and Davis, Sean and Gatto, Laurent and Girke, Thomas and Gottardo, Raphael and Hahne, Florian and Hansen, Kasper D. and Irizarry, Rafael A. and Lawrence, Michael and Love, Michael I. and {MacDonald}, James and Obenchain, Valerie and Oleś, Andrzej K. and Pagès, Hervé and Reyes, Alejandro and Shannon, Paul and Smyth, Gordon K. and Tenenbaum, Dan and Waldron, Levi and Morgan, Martin},
	urldate = {2018-06-12},
	date = {2015-02},
	langid = {english},
	file = {Full Text PDF:C\:\\Users\\Ivan\\Zotero\\storage\\3SZNLEBD\\Huber et al. - 2015 - Orchestrating high-throughput genomic analysis wit.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\YZDGQ6IZ\\nmeth.html:text/html}
}

@article{mccarthy_scater:_2017,
	title = {Scater: pre-processing, quality control, normalization and visualization of single-cell {RNA}-seq data in R},
	volume = {33},
	issn = {1367-4803},
	url = {https://academic.oup.com/bioinformatics/article/33/8/1179/2907823},
	doi = {10.1093/bioinformatics/btw777},
	shorttitle = {Scater},
	abstract = {Abstract.  Motivation: Single-cell {RNA} sequencing ({scRNA}-seq) is increasingly used to study gene expression at the level of individual cells. However, preparing},
	pages = {1179--1186},
	number = {8},
	journaltitle = {Bioinformatics},
	shortjournal = {Bioinformatics},
	author = {{McCarthy}, Davis J. and Campbell, Kieran R. and Lun, Aaron T. L. and Wills, Quin F.},
	urldate = {2018-07-10},
	date = {2017-04-15},
	langid = {english},
	file = {McCarthy et al. - 2017 - Scater pre-processing, quality control, normaliza.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\BQNUWZHB\\McCarthy et al. - 2017 - Scater pre-processing, quality control, normaliza.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\8BIDNQ7J\\2907823.html:text/html}
}

@article{scialdone_computational_2015,
	title = {Computational assignment of cell-cycle stage from single-cell transcriptome data},
	volume = {85},
	issn = {1046-2023},
	url = {http://www.sciencedirect.com/science/article/pii/S1046202315300098},
	doi = {10.1016/j.ymeth.2015.06.021},
	series = {Inferring Gene Regulatory Interactions from Quantitative High-Throughput Measurements},
	abstract = {The transcriptome of single cells can reveal important information about cellular states and heterogeneity within populations of cells. Recently, single-cell {RNA}-sequencing has facilitated expression profiling of large numbers of single cells in parallel. To fully exploit these data, it is critical that suitable computational approaches are developed. One key challenge, especially pertinent when considering dividing populations of cells, is to understand the cell-cycle stage of each captured cell. Here we describe and compare five established supervised machine learning methods and a custom-built predictor for allocating cells to their cell-cycle stage on the basis of their transcriptome. In particular, we assess the impact of different normalisation strategies and the usage of prior knowledge on the predictive power of the classifiers. We tested the methods on previously published datasets and found that a {PCA}-based approach and the custom predictor performed best. Moreover, our analysis shows that the performance depends strongly on normalisation and the usage of prior knowledge. Only by leveraging prior knowledge in form of cell-cycle annotated genes and by preprocessing the data using a rank-based normalisation, is it possible to robustly capture the transcriptional cell-cycle signature across different cell types, organisms and experimental protocols.},
	pages = {54--61},
	journaltitle = {Methods},
	shortjournal = {Methods},
	author = {Scialdone, Antonio and Natarajan, Kedar N. and Saraiva, Luis R. and Proserpio, Valentina and Teichmann, Sarah A. and Stegle, Oliver and Marioni, John C. and Buettner, Florian},
	urldate = {2018-08-29},
	date = {2015-09-01},
	file = {Scialdone et al. - 2015 - Computational assignment of cell-cycle stage from .pdf:C\:\\Users\\Ivan\\Zotero\\storage\\N7RL5BYH\\Scialdone et al. - 2015 - Computational assignment of cell-cycle stage from .pdf:application/pdf;ScienceDirect Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\8H6JIWLG\\S1046202315300098.html:text/html}
}

@article{ilicic_classification_2016,
	title = {Classification of low quality cells from single-cell {RNA}-seq data},
	volume = {17},
	issn = {1474-760X},
	url = {https://doi.org/10.1186/s13059-016-0888-1},
	doi = {10.1186/s13059-016-0888-1},
	abstract = {Single-cell {RNA} sequencing ({scRNA}-seq) has broad applications across biomedical research. One of the key challenges is to ensure that only single, live cells are included in downstream analysis, as the inclusion of compromised cells inevitably affects data interpretation. Here, we present a generic approach for processing {scRNA}-seq data and detecting low quality cells, using a curated set of over 20 biological and technical features. Our approach improves classification accuracy by over 30 \% compared to traditional methods when tested on over 5,000 cells, including {CD}4+ T cells, bone marrow dendritic cells, and mouse embryonic stem cells.},
	pages = {29},
	number = {1},
	journaltitle = {Genome Biology},
	shortjournal = {Genome Biology},
	author = {Ilicic, Tomislav and Kim, Jong Kyoung and Kolodziejczyk, Aleksandra A. and Bagger, Frederik Otzen and {McCarthy}, Davis James and Marioni, John C. and Teichmann, Sarah A.},
	urldate = {2018-08-29},
	date = {2016-02-17},
	file = {Ilicic et al. - 2016 - Classification of low quality cells from single-ce.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\LSPIBB2B\\Ilicic et al. - 2016 - Classification of low quality cells from single-ce.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\SX7QUUBY\\s13059-016-0888-1.html:text/html}
}

@article{islam_quantitative_2014,
	title = {Quantitative single-cell {RNA}-seq with unique molecular identifiers},
	volume = {11},
	rights = {2013 Nature Publishing Group},
	issn = {1548-7105},
	url = {https://www.nature.com/articles/nmeth.2772},
	doi = {10.1038/nmeth.2772},
	abstract = {Single-cell {RNA} sequencing ({RNA}-seq) is a powerful tool to reveal cellular heterogeneity, discover new cell types and characterize tumor microevolution. However, losses in {cDNA} synthesis and bias in {cDNA} amplification lead to severe quantitative errors. We show that molecular labels—random sequences that label individual molecules—can nearly eliminate amplification noise, and that microfluidic sample preparation and optimized reagents produce a fivefold improvement in {mRNA} capture efficiency.},
	pages = {163--166},
	number = {2},
	journaltitle = {Nature Methods},
	author = {Islam, Saiful and Zeisel, Amit and Joost, Simon and Manno, Gioele La and Zajac, Pawel and Kasper, Maria and Lönnerberg, Peter and Linnarsson, Sten},
	urldate = {2018-08-29},
	date = {2014-02},
	langid = {english},
	file = {Islam et al. - 2014 - Quantitative single-cell RNA-seq with unique molec.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\G7QQY2ZC\\Islam et al. - 2014 - Quantitative single-cell RNA-seq with unique molec.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\DQDNTX79\\nmeth.html:text/html}
}

@article{oliveri_global_2008,
	title = {Global regulatory logic for specification of an embryonic cell lineage},
	volume = {105},
	rights = {© 2008 by The National Academy of Sciences of the {USA}},
	issn = {0027-8424, 1091-6490},
	url = {http://www.pnas.org/content/105/16/5955},
	doi = {10.1073/pnas.0711220105},
	abstract = {Explanation of a process of development must ultimately be couched in the terms of the genomic regulatory code. Specification of an embryonic cell lineage is driven by a network of interactions among genes encoding transcription factors. Here, we present the gene regulatory network ({GRN}) that directs the specification of the skeletogenic micromere lineage of the sea urchin embryo. The {GRN} now includes all regulatory genes expressed in this lineage up to late blastula stage, as identified in a genomewide survey. The architecture of the {GRN} was established by a large-scale perturbation analysis in which the expression of each gene in the {GRN} was cut off by use of morpholinos, and the effects on all other genes were measured quantitatively. Several cis-regulatory analyses provided additional evidence. The explanatory power of the {GRN} suffices to provide a causal explanation for all observable developmental functions of the micromere lineage during the specification period. These functions are: (i) initial acquisition of identity through transcriptional interpretation of localized maternal cues; (ii) activation of specific regulatory genes by use of a double negative gate; (iii) dynamic stabilization of the regulatory state by activation of a feedback subcircuit; (iv) exclusion of alternative regulatory states; (v) presentation of a signal required by the micromeres themselves and of two different signals required for development of adjacent endomesodermal lineages; and (vi) lineage-specific activation of batteries of skeletogenic genes. The {GRN} precisely predicts gene expression responses and provides a coherent explanation of the biology of specification.},
	pages = {5955--5962},
	number = {16},
	journaltitle = {Proceedings of the National Academy of Sciences},
	shortjournal = {{PNAS}},
	author = {Oliveri, Paola and Tu, Qiang and Davidson, Eric H.},
	urldate = {2018-08-29},
	date = {2008-04-22},
	langid = {english},
	pmid = {18413610},
	file = {Oliveri et al. - 2008 - Global regulatory logic for specification of an em.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\C8G4I5X5\\Oliveri et al. - 2008 - Global regulatory logic for specification of an em.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\JIYV9A9J\\5955.html:text/html}
}

@article{chambers_functional_2003,
	title = {Functional Expression Cloning of Nanog, a Pluripotency Sustaining Factor in Embryonic Stem Cells},
	volume = {113},
	issn = {0092-8674, 1097-4172},
	url = {https://www.cell.com/cell/abstract/S0092-8674(03)00392-1},
	doi = {10.1016/S0092-8674(03)00392-1},
	pages = {643--655},
	number = {5},
	journaltitle = {Cell},
	shortjournal = {Cell},
	author = {Chambers, Ian and Colby, Douglas and Robertson, Morag and Nichols, Jennifer and Lee, Sonia and Tweedie, Susan and Smith, Austin},
	urldate = {2018-08-29},
	date = {2003-05-30},
	pmid = {12787505},
	file = {Chambers et al. - 2003 - Functional Expression Cloning of Nanog, a Pluripot.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\8IIL6LSR\\Chambers et al. - 2003 - Functional Expression Cloning of Nanog, a Pluripot.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\NLI2BRXH\\S0092-8674(03)00392-1.html:text/html}
}

@article{ewels_multiqc:_2016,
	title = {{MultiQC}: summarize analysis results for multiple tools and samples in a single report},
	volume = {32},
	issn = {1367-4803},
	url = {https://academic.oup.com/bioinformatics/article/32/19/3047/2196507},
	doi = {10.1093/bioinformatics/btw354},
	shorttitle = {{MultiQC}},
	abstract = {Abstract.  Motivation: Fast and accurate quality control is essential for studies involving next-generation sequencing data. Whilst numerous tools exist to quan},
	pages = {3047--3048},
	number = {19},
	journaltitle = {Bioinformatics},
	shortjournal = {Bioinformatics},
	author = {Ewels, Philip and Magnusson, Måns and Lundin, Sverker and Käller, Max},
	urldate = {2018-08-29},
	date = {2016-10-01},
	langid = {english},
	file = {Ewels et al. - 2016 - MultiQC summarize analysis results for multiple t.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\8SLYCYD4\\Ewels et al. - 2016 - MultiQC summarize analysis results for multiple t.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\LWKET2LY\\2196507.html:text/html}
}

@online{andrews_fastqc_2018,
	title = {{FastQC} - A Quality Control tool for High Throughput Sequence Data},
	url = {http://www.bioinformatics.babraham.ac.uk/projects/fastqc/},
	titleaddon = {Babraham Bioinformatics},
	author = {Andrews, Simon},
	urldate = {2018-04-24},
	date = {2018-01-10},
	file = {FastQC_Manual.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\YEYXFTLN\\FastQC_Manual.pdf:application/pdf}
}

@article{bray_near-optimal_2016,
	title = {Near-optimal probabilistic {RNA}-seq quantification},
	volume = {34},
	rights = {2016 Nature Publishing Group},
	issn = {1546-1696},
	url = {https://www.nature.com/articles/nbt.3519},
	doi = {10.1038/nbt.3519},
	abstract = {We present kallisto, an {RNA}-seq quantification program that is two orders of magnitude faster than previous approaches and achieves similar accuracy. Kallisto pseudoaligns reads to a reference, producing a list of transcripts that are compatible with each read while avoiding alignment of individual bases. We use kallisto to analyze 30 million unaligned paired-end {RNA}-seq reads in {\textless}10 min on a standard laptop computer. This removes a major computational bottleneck in {RNA}-seq analysis.},
	pages = {525--527},
	number = {5},
	journaltitle = {Nature Biotechnology},
	author = {Bray, Nicolas L. and Pimentel, Harold and Melsted, Páll and Pachter, Lior},
	urldate = {2018-06-23},
	date = {2016-05},
	langid = {english},
	file = {Full Text PDF:C\:\\Users\\Ivan\\Zotero\\storage\\JN7EYK6Z\\Bray et al. - 2016 - Near-optimal probabilistic RNA-seq quantification.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\5CYAP3HC\\nbt.html:text/html}
}

@article{mohammed_single-cell_2017,
	title = {Single-Cell Landscape of Transcriptional Heterogeneity and Cell Fate Decisions during Mouse Early Gastrulation},
	volume = {20},
	issn = {2211-1247},
	url = {http://www.sciencedirect.com/science/article/pii/S2211124717309610},
	doi = {10.1016/j.celrep.2017.07.009},
	abstract = {Summary
The mouse inner cell mass ({ICM}) segregates into the epiblast and primitive endoderm ({PrE}) lineages coincident with implantation of the embryo. The epiblast subsequently undergoes considerable expansion of cell numbers prior to gastrulation. To investigate underlying regulatory principles, we performed systematic single-cell {RNA} sequencing (seq) of conceptuses from E3.5 to E6.5. The epiblast shows reactivation and subsequent inactivation of the X chromosome, with Zfp57 expression associated with reactivation and inactivation together with other candidate regulators. At E6.5, the transition from epiblast to primitive streak is linked with decreased expression of polycomb subunits, suggesting a key regulatory role. Notably, our analyses suggest elevated transcriptional noise at E3.5 and within the non-committed epiblast at E6.5, coinciding with exit from pluripotency. By contrast, E6.5 primitive streak cells became highly synchronized and exhibit a shortened G1 cell-cycle phase, consistent with accelerated proliferation. Our study systematically charts transcriptional noise and uncovers molecular processes associated with early lineage decisions.},
	pages = {1215--1228},
	number = {5},
	journaltitle = {Cell Reports},
	shortjournal = {Cell Reports},
	author = {Mohammed, Hisham and Hernando-Herraez, Irene and Savino, Aurora and Scialdone, Antonio and Macaulay, Iain and Mulas, Carla and Chandra, Tamir and Voet, Thierry and Dean, Wendy and Nichols, Jennifer and Marioni, John C. and Reik, Wolf},
	urldate = {2018-03-20},
	date = {2017-08-01},
	keywords = {read, reading, to read},
	file = {1-s2.0-S2211124717309610-mmc2.xlsx:C\:\\Users\\Ivan\\Zotero\\storage\\L8LLTLTP\\1-s2.0-S2211124717309610-mmc2.xlsx:application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;1-s2.0-S2211124717309610-mmc3.xlsx:C\:\\Users\\Ivan\\Zotero\\storage\\EGVSXR7T\\1-s2.0-S2211124717309610-mmc3.xlsx:application/msword;1-s2.0-S2211124717309610-mmc4.xlsx:C\:\\Users\\Ivan\\Zotero\\storage\\DZ5VTMEL\\1-s2.0-S2211124717309610-mmc4.xlsx:application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;1-s2.0-S2211124717309610-mmc5.xlsx:C\:\\Users\\Ivan\\Zotero\\storage\\4Q9XNKEP\\1-s2.0-S2211124717309610-mmc5.xlsx:application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;1-s2.0-S2211124717309610-mmc6.xlsx:C\:\\Users\\Ivan\\Zotero\\storage\\QBEKKE27\\1-s2.0-S2211124717309610-mmc6.xlsx:application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;GSE100597_count_table_QC_filtered.txt:C\:\\Users\\Ivan\\Zotero\\storage\\S7LGT2KK\\GSE100597_count_table_QC_filtered.txt:text/plain;Mohammed et al. - 2017 - Single-Cell Landscape of Transcriptional Heterogen.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\QGS6UTJM\\Mohammed et al. - 2017 - Single-Cell Landscape of Transcriptional Heterogen.pdf:application/pdf;ScienceDirect Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\29JIKGBY\\S2211124717309610.html:text/html}
}

@article{ptashne_transcription_2001,
	title = {Transcription initiation: imposing specificity by localization},
	volume = {37},
	issn = {0071-1365},
	shorttitle = {Transcription initiation},
	abstract = {Many crucial cellular enzymes--including {RNA} polymerases, kinases, phosphatases, proteases, acetylaters, etc.--have multiple potential substrates. Regulation entails substrate selection, a process effected by a mechanism we call regulated localization. This formulation is particularly well illustrated by the mechanisms of gene regulation. Analysis of these mechanisms reveals that regulated localization requires simple molecular interactions. These molecular interactions readily lend themselves to combinatorial control. This system of regulation is highly 'evolvable'. Its use accounts, at least in part, for the nature of many of the complexities observed in biological systems.},
	pages = {1--15},
	journaltitle = {Essays in Biochemistry},
	shortjournal = {Essays Biochem.},
	author = {Ptashne, M. and Gann, A.},
	date = {2001},
	pmid = {11758451}
}

@article{boyer_core_2005,
	title = {Core Transcriptional Regulatory Circuitry in Human Embryonic Stem Cells},
	volume = {122},
	issn = {0092-8674, 1097-4172},
	url = {http://www.cell.com/cell/abstract/S0092-8674(05)00825-1},
	doi = {10.1016/j.cell.2005.08.020},
	pages = {947--956},
	number = {6},
	journaltitle = {Cell},
	shortjournal = {Cell},
	author = {Boyer, Laurie A. and Lee, Tong Ihn and Cole, Megan F. and Johnstone, Sarah E. and Levine, Stuart S. and Zucker, Jacob P. and Guenther, Matthew G. and Kumar, Roshan M. and Murray, Heather L. and Jenner, Richard G. and Gifford, David K. and Melton, Douglas A. and Jaenisch, Rudolf and Young, Richard A.},
	urldate = {2018-04-23},
	date = {2005-09-23},
	file = {Boyer et al. - 2005 - Core Transcriptional Regulatory Circuitry in Human.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\B2G4KGUB\\Boyer et al. - 2005 - Core Transcriptional Regulatory Circuitry in Human.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\YBMYTQ5N\\S0092-8674(05)00825-1.html:text/html}
}

@article{munoz_descalzo_silvia_correlations_2012,
	title = {Correlations Between the Levels of Oct4 and Nanog as a Signature for Naïve Pluripotency in Mouse Embryonic Stem Cells},
	volume = {30},
	issn = {1066-5099},
	url = {http://stemcellsjournals.onlinelibrary.wiley.com/doi/abs/10.1002/stem.1230},
	doi = {10.1002/stem.1230},
	abstract = {Abstract The pluripotent state is traditionally associated with large absolute levels of certain transcription factors such as Nanog and Oct4. Here, we present experimental observations using quantitative immunofluorescence that pluripotency in mouse embryonic stem cells ({mESCs}) is established by specific ratios between Oct4 and Nanog. When cells are grown in 2i conditions, they exhibit uniform levels of pluripotency and this is associated with a high correlation between the levels of Oct4 and Nanog in individual cells. The correlation is lost when cells differentiate. Our results suggest that the correlation between these two factors and the distribution of Oct4/Nanog ratios can be used as quantifiers to distinguish between three subpopulations in an {mESC} culture: pluripotent, lineage?primed, and differentiating cells. When we apply these quantifiers to cells with lower levels of Nanog or mutant for ??Catenin or Tcf3, the results suggest that these cells exhibit higher probability of differentiation. {STEM} {CELLS} 2012;30:2683?2691},
	pages = {2683--2691},
	number = {12},
	journaltitle = {{STEM} {CELLS}},
	shortjournal = {{STEM} {CELLS}},
	author = {{Muñoz Descalzo Silvia} and {RuÉ Pau} and {Garcia‐Ojalvo Jordi} and {Arias Alfonso Martinez}},
	urldate = {2018-04-23},
	date = {2012-11-27},
	file = {Muñoz Descalzo Silvia et al. - 2012 - Correlations Between the Levels of Oct4 and Nanog .pdf:C\:\\Users\\Ivan\\Zotero\\storage\\ILTLUQ5F\\Muñoz Descalzo Silvia et al. - 2012 - Correlations Between the Levels of Oct4 and Nanog .pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\G8CQLUWK\\stem.html:text/html}
}

@article{kalkan_tracking_2017,
	title = {Tracking the embryonic stem cell transition from ground state pluripotency},
	volume = {144},
	rights = {© 2017. Published by The Company of Biologists Ltd. http://creativecommons.org/licenses/by/3.0This is an Open Access article distributed under the terms of the Creative Commons Attribution License (http://creativecommons.org/licenses/by/3.0), which permits unrestricted use, distribution and reproduction in any medium provided that the original work is properly attributed.},
	issn = {0950-1991, 1477-9129},
	url = {http://dev.biologists.org/content/144/7/1221},
	doi = {10.1242/dev.142711},
	abstract = {Skip to Next Section
Mouse embryonic stem ({ES}) cells are locked into self-renewal by shielding from inductive cues. Release from this ground state in minimal conditions offers a system for delineating developmental progression from naïve pluripotency. Here, we examine the initial transition process. The {ES} cell population behaves asynchronously. We therefore exploited a short-half-life Rex1::{GFP} reporter to isolate cells either side of exit from naïve status. Extinction of {ES} cell identity in single cells is acute. It occurs only after near-complete elimination of naïve pluripotency factors, but precedes appearance of lineage specification markers. Cells newly departed from the {ES} cell state display features of early post-implantation epiblast and are distinct from primed epiblast. They also exhibit a genome-wide increase in {DNA} methylation, intermediate between early and late epiblast. These findings are consistent with the proposition that naïve cells transition to a distinct formative phase of pluripotency preparatory to lineage priming.},
	pages = {1221--1234},
	number = {7},
	journaltitle = {Development},
	shortjournal = {Development},
	author = {Kalkan, Tüzer and Olova, Nelly and Roode, Mila and Mulas, Carla and Lee, Heather J. and Nett, Isabelle and Marks, Hendrik and Walker, Rachael and Stunnenberg, Hendrik G. and Lilley, Kathryn S. and Nichols, Jennifer and Reik, Wolf and Bertone, Paul and Smith, Austin},
	urldate = {2018-04-23},
	date = {2017-04-01},
	langid = {english},
	pmid = {28174249},
	file = {Kalkan et al. - 2017 - Tracking the embryonic stem cell transition from g.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\IEAVU6QS\\Kalkan et al. - 2017 - Tracking the embryonic stem cell transition from g.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\MZ4CR573\\1221.html:text/html}
}

@article{festuccia_esrrb_2012,
	title = {Esrrb Is a Direct Nanog Target Gene that Can Substitute for Nanog Function in Pluripotent Cells},
	volume = {11},
	issn = {1934-5909},
	url = {http://www.sciencedirect.com/science/article/pii/S1934590912004808},
	doi = {10.1016/j.stem.2012.08.002},
	abstract = {Summary
Embryonic stem cell ({ESC}) self-renewal efficiency is determined by the level of Nanog expression. However, the mechanisms by which Nanog functions remain unclear, and in particular, direct Nanog target genes are uncharacterized. Here we investigate {ESCs} expressing different Nanog levels and Nanog−/− cells with distinct functionally inducible Nanog proteins to identify Nanog-responsive genes. Surprisingly, these constitute a minor fraction of genes that Nanog binds. Prominent among Nanog-reponsive genes is Estrogen-related receptor b (Esrrb). Nanog binds directly to Esrrb, enhances binding of {RNAPolII}, and stimulates Esrrb transcription. Overexpression of Esrrb in {ESCs} maintains cytokine-independent self-renewal and pluripotency. Remarkably, this activity is retained in Nanog−/− {ESCs}. Moreover, Esrrb can reprogram Nanog−/− {EpiSCs} and can rescue stalled reprogramming in Nanog−/− pre-{iPSCs}. Finally, Esrrb deletion abolishes the defining ability of Nanog to confer {LIF}-independent {ESC} self-renewal. These findings are consistent with the functional placement of Esrrb downstream of Nanog.},
	pages = {477--490},
	number = {4},
	journaltitle = {Cell Stem Cell},
	shortjournal = {Cell Stem Cell},
	author = {Festuccia, Nicola and Osorno, Rodrigo and Halbritter, Florian and Karwacki-Neisius, Violetta and Navarro, Pablo and Colby, Douglas and Wong, Frederick and Yates, Adam and Tomlinson, Simon R. and Chambers, Ian},
	urldate = {2018-04-23},
	date = {2012-10-05},
	keywords = {reading},
	file = {Festuccia et al. - 2012 - Esrrb Is a Direct Nanog Target Gene that Can Subst.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\8NWLKLUK\\Festuccia et al. - 2012 - Esrrb Is a Direct Nanog Target Gene that Can Subst.pdf:application/pdf;ScienceDirect Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\FD2JFVGB\\S1934590912004808.html:text/html}
}

@article{chambers_nanog_2007,
	title = {Nanog safeguards pluripotency and mediates germline development},
	volume = {450},
	rights = {2007 Nature Publishing Group},
	issn = {1476-4687},
	url = {https://www.nature.com/articles/nature06403},
	doi = {10.1038/nature06403},
	abstract = {Nanog is a divergent homeodomain protein found in mammalian pluripotent cells and developing germ cells1,2. Deletion of Nanog causes early embryonic lethality2, whereas constitutive expression enables autonomous self-renewal of embryonic stem cells1. Nanog is accordingly considered a core element of the pluripotent transcriptional network3,4,5,6,7. However, here we report that Nanog fluctuates in mouse embryonic stem cells. Transient downregulation of Nanog appears to predispose cells towards differentiation but does not mark commitment. By genetic deletion we show that, although they are prone to differentiate, embryonic stem cells can self-renew indefinitely in the permanent absence of Nanog. Expanded Nanog null cells colonize embryonic germ layers and exhibit multilineage differentiation both in fetal and adult chimaeras. Although they are also recruited to the germ line, primordial germ cells lacking Nanog fail to mature on reaching the genital ridge. This defect is rescued by repair of the mutant allele. Thus Nanog is dispensible for expression of somatic pluripotency but is specifically required for formation of germ cells. Nanog therefore acts primarily in construction of inner cell mass and germ cell states rather than in the housekeeping machinery of pluripotency. We surmise that Nanog stabilizes embryonic stem cells in culture by resisting or reversing alternative gene expression states.},
	pages = {1230--1234},
	number = {7173},
	journaltitle = {Nature},
	author = {Chambers, Ian and Silva, Jose and Colby, Douglas and Nichols, Jennifer and Nijmeijer, Bianca and Robertson, Morag and Vrana, Jan and Jones, Ken and Grotewold, Lars and Smith, Austin},
	urldate = {2018-04-22},
	date = {2007-12},
	langid = {english},
	keywords = {Nanog},
	file = {Chambers et al. - 2007 - Nanog safeguards pluripotency and mediates germlin.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\YNE4H6YA\\Chambers et al. - 2007 - Nanog safeguards pluripotency and mediates germlin.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\FYFF6K3V\\nature06403.html:text/html}
}

@article{ibarra-soria_defining_2018,
	title = {Defining murine organogenesis at single-cell resolution reveals a role for the leukotriene pathway in regulating blood progenitor formation},
	rights = {2018 © The Authors 2017, under exclusive licence to Macmillan Publishers Ltd., part of Springer Nature},
	issn = {1476-4679},
	url = {https://www.nature.com/articles/s41556-017-0013-z},
	doi = {10.1038/s41556-017-0013-z},
	abstract = {{\textless}p{\textgreater}Ibarra-Soria et al. study cellular diversity, transcriptional signatures, lineage specification and somitogenesis on a single-cell level in E8.25 mouse embryos, and reveal the regulation of blood progenitor formation by the leukotriene pathway.{\textless}/p{\textgreater}},
	pages = {1},
	journaltitle = {Nature Cell Biology},
	author = {Ibarra-Soria, Ximena and Jawaid, Wajid and Pijuan-Sala, Blanca and Ladopoulos, Vasileios and Scialdone, Antonio and Jörg, David J. and Tyser, Richard C. V. and Calero-Nieto, Fernando J. and Mulas, Carla and Nichols, Jennifer and Vallier, Ludovic and Srinivas, Shankar and Simons, Benjamin D. and Göttgens, Berthold and Marioni, John C.},
	urldate = {2018-01-22},
	date = {2018-01-08},
	file = {Ibarra-Soria et al. - 2018 - Defining murine organogenesis at single-cell resol.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\I3G7T6VR\\Ibarra-Soria et al. - 2018 - Defining murine organogenesis at single-cell resol.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\A9JABH5W\\s41556-017-0013-z.html:text/html}
}

@article{ferraris_combinatorial_2011,
	title = {Combinatorial binding of transcription factors in the pluripotency control regions of the genome},
	volume = {21},
	issn = {1088-9051, 1549-5469},
	url = {http://genome.cshlp.org/content/21/7/1055},
	doi = {10.1101/gr.115824.110},
	abstract = {The pluripotency control regions ({PluCRs}) are defined as genomic regions that are bound by {POU}5F1, {SOX}2, and {NANOG} in vivo. We utilized a high-throughput binding assay to record more than 270,000 different {DNA}/protein binding measurements along incrementally tiled windows of {DNA} within these {PluCRs}. This high-resolution binding map is then used to systematically define the context of {POU} factor binding, and reveals patterns of cooperativity and competition in the pluripotency network. The most prominent pattern is a pervasive binding competition between {POU}5F1 and the forkhead transcription factors. Like many transcription factors, {POU}5F1 is co-expressed with a paralog, {POU}2F1, that shares an apparently identical binding specificity. By analyzing thousands of binding measurements, we discover context effects that discriminate {POU}2F1 from {POU}5F1 binding. Proximal {NANOG} binding promotes {POU}5F1 binding, whereas nearby {SOX}2 binding favors {POU}2F1. We demonstrate by cross-species comparison and by chromatin immunoprecipitation ({ChIP}) that the contextual sequence determinants learned in vitro are sufficient to predict {POU}2F1 binding in vivo.},
	pages = {1055--1064},
	number = {7},
	journaltitle = {Genome Research},
	shortjournal = {Genome Res.},
	author = {Ferraris, Luciana and Stewart, Allan P. and Kang, Jinsuk and {DeSimone}, Alec M. and Gemberling, Matthew and Tantin, Dean and Fairbrother, William G.},
	urldate = {2018-06-07},
	date = {2011-07-01},
	langid = {english},
	pmid = {21527551},
	keywords = {by Natashia},
	file = {Full Text PDF:C\:\\Users\\Ivan\\Zotero\\storage\\D7UL69US\\Ferraris et al. - 2011 - Combinatorial binding of transcription factors in .pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\Y55ALV9B\\1055.html:text/html}
}

@article{chen_integration_2008,
	title = {Integration of External Signaling Pathways with the Core Transcriptional Network in Embryonic Stem Cells},
	volume = {133},
	issn = {0092-8674},
	url = {http://www.sciencedirect.com/science/article/pii/S009286740800617X},
	doi = {10.1016/j.cell.2008.04.043},
	abstract = {Summary
Transcription factors ({TFs}) and their specific interactions with targets are crucial for specifying gene-expression programs. To gain insights into the transcriptional regulatory networks in embryonic stem ({ES}) cells, we use chromatin immunoprecipitation coupled with ultra-high-throughput {DNA} sequencing ({ChIP}-seq) to map the locations of 13 sequence-specific {TFs} (Nanog, Oct4, {STAT}3, Smad1, Sox2, Zfx, c-Myc, n-Myc, Klf4, Esrrb, Tcfcp2l1, E2f1, and {CTCF}) and 2 transcription regulators (p300 and Suz12). These factors are known to play different roles in {ES}-cell biology as components of the {LIF} and {BMP} signaling pathways, self-renewal regulators, and key reprogramming factors. Our study provides insights into the integration of the signaling pathways into the {ES}-cell-specific transcription circuitries. Intriguingly, we find specific genomic regions extensively targeted by different {TFs}. Collectively, the comprehensive mapping of {TF}-binding sites identifies important features of the transcriptional regulatory networks that define {ES}-cell identity.},
	pages = {1106--1117},
	number = {6},
	journaltitle = {Cell},
	shortjournal = {Cell},
	author = {Chen, Xi and Xu, Han and Yuan, Ping and Fang, Fang and Huss, Mikael and Vega, Vinsensius B. and Wong, Eleanor and Orlov, Yuriy L. and Zhang, Weiwei and Jiang, Jianming and Loh, Yuin-Han and Yeo, Hock Chuan and Yeo, Zhen Xuan and Narang, Vipin and Govindarajan, Kunde Ramamoorthy and Leong, Bernard and Shahab, Atif and Ruan, Yijun and Bourque, Guillaume and Sung, Wing-Kin and Clarke, Neil D. and Wei, Chia-Lin and Ng, Huck-Hui},
	urldate = {2018-06-07},
	date = {2008-06-13},
	keywords = {by Natashia, reading},
	file = {ScienceDirect Full Text PDF:C\:\\Users\\Ivan\\Zotero\\storage\\HH6Q3KLY\\Chen et al. - 2008 - Integration of External Signaling Pathways with th.pdf:application/pdf;ScienceDirect Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\HAMS97IH\\S009286740800617X.html:text/html}
}

@article{loh_oct4_2006,
	title = {The Oct4 and Nanog transcription network regulates pluripotency in mouse embryonic stem cells},
	volume = {38},
	rights = {2006 Nature Publishing Group},
	issn = {1546-1718},
	url = {https://www.nature.com/articles/ng1760},
	doi = {10.1038/ng1760},
	abstract = {Oct4 and Nanog are transcription factors required to maintain the pluripotency and self-renewal of embryonic stem ({ES}) cells. Using the chromatin immunoprecipitation paired-end ditags method, we mapped the binding sites of these factors in the mouse {ES} cell genome. We identified 1,083 and 3,006 high-confidence binding sites for Oct4 and Nanog, respectively. Comparative location analyses indicated that Oct4 and Nanog overlap substantially in their targets, and they are bound to genes in different configurations. Using de novo motif discovery algorithms, we defined the cis-acting elements mediating their respective binding to genomic sites. By integrating {RNA} interference–mediated depletion of Oct4 and Nanog with microarray expression profiling, we demonstrated that these factors can activate or suppress transcription. We further showed that common core downstream targets are important to keep {ES} cells from differentiating. The emerging picture is one in which Oct4 and Nanog control a cascade of pathways that are intricately connected to govern pluripotency, self-renewal, genome surveillance and cell fate determination.},
	pages = {431--440},
	number = {4},
	journaltitle = {Nature Genetics},
	author = {Loh, Yuin-Han and Wu, Qiang and Chew, Joon-Lin and Vega, Vinsensius B. and Zhang, Weiwei and Chen, Xi and Bourque, Guillaume and George, Joshy and Leong, Bernard and Liu, Jun and Wong, Kee-Yew and Sung, Ken W. and Lee, Charlie W. H. and Zhao, Xiao-Dong and Chiu, Kuo-Ping and Lipovich, Leonard and Kuznetsov, Vladimir A. and Robson, Paul and Stanton, Lawrence W. and Wei, Chia-Lin and Ruan, Yijun and Lim, Bing and Ng, Huck-Hui},
	urldate = {2018-06-07},
	date = {2006-04},
	langid = {english},
	keywords = {by Natashia},
	file = {Full Text PDF:C\:\\Users\\Ivan\\Zotero\\storage\\URTPS93V\\Loh et al. - 2006 - The Oct4 and Nanog transcription network regulates.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\L7I2G36N\\ng1760.html:text/html}
}

@article{gagliardi_direct_2013,
	title = {A direct physical interaction between Nanog and Sox2 regulates embryonic stem cell self‐renewal},
	volume = {32},
	rights = {Copyright © 2013 European Molecular Biology Organization. This is an open access article under the terms of the Creative Commons Attribution License, which permits use, distribution and reproduction in any medium, provided the original work is properly cited.},
	issn = {0261-4189, 1460-2075},
	url = {http://emboj.embopress.org/content/32/16/2231},
	doi = {10.1038/emboj.2013.161},
	abstract = {Embryonic stem ({ES}) cell self‐renewal efficiency is determined by the Nanog protein level. However, the protein partners of Nanog that function to direct self‐renewal are unclear. Here, we identify a Nanog interactome of over 130 proteins including transcription factors, chromatin modifying complexes, phosphorylation and ubiquitination enzymes, basal transcriptional machinery members, and {RNA} processing factors. Sox2 was identified as a robust interacting partner of Nanog. The purified Nanog–Sox2 complex identified a {DNA} recognition sequence present in multiple overlapping Nanog/Sox2 {ChIP}‐Seq data sets. The Nanog tryptophan repeat region is necessary and sufficient for interaction with Sox2, with tryptophan residues required. In Sox2, tyrosine to alanine mutations within a triple‐repeat motif (S X T/S Y) abrogates the Nanog–Sox2 interaction, alters expression of genes associated with the Nanog‐Sox2 cognate sequence, and reduces the ability of Sox2 to rescue {ES} cell differentiation induced by endogenous Sox2 deletion. Substitution of the tyrosines with phenylalanine rescues both the Sox2–Nanog interaction and efficient self‐renewal. These results suggest that aromatic stacking of Nanog tryptophans and Sox2 tyrosines mediates an interaction central to {ES} cell self‐renewal.},
	pages = {2231--2247},
	number = {16},
	journaltitle = {The {EMBO} Journal},
	author = {Gagliardi, Alessia and Mullin, Nicholas P. and Tan, Zi Ying and Colby, Douglas and Kousa, Anastasia I. and Halbritter, Florian and Weiss, Jason T. and Felker, Anastasia and Bezstarosti, Karel and Favaro, Rebecca and Demmers, Jeroen and Nicolis, Silvia K. and Tomlinson, Simon R. and Poot, Raymond A. and Chambers, Ian},
	urldate = {2018-06-07},
	date = {2013-08-14},
	langid = {english},
	pmid = {23892456},
	keywords = {by Natashia},
	file = {Full Text PDF:C\:\\Users\\Ivan\\Zotero\\storage\\QAJIMJEK\\Gagliardi et al. - 2013 - A direct physical interaction between Nanog and So.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\AUFMXB7Q\\2231.html:text/html}
}

@article{bedzhov_developmental_2014,
	title = {Developmental plasticity, cell fate specification and morphogenesis in the early mouse embryo},
	volume = {369},
	rights = {. © 2014 The Authors. Published by the Royal Society under the terms of the Creative Commons Attribution License http://creativecommons.org/licenses/by/3.0/, which permits unrestricted use, provided the original author and source are credited.},
	issn = {0962-8436, 1471-2970},
	url = {http://rstb.royalsocietypublishing.org/content/369/1657/20130538},
	doi = {10.1098/rstb.2013.0538},
	abstract = {A critical point in mammalian development is when the early embryo implants into its mother's uterus. This event has historically been difficult to study due to the fact that it occurs within the maternal tissue and therefore is hidden from view. In this review, we discuss how the mouse embryo is prepared for implantation and the molecular mechanisms involved in directing and coordinating this crucial event. Prior to implantation, the cells of the embryo are specified as precursors of future embryonic and extra-embryonic lineages. These preimplantation cell fate decisions rely on a combination of factors including cell polarity, position and cell–cell signalling and are influenced by the heterogeneity between early embryo cells. At the point of implantation, signalling events between the embryo and mother, and between the embryonic and extraembryonic compartments of the embryo itself, orchestrate a total reorganization of the embryo, coupled with a burst of cell proliferation. New developments in embryo culture and imaging techniques have recently revealed the growth and morphogenesis of the embryo at the time of implantation, leading to a new model for the blastocyst to egg cylinder transition. In this model, pluripotent cells that will give rise to the fetus self-organize into a polarized three-dimensional rosette-like structure that initiates egg cylinder formation.},
	pages = {20130538},
	number = {1657},
	journaltitle = {Phil. Trans. R. Soc. B},
	shortjournal = {Phil. Trans. R. Soc. B},
	author = {Bedzhov, Ivan and Graham, Sarah J. L. and Leung, Chuen Yan and Zernicka-Goetz, Magdalena},
	urldate = {2018-07-05},
	date = {2014-12-05},
	langid = {english},
	pmid = {25349447},
	file = {Bedzhov et al. - 2014 - Developmental plasticity, cell fate specification .pdf:C\:\\Users\\Ivan\\Zotero\\storage\\PWEWZQ86\\Bedzhov et al. - 2014 - Developmental plasticity, cell fate specification .pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\6DPBGWNC\\20130538.html:text/html}
}

@article{takaoka_cell_2012,
	title = {Cell fate decisions and axis determination in the early mouse embryo},
	volume = {139},
	rights = {© 2012.},
	issn = {0950-1991, 1477-9129},
	url = {http://dev.biologists.org/content/139/1/3},
	doi = {10.1242/dev.060095},
	abstract = {{\textless}p{\textgreater}The mouse embryo generates multiple cell lineages, as well as its future body axes in the early phase of its development. The early cell fate decisions lead to the generation of three lineages in the pre-implantation embryo: the epiblast, the primitive endoderm and the trophectoderm. Shortly after implantation, the anterior-posterior axis is firmly established. Recent studies have provided a better understanding of how the earliest cell fate decisions are regulated in the pre-implantation embryo, and how and when the body axes are established in the pregastrulation embryo. In this review, we address the timing of the first cell fate decisions and of the establishment of embryonic polarity, and we ask how far back one can trace their origins.{\textless}/p{\textgreater}},
	pages = {3--14},
	number = {1},
	journaltitle = {Development},
	shortjournal = {Development},
	author = {Takaoka, Katsuyoshi and Hamada, Hiroshi},
	urldate = {2018-07-19},
	date = {2012-01-01},
	langid = {english},
	pmid = {22147950},
	file = {Full Text PDF:C\:\\Users\\Ivan\\Zotero\\storage\\TV4DCJ6F\\Takaoka and Hamada - 2012 - Cell fate decisions and axis determination in the .pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\822VZ2ZU\\3.html:text/html}
}

@article{davidson_pluripotent_2015,
	title = {The pluripotent state in mouse and human},
	volume = {142},
	issn = {0950-1991, 1477-9129},
	url = {http://dev.biologists.org/cgi/doi/10.1242/dev.116061},
	doi = {10.1242/dev.116061},
	abstract = {In the mouse, naïve pluripotent stem cells ({PSCs}) are thought to represent the cell culture equivalent of the late epiblast in the preimplantation embryo, with which they share a unique defining set of features. Recent studies have focused on the identification and propagation of a similar cell state in human. Although the capture of an exact human equivalent of the mouse naïve {PSC} remains an elusive goal, comparative studies spurred on by this quest are lighting the path to a deeper understanding of pluripotent state regulation in early mammalian development.},
	pages = {3090--3099},
	number = {18},
	journaltitle = {Development},
	author = {Davidson, K. C. and Mason, E. A. and Pera, M. F.},
	urldate = {2018-08-03},
	date = {2015-09-15},
	langid = {english},
	file = {Davidson et al. - 2015 - The pluripotent state in mouse and human.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\I5WM4A8M\\Davidson et al. - 2015 - The pluripotent state in mouse and human.pdf:application/pdf}
}

@article{moore_primitive_2014,
	title = {The Primitive Endoderm Segregates from the Epiblast in 1 Integrin-Deficient Early Mouse Embryos},
	volume = {34},
	issn = {0270-7306},
	url = {http://mcb.asm.org/cgi/doi/10.1128/MCB.00937-13},
	doi = {10.1128/MCB.00937-13},
	pages = {560--572},
	number = {3},
	journaltitle = {Molecular and Cellular Biology},
	author = {Moore, R. and Tao, W. and Smith, E. R. and Xu, X.-X.},
	urldate = {2018-08-03},
	date = {2014-02-01},
	langid = {english},
	file = {Moore et al. - 2014 - The Primitive Endoderm Segregates from the Epiblas.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\7VBA8PKF\\Moore et al. - 2014 - The Primitive Endoderm Segregates from the Epiblas.pdf:application/pdf}
}

@article{aksoy_oct4_2013,
	title = {Oct4 switches partnering from Sox2 to Sox17 to reinterpret the enhancer code and specify endoderm},
	volume = {32},
	rights = {Copyright © 2013 European Molecular Biology Organization},
	issn = {0261-4189, 1460-2075},
	url = {http://emboj.embopress.org/content/32/7/938},
	doi = {10.1038/emboj.2013.31},
	abstract = {How regulatory information is encoded in the genome is poorly understood and poses a challenge when studying biological processes. We demonstrate here that genomic redistribution of Oct4 by alternative partnering with Sox2 and Sox17 is a fundamental regulatory event of endodermal specification. We show that Sox17 partners with Oct4 and binds to a unique ‘compressed’ Sox/Oct motif that earmarks endodermal genes. This is in contrast to the pluripotent state where Oct4 selectively partners with Sox2 at ‘canonical’ binding sites. The distinct selection of binding sites by alternative Sox/Oct partnering is underscored by our demonstration that rationally point‐mutated Sox17 partners with Oct4 on pluripotency genes earmarked by the canonical Sox/Oct motif. In an endodermal differentiation assay, we demonstrate that the compressed motif is required for proper expression of endodermal genes. Evidently, Oct4 drives alternative developmental programs by switching Sox partners that affects enhancer selection, leading to either an endodermal or pluripotent cell fate. This work provides insights in understanding cell fate transcriptional regulation by highlighting the direct link between the {DNA} sequence of an enhancer and a developmental outcome.
There is a Have you seen? (April 2013) associated with this Article.},
	pages = {938--953},
	number = {7},
	journaltitle = {The {EMBO} Journal},
	author = {Aksoy, Irene and Jauch, Ralf and Chen, Jiaxuan and Dyla, Mateusz and Divakar, Ushashree and Bogu, Gireesh K. and Teo, Roy and Ng, Calista Keow Leng and Herath, Wishva and Lili, Sun and Hutchins, Andrew P. and Robson, Paul and Kolatkar, Prasanna R. and Stanton, Lawrence W.},
	urldate = {2018-08-20},
	date = {2013-04-03},
	langid = {english},
	pmid = {23474895},
	file = {Aksoy et al. - 2013 - Oct4 switches partnering from Sox2 to Sox17 to rei.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\KKX95V86\\Aksoy et al. - 2013 - Oct4 switches partnering from Sox2 to Sox17 to rei.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\L6MHELN5\\938.html:text/html}
}

@article{tam_gene_2007,
	title = {Gene function in mouse embryogenesis: get set for gastrulation},
	volume = {8},
	rights = {2007 Nature Publishing Group},
	issn = {1471-0064},
	url = {http://www.nature.com/articles/nrg2084},
	doi = {10.1038/nrg2084},
	shorttitle = {Gene function in mouse embryogenesis},
	abstract = {During early mouse embryogenesis, temporal and spatial regulation of gene expression and cell signalling influences lineage specification, embryonic polarity, the patterning of tissue progenitors and the morphogenetic movement of cells and tissues. Uniquely in mammals, the extraembryonic tissues are the source of signals for lineage specification and tissue patterning. Here we discuss recent discoveries about the lead up to gastrulation, including early manifestations of asymmetry, coordination of cell and tissue movements and the interactions of transcription factors and signalling activity for lineage allocation and germ-layer specification.},
	pages = {368--381},
	number = {5},
	journaltitle = {Nature Reviews Genetics},
	author = {Tam, Patrick P. L. and Loebel, David A. F.},
	urldate = {2018-08-20},
	date = {2007-05},
	langid = {english},
	file = {Full Text PDF:C\:\\Users\\Ivan\\Zotero\\storage\\URPHNV8R\\Tam and Loebel - 2007 - Gene function in mouse embryogenesis get set for .pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\9Y6TVDJ8\\nrg2084.html:text/html}
}

@article{oliveri_global_2008,
	title = {Global regulatory logic for specification of an embryonic cell lineage},
	volume = {105},
	rights = {© 2008 by The National Academy of Sciences of the {USA}},
	issn = {0027-8424, 1091-6490},
	url = {http://www.pnas.org/content/105/16/5955},
	doi = {10.1073/pnas.0711220105},
	abstract = {Explanation of a process of development must ultimately be couched in the terms of the genomic regulatory code. Specification of an embryonic cell lineage is driven by a network of interactions among genes encoding transcription factors. Here, we present the gene regulatory network ({GRN}) that directs the specification of the skeletogenic micromere lineage of the sea urchin embryo. The {GRN} now includes all regulatory genes expressed in this lineage up to late blastula stage, as identified in a genomewide survey. The architecture of the {GRN} was established by a large-scale perturbation analysis in which the expression of each gene in the {GRN} was cut off by use of morpholinos, and the effects on all other genes were measured quantitatively. Several cis-regulatory analyses provided additional evidence. The explanatory power of the {GRN} suffices to provide a causal explanation for all observable developmental functions of the micromere lineage during the specification period. These functions are: (i) initial acquisition of identity through transcriptional interpretation of localized maternal cues; (ii) activation of specific regulatory genes by use of a double negative gate; (iii) dynamic stabilization of the regulatory state by activation of a feedback subcircuit; (iv) exclusion of alternative regulatory states; (v) presentation of a signal required by the micromeres themselves and of two different signals required for development of adjacent endomesodermal lineages; and (vi) lineage-specific activation of batteries of skeletogenic genes. The {GRN} precisely predicts gene expression responses and provides a coherent explanation of the biology of specification.},
	pages = {5955--5962},
	number = {16},
	journaltitle = {Proceedings of the National Academy of Sciences},
	shortjournal = {{PNAS}},
	author = {Oliveri, Paola and Tu, Qiang and Davidson, Eric H.},
	urldate = {2018-08-28},
	date = {2008-04-22},
	langid = {english},
	pmid = {18413610},
	file = {Oliveri et al. - 2008 - Global regulatory logic for specification of an em.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\7VDT2S7L\\Oliveri et al. - 2008 - Global regulatory logic for specification of an em.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\5R6S2YZ9\\5955.html:text/html}
}

@article{marson_connecting_2008,
	title = {Connecting {microRNA} Genes to the Core Transcriptional Regulatory Circuitry of Embryonic Stem Cells},
	volume = {134},
	issn = {0092-8674, 1097-4172},
	url = {https://www.cell.com/cell/abstract/S0092-8674(08)00938-0},
	doi = {10.1016/j.cell.2008.07.020},
	pages = {521--533},
	number = {3},
	journaltitle = {Cell},
	shortjournal = {Cell},
	author = {Marson, Alexander and Levine, Stuart S. and Cole, Megan F. and Frampton, Garrett M. and Brambrink, Tobias and Johnstone, Sarah and Guenther, Matthew G. and Johnston, Wendy K. and Wernig, Marius and Newman, Jamie and Calabrese, J. Mauro and Dennis, Lucas M. and Volkert, Thomas L. and Gupta, Sumeet and Love, Jennifer and Hannett, Nancy and Sharp, Phillip A. and Bartel, David P. and Jaenisch, Rudolf and Young, Richard A.},
	urldate = {2018-08-28},
	date = {2008-08-08},
	pmid = {18692474},
	file = {Marson et al. - 2008 - Connecting microRNA Genes to the Core Transcriptio.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\ZVTZT96C\\Marson et al. - 2008 - Connecting microRNA Genes to the Core Transcriptio.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\8WYTBKI6\\S0092-8674(08)00938-0.html:text/html}
}

@article{silva_nanog_2009,
	title = {Nanog Is the Gateway to the Pluripotent Ground State},
	volume = {138},
	issn = {0092-8674, 1097-4172},
	url = {https://www.cell.com/cell/abstract/S0092-8674(09)00969-6},
	doi = {10.1016/j.cell.2009.07.039},
	pages = {722--737},
	number = {4},
	journaltitle = {Cell},
	shortjournal = {Cell},
	author = {Silva, Jose and Nichols, Jennifer and Theunissen, Thorold W. and Guo, Ge and Oosten, Anouk L. van and Barrandon, Ornella and Wray, Jason and Yamanaka, Shinya and Chambers, Ian and Smith, Austin},
	urldate = {2018-08-28},
	date = {2009-08-21},
	pmid = {19703398},
	file = {Silva et al. - 2009 - Nanog Is the Gateway to the Pluripotent Ground Sta.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\UUXZQTXZ\\Silva et al. - 2009 - Nanog Is the Gateway to the Pluripotent Ground Sta.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\ZHTII66V\\S0092-8674(09)00969-6.html:text/html}
}

@article{wang_protein_2006,
	title = {A protein interaction network for pluripotency of embryonic stem cells},
	volume = {444},
	rights = {2006 Nature Publishing Group},
	issn = {1476-4687},
	url = {https://www.nature.com/articles/nature05284},
	doi = {10.1038/nature05284},
	abstract = {Embryonic stem ({ES}) cells are pluripotent1,2 and of therapeutic potential in regenerative medicine3,4. Understanding pluripotency at the molecular level should illuminate fundamental properties of stem cells and the process of cellular reprogramming. Through cell fusion the embryonic cell phenotype can be imposed on somatic cells, a process promoted by the homeodomain protein Nanog5, which is central to the maintenance of {ES} cell pluripotency6,7. Nanog is thought to function in concert with other factors such as Oct4 (ref. 8) and Sox2 (ref. 9) to establish {ES} cell identity. Here we explore the protein network in which Nanog operates in mouse {ES} cells. Using affinity purification of Nanog under native conditions followed by mass spectrometry, we have identified physically associated proteins. In an iterative fashion we also identified partners of several Nanog-associated proteins (including Oct4), validated the functional relevance of selected newly identified components and constructed a protein interaction network. The network is highly enriched for nuclear factors that are individually critical for maintenance of the {ES} cell state and co-regulated on differentiation. The network is linked to multiple co-repressor pathways and is composed of numerous proteins whose encoding genes are putative direct transcriptional targets of its members. This tight protein network seems to function as a cellular module dedicated to pluripotency.},
	pages = {364--368},
	number = {7117},
	journaltitle = {Nature},
	author = {Wang, Jianlong and Rao, Sridhar and Chu, Jianlin and Shen, Xiaohua and Levasseur, Dana N. and Theunissen, Thorold W. and Orkin, Stuart H.},
	urldate = {2018-08-28},
	date = {2006-11},
	langid = {english},
	file = {Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\AX8SJQ7Y\\nature05284.html:text/html;Wang et al. - 2006 - A protein interaction network for pluripotency of .pdf:C\:\\Users\\Ivan\\Zotero\\storage\\LRE5WLJP\\Wang et al. - 2006 - A protein interaction network for pluripotency of .pdf:application/pdf}
}

@article{hall_oct4_2009,
	title = {Oct4 and {LIF}/Stat3 Additively Induce Krüppel Factors to Sustain Embryonic Stem Cell Self-Renewal},
	volume = {5},
	issn = {1934-5909},
	url = {http://www.sciencedirect.com/science/article/pii/S1934590909005761},
	doi = {10.1016/j.stem.2009.11.003},
	abstract = {Summary
Embryonic stem cell ({ESC}) pluripotency is dependent on an intrinsic gene regulatory network centered on Oct4. Propagation of the pluripotent state is stimulated by the cytokine leukemia inhibitory factor ({LIF}) acting through the transcriptional regulator Stat3. Here, we show that this extrinsic stimulus converges with the intrinsic circuitry in Krüppel-factor activation. Oct4 primarily induces Klf2 while {LIF}/Stat3 selectively enhances Klf4 expression. Overexpression of either factor reduces {LIF} dependence, but with quantitative and qualitative differences. Unlike Klf4, Klf2 increases {ESC} clonogenicity, maintains undifferentiated {ESCs} in the genetic absence of Stat3, and confers resistance to {BMP}-induced differentiation. {ESCs} expanded with Klf2 remain capable of contributing to adult chimeras. Postimplantation-embryo-derived {EpiSCs} lack both Klf2 and Klf4 and expression of either can reinstate naive pluripotency. These findings indicate that Oct4 and Stat3 intersect in directing expression of Klf transcriptional regulators with overlapping properties that additively reinforce ground-state {ESC} pluripotency, identity, and self-renewal.},
	pages = {597--609},
	number = {6},
	journaltitle = {Cell Stem Cell},
	shortjournal = {Cell Stem Cell},
	author = {Hall, John and Guo, Ge and Wray, Jason and Eyres, Isobel and Nichols, Jennifer and Grotewold, Lars and Morfopoulou, Sofia and Humphreys, Peter and Mansfield, William and Walker, Rachael and Tomlinson, Simon and Smith, Austin},
	urldate = {2018-08-28},
	date = {2009-12-04},
	file = {Hall et al. - 2009 - Oct4 and LIFStat3 Additively Induce Krüppel Facto.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\68ZZQE56\\Hall et al. - 2009 - Oct4 and LIFStat3 Additively Induce Krüppel Facto.pdf:application/pdf;ScienceDirect Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\853VM6X6\\S1934590909005761.html:text/html}
}

@online{love_tximport_2018,
	title = {tximport - Import and summarize transcript-level estimates for transcript- and gene-level analysis},
	url = {http://bioconductor.org/packages/tximport/},
	abstract = {Imports transcript-level abundance, estimated counts and transcript lengths, and summarizes into matrices for use with downstream gene-level analysis packages. Average transcript length, weighted by sample-specific transcript abundance estimates, is provided as a matrix which can be used as an offset for different expression of gene-level counts.},
	titleaddon = {Bioconductor},
	author = {Love, Michael},
	urldate = {2018-06-26},
	date = {2018-06-25},
	langid = {american},
	keywords = {Bioconductor, R},
	file = {Importing transcript abundance datasets with tximport:C\:\\Users\\Ivan\\Zotero\\storage\\2PNHLG2I\\tximport.html:text/html;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\AZ78AHWJ\\tximport.html:text/html;tximport.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\TKKMB5DG\\tximport.pdf:application/pdf}
}

@online{mccarthy_scater_2018,
	title = {scater - Single-Cell Analysis Toolkit for Gene Expression Data in R},
	url = {http://bioconductor.org/packages/scater/},
	abstract = {A collection of tools for doing various analyses of single-cell {RNA}-seq gene expression data, with a focus on quality control.},
	titleaddon = {Bioconductor},
	author = {{McCarthy}, Davis},
	urldate = {2018-06-26},
	date = {2018-06-25},
	langid = {american},
	keywords = {Bioconductor, R},
	file = {scater.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\S8CW24V4\\scater.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\RIUH7JUR\\scater.html:text/html}
}

@online{rainer_ensembldb_2018,
	title = {ensembldb - Utilities to create and use Ensembl-based annotation databases},
	url = {http://bioconductor.org/packages/ensembldb/},
	abstract = {The package provides functions to create and use transcript centric annotation databases/packages. The annotation for the databases are directly fetched from Ensembl using their Perl {API}. The functionality and data is similar to that of the {TxDb} packages from the {GenomicFeatures} package, but, in addition to retrieve all gene/transcript models and annotations from the database, ensembldb provides a filter framework allowing to retrieve annotations for specific entries like genes encoded on a chromosome region or transcript models of {lincRNA} genes. {EnsDb} databases built with ensembldb contain also protein annotations and mappings between proteins and their encoding transcripts. Finally, ensembldb provides functions to map between genomic, transcript and protein coordinates.},
	titleaddon = {Bioconductor},
	author = {Rainer, Johannes},
	urldate = {2018-07-03},
	date = {2018-07-02},
	langid = {american},
	keywords = {Bioconductor, R},
	file = {ensembldb.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\QG3KRRBK\\ensembldb.pdf:application/pdf;Generating an using Ensembl based annotation packages:C\:\\Users\\Ivan\\Zotero\\storage\\FPKHRHB2\\ensembldb.html:text/html;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\DY7NAGPU\\ensembldb.html:text/html}
}

@online{rainer_ensdb.mmusculus.v79_2018,
	title = {{EnsDb}.Mmusculus.v79 - Ensembl based annotation package},
	url = {http://bioconductor.org/packages/EnsDb.Mmusculus.v79/},
	abstract = {Exposes an annotation databases generated from Ensembl.},
	titleaddon = {Bioconductor},
	author = {Rainer, Johannes},
	urldate = {2018-07-03},
	date = {2018-07-03},
	langid = {american},
	keywords = {Bioconductor, R},
	file = {Rainer - This data package was made from resources at Ensem.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\F7UKD86I\\Rainer - This data package was made from resources at Ensem.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\FT7NDCYF\\EnsDb.Mmusculus.v79.html:text/html}
}

@online{morgan_summarizedexperiment_2018,
	title = {{SummarizedExperiment} - {SummarizedExperiment} container},
	url = {http://bioconductor.org/packages/SummarizedExperiment/},
	abstract = {The {SummarizedExperiment} container contains one or more assays, each represented by a matrix-like object of numeric or other mode. The rows typically represent genomic ranges of interest and the columns represent samples.},
	titleaddon = {Bioconductor},
	author = {Morgan, Martin and Obenchain, Valerie and Hester, Jim and Pagès, Hervé},
	urldate = {2018-07-04},
	date = {2018-07-02},
	langid = {american},
	keywords = {Bioconductor, R},
	file = {Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\IS3693TR\\SummarizedExperiment.html:text/html;SummarizedExperiment for Coordinating Experimental Assays, Samples, and Regions of Interest:C\:\\Users\\Ivan\\Zotero\\storage\\FMU2MARR\\SummarizedExperiment.html:text/html;SummarizedExperiment.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\D8VT7KZZ\\SummarizedExperiment.pdf:application/pdf}
}

@online{pages_annotationdbi_2018,
	title = {{AnnotationDbi} - Annotation Database Interface},
	url = {http://bioconductor.org/packages/AnnotationDbi/},
	abstract = {Provides user interface and database connection code for annotation data packages using {SQLite} data storage.},
	titleaddon = {Bioconductor},
	author = {Pagès, Hervé and Carlson, Marc and Falcon, Seth and Li, Nianhua},
	urldate = {2018-07-05},
	date = {2018-07-03},
	langid = {american},
	keywords = {Bioconductor, R},
	file = {AnnotationDbi.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\A9QZ5VRY\\AnnotationDbi.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\NLB4CL5F\\AnnotationDbi.html:text/html}
}

@online{lun_simplesinglecell_2018,
	title = {{simpleSingleCell} - A step-by-step workflow for low-level analysis of single-cell {RNA}-seq data with Bioconductor},
	url = {http://bioconductor.org/packages/simpleSingleCell/},
	abstract = {This workflow implements a low-level analysis pipeline for {scRNA}-seq data using scran, scater and other Bioconductor packages. It describes how to perform quality control on the libraries, normalization of cell-specific biases, basic data exploration and cell cycle phase identification. Procedures to detect highly variable genes, significantly correlated genes and subpopulation-specific marker genes are also shown. These analyses are demonstrated on a range of publicly available {scRNA}-seq data sets.},
	titleaddon = {Bioconductor},
	author = {Lun, Aaron T. L. and {McCarthy}, Davis J. and Marioni, John C.},
	urldate = {2018-07-05},
	date = {2018-05-25},
	langid = {american},
	keywords = {Bioconductor, R},
	file = {Analyzing single-cell RNA sequencing data from droplet-based protocols:C\:\\Users\\Ivan\\Zotero\\storage\\CH46U3PS\\work-3-tenx.html:text/html;Analyzing single-cell RNA-seq data containing read counts:C\:\\Users\\Ivan\\Zotero\\storage\\CW9KU6W5\\work-1-reads.html:text/html;Analyzing single-cell RNA-seq data containing UMI counts:C\:\\Users\\Ivan\\Zotero\\storage\\ZXC7PXU4\\work-2-umis.html:text/html;Correcting batch effects in single-cell RNA-seq data:C\:\\Users\\Ivan\\Zotero\\storage\\6FIHGWH6\\work-5-mnn.html:text/html;Further strategies for analyzing single-cell RNA-seq data:C\:\\Users\\Ivan\\Zotero\\storage\\8RDA2GIM\\work-4-misc.html:text/html;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\N8G47PY9\\simpleSingleCell.html:text/html;Workflows for analyzing single-cell RNA-seq data with R/Bioconductor:C\:\\Users\\Ivan\\Zotero\\storage\\SN9FYRHD\\work-0-intro.html:text/html}
}

@online{durinck_biomart_2018,
	title = {{biomaRt} - Interface to {BioMart} databases (e.g. Ensembl, {COSMIC}, Wormbase and Gramene)},
	url = {http://bioconductor.org/packages/biomaRt/},
	abstract = {In recent years a wealth of biological data has become available in public data repositories. Easy access to these valuable data resources and firm integration with data analysis is needed for comprehensive bioinformatics data analysis. {biomaRt} provides an interface to a growing collection of databases implementing the {BioMart} software suite (http://www.biomart.org). The package enables retrieval of large amounts of data in a uniform way without the need to know the underlying database schemas or write complex {SQL} queries. Examples of {BioMart} databases are Ensembl, {COSMIC}, Uniprot, {HGNC}, Gramene, Wormbase and {dbSNP} mapped to Ensembl. These major databases give {biomaRt} users direct access to a diverse set of data and enable a wide range of powerful online queries from gene annotation to database mining.},
	titleaddon = {Bioconductor},
	author = {Durinck, Steffen and Huber, Wolfgang and Davis, Sean and Pepin, Francois and Buffalo, Vince S and Smith, Mike},
	urldate = {2018-07-07},
	date = {2018-07-05},
	langid = {american},
	keywords = {Bioconductor, R},
	file = {biomaRt.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\J49AGE2K\\biomaRt.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\2LKDXITR\\biomaRt.html:text/html}
}

@online{krueger_trim_2018,
	title = {Trim Galore!},
	url = {https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/},
	abstract = {A wrapper around Cutadapt and {FastQC} to consistently apply adapter and quality trimming to {FastQ} files, with extra functionality for {RRBS} data},
	titleaddon = {Babraham Bioinformatics},
	author = {Krueger, Felix},
	urldate = {2018-08-29},
	date = {2018-06-28},
	file = {Babraham Bioinformatics - Trim Galore!:C\:\\Users\\Ivan\\Zotero\\storage\\424Y6HHW\\trim_galore.html:text/html}
}

@article{gruning_bioconda:_2018,
	title = {Bioconda: sustainable and comprehensive software distribution for the life sciences},
	volume = {15},
	rights = {2018 The Author(s)},
	issn = {1548-7105},
	url = {https://www.nature.com/articles/s41592-018-0046-7},
	doi = {10.1038/s41592-018-0046-7},
	shorttitle = {Bioconda},
	abstract = {Correspondence},
	pages = {475--476},
	number = {7},
	journaltitle = {Nature Methods},
	author = {Grüning, Björn and Dale, Ryan and Sjödin, Andreas and Chapman, Brad A. and Rowe, Jillian and Tomkins-Tinch, Christopher H. and Valieris, Renan and Köster, Johannes},
	urldate = {2018-08-29},
	date = {2018-07},
	langid = {english},
	file = {Grüning et al. - 2018 - Bioconda sustainable and comprehensive software d.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\L22D2A4G\\Grüning et al. - 2018 - Bioconda sustainable and comprehensive software d.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\XFDAMEEP\\s41592-018-0046-7.html:text/html}
}

@article{ritchie_limma_2015,
	title = {limma powers differential expression analyses for {RNA}-sequencing and microarray studies},
	volume = {43},
	issn = {0305-1048},
	url = {https://academic.oup.com/nar/article/43/7/e47/2414268},
	doi = {10.1093/nar/gkv007},
	abstract = {Abstract.  limma is an R/Bioconductor software package that provides an integrated solution for analysing data from gene expression experiments. It contains ric},
	pages = {e47--e47},
	number = {7},
	journaltitle = {Nucleic Acids Research},
	shortjournal = {Nucleic Acids Res},
	author = {Ritchie, Matthew E. and Phipson, Belinda and Wu, Di and Hu, Yifang and Law, Charity W. and Shi, Wei and Smyth, Gordon K.},
	urldate = {2018-08-29},
	date = {2015-04-20},
	langid = {english},
	file = {Ritchie et al. - 2015 - limma powers differential expression analyses for .pdf:C\:\\Users\\Ivan\\Zotero\\storage\\TYMEVUXM\\Ritchie et al. - 2015 - limma powers differential expression analyses for .pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\3C89HCIF\\2414268.html:text/html}
}

@online{risso_singlecellexperiment_2018,
	title = {{SingleCellExperiment} - S4 Classes for Single Cell Data},
	url = {http://bioconductor.org/packages/SingleCellExperiment/},
	abstract = {Defines a S4 class for storing data from single-cell experiments. This includes specialized methods to store and retrieve spike-in information, dimensionality reduction coordinates and size factors for each cell, along with the usual metadata for genes and libraries.},
	titleaddon = {Bioconductor},
	author = {Risso, Davide},
	urldate = {2018-07-12},
	date = {2018-07-11},
	langid = {american},
	file = {SingleCellExperiment.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\PY2M5BBQ\\SingleCellExperiment.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\TN37WP22\\SingleCellExperiment.html:text/html}
}

@online{lun_scran_2018,
	title = {scran - Methods for Single-Cell {RNA}-Seq Data Analysis},
	url = {http://bioconductor.org/packages/scran/},
	abstract = {Implements functions for low-level analyses of single-cell {RNA}-seq data. Methods are provided for normalization of cell-specific biases, assignment of cell cycle phase, detection of highly variable and significantly correlated genes, correction of batch effects, identification of marker genes, and other common tasks in single-cell analysis workflows.},
	titleaddon = {Bioconductor},
	author = {Lun, Aaron},
	urldate = {2018-07-18},
	date = {2018-07-17},
	langid = {american},
	file = {scran.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\BBHTU2ZZ\\scran.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\CQDUNEE2\\scran.html:text/html;Using scran to perform basic analyses of single-cell RNA-seq data:C\:\\Users\\Ivan\\Zotero\\storage\\SEUBNXGM\\scran.html:text/html}
}

@online{langfelder_dynamictreecut:_2016,
	title = {{dynamicTreeCut}: Methods for Detection of Clusters in Hierarchical Clustering Dendrograms},
	rights = {{GPL}-2 {\textbar} {GPL}-3 [expanded from: {GPL} (≥ 2)]},
	url = {https://CRAN.R-project.org/package=dynamicTreeCut},
	shorttitle = {{dynamicTreeCut}},
	abstract = {Contains methods for detection of clusters in hierarchical clustering dendrograms.},
	titleaddon = {{CRAN}},
	author = {Langfelder, Peter and Zhang, Bin},
	urldate = {2018-08-20},
	date = {2016-03-11},
	file = {dynamicTreeCut.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\9CLVSSL4\\dynamicTreeCut.pdf:application/pdf}
}

@online{michael_lawrence_rtracklayer_2018,
	title = {rtracklayer - R interface to genome annotation files and the {UCSC} genome browser},
	url = {http://bioconductor.org/packages/rtracklayer/},
	abstract = {Extensible framework for interacting with multiple genome browsers (currently {UCSC} built-in) and manipulating annotation tracks in various formats (currently {GFF}, {BED}, {bedGraph}, {BED}15, {WIG}, {BigWig} and 2bit built-in). The user may export/import tracks to/from the supported browsers, as well as query and modify the browser state, such as the current viewport.},
	titleaddon = {Bioconductor},
	author = {Michael Lawrence},
	urldate = {2018-08-21},
	date = {2018-08-10},
	langid = {american},
	file = {rtracklayer_manual.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\YWYMWZRD\\rtracklayer_manual.pdf:application/pdf;rtracklayer_vignette.pdf:C\:\\Users\\Ivan\\Zotero\\storage\\LF6HKTX4\\rtracklayer_vignette.pdf:application/pdf;Snapshot:C\:\\Users\\Ivan\\Zotero\\storage\\3L9GTPDN\\rtracklayer.html:text/html}
}
