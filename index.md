# CRUK Bioinformatics Autumn School 2017: Functional Genomics

## 18th - 22nd September 2017: Craik-Marshall Room, Downing Site, University of Cambridge


# Overview

Functional genomics looks at the dynamic aspects of how the genome functions within cells,
particularly in the form of gene expression (transcription) and gene regulation. This workshop surveys
current methods for functional genomics using high-throughput technologies. 

High-throughput technologies such as next generation sequencing (NGS) can routinely produce massive amounts of data. However, such datasets pose new challenges in the way the data have to be analyzed, annotated and interpreted which are not trivial and are daunting to the wet-lab biologist. This course covers state-of-the-art and best-practice tools for bulk RNA-seq and ChIP-seq data analysis, and will also introduce approaches for analysing data arising from single-cell RNA-seq studies.

# Audience

Enthusiastic and motivated wet-lab biologists who want to gain more of an understanding of NGS data and eventually progress to analysing their own data

# Pre-requisites

The course will include a great deal of hands-on work in R and at the command line. In order for you to make the most of the course we strongly recommend that you take an introductory course, or have sufficient experience in the following areas:

- R
 - Unix
 - Introductory statistics

More specific requirements and references can be found [here](http://www.cruk.cam.ac.uk/bioinformatics-summer-school-prerequisites)


# Instructors

- [Rory Stark (CRUK CI)](http://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
- [Mark Dunning (CRUK CI)](http://markdunning.github.io/)
- [Shamith Samarajiwa (MRC CU)](http://www.mrc-cu.cam.ac.uk/research/Shamith-Samarajiwa-folder)
- [Dora Bihary (MRC CU)](http://www.mrc-cu.cam.ac.uk/research/Shamith-Samarajiwa-folder)
- [Ashley Sawle (CRUK CI)](http://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
- [Alistair Martin (CRUK CI)](http://www.cruk.cam.ac.uk/research-groups/caldas-group)
- [Stephane Ballereau (CRUK CI)](http://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
- [Davis McCarthy (EBI)](https://sites.google.com/site/davismcc/home)


# Aims
During this course you will learn about:-

- To provide an understanding of how aligned sequencing reads, genome sequences and genomic regions are represented in R.
- To encourage confidence in reading sequencing reads into R, performing quality assessment and executing standard pipelines for (bulk) RNA-Seq and ChIP-Seq analysis 
- Analysis of transcription factor (TF) and epigenomic (histone mark) ChIP-seq data 
- Recent advances in single-cell sequencing

# Objectives
After the course you should be able to:-

- Know what tools are available in Bioconductor for HTS analysis and understand the basic object-types that are utilised.
- Process and quality control short read sequencing data 
- Given a set of gene identifiers, find out whereabouts in the genome they are located, and vice-versa 
- Produce a list of differentially expressed genes from an RNA-Seq experiment.
- Import a set of ChIP-Seq peaks and investigate their biological context.
- Appreciate the differences between bulk and single-cell RNA-seq analyses, and why the same methodologies might not be applicable

# Materials

# Day 1 (September 18th)

- 09:30 [Course Introduction](Introduction/Session1-intro.html)
- 09:30 - 10:30; [Introduction to Functional Genomics](Introduction/Functional Genomics Overview.pdf)
- 10:30 - 12:30; [Introduction (Recap) of R and Bioconductor](Introduction/bioc-intro.nb.html)
  + [source](Introduction/bioc-intro.Rmd)
- 12:30 - 13:30; LUNCH
- 13:30 - 14:30 [Principles of Experimental Design](Introduction/Experimental Design.pdf)
- 14:30 - 17:00;
    + Data Processing for Next Generation Sequencing
    + Alignment to reference genomes and QC 
    
# Day 2 (September 19th)

- 09:30 - 10:30; [Introduction to RNA-seq](RNASeq/slides/rnaSeq_Sept2017.pdf)
- 10:30 - 12:30; 
  + [Counting](RNASeq/count.nb.html)
  + [Importing and QC of RNA-seq data](RNASeq/rna-seq-preprocessing.nb.html)
  + [source file](RNASeq/slides/rna-seq-preprocessing.nb.Rmd)
- 12:30 - 13:30; LUNCH
- 13:30 - 17:00; Differential Expression
  + Slides
  + [Practical](DifferentialExpression/rna-seq-de.nb.html)
  + [Source](DifferentialExpression/rna-seq.Rmd)

# Day 3 (September 20th)

- 09:30 - 11:00; [Annotation and Visualisation of Differential Expression](RNASeq/rna-seq-annotation-visualisation.nb.html)
  + [source file](RNASeq/rna-seq-annotation-visualisation.nb.Rmd)
- 11:00 - 12:30; [Gene set analysis and Gene Ontology testing](RNASeq/rna-seq-gene-set-testing.nb.html)
  + [source file](RNASeq/rna-seq-gene-set-testing.Rmd)
- 12:30 - 13:30; LUNCH
- 13:30 - 17:00; 
  + Introduction to ChIP-seq
  + Quality control methods for ChIP-seq 
  + Introduction to BedTools, GenomicRanges et al 

# Day 4 (September 21st)

- 09:30 - 17:00; 
  + Downloading public ChIP-seq data
  + Downstream analysis of public ChIP-seq datasets
  + Differential binding analysis
  + Identifying direct targets of Transcription Factors
  + Introduction to the analysis of long distance interactions (Hi-C)
  + Introduction to the analysis of specialised interaction data
  
# Day 5 (September 22nd)

- 09:30 - 12:30; [Introduction to single-cell sequencing](SingleCell/index.html)

## Data

- Mouse mammary data (counts): [https://figshare.com/s/1d788fd384d33e913a2a](https://figshare.com/s/1d788fd384d33e913a2a)
