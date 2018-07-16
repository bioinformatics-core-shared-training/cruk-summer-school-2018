# CRUK Bioinformatics Summer School 2018: Functional Genomics

## 23rd - 27th July 2018: Craik-Marshall Room, Downing Site, University of Cambridge

# Overview

Functional genomics looks at the dynamic aspects of how the genome functions within cells,
particularly in the form of gene expression (transcription) and gene regulation. This workshop surveys
current methods for functional genomics using high-throughput technologies. 

High-throughput technologies such as next generation sequencing (NGS) can routinely produce massive amounts of data. However, such datasets pose new challenges in the way the data have to be analyzed, annotated and interpreted which are not trivial and are daunting to the wet-lab biologist. This course covers state-of-the-art and best-practice tools for bulk RNA-seq and ChIP-seq data analysis, and will also introduce approaches for analysing data arising from single-cell RNA-seq studies.

# Audience

Enthusiastic and motivated wet-lab biologists who want to gain more of an understanding of NGS data and eventually progress to analysing their own data

# Pre-requisites

**The course will include a great deal of hands-on work in R and at the command line. In order for you to make the most of the course we strongly recommend that you take an introductory course, or have sufficient experience in the following areas:**

- R
- Unix
- Introductory statistics

**More specific requirements and references can be found [here](http://www.cruk.cam.ac.uk/bioinformatics-summer-school-prerequisites)**


# Instructors

- [Mark Fernandes (CRUK CI)](http://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
- [Shamith Samarajiwa (MRC CU)](http://www.mrc-cu.cam.ac.uk/research/Shamith-Samarajiwa-folder)
- [Dora Bihary (MRC CU)](http://www.mrc-cu.cam.ac.uk/research/Shamith-Samarajiwa-folder)
- [Ashley Sawle (CRUK CI)](http://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
- [Abigail Edwards (CRUK CI)](http://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
- [Alistair Martin (CRUK CI)](http://www.cruk.cam.ac.uk/research-groups/caldas-group)
- [Stephane Ballereau (CRUK CI)](http://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
- [Michael Morgan (CRUK CI)](http://www.cruk.cam.ac.uk/).  


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

# Day 0 (July 22nd )
18:00 - ..
Informal get-together at The Grain and Hop Store (close to accommodation in Downing College)
Join us for a drink and dinner (self-paying), and to meet your colleagues for the next few days
http://www.grainandhopstore-cambridge.co.uk/

# Materials
Note Training Room in Craik-Marshall building (1st Floor) will be open from 9am.
# Day 1 (July 23rd)

- 09:30 [Course Introduction](Introduction/Session1-intro.html)
- 09:30 - 10:30; [Introduction to Functional Genomics](Introduction/Functional_Genomics_Overview.pdf)
- 10:30 - 12:30; [Introduction (Recap) of R and Bioconductor](Introduction/bioc-intro.nb.html)
  + [source](Introduction/bioc-intro.Rmd)
- 12:30 - 13:30; LUNCH
- 13:30 - 14:30 [Principles of Experimental Design](Introduction/Experimental_Design.pdf)
- 14:30 - 17:00;
    + [File formats (lecture)](Introduction/SS_DB/Materials/Lectures/Lecture1_fileFormats_DB.pdf)
    + [Quality Control and artefacts (lecture)](Introduction/SS_DB/Materials/Lectures/Lecture2_qualityControl_artefactRemoval_DB.pdf)
    + [Quality Control and artefacts (practical)](Introduction/SS_DB/Materials/Practicals/Practical1_qualityControl_artefactRemoval_DB.pdf)
    + [Short Read Alignment (lecture)](Introduction/SS_DB/Materials/Lectures/Lecture3_ShortRead_Alignment_SS.pdf)
    + [Short Read Alignment (practical)](Introduction/SS_DB/Materials/Practicals/Practical2_Sequence_Alignment_SS.html)
    
# Day 2 (July 24th)

- 09:00 - 09:30; [Introduction to RNA-seq](RNASeq/slides/rnaSeq_Sept2017.pdf)
- 09:30 - 11:00; 
  + [Counting](RNASeq/count.nb.html)
  + [Importing and QC of RNA-seq data](RNASeq/rna-seq-preprocessing.nb.html)
  + [source file](RNASeq/slides/rna-seq-preprocessing.nb.Rmd)
- 11:00 - 12:30 Linear models & differential expression
- 12:30 - 13:30; LUNCH
- 13:30 - 15:00; Linear models & differential expression
- 15:00 - 17:00 Differential expression analysis with DESeq2
  + [Slides](RNASeq/slides/LinearModels.pdf)
  + [Practical](DifferentialExpression/rna-seq-de.nb.html)
  + [Source](DifferentialExpression/rna-seq.Rmd)

# Day 3 (July 25th)

- 09:30 - 11:00; [Annotation and Visualisation of Differential Expression](RNASeq/rna-seq-annotation-visualisation.nb.html)
  + [source file](RNASeq/rna-seq-annotation-visualisation.nb.Rmd)
- 11:00 - 12:30; [Gene set analysis and Gene Ontology testing](RNASeq/rna-seq-gene-set-testing.nb.html)
  + [source file](RNASeq/rna-seq-gene-set-testing.Rmd)
- 12:30 - 13:30; LUNCH
- 13:30 - 16:30; [Introduction to single-cell sequencing](SingleCell/index.html) NB We do not have sufficient time
to teach this entire course in half a day. However, some concepts are covered in the Bulk RNASeq course and we provide
access to the full materials. We will teach topics such as PCA which should be of interest even to those not interested in single-cell work.

# Day 4 (July 26th)

- 09:30 - 10:00; [Introduction to ChIP-Seq (lecture)](ChIP/Materials/Lectures/Lecture4_Introduction_to_ChIP-seq_and_ATAC-seq_SS.pdf). 
- 10:00 - 11:30; Peak calling and Visualisation  
  + [Peak Calling (lecture) ](ChIP/Materials/Lectures/Lecture5_Peak_Calling_SS.pdf). 
  + [Peak Calling (practical) ](ChIP/Materials/Practicals/Prctical4_PeakCalling_SS.pdf). 
  + [Viewing ChIP-seq results in IGV (practical) ](ChIP/Materials/Practicals/Practical3_IGV_DB.pdf)
- 11:30 - 12:30;
  + [Quality control methods for ChIP-seq (lecture)](ChIP/Materials/Lectures/Lecture6_chipqc_DB.pdf)
  + [Quality control methods for ChIP-seq (practical)](ChIP/Materials/Practicals/Practical5_chipqc_DB.pdf)
- 12:30 - 13:30; LUNCH
- 13:30 -14:00;
  + [Finish off: Quality control methods for ChIP-seq (practical)](ChIP/Materials/Practicals/Practical5_chipqc_DB.pdf)
- 14:00 - 15:30; 
 + [Useful software for the analysis of genomic data (lecture)](ChIP/Materials/Lectures/Lecture7_Useful_software_utilities_for_computational_genomics_SS.pdf)
 + [Useful software for the analysis of genomic data (practical)](ChIP/Materials/Practicals/Practical6_Useful_Utilities_for_Genomics.pdf)
- 15:30 - 17:00;
  + [Downstream analysis of ChIP-seq (lecture)](ChIP/Materials/Lectures/Lecture_8_and_9_Downstream_Analysis_of_ChIPseq_SS.pdf)
  + [Downstream analysis of ChIP-seq (practical)](ChIP/Materials/Practicals/Practical7_DownStreamAnalysis.pdf)
  
# Day 5 (July 27th)

- 09:30 - 12:30; 
- 09:30 - 10:30;  
  + [Identifying direct targets of Transcription Factors (practical)](ChIP/Materials/Practicals/Practical8_Rcade_SS.pdf)
- 10:30 - 11:30;  
  + [Differential binding analysis (lecture)](ChIP/Materials/Lectures/Lecture10_Differential_binding.pdf)
  + [Differential binding analysis (practical)](ChIP/Materials/Practicals/Practical9_diffbind_DB.pdf)
- 11:30 - 12:30; 
  + [Introduction to Epigenomics and Chromatin Interactions (lecture)](ChIP/Materials/Lectures/Lecture11_Intro_to_Epigenomics_SS.pdf)
- 12:30 - 13:30; LUNCH
  
<!--
## Data
- Mouse mammary data (counts): [https://figshare.com/s/1d788fd384d33e913a2a](https://figshare.com/s/1d788fd384d33e913a2a)
-->
