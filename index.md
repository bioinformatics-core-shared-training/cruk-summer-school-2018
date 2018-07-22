# CRUK Bioinformatics Summer School 2018: Functional Genomics

**23rd - 27th July 2018: Craik-Marshall Room, Downing Site, University of Cambridge**

## Overview

Functional genomics looks at the dynamic aspects of how the genome functions within cells,
particularly in the form of gene expression (transcription) and gene regulation. This workshop surveys
current methods for functional genomics using high-throughput technologies. 

High-throughput technologies such as next generation sequencing (NGS) can routinely produce massive amounts of data. However, such datasets pose new challenges in the way the data have to be analyzed, annotated and interpreted which are not trivial and are daunting to the wet-lab biologist. This course covers state-of-the-art and best-practice tools for bulk RNA-seq and ChIP-seq data analysis, and will also introduce approaches for analysing data arising from single-cell RNA-seq studies.

### Audience

Enthusiastic and motivated wet-lab biologists who want to gain more of an understanding of NGS data and eventually progress to analysing their own data

### Pre-requisites

**The course will include a great deal of hands-on work in R and at the command line. In order for you to make the most of the course we strongly recommend that you take an introductory course, or have sufficient experience in the following areas:**

- R
- Unix
- Introductory statistics

**More specific requirements and references can be found [here](http://www.cruk.cam.ac.uk/bioinformatics-summer-school-prerequisites)**


### Instructors

- [Mark Fernandes (CRUK CI)](http://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
- [Shamith Samarajiwa (MRC CU)](http://www.mrc-cu.cam.ac.uk/research/Shamith-Samarajiwa-folder)
- [Dora Bihary (MRC CU)](http://www.mrc-cu.cam.ac.uk/research/Shamith-Samarajiwa-folder)
- [Ashley Sawle (CRUK CI)](http://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
- [Abigail Edwards (CRUK CI)](http://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
- [Alistair Martin (CRUK CI)](http://www.cruk.cam.ac.uk/research-groups/caldas-group)
- [Stephane Ballereau (CRUK CI)](http://www.cruk.cam.ac.uk/core-facilities/bioinformatics-core)
- [Michael Morgan (CRUK CI)](http://www.cruk.cam.ac.uk/).  


### Aims
During this course you will learn about:-

- To provide an understanding of how aligned sequencing reads, genome sequences and genomic regions are represented in R.
- To encourage confidence in reading sequencing reads into R, performing quality assessment and executing standard pipelines for (bulk) RNA-Seq and ChIP-Seq analysis 
- Analysis of transcription factor (TF) and epigenomic (histone mark) ChIP-seq data 
- Recent advances in single-cell sequencing

### Objectives
After the course you should be able to:-

- Know what tools are available in Bioconductor for HTS analysis and understand the basic object-types that are utilised.
- Process and quality control short read sequencing data 
- Given a set of gene identifiers, find out whereabouts in the genome they are located, and vice-versa 
- Produce a list of differentially expressed genes from an RNA-Seq experiment.
- Import a set of ChIP-Seq peaks and investigate their biological context.
- Appreciate the differences between bulk and single-cell RNA-seq analyses, and why the same methodologies might not be applicable

# Day 0 (July 22nd )

**SOCIAL
18:00 - ..
Informal get-together at The Grain and Hop Store (close to accommodation in Downing College)
Join us for a drink and dinner (self-paying), and to meet your colleagues for the next few days
http://www.grainandhopstore-cambridge.co.uk/ 
Note that the Training Room in Craik-Marshall building (1st Floor) will be open from 9am. **

# Day 1 (July 23rd)

- 09:30 [Course Introduction](Introduction/Session1-intro.html)
- 09:30 - 10:30; [Introduction to Functional Genomics](Introduction/Functional_Genomics_Overview.pdf)
- 10:30 - 12:30; [Introduction (Recap) of R and Bioconductor](Introduction/bioc-intro.nb.html)
  + [source](Introduction/bioc-intro.Rmd)
- 12:30 - 13:30; LUNCH
- 13:30 - 14:30 [Principles of Experimental Design](Introduction/Experimental_Design.pdf)
- 14:30 - 17:00;
- Data processing for Next Generation Sequencing
  + Lecture 1: [Introduction to next generation sequencing](Introduction/SS_DB/Materials/Lectures/L1_Introduction_to_Next_Generation_Sequencing_SS.pdf) (2.30- 2.45pm)
  + Lecture 2: [Brief introduction to file formats](Introduction/SS_DB/Materials/Lectures/L2_fileFormats_DB.pdf) (2.45- 3.00pm)
  + Lecture 3: [Quality control and artefact removal](Introduction/SS_DB/Materials/Lectures/L3_qualityControl_artefactRemoval_DB.pdf) (3.00- 3.45pm)
  + Practical 1: [learn to use FastQC and Cutadapt](ChIPSeq/Practicals/Practical1_fastQC_DB.html) (20 min) on a sample dataset
  + Lecture 4: [Short read alignment and Quality Control](Introduction/SS_DB/Materials/Lectures/L4_Short_Read_Alignment_to_a_Reference-Genome_SS.pdf) (3.45-5.00pm)
  + Practical 2: [Alignment of a ChIP-seq dataset to a reference genome using BWA OR Bowtie2 and a RNA-seq dataset to STAR](ChIPSeq/Practicals/Practical2_alignment_SS.html (45 min)
    
# Day 2 (July 24th)

**Please note we use several Rstudio Notebook html files as material for the RNAseq course. To obtain the source code
(.Rmd file) you can simply click on the code button in the top right-hand corner).**

- 09:00 - 09:30; [Introduction to RNA-seq](RNASeq2018/html/00_Introduction_to_RNAseq_Analysis.html)
- 09:30 - 11:00; 
  + [Counting](RNASeq2018/html/01_Read_Counts_with_Subread.html)
  + [Importing and QC of RNA-seq data](RNASeq2018/html/02_Preprocessing_Data.nb.html)
- 11:00 - 12:30 Linear models & differential expression
  + [Slides](RNASeq2018/slides/LinearModels.pdf)
  + [Linear models html nb](RNASeq2018/html/03_Linear_Models.nb.html)
- 12:30 - 13:30; LUNCH
- 13:30 - 15:00; Linear models & differential expression
- 15:00 - 17:00 
  + [Differential expression analysis with DESeq2](RNASeq2018/html/04_DE_analysis_with_DESeq2.nb.html)

# Day 3 (July 25th)

- 09:30 - 11:00; [Annotation and Visualisation of Differential Expression](RNASeq2018/html/05_Annotation_and_Visualisation.nb.html)
- 11:00 - 12:30; [Gene set analysis and Gene Ontology testing](RNASeq2018/html/06_Gene_set_testing.nb.html)

- 12:30 - 13:30; LUNCH
- 13:30 - 16:30; [Introduction to single-cell sequencing](SingleCell/index.html) NB We do not have sufficient time
to teach this entire course in half a day. However, some concepts are covered in the Bulk RNASeq course and we provide
access to the full materials. We will teach topics such as PCA which should be of interest even to those not interested in single-cell work.

**SOCIAL: Punting trip - leave from Mill Lane punting site at 18:00 (~10 min walk from Craik-Marshall)
[Google Map](https://www.google.com/maps/dir/Craik-Marshall+Building,+Cambridge/Scudamore's+Mill+Lane+Punting+Station,+Mill+Ln,+Cambridge+CB2+1RS/@52.2021771,0.1169979,17z/data=!3m1!4b1!4m14!4m13!1m5!1m1!1s0x47d8709788f7fb6b:0x117c3858ab077fb0!2m2!1d0.1224918!2d52.2019233!1m5!1m1!1s0x47d870a335f0f833:0xdb33dae59780b590!2m2!1d0.115892!2d52.201511!3e3) **

# Day 4 (July 26th)

- ChIP-seq data analysis
    + Lecture 5: [Introduction to ChIP-seq](ChIPSeq/Lectures//L5_Introduction_to_ChIP-seq_SS.pdf) (9.30-10.00pm)
    + Lecture 6: [Peak Calling](ChIPSeq/Lectures/L6_Peak_Calling_SS.pdf) (10.00-11.00pm)
    + Practical 3: [Peak calling using MACS2](ChIPSeq/Practicals/Practical3_peakcalling_SS.html) (30 min)
    + Lecture 7: [Differential binding analysis](ChIPSeq/Lectures/L7_Differential_binding_analysis_DB.pdf) (11.00-12.30pm)
    + Practical 4: [THOR (and Diffbind)](ChIPSeq/Practicals/Practical4_differentialBinding_DB.html) (20 min)
    + Lecture 8: [Quality control methods for ChIP-seq](ChIPSeq/Lectures/L8_Quality_control_methods_for_ChIP-seq_DB.pdf) (1 hr)
    + Practical 5: [ChIPQC package](ChIPSeq/Practicals/Practical5_ChIPQC_DB.html) (30 min)
    + Practical 6: [Integrative Genome Viewer](ChIPSeq/Practicals/Practical6_IGV_DB.html) (30 min)
    LUNCH (12.30-1.30pm)
    + Lecture 9: [Downstream analysis of ChIP-seq](ChIPSeq/Lectures/L9_Downstream_Analysis_of_ChIP-seq_data_SS.pdf) (1.30-3.15pm)
    + Practical 7: [Downstream analysis of ChIP-seq](ChIPSeq/Practicals/Practical7_Downstream_Analysis_of_ChIP-seq_SS.html) (30 min)
    + Practical 8: [Identifying direct targets of transcription factors with Rcade](ChIPSeq/Practicals/Practical8_Rcade_SS.pdf) (30 min)
    + Lecture 10:: [Useful software utilities for the analysis of genomic data](ChIPSeq/Lectures/L10_Software_utilities_for_computational_genomics_SS.pdf) (4.30-5.00pm)
 
**SOCIAL: Summer School evening meal & Reception at the nearby Downing College at 18:00 to 22:30.. If you do not wish to attend this meal (free to attendees) then please let us know ASAP. http://www.dow.cam.ac.uk  Downing College, Regent Street, Cambridge, CB2 1DQ**

# Day 5 (July 27th)

- 09:30 - 12:30; 
- ATAC-seq and Epigenomics
  +  Practical 9: [Useful software utilities for the analysis of genomic data](ChIPSeq/Practicals/Practical9_Useful_software_utilities_for_the_analysis_of_genomic_data_SS.html) (9.30-10.30am)
  +  Lecture 11 [ATAC-seq data analysis](ChIPSeq/Lectures/L11_Introduction to ATAC-seq_SS.pdf) (10.30-11.30am)
  +  Practical 10: [ATAC-seq analysis](ChIPSeq/Practicals/) (30 min)
  +  Lecture 12 [Introduction to Epigenomics and Chromatin Interactions](ChIPSeq/Lectures/L12_Introduction_to_Epigenomics_SS.pdf) (11.30-12.30)

- 12:30 - 13:30; LUNCH
- Safe Journey home (Thank you for participating in the Summer School)
  
<!--
## Data
- Mouse mammary data (counts): [https://figshare.com/s/1d788fd384d33e913a2a](https://figshare.com/s/1d788fd384d33e913a2a)
-->
