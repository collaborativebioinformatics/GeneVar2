# GeneVar2: Effortless SV annotation and interpretation 

![]()![VectorEPS_ByTailorBrands_v2](https://user-images.githubusercontent.com/41301333/136847583-fa82b8ec-6762-461f-be20-b4fec6d23561.jpg)


## Contributors

1. Tim Hefferon `(Leader & Liaison)`
2. Ahmad Al Khleifat `(Data Guru and Writer)`
3. Rupesh Kesharwani `(Sysadmin and code developer)`
4. Divya Kalra `(Sysadmin and code developer)`
5. Kimberly Walker `(Writer)`
6. Priya Lakra `(Code developer)`
7. Jianzhi(Quentin) Yang `(QC checker)` 
8. Jean Monlong `(Sysadmin and code developer)`
9. David Henke `(Writer and Guide)`
10. Weiyu Zhou (`App developer`)

## Motivation

Structural variants represent a major difference between individuals in health and disease. For those outside the immediate field of genetics, a group that includes researchers and clinicians the interpretation of findings is particularly challenging.  GeneVar2 is an extremely fast and computationally efficient platform for analysis, visualization, and interpretation of structural variation data. It is designed to provide a powerful and easy-to-use tool for applications in biomedical research and diagnostic medicine, at no computational cost. Its comprehensive approach brings the analyses of structural variation within the reach of non-specialist laboratories, and those from centres with limited funding available.


## Goals

[GeneVar](https://github.com/collaborativebioinformatics/GeneVar) is an open access, gene centric data browser for SV analysis. GeneVar takes as input a gene name or ID and produces a report that informs the user of all SVs overlapping the gene and any non-coding regulatory elements affecting expression of the gene. [Clinical SV](https://github.com/collaborativebioinformatics/clinical_SVs) is an open access software that can annotate vcf files with clinically relavant information as well as provide useful visualizations such as disease ontology plots.

GeneVar2 is the integration of these two apps which work together to facilitate reporting of structural variations data. GeneVar2 tool is intended to have a clinical focus, informing the interpretation of SV pertaining to a gene name. In addition, GeneVar2 gives the user the option to upload genotyping data and produces a report, file, and genome browser session that informs the user of all structural variants overlapping the gene, including any non-coding regulatory elements affecting expression of the gene.


## Description

The aim of this project is to merge the functionality of GeneVar and Clinical_SV into one new application, GeneVar2. GeneVar is a gene centric data browser which is great to review a small list of genes individually for the browser allows for in-depth analysis of SVs that overlap a gene of interest. However, many users will typically have variant caller files (VCFs) as output from analysis pipelines.  To better accomodate this use-case, we are combining GeneVar with Clinical_SV which already encompasses the ability to upload and annotat SV vcfs.  In addition, Clinical_SV produces helpful visualizations of Disease Ontology and enrichment pathway analysis based on SV types.


## Overview Diagram

![](GeneVar2_workflow_v2.png)


## How it works

GeneVar2 is a web page application.  Users have two ways of interacting with the tool depending on what their input is.

1) Enter individual genes: After entering the gene name (HGNC, Ensembl gene (ENSG), or transcript (ENST) identifier) in the search box on the homepage, you will be directed to the gene-specific page containing:

    - Gene-level summary with number of SVs, number of clinival SVs or SVs overlapping clinical SNVs.
    - Links to the gene's page on OMIM, GTEx, gnomAD.
    - A dynamic table with the annotated variants overlapping the gene.
    - A graph with the distribution of the allele frequency for variants matched with gnomAD-SV (50% reciprocal overlap).

    The profile of the SV to consider, such as type and size range, can be specified on the side bar. 
    Each column in the dynamic table can be "searched" into or reorder dynamically. 
    All data used by the app will be available for download in tab-delimited files. 
    By default, allele frequency is reported based on gnomAD genomes and exomes.

2) Upload a vcf file: After users upload their own structural variant vcf file, the application will annotate each SV with the following annotation:

   - Allele Frequency: For variants found in gnomAD-SV, allele frequency for available populations will be annotated.
   - ClinVar Information: Pathogenic SVs, SNVs and Indels from ClinVar that overlap with called SVs will be annotated.

In addition a number of tables are provided to download including Disease Ontology details table and Clinical rank table.

Finally, a number of plots are also provided to download including KEGG pathway, Disease Ontology, and Disease associations from DisGeNET.


## Installation and Quick Start
For users interested in annotating their own vcf files without using the web application, the following R scripts are available.  First install the required packages in R.

```r
install.packages("easypackages")
if(!"easypackages" %in% row.names(installed.packages())){
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  library(easypackages, character.only = TRUE, quietly = TRUE)
}
pkgs=c("shiny","tippy","shinythemes", "tidyverse", "tidygraph", "clusterProfiler","org.Hs.eg.db","DOSE","ggnewscale","cowplot","tidyverse","plyr","ReactomePA","reactome.db","reactome.db", 
"KEGG.db","enrichplot","dplyr","GenomicRanges", "rtracklayer", "VariantAnnotation", "tidyr")
suppressWarnings(suppressMessages(easypackages::packages(pkgs, prompt = FALSE)))
```

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install('jmonlong/sveval')
```

#### Note: For any issues encountered installing these packages, please review: https://www.bioconductor.org/packages/release/bioc

### SV calling

If users first need to call SVs on their samples, the developers recommend Parliament2.  Parliament2 runs a combination of tools to generate structural variant calls on whole-genome sequencing data. It can run the following callers: Breakdancer, Breakseq2, CNVnator, Delly2, Manta, and Lumpy. Because of synergies in how the programs use computational resources, these are all run in parallel. Parliament2 will produce the outputs of each of the tools for subsequent investigation.  See https://github.com/fritzsedlazeck/parliament2 for further details.

### Annotation of SVs in R

The different modules of the annotation are written as functions and saved in separate files.
Then the master annotation script can read a VCF, *source* these functions and use them to annotate the SVs. 
See the current master annotation script [`annotate_vcf.R`](R/annotate_vcf.R) and the different scripts *source*d inside.


```r
Rscript annotate_vcf.R test.input.vcf annotation_data.RData test.output.vcf test.output.csv
```

### Gene/Diesease ontology and Pathways analysis and plotting

This module supports the enrichment analysis of Disease Ontology (DO) (Schriml et al. 2011), Network of Cancer Gene (A. et al. 2016) and Disease Gene Network (DisGeNET) (Janet et al. 2015). In addition, several visualization methods were provided by enrichplot to help interpreting enrichment and disease ontology results.


```r
Rscript GeneAnnotationFromCSV.R test.output.vcf 0.05 (pvalueCutoff) DUP (svtype) chr12 (Chromosome)
```


## Test data

test.input.vcf

## Results
The following are example result files that can be generated by either the R scripts described above or by the web application.

1A. Annotated SV VCF (test.output.vcf) 

1B. CSV (test.output.csv) has several fields of information, explained in the table below, including clinically relevant ID and SV RANK

| name       | description                                                           |
|------------|-----------------------------------------------------------------------|
| gene       | names of genes overlapped, separated by `\|`                          |
| variant_id | SV ID                                                                 |
| chr        | chromosome name                                                       |
| start      | start position                                                        |
| end        | end position                                                          |
| size       | size of the SV in bp                                                  |
| frequency  | allele frequency                                                      |
| svtype     | type of SV. E.g. DEL, DUP, INS, ...                                   |
| clinsv     | dbVar accession IDs of matching known clinical SVs (separated by `\|` |
| clinrk     | clinical importance rank, for example to select top 5 SVs             |


![](Results-Table1.png)


2. The results obtained from previous step are used to generate the gene and disease ontology analysis. The user can specify a chr and cvtype to plot.

## An example of chr12 and DUP svtype


![](genesDiseaseOntology.png)


![](PathwaysCompare.png)


![](Table-results.png)


![freq_plot_perchr_plot](https://user-images.githubusercontent.com/41301333/137140949-aa41ab63-f3ed-4140-bd33-67c97fecd6f5.jpeg)



## Citation

GeneVar2 is available on GitHub (https://github.com/collaborativebioinformatics/GeneVar2). The repository provides detailed instructions for tool usage and installation. 


## References
- GeneVar: https://github.com/collaborativebioinformatics/GeneVar
- Clinical SV: https://github.com/collaborativebioinformatics/clinical_SVs 
- Parliament2: https://github.com/fritzsedlazeck/parliament2
