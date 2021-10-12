# GeneVar2: Effortless SV annotation and interpretation 

![]()![VectorEPS_ByTailorBrands_v2](https://user-images.githubusercontent.com/41301333/136847583-fa82b8ec-6762-461f-be20-b4fec6d23561.jpg)


## Contributors

1. Tim Hefferon `(Leader & Liaison)`
2. Ahmad Al Khleifat `(Writer and Guide)`
3. Rupesh Kesharwani `(Sysadmin and code developer)`
4. Divya Kalra `(Sysadmin and code developer)`
5. Kimberly Walker `(Writer)`
6. Priya Lakra `(Code developer)`
7. Jianzhi(Quentin) Yang `(QC checker)` 
8. Jean Monlong `(Sysadmin and code developer)`
9. David Henke `(Writer and Guide)`
10. Weiyu Zhou (`App developer`)

## Goals

GeneVar is an open access, gene centric data browser for SV analysis. GeneVar takes as input a gene name or ID and produces a report that informs the user of all SVs overlapping the gene and any non-coding regulatory elements affecting expression of the gene. Clinical_SV is an open access software that can annotate vcf files with clinically relavant information as well as provide useful visualizations such as disease ontology plots.

GeneVar2 is the integration of these two apps which work together to facilitate reporting of structural variations data. GeneVar2 tool is intended to have a clinical focus, informing the interpretation of SV pertaining to a gene name. In addition, GeneVar2 gives the user the option to upload genotyping data and produces a report, file, and genome browser session that informs the user of all structural variants overlapping the gene, including any non-coding regulatory elements affecting expression of the gene.


## Description

The aim of this project is to merge the functionality of GeneVar and Clinical_SV into one new application, GeneVar2. GeneVar is a gene centric data browser which is great to review a small list of genes individually for the browser allows for in-depth analysis of SVs that overlap a gene of interest. However, many users will typically have variant caller files (VCFs) as output from analysis pipelines.  To better accomodate this use-case, we are combining GeneVar with Clinical_SV which already encompasses the ability to upload and annotat SV vcfs.  In addition, Clinical_SV produces helpful visualizations of Disease Ontology and enrichment pathway analysis based on SV types.


## Overview Diagram

![](GeneVar2_workflow_v2.png)


## How it works

GeneVar2 is a web page application.  Users have two ways of interacting with the tool depending on what their input is.

1) Enter individual genes: After entering the gene name (HGNC, Ensembl gene (ENSG), or transcript (ENST) identifier) in the search box on the homepage, you will be directed to the gene-specific page containing:

    Gene-level summary with number of SVs, number of clinival SVs or SVs overlapping clinical SNVs.
    Links to the gene's page on OMIM, GTEx, gnomAD.
    A dynamic table with the annotated variants overlapping the gene.
    A graph with the distribution of the allele frequency for variants matched with gnomAD-SV (50% reciprocal overlap).

The profile of the SV to consider, such as type and size range, can be specified on the side bar. Each column in the dynamic table can be "searched" into or reorder dynamically. All data used by the app will be available for download in tab-delimited files. By default, allele frequency is reported based on gnomAD genomes and exomes.

2) Upload a vcf file:After uploading their own structural variant vcf file, the application will annotate each SV with the following annotation:

A) Allele Frequency: For variants found in gnomAD-SV, allele frequency for available populations will be annotated.
B) ClinVar Information: Pathogenic SVs, SNVs and Indels from ClinVar that overlap with called SVs will be annotated.

In addition a number of tables will be provided for download.  These tables include:
need a list of tables so I can describe them here.

Additionally the following plots will be available for download:
Need to add verbiage explaining graphs


## Installation

```r
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

## Note: Please have a look (https://www.bioconductor.org/packages/release/bioc) for any kind of issues.


## Quick Start

### SV calling

If needed SVs can be called using parliament2 (https://github.com/dnanexus/parliament2).

### Annotation of SVs in R

The different modules of the annotation are written as functions and saved in separate files.
Then the master annotation script can read a VCF, *source* these functions and use them to annotate the SVs. 
See the current master annotation script [`annotate_vcf.R`](R/annotate_vcf.R) and the different scripts *source*d inside.


```
Rscript annotate_vcf.R test.input.vcf annotation_data.RData test.output.vcf test.output.csv
```

### Gene/Diesease ontology and Pathways analysis and plotting

This module supports the enrichment analysis of Disease Ontology (DO) (Schriml et al. 2011), Network of Cancer Gene (A. et al. 2016) and Disease Gene Network (DisGeNET) (Janet et al. 2015). In addition, several visualization methods were provided by enrichplot to help interpreting enrichment and disease ontology results.


```
Rscript GeneAnnotationFromCSV.R test.output.vcf 0.05 (pvalueCutoff) DUP (svtype) chr12 (Chromosome)
```


## Test data

test.input.vcf

## Results

1A. Annotated SV VCF (test.output.vcf) 

1B. CSV (test.output.csv) with clinically relevant ID and RANK

The table contains the following columns.

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


2. The results obtained from previous step used to perform the gene and disease ontology analysis. The user can choose chr and cvtype to plot these.

## An example of chr12 and DUP svtype


![](genesDiseaseOntology.png)


![](PathwaysCompare.png)


![](Table-results.png)


## Citation

GeneVar2 is available on GitHub (https://github.com/collaborativebioinformatics/GeneVar2). The repository provides detailed instructions for tool usage and installation. 


## References

