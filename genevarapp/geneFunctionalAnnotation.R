#! /usr/bin/env Rscript

###################################################################################
# File: geneFunctionalAnnotation.R
# Aim: Take the list of genes and annoate them based on human disease ontology & diesease gene-network and pathways like Reactome and KEGG  
# Author: Rupesh Kesharwani
# Last update: June 3, 2021
# Copyright (c) 2021 Kesharwani RK
###################################################################################

options(stringsAsFactors = F)
suppressWarnings(suppressMessages(library(clusterProfiler, quietly = T)))
suppressWarnings(suppressMessages(library(pathview, quietly = T)))
suppressWarnings(suppressMessages(library(org.Hs.eg.db, quietly = T)))
suppressWarnings(suppressMessages(library(enrichplot, quietly = T)))
suppressWarnings(suppressMessages(library(DOSE, quietly = T)))
suppressWarnings(suppressMessages(library(ggnewscale, quietly = T)))
suppressWarnings(suppressMessages(library(cowplot, quietly = T)))
suppressWarnings(suppressMessages(library(tidyverse, quietly = T)))
suppressWarnings(suppressMessages(library(plyr)))
suppressWarnings(suppressMessages(library(ReactomePA)))
suppressWarnings(suppressMessages(library(reactome.db)))
suppressWarnings(suppressMessages(library(KEGG.db)))

args = commandArgs(TRUE)

if(length(args) < 3)
{
 stop("#ERROR! please supply all inputs.
\nUSAGE: Rscript ./geneFunctionalAnnotation.R <listofgene> <pvalueCutoff> <GeneIDtype>\n")
}

listofgene = args[1] #list of genes
pvalueCutoff = as.numeric(args[2]) #the minimum p value for enriched GO and pathways
GeneIDtype = as.character(args[3]) #ID should be one of: SYMBOL,REFSEQ,ENSEMBL (## case sensitive)

# listofgene = "listofENSEMBLID.txt"
# pvalueCutoff = 0.05
# GeneIDtype = "ENSEMBL"

## input list of genes
listofgene <- read.table(listofgene, header = F)
listofgene <- listofgene[,1]

## Covert list of genes to entrezID
genelist2convertID <- function(listofgene, fromType, toType){
  id <- suppressMessages(suppressWarnings(bitr(listofgene, fromType=fromType, toType=toType, OrgDb="org.Hs.eg.db", drop = TRUE)))
  ### convert as charecter vector
  id2 = as.character(id[,2])
  ## to strip the NA values
  x <- id2[!is.na(id2)] ## entrez id
  x<- unique(x)
  return(x)
}

listEntrezID<-genelist2convertID(listofgene = listofgene, fromType = GeneIDtype, toType = 'ENTREZID')

## function to create GO/DO enriched output (txt)
print.go.out <- function(goResults, filename) {
  go_rich <- as.data.frame(goResults)
  go_empty <- is.data.frame(go_rich) && nrow(go_rich)==0
  if(go_empty==TRUE) {
    message("Seems some of Your ontology enrichment (DO/pathways) results are empty, please try with lower pvalue and qvalue Cutoff.\n")
  } else {
    ## enriched GO+Description+GenERatio+pvalue
    go <- suppressWarnings(as.data.frame(goResults)[,c(1:3,5:7)])
    #   #order by padj
    go <- go[with(go, order(p.adjust)), ]
    colnames(go)[1] <- "OntologyID"
    ##rownames(go) <- NULL
    write.table(as.data.frame(go), paste(filename,"Description.tsv", sep = "_"),quote=F, sep = "\t",row.names = F)
    ## enriched GO and their genes
    ## unlist object
    GO.list <- ldply(goResults@geneSets, data.frame)
    colnames(GO.list) <- c("EnrichedDO","AffliatedEntrezGeneID")
    write.table(GO.list, paste(filename,"enrichedGenes.tsv", sep = "_"),quote=F,row.names=F,sep = "\t")
  }
}

## gene Enrichment analysis (only Disease ontology as this is clinical data)
enrichedGene <- function (listEntrezID, showCategory=20){
	## enrichent analysis
  endgn <- DOSE::enrichDGN(listEntrezID, readable = T, pvalueCutoff = pvalueCutoff)
	enDO <- DOSE::enrichDO(listEntrezID, readable = T, pvalueCutoff = pvalueCutoff, qvalueCutoff = 0.1)
	## pathway analysis
	reactomePath <- enrichPathway(listEntrezID, pvalueCutoff = pvalueCutoff, readable = T)
	keggpath <- enrichKEGG(gene= listEntrezID, pvalueCutoff = pvalueCutoff, use_internal_data = T)
	## enrich plots
	bar_plot <- barplot(enDO, showCategory=showCategory) + ggtitle("Higher level of Disease Ontology")
	dot_plot <- enrichplot::dotplot(endgn, showCategory=showCategory, font.size = 8) + ggtitle("Disease associations from DisGeNET")
	edox <- DOSE::setReadable(endgn, 'org.Hs.eg.db', 'ENTREZID')
	cnet_plot <- cnetplot(edox, colorEdge = TRUE, categorySize="pvalue") + ggtitle("Top 5 Category")
	cnet_plot <- cnet_plot+theme(plot.title = element_text(hjust=0.7))
	anno_plot <- cowplot::plot_grid(bar_plot, dot_plot, cnet_plot, label_size = 8, ncol=3)
	## pathway plots
	ed1p <- pairwise_termsim(reactomePath)
	ed1 <- enrichplot::emapplot(ed1p)
	ed2p <- pairwise_termsim(keggpath)
	ed2 <- enrichplot::emapplot(ed2p)
	path_plot <- cowplot::plot_grid(ed2, ed1, ncol=2, labels=c("KEGGPathwasy", "ReactomePathways"))
	reactomeplot <-enrichplot::dotplot(reactomePath, showCategory=20) + ggtitle("ReactomePathway")
	keggplot <-enrichplot::dotplot(keggpath, showCategory=20) + ggtitle("KEGGPathway")
	comp_path <-cowplot::plot_grid(keggplot,reactomeplot, ncol=2)
	results <- list(endgn, enDO, reactomePath, keggpath, anno_plot, path_plot, comp_path)
	return(results)
}

## run enrich and plot
p <- enrichedGene(listEntrezID=listEntrezID)
## plot disease ontology
ggsave('genesDiseaseOntology.png', width = 16, height = 10, plot=p[[5]])
ggsave('PathwayNetwork.png', width = 16, height = 10, plot=p[[6]])
ggsave('PathwaysCompare.png', width = 16, height = 10, plot=p[[7]])
## print out
print.go.out(goResults = p[[1]], filename = "DGN")
print.go.out(goResults = p[[2]], filename = "DO")
print.go.out(goResults = p[[3]], filename = "ReactomePathway")
print.go.out(goResults = p[[4]], filename = "KEGGPathway")

##END
