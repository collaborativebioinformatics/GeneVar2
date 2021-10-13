#! /usr/bin/env Rscript

###################################################################################
# File: GeneAnnotationFromCSV.R
# Aim: Take the list of genes and annotate them based on human disease ontology & disease gene-network and pathways like Reactome and KEGG  
# Author: Rupesh Kesharwani
# Last update: Oct 12, 2021
# Copyright (c) 2021 Kesharwani RK
###################################################################################

##########################
##    Hackathon 2021    ##
##    2021 Oct 10-13    ##
##      GeneVar2        ##
##      SV vcf/csv      ## 
##########################

args = commandArgs(TRUE)

if(length(args) < 4)
{
  stop("#ERROR! please supply all inputs.
\nUSAGE: Rscript ./GeneAnnotationFromCSV.R <CSV output from annotate_vcf.R> <pvalueCutoff> <svtype> <chrom>\n")
}

csv_output = args[1] #csv file output from annotate_vcf.R (Clinical_SV package)
pvalueCutoff = as.numeric(args[2]) #the minimum p value for enriched GO and pathways
svtype = args[3] ## user defined SV type
chrom = args[4] ## user defined  SV chromosome

options(stringsAsFactors = F)
cat("Loading packages...\n")
suppressWarnings(suppressMessages(library(clusterProfiler, quietly = T)))
suppressWarnings(suppressMessages(library(pathview, quietly = T)))
suppressWarnings(suppressMessages(library(org.Hs.eg.db, quietly = T)))
suppressWarnings(suppressMessages(library(enrichplot, quietly = T)))
suppressWarnings(suppressMessages(library(DOSE, quietly = T)))
suppressWarnings(suppressMessages(library(ggnewscale, quietly = T)))
suppressWarnings(suppressMessages(library(cowplot, quietly = T)))
suppressWarnings(suppressMessages(library(tidyverse, quietly = T)))
suppressWarnings(suppressMessages(library(plyr, quietly = T)))
suppressWarnings(suppressMessages(library(ReactomePA, quietly = T)))
suppressWarnings(suppressMessages(library(reactome.db, quietly = T)))
suppressWarnings(suppressMessages(library(KEGG.db, quietly = T)))
suppressWarnings(suppressMessages(library(dplyr, quietly = T)))
suppressWarnings(suppressMessages(library(tidyr, quietly = T)))
cat("done.\n")

# ## Example:
# csv_output = "test.output.csv"
# pvalueCutoff = 0.1
# svtype = "DEL"
# chrom = "chr7"

GeneIDtype = "SYMBOL" #ID should be one of: SYMBOL,REFSEQ,ENSEMBL (## case sensitive)

cat("#####################################\n")
cat("Analysis started..\n")
date()
cat("#####################################\n")

cat("Your inputs...\n")
cat("CSV as input:", csv_output, "\n")
cat("P-value cutoff:", pvalueCutoff, "\n")
cat("svtype:", svtype, "\n")
cat("chromosome:", chrom, "\n")


## input CSV file
outcvs <- read.csv(csv_output, header = T)
## select only gene, chr and svtype
myout <- outcvs %>% dplyr::select(gene, chr, svtype)
## filter by chr as well as by svtype
user_filt <- dplyr::filter(myout, svtype == svtype & chr == chrom) %>% select(gene)
user_genes <- user_filt %>% tidyr::separate_rows(gene, sep = '\\|')
listgene<-as.data.frame(user_genes$gene)
colnames(listgene) <- NULL
listofgene <- listgene[,1]

## Covert list of genes to entrezID
genelist2convertID <- function(listofgene, fromType, toType){
  id <- suppressMessages(suppressWarnings(bitr(listofgene, fromType=fromType, 
                                               toType=toType, OrgDb="org.Hs.eg.db", 
                                               drop = TRUE)))
  ### convert as charecter vector
  id2 = as.character(id[,2])
  ## to strip the NA values
  x <- id2[!is.na(id2)] ## entrez id
  x<- unique(x)
  return(x)
}

cat("converting list of gene SYMBOL to ENTREZ ID...\n")
listEntrezID<-genelist2convertID(listofgene = listofgene, fromType = GeneIDtype, 
                                 toType = 'ENTREZID')

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

# ## gene Enrichment analysis (only Disease ontology as this is clinical data)
# enrichedGene <- function (listEntrezID, showCategory=20){
#   cat("## enrichent analysis...\n")
#   endgn <- DOSE::enrichDGN(listEntrezID, readable = T, pvalueCutoff = pvalueCutoff)
#   enDO <- DOSE::enrichDO(listEntrezID, readable = T, pvalueCutoff = pvalueCutoff, qvalueCutoff = 0.1)
#   cat("## pathway analysis...\n")
#   reactomePath <- enrichPathway(listEntrezID, pvalueCutoff = pvalueCutoff, readable = T)
#   keggpath <- enrichKEGG(gene= listEntrezID, pvalueCutoff = pvalueCutoff, use_internal_data = T)
#   cat("## enrich plots...\n")
#   bar_plot <- barplot(enDO, showCategory=showCategory) + ggtitle("Higher level of Disease Ontology")
#   dot_plot <- enrichplot::dotplot(endgn, showCategory=showCategory, font.size = 8) + ggtitle("Disease associations from DisGeNET")
#   edox <- DOSE::setReadable(endgn, 'org.Hs.eg.db', 'ENTREZID')
#   cnet_plot <- cnetplot(edox, colorEdge = TRUE, categorySize="pvalue") + ggtitle("Top 5 Category")
#   cnet_plot <- cnet_plot+theme(plot.title = element_text(hjust=0.7))
#   anno_plot <- cowplot::plot_grid(bar_plot, dot_plot, cnet_plot, label_size = 8, ncol=3)
#   cat("## pathway plots...\n")
#   ed1p <- pairwise_termsim(reactomePath)
#   ed1 <- enrichplot::emapplot(ed1p)
#   ed2p <- pairwise_termsim(keggpath)
#   ed2 <- enrichplot::emapplot(ed2p)
#   path_plot <- cowplot::plot_grid(ed2, ed1, ncol=2, labels=c("KEGGPathways", "ReactomePathways"))
#   reactomeplot <-enrichplot::dotplot(reactomePath, showCategory=20) + ggtitle("ReactomePathway")
#   keggplot <-enrichplot::dotplot(keggpath, showCategory=20) + ggtitle("KEGGPathway")
#   comp_path <-cowplot::plot_grid(keggplot,reactomeplot, ncol=2)
#   results <- list(endgn, enDO, reactomePath, keggpath, anno_plot, path_plot, comp_path)
#   return(results)
# }

results <- list()
## gene Enrichment analysis (only Disease ontology as this is clinical data)
enrichedGene <- function (listEntrezID, showCategory=20, chrom, svtype){
  ## enrichent analysis
  enDO <- DOSE::enrichDO(listEntrezID, readable = T, pvalueCutoff = pvalueCutoff, 
                         qvalueCutoff = 0.1, minGSSize = 1)
  results[[1]] <- enDO
  endgn <- DOSE::enrichDGN(listEntrezID, readable = T, pvalueCutoff = pvalueCutoff, 
                           minGSSize = 1, qvalueCutoff = 0.1)
  results[[2]] <- endgn
  ## enrich plots1
  if(as.double(summary(as.data.frame(enDO)$geneID)[1]) > 0) {
    bar_plot <- barplot(enDO, showCategory=showCategory) + 
      ggtitle("Higher level of Disease Ontology") 
    ggsave(paste(chrom, svtype,'genesDiseaseOntology_bar_plot.png', sep="_"), width = 16, 
           height = 10, bar_plot)
  } else {
    cat("Note: bar_plot can't be plotted as no enrichDO found. \n")
  }
  ## enrich plots2
  if(as.double(summary(as.data.frame(endgn)$geneID)[1]) > 0) 
     {
       dot_plot <- enrichplot::dotplot(endgn, showCategory=showCategory, font.size = 6) + ggtitle("Disease associations from DisGeNET")
       edox <- DOSE::setReadable(endgn, 'org.Hs.eg.db', 'ENTREZID')
       #cnet_plot <- cnetplot(edox, colorEdge = TRUE, categorySize="pvalue") + ggtitle("Top 5 Category")
       #cnet_plot <- cnet_plot+theme(plot.title = element_text(hjust=0.7))
       cnet_plot <- cnetplot(edox, color_category='firebrick', color_gene='steelblue', 
        categorySize="pvalue") + ggtitle("Top 5 Category")
       ggsave(paste(chrom, svtype,'genesDiseaseOntology_dot_plot.png', sep="_"), width = 14, 
              height = 8 , dot_plot)
       ggsave(paste(chrom, svtype,'genesDiseaseOntology_cnet_plot.png', sep="_"), width = 14, 
              height = 8 , cnet_plot)
  } else {
    cat("Note: dot_plot can't be plotted as no enrichDGN found. \n")
  }
  
  ## pathway analysis
  reactomePath <- enrichPathway(listEntrezID, pvalueCutoff = pvalueCutoff, 
                                readable = T, qvalueCutoff = 0.1, minGSSize = 1)
  results[[3]] <- reactomePath
  keggpath <- enrichKEGG(gene= listEntrezID, pvalueCutoff = pvalueCutoff, 
                         use_internal_data = T, minGSSize = 1, qvalueCutoff = 0.1)
  results[[4]] <- keggpath
  ## pathway plots1
  if(as.double(summary(as.data.frame(reactomePath)$geneID)[1]) > 0) {
    reactomeplot <-enrichplot::dotplot(reactomePath, showCategory=showCategory) + 
      ggtitle("ReactomePathway")
    ggsave(paste(chrom, svtype,'reactomeplot_dot_plot.png', sep="_"), width = 14, 
           height = 8, reactomeplot)
  } else {
    cat("Note: reactome_plot can't be plotted as no enrichPathway found. \n")
  }
  
  ## pathway plots2
  if(as.double(summary(as.data.frame(keggpath)$geneID)[1]) > 0) {
    keggplot <-enrichplot::dotplot(keggpath, showCategory=showCategory) + 
      ggtitle("KEGGPathway")
    ggsave(paste(chrom, svtype,'keggplot_dot_plot.png', sep="_"), width = 14, 
           height = 8, keggplot)
  } else {
    cat("Note: kegg_plot can't be plotted as no enrichKEGG found. \n")
  }
  #results <- list(endgn, enDO, reactomePath, keggpath, anno_plot, path_plot, comp_path)
  return(results)
}


cat("Running enrichment and ontology...\n")

## run enrich and plot
p <-enrichedGene(listEntrezID = listEntrezID, showCategory = 10, chrom = chrom, 
             svtype = svtype)


# cat("Plotting...\n")
# ## plot disease ontology
# ggsave(paste(chrom,svtype,'genesDiseaseOntology.png',sep="_"), width = 16, height = 10, plot=p[[5]])
# ggsave(paste(chrom,svtype,'PathwayNetwork.png',sep="_"), width = 16, height = 10, plot=p[[6]])
# ggsave(paste(chrom,svtype,'PathwaysCompare.png', sep="_"), width = 16, height = 10, plot=p[[7]])

cat("Generating results output...\n")
## print out
print.go.out(goResults = p[[1]], filename = paste(chrom,svtype,"DO",sep="_"))
print.go.out(goResults = p[[2]], filename = paste(chrom,svtype,"DGN", sep="_"))
print.go.out(goResults = p[[3]], filename = paste(chrom,svtype,"ReactomePathway", sep="_"))
print.go.out(goResults = p[[4]], filename = paste(chrom,svtype,"KEGGPathway", sep="_"))
##END
save(listEntrezID, p, print.go.out, enrichedGene, file = paste(chrom,svtype,"GeneOntology.RData", sep="_"))
cat("#Done, the results are generated in current dir.\n")

