library(VariantAnnotation)
library(dplyr)
library(GenomicRanges)

annotate_genes <- function(vcf.o, genc){
  ## Let's say we are only interested in SVs affecting coding regions
  genc = subset(genc, type=='CDS')

  gene.ol = suppressWarnings(findOverlaps(vcf.o, genc)) %>% as.data.frame %>%
    mutate(gene_id=genc$gene_id[subjectHits],
           gene_name=genc$gene_name[subjectHits],
           type=genc$type[subjectHits]) %>%
    group_by(queryHits) %>%
    summarize(gene_name=paste(unique(sort(gene_name)), collapse='|'))

  ## IDEA: we could only keep SVs that overlap a gene to reduce the amount
  ## of variants to analyze in the next modules
  vcf.o = vcf.o[gene.ol$queryHits]
  
  ## IDEA: maybe it will be easier later if we duplicate SVs that overlap multiple genes
  ## i.e. one record for each SV-gene pair. Then the GENE field will contain only one gene
  ## name which might be easier to filter later.
  ## -> for now the new GENE field will have a list of all the genes overlapped, separated by '|'
  
  ## add field to VCF object
  hh = S4Vectors::DataFrame(Number='1', Type='String', Description='Overlapped gene(s)')
  rownames(hh) = 'GENE'
  info(header(vcf.o)) = rbind(info(header(vcf.o)), hh)
  info(vcf.o)$GENE = gene.ol$gene_name
  
  return(vcf.o)
}
