library(VariantAnnotation)
library(dplyr)
library(GenomicRanges)
library(sveval)

## USES VCF FIELDS: SVLEN, SVTYPE

annotate_known_clinical_SVs <- function(vcf.o, clinsv){

  svs.gr = rowRanges(vcf.o)
  svs.gr$type = info(vcf.o)$SVTYPE
  svs.gr$size = unlist(lapply(info(vcf.o)$SVLEN, function(x) abs(x[1])))
  names(svs.gr) = NULL

  clins = rep(NA, length(vcf.o))

  clsig = c('Pathogenic', 'Pathogenic/Likely pathogenic', 'Likely pathogenic')
  clinsv = subset(clinsv, clinical_significance %in% clsig)
  
  clin.ol =  svOverlap(svs.gr, clinsv, min.ol=.5, max.ins.dist=100) %>% as.data.frame
  if(nrow(clin.ol)>0){
    clin.ol = clin.ol %>% mutate(dbvar_id=clinsv$dbvar_id[subjectHits]) %>%
      group_by(queryHits) %>%
      summarize(dbvar_id=paste(unique(sort(dbvar_id)), collapse='|'))
    clins[clin.ol$queryHits] = clin.ol$dbvar_id
  }

  ## add field to VCF object
  hh = S4Vectors::DataFrame(Number='1', Type='String', Description='IDs of matching known clinical SVs')
  rownames(hh) = 'CLINSV'
  info(header(vcf.o)) = rbind(info(header(vcf.o)), hh)
  info(vcf.o)$CLINSV = clins
    
  return(vcf.o)
}
