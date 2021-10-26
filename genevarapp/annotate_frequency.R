library(VariantAnnotation)
library(dplyr)
library(GenomicRanges)
library(sveval)

## USES VCF FIELDS: SVLEN, SVTYPE

annotate_frequency <- function(vcf.o, gnomad){

  svs.gr = rowRanges(vcf.o)
  svs.gr$type = info(vcf.o)$SVTYPE
  svs.gr$size = unlist(lapply(info(vcf.o)$SVLEN, function(x) abs(x[1])))
  names(svs.gr) = NULL

  af.ol = suppressWarnings(svOverlap(svs.gr, gnomad, min.ol=.1, max.ins.dist=100)) %>%
    as.data.frame %>%
    mutate(freq=gnomad$AF[subjectHits]) %>%
    group_by(queryHits) %>%
    summarize(freq=max(freq))

  afs = rep(0, length(vcf.o))
  afs[af.ol$queryHits] = af.ol$freq
  
  ## IDEA: we could only keep SVs that are rare to reduce the amount
  ## of variants to analyze in the next modules
  
  ## add field to VCF object
  hh = S4Vectors::DataFrame(Number='1', Type='Float', Description='Allele frequency')
  rownames(hh) = 'AF'
  info(header(vcf.o)) = rbind(info(header(vcf.o)), hh)
  info(vcf.o)$AF = afs
    
  return(vcf.o)
}
