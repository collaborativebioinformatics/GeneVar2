library(dplyr)

annotate_clinical_score <- function(vcf.o, genc){
  clinscore = rep(0, length(vcf.o))

  ## +1 if known clinical SV
  clinscore = clinscore + ifelse(is.na(info(vcf.o)$CLINSV), 0, 1)
  
  ## +1 if never been seen in gnomAD-SV
  clinscore = clinscore + ifelse(info(vcf.o)$AF>0, 0, 1)

  ## +1 if gene is predicted to be loss-of-function intolerant
  sv.gene.pli = unlist(lapply(strsplit(info(vcf.o)$GENE, '\\|'),
                              function(genes) any(genes %in% gene.pli)))
  clinscore = clinscore + ifelse(sv.gene.pli, 1, 0)

  ## +1 if gene a known OMIM disease gene
  sv.gene.omim = unlist(lapply(strsplit(info(vcf.o)$GENE, '\\|'),
                              function(genes) any(genes %in% gene.omim)))
  clinscore = clinscore + ifelse(sv.gene.pli, 1, 0)

  nb.clinsv = unlist(lapply(strsplit(info(vcf.o)$CLINSV, '\\|'), length))
  nb.clinsv = ifelse(is.na(info(vcf.o)$CLINSV), 0, nb.clinsv)

  nb.gene = unlist(lapply(strsplit(info(vcf.o)$GENE, '\\|'), length))
  nb.gene = ifelse(is.na(info(vcf.o)$GENE), 0, nb.gene)

  df = tibble(idx=1:length(vcf.o), score=clinscore,
              nb.clinsv=nb.clinsv, nb.gene=nb.gene) %>%
    arrange(desc(clinscore), desc(nb.clinsv), desc(nb.gene)) %>%
    mutate(rk=1:n())
  
  clinrk = rep(NA, length(vcf.o))
  clinrk[df$idx] = df$rk
  
  ## add field to VCF object
  hh = S4Vectors::DataFrame(Number='1', Type='Integer', Description='Clinical rank')
  rownames(hh) = 'CLINRK'
  info(header(vcf.o)) = rbind(info(header(vcf.o)), hh)
  info(vcf.o)$CLINRK = clinrk
  
  return(vcf.o)
}
