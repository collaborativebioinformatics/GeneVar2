message('Loading packages...')
suppressWarnings(suppressMessages(library(rtracklayer)))
suppressWarnings(suppressMessages(library(GenomicRanges)))
suppressWarnings(suppressMessages(library(dplyr)))

## arguments: output filename
args = commandArgs(TRUE)
out.rdata = 'annotation_data.RData'
if(length(args) > 0){
  out.rdata = args[1]
}

##
message('Gene annotation from gencode...')
##
if(!file.exists('gencode.v35.annotation.gff3.gz')){
  download.file('http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gff3.gz', 'gencode.v35.annotation.gff3.gz', method='wget')
}
genc.all = import('gencode.v35.annotation.gff3.gz')

## subset to information of interest
genc = subset(genc.all, type %in% c('gene', 'CDS', 'five_prime_UTR', 'three_prime_UTR'))
mcols(genc) = mcols(genc)[,c('type', 'gene_name', 'gene_type', 'gene_id')]
genc = genc %>% as.data.frame %>% 
  select(seqnames, start, end, strand, type, gene_name, gene_id, gene_type) %>%
  mutate(type=as.character(type)) %>%
  unique %>% 
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)

##
message('gnomAD-SV SV catalog with frequencies...')
##
## gnomad_v2.1_sv.sites.lifted.tsv.gz has been saved in our dnanexus project space
gnomad = read.table('gnomad_v2.1_sv.sites.lifted.tsv.gz', as.is=TRUE, header=TRUE)
## match frequency to variant
matchAF <- function(df){
  if(nrow(df) == 1){
    ## if only one allele, no need to split the AF field
    return(df)
  }
  afs = unlist(strsplit(df$AF[1], split=','))
  df$AF = afs
  return(df)
}
gnomad = gnomad %>% group_by(seqnames, start, end, type, size, ref) %>%
  do(matchAF(.)) %>% mutate(AF=as.numeric(AF))
## convert CNV to DEL/DUP
gnomad = gnomad %>% mutate(type=ifelse(grepl('<CN=', alt), 'DUP', type),
                           type=ifelse(alt %in% c('<CN=0>', '<CN=1>'), 'DEL', type)) %>%
  filter(alt!='<CN=2>')
gnomad = gnomad %>% select(-qual) %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE)

##
message('Known clinical SVs...')
##
if(!file.exists('nstd102.GRCh38.variant_call.tsv.gz')){
  download.file('http://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/tsv/nstd102.GRCh38.variant_call.tsv.gz', 'nstd102.GRCh38.variant_call.tsv.gz', method='wget')
}
clinsv = read.table('nstd102.GRCh38.variant_call.tsv.gz', as.is=TRUE, skip=1, quote='', comment='', sep='\t', header=TRUE)
clinsv = clinsv[, c(1,3,8,10:16,37)] %>%
  rename(dbvar_id=X.variant_call_accession) %>%
  mutate(chr=paste0('chr', chr), 
         type=ifelse(grepl('gain', variant_call_type), 'DUP', NA),
         type=ifelse(grepl('loss', variant_call_type), 'DEL', type),
         type=ifelse(grepl('duplication', variant_call_type), 'DUP', type),
         type=ifelse(variant_call_type=='deletion', 'DEL', type),
         type=ifelse(variant_call_type=='insertion', 'INS', type),
         type=ifelse(variant_call_type=='inversion', 'INV', type),
         inner_start=ifelse(is.na(inner_start), outer_start, inner_start),
         inner_stop=ifelse(is.na(inner_stop), outer_stop, inner_stop),
         outer_start=ifelse(is.na(outer_start), inner_start, outer_start),
         outer_stop=ifelse(is.na(outer_stop), inner_stop, outer_stop),
         start=ifelse(is.na(start), (inner_start+outer_start)/2, start),
         end=ifelse(is.na(stop), (inner_stop+outer_stop)/2, stop),
         size=ifelse(!is.na(insertion_length), insertion_length, end-start)) %>% 
  select(chr, start, end, type, size, dbvar_id, clinical_significance) %>% 
  filter(!is.na(type), !is.na(start), size>=50) %>%
  makeGRangesFromDataFrame(keep.extra.columns=TRUE)

##
message('Gene with high probability of loss-of-function intolerance...')
##
if(!file.exists('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')){
  download.file('https://azureopendatastorage.blob.core.windows.net/gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', 'gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', method='wget')
}
pli = read.table('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', header=TRUE, as.is=TRUE, sep='\t')
gene.pli = pli %>% filter(pLI>.9) %>% .$gene %>% unique

##
message('Gene associated with phenotype in OMIM...')
##
if(!file.exists('gene_list_terms.txt.gz')){
  download.file('https://maayanlab.cloud/static/hdfs/harmonizome/data/omim/gene_list_terms.txt.gz', 'gene_list_terms.txt.gz', method='wget')
}
omim = read.table('gene_list_terms.txt.gz', header=TRUE, as.is=TRUE, sep='\t')
gene.omim = unique(omim$GeneSym)

##
## save in one file
##
message('Save annotation to ', out.rdata, '...')
save(genc, gnomad, clinsv, gene.pli, gene.omim, file=out.rdata)


