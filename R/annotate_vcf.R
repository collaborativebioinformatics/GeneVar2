args = commandArgs(TRUE)
if(length(args) == 0){
  ## if no arguments given, download and use a test data
  if(!file.exists('HG002_SVs_Tier1_v0.6.vcf.gz')){
    download.file('ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz', 'HG002_SVs_Tier1_v0.6.vcf.gz')
  }
  in.vcf = 'HG002_SVs_Tier1_v0.6.vcf.gz'
  annot.rdata = 'annotation_data.RData'
  out.vcf = 'clinical-sv-annotated.vcf'
  out.csv = 'clinical-sv-table.csv'
} else {
  in.vcf = args[1]
  annot.rdata = args[2]
  out.vcf = args[3]
  out.csv = args[4]
}

## annotation data contains pre-formatted annotation (e.g. genes, clinical svs)
## see prepare_annotation_data.R
load(annot.rdata)

## read VCF
suppressWarnings(suppressMessages(library(sveval)))
suppressWarnings(suppressMessages(library(GenomicRanges)))
vcf.o = readSVvcf(in.vcf, out.fmt='vcf', keep.ids=TRUE)

## make sure chromosomes are in the form 'chrX'
if(all(!grepl('chr', seqlevels(vcf.o)))){
  seqlevels(vcf.o) = paste0('chr', seqlevels(vcf.o))
}

## if the scripts are not in the current directory, try in a 'R' folder (e.g. when running from the repo's root)
scripts.dir = ''
if(!file.exists('annotate_genes.R')){
  scripts.dir = 'R/'
}

## annotate gene overlapped by SVs
source(paste0(scripts.dir, 'annotate_genes.R'))
vcf.o = annotate_genes(vcf.o, genc)

## annotate frequency
source(paste0(scripts.dir, 'annotate_frequency.R'))
vcf.o = annotate_frequency(vcf.o, gnomad)

## annotate known clinical SVs
source(paste0(scripts.dir, 'annotate_known_clinical_SVs.R'))
vcf.o = annotate_known_clinical_SVs(vcf.o, clinsv)

## clinical ranks, to order the SVs and select top 5 for example
source(paste0(scripts.dir, 'annotate_clinical_score.R'))
vcf.o = annotate_clinical_score(vcf.o)

## write annotated VCF
writeVcf(vcf.o, file=out.vcf)

## write tables
svs = tibble(gene=info(vcf.o)$GENE,
             variant_id=names(vcf.o),
             chr=as.character(seqnames(vcf.o)),
             start=start(vcf.o),
             end=end(vcf.o),
             size=abs(unlist(lapply(info(vcf.o)$SVLEN, '[', 1))),
             frequency=info(vcf.o)$AF,
             svtype=info(vcf.o)$SVTYPE,
             clinsv=info(vcf.o)$CLINSV,
             clinrk=info(vcf.o)$CLINRK) %>%
  arrange(clinrk)

write.table(svs, file=out.csv, sep=',', quote=TRUE, row.names=FALSE)
message("Output written in:", out.csv)
