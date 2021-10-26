library(shiny)
library(sveval)
library(VariantAnnotation)
library(reticulate)
## clinical sv --------------------------
annot.rdata <- 'annotation_data.RData'
out.vcf <- 'clinical-sv-annotated.vcf'
out.csv <- 'clinical-sv-table.csv'
out.dir <- 'input_location.txt'
load(annot.rdata)
## --------------------------------------
#### launch app locally
runApp()

#### launch app on UCSC server
runApp(port=3457, host='0.0.0.0')
