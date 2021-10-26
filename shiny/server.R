library(dplyr)
library(shiny)
library(DT)
library(ggplot2)
library(shinydashboard)
library(Rsamtools)
library(GenomicRanges)
library(stringdist)

message('local files')
message(list.files())

message('load genes')
genes.df = read.table('gencode.tsv.gz', as.is=TRUE, header=TRUE)
pli.df = read.table('pli.gene.tsv.gz', as.is=TRUE, header=TRUE)

## list all genes
genes = unique(c(genes.df$gene_id, genes.df$gene_name, genes.df$transcript_id))
## genes = unique(c('ENST00000400454.6', genes))

vars.tbx <- TabixFile('dbvar38.ann.tsv.gz', index='dbvar38.ann.tsv.gz.tbi')

getVars <- function(chr, start, end){
  param <- GRanges(chr, IRanges(start, end))
  res <- scanTabix(vars.tbx, param=param)
  vars.df = read.csv(textConnection(res[[1]]), sep="\t", header=FALSE)
  colnames(vars.df) = c('chr', 'start', 'end', 'variant_id', 'type', 'af', 'clinical_sv', 'clinical_snv')
  vars.df %>% mutate(size=end-start)
}
## vars.df = getVars(1, 0, 1e5)

overlapVarsGenes <- function(vars.df, genes.df){
  vars.gr = makeGRangesFromDataFrame(vars.df)
  genes.gr = makeGRangesFromDataFrame(genes.df)
  var.gene.df = findOverlaps(vars.gr, genes.gr) %>% as.data.frame %>%
    mutate(variant_id=vars.df$variant_id[queryHits],
           gene_id=genes.df$gene_id[subjectHits],
           gene_name=genes.df$gene_name[subjectHits],
           transcript_id=genes.df$transcript_id[subjectHits],
           exon_number=genes.df$exon_number[subjectHits],
           type=genes.df$type[subjectHits]) %>%
    select(-queryHits, -subjectHits) %>% unique
  var.gene.df = var.gene.df %>%
    group_by(variant_id, gene_id, gene_name, transcript_id, type) %>%
    summarize(exon_number=paste(sort(unique(exon_number)), collapse=';'))
  gene.var = var.gene.df %>%
      group_by(variant_id, type) %>%
    summarize(exon_number=ifelse(any(exon_number!=''),
                                 paste0(sort(unique(exon_number)), collapse='|'),
                                 '')) %>%
      mutate(elt_type=ifelse(exon_number!='', paste0(type, '(', exon_number, ')'), type)) %>%
      select(-type)
  merge(gene.var, vars.df) %>%
    group_by(chr, start, end, variant_id, size, type, af, clinical_sv, clinical_snv) %>%
    summarize(gene_impact=paste(sort(unique(elt_type)), collapse=';')) %>% ungroup
}

svtypes = c('DEL', 'DUP', 'INV', 'INS')
svsize.max = 1e9

## function to add links to the variant table
dtify <- function(df){
  if(nrow(df) == 0) return()
  df = df %>% mutate(type=factor(type), coord=paste0(chr, ':', start, '-', end)) %>%
    dplyr::select(-chr, -start, -end) %>%
    dplyr::select(variant_id, coord, type, size, af, everything()) %>%
    mutate(clinical_snv=as.character(clinical_snv),
           variant_id=paste0('<a href="https://www.ncbi.nlm.nih.gov/dbvar/variants/',
                             variant_id, '" target="_blank">', variant_id, '</a>'))
}



## server side of the app
server <- function(input, output, session) {
  ## reactive conductor to extract gene name for a search
  geneName <- reactive({
    message('Gene searched: ', input$gene_search)
    gene_name = input$gene_search
    ## if no genes selected, select first gene name
    if(gene_name == ''){
      gene_name = 'PCSK9'
    }
    if(all(gene_name != genes.df$gene_name)){
      gene.var = genes.df %>%
        filter(gene_id==input$gene_search | gene_name==input$gene_search | transcript_id==input$gene_search)
      gene_name = head(gene.var$gene_name, 1)
    }
    if(length(gene_name) == 0) gene_name = ''
    return(gene_name)
  })
  ## reactive conductor to apply the filtering only once for all elements that need it
  selVars <- reactive({
    genen = geneName()
    message('Gene: ', genen)
    ## find gene
    gene.sel = genes.df %>% filter(gene_id==genen | gene_name==genen | transcript_id==genen)
    if(nrow(gene.sel)==0){
      message('no gene, return')
      return(tibble())
    }
    ## get variants for the gene's region
    vars.df = getVars(gene.sel$chr[1], min(gene.sel$start), max(gene.sel$end))
    message(nrow(vars.df), ' variants')
    ## overlap with genes
    vars.sel = overlapVarsGenes(vars.df, gene.sel)
    ## filter variants
    vars.sel = vars.sel %>% filter(type %in% input$svtypes,
                                   size >= input$size.min,
                                   size <= input$size.max) %>%
      arrange(desc(clinical_sv))
  })
  ## Text
  output$title = renderText({
    genen = geneName()
    if(genen == ''){
      ## agrep(input$gene_search, genes, value=TRUE)
      hints = genes[head(order(stringdist(input$gene_search, genes)))]
      genen = paste('Gene', input$gene_search, 'not found. Close matches:', paste(hints, collapse=' '))
    }
    paste0('<h1>', genen, '</h1>')
  })
  output$omim_url = renderText({
    return(as.character(a('OMIM', href=paste0('https://www.genenames.org/tools/search/#!/genes?query=', geneName()), target='_blank')))
  })
  output$gtex_url = renderText({
    return(as.character(a('GTEx', href=paste0('https://gtexportal.org/home/gene/', geneName()), target='_blank')))
  })
  output$gnomad_url = renderText({
    return(as.character(a('gnomAD', href=paste0('https://gnomad.broadinstitute.org/gene/', geneName()), target='_blank')))
  })
  ## boxes
  output$sv_box <- renderInfoBox({
    infoBox("SVs", nrow(selVars()), icon=icon("dna"), color="blue")
  })
  output$path_sv_box <- renderInfoBox({
    infoBox("Clinical SVs", sum(grepl('Pathogenic', selVars()$clinical_sv)),
            icon=icon("stethoscope"), color="red")
  })
  output$path_snv_box <- renderInfoBox({
    infoBox("Overlap Clinical SNVs", sum(grepl('Pathogenic', selVars()$clinical_snv)),
            icon=icon("stethoscope"), color="red")
  })
  ## dynamic tables
  output$vars_table <- renderDataTable(
    datatable(dtify(selVars()),
    filter='top',
    rownames=FALSE,
    escape=FALSE,
    options=list(pageLength=15, searching=TRUE)))
  ## Graph
  output$af_plot = renderPlot({
    df = selVars()
    if(nrow(df)==0) return()
    ggplot(selVars(), aes(x=af)) + geom_histogram() + theme_bw() + xlab('allele frequency') + xlim(-.1,1.1)
  })

  ## clinical sv ----------------------

  observe({
    shinyjs::toggleState("submit", !is.null(input$file) && input$file != "")
  })

  output$text1 <- renderText({
    paste0(input$file)
  })
  output$text2 <- renderText({
    readSVvcf(paste0(input$file$datapath), out.fmt='vcf', keep.ids=TRUE)
  })


  dataset <- reactive({
    req(input$file$datapath)
    ## read VCF
    suppressWarnings(suppressMessages(library(sveval)))
    suppressWarnings(suppressMessages(library(GenomicRanges)))
    vcf.o = readSVvcf(input$file$datapath, out.fmt='vcf', keep.ids=TRUE)

    ## make sure chromosomes are in the form 'chrX'
    if(all(!grepl('chr', seqlevels(vcf.o)))){
      seqlevels(vcf.o) = paste0('chr', seqlevels(vcf.o))
    }

    ## annotate gene overlapped by SVs
    source('annotate_genes.R')
    vcf.o = annotate_genes(vcf.o, genc)

    ## annotate frequency
    source('annotate_frequency.R')
    vcf.o = annotate_frequency(vcf.o, gnomad)

    ## annotate known clinical SVs
    source('annotate_known_clinical_SVs.R')
    vcf.o = annotate_known_clinical_SVs(vcf.o, clinsv)

    ## clinical ranks, to order the SVs and select top 5 for example
    source('annotate_clinical_score.R')
    vcf.o = annotate_clinical_score(vcf.o,genc)

    ## write annotated VCF
    writeVcf(vcf.o, file=out.vcf)

    ## write tables
    svs <- tibble(gene=info(vcf.o)$GENE,
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
    file.create("input_location.txt")
    fileConn<-base::file(out.dir)
    writeLines(input$file$datapath, "input_location.txt")
    close(fileConn)
    source_python('retrieve.py')
    return(svs)
  })

  output$newvcf <- renderDataTable({
    dataset()
  })
  vcf.o <- readSVvcf("clinical-sv-annotated.vcf", out.fmt='vcf', keep.ids=TRUE)
  # vcf.o <- readVcf("clinical-sv-annotated.vcf")
  # vr.vcf <- makeVRangesFromGRanges(vcf.o)

  observeEvent(input$submit, {
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Annotated")
  })

  output$downloadvcf <- downloadHandler(
    filename = function() {
      paste(input$file, '_annotated', '.vcf', sep='')
    },
    content = function(con) {
      writeVcf(vcf.o, con)
    }
  )

  output$downloadcsv <- downloadHandler(
    filename = function() {
      paste(input$file, '_annotated', '.csv', sep='')
    },
    content = function(con) {
      write.csv(dataset(),con)
    }
  )

  output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$file, '.png', sep='')},
    content = function(file) {
      file.copy("output.png", file, overwrite=TRUE)
    }
  )
  ##---------------------------
}
