##########################
##    Hackathon 2021    ##
##    2021 Oct 10-13    ##
##      GeneVar2        ##
##      SV vcf/csv      ## 
##      Plots in R      ##
##########################

#load libraries
lapply(c("vcfR","ggplot2","scales","ggpubr"),require, character.only = TRUE)
chr_names <-  c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                "chr12","chr13","chr14","chr15","chr16","chr17","chr18",
                "chr19","chr20","chr21","chr22","chrX","chrY")

#call in VCF
in_filename <- "NA19461.clinicalsv.csv";in_filename <- "NA19461.clinicalsv.vcf"
dir_data <- "~/"
if(length(grep('vcf',in_filename))){
  in_data <- read.vcfR(paste0(dir_data,in_filename))
  data_fix <- data.frame(getFIX(in_data))
  try(data_fix$svtype <- gsub('Manta','',do.call(rbind,strsplit(data_fix$ID,":"))[,1]),silent = T)
  data_fix <- data_fix[,c("CHROM","POS","svtype")]
  meta_info <- in_data@meta[grep('contig',in_data@meta)]
  meta_chr <- unlist(lapply(chr_names,function(x){
    meta_info[grep(paste0("##contig=<ID=",x,","),meta_info)]
  }))
  sum_GRChX <- data.frame(CHROM=chr_names,
                          LEN=as.numeric(gsub('>','',gsub('length=','',do.call(rbind,strsplit(meta_chr,','))[,2]))))
}else if(length(grep('csv',in_filename))){
  in_data <- read.csv(paste0(dir_data,in_filename))
  data_fix <- in_data[,c("chr","start","svtype")]
  colnames(data_fix) <- c("CHROM","POS","svtype")
  # GRCh38 data pulled from: https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/#/st
  sum_GRChX <- data.frame(CHROM=chr_names,
                          LEN=c(249698942,242508799,198450956,190424264,181630948,170805979,159345973,145138636,138688728,133797422,135186938,133275309,114364328,108136338,102439437,92211104,83836422,80373285,58617616,64444167,46709983,51857516,156040895,57264655))
}

sum_GRChX$ToTLEN <- unlist(lapply(1:nrow(sum_GRChX),function(x) sum(sum_GRChX$LEN[1:x])))
sum_GRChX$StartPos <- c(1,1+sum_GRChX$ToTLEN[-nrow(sum_GRChX)])

#manipulate data for graph
data_fix$CHROM <- factor(data_fix$CHROM, ordered = TRUE, 
                         levels = chr_names)
# Plot: Allele Freq ####
allel_freq <- data.frame(table(data_fix$CHROM,data_fix$svtype));colnames(allel_freq) <- c("CHROM","InDel","Freq")
allel_freq$InDel <- toupper(allel_freq$InDel)
alFr_title <- NULL

plt_allelFreq <- ggplot(allel_freq,aes(x=InDel,y=Freq ,fill=InDel))+geom_bar(stat = 'identity')+facet_wrap(.~CHROM,nrow = 1)+
  labs(x=NULL,y="Count",title=alFr_title)+
  scale_fill_manual(values=c("red","darkgreen",'blue',"purple"),breaks=c("DEL","DUP","INS","BND"),labels=c("DEL - Deletion","DUP -Duplication","INS - Insertion","BND - Breakends"))+
  guides(fill=guide_legend(title="VCF Call"))+
  scale_y_continuous(expand = c(0, 0),labels = label_number(accuracy = 1),breaks=pretty_breaks()) +
  theme_bw()+theme(strip.background =element_rect(fill="white"),
                   # axis.line = element_line(colour = "black"),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

# plot:: SVType FREQ ####
plt_svtype_fun <- function(SVType="DEL",include_title=T,bin_wid=1000000){
  if(length(SVType)==1){
    if(toupper(SVType) %in% c("DEL","DUP","INS","BND")){
      in_del <- data_fix[which(data_fix$svtype==SVType),]
      in_del$genoPOS <- apply(in_del,1,function(x)sum_GRChX$StartPos[which(sum_GRChX$CHROM==as.character(x[1]))]+as.numeric(x[2]))
      # options(scipen = 100)
      alFr_title <- if(include_title) paste("Structural Variant:",SVType) else NULL
      plt_del <- ggplot(in_del,aes(x=genoPOS))+annotate('rect',xmin=sum_GRChX$StartPos[seq(2,nrow(sum_GRChX),2)],
                                                        xmax=sum_GRChX$ToTLEN[seq(2,nrow(sum_GRChX),2)],ymin=0,ymax=2,alpha=0.2)+
        geom_histogram(binwidth = bin_wid,aes(fill = CHROM))+
        scale_x_continuous(limits = c(0,max(sum_GRChX$ToTLEN)),expand = c(0, 0),name=NULL, breaks=sum_GRChX$StartPos, labels=sum_GRChX$CHROM)+
        labs(y="Counts",title=alFr_title)+
        scale_fill_manual(values=rep(c("red","black"),length(sum_GRChX$CHROM)/2),breaks=sum_GRChX$CHROM,labels=NULL,guide='none')+
        scale_y_continuous(expand = c(0, 0),breaks= function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
        theme_bw()+theme(
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 
      plt_del
    }}else{
      if(all(toupper(SVType) %in% c("DEL","DUP","INS","BND"))){
        in_del <- data_fix[data_fix$svtype%in% SVType,]
        in_del$genoPOS <- apply(in_del,1,function(x)sum_GRChX$StartPos[which(sum_GRChX$CHROM==as.character(x[1]))]+as.numeric(x[2]))
        # options(scipen = 100)
        alFr_title <- paste("Structural Variant:",SVType);names(alFr_title) <- SVType
        plt_del <- ggplot(in_del,aes(x=genoPOS))+annotate('rect',xmin=sum_GRChX$StartPos[seq(2,nrow(sum_GRChX),2)],
                                                          xmax=sum_GRChX$ToTLEN[seq(2,nrow(sum_GRChX),2)],ymin=0,ymax=2,alpha=0.2)+
          geom_histogram(binwidth = bin_wid,aes(fill = CHROM))+
          
          facet_wrap(svtype~.,ncol = 1,labeller = labeller(svtype=alFr_title))+
          scale_x_continuous(limits = c(0,max(sum_GRChX$ToTLEN)),expand = c(0, 0),name=NULL, breaks=sum_GRChX$StartPos, labels=sum_GRChX$CHROM)+
          labs(y="Counts",title=NULL)+
          scale_fill_manual(values=rep(c("red","black"),length(sum_GRChX$CHROM)/2),breaks=sum_GRChX$CHROM,labels=NULL,guide='none')+
          scale_y_continuous(expand = c(0, 0),breaks= function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))) +
          theme_bw()+theme(
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) 
        plt_del
      }
    }
}
plt_del <- plt_svtype_fun("DEL")
plt_dup <- plt_svtype_fun("DUP")
plt_ins <- plt_svtype_fun("INS")
plt_bnd <- plt_svtype_fun("BND")
plt_all <- plt_svtype_fun(c("DEL","DUP","INS","BND"),bin_wid=1000000)

# save plots ####
ggsave(plt_allelFreq,file='./foo_freq.jpeg')

# Individual Structural Variants 
ggsave("foo.pdf",plt_del,width = 8,height = 3)
# library("gridExtra")
# ggsave("foo4.pdf", arrangeGrob(plt_del,plt_dup,plt_ins,plt_bnd,ncol=1),width = 8,height = 3)
# Freq of Structural Variants
# all Structural Variants 
ggsave(plt_all,file='./foo_grid.jpeg',width = 9,height = 3)

