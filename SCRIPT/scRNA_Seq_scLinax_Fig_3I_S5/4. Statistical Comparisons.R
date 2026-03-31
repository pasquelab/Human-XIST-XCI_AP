library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(scLinaX)
library(reshape2)

#this analysis is the same as done in script 3, but adds extra bits at the end requested by the reviewers.
#mainly this adds statistical comparions between the allelic ratios of different days and cell types, 
#and limits it so that each comparison only includes genes that are informative in both days/cell types
  
  #defining colors that I will use for my graphs
  XCIcolors <-scale_color_manual(values=c("escapes from XCI"="#9AB1DE", 
                                          "no XCI status known"="grey", 
                                          "PAR"="dark blue",
                                          "subject to XCI"="#FF9777", 
                                          "variably escapes XCI"="#FFE543",
                                          "autosomal" = "black",
                                          "escape"="#4863a3",
                                          "inactive"="#9e4e36",
                                          "variable"="#b59d0e",
                                          "Unknown"="dark grey"))
  
  XCIfill <-scale_fill_manual(values=c("escapes from XCI"="#8AA1CE", 
                                       "no XCI status known"="light grey", 
                                       "PAR"="dark blue",
                                       "subject to XCI"="#FC8767", 
                                       "variably escapes XCI"="#FFD533",
                                       "autosomal" = "dark grey",
                                       "escape"="#8AA1CE",
                                       "inactive"="#FC8767",
                                       "variable"="#FFD533",
                                       "Unknown"="light grey"))
  

  
  
  data("XCI_ref")
  XCI_ref$XCI_status[which(XCI_ref$Gene=="XIST")] <- "inactive"

    #assuming that you only have one file formatted for sclinax in your folder. Can manually set this to your filename.
    inputs <- read.table(list.files("./", pattern="*.formattedForscLinax.tsv")[1], h=T, sep="\t", quote="")
  
  
  
  #adding a filter to remove SNPs that are rarely seen in gnomad and might not be real SNPs
  QCfull<-read.table(list.files("./", pattern="hg38_multianno.txt")[1], h=T, sep="\t")
  inputs <- inputs[which(inputs$POS%in%QCfull$Start[which(QCfull$gnomad41_genome_AF>0.01)]),]
  
  inputs$Sample_ID <- "SigmaWT"

  #make sure that the cell barcodes in the annotation file match the cell barcodes in the input file
  EXMCannotation <- read.table("timecourse_annotation.tsv", h=T)
  EXMCannotation$cell_barcode <- gsub(EXMCannotation$cell_barcode, pattern = "[.]", replacement = "-")
  EXMCannotationJustDay <- EXMCannotation
  EXMCannotationJustDay$Annotation <- unlist(lapply(EXMCannotation$Annotation, function(x){return(strsplit(x, split="_")[[1]][1])}))  

  
  
  QCref <- run_RefGeneQC(inputs[which(inputs$cell_barcode%in%EXMCannotationJustDay$cell_barcode[which(EXMCannotationJustDay$Annotation%in%c("D8", "D13", "D18", "D70"))]),], 
                         XCI_ref, SAMPLE_NUM_THR=1)

  
  
  sclinaxResult <- run_scLinaX(ASE_df=inputs,XCI_ref=XCI_ref,QCREF = QCref,
                               Inactive_Gene_ratio_THR=0.15,SNP_DETECTION_DP=10,SNP_DETECTION_MAF=0.05,QC_total_allele_THR=10,
                               HE_allele_cell_number_THR=10,REMOVE_ESCAPE=TRUE,PVAL_THR=0.01,RHO_THR=0.5)

  
  #for the supplemental version the cells are split just by day
  combinedScLinaXSummarySplitByDay <- summarize_scLinaX(sclinaxResult, QC_total_allele_THR=10, Annotation=EXMCannotationJustDay)
  combinedScLinaXSummarySplitByDay$Day <- factor(combinedScLinaXSummarySplitByDay$Annotation, levels=c("D0", "D1", "D2", "D4", "D8", "D13", "D18", "D70"))
  combinedScLinaXSummarySplitByDay$XCI_status[which(combinedScLinaXSummarySplitByDay$XCI_status=="escape")] <- "Escapes from XCI"
  combinedScLinaXSummarySplitByDay$XCI_status[which(combinedScLinaXSummarySplitByDay$XCI_status=="inactive")] <- "Subject to XCI"  


  #For the main figure version the cells are split by cell type and day
  #I am combinging some of the cell annotations from Pham, Panda, et al. to make interpreting the output easier
  EXMCannotationCombineCellTypes <- EXMCannotation
  EXMCannotationCombineCellTypes$Annotation <- as.character(EXMCannotationCombineCellTypes$Annotation)
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D0_8CLC", "D0_naive"))]<- "D0_naive"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D0_TSC"))]<- "D0_TSC"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D1_8CLC", "D1_naive", "D1_intermediate.1","D1_intermediate.3"))]<- "D1_Epiblast"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D1_TSC"))]<- "D1_TSC"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D13_early.EXMC"))]<- "D13_EXMC"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D13_8CLC", "D13_naive", "D13_intermediate.1","D13_intermediate.2","D13_intermediate.3"))]<- "D13_Epiblast"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D13_TSC"))]<- "D13_TSC"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D18_early.EXMC"))]<- "D18_EXMC"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D18_TSC"))]<- "D18_TSC"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D2_8CLC", "D2_naive", "D2_intermediate.1","D2_intermediate.2","D2_intermediate.3"))]<- "D2_Epiblast"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D2_TSC"))]<- "D2_TSC"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D4_8CLC", "D4_naive", "D4_intermediate.1","D4_intermediate.2","D4_intermediate.3"))]<- "D4_Epiblast"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D4_TSC"))]<- "D4_TSC"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D70_late.EXMC"))]<- "D70_EXMC"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D8_early.EXMC"))]<- "D8_EXMC"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D8_8CLC", "D8_naive", "D8_intermediate.1","D8_intermediate.2","D8_intermediate.3"))]<- "D8_Epiblast"
  EXMCannotationCombineCellTypes$Annotation[which(EXMCannotation$Annotation%in%c("D8_TSC"))]<- "D8_TSC"  

  
  
  
  
  combinedScLinaXSummary <- summarize_scLinaX(sclinaxResult, QC_total_allele_THR=10, Annotation=EXMCannotationCombineCellTypes)
  combinedScLinaXSummary$Day <- unlist(lapply(combinedScLinaXSummary$Annotation, function(x){return(strsplit(x, split="_")[[1]][1])}))
  combinedScLinaXSummary$CellType <- unlist(lapply(combinedScLinaXSummary$Annotation, function(x){return(strsplit(x, split="_")[[1]][2])}))
  combinedScLinaXSummary$XCI_status[which(combinedScLinaXSummary$XCI_status=="escape")] <- "Escapes from XCI"
  combinedScLinaXSummary$XCI_status[which(combinedScLinaXSummary$XCI_status=="inactive")] <- "Subject to XCI"
  combinedScLinaXSummary$Day <- factor(combinedScLinaXSummary$Day, levels=c("D0", "D1", "D2", "D4", "D8", "D13", "D18", "D70"))

  
  
  
  ##stats table doing t.tests to compare using only genes informative in both cell types
  
  statsTable <- matrix(nrow=length(unique(combinedScLinaXSummary$Annotation)), ncol=length(unique(combinedScLinaXSummary$Annotation)))
  colnames(statsTable) <- unique(combinedScLinaXSummary$Annotation)
  rownames(statsTable) <- unique(unique(combinedScLinaXSummary$Annotation))
  
  statsTableSharedInformativeGenes <- statsTable

  statsTableMeta <- data.frame(anno=colnames(statsTable), day=NA, cellType=NA)
  statsTableMeta$day <- unlist(lapply(statsTableMeta$anno, function(x){return(strsplit(x, split="_")[[1]][1])}))
  statsTableMeta$cellType <- unlist(lapply(statsTableMeta$anno, function(x){return(strsplit(x, split="_")[[1]][2])}))
  
  for(i in 1:nrow(statsTableMeta)){
    for(j in 1:nrow(statsTableMeta)){
      informativeGenes <- combinedScLinaXSummary$Gene[which(combinedScLinaXSummary$Gene%in%combinedScLinaXSummary$Gene[which(combinedScLinaXSummary$Annotation==statsTableMeta$anno[i])] &
                                                              combinedScLinaXSummary$Annotation==statsTableMeta$anno[j]  &
                                                              combinedScLinaXSummary$XCI_status=="Subject to XCI")]
      statsTableSharedInformativeGenes[i,j] <- length(informativeGenes)
      
      if(length(which(combinedScLinaXSummary$Gene%in%informativeGenes &
                      combinedScLinaXSummary$Annotation==statsTableMeta$anno[i]))>1 &
         length(which(combinedScLinaXSummary$Gene%in%informativeGenes &
                      combinedScLinaXSummary$Annotation==statsTableMeta$anno[j]))>1 )
      statsTable[i,j] <- t.test(combinedScLinaXSummary$minor_allele_ratio[which(combinedScLinaXSummary$Gene%in%informativeGenes &
                                                                                  combinedScLinaXSummary$Annotation==statsTableMeta$anno[i])],
                                combinedScLinaXSummary$minor_allele_ratio[which(combinedScLinaXSummary$Gene%in%informativeGenes &
                                                                                  combinedScLinaXSummary$Annotation==statsTableMeta$anno[j])])$p.value
            
    }
    
  }
  
  
  statsTable
write.table(statsTable, file="t.testComparisons.scLinaX.tsv", sep="\t", col.names = T, row.names=T, quote=F)  
write.table(statsTableSharedInformativeGenes, file="t.testComparisons.sharedInformativeGenes.scLinaX.tsv", sep="\t", col.names = T, row.names=T, quote=F)  


statsTableMelt <- data.frame(melt(statsTable), ninformativeGenes =melt(statsTableSharedInformativeGenes)$value)
colnames(statsTableMelt)[3] <- "p-values"

head(statsTableMelt)
write.table(statsTableMelt, file="t.testComparisons.melt.scLinaX.tsv", sep="\t", col.names = T, row.names=F, quote=F)  

write.table(combinedScLinaXSummary, file="scLinaXResults.tsv", sep="\t", col.names=T, row.names=F, quote=T)

combinedScLinaXSummary <- read.table("scLinaXResults.tsv", h=T)

ggplot(combinedScLinaXSummary, aes(x=CellType, y=Reference_Cell_Count))+geom_violin(scale="width")+
  theme_bw()+ylab("Number of cells used for scLinaX")+xlab("Cell Type")
