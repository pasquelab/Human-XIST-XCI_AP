
setwd('C:/Users/gael/Desktop/test_jeanne/bulk_RNAseq_Leo/allelic_analysis/R_analysis_like_bradley/')

library(ggplot2)
library(dplyr)
library(tibble)
library(gridExtra)
library(viridis)

outFilePrefix <- "C:/Users/gael/Desktop/test_jeanne/bulk_RNAseq_Leo/allelic_analysis/R_analysis_like_bradley/output_figs"



# # Uploading depthperSNP files
# setwd("X:/jeanne_test/bulk_RNAseq_Leo/allelic_analysis_bradley_method/bam/SNPsplit/depthPerSNP/")
# refFiles <- list.files(pattern = "*.genome1.depthPERSNP")
# altFiles <- list.files(pattern = "*.genome2.depthPERSNP")
# 
# print("checking length of ref and alt files")
# length(refFiles)==length(altFiles)
# 
# combinedFiles <- vector(mode="list")
# for(i in 1:length(refFiles)){
#   combinedFiles[[i]] <- vector(mode="list")
#   combinedFiles[[i]]$ref <- unique(read.table(refFiles[i], h=F)[,1:5])
#   colnames(combinedFiles[[i]]$ref) <- c("chr", "start", "stop", "length", "refDepth")
#   combinedFiles[[i]]$alt <- unique(read.table(altFiles[i], h=F)[,1:5])
#   colnames(combinedFiles[[i]]$alt) <- c("chr", "start", "stop", "length", "altDepth")
# }
# 
# 
# 
# # combining them into one file
# allReadTotal <- matrix(nrow=nrow(combinedFiles[[1]]$ref), ncol= 2*length(combinedFiles))
# colnames(allReadTotal) <- c(refFiles, altFiles)
# 
# for(i in 1:length(combinedFiles)){
#   allReadTotal[,i]<- combinedFiles[[i]]$ref$refDepth
#   allReadTotal[,i+length(combinedFiles)] <- combinedFiles[[i]]$alt$altDepth
# }
# 
# allReadTotalStack <- data.frame(chr = combinedFiles[[i]]$ref$chr, SNP = combinedFiles[[1]]$ref$start, stack(as.data.frame(allReadTotal[,1:length(combinedFiles)])), 
#                                 stack(as.data.frame(allReadTotal[,(1+length(combinedFiles)):(2*length(combinedFiles))])))
# 
# levels(allReadTotalStack$ind) <- c('KD dox -', 'KD dox +', 'WT dox -', 'WT dox +')
# 
# 
# 
# 
# 
# ##### associating snp to gene data (name, status) #####
# setwd('C:/Users/gael/Desktop/test_jeanne/bulk_RNAseq_Leo/allelic_analysis/R_analysis_like_bradley/')
# 
XCIstatuses <- read.table("newcalls.xci.txt", h=T)# this is from my 2015 paper (Overallscore), newscore is newer calls I made but never published. just use the published ones
## exons <- read.table("hg38.exonsWgene.tsv", h=F)#this is chrX only
allExons <- read.table("hg38.allChrs.exonsWgenes.tsv", h=T)
# 
# #could rewrite getStatus to have gene as input instead of Xpos
# getGene <- function(Xpos){
#   #Xpos is a single entry from V2 from my vcf, but filtered to only be chrX as that is all I have the exons for
#   if(any(exons$V2<Xpos & exons$V3>Xpos)){
#     return(paste(unique(exons$V4[which(exons$V2<Xpos & exons$V3>Xpos)]), collapse=","))
#   } else return("no gene found")
# }
# 
# getStatus <- function(Xpos){
#   #Xpos is a single entry from V2 from my vcf, but filtered to only be chrX as that is all I have the exons for
#   
#   gene <- exons$V4[which(exons$V2<Xpos & exons$V3>Xpos)]
#   status <- XCIstatuses$Overallscore[which(XCIstatuses$hg19.kgXref.geneSymbol%in%gene)]
#   
#   if(any(c("Discordant", "VE", "MostlyVE")%in%status) | 
#      (any(c("E","MostlyE")%in%status) & any(c( "S", "MostlyS")%in%status))){
#     return("variably escapes XCI")
#   } else if(any(c("S", "MostlyS")%in%status)) {
#     return("subject to XCI")
#   } else if(any(c("E", "MostlyE")%in%status)){
#     return("escapes from XCI")
#   } else if("PAR"%in%status){
#     return("PAR")
#   } else{ return("no XCI status known")}
# }
# 
#
# allReadTotalStack$gene <- NA
# allReadTotalStack$status <- NA
# 
# allReadTotalStack$gene[which(allReadTotalStack$chr=="chrX")] <- unlist(lapply(allReadTotalStack$SNP[which(allReadTotalStack$chr=="chrX")], getGene))
# allReadTotalStack$status[which(allReadTotalStack$chr=="chrX")] <- unlist(lapply(allReadTotalStack$SNP[which(allReadTotalStack$chr=="chrX")], getStatus))
# 
# # Expand to autosomes
# getGeneIncludeAutosomes <- function(SNPpos, chrExons){
# 
#   if(any(chrExons$start<SNPpos & chrExons$stop>SNPpos)){
#     return(paste(unique(chrExons$gene[which(chrExons$start<SNPpos & chrExons$stop>SNPpos)]), collapse=","))
#   } else return("no gene found")
# 
# }
# # takes forever
# for(chr in unique(allReadTotalStack$chr)){
#   print(chr)
#   allReadTotalStack$gene[which(allReadTotalStack$chr==chr)] <- unlist(lapply(allReadTotalStack[which(allReadTotalStack$chr==chr),2], getGeneIncludeAutosomes, chrExons = allExons[which(allExons$chr==chr),]))
#   print(paste0(chr, " done. Next chr"))
# }
# 
# # Jeanne:
# allReadTotalStack$status[which(allReadTotalStack$chr!="chrX" & allReadTotalStack$gene != "no gene found")] <- "autosomal" 
# 
# 
# # setwd('C:/Users/gael/Desktop/test_jeanne/bulk_RNAseq_Leo/allelic_analysis/R_analysis_like_bradley/')
# # write.table(allReadTotalStack, "allReadTotalStack_bulk-Leo.csv", sep = ",", col.names = T, row.names = F, quote = F)
# 
# 
# 
# 
# ##### Ratios #####
# 
# getMinMaxRatio<- function(reads){ # adapting this from my previous papers XCI calling
#   
#   if(length(which(!is.na(reads)))==3){
#     minX <- median(reads,  na.rm=T) #using median instead of min because I have 3 alleles and don't want the 3rd one.
#     maxX <- max(reads,  na.rm=T)
#   } else if(length(which(!is.na(reads)))==2){
#     minX <- min(reads,  na.rm=T) 
#     maxX <- max(reads,  na.rm=T)
#   } else{print("Error in getMinMaxRatio. Reads input should have 2 or 3 non-NA integers")}
#   ratio=minX/maxX
#   
#   return(data.frame(
#     ratio=ratio,
#     minR=ratio-1.96*sqrt((1-ratio)*ratio)/(minX+maxX),
#     maxR=ratio+1.96*sqrt((1-ratio)*ratio)/(minX+maxX)
#   ))
# }
# 
# 
# makeXCIstatusCalls <- function(minMax){#input is the min and max (list of 2 ?)
#   
#   if(is.na(minMax[2])){
#     return(NA)
#   }
#   
#   if(minMax[2]<0.10){
#     return("subject to XCI")
#   } else if(minMax[1]>0.1){
#     return("escapes from XCI")
#   } else if(minMax[1]<0.1 & minMax[2]>0.1){
#     return("threshold")
#   } else return(NA)
#   
# }
# 
# head(allReadTotalStack[which(allReadTotalStack$chr=="chrX"),])
# 
# allReadTotalStack$XminXmaxRatio <- NA
# allReadTotalStack$minError <- NA
# allReadTotalStack$maxError <- NA
# 
# # chrX only 
# # allReadTotalStack[which(allReadTotalStack$chr=="chrX"),9:11] <-matrix(unlist(apply(allReadTotalStack[which(allReadTotalStack$chr=="chrX"),c(3,5)], 1, getMinMaxRatio)), ncol = 3, byrow = T)
# # colnames(allReadTotalStack)[9:11] <- c("XminXmaxRatio", "minError","maxError")
# 
# 
# # chrX + autosomes (Jeanne)
# #allReadTotalStack[,9:11] <-matrix(unlist(apply(allReadTotalStack[,c(3,5)], 1, getMinMaxRatio)), ncol = 3, byrow = T) # too long
# allReadTotalStack[which(allReadTotalStack$gene != "no gene found"),9:11] <- matrix(unlist(apply(allReadTotalStack[which(allReadTotalStack$gene != "no gene found"),c(3,5)], 1, getMinMaxRatio)), ncol = 3, byrow = T)
# colnames(allReadTotalStack)[9:11] <- c("XminXmaxRatio", "minError","maxError")
# 
# allReadTotalStack$hereStatus <- NA
# allReadTotalStack$hereStatus[which(allReadTotalStack$chr=="chrX" & allReadTotalStack$gene != "no gene found")] <- apply(allReadTotalStack[which(allReadTotalStack$chr=="chrX" & allReadTotalStack$gene != "no gene found"),10:11], 1, makeXCIstatusCalls)
# 
# # Filter read depth >= 30
# allReadTotalStack$depthFilter <- apply(allReadTotalStack[,c(3,5)], 1, sum)>29
# # add a column with depth info
# allReadTotalStack$depth <- apply(allReadTotalStack[,c(3,5)], 1, sum)
# 
# # Saving the table (because very long to generate again): 
# setwd('C:/Users/gael/Desktop/test_jeanne/bulk_RNAseq_Leo/allelic_analysis/R_analysis_like_bradley/')
# # write.table(allReadTotalStack, "allReadTotalStack_bulk-Leo.csv", sep = ",", col.names = T, row.names = F, quote = F)
# allReadTotalStack <- read.csv("allReadTotalStack_bulk-Leo.csv", header = T)
# 
# 
# ##### Per gene #####
# #allReadTotalStack <- allReadTotalStack[which(allReadTotalStack$chr=="chrX" & allReadTotalStack$ind!="Undetermined"),] # just chrX
# 
# allelicRatiosPerGene <- data.frame(unique(allReadTotalStack[order(allReadTotalStack$SNP),7:8]))[-1,] # remove the first line corresponding to 'no gene found'
# allelicRatiosPerGene <- cbind(allelicRatiosPerGene, matrix(nrow = nrow(allelicRatiosPerGene), ncol=4))
# colnames(allelicRatiosPerGene)[3:6] <- levels(allReadTotalStack$ind)[1:4]
# 
# getMedianperGene <- function(geneName, sampleName, minSNP=1, minDepth=29){
#   if(length(which(allReadTotalStack$gene==geneName & allReadTotalStack$ind==sampleName & allReadTotalStack$depth>minDepth))>minSNP-1){
#     return(median(allReadTotalStack$XminXmaxRatio[which(allReadTotalStack$gene==geneName & allReadTotalStack$ind==sampleName& allReadTotalStack$depth>minDepth)], na.rm=T))   
#   } else return(NA)
# }
# 
# #for(r in which(allelicRatiosPerGene$gene%in%exons$V4[which(exons$V1=="chrX")])){###use this line for only X chromosome
# for(r in 1:nrow(allelicRatiosPerGene)){ ### use this line to process all autosomes
#   for(c in 3:6){
#     allelicRatiosPerGene[r,c] <- getMedianperGene(allelicRatiosPerGene$gene[r], colnames(allelicRatiosPerGene)[c], minDepth=29)
#     
#   }
# }
# 
# summary(allelicRatiosPerGene) #should be between 0 and 1
# 
# # adding column with chr
# allelicRatiosPerGene$chr <- NA
# for(chr in unique(allReadTotalStack$chr)){
#   allelicRatiosPerGene$chr[which(allelicRatiosPerGene$gene%in%allReadTotalStack$gene[which(allReadTotalStack$chr==chr)])] <- chr   
# }



# saving table because very long to generate: 
# write.table(allelicRatiosPerGene, "allelicRatiosPerGene.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
# opening the saved table: 
allelicRatiosPerGene <- read.table("allelicRatiosPerGene.tsv", header = T, sep = "\t")
colnames(allelicRatiosPerGene)[3:6] <- unlist(c('KD dox -', 'KD dox +', 'WT dox -', 'WT dox +'))

# from wide to long df
allelicRatiosPerGeneStack <- data.frame(allelicRatiosPerGene[,1:2], chr = allelicRatiosPerGene$chr, stack(allelicRatiosPerGene[,3:6]))

# orders for plot via factors
allelicRatiosPerGeneStack$statusOrdered <- factor(allelicRatiosPerGeneStack$status, levels=c("autosomal", "escapes from XCI", "subject to XCI", "no XCI status known", "variably escapes XCI"))
allelicRatiosPerGeneStack$cellTypeOrdered <- factor(allelicRatiosPerGeneStack$ind, levels=c('WT dox -', 'WT dox +','KD dox -', 'KD dox +'))


##### Plot Bradley #####
# plot 1 : median of ratio per gene 
#options(repr.plot.width=20)
#pdf(file=paste(c(outFilePrefix, "Xmin/XmaxRatio.medianSNPperGene.wLegend.pdf"), collapse=""), height=4, width=10)
ggplot(allelicRatiosPerGeneStack[which(allelicRatiosPerGeneStack$status%in%c("autosomal", "escapes from XCI", "subject to XCI", "no XCI status known")),], aes(x=statusOrdered, y=values, fill=cellTypeOrdered))+
  geom_boxplot(size=0.5, outlier.size = 0.5) +
  theme_linedraw() +
  ylab("Median SNP Xmin/Xmax ratio per gene") +
  xlab("Published XCI status") + guides(fill=guide_legend(title = "Sample")) +
  scale_fill_manual(values=c('WT dox -'='#8CCEBA', 'WT dox +'='#EBAE80','KD dox -'="#B9B7D8", 'KD dox +'="#F293C4")) +
  theme(panel.grid= element_blank(), panel.grid.major.y=element_blank(), legend.position="bottom")
#dev.off()



# plot 2: violin plot mean allelic ratio per gene 
## CAUTION: JUST CHRX WANTED HERE
#options(repr.plot.width=5)
#pdf(file=paste(c(outFilePrefix, "violinplot.minmaxRatio.medianSNPperGene.allgenes.pdf"), collapse=""), height=5, width=4)
ggplot(allelicRatiosPerGeneStack[allelicRatiosPerGeneStack$chr=="chrX",],
       aes(x=cellTypeOrdered, y=values*100, fill=cellTypeOrdered)) +
  geom_violin(scale="width", alpha=0.5) +
  theme_linedraw() +
  ylab("Average allelic ratio (min/max) per gene") +
  geom_boxplot(position=position_dodge(0.9), width=0.2, alpha=1) +
  xlab("") +
  guides(fill=guide_legend(title = "Gene allelic expression status in primed cells")) +
  scale_fill_manual(values=c('WT dox -'='#8CCEBA', 'WT dox +'='#EBAE80','KD dox -'="#B9B7D8", 'KD dox +'="#F293C4")) +
  geom_hline(yintercept=c(25), linetype = 2, color = 'red') +
  theme(text = element_text(size=16), legend.position="none", panel.grid.major = element_blank(), panel.grid.minor=element_blank())
#dev.off()


##### Plot Jeanne #####
# plot 3: simple bar chart with proportion of mono and bi allelic genes + nb of informative genes 
# adding allelism status info
allelicRatiosPerGeneStack$allelism <- if_else(allelicRatiosPerGeneStack$values >= 0.25, "bi", "mono", NA)

# counts genes 
counts_all_genes <- allelicRatiosPerGeneStack %>%
  filter(values != 'NA', chr == 'chrX') %>%
  group_by(cellTypeOrdered) %>%
  summarize(number_of_genes = n(), nb_mono = n()) 
counts_biall_genes <- allelicRatiosPerGeneStack %>%
  filter(values != 'NA', chr == 'chrX') %>%
  group_by(cellTypeOrdered) %>% 
  filter(allelism == 'bi') %>%
  summarize(nb_bi_genes = n())
counts_mono_genes <- allelicRatiosPerGeneStack %>%
  filter(values != 'NA', chr == 'chrX') %>%
  group_by(cellTypeOrdered) %>% 
  filter(allelism == 'mono') %>%
  summarize(nb_mono_genes = n())
counts <- data.frame(cell_type = counts_all_genes$cellTypeOrdered, nb_genes_tot = counts_all_genes$number_of_genes, nb_bi_genes = counts_biall_genes$nb_bi_genes, nb_mono_genes = count_mono_genes$nb_mono_genes)
counts$percent_bi <- counts$nb_bi_genes / counts$nb_genes_tot * 100
counts$percent_mono <- counts$nb_mono_genes / counts$nb_genes_tot * 100


# barplot
allelicRatiosPerGeneStack %>% 
  filter(values != 'NA', chr == 'chrX') %>%
  ggplot( aes(x = cellTypeOrdered, fill = allelism)) +
  geom_bar(stat = "count", , position = "fill") +
  labs(x = "",
       y = "Proportion of genes",
       fill = "Allelism") +
  theme_classic() 
#scale_y_continuous(limits = c(0, 1))
#geom_text(aes(x = Player, y = EFG,label=EFG),vjust=-0.25)+


# plot 4: map of biallelic genes in KD dox + that are not biallelic in KD dox - (subjected to xist-xci?)
list_gene_bi_KD <- allelicRatiosPerGeneStack[allelicRatiosPerGeneStack$cellTypeOrdered == 'KD dox +' & !is.na(allelicRatiosPerGeneStack$values) & allelicRatiosPerGeneStack$allelism == "bi" & allelicRatiosPerGeneStack$chr == 'chrX' ,]$gene 
list_gene_bi_KDnodox <- allelicRatiosPerGeneStack[allelicRatiosPerGeneStack$cellTypeOrdered == 'KD dox -' & !is.na(allelicRatiosPerGeneStack$values) & allelicRatiosPerGeneStack$allelism == "bi" & allelicRatiosPerGeneStack$chr == 'chrX' ,]$gene 
gene_list <- list_gene_bi_KD[!(list_gene_bi_KD %in% list_gene_bi_KDnodox)]
positions <- c()
for (gene_name in gene_list){
  start_pos <- min(allExons %>% filter(gene == gene_name) %>% select(start)) # takes the first start position int he list of all start positions of all exons from the gene
  positions <- as.vector(c(positions, start_pos))
}
length(positions) == length(gene_list)
gene_df <- data.frame(gene = gene_list, position=positions)

ggplot(data = gene_df, aes(x = position)) +
  # chrX extremities in black
  geom_segment(x = 0, y = 0, xend = 0, yend = 0.2, color = 'black') +
  geom_segment(x = 156040895, y = 0, xend = 156040895, yend = 0.2, color = 'black' ) +
  # centromere in blue
  geom_segment(x = 58100000, y = 0, xend = 58100000, yend = 0.2, color = 'blue' ) +
  geom_segment(x = 63800000, y = 0, xend = 63800000, yend = 0.2, color = 'blue' ) +
  # genes in red 
  geom_segment(aes(x = position, xend = position, y = 0, yend = 0.1), color = "red", linewidth = 0.05) +  # Draw red lines
  # scaling the axis
  ylim(0,0.2) +
  scale_x_continuous(limits = c(0, 156040895)) +  # Set x-axis limits
  theme_void()  # Remove unnecessary plot elements
