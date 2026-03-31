## Figure Sup.7D: expression heatmap

Input data files:
- trimmed bulk RNAseq fastq files (Naive, TSC, EXMC, Primed, and WT dox +/-, XIST iKD +/- dox)

### Transcript abundance using kallisto

Reads were pseudo-aligned using kallisto/0.46.2 [(Bray et al., 2016)](https://github.com/pachterlab/kallisto) to quantify transcripts abundance in TPM (Transcripts Per Kilobase Million).

First, kallisto index was generated as follow: 
1) transcriptome fasta files of cDNA and ncRNA were downloaded from: https://www.ensembl.org/info/data/ftp/index.html, and concatenated.
2) The resulting fasta file was indexed using: `kallisto index -i kallisto_transcripts.idx Homo_sapiens.GRCh38.cdna.ncrna.fa.gz`

Then reads were pseudo-aligned using this bash script:
```
#!/bin/bash

#SBATCH --partition=ipop-up
#SBATCH --job-name=kallisto
#SBATCH --mem=32G 
#SBATCH --cpus-per-task=8

module purge
module load kallisto/0.46.2

index=/path/to/kallisto/index/kallisto_transcripts.idx

kallisto quant -i $index -o ./D1688T25/ --rf-stranded -b 50 --threads 8 ../fastq/trimmed/D1688T25_trimmed_R1.fastq.gz ../fastq/trimmed/D1688T25_trimmed_R2.fastq.gz
kallisto quant -i $index -o ./D1688T26/ --rf-stranded -b 50 --threads 8 ../fastq/trimmed/D1688T26_trimmed_R1.fastq.gz ../fastq/trimmed/D1688T26_trimmed_R2.fastq.gz
kallisto quant -i $index -o ./D1688T27/ --rf-stranded -b 50 --threads 8 ../fastq/trimmed/D1688T27_trimmed_R1.fastq.gz ../fastq/trimmed/D1688T27_trimmed_R2.fastq.gz
kallisto quant -i $index -o ./D1688T28/ --rf-stranded -b 50 --threads 8 ../fastq/trimmed/D1688T28_trimmed_R1.fastq.gz ../fastq/trimmed/D1688T28_trimmed_R2.fastq.gz
kallisto quant -i $index -o ./EXMC/ --rf-stranded -b 50 --threads 8 ../fastq/trimmed/E_R1_val_1.fq.gz ../fastq/trimmed/E_R1_val_2.fq.gz
kallisto quant -i $index -o ./TSC/ --rf-stranded -b 50 --threads 8 ../fastq/trimmed/T4-ACE_R1_val_1.fq.gz ../fastq/trimmed/T4-ACE_R1_val_2.fq.gz
kallisto quant -i $index -o ./Naive/ --rf-stranded -b 50 --threads 8 ../fastq/trimmed/Naive_R1_val_1.fq.gz ../fastq/trimmed/Naive_R1_val_2.fq.gz
kallisto quant -i $index -o ./Primed/ --rf-stranded -b 50 --threads 8 ../fastq/trimmed/Primed_R1_val_1.fq.gz ../fastq/trimmed/Primed_R1_val_2.fq.gz
```

### Heatmap (R)

The output abundance.tsv files were renamed (adding sample name, to differenciate them from one another), gathered in a `kallisto_tpm_tables/` folder and then processed on R/4.3.2 to generate a heatmap using the following R script. This script requires a gtf file (here gencode.v31.annotation.gtf):

```
library(ggplot2)
library(dplyr)
library(tibble)
library(gridExtra)
library(viridis)
library(stringr)

# uploading data
samples <- c("Naive", "D1688T27", "D1688T28", "D1688T25", "D1688T26", "TSC", "EXMC", "Primed")

sample_data <- list()
for(sample_name in samples) {
  # Construct the file path for the sample
  file_path <- paste0("kallisto_tpm_tables/abundance_", sample_name, ".tsv")
  # Read the file corresponding to the sample, save it as a df in a list
  sample_data[[sample_name]] <- read.table(file = file_path, header = TRUE)
}
# Extract tpm from each sample and combine into tpm_matrix
tpm_matrix <- data.frame(transcript_id = sample_data[[samples[1]]]$target_id)
for(sample_name in samples) {
  tpm_matrix[[sample_name]] <- sample_data[[sample_name]]$tpm
}


#  Adding gene_name info
# gtf preparation to have a column with transcript ID and the corresponding gene_name -> takes several minutes so was downloaded once generated
gtf <- read.table("X:/ReferenceGenomes/GRCh38/GTF/gencode.v31.annotation.gtf", header = F, sep = "\t")
gtf <- gtf[gtf$V3 == 'transcript',]
# write.table(gtf, "gtf_trancripts.tsv", sep = "\t", col.names = F, row.names = F)
# gtf <- read.table("gtf_trancripts.tsv", sep = '\t', header = F)

# extract gene_name, gene_id and transcript ID info from gtf last column, into new columns, to link trID with gene name in kallisto output
extract_info <- function(last_column_str) {
  gene_name <- str_extract(last_column_str, "gene_name\\s([^;]+);")
  gene_id <- str_extract(last_column_str, "gene_id\\s([^;]+);")
  transcript_id <- str_extract(last_column_str, "transcript_id\\s([^;]+);")
  # Remove "gene_name" and final semicolon from each extraction
  gene_name <- gsub("gene_name\\s|;$", "", gene_name)
  gene_id <- gsub("gene_id\\s|;$", "", gene_id)
  transcript_id <- gsub("transcript_id\\s|;$", "", transcript_id)
  
  return(c(gene_name, gene_id, transcript_id))
}

gtf<- cbind(gtf, t(sapply(gtf$V9, extract_info)))
colnames(gtf)[10:12] <- c("gene_name", "gene_id", "transcript_id")

# Add gene_name info in tpm_matrix
tpm_matrix <- merge(tpm_matrix, gtf[,c("transcript_id", "gene_name")], by = "transcript_id", all.x = TRUE)
# Stack the tpm_matrix into a long dataframe
tpm_data_long <- data.frame(transcript_id = tpm_matrix$transcript_id, gene_name = tpm_matrix$gene_name, stack(tpm_matrix[, - c(1,10)]))
colnames(tpm_data_long)[3:4] <- c("tpm", "sample")
# Take the sum of transcript counts corresponding to the same gene (but different transcripts)
tpm_data_long <- tpm_data_long %>%
  group_by(sample, gene_name) %>%
  summarize(tpm = sum(tpm)) 

# Genes in categories (B. Balaton)
feats.core <- c("POU5F1", "SOX2","NANOG")
feats.primed <- c("OTX2", "ZIC2", "CD24", "DUSP6", "TCF4")
feats.naive <- c("KLF17", "KLF4", "SUSD2", "DNMT3L", "DPPA5", "TFCP2L1")
feats.TE <- c("GATA2","GATA3","ITGA6","TP63","KRT7","KRT18","HAND1","NR2F2")
feats.EXMC <- c("LUM","NID2", "FOXF1", "VIM", "POSTN", "ANXA1")#,  "PITX1")#, "BST2", "DCN")
feats.PrE <- c("SOX17", "GATA4", "GATA6", "FOXA2", "PDGFRA","CDH2")
feats.EVT <- c("HLA-G", "MMP2", "ITGA5")
feats.STB <- c("CGA", "CGB3", "SDC1")
feats.amnion<- c("WNT6", "GABRP", "ISL1", "HEY1", "CDH10", "CTSV", "TPM1")
feats.mesoderm <- c("MIXL1", "MESP1", "ZIC3", "CDX1", "CDX2", "CDX4", "EOMES")
feats.XCI <- c("XIST", "XACT")

combinedFeats <- list(core=feats.core, primed=feats.primed, naive=feats.naive, trophectoderm=feats.TE, EXMC=feats.EXMC, PE=feats.PrE, EVT=feats.EVT, STB=feats.STB, amnion=feats.amnion, mesoderm=feats.mesoderm, XCI=feats.XCI)

# ordering genes
temp <- as.factor(unlist(combinedFeats))
temp <- factor(temp, levels=temp)

tpm_data_long$geneOrder <- NA
for(i in 1:length(temp)){
  tpm_data_long$geneOrder[which(tpm_data_long$gene_name==temp[i])] <- i 
}



##### Plot #####

# order of samples 
tpm_data_long$sample <- factor(tpm_data_long$sample, levels = samples, labels = c('Naive', 'WT dox -', 'WT dox +', 'KD dox -', 'KD dox +', 'TSC', 'EXMC', 'Primed'))
# selecting only genes to plot
tpm_data_genesubset <- tpm_data_long[tpm_data_long$gene_name %in% c(feats.core, feats.naive, feats.primed, feats.TE, feats.EXMC),]

# plot heatmap

ggplot(tpm_data_genesubset) +
  geom_rect(aes(xmin=as.numeric(sample)-0.5, xmax=as.numeric(sample)+0.5, ymin=geneOrder-0.5, ymax=geneOrder+0.5, fill = tpm), linejoin = "bevel") +
  #display sample names
  geom_text(aes(x=as.numeric(sample), y=-1, label=sample)) +
  theme_linedraw() +
  scale_fill_viridis(option="A", trans = "log1p", breaks = c(0,10,25,50,100,250,500,750,1000)) +
  # display gene names 
  geom_text(aes(x=0, y=geneOrder, label=gene_name)) + 
  xlab("sample") +
  ylab("gene") +
  guides(fill=guide_colorbar(title="log(TPM + 1)")) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(plot.caption = element_text(size = 10))
```
