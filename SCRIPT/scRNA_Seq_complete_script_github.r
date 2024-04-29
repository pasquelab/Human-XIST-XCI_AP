library(Seurat)
library(dplyr)
library(tibble)
library(patchwork)
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(edgeR)

setwd ('/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/')

Day_0_Panda_Carrillo.data <- Read10X(data.dir = "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/Sequencing_files_matrix/Day-0/filtered/")

Day_0_Panda_Carrillo.df <- as.data.frame(Read10X(data.dir = "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/Sequencing_files_matrix/Day-0/filtered/"))

head(Day_0_Panda_Carrillo.df)
dim(Day_0_Panda_Carrillo.df)
tail(Day_0_Panda_Carrillo.df)

outputFolderPath <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filePath <- file.path(outputFolderPath, "Day_0_Panda_Carrillo.rds")
saveRDS(Day_0_Panda_Carrillo.df, file = filePath)

so_Day_0_Panda_Carrillo <- CreateSeuratObject(counts = Day_0_Panda_Carrillo.df, project = "Day_0_Panda_Carrillo", min.cells = 0, min.features = 0)
so_Day_0_Panda_Carrillo

n_last <- 10 
so_Day_0_Panda_Carrillo@meta.data$conditions<- substr(row.names(so_Day_0_Panda_Carrillo@meta.data), nchar(row.names(so_Day_0_Panda_Carrillo@meta.data)) - n_last + 1, nchar(row.names(so_Day_0_Panda_Carrillo@meta.data)))

lookup_table <- data.frame(Original = c("ACTTTAGG-1", "AACGGGAA-1", "AGTAGGCT-1", "ATGTTGAC-1"),
                            Replacement = c("WT_DOX-", "XIST_KD_DOX-", "WT_DOX+", "XIST_KD_DOX+"))

so_Day_0_Panda_Carrillo@meta.data$conditions <- lookup_table$Replacement[match(so_Day_0_Panda_Carrillo@meta.data$conditions, lookup_table$Original)]

so_Day_0_Panda_Carrillo@meta.data <- cbind(so_Day_0_Panda_Carrillo@meta.data, source = rep("DAY-0", nrow(so_Day_0_Panda_Carrillo@meta.data)))

outputFolderPath <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filePath <- file.path(outputFolderPath, "so_Day_0_Panda_Carrillo.rds")
saveRDS(so_Day_0_Panda_Carrillo, file = filePath)

Day_4_Panda_Carrillo.data <- Read10X(data.dir = "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/Sequencing_files_matrix/Day-4/filtered/")

Day_4_Panda_Carrillo.df <- as.data.frame(Read10X(data.dir = "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/Sequencing_files_matrix/Day-4/filtered/"))

outputFolderPath <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filePath <- file.path(outputFolderPath, "Day_4_Panda_Carrillo.rds")
saveRDS(Day_4_Panda_Carrillo.df, file = filePath)

Day_4_Panda_Carrillo.df <- readRDS("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/Day_4_Panda_Carrillo.rds")

so_Day_4_Panda_Carrillo <- CreateSeuratObject(counts = Day_4_Panda_Carrillo.df, project = "Day_4_Panda_Carrillo", min.cells = 0, min.features = 0)
so_Day_4_Panda_Carrillo

so_Day_4_Panda_Carrillo <- CreateSeuratObject(counts = Day_4_Panda_Carrillo.df, project = "Day_4_Panda_Carrillo", min.cells = 0, min.features = 0)
so_Day_4_Panda_Carrillo

lookup_table <- data.frame(Original = c("ACTTTAGG-1", "AACGGGAA-1", "AGTAGGCT-1", "ATGTTGAC-1"),
                            Replacement = c("WT_DOX-", "XIST_KD_DOX-", "WT_DOX+", "XIST_KD_DOX+"))

so_Day_4_Panda_Carrillo@meta.data$conditions <- lookup_table$Replacement[match(so_Day_4_Panda_Carrillo@meta.data$conditions, lookup_table$Original)]

so_Day_4_Panda_Carrillo@meta.data <- cbind(so_Day_4_Panda_Carrillo@meta.data, source = rep("DAY-4", nrow(so_Day_4_Panda_Carrillo@meta.data)))

outputFolderPath <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filePath <- file.path(outputFolderPath, "so_Day_4_Panda_Carrillo.rds")
saveRDS(so_Day_4_Panda_Carrillo, file = filePath)

Day_8_Panda_Carrillo.data <- Read10X(data.dir = "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/Sequencing_files_matrix/Day-8/filtered/")

Day_8_Panda_Carrillo.df <- as.data.frame(Read10X(data.dir = "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/Sequencing_files_matrix/Day-8/filtered/"))

outputFolderPath <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filePath <- file.path(outputFolderPath, "Day_8_Panda_Carrillo.rds")
saveRDS(Day_8_Panda_Carrillo.df, file = filePath)

so_Day_8_Panda_Carrillo <- CreateSeuratObject(counts = Day_8_Panda_Carrillo.df, project = "Day_8_Panda_Carrillo", min.cells = 0, min.features = 0)
so_Day_8_Panda_Carrillo

n_last <- 10 
so_Day_8_Panda_Carrillo@meta.data$conditions<- substr(row.names(so_Day_8_Panda_Carrillo@meta.data), nchar(row.names(so_Day_8_Panda_Carrillo@meta.data)) - n_last + 1, nchar(row.names(so_Day_8_Panda_Carrillo@meta.data)))

lookup_table <- data.frame(Original = c("ACTTTAGG-1", "AACGGGAA-1", "AGTAGGCT-1", "ATGTTGAC-1"),
                            Replacement = c("WT_DOX-", "XIST_KD_DOX-", "WT_DOX+", "XIST_KD_DOX+"))
so_Day_8_Panda_Carrillo@meta.data$conditions <- lookup_table$Replacement[match(so_Day_8_Panda_Carrillo@meta.data$conditions, lookup_table$Original)]

so_Day_8_Panda_Carrillo@meta.data <- cbind(so_Day_8_Panda_Carrillo@meta.data, source = rep("DAY-8", nrow(so_Day_8_Panda_Carrillo@meta.data)))

outputFolderPath <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filePath <- file.path(outputFolderPath, "so_Day_8_Panda_Carrillo.rds")
saveRDS(so_Day_8_Panda_Carrillo, file = filePath)

Day_0_Panda_Carrillo <- readRDS("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/Day_0_Panda_Carrillo.rds")

Day_4_Panda_Carrillo <- readRDS("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/Day_4_Panda_Carrillo.rds")

Day_8_Panda_Carrillo <- readRDS("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/Day_8_Panda_Carrillo.rds")

Day_0_Panda_Carrillo$genes <- rownames(Day_0_Panda_Carrillo)
rownames(Day_0_Panda_Carrillo) <- NULL
Day_4_Panda_Carrillo$genes <- rownames(Day_4_Panda_Carrillo)
rownames(Day_4_Panda_Carrillo) <- NULL
Day_8_Panda_Carrillo$genes <- rownames(Day_8_Panda_Carrillo)
rownames(Day_8_Panda_Carrillo) <- NULL

Day_0_4_Panda_Carrillo <- merge(Day_0_Panda_Carrillo, Day_4_Panda_Carrillo,
                          by = 'genes', all = TRUE)

Day_0_4_8_Panda_Carrillo <- merge(Day_0_4_Panda_Carrillo, Day_8_Panda_Carrillo,
                          by = 'genes', all = TRUE)

row.names(Day_0_4_8_Panda_Carrillo) <- Day_0_4_8_Panda_Carrillo$genes

column_to_remove <- "genes"  
Day_0_4_8_Panda_Carrillo <- Day_0_4_8_Panda_Carrillo[, -which(colnames(Day_0_4_8_Panda_Carrillo) == column_to_remove)]


Day_0_4_8_Panda_Carrillo[is.na(Day_0_4_8_Panda_Carrillo)] <- 0

outputFolderPath <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filePath <- file.path(outputFolderPath, "Day_0_4_8_Panda_Carrillo.rds")
saveRDS(Day_0_4_8_Panda_Carrillo, file = filePath)

so_Day_0_4_8_Panda_Carrillo <- CreateSeuratObject(counts = Day_0_4_8_Panda_Carrillo, project = "Day_0_4_8_Panda_Carrillo", min.cells = 0, min.features = 0)
so_Day_0_4_8_Panda_Carrillo

outputFolderPath <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filePath <- file.path(outputFolderPath, "so_Day_0_4_8_Panda_Carrillo.rds")
saveRDS(so_Day_0_4_8_Panda_Carrillo, file = filePath)

so_Day_0_4_8_Panda_Carrillo[["percent.mt"]] <- PercentageFeatureSet(so_Day_0_4_8_Panda_Carrillo, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(so_Day_0_4_8_Panda_Carrillo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

so_Day_0_4_8_Panda_Carrillo_filtered <- subset(so_Day_0_4_8_Panda_Carrillo, subset = nFeature_RNA > 1500 & nCount_RNA > 2200 & percent.mt < 25 )

so_Day_0_Panda_Carrillo <- readRDS("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/so_Day_0_Panda_Carrillo.rds")

so_Day_4_Panda_Carrillo <- readRDS("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/so_Day_4_Panda_Carrillo.rds")

so_Day_8_Panda_Carrillo <- readRDS("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/so_Day_8_Panda_Carrillo.rds")

# Subsetting the original matrix based on row names and columns 4 and 5
so_Day_0_Panda_Carrillo_subset_metadata <- so_Day_0_Panda_Carrillo@meta.data[, c(4, 5)]
so_Day_0_Panda_Carrillo_subset_metadata

outputFolderPath <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filePath <- file.path(outputFolderPath, "so_Day_0_Panda_Carrillo_subset_metadata.rds")
saveRDS(so_Day_0_Panda_Carrillo_subset_metadata, file = filePath)

# Subsetting the original matrix based on row names and columns 4 and 5
so_Day_4_Panda_Carrillo_subset_metadata <- so_Day_4_Panda_Carrillo@meta.data[, c(4, 5)]
so_Day_4_Panda_Carrillo_subset_metadata

outputFolderPath <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filePath <- file.path(outputFolderPath, "so_Day_4_Panda_Carrillo_subset_metadata.rds")
saveRDS(so_Day_4_Panda_Carrillo_subset_metadata, file = filePath)

# Subsetting the original matrix based on row names and columns 4 and 5
so_Day_8_Panda_Carrillo_subset_metadata <- so_Day_8_Panda_Carrillo@meta.data[, c(4, 5)]
so_Day_8_Panda_Carrillo_subset_metadata

outputFolderPath <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filePath <- file.path(outputFolderPath, "so_Day_8_Panda_Carrillo_subset_metadata.rds")
saveRDS(so_Day_8_Panda_Carrillo_subset_metadata, file = filePath)

so_Day_0_4_8_Panda_Carrillo_subset_metadata <- rbind(so_Day_0_Panda_Carrillo_subset_metadata, so_Day_4_Panda_Carrillo_subset_metadata, so_Day_8_Panda_Carrillo_subset_metadata)

head(so_Day_0_4_8_Panda_Carrillo_subset_metadata)
dim(so_Day_0_4_8_Panda_Carrillo_subset_metadata)
tail(so_Day_0_4_8_Panda_Carrillo_subset_metadata)

so_Day_0_4_8_Panda_Carrillo_norm@meta.data <- cbind(so_Day_0_4_8_Panda_Carrillo_norm@meta.data, so_Day_0_4_8_Panda_Carrillo_subset_metadata[, c("conditions", "source")])

so_Day_0_4_8_Panda_Carrillo_subset_metadata <- readRDS("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/so_Day_0_4_8_Panda_Carrillo_subset_metadata.rds")

# Create a new column by combining "conditions" and "source"
so_Day_0_4_8_Panda_Carrillo_norm@meta.data$combined <- paste(so_Day_0_4_8_Panda_Carrillo_norm@meta.data$source, "_", so_Day_0_4_8_Panda_Carrillo_norm@meta.data$conditions, sep = "")

so_Day_0_4_8_Panda_Carrillo_filtered_norm <- readRDS("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/so_Day_0_4_8_Panda_Carrillo_filtered_norm_2.rds")

# Columns to remove
columns_to_remove <- c(8)
# Remove specific columns
so_Day_0_4_8_Panda_Carrillo_filtered_norm@meta.data <- so_Day_0_4_8_Panda_Carrillo_filtered_norm@meta.data[, -columns_to_remove]
head(so_Day_0_4_8_Panda_Carrillo_filtered_norm@meta.data)

outputFolderPath <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filePath <- file.path(outputFolderPath, "so_Day_0_4_8_Panda_Carrillo_norm.rds")
saveRDS(so_Day_0_4_8_Panda_Carrillo_norm, file = filePath)

so_Day_0_4_8_Panda_Carrillo_filtered <- readRDS("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/so_Day_0_4_8_Panda_Carrillo_filtered.rds")

so_Day_0_4_8_Panda_Carrillo_filtered_norm <- NormalizeData(so_Day_0_4_8_Panda_Carrillo_filtered, normalization.method = "LogNormalize", scale.factor = 10000)

so_Day_0_4_8_Panda_Carrillo_filtered_norm <- FindVariableFeatures(so_Day_0_4_8_Panda_Carrillo_filtered_norm, selection.method = "vst", nfeatures = 5000)

top10 <- head(VariableFeatures(so_Day_0_4_8_Panda_Carrillo_filtered_norm), 10)
print(top10)

plot1 <- VariableFeaturePlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm)
plot1

output_path <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filename <- "plot1_d_0_4_8_filtered.pdf"

ggsave(file.path(output_path, filename), plot = plot1, height = 8, width = 8)

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

output_path <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filename <- "plot2_d_0_4_8_filtered.pdf"

ggsave(file.path(output_path, filename), plot = plot2, height = 8, width = 8)

all.genes <- rownames(so_Day_0_4_8_Panda_Carrillo_filtered_norm)
so_Day_0_4_8_Panda_Carrillo_filtered_norm <- ScaleData(so_Day_0_4_8_Panda_Carrillo_filtered_norm, features = all.genes)


so_Day_0_4_8_Panda_Carrillo_filtered_norm <- RunPCA(so_Day_0_4_8_Panda_Carrillo_filtered_norm, features = VariableFeatures(object = so_Day_0_4_8_Panda_Carrillo_filtered_norm), npcs=60,nfeatures.print=10,ndims.print=1:5,seed.use=13)

print(so_Day_0_4_8_Panda_Carrillo_filtered_norm[["pca"]], dims = 1:5, nfeatures = 5)

VIZDIMPLOT<- VizDimLoadings(so_Day_0_4_8_Panda_Carrillo_filtered_norm, dims = 1:2, reduction = "pca")
VIZDIMPLOT

output_path <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filename <- "VIZDIMPLOT_filtered.pdf"

ggsave(file.path(output_path, filename), plot = VIZDIMPLOT, height = 8, width = 8)

DIMPLOT_D0_4_8<- DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "pca")
DIMPLOT_D0_4_8

output_path <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filename <- "DIMPLOT_D0_4_8_filtered.pdf"

ggsave(file.path(output_path, filename), plot = DIMPLOT_D0_4_8, height = 8, width = 8)

ELBOWPLOT_D_0_4_8<- ElbowPlot(
 so_Day_0_4_8_Panda_Carrillo_filtered_norm,
  ndims = 60)
ELBOWPLOT_D_0_4_8

output_path <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filename <- "ELBOWPLOT_D_0_4_8_filtered.pdf"

ggsave(file.path(output_path, filename), plot = ELBOWPLOT_D_0_4_8, height = 8, width = 8)

outputFolderPath <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/"
filePath <- file.path(outputFolderPath, "so_Day_0_4_8_Panda_Carrillo_filtered_norm.rds")
saveRDS(so_Day_0_4_8_Panda_Carrillo_filtered_norm, file = filePath)

DimHeatmap(so_Day_0_4_8_Panda_Carrillo_filtered_norm, dims = 1:10, cells = 500, balanced = TRUE)

use.pcs = 1:10

so_Day_0_4_8_Panda_Carrillo_filtered_norm <- FindNeighbors(so_Day_0_4_8_Panda_Carrillo_filtered_norm, dims = use.pcs)
so_Day_0_4_8_Panda_Carrillo_filtered_norm <- FindClusters(so_Day_0_4_8_Panda_Carrillo_filtered_norm, resolution = 0.5)

so_Day_0_4_8_Panda_Carrillo_filtered_norm <- RunUMAP(so_Day_0_4_8_Panda_Carrillo_filtered_norm, dims = 1:10)

options(repr.plot.width = 11, repr.plot.height = 10, repr.plot.res = 150)
UMAP_Day_0_4_8_Panda_Carrillo_combined_conditions_filtered <- DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", label=T, group.by = "combined", raster=FALSE)
UMAP_Day_0_4_8_Panda_Carrillo_combined_conditions_filtered

output_path <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/OUTPUT/filtered_full/Final_assay_after_removing_vimcells"
filename <- "UMAP_Day_0_4_8_Panda_Carrillo_cluster_dims_10_res_point5_seurat_clusters.pdf"

ggsave(file.path(output_path, filename), plot = UMAP_Day_0_4_8_Panda_Carrillo_combined_conditions_filtered, height = 8, width = 8)

options(repr.plot.width = 36, repr.plot.height = 36, repr.plot.res = 150)
Featureplot_umap_geneexpression <- FeaturePlot(cols = c('grey',(RColorBrewer::brewer.pal(9,'Reds'))),
  object = so_Day_0_4_8_Panda_Carrillo_filtered_norm,
  features = c('GATA3', 'GATA2','HAND1', 'KRT7', 'KRT19', 'NR2F2',
               'SPIC', 'KLF4','KLF17','VENTX','SOX2', 'SUSD2', 'DNMT3L', 'TFCP2L1', 'ENPEP', 'TACSTD2',
               'ZIC2', 'OTX2', 'DUSP6',
              'VIM', 'POSTN', 'LUM', 'NID2', 'DCN', 'ANXA1', 'GATA4', 'BST2', 'COL1A1', 'CDH2', 'CDX2', 'COL6A1', 'LAMC1', 'GATA6'),
  pt.size = .5,
  max.cutoff = 'q75',
  ncol = 3,
                                               reduction = "umap",
)
Featureplot_umap_geneexpression

so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster<- as.character(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters)

so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster[which(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters=="0" )]<- "Naive"
so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster[which(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters=="1" )]<- "EXMC"
so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster[which(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters=="2" )]<- "Naive"
so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster[which(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters=="3" )]<- "Naive"
so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster[which(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters=="4" )]<- "Epiblast_interm_1"
so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster[which(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters=="5" )]<- "Epiblast_interm_3"
so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster[which(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters=="6" )]<- "Epiblast_interm_2"
so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster[which(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters=="7" )]<- "Epiblast_interm_2"
so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster[which(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters=="8" )]<- "Naive"
so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster[which(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters=="9" )]<- "Epiblast_interm_4"
so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster[which(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters=="10" )]<- "EXMC"
so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster[which(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters=="11" )]<- "Differ_Naive"
so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster[which(so_Day_0_4_8_Panda_Carrillo_filtered_norm$seurat_clusters=="12" )]<- "TSC"

Idents(so_Day_0_4_8_Panda_Carrillo_filtered_norm) <- "anno.cluster"

options(repr.plot.width = 8, repr.plot.height = 8, repr.plot.res = 150)
cluster_colors <- c("#522D68", "#F26522", "#DF3383", "#2D419A", "#00B1F0", "#005FAD", "#0089CE", "#39B54A")
UMAP_so_Day_0_4_8_Panda_Carrillo_filtered_norm_cluster <- DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", cols = cluster_colors, label=T, raster=FALSE)
UMAP_so_Day_0_4_8_Panda_Carrillo_filtered_norm_cluster

output_path <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/OUTPUT/filtered_full/Final_assay_after_removing_vimcells"
filename <- "UMAP_Day_0_4_8_Panda_Carrillo_cluster_dims_10_res_point5_annot_cluster.pdf"

ggsave(file.path(output_path, filename), plot = UMAP_so_Day_0_4_8_Panda_Carrillo_filtered_norm_cluster, height = 8, width = 8)

library(viridis)

viridis_palette <- viridis(9)

options(repr.plot.width = 14, repr.plot.height = 6, repr.plot.res = 150)
Featureplot_umap_geneexpression <- FeaturePlot(
  cols = c('grey', viridis_palette),
  object = so_Day_0_4_8_Panda_Carrillo_filtered_norm,
  features = c('KLF4', 'DNMT3L', 'GATA3', 'GATA2', 'LUM', 'NID2'),
  pt.size = .8,
  max.cutoff = 'q75',
  ncol = 3,
  reduction = "umap"
)
Featureplot_umap_geneexpression

output_path <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/"
filename <- "Featureplot_umap_geneexpression_dims_10_res_point5_final_selected gene.pdf"

ggsave(file.path(output_path, filename), plot = Featureplot_umap_geneexpression, height = 18, width = 20)

signature_markers<-c("ITGA2","EPCAM","VGLL1", "ADAP2","KRT7", "KRT18", "CLDN4","ARID3A","FGF4", "HAND1", "GATA2", "GATA3","DAB2","KRT19", "TFAP2A", "ENPEP", "TACSTD2", #TE
                        "DNMT3L","SUSD2", "SPIC", "VENTX", "KLF17", "KLF4", #Naive
                        "NANOG", "POU5F1", "ZIC2", "SOX2", # corepluripotency
                    "BST2", "VIM", "GATA4", "LUM", "ANXA1", "COL3A1", "SNAI2", "NID2", "DCN", "GATA6", "FN1") # EXMC

options(repr.plot.width=15, repr.plot.height=8)
dotplot_filtered<- DotPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm_2,features = signature_markers,
               dot.scale = 8) + 
    RotatedAxis() +
    scale_colour_viridis()
dotplot_filtered

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/dotplot_dims10_res_points5_final_anno_cluster.pdf", width = 15, height = 8)
plot(dotplot_filtered)
dev.off()

library(gridExtra)
library(ggplot2)
library(wacolors)

options(repr.plot.width=25, repr.plot.height=5)
number_perSample<- table(so_Day_0_4_8_Panda_Carrillo_filtered_norm$anno.cluster, so_Day_0_4_8_Panda_Carrillo_filtered_norm$combined)
grid.table(number_perSample)

DF_Count <- so_Day_0_4_8_Panda_Carrillo_filtered_norm@meta.data %>%group_by(combined,anno.cluster) %>%
  count() %>%
  ungroup() %>%
  group_by(combined) %>%
  mutate(Freq = n/sum(n)*100)
options(repr.plot.width=20, repr.plot.height=10)
p<- ggplot(DF_Count, aes(x = combined, y = Freq, fill = anno.cluster))+
  geom_col()+
  geom_text(aes(label = paste(round(Freq, 1),"%")),position = position_stack(vjust = 0.5))+
    scale_fill_wa_d(wacolors$rainier)+ 
    theme_classic()
p

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Day_0_4_8_Panda_Carrillo_table_cell_type_percentage_filtered_2_final.pdf", width = 20, height = 10)
plot(p)
dev.off()

so_Day_0_4_8_Panda_Carrillo_filtered_norm<- readRDS("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/so_Day_0_4_8_Panda_Carrillo_filtered_norm.rds")

table(so_Day_0_4_8_Panda_Carrillo_filtered_norm_2@meta.data$combined)

Idents(so_Day_0_4_8_Panda_Carrillo_filtered_norm_2) <- "combined"
DEGs_DAY8_XIST_DOX_vs_noDOX <- FindMarkers(so_Day_0_4_8_Panda_Carrillo_filtered_norm_2, ident.1 = "DAY-8_XIST_KD_DOX+", ident.2 = "DAY-8_XIST_KD_DOX-", verbose = FALSE)
head(DEGs_DAY8_XIST_DOX_vs_noDOX)
tail(DEGs_DAY8_XIST_DOX_vs_noDOX)

dim(DEGs_DAY8_XIST_DOX_vs_noDOX)

row_names <- rownames(DEGs_DAY8_XIST_DOX_vs_noDOX)
DEGs_DAY8_XIST_DOX_vs_noDOX$genes <- row_names

write.csv(DEGs_DAY8_XIST_KD_DOX_vs_noDOX, "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/all_DEGs_DAY8_WT_DOX_vs_noDOX_final.csv")

summary(DEGs_DAY8_XIST_KD_DOX_vs_noDOX)

library(ggplot2)
library(ggrepel)
library(dplyr)

topDEGs_DEGs_DAY8_XIST_DOX_vs_noDOX<- DEGs_DAY8_XIST_DOX_vs_noDOX %>%
                dplyr::filter(p_val_adj < 0.05) %>%
                arrange(desc(abs(avg_log2FC)))
        
head(topDEGs_DAY8_XIST_DOX_vs_noDOX)
dim(topDEGs_DAY8_XIST_DOX_vs_noDOX)

write.csv(topDEGs_DAY8_XIST_DOX_vs_noDOX, "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/topDEGs_DAY8_XIST_DOX_vs_noDOX_final.csv")

library(ggrepel)

options(repr.plot.height = 8, repr.plot.width = 8, ggrepel.max.overlaps = Inf)

p_DEGs_DAY8_XIST_DOX_vs_noDOX<-DEGs_DAY8_XIST_DOX_vs_noDOX %>%
  ggplot(aes(avg_log2FC, -log10(p_val_adj))) +
  geom_point(size = 0.5, aes(color = case_when(
    avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "dodgerblue4",
    avg_log2FC <= -0.5 & p_val_adj < 0.05 ~ "tomato",
    TRUE ~ "grey70"
  ))) +
  geom_text_repel(
    aes(label = genes, color = case_when(
      avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "dodgerblue4",
      avg_log2FC <= -0.5 & p_val_adj < 0.05 ~ "tomato",
      TRUE ~ "grey70"
    )),
    data = dplyr::filter(DEGs_DAY8_XIST_DOX_vs_noDOX, genes %in% topDEGs_DAY8_XIST_DOX_vs_noDOX$genes),
    segment.size = 0.2,
    point.padding = unit(0.3, "lines"),
    box.padding = unit(0.6, "lines")
  ) +
  geom_hline(yintercept = -log(0.05, base=100), linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey70") +
  scale_color_identity() +  # To use the specified colors
  ggtitle("Differentially Expressed Genes XISTKD_DOX_VS_NODOX") +
  theme_bw() +
coord_cartesian(xlim=c(-1,1), ylim=c(0,300))
#ylim( c(0, 30000000000000000000)) # if you want a specific scale to your plot

p_DEGs_DAY8_XIST_DOX_vs_noDOX


pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/volcano_p_DEGs_DAY8_XIST_DOX_vs_noDOX.pdf", width = 8, height = 8)
plot(p_DEGs_DAY8_XIST_DOX_vs_noDOX)
dev.off()

Idents(so_Day_0_4_8_Panda_Carrillo_filtered_norm_2) <- "combined"

options(repr.plot.width = 11, repr.plot.height = 10, repr.plot.res = 150)


UMAP_Day_0_4_8_Panda_Carrillo_combined_conditions_filtered <- DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm_2, reduction = "umap", label=T, group.by = "combined", raster=FALSE)
UMAP_Day_0_4_8_Panda_Carrillo_combined_conditions_filtered

output_path <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/OUTPUT/filtered_full/Final_assay_after_removing_vimcells"
filename <- "UMAP_Day_0_4_8_Panda_Carrillo_cluster_dims_10_res_point5_combined_conditions.pdf"

ggsave(file.path(output_path, filename), plot = UMAP_Day_0_4_8_Panda_Carrillo_combined_conditions_filtered, height = 8, width = 8)

options(repr.plot.width = 11, repr.plot.height = 10, repr.plot.res = 150)


UMAP_Day_0_4_8_Panda_Carrillo_combined_conditions_filtered <- DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm_2, reduction = "umap", label=T, group.by = "combined", raster=FALSE)
UMAP_Day_0_4_8_Panda_Carrillo_combined_conditions_filtered

output_path <- "/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/OUTPUT/filtered_full/Final_assay_after_removing_vimcells"
filename <- "UMAP_Day_0_4_8_Panda_Carrillo_cluster_dims_10_res_point5_combined_conditions.pdf"

ggsave(file.path(output_path, filename), plot = UMAP_Day_0_4_8_Panda_Carrillo_combined_conditions_filtered, height = 8, width = 8)

Idents(so_Day_0_4_8_Panda_Carrillo_filtered_norm_2) <- "combined"

DAY_0_WT_DOX_pos <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("DAY-0_WT_DOX+"))
DAY_0_WT_DOX_neg<- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("DAY-0_WT_DOX-"))
DAY_0_XIST_KD_DOX_pos<- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("DAY-0_XIST_KD_DOX+"))
DAY_0_XIST_KD_DOX_neg<- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("DAY-0_XIST_KD_DOX-"))
DAY_4_WT_DOX_neg <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("DAY-4_WT_DOX-"))
DAY_4_WT_DOX_pos <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("DAY-4_WT_DOX+"))
DAY_4_XIST_KD_DOX_neg <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("DAY-4_XIST_KD_DOX-"))
DAY_4_XIST_KD_DOX_pos <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("DAY-4_XIST_KD_DOX+"))
DAY_8_WT_DOX_neg <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("DAY-8_WT_DOX-"))
DAY_8_WT_DOX_pos <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("DAY-8_WT_DOX+"))
DAY_8_XIST_KD_DOX_neg <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("DAY-8_XIST_KD_DOX-"))
DAY_8_XIST_KD_DOX_pos <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("DAY-8_XIST_KD_DOX+"))

options(repr.plot.width=10, repr.plot.height=8)
pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_DAY_0_WT_DOX_neg.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(DAY_0_WT_DOX_neg), cols.highlight = c("#F8766D") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","DAY_0_WT_DOX_neg"), values = c("grey80","#F8766D")) 
dev.off()

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_DAY_0_WT_DOX_pos.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(DAY_0_WT_DOX_pos), cols.highlight = c("#DE8C00") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","DAY_0_WT_DOX_pos"), values = c("grey80","#DE8C00")) 
dev.off()

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_DAY_0_XIST_KD_DOX_neg.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(DAY_0_XIST_KD_DOX_neg), cols.highlight = c("#B79F00") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","DAY_0_XIST_KD_DOX_neg"), values = c("grey80","#B79F00")) 
dev.off()

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_DAY_0_XIST_KD_DOX_pos.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(DAY_0_XIST_KD_DOX_pos), cols.highlight = c("#7CAE00") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","DAY_0_XIST_KD_DOX_pos"), values = c("grey80","#7CAE00")) 
dev.off()

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_DAY_4_WT_DOX_neg.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(DAY_4_WT_DOX_neg), cols.highlight = c("#00BA38") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","DAY_4_WT_DOX_neg"), values = c("grey80","#00BA38")) 
dev.off()

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_DAY_4_WT_DOX_pos.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(DAY_4_WT_DOX_pos), cols.highlight = c("#00C08B") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","DAY_4_WT_DOX_pos"), values = c("grey80","#00C08B")) 
dev.off()

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_DAY_4_XIST_KD_DOX_neg.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(DAY_4_XIST_KD_DOX_neg), cols.highlight = c("#00BFC4") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","DAY_4_XIST_KD_DOX_neg"), values = c("grey80","#00BFC4")) 
dev.off()

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_DAY_4_XIST_KD_DOX_pos.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(DAY_4_XIST_KD_DOX_pos), cols.highlight = c("#00B4F0") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","DAY_4_XIST_KD_DOX_pos"), values = c("grey80","#00B4F0")) 
dev.off()

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_DAY_8_WT_DOX_neg.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(DAY_8_WT_DOX_neg), cols.highlight = c("#619CFF") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","DAY_8_WT_DOX_neg"), values = c("grey80","#619CFF")) 
dev.off()

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_DAY_8_WT_DOX_pos.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(DAY_8_WT_DOX_pos), cols.highlight = c("#C77CFF") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","DAY_8_WT_DOX_pos"), values = c("grey80","#C77CFF")) 
dev.off()

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_DAY_8_XIST_KD_DOX_neg.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(DAY_8_XIST_KD_DOX_neg), cols.highlight = c("#F564E3") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","DAY_8_XIST_KD_DOX_neg"), values = c("grey80","#F564E3")) 
dev.off()

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_DAY_8_XIST_KD_DOX_pos.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(DAY_8_XIST_KD_DOX_pos), cols.highlight = c("#FF64B0") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","DAY_8_XIST_KD_DOX_pos"), values = c("grey80","#FF64B0")) 
dev.off()

Idents(so_Day_0_4_8_Panda_Carrillo_filtered_norm_2) <- "anno.cluster"

Differ_Naive <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("Differ_Naive"))
Epiblast_interm_1<- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("Epiblast_interm_1"))
Epiblast_interm_2<- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("Epiblast_interm_2"))
Epiblast_interm_3<- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("Epiblast_interm_3"))
Epiblast_interm_4 <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("Epiblast_interm_4"))
EXMC <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("EXMC"))
Naive <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("Naive"))
TSC  <- WhichCells(so_Day_0_4_8_Panda_Carrillo_filtered_norm, idents = c("TSC"))

options(repr.plot.width=10, repr.plot.height=8)
pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_Naive_color.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(Naive), cols.highlight = c("#522D68") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","Naive"), values = c("grey80","#522D68")) 
dev.off()
options(repr.plot.width=10, repr.plot.height=8)
pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_Differ_Naive_color.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(Differ_Naive), cols.highlight = c("#F26522") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","Differ_Naive"), values = c("grey80","#F26522")) 
dev.off()
options(repr.plot.width=10, repr.plot.height=8)
pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_Epi_interm1_color.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(Epiblast_interm_1), cols.highlight = c("#2D419A") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","Epiblast_interm_1"), values = c("grey80","#2D419A")) 
dev.off()
options(repr.plot.width=10, repr.plot.height=8)
pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_Epiblast_interm_2_color.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(Epiblast_interm_2), cols.highlight = c("#005FAD") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","Epiblast_interm_2"), values = c("grey80","#005FAD")) 
dev.off()
options(repr.plot.width=10, repr.plot.height=8)
pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_Epiblast_interm_3_color.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(Epiblast_interm_3), cols.highlight = c("#0089CE") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","Epiblast_interm_3"), values = c("grey80","#0089CE")) 
dev.off()
options(repr.plot.width=10, repr.plot.height=8)
pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_Epiblast_interm_4_color.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(Epiblast_interm_4), cols.highlight = c("#00B1F0") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","Epiblast_interm_4"), values = c("grey80","#00B1F0")) 
dev.off()
options(repr.plot.width=10, repr.plot.height=8)
pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_TSC_color.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(TSC), cols.highlight = c("#DF3383") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","TSC"), values = c("grey80","#DF3383")) 
dev.off()
options(repr.plot.width=10, repr.plot.height=8)
pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2024_XIST_KD/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/Projection_UMAP/projection_EXMC_color.pdf", width=6, height=4)
DimPlot(so_Day_0_4_8_Panda_Carrillo_filtered_norm, reduction = "umap", raster=FALSE, cells.highlight= list(EXMC), cols.highlight = c("#39B54A") , pt.size =2.5, label.size = 3, cols= "grey") + 
  scale_color_manual(labels = c("unselected","EXMC"), values = c("grey80","#39B54A")) 
dev.off()



gene_list <-  c( "POU5F1","SOX2", "NANOG", #core
               "OTX2", "ZIC2", "CD24","DUSP6", "TCF4", #primed
                "KLF17", "KLF4", "SUSD2", "SPIC", "VENTX", "DNMT3L","DPPA5", "TFCP2L1",   #naive
               "GATA2","GATA3", "ITGA6", "TP63", "KRT7", "KRT18", "HAND1","NR2F2",#trophectoderm
             "LUM", "NID2", "FOXF1", "HAND1", "VIM", "POSTN","ANXA1", "PITX1", # EXM
                "SOX17", "GATA4", "GATA6","FOXA2", "PDGFRA", "CDH2",#prE
                "HLA-G", "MMP2",  #EVT
                "CGA", "CGB3", "SDC1", "CK7", #ST 
                "WNT6", "GABRP", "ISL1", "HEY1", "HAND1", "CDH10","CTSV","TPM1", #Amnion  
                "MIXL1", "MESP1", "EPHA4", "ZIC3", "GSC", "TBXT", "FOXF1", "HAND1", "CDX1", "CDX2", "CDX4","EOMES", #Mesoderm
                "PRAMEF1", "MBD3L2", "MBD3L3", "H3Y1", "TPRX1", "DPPA3", "CCNA1", "TRIM43", "TRIM49", #8CLC
                "DPPA2", "GDF3", "ZNF208", "ZNF528", "ZNF560", "ZNF76", "ZNF728", "ZNF729", "ZNF98"#Formative
            )

options(repr.plot.width=25
    , repr.plot.height=15)
heatmap<- DoHeatmap(so_Day_0_4_8_Panda_Carrillo_filtered_norm, group.by="anno.cluster", features = gene_list, size = 5)+ scale_fill_viridis(option = "A")

pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/heatmap_DAY_0_4_8_PANDA_CARRILLO_filtered_2_anno_cluster_final.pdf", width = 60, height = 30)
plot(heatmap)
dev.off()

library(ggplot2)
library(Seurat)
library(reshape2)
library(umap)
library(dplyr)
library(gridExtra)
library(tidyr)
library(RColorBrewer)
library(viridis)

setwd("/staging/leuven/stg_00041/Bradley/XCI.Robjects/")

bigSeurat <- readRDS("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/OUTPUT/filtered_full/Final_assay_after_removing_vimcells/so_Day_0_4_8_Panda_Carrillo_filtered_norm_2_final.rds")

summary(bigSeurat@meta.data)
unique(bigSeurat@meta.data$combined)

utils:::format.object_size(object.size(bigSeurat), "auto")

UCSCgeneAnno <- read.table("../genesByChr.hg38.tsv")
colnames(UCSCgeneAnno) <- c("chr", "gene")

#adding a bit to give pseudo-autosomal genes their own classification. just doing a basic of genes mentioned on the X and Y
PARgenes <- unique(UCSCgeneAnno$gene)[which(
    unique(UCSCgeneAnno$gene)%in%UCSCgeneAnno$gene[which(UCSCgeneAnno$chr=="chrX")] &
    unique(UCSCgeneAnno$gene)%in%UCSCgeneAnno$gene[which(UCSCgeneAnno$chr=="chrY")]
)]

PARgenes
UCSCgeneAnno$chr[which(UCSCgeneAnno$gene%in%PARgenes)] <- "PAR"

chrAnno <- data.frame(gene=rownames(bigSeurat@assays$RNA@counts), chr=NA)

for(i in unique(UCSCgeneAnno$chr)){
 chrAnno$chr[which(chrAnno$gene%in%UCSCgeneAnno$gene[which(UCSCgeneAnno$chr==i)])]<- i   
}

unique(bigSeurat@meta.data$combined)

#initiating new column names. this is needed when you subset in the next cells
bigSeurat@meta.data$chrXmeans <- NA
bigSeurat@meta.data$chrYmeans<- NA
bigSeurat@meta.data$chr2means<- NA
bigSeurat@meta.data$chr7means<- NA
bigSeurat@meta.data$chrAllmeans<- NA

#the sparse matrix is too big to run means on because the function makes a dense matrix which takes up too much space and crashes the session.
#I can get around this by running fewer chromosomes and only including one condition at a time
for(condition in unique(bigSeurat@meta.data$combined)){
    #subsetting 1 chromosome at a time and 1 condition at a time
    bigSeurat@meta.data$chrXmeans[which(bigSeurat@meta.data$combined==condition)] <- as.numeric(apply(bigSeurat@assays$RNA@counts[which(chrAnno$chr%in%"chrX"),
                                                                                          which(bigSeurat@meta.data$combined==condition)], 2, mean))
    bigSeurat@meta.data$chrYmeans[which(bigSeurat@meta.data$combined==condition)] <- as.numeric(apply(bigSeurat@assays$RNA@counts[which(chrAnno$chr%in%"chrY"),
                                                                                          which(bigSeurat@meta.data$combined==condition)], 2, mean))
    bigSeurat@meta.data$chr7means[which(bigSeurat@meta.data$combined==condition)] <- as.numeric(apply(bigSeurat@assays$RNA@counts[which(chrAnno$chr%in%"chr7"),
                                                                                          which(bigSeurat@meta.data$combined==condition)], 2, mean))
    bigSeurat@meta.data$chr2means[which(bigSeurat@meta.data$combined==condition)] <- as.numeric(apply(bigSeurat@assays$RNA@counts[which(chrAnno$chr%in%"chr2"),
                                                                                          which(bigSeurat@meta.data$combined==condition)], 2, mean))
    
    temp <- as.numeric(apply(bigSeurat@assays$RNA@counts[which(chrAnno$chr%in%
            c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9",
                "chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17",
                 "chr18","chr19","chr20","chr21","chr22")),
                  which(bigSeurat@meta.data$combined==condition)], 2, mean))
    

  bigSeurat@meta.data$chrAllmeans[which(bigSeurat@meta.data$combined==condition)] <- temp

}

head(bigSeurat@meta.data)

options(repr.plot.width=110/6, repr.plot.height=40/6)
pdf("/data/leuven/343/vsc34313/jupyter_notebooks/XCI/bulkRNAseq/figs.bulkRNAseq4Amitesh/scRNAseq/forAmitesh.XISTKDscRNAseq.XAratio.pdf", width=8, height=4)
ggplot(bigSeurat@meta.data, aes(x=combined, y=chrXmeans/chrAllmeans, fill=conditions))+
geom_boxplot(width=0.3, size=0.5, outlier.size=0.5)+
geom_violin(alpha=0.5,scale="width", linewidth=0.5)+
theme_linedraw()+
theme(axis.text=element_text(size=6),
    axis.text.x=element_text(angle=40, hjust=1, vjust=1), legend.position="none", panel.grid.major.x  = element_blank(), panel.grid.minor=element_blank()
     )+
xlab("")+ylab("mean chrX to autosomal expression")+geom_vline(xintercept = c(4.5,8.5))
dev.off()

Idents(bigSeurat) <- "combined" 
so_sub<- subset(x = bigSeurat, idents = c("DAY-8_WT_DOX+", "DAY-8_WT_DOX-", "DAY-8_XIST_KD_DOX+", "DAY-8_XIST_KD_DOX-"))

Idents(so_sub) <- "anno.cluster"

so_sub_anno<- subset(x = so_sub, idents = c("TSC", "EXMC"))

mutate(so_sub_anno@meta.data, test=paste0(anno.cluster, "_", conditions))

options(repr.plot.width=10, repr.plot.height=8)
pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/XCI_analysis/forAmitesh.XISTKDscRNAseq.X2ratio_TSC_EXMC.pdf", width=10, height=8)
p<- ggplot(so_sub_anno@meta.data, aes(x=anno.cluster, y=chrXmeans/chrAllmeans, fill=conditions))+geom_hline(yintercept=1, color="red")+
geom_boxplot()+
theme_linedraw()+
theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), strip.background = element_rect(fill="light grey"), strip.text = element_text(color="black"))
dev.off()
+
facet_grid(conditions~source)
p

options(repr.plot.width=10, repr.plot.height=8)
pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/XCI_analysis/forAmitesh.XISTKDscRNAseq.X2ratio_TSC_EXMC_withviolin.pdf", width=10, height=8)
ggplot(dplyr::mutate(so_sub_anno@meta.data, test=paste0(anno.cluster, "_", conditions)), aes(x=test, y=chrXmeans/chrAllmeans, fill=test))+
geom_boxplot(width=0.3, size=0.5, linewidth=0.5, outlier.size=1)+geom_hline(yintercept=1, color="red")+
geom_violin(alpha=0.5,scale="width", linewidth=0.5)+
theme_linedraw()+
theme(axis.text=element_text(size=6),
    axis.text.x=element_text(angle=40, hjust=1, vjust=1), legend.position="none", panel.grid.major.x  = element_blank(), panel.grid.minor=element_blank())+
xlab("")+ylab("mean chrX to chr2 expression")+geom_vline(xintercept = c(4.5,8.5))+coord_cartesian(ylim = c(0,3))+
scale_fill_brewer(palette = "Dark2")
dev.off()

options(repr.plot.width=8, repr.plot.height=5)
pdf("/lustre1/project/stg_00041/Amitesh/XCI_PAPER_2023/XCI_analysis/XISTKDscRNAseq.X2Aratio_all_sample.pdf", width=8, height=4)
ggplot(bigSeurat@meta.data, aes(x=combined, y=chrXmeans/chrAllmeans, fill=conditions))+
geom_boxplot(width=0.3, size=0.5, linewidth=0.5, outlier.size=1)+
geom_violin(alpha=0.5,scale="width", linewidth=0.5)+
theme_linedraw()+
theme(axis.text=element_text(size=6),
    axis.text.x=element_text(angle=40, hjust=1, vjust=1), legend.position="none", panel.grid.major.x  = element_blank(), panel.grid.minor=element_blank())+
xlab("")+ylab("mean chrX to autosome expression")+geom_vline(xintercept = c(4.5,8.5))+coord_cartesian(ylim = c(0,3))+
scale_fill_brewer(palette = "Dark2")
dev.off()
