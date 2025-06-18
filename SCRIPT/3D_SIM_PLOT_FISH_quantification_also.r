library(Seurat)
library(dplyr)
library(tibble)
library(patchwork)
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(edgeR)
library(dplyr)
library(readxl)
library(ggplot2)



library(reshape2)
library(data.table)

library(openxlsx)

setwd ('/lustre1/project/stg_00041/Amitesh')

SIM_data <- read.xlsx("/lustre1/project/stg_00041/Amitesh/3D_SIM/Pooled_data.xlsx", sheet = 1)
head(SIM_data)

SIM_data <- SIM_data %>% rename(image="Original.Image.Name", cell_type="Cell.types")
head(SIM_data)

cells_per_image_per_celltype <- SIM_data %>%
    group_by(image, CellID, cell_type) %>%
    summarize(foci_out = case_when(NearestNeighborDistance_5 > 4 ~ 1,
                                   TRUE ~ 0)) %>%
    unique() %>%
    filter(foci_out == 1) %>%
    ungroup() %>%
    group_by(image, cell_type) %>%
    summarize(num_cells = length(foci_out))

total_cells_per_celltype <- SIM_data %>%
    select(image, CellID, cell_type) %>%
    unique() %>%
    group_by(cell_type) %>%
    summarize(total_cells=length(CellID))

cells_per_image_per_celltype <- cells_per_image_per_celltype %>% group_by(cell_type) %>% summarize(num_cells_total = sum(num_cells))

cells_per_image_per_celltype <- left_join(cells_per_image_per_celltype, total_cells_per_celltype)

head(cells_per_image_per_celltype)

p<- ggplot(cells_per_image_per_celltype, aes(x=cell_type, y=(num_cells_total/total_cells) * 100, fill=cell_type)) + ylim(0, 100) +
    geom_bar(stat="identity")
p

library(ggplot2)
library(dplyr)

# Compute the number of cells per image per cell type
cells_per_image_per_celltype <- SIM_data %>%
    group_by(image, CellID, cell_type) %>%
    summarize(foci_out = case_when(NearestNeighborDistance_5 > 5 ~ 1, TRUE ~ 0)) %>%
    unique() %>%
    filter(foci_out == 1) %>%
    ungroup() %>%
    group_by(image, cell_type) %>%
    summarize(num_cells = length(foci_out))

# Compute the total number of cells per cell type
total_cells_per_celltype <- SIM_data %>%
    select(image, CellID, cell_type) %>%
    unique() %>%
    group_by(cell_type) %>%
    summarize(total_cells = length(CellID))

# Aggregate data
cells_per_image_per_celltype <- cells_per_image_per_celltype %>%
    group_by(cell_type) %>%
    summarize(num_cells_total = sum(num_cells))

# Merge data
cells_per_image_per_celltype <- left_join(cells_per_image_per_celltype, total_cells_per_celltype)

# Compute percentage
cells_per_image_per_celltype <- cells_per_image_per_celltype %>%
    mutate(percent_cells = (num_cells_total / total_cells) * 100)

# Plot with box plot and overlaid violin plot
p <- ggplot(cells_per_image_per_celltype, aes(x = cell_type, y = percent_cells, fill = cell_type)) +
    geom_violin(alpha = 0.4) +  # Overlay violin plot with transparency
    geom_boxplot(width = 0.2, outlier.shape = NA) +  # Add boxplot, remove outliers to avoid clutter
    ylim(0, 100) +
    theme_minimal() +
    labs(y = "Percentage of Cells (%)", x = "Cell Type") +
    theme(legend.position = "none")

p


pdf("/lustre1/project/stg_00041/Amitesh/3D_SIM/output/NND_5_THRESHOLD_5.pdf", width = 10, height = 10)
plot(p)
dev.off()

test_data_h9naive_exmc <- c(cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESCs Naive"),]$num_cells_total,
               cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESCs Naive"),]$total_cells - cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESCs Naive"),]$num_cells_total,
               cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESC EXMCs"),]$num_cells_total,
               cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESC EXMCs"),]$total_cells - cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESC EXMCs"),]$num_cells_total)
test_data_h9naive_exmc

test_data_h9naive_tsc <- c(cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESCs Naive"),]$num_cells_total,
               cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESCs Naive"),]$total_cells - cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESCs Naive"),]$num_cells_total,
               cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESC TSCs"),]$num_cells_total,
               cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESC TSCs"),]$total_cells - cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESC TSCs"),]$num_cells_total)
test_data_h9naive_tsc

test_data_tsc_exmc <- c(cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESC TSCs"),]$num_cells_total,
               cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESC TSCs"),]$total_cells - cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESC TSCs"),]$num_cells_total,
               cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESC EXMCs"),]$num_cells_total,
               cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESC EXMCs"),]$total_cells - cells_per_image_per_celltype[which(cells_per_image_per_celltype$cell_type == "H9 ESC EXMCs"),]$num_cells_total)
test_data_tsc_exmc

enrichment <- function(data){
    out <- list()
    df <- as.table(rbind(c(data[[1]], data[[2]]), c(data[[3]], data[[4]])))
    dimnames(df) <- list(genes = c("naive", "not_naive"), test = c("Yes", "No"))
    out[[1]] <- fisher.test(df)
    out[[2]] <- chisq.test(df)
    out[[3]] <- as.data.table(df) %>% 
        dcast(genes~test) %>% 
        mutate(total=Yes+No, yes_perc=Yes/(Yes+No), no_perc=No/(Yes+No)) %>%
        select(genes, yes_perc, no_perc) %>%
        melt(id.vars = "genes", ) %>%
        select(-variable) %>%
       mutate(genes_factor=factor(.$genes, levels=c("X-linked genes", "Upregulated genes")))
    return(out)
}

enrichment_barplot <- function(df, title, lim){
    out <- ggplot(df, aes(x=genes_factor, 
                          fill=genes_factor,
                          y=value*100)) +
        geom_bar(stat="identity", color="black", width=1) +
        facet_wrap(~genes_factor, strip.position="bottom", scales="free_x") +
        scale_fill_manual(values=c("#77AADB", "#ee8866")) +
        ggtitle(title) +
        ylab(NULL) +
        xlab(NULL) +
        coord_cartesian(clip="off") +
        scale_x_discrete(labels=c("n=255", "n=72")) +
        scale_y_continuous(labels = function(x){paste0(x, "%")}, expand=c(0,0), limits = c(0, lim)) +
        theme(panel.grid=element_blank(),
              panel.background=element_blank(),
              axis.title.y=element_text(color="black", size=7),
              axis.text.y=element_text(color="black", size=7),
              axis.text.x=element_text(color="black", size=7, margin=margin(t=0, b=0)),
              legend.title=element_text(size=7, hjust=0.5),
              axis.ticks=element_line(color="black"),
              axis.ticks.x=element_blank(),
              axis.line.x=element_line(color="black"),
              axis.line.y=element_line(color="black"),
              legend.text=element_text(size=7),
              legend.position="none",
              plot.title=element_text(size=7, hjust=0.5),
              strip.background = element_blank(),
              strip.placement = "outside",
              strip.text.x = element_text(color="black", angle=90, size=7, hjust=1, vjust=0.5, margin=margin(t=0)),
              panel.spacing = unit(0, "lines"))    
    return(out)
}           

enrichment(test_data_h9naive_exmc)
enrichment(test_data_h9naive_tsc)
enrichment(test_data_tsc_exmc)

library(tidyverse)

cells_per_image_per_celltype <- cells_per_image_per_celltype %>%
    mutate(non_spreaded = total_cells - num_cells_total) %>%
    pivot_longer(cols = c(num_cells_total, non_spreaded), names_to = "XIST_cloud", values_to = "cell_count") %>%
    mutate(XIST_cloud = ifelse(XIST_cloud == "num_cells_total", "Spreaded", "Non-Spreaded"))

p1<- ggplot(cells_per_image_per_celltype, aes(x=cell_type, y=(cell_count/total_cells) * 100, fill=XIST_cloud)) +
    geom_bar(stat="identity") +
    ylim(0, 100) +
    ylab("Percentage of Cells") +
    xlab("Cell Type") +
    scale_fill_manual(values=c("Spreaded"="red", "Non-Spreaded"="gray")) +
    theme_minimal()
p1


pdf("/lustre1/project/stg_00041/Amitesh/3D_SIM/output/NND_5_THRESHOLD_5_withandwithoutXIST_spreaded_nuclei.pdf", width = 10, height = 10)
plot(p1)
dev.off()

SIM_data <- SIM_data %>% rename(number_of_foci_per_cell="Number.of.foci/nucleus")

unique_foci_data <- SIM_data %>%
  group_by(cell_type) %>%
  summarise(unique_foci = unique(`number_of_foci_per_cell`)) %>%
  ungroup()
unique_foci_data

# Create the box plot and violin plot overlay
p3<- ggplot(unique_foci_data, aes(x = cell_type, y = unique_foci, fill = cell_type)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot with some transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +  
geom_point(position = position_dodge(width = 0.9), shape = 21, color = "black", fill = "black", size = 2) + # Box plot with smaller width
  theme_minimal() +
  labs(title = "number of foci per cell",
       x = "Cell Type",
       y = "Unique Number of Foci per Nucleus") +
  theme(legend.position = "none")
p3


pdf("/lustre1/project/stg_00041/Amitesh/3D_SIM/output/number_of_foci_per_celltype.pdf", width = 10, height = 10)
plot(p3)
dev.off()

unique_foci_volume <- SIM_data %>%
  group_by(cell_type) %>%
  summarise(unique_volume = unique(`Foci.Volume`)) %>%
  ungroup()
unique_foci_volume

# Create the box plot and violin plot overlay
p4<- ggplot(unique_foci_volume, aes(x = cell_type, y = unique_volume, fill = cell_type)) +
  geom_violin(trim = FALSE, alpha = 0.5) +  # Violin plot with some transparency
  geom_boxplot(width = 0.2, color = "black", alpha = 0.8) +  
geom_point(position = position_dodge(width = 0.9), shape = 21, color = "black", fill = "black", size = 2) + # Box plot with smaller width
  theme_minimal() +
  labs(title = "volume of XIST",
       x = "Cell Type",
       y = "volume") +
  theme(legend.position = "none")
p4





ggplot(SIM_data, aes(x=cell_type, y=NearestNeighborDistance_5, color=cell_type)) +
    geom_jitter(alpha=0.4, width=0.2) +
    geom_hline(yintercept=5, linetype="dashed", color="red", size=1) +  # Add horizontal line at y=4
    theme_minimal() +
    ylab("Nearest Neighbor Distance (XIST Spread)")

table(SIM_data$NewCellID)

head(SIM_data)

# Load necessary libraries
library(ggplot2)
library(dplyr)

distance_cols <- c("NearestNeighborDistance_1", 
                   "NearestNeighborDistance_2", 
                   "NearestNeighborDistance_3", 
                   "NearestNeighborDistance_4", 
                   "NearestNeighborDistance_5")


# Define the cell types and distances for plotting
cell_types <- unique(SIM_data$cell_type)
plot_list <- list()  # To store plots for different distances

# Loop over NearestNeighborDistance columns to create histograms for each
for (dist in distance_cols) {
  p <- ggplot(SIM_data, aes_string(x = dist, fill = "cell_type")) +
    geom_histogram(binwidth = 0.1, alpha = 0.7, position = "identity") +
    facet_wrap(~ cell_type, nrow = 2) +
    labs(
      title = paste("Histogram of", dist),
      x = dist,
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    )
  
  # Add to the plot list
  plot_list[[dist]] <- p
  
  # Optionally display the plot during the loop
  print(p)
}

# Save the plots (optional)
#for (dist in names(plot_list)) {
 # ggsave(
   # filename = paste0("Histogram_", dist, ".png"),
 ##   plot = plot_list[[dist]],
   ## width = 8,
   # height = 6,
   # dpi = 300
  #)
#}



SIM_data_coordinate <- read.xlsx("/lustre1/project/stg_00041/Amitesh/3D_SIM/Vesicles_Type_A_Position.xlsx", sheet = 1)
head(SIM_data_coordinate)

library(ggplot2)
library(dplyr)

# Step 1: Count the number of foci per cell where NearestNeighborDistance_5 > 4
foci_per_cell <- SIM_data %>%
  filter(NearestNeighborDistance_5 > 4) %>%
  group_by(image, CellID, cell_type) %>%
  summarise(num_foci = n(), .groups = "drop")  # Count foci per cell

# Step 2: Create a bin plot (histogram) for each condition
p <- ggplot(foci_per_cell, aes(x = num_foci)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black") +  # Adjust bin width if needed
  facet_wrap(~cell_type, scales = "free_y") +  # Separate plots for each condition
  labs(x = "Number of Foci with NearestNeighborDistance_5 > 4",
       y = "Number of Cells",
       title = "Distribution of Foci Per Cell Across Conditions") +
  theme_minimal()

# Display plot
print(p)


library(ggplot2)
library(dplyr)

# Step 1: Create the bin plot for all foci and all cells within each condition
p <- ggplot(SIM_data, aes(x = NearestNeighborDistance_5)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "black", alpha = 0.7) +  # Adjust bin width if needed
  facet_wrap(~cell_type, scales = "free_y") +  # Separate plots for each condition
  geom_vline(xintercept = 4, linetype = "dashed", color = "red", size = 1) +  # Mark threshold
  labs(x = "Nearest Neighbor Distance 5",
       y = "Count of Foci",
       title = "Distribution of Nearest Neighbor Distance 5 Across Conditions") +
  theme_minimal()

# Display plot
print(p)

library(ggplot2)
library(dplyr)

# Step 1: Create the bin plot for all foci and all cells within each condition
p <- ggplot(SIM_data, aes(x = NearestNeighborDistance_5)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "black", alpha = 0.7) +  # Adjust bin width if needed
  facet_wrap(~cell_type, scales = "free_y") +  # Separate plots for each condition
  geom_vline(xintercept = 4, linetype = "dashed", color = "red", size = 1) +  # Mark threshold
  labs(x = "Nearest Neighbor Distance 5",
       y = "Count of Foci",
       title = "Distribution of Nearest Neighbor Distance 5 Across Conditions") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 15)) +  # Limit y-axis to 0-100 range
  xlim(5, NA)  # Focus on values beyond 5 for the x-axis

# Display plot
print(p)

library(ggplot2)
library(dplyr)

# Step 1: Aggregate the data by cell and count the number of cells for each foci count per cell
SIM_data_per_cell <- SIM_data %>%
  group_by(cell_type, number_of_foci_per_cell) %>%
  summarise(cell_count = n())  # Count number of cells for each number of foci per cell

# Step 2: Create the bin plot for the number of cells per foci count within each condition
p <- ggplot(SIM_data_per_cell, aes(x = number_of_foci_per_cell, y = cell_count)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", alpha = 0.7) +  # Bar plot for cell count
  facet_wrap(~cell_type, scales = "free_y") +  # Separate plots for each condition
  geom_vline(xintercept = 4, linetype = "dashed", color = "red", size = 1) +  # Mark threshold
  labs(x = "Number of Foci per Cell",
       y = "Count of Cells",
       title = "Distribution of Number of Foci per Cell Across Conditions") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 100)) +  # Limit y-axis to 0-100 range
  xlim(5, NA)  # Focus on values beyond 5 for the x-axis

# Display plot
print(p)

library(ggplot2)
library(dplyr)
library(tidyr)

# Binning the NearestNeighborDistance_5 into defined bins
SIM_data <- SIM_data %>%
  mutate(distance_bin = cut(NearestNeighborDistance_5, 
                            breaks = seq(0, max(NearestNeighborDistance_5, na.rm = TRUE), by = 1), 
                            include.lowest = TRUE, 
                            labels = paste(seq(0, max(NearestNeighborDistance_5, na.rm = TRUE) - 1), 
                                           seq(1, max(NearestNeighborDistance_5, na.rm = TRUE)))))

# Initialize list to store plots
p_list <- list()

# Get unique cell types (which we treat as conditions)
cell_types <- unique(SIM_data$cell_type)

# Loop over each cell_type and create the plot
for (cell_type in cell_types) {
  p <- SIM_data %>%
    filter(cell_type == !!cell_type) %>%
    ggplot(aes(x = distance_bin, fill = as.factor(CellID))) +
    geom_bar(position = "dodge", stat = "count", show.legend = FALSE) +
    labs(title = paste("Foci in Cell Type:", cell_type), 
         x = "Nearest Neighbor Distance (binned)", 
         y = "Number of Foci") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Set3") +
    facet_wrap(~ cell_type)  # To create a facet for each cell type
  p_list[[cell_type]] <- p
}

# Displaying the first plot (change index to view others)
print(p_list[[3]])

# You can save the plot if needed:
# ggsave("cell_type_plot.png", plot = p_list[[1]], width = 8, height = 6)



pdf("/lustre1/project/stg_00041/Amitesh/3D_SIM/output/bin_plot_naive.pdf", width = 10, height = 10)
plot(p_list[[3]])
dev.off()

library(ggplot2)
library(dplyr)

# Create bins for NearestNeighborDistance_5 (e.g., increments of 1)
SIM_data <- SIM_data %>%
  mutate(distance_bin = cut(NearestNeighborDistance_5, 
                            breaks = seq(0, max(NearestNeighborDistance_5, na.rm = TRUE), by = 1), 
                            include.lowest = TRUE, 
                            labels = paste(seq(0, max(NearestNeighborDistance_5, na.rm = TRUE) - 1), 
                                           seq(1, max(NearestNeighborDistance_5, na.rm = TRUE)))))

# Count the number of cells in each bin for each cell_type
cells_per_bin_per_celltype <- SIM_data %>%
  group_by(cell_type, distance_bin) %>%
  summarize(num_cells = n_distinct(CellID)) %>%
  ungroup()

# Create the plot
p <- ggplot(cells_per_bin_per_celltype, aes(x = distance_bin, y = num_cells, fill = cell_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of Cells per Distance Bin per Cell Type", 
       x = "Nearest Neighbor Distance (binned)", 
       y = "Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~ cell_type, scales = "free_y")  # Facet by cell_type

# Display the plot
print(p)


library(ggplot2)
library(dplyr)

# Create bins for NearestNeighborDistance_5 with step size 2
SIM_data <- SIM_data %>%
  mutate(distance_bin = cut(NearestNeighborDistance_5, 
                            breaks = seq(0, max(NearestNeighborDistance_5, na.rm = TRUE), by = 2), 
                            include.lowest = TRUE, 
                            right = FALSE))  # Right=FALSE makes intervals inclusive of the left edge

# Count the number of cells in each bin for each cell_type
cells_per_bin_per_celltype <- SIM_data %>%
  group_by(cell_type, distance_bin) %>%
  summarize(num_cells = n_distinct(CellID)) %>%
  ungroup()

# Create the plot
p <- ggplot(cells_per_bin_per_celltype, aes(x = distance_bin, y = num_cells, fill = cell_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of Cells per Distance Bin per Cell Type", 
       x = "Nearest Neighbor Distance (binned)", 
       y = "Number of Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~ cell_type, scales = "free_y")  # Facet by cell_type

# Display the plot
print(p)

library(ggplot2)
library(dplyr)

# Extract NearestNeighborDistance_5 values per cell type
distance_distribution <- SIM_data %>%
    select(image, CellID, cell_type, NearestNeighborDistance_5) %>%
    unique()

# Plot histogram
p <- ggplot(distance_distribution, aes(x=NearestNeighborDistance_5, fill=cell_type)) +
    geom_histogram(binwidth=1, position="dodge", alpha=0.7) +
    labs(x="Nearest Neighbor Distance", y="Cell Count", title="Distribution of Nearest Neighbor Distance") +
    theme_minimal()
p


# Bin distances into ranges
distance_distribution_binned <- SIM_data %>%
    mutate(distance_bin = cut(NearestNeighborDistance_5, breaks=seq(0, max(NearestNeighborDistance_5, na.rm=TRUE), by=2))) %>%
    group_by(cell_type, distance_bin) %>%
    summarize(count = n(), .groups = 'drop')

# Plot binned distribution
p <- ggplot(distance_distribution_binned, aes(x=distance_bin, y=count, fill=cell_type)) +
    geom_bar(stat="identity", position="dodge") +
    labs(x="Distance Bin", y="Cell Count", title="Binned Nearest Neighbor Distance Distribution") +
    theme_minimal()
p


library(ggplot2)
library(dplyr)

# Remove duplicates
distance_distribution <- SIM_data %>%
    select(image, CellID, cell_type, NearestNeighborDistance_5) %>%
    distinct()

# Plot histogram with proper bin size
p <- ggplot(distance_distribution, aes(x=NearestNeighborDistance_5, fill=cell_type)) +
    geom_histogram(binwidth=2, position="identity", alpha=0.6, color="black") +
    labs(x="Nearest Neighbor Distance", y="Cell Count", title="Distribution of Nearest Neighbor Distance") +
    theme_minimal()
p






# Load the required libraries
library(ggplot2)
library(reshape2)

# Data (replace with your actual data)
D0 <- c(9.237011069, 6.674485418, 15.68937481, 8.012633062, 6.880354286, 5.780581121, 12.88847468, 37.48438685, 9.556649575, 12.99140912, 12.27086808, 12.05958161, 8.413535595, 11.75619591, 51.56473393, 15.09343861, 13.46815807, 14.4324912, 17.87808594, 22.84060919, 31.87175139, 19.48711367, 7.384191253, 14.31872156, 17.7534811, 8.012633062, 13.8528078, 16.41533345, 16.63203753, 24.35753769, 13.00766192, 12.36296731, 19.5629601, 19.14038716, 15.44558273, 7.178322385, 12.69885862, 13.16477237, 12.00540559, 12.23836247, 10.43430107, 12.48757215, 12.6501002, 17.39591938, 18.81533105)
D4 <- c(5.08712809, 3.234308275, 6.625727001, 6.896607091, 4.75665438, 10.67809315, 2.676295289, 2.291645561, 3.185549858, 5.726405103, 1.376070857, 5.655976279, 4.729566371, 13.86906061, 3.147626646, 9.350780707, 10.42346586, 4.06320135, 5.325502569, 3.174714655, 7.286674421, 3.174714655, 3.889838092, 3.169297053, 4.052366146, 11.65326147, 6.446946142, 3.174714655, 16.41533345, 14.82797613, 10.17967378, 7.216245598, 2.042435879, 8.245589939, 14.17786391, 4.111959766, 5.921438767, 4.854171213, 9.502473557, 5.666811483, 4.317828635, 7.091640756, 9.182835051, 7.302927226, 8.353941975, 4.464103883, 7.925951433, 5.504283429, 2.459591217, 1.262301219, 5.910603564, 5.108798497, 13.34897084, 9.128659033, 8.483964419)
D8 <- c(9.675836815, 7.248751208, 7.91511623, 8.711503694, 9.15032944, 1.072685156, 0.93182751, 5.910603564, 5.921438767, 10.60766432, 5.856427546, 15.16928504, 20.4568644, 13.61985093, 8.678998083, 4.458686281, 2.069523888, 4.361169449, 7.210827996, 3.337242709, 2.351239181, 6.641979807, 7.113311163, 5.292996959, 3.46184755, 5.645141076, 3.543111577, 5.406766596, 7.438367271, 4.54536791, 4.745819177, 1.988259861, 1.961171852, 2.822570538, 3.618958002, 4.561620716, 8.695250889, 6.148978043, 2.535437642, 3.981937323, 6.300670893, 2.085776693, 3.396836329, 5.677646686, 5.125051303, 7.888028221, 9.870870479, 3.895255694, 3.212637867, 4.231147006, 5.363425782, 5.921438767, 6.029790803)
max_length <- max(length(D0), length(D4), length(D8))

D0 <- c(D0, rep(NA, max_length - length(D0)))
D4 <- c(D4, rep(NA, max_length - length(D4)))
D8 <- c(D8, rep(NA, max_length - length(D8)))

# Create a data frame
df <- data.frame(
  D0 = D0,
  D4 = D4,
    D8 = D8
)

df

df_melt <- melt(df, variable.name = "Condition", value.name = "XIST_VOULME")

df_melt

# Create the violin plot with the box plot overlay and individual data points
p2 <- ggplot(df_melt, aes(x = Condition, y = XIST_VOULME, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5) +   # Violin plot (trim = FALSE shows full distribution)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +  # Overlay box plot
  geom_point(position = position_dodge(width = 0.9), shape = 21, color = "black", fill = "black", size = 2) +  # Add individual data points
  labs(title = "POLA1 Max Feret Diameter Across Conditions (D0, D4, D8)", 
       x = "Condition", y = "XIST_VOULME") +
  theme_minimal() +
  scale_fill_manual(values = c("D0" = "blue", "D4" = "green", "D8" = "red"))

# Print the plot
p2

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "XIST_volume_dynamics_timecourse.pdf"

ggsave(file.path(output_path, filename), plot = p2, height = 8, width = 8)

library(ggplot2)
library(ggpubr)

# Define pairwise comparisons
my_comparisons <- list(c("D0", "D4"), c("D0", "D8"), c("D4", "D8"))

# Create the violin + box + point + stat plot
p2 <- ggplot(df_melt, aes(x = Condition, y = XIST_VOULME, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.9), shape = 21, color = "black", fill = "black", size = 2) +
  labs(title = "POLA1 Volume Across Conditions (D0, D4, D8)", 
       x = "Condition", y = "XIST_VOULME") +
  theme_minimal() +
  scale_fill_manual(values = c("D0" = "blue", "D4" = "green", "D8" = "red")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format")
p2

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "XIST_volume_dynamics_timecourse_with_wilcox_stat.pdf"

ggsave(file.path(output_path, filename), plot = p2, height = 8, width = 8)

# Load the required libraries
library(ggplot2)
library(reshape2)

# Data (replace with your actual data)
D0 <- c(4.590891076, 3.781442063, 7.070439745, 4.320861733, 4.503720811, 3.979324347, 7.037331966, 17.1028569, 5.889675417, 5.385685415, 5.07822235, 4.970244825, 5.518926312, 4.50210967, 21.37216415, 6.204015188, 5.117665188, 5.580997751, 9.204172016, 7.106665452, 11.46862358, 8.540937283, 4.627193532, 7.861871302, 5.824761668, 7.372796121, 5.030951314, 5.919995609, 5.820321815, 9.310975989, 5.143840499, 5.683707569, 5.539095075, 9.317246737, 4.641402314, 4.490173604, 7.462558665, 5.198874624, 5.79744163, 5.954271059, 4.026374954, 4.231073981, 4.673970597, 6.013767404, 6.391860318)
D4 <- c(4.059570424, 3.753100675, 3.559236278, 4.378614503, 3.222902051, 3.84891975, 2.645190032, 3.254348431, 2.859203226, 3.296063524, 2.383834248, 3.884039737, 3.802559819, 4.946836292, 3.279854371, 3.95792789, 5.065349161, 2.790156942, 3.107511275, 3.125269307, 5.384355345, 3.437636471, 2.68320596, 3.781737453, 2.715319039, 4.38504547, 3.231806272, 2.514874594, 5.246139039, 6.817184031, 4.674070146, 3.533009276, 2.095248751, 3.489728407, 7.701766296, 4.562257497, 3.617733456, 3.053619856, 3.905816455, 3.233849308, 3.164631005, 3.972703811, 4.3777643, 3.503108563, 4.066917785, 3.278975944, 3.739595351, 2.982245744, 2.888361934, 1.808236784, 3.771491132, 12.49843581, 3.986399511, 5.425601751, 6.602334012, 4.035303664)
D8 <- c(5.652530667, 3.394153996, 3.695764557, 4.439597467, 4.149325954, 1.660554512, 1.869734779, 3.456254322, 3.5583475, 14.67460854, 4.586493212, 6.155387438, 6.911607545, 5.352617023, 4.878834885, 4.052989801, 2.769043481, 2.095248751, 3.065384387, 4.25273614, 2.655191723, 2.424799566, 3.434766985, 7.404215547, 3.901146624, 4.91355923, 2.790156942, 3.473138086, 3.16398444, 6.928716526, 6.72711895, 4.78653528, 3.338677823, 3.158070035, 3.497183147, 3.793597377, 4.256933104, 3.067568379, 3.851649465, 3.453265569, 3.337786689, 5.468476145, 3.231403966, 2.544983405, 4.348047585, 3.189105448, 3.030599001, 6.264109037, 4.021728991, 3.180257027, 4.773596494, 3.758944422, 3.29603498, 3.642431638, 2.893606311)
max_length <- max(length(D0), length(D4), length(D8))

D0 <- c(D0, rep(NA, max_length - length(D0)))
D4 <- c(D4, rep(NA, max_length - length(D4)))
D8 <- c(D8, rep(NA, max_length - length(D8)))

# Create a data frame
df_1 <- data.frame(
  D0 = D0,
  D4 = D4,
    D8 = D8
)

df_1

df_melt_1 <- melt(df_1, variable.name = "Condition", value.name = "XIST_max_feret_diameter")

df_melt_1


# Create the violin plot with the box plot overlay and individual data points
p3 <- ggplot(df_melt_1, aes(x = Condition, y = XIST_max_feret_diameter, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5) +   # Violin plot (trim = FALSE shows full distribution)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +  # Overlay box plot
  geom_point(position = position_dodge(width = 0.9), shape = 21, color = "black", fill = "black", size = 2) +  # Add individual data points
  labs(title = "POLA1 Max Feret Diameter Across Conditions (D0, D4, D8)", 
       x = "Condition", y = "XIST_max_feret_diameter") +
  theme_minimal() +
  scale_fill_manual(values = c("D0" = "blue", "D4" = "green", "D8" = "red"))


# Print the plot
p3

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "XIST_max_feret_diameter_dynamics_timecourse.pdf"

ggsave(file.path(output_path, filename), plot = p3, height = 8, width = 8)

library(ggpubr)

# Define comparisons
my_comparisons <- list(c("D0", "D4"), c("D0", "D8"), c("D4", "D8"))

# Create plot
p3 <- ggplot(df_melt_1, aes(x = Condition, y = XIST_max_feret_diameter, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.9), shape = 21, color = "black", fill = "black", size = 2) +
  labs(title = "POLA1 Max Feret Diameter Across Conditions (D0, D4, D8)", 
       x = "Condition", y = "XIST_max_feret_diameter") +
  theme_minimal() +
  scale_fill_manual(values = c("D0" = "blue", "D4" = "green", "D8" = "red")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format")
p3

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "XIST_max_feret_diameter_dynamics_timecourse_with_wilcox_test.pdf"

ggsave(file.path(output_path, filename), plot = p3, height = 8, width = 8)

# Load required libraries
library(ggplot2)
library(dplyr)

# Combine all your data into vectors
POLA1_Volume <- c(
  # Day 0
  0.55644, 0.907, 1.025, 0.14643, 0.64, 0.322, 0.34168, 0.32215, 0.12, 0.312, 0.16, 1.2, 0.31, 0.41001, 0.2, 0.156,
  # Day 4
  0.19524, 0.62478, 0.556, 0.33191, 0.70288, 0.1757, 0.68335, 0.8786, 0.73216, 0.8786, 0.45882, 0.56621, 0.54668,
  0.21477, 0.41001, 0.45882, 0.83955, 0.62478, 0.94693, 0.31239, 0.56621
)

Condition <- rep("XIST_positive", times = length(POLA1_Volume))

Timepoint <- c(rep("Day 0", 16), rep("Day 4", 21))

# Repeat for XIST-negative
POLA1_Volume_neg <- c(
  # Day 0
  1.0348, 1.1227, 1.2398, 0.8395, 0.663, 0.429, 1.425, 1.0348, 0.4393, 0.751, 0.4, 1.8, 0.95669, 1.4643, 1.161, 0.741,
  # Day 4
  0.68335, 1.3374, 0.97622, 1.454, 1.4448, 1.054, 2.6163, 1.2593, 0.75169, 1.2593, 0.56621, 1.1129, 0.55644,
  0.40025, 0.5174, 0.77121, 1.5131, 1.5717, 1.1129, 0.47835, 0.90788
)

Condition_neg <- rep("XIST_negative", times = length(POLA1_Volume_neg))
Timepoint_neg <- c(rep("Day 0", 16), rep("Day 4", 21))

# Combine both into one data frame
df <- data.frame(
  POLA1_Volume = c(POLA1_Volume, POLA1_Volume_neg),
  Condition = c(Condition, Condition_neg),
  Timepoint = c(Timepoint, Timepoint_neg)
)

df

library(ggplot2)
library(ggpubr)

p <- ggplot(df, aes(x = Condition, y = POLA1_Volume, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.04, size = 1.5, alpha = 0.6) +
  facet_wrap(~Timepoint) +
  labs(title = "POLA1 Volume Comparison by XIST Status and Timepoint",
       y = "POLA1 Volume", x = "Condition") +
  scale_fill_manual(values = c("XIST_positive" = "red", "XIST_negative" = "blue"))

p1 <- p + stat_compare_means(method = "wilcox.test", label = "p.format",
                       comparisons = list(c("XIST_negative", "XIST_positive")),
                       label.y = c(3.2, 3.2))

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "POLA1_Volumeax_feret_diameter_dynamics_timecourse.pdf"

ggsave(file.path(output_path, filename), plot = p1, height = 8, width = 8)

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/ggpubr")



# Load necessary libraries
library(ggplot2)
library(dplyr)

# Example data
data <- data.frame(
  CellType = rep(c("Naive", "TSC", "EXMC"), each = 3),
  Experiment = rep(1:3, times = 3),
  Monoallelic = c(16.34, 9.43, 6.06, 87.72, 85.23, 87.23, 97.78, 96.36, 97.73),
  Biallelic = c(83.66, 71.70, 93.94, 12.28, 13.64, 13.83, 2.22, 5.45, 2.27)
)

# Reshape data for plotting
data_long <- data %>%
  pivot_longer(cols = c("Monoallelic", "Biallelic"), names_to = "ExpressionType", values_to = "Proportion")

# Plot
ggplot(data_long, aes(x = CellType, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ Experiment) +  # Facet by experiment
  labs(title = "Proportion of Monoallelic vs. Biallelic Expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("Monoallelic" = "skyblue", "Biallelic" = "orange")) +
  theme_minimal()

# Install tidyr if you haven't already
install.packages("tidyr")

# Load the tidyr package
library(tidyr)

# Example data
data <- data.frame(
  CellType = rep(c("Naive", "TSC", "EXMC"), each = 3),
  Experiment = rep(1:3, times = 3),
  Monoallelic = c(97.91666667, 100, 100, 100, 100, 100, 99.43181818, 100, 100),
  Biallelic = c(2.083333333, 0, 0, 0, 0, 0, 0.5681818182, 0, 0)
)
# Reshape data to long format
library(tidyr)
data_long <- data %>%
  pivot_longer(cols = c("Monoallelic", "Biallelic"), names_to = "ExpressionType", values_to = "Proportion")

# Load necessary libraries
library(ggplot2)
# Custom color palette
bi_color <- "#94a096"
mono_color <- "#ce4751"
# Plot with stacked bars and dots for individual experiments
p<- ggplot(data_long, aes(x = CellType, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +
  geom_point(aes(x = CellType, y = Proportion, color = ExpressionType),
             position = position_jitter(width = 0.2), size = 2, shape = 21, stroke = 1.5) +
  labs(title = "Proportion of Monoallelic vs. Biallelic Expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("Monoallelic" =  mono_color, "Biallelic" = bi_color)) +
  scale_color_manual(values = c("Monoallelic" =  mono_color, "Biallelic" = bi_color)) +
  theme_minimal() +
  theme(legend.position = "top")
p

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "XIST_with_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p, height = 8, width = 8)

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Create the new data with percentages
data <- data.frame(
  Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP3', 'EXP3', 'EXP3'),
  Cell_type = c('NAIVE', 'TSC', 'EXMC', 'NAIVE', 'TSC', 'EXMC', 'NAIVE', 'TSC', 'EXMC'),
  Percent_Mono = c(97.91666667, 100, 100, 100, 100, 100, 99.43181818, 100, 100),
  Percent_Bi = c(2.083333333, 0, 0, 0, 0, 0, 0.5681818182, 0, 0)
)

# Calculate mean and standard deviation for each cell type across experiments
summary_data <- data %>%
  group_by(Cell_type) %>%
  summarise(
    Mono_mean = mean(Percent_Mono),
    Mono_sd = sd(Percent_Mono),
    Bi_mean = mean(Percent_Bi),
    Bi_sd = sd(Percent_Bi)
  )

# Custom color palette
mono_color <- "#ce4751"
bi_color <- "#94a096"

# Create a stacked bar plot
p2<- ggplot(summary_data, aes(x = Cell_type)) +
  geom_bar(aes(y = Mono_mean, fill = 'Monoallelic'), stat = "identity", width = 0.7) +
  geom_bar(aes(y = Bi_mean, fill = 'Biallelic'), stat = "identity", width = 0.7, position = 'stack') +
  geom_errorbar(aes(ymin = Mono_mean - Mono_sd, ymax = Mono_mean + Mono_sd), width = 0.2) +
  geom_errorbar(aes(ymin = Mono_mean + Bi_mean - Bi_sd, ymax = Mono_mean + Bi_mean + Bi_sd), width = 0.2) +
  scale_fill_manual(values = c('Monoallelic' = mono_color, 'Biallelic' = bi_color)) +
  labs(title = "Mean Percentage of Monoallelic and Biallelic Expression Across Cell Types",
       x = "Cell Type",
       y = "Mean Percentage") +
  theme_minimal()
p2

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "XIST_with_individual_stdev_error_bar.pdf"

ggsave(file.path(output_path, filename), plot = p2, height = 8, width = 8)

# Example data
data_THOC2 <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP3', 'EXP3', 'EXP3'),
    Cell_type = c('NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC'),
  Monoallelic = c(21.590909, 80.597015, 79.130435, 27.737226, 57.5, 73.043478, 19.91342, 69.791667, 74.452555),
  Biallelic = c(78.40909091, 19.40298507, 20.86956522, 72.26277372, 42.5, 26.95652174, 80.08658009, 30.20833333, 25.54744526)
)

# Reshape data to long format
library(tidyr)
data_long <- data_THOC2 %>%
  pivot_longer(cols = c("Biallelic", "Monoallelic"), names_to = "ExpressionType", values_to = "Proportion")

print(data_long)

# Load necessary libraries
library(ggplot2)
# Custom color palette
bi_color <- "#94a096"
mono_color <- "#ce4751"
# Plot with stacked bars and dots for individual experiments
p_THOC2<- ggplot(data_long, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +
  geom_point(aes(x = Cell_type, y = Proportion, color = ExpressionType),
             position = position_jitter(width = 0.2), size = 2, shape = 21, stroke = 1.5) +
  labs(title = "Proportion of Monoallelic vs. Biallelic Expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("Monoallelic" =  mono_color, "Biallelic" = bi_color)) +
  scale_color_manual(values = c("Monoallelic" =  mono_color, "Biallelic" = bi_color)) +
  theme_minimal() +
  theme(legend.position = "top")
p_THOC2

# Load necessary libraries
library(ggplot2)

# Custom color palette for bars
bi_color <- "#94a096"
mono_color <- "#ce4751"

# New color palette for dots (distinct but corresponding)
bi_dot_color <- "#5b6d7c"  # Darker shade of bi_color
mono_dot_color <- "#b03038"  # Darker shade of mono_color

# Plot with stacked bars and bold dots for individual experiments
p_THOC2 <- ggplot(data_long, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +  # Bars remain the same
  geom_point(aes(x = Cell_type, y = Proportion, color = ExpressionType),  # Different colors for dots
             position = position_jitter(width = 0.2, height = 0), 
             size = 3.5, shape = 21, stroke = 1.8) +  # Bold dots with larger size and stroke
  labs(title = "Proportion of Monoallelic vs. Biallelic Expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("Monoallelic" = mono_color, "Biallelic" = bi_color)) +  # Keep bar colors
  scale_color_manual(values = c("Monoallelic" = mono_dot_color, "Biallelic" = bi_dot_color)) +  # Different dot colors
  theme_minimal() +
  theme(legend.position = "top")

p_THOC2

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "THOC2_with_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p_THOC2, height = 8, width = 8)

# Example data
data_HUWE1 <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP3', 'EXP3', 'EXP3'),
    Cell_type = c('NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC'),
  Monoallelic = c(19.13580247, 90.26548673, 64.56692913, 15.46391753, 87.35632184, 66.12903226, 14.87603306, 77.5862069, 67.32673267),
  Biallelic = c(80.86419753, 9.734513274, 35.43307087, 84.53608247, 12.64367816, 33.87096774, 85.12396694, 22.4137931, 32.67326733)
)

# Reshape data to long format
library(tidyr)
data_long_HUWE1 <- data_HUWE1 %>%
  pivot_longer(cols = c("Biallelic", "Monoallelic"), names_to = "ExpressionType", values_to = "Proportion")

# Load necessary libraries
library(ggplot2)

# Custom color palette for bars
bi_color <- "#94a096"
mono_color <- "#ce4751"

# New color palette for dots (distinct but corresponding)
bi_dot_color <- "#5b6d7c"  # Darker shade of bi_color
mono_dot_color <- "#b03038"  # Darker shade of mono_color

# Plot with stacked bars and bold dots for individual experiments
p_HUWE1 <- ggplot(data_long_HUWE1, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +  # Bars remain the same
  geom_point(aes(x = Cell_type, y = Proportion, color = ExpressionType),  # Different colors for dots
             position = position_jitter(width = 0.2, height = 0), 
             size = 3.5, shape = 21, stroke = 1.8) +  # Bold dots with larger size and stroke
  labs(title = "Proportion of Monoallelic vs. Biallelic Expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("Monoallelic" = mono_color, "Biallelic" = bi_color)) +  # Keep bar colors
  scale_color_manual(values = c("Monoallelic" = mono_dot_color, "Biallelic" = bi_dot_color)) +  # Different dot colors
  theme_minimal() +
  theme(legend.position = "top")

p_HUWE1

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "HUWE1_with_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p_HUWE1, height = 8, width = 8)



# Example data
data_POLA1 <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP3', 'EXP3', 'EXP3'),
    Cell_type = c('NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC'),
  Monoallelic = c(15.43209877, 96.49122807, 86.75496689, 31.2, 94.64285714, 86.20689655, 12.08791209, 97.72727273, 86.31578947),
  Biallelic = c(84.56790123, 3.50877193, 13.24503311, 68.8, 5.357142857, 13.79310345, 87.91208791, 2.272727273, 13.68421053)
)

# Reshape data to long format
library(tidyr)
data_long_POLA1 <- data_POLA1 %>%
  pivot_longer(cols = c("Biallelic", "Monoallelic"), names_to = "ExpressionType", values_to = "Proportion")

# Load necessary libraries
library(ggplot2)

# Custom color palette for bars
bi_color <- "#94a096"
mono_color <- "#ce4751"

# New color palette for dots (distinct but corresponding)
bi_dot_color <- "#5b6d7c"  # Darker shade of bi_color
mono_dot_color <- "#b03038"  # Darker shade of mono_color

# Plot with stacked bars and bold dots for individual experiments
p_POLA1 <- ggplot(data_long_POLA1, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +  # Bars remain the same
  geom_point(aes(x = Cell_type, y = Proportion, color = ExpressionType),  # Different colors for dots
             position = position_jitter(width = 0.2, height = 0), 
             size = 3.5, shape = 21, stroke = 1.8) +  # Bold dots with larger size and stroke
  labs(title = "Proportion of Monoallelic vs. Biallelic Expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("Monoallelic" = mono_color, "Biallelic" = bi_color)) +  # Keep bar colors
  scale_color_manual(values = c("Monoallelic" = mono_dot_color, "Biallelic" = bi_dot_color)) +  # Different dot colors
  theme_minimal() +
  theme(legend.position = "top")

p_POLA1

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "POLA1_with_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p_POLA1, height = 8, width = 8)



# Example data
data_XIST_POLA1 <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP3', 'EXP3', 'EXP3'),
    Cell_type = c('NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC'),
  Monoallelic = c(10.22727273, 89.10891089, 87.37864078, 20.45454545, 98.82352941, 97.22222222, 12.30769231, 95.77464789, 90.51724138),
  Biallelic = c(89.77272727, 10.89108911, 12.62135922, 79.54545455, 1.176470588, 2.777777778, 87.69230769, 4.225352113, 9.482758621)
)

# Reshape data to long format
library(tidyr)
data_long_XIST_POLA1 <- data_XIST_POLA1 %>%
  pivot_longer(cols = c("Biallelic", "Monoallelic"), names_to = "ExpressionType", values_to = "Proportion")

# Load necessary libraries
library(ggplot2)

# Custom color palette for bars
bi_color <- "#94a096"
mono_color <- "#ce4751"

# New color palette for dots (distinct but corresponding)
bi_dot_color <- "#5b6d7c"  # Darker shade of bi_color
mono_dot_color <- "#b03038"  # Darker shade of mono_color

# Plot with stacked bars and bold dots for individual experiments
p_XIST_POLA1 <- ggplot(data_long_XIST_POLA1, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +  # Bars remain the same
  geom_point(aes(x = Cell_type, y = Proportion, color = ExpressionType),  # Different colors for dots
             position = position_jitter(width = 0.2, height = 0), 
             size = 3.5, shape = 21, stroke = 1.8) +  # Bold dots with larger size and stroke
  labs(title = "Proportion of Monoallelic vs. Biallelic Expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("Monoallelic" = mono_color, "Biallelic" = bi_color)) +  # Keep bar colors
  scale_color_manual(values = c("Monoallelic" = mono_dot_color, "Biallelic" = bi_dot_color)) +  # Different dot colors
  theme_minimal() +
  theme(legend.position = "top")

p_XIST_POLA1

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "XIST_POLA1_with_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p_XIST_POLA1, height = 8, width = 8)



# Example data
data_XIST <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP3', 'EXP3', 'EXP3'),
    Cell_type = c('NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC'),
  Monoallelic = c(97.91666667, 100, 100, 100, 100, 100, 99.43181818, 100, 100),
  Biallelic = c(2.083333333, 0, 0, 0, 0, 0, 0.5681818182, 0, 0)
)

# Reshape data to long format
library(tidyr)
data_long_XIST <- data_XIST %>%
  pivot_longer(cols = c("Biallelic", "Monoallelic"), names_to = "ExpressionType", values_to = "Proportion")

# Load necessary libraries
library(ggplot2)

# Custom color palette for bars
bi_color <- "#94a096"
mono_color <- "#ce4751"

# New color palette for dots (distinct but corresponding)
bi_dot_color <- "#5b6d7c"  # Darker shade of bi_color
mono_dot_color <- "#b03038"  # Darker shade of mono_color

# Plot with stacked bars and bold dots for individual experiments
p_XIST <- ggplot(data_long_XIST, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +  # Bars remain the same
  geom_point(aes(x = Cell_type, y = Proportion, color = ExpressionType),  # Different colors for dots
             position = position_jitter(width = 0.2, height = 0), 
             size = 3.5, shape = 21, stroke = 1.8) +  # Bold dots with larger size and stroke
  labs(title = "Proportion of Monoallelic vs. Biallelic Expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("Monoallelic" = mono_color, "Biallelic" = bi_color)) +  # Keep bar colors
  scale_color_manual(values = c("Monoallelic" = mono_dot_color, "Biallelic" = bi_dot_color)) +  # Different dot colors
  theme_minimal() +
  theme(legend.position = "top")

p_XIST

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "1_with_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p_XIST, height = 8, width = 8)

# Example data
data_ATRX <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP3', 'EXP3', 'EXP3'),
    Cell_type = c('NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC'),
  Monoallelic = c(16.66666667, 78.83211679, 74.07407407, 20, 78.26086957, 82.65306122, 7.407407407, 75.34246575, 78.04878049),
  Biallelic = c(83.33333333, 21.16788321, 25.92592593, 80, 21.73913043, 17.34693878, 92.59259259, 24.65753425, 21.95121951)
)

# Reshape data to long format
library(tidyr)
data_long_ATRX <- data_ATRX %>%
  pivot_longer(cols = c("Biallelic", "Monoallelic"), names_to = "ExpressionType", values_to = "Proportion")

# Load necessary libraries
library(ggplot2)

# Custom color palette for bars
bi_color <- "#94a096"
mono_color <- "#ce4751"

# New color palette for dots (distinct but corresponding)
bi_dot_color <- "#5b6d7c"  # Darker shade of bi_color
mono_dot_color <- "#b03038"  # Darker shade of mono_color

# Plot with stacked bars and bold dots for individual experiments
p_ATRX <- ggplot(data_long_ATRX, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +  # Bars remain the same
  geom_point(aes(x = Cell_type, y = Proportion, color = ExpressionType),  # Different colors for dots
             position = position_jitter(width = 0.2, height = 0), 
             size = 3.5, shape = 21, stroke = 1.8) +  # Bold dots with larger size and stroke
  labs(title = "Proportion of Monoallelic vs. Biallelic Expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("Monoallelic" = mono_color, "Biallelic" = bi_color)) +  # Keep bar colors
  scale_color_manual(values = c("Monoallelic" = mono_dot_color, "Biallelic" = bi_dot_color)) +  # Different dot colors
  theme_minimal() +
  theme(legend.position = "top")

p_ATRX

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "ATRX_with_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p_ATRX, height = 8, width = 8)



# Example data
data_GAT3_POLA1 <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP3', 'EXP3', 'EXP3', 'EXP3', 'EXP3'),
    Cell_type = c('D8', 'D10', 'D12', 'D14', 'D16', 'D8', 'D10', 'D12', 'D14', 'D16', 'D8', 'D10', 'D12', 'D14', 'D16'),
  Monoallelic = c(60.60606061, 64.54545455, 85, 85.71428571, 75.2688172, 79.3814433, 84.81012658, 71.92982456, 80.46875, 89.09090909, 66.66666667, 85.96491228, 81.72043011, 83.33333333, 85.45454545),
  Biallelic = c(39.39393939, 35.45454545, 15, 14.28571429, 24.7311828, 20.6185567, 15.18987342, 28.07017544, 19.53125, 10.90909091, 33.33333333, 14.03508772, 18.27956989, 16.66666667, 14.54545455)
)

# Reshape data to long format
library(tidyr)
data_long_GAT3_POLA1 <- data_GAT3_POLA1 %>%
  pivot_longer(cols = c("Biallelic", "Monoallelic"), names_to = "ExpressionType", values_to = "Proportion")

#view(data_long_GAT3_POLA1)
#tail(data_long_GAT3_POLA1)
#dim(data_long_GAT3_POLA1)
head(data_long_GAT3_POLA1, 30)

# Set factor levels for Cell_type in the desired order
data_long_GAT3_POLA1$Cell_type <- factor(data_long_GAT3_POLA1$Cell_type, levels = c("D8", "D10", "D12", "D14", "D16"))

# Load necessary libraries
library(ggplot2)

# Custom color palette for bars
bi_color <- "#94a096"
mono_color <- "#ce4751"

# New color palette for dots (distinct but corresponding)
bi_dot_color <- "#5b6d7c"  # Darker shade of bi_color
mono_dot_color <- "#b03038"  # Darker shade of mono_color

# Plot with stacked bars and bold dots for individual experiments
p_GAT3_POLA1 <- ggplot(data_long_GAT3_POLA1, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +  # Bars remain the same
  geom_point(aes(x = Cell_type, y = Proportion, color = ExpressionType),  # Different colors for dots
             position = position_jitter(width = 0.2, height = 0), 
             size = 3.5, shape = 21, stroke = 1.8) +  # Bold dots with larger size and stroke
  labs(title = "Proportion of Monoallelic vs. Biallelic Expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("Monoallelic" = mono_color, "Biallelic" = bi_color)) +  # Keep bar colors
  scale_color_manual(values = c("Monoallelic" = mono_dot_color, "Biallelic" = bi_dot_color)) +  # Different dot colors
  theme_minimal() +
  theme(legend.position = "top")

p_GAT3_POLA1

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "GATA3_POLA1_with_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p_GAT3_POLA1, height = 8, width = 8)

# Example data
data_VIM_POLA1 <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP3', 'EXP3', 'EXP3', 'EXP3', 'EXP3'),
    Cell_type = c('D8', 'D10', 'D12', 'D14', 'D16', 'D8', 'D10', 'D12', 'D14', 'D16', 'D8', 'D10', 'D12', 'D14', 'D16'),
  Monoallelic = c(55.93220339, 67.5, 72.36842105, 80.48780488, 81.25, 79.31034483, 81.3559322, 71.69811321, 73.52941176, 80, 63.76811594, 79.54545455, 82.92682927, 94.33962264, 78.20512821),
  Biallelic = c(44.06779661, 32.5, 27.63157895, 19.51219512, 18.75, 20.68965517, 18.6440678, 28.30188679, 26.47058824, 20, 36.23188406, 20.45454545, 17.07317073, 5.660377358, 21.79487179)
)
# Reshape data to long format
library(tidyr)
data_long_VIM_POLA1 <- data_VIM_POLA1 %>%
  pivot_longer(cols = c("Biallelic", "Monoallelic"), names_to = "ExpressionType", values_to = "Proportion")

#view(data_long_GAT3_POLA1)
#tail(data_long_GAT3_POLA1)
#dim(data_long_GAT3_POLA1)
head(data_long_VIM_POLA1, 30)

# Set factor levels for Cell_type in the desired order
data_long_VIM_POLA1$Cell_type <- factor(data_long_VIM_POLA1$Cell_type, levels = c("D8", "D10", "D12", "D14", "D16"))

# Load necessary libraries
library(ggplot2)

# Custom color palette for bars
bi_color <- "#94a096"
mono_color <- "#ce4751"

# New color palette for dots (distinct but corresponding)
bi_dot_color <- "#5b6d7c"  # Darker shade of bi_color
mono_dot_color <- "#b03038"  # Darker shade of mono_color

# Plot with stacked bars and bold dots for individual experiments
p_VIM_POLA1 <- ggplot(data_long_VIM_POLA1, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +  # Bars remain the same
  geom_point(aes(x = Cell_type, y = Proportion, color = ExpressionType),  # Different colors for dots
             position = position_jitter(width = 0.2, height = 0), 
             size = 3.5, shape = 21, stroke = 1.8) +  # Bold dots with larger size and stroke
  labs(title = "Proportion of Monoallelic vs. Biallelic Expression IN VIM-POLA1",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("Monoallelic" = mono_color, "Biallelic" = bi_color)) +  # Keep bar colors
  scale_color_manual(values = c("Monoallelic" = mono_dot_color, "Biallelic" = bi_dot_color)) +  # Different dot colors
  theme_minimal() +
  theme(legend.position = "top")

p_VIM_POLA1

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "VIM_POLA1_with_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p_VIM_POLA1, height = 8, width = 8)

# Example data
data_POLA1 <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP3', 'EXP3', 'EXP3', 'EXP3', 'EXP3', 'EXP3', 'EXP3'),
    Cell_type = c('D0', 'D4', 'D8', 'D10', 'D12', 'D14', 'D16', 'D0', 'D4', 'D8', 'D10', 'D12', 'D14', 'D16', 'D0', 'D4', 'D8', 'D10', 'D12', 'D14', 'D16'),
  Monoallelic = c(23.65591398, 21.72131148, 53.65853659, 75.20661157, 75.2, 78.47222222, 81.64556962, 18.24324324, 15.92920354, 69.8630137, 75.47169811, 70.50359712, 77.48691099, 76.92307692, 34.375, 16.21621622, 66.10169492, 89.70588235, 79.61165049, 94.14893617, 83.16831683),
  Biallelic = c(76.34408602, 78.27868852, 46.34146341, 24.79338843, 24.8, 21.52777778, 18.35443038, 81.75675676, 84.07079646, 30.1369863, 24.52830189, 29.49640288, 22.51308901, 23.07692308, 65.625, 83.78378378, 33.89830508, 10.29411765, 20.38834951, 5.85106383, 16.83168317)
    )
# Reshape data to long format
library(tidyr)
data_long_POLA1 <- data_POLA1 %>%
  pivot_longer(cols = c("Biallelic", "Monoallelic"), names_to = "ExpressionType", values_to = "Proportion")

#view(data_long_GAT3_POLA1)
#tail(data_long_GAT3_POLA1)
dim(data_long_POLA1)
head(data_long_POLA1, 50)

# Set factor levels for Cell_type in the desired order
data_long_POLA1$Cell_type <- factor(data_long_POLA1$Cell_type, levels = c("D0", "D4", "D8", "D10", "D12", "D14", "D16"))

# Load necessary libraries
library(ggplot2)

# Custom color palette for bars
bi_color <- "#94a096"
mono_color <- "#ce4751"

# New color palette for dots (distinct but corresponding)
bi_dot_color <- "#5b6d7c"  # Darker shade of bi_color
mono_dot_color <- "#b03038"  # Darker shade of mono_color

# Plot with stacked bars and bold dots for individual experiments
p_POLA1 <- ggplot(data_long_POLA1, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +  # Bars remain the same
  geom_point(aes(x = Cell_type, y = Proportion, color = ExpressionType),  # Different colors for dots
             position = position_jitter(width = 0.2, height = 0), 
             size = 3.5, shape = 21, stroke = 1.8) +  # Bold dots with larger size and stroke
  labs(title = "Proportion of Monoallelic vs. Biallelic Expression IN POLA1 TIMECOURSE",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("Monoallelic" = mono_color, "Biallelic" = bi_color)) +  # Keep bar colors
  scale_color_manual(values = c("Monoallelic" = mono_dot_color, "Biallelic" = bi_dot_color)) +  # Different dot colors
  theme_minimal() +
  theme(legend.position = "top")

p_POLA1

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "POLA1_TIMECOURSE_with_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p_POLA1, height = 8, width = 8)

# Data (replace with your actual data)
XIST_negative <- c(1.6092811715847157, 0.7377378854006338, 1.019814135700876, 0.8064487668840261, 0.4773598082004101, 1.258494039801081, 1.4429285111512395, 0.9221723567507922, 0.6798760904672507, 0.0180818109166822, 0.4339634620003728, 0.8064487668840261, 1.0595941197175769, 0.8715432861840821, 0.5417601799930056, 0.7963874645897183, 0.7530466501902778, 1.2027075995844725, 0.96433312038755, 0.5146721709933554, 1.814896602976569, 0.9101571023882494, 2.194128728971673, 1.0022563329870604, 1.0076739347869905, 0.9209923059881095, 0.6934530303910472, 1.051014749186431, 0.8343106771892287, 0.482166560193775)

XIST_positive <- c(0.48459253256708296, 0.4014162023503448, 0.4701270838337372, 0.42311437545036346, 0.30377442340026095, 0.5279888787671202, 0.6581779173672321, 0.6292470199005405, 0.06871088148339236, 0.00361636218333644, 0.04339634620003728, 0.0542454327500466, 0.3761016670669898, 0.34717076960029825, 0.28171529359636294, 0.08126402699895084, 0.4604961529940548, 0.7476290483903478, 0.7042882339909073, 0.28171529359636294, 0.7205410393906975, 0.4388257457943346, 0.666365021391397, 0.2762976917964329, 0.5688481889926559, 0.5688481889926559, 0.4388257457943346, 0.2708800899965028, 0.3033857007960832, 0.23837447919692248)


# Base R scatter plot
plot(XIST_negative, type = "b", col = "blue", pch = 19, 
     xlab = "Nucleus Number", ylab = "POLA1 Volume", 
     main = "POLA1 Volume: XIST-negative vs XIST-positive")
lines(XIST_positive, type = "b", col = "red", pch = 19)

# Add legend
legend("topright", legend = c("XIST-negative", "XIST-positive"), 
       col = c("blue", "red"), lty = 1, pch = 19)

# Load ggplot2
library(ggplot2)

# Create a data frame
df <- data.frame(
  Nucleus = 1:30,
  XIST_negative = XIST_negative,
  XIST_positive = XIST_positive
)

# Melt the data for ggplot (long format)
library(reshape2)
df_melt <- melt(df, id.vars = "Nucleus", 
                variable.name = "Condition", 
                value.name = "POLA1_Volume")

# Create the plot
ggplot(df_melt, aes(x = Nucleus, y = POLA1_Volume, color = Condition, group = Condition)) +
  geom_line() +
  geom_point() +
  labs(title = "POLA1 Volume: XIST-negative vs XIST-positive", 
       x = "Nucleus Number", y = "POLA1 Volume") +
  theme_minimal() +
  scale_color_manual(values = c("XIST_negative" = "blue", "XIST_positive" = "red"))

# Data (replace with your actual data)
XIST_negative <- c(1.6092811715847157, 0.7377378854006338, 1.019814135700876, 0.8064487668840261, 0.4773598082004101, 1.258494039801081, 1.4429285111512395, 0.9221723567507922, 0.6798760904672507, 0.0180818109166822, 0.4339634620003728, 0.8064487668840261, 1.0595941197175769, 0.8715432861840821, 0.5417601799930056, 0.7963874645897183, 0.7530466501902778, 1.2027075995844725, 0.96433312038755, 0.5146721709933554, 1.814896602976569, 0.9101571023882494, 2.194128728971673, 1.0022563329870604, 1.0076739347869905, 0.9209923059881095, 0.6934530303910472, 1.051014749186431, 0.8343106771892287, 0.482166560193775)

XIST_positive <- c(0.48459253256708296, 0.4014162023503448, 0.4701270838337372, 0.42311437545036346, 0.30377442340026095, 0.5279888787671202, 0.6581779173672321, 0.6292470199005405, 0.06871088148339236, 0.00361636218333644, 0.04339634620003728, 0.0542454327500466, 0.3761016670669898, 0.34717076960029825, 0.28171529359636294, 0.08126402699895084, 0.4604961529940548, 0.7476290483903478, 0.7042882339909073, 0.28171529359636294, 0.7205410393906975, 0.4388257457943346, 0.666365021391397, 0.2762976917964329, 0.5688481889926559, 0.5688481889926559, 0.4388257457943346, 0.2708800899965028, 0.3033857007960832, 0.23837447919692248)


# Combine the data into a data frame
df <- data.frame(XIST_negative, XIST_positive)

# Create a box plot
boxplot(df, main = "POLA1 Volume Comparison", 
        names = c("XIST-negative", "XIST-positive"),
        col = c("blue", "red"),
        ylab = "POLA1 Volume")


# Load ggplot2
library(ggplot2)

# Convert the data into long format using melt from reshape2 package
library(reshape2)
df_melt <- melt(df, variable.name = "Condition", value.name = "POLA1_Volume")

# Create the box plot
ggplot(df_melt, aes(x = Condition, y = POLA1_Volume, fill = Condition)) +
  geom_boxplot() +
  labs(title = "POLA1 Volume: XIST-negative vs XIST-positive", 
       x = "Condition", y = "POLA1 Volume") +
  theme_minimal() +
  scale_fill_manual(values = c("XIST_negative" = "blue", "XIST_positive" = "red"))


# Load the required libraries
library(ggplot2)
library(reshape2)

# Data (replace with your actual data)
XIST_negative <- c(1.6092811715847157, 0.7377378854006338, 1.019814135700876, 0.8064487668840261, 0.4773598082004101, 1.258494039801081, 1.4429285111512395, 0.9221723567507922, 0.6798760904672507, 0.0180818109166822, 0.4339634620003728, 0.8064487668840261, 1.0595941197175769, 0.8715432861840821, 0.5417601799930056, 0.7963874645897183, 0.7530466501902778, 1.2027075995844725, 0.96433312038755, 0.5146721709933554, 1.814896602976569, 0.9101571023882494, 2.194128728971673, 1.0022563329870604, 1.0076739347869905, 0.9209923059881095, 0.6934530303910472, 1.051014749186431, 0.8343106771892287, 0.482166560193775)

XIST_positive <- c(0.48459253256708296, 0.4014162023503448, 0.4701270838337372, 0.42311437545036346, 0.30377442340026095, 0.5279888787671202, 0.6581779173672321, 0.6292470199005405, 0.06871088148339236, 0.00361636218333644, 0.04339634620003728, 0.0542454327500466, 0.3761016670669898, 0.34717076960029825, 0.28171529359636294, 0.08126402699895084, 0.4604961529940548, 0.7476290483903478, 0.7042882339909073, 0.28171529359636294, 0.7205410393906975, 0.4388257457943346, 0.666365021391397, 0.2762976917964329, 0.5688481889926559, 0.5688481889926559, 0.4388257457943346, 0.2708800899965028, 0.3033857007960832, 0.23837447919692248)


# Create a data frame
df <- data.frame(
  XIST_negative = XIST_negative,
  XIST_positive = XIST_positive
)

# Reshape the data into long format
df_melt <- melt(df, variable.name = "Condition", value.name = "POLA1_Volume")

# Create the violin plot with the box plot overlay
ggplot(df_melt, aes(x = Condition, y = POLA1_Volume, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5) +   # Violin plot (trim = FALSE shows full distribution)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +  # Overlay box plot
  labs(title = "POLA1 Volume: XIST-negative vs XIST-positive", 
       x = "Condition", y = "POLA1 Volume") +
  theme_minimal() +
  scale_fill_manual(values = c("XIST_negative" = "blue", "XIST_positive" = "red"))

# Load the required libraries
library(ggplot2)
library(reshape2)

# Data (replace with your actual data)
XIST_negative <- c(1.6092811715847157, 0.7377378854006338, 1.019814135700876, 0.8064487668840261, 0.4773598082004101, 1.258494039801081, 1.4429285111512395, 0.9221723567507922, 0.6798760904672507, 0.0180818109166822, 0.4339634620003728, 0.8064487668840261, 1.0595941197175769, 0.8715432861840821, 0.5417601799930056, 0.7963874645897183, 0.7530466501902778, 1.2027075995844725, 0.96433312038755, 0.5146721709933554, 1.814896602976569, 0.9101571023882494, 2.194128728971673, 1.0022563329870604, 1.0076739347869905, 0.9209923059881095, 0.6934530303910472, 1.051014749186431, 0.8343106771892287, 0.482166560193775)

XIST_positive <- c(0.48459253256708296, 0.4014162023503448, 0.4701270838337372, 0.42311437545036346, 0.30377442340026095, 0.5279888787671202, 0.6581779173672321, 0.6292470199005405, 0.06871088148339236, 0.00361636218333644, 0.04339634620003728, 0.0542454327500466, 0.3761016670669898, 0.34717076960029825, 0.28171529359636294, 0.08126402699895084, 0.4604961529940548, 0.7476290483903478, 0.7042882339909073, 0.28171529359636294, 0.7205410393906975, 0.4388257457943346, 0.666365021391397, 0.2762976917964329, 0.5688481889926559, 0.5688481889926559, 0.4388257457943346, 0.2708800899965028, 0.3033857007960832, 0.23837447919692248)

# Create a data frame
df <- data.frame(
  XIST_negative = XIST_negative,
  XIST_positive = XIST_positive
)

# Reshape the data into long format
df_melt <- melt(df, variable.name = "Condition", value.name = "POLA1_Volume")

# Create the violin plot with the box plot overlay and individual data points
ggplot(df_melt, aes(x = Condition, y = POLA1_Volume, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5) +   # Violin plot (trim = FALSE shows full distribution)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +  # Overlay box plot
  geom_jitter(width = 0.2, alpha = 0.6, shape = 21, color = "black", fill = NA) +  # Add individual data points
  labs(title = "POLA1 Volume: XIST-negative vs XIST-positive", 
       x = "Condition", y = "POLA1 Volume") +
  theme_minimal() +
  scale_fill_manual(values = c("XIST_negative" = "blue", "XIST_positive" = "red"))


# Load the required libraries
library(ggplot2)
library(reshape2)

# Data (replace with your actual data)
XIST_negative <- c(1.6092811715847157, 0.7377378854006338, 1.019814135700876, 0.8064487668840261, 0.4773598082004101, 1.258494039801081, 1.4429285111512395, 0.9221723567507922, 0.6798760904672507, 0.0180818109166822, 0.4339634620003728, 0.8064487668840261, 1.0595941197175769, 0.8715432861840821, 0.5417601799930056, 0.7963874645897183, 0.7530466501902778, 1.2027075995844725, 0.96433312038755, 0.5146721709933554, 1.814896602976569, 0.9101571023882494, 2.194128728971673, 1.0022563329870604, 1.0076739347869905, 0.9209923059881095, 0.6934530303910472, 1.051014749186431, 0.8343106771892287, 0.482166560193775)

XIST_positive <- c(0.48459253256708296, 0.4014162023503448, 0.4701270838337372, 0.42311437545036346, 0.30377442340026095, 0.5279888787671202, 0.6581779173672321, 0.6292470199005405, 0.06871088148339236, 0.00361636218333644, 0.04339634620003728, 0.0542454327500466, 0.3761016670669898, 0.34717076960029825, 0.28171529359636294, 0.08126402699895084, 0.4604961529940548, 0.7476290483903478, 0.7042882339909073, 0.28171529359636294, 0.7205410393906975, 0.4388257457943346, 0.666365021391397, 0.2762976917964329, 0.5688481889926559, 0.5688481889926559, 0.4388257457943346, 0.2708800899965028, 0.3033857007960832, 0.23837447919692248)

# Create a data frame
df <- data.frame(
  XIST_negative = XIST_negative,
  XIST_positive = XIST_positive
)

# Reshape the data into long format
df_melt <- melt(df, variable.name = "Condition", value.name = "POLA1_Volume")

# Create the violin plot with the box plot overlay and individual data points
p3<- ggplot(df_melt, aes(x = Condition, y = POLA1_Volume, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5) +   # Violin plot (trim = FALSE shows full distribution)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +  # Overlay box plot
  geom_point(position = position_dodge(width = 0.9), shape = 21, color = "black", fill = "black", size = 2) +  # Add individual data points
  labs(title = "POLA1 Volume: XIST-negative vs XIST-positive", 
       x = "Condition", y = "POLA1 Volume") +
  theme_minimal() +
  scale_fill_manual(values = c("XIST_negative" = "blue", "XIST_positive" = "red"))
p3

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "POLA1_volume quantification_naive.pdf"

ggsave(file.path(output_path, filename), plot = p3, height = 8, width = 8)

# Load the required libraries
library(ggplot2)
library(reshape2)

# Sample Data (Add your actual cell identifiers if available)
cell_id <- 1:30  # Example cell identifiers
XIST_negative <- c(1.6092811715847157, 0.7377378854006338, 1.019814135700876, 0.8064487668840261, 0.4773598082004101, 1.258494039801081, 1.4429285111512395, 0.9221723567507922, 0.6798760904672507, 0.0180818109166822, 0.4339634620003728, 0.8064487668840261, 1.0595941197175769, 0.8715432861840821, 0.5417601799930056, 0.7963874645897183, 0.7530466501902778, 1.2027075995844725, 0.96433312038755, 0.5146721709933554, 1.814896602976569, 0.9101571023882494, 2.194128728971673, 1.0022563329870604, 1.0076739347869905, 0.9209923059881095, 0.6934530303910472, 1.051014749186431, 0.8343106771892287, 0.482166560193775)
XIST_positive <- c(0.48459253256708296, 0.4014162023503448, 0.4701270838337372, 0.42311437545036346, 0.30377442340026095, 0.5279888787671202, 0.6581779173672321, 0.6292470199005405, 0.06871088148339236, 0.00361636218333644, 0.04339634620003728, 0.0542454327500466, 0.3761016670669898, 0.34717076960029825, 0.28171529359636294, 0.08126402699895084, 0.4604961529940548, 0.7476290483903478, 0.7042882339909073, 0.28171529359636294, 0.7205410393906975, 0.4388257457943346, 0.666365021391397, 0.2762976917964329, 0.5688481889926559, 0.5688481889926559, 0.4388257457943346, 0.2708800899965028, 0.3033857007960832, 0.23837447919692248)

# Create a data frame with cell IDs
df <- data.frame(
  Cell_ID = cell_id,
  XIST_negative = XIST_negative,
  XIST_positive = XIST_positive
)

# Reshape the data into long format for plotting
df_melt <- melt(df, id.vars = "Cell_ID", variable.name = "Condition", value.name = "POLA1_Volume")

# Create the bar plot
p4<- ggplot(df_melt, aes(x = Cell_ID, y = POLA1_Volume, fill = Condition)) +
  geom_col(position = "dodge") +  # Use geom_col() for a bar plot
  labs(title = "POLA1 Volume: XIST-negative vs XIST-positive by Cell", 
       x = "Cell ID", y = "POLA1 Volume") +
  theme_minimal() +
  scale_fill_manual(values = c("XIST_negative" = "blue", "XIST_positive" = "red"))
p4

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "POLA1_volume quantification_naive_bar_plot.pdf"

ggsave(file.path(output_path, filename), plot = p4, height = 8, width = 8)

# Load the required libraries
library(ggplot2)
library(reshape2)

# Data (replace with your actual data)
XIST_negative <- c(0.768566494, 0.647619048, 0.684466019, 0.655882353, 0.611111111, 0.704453441, 0.686746988, 0.594405594, 0.90821256, 0.833333333, 0.909090909, 0.93697479, 0.738035264, 0.715133531, 0.657894737, 0.907407407, 0.620535714, 0.616666667, 0.577922078, 0.646258503, 0.715811966, 0.674698795, 0.767045455, 0.783898305, 0.639175258, 0.618181818, 0.612440191, 0.795081967, 0.733333333, 0.669172932)
XIST_positive <- c(0.231433506, 0.352380952, 0.315533981, 0.344117647, 0.388888889, 0.295546559, 0.313253012, 0.405594406, 0.09178744, 0.166666667, 0.090909091, 0.06302521, 0.261964736, 0.284866469, 0.342105263, 0.092592593, 0.379464286, 0.383333333, 0.422077922, 0.353741497, 0.284188034, 0.325301205, 0.232954545, 0.216101695, 0.360824742, 0.381818182, 0.387559809, 0.204918033, 0.266666667, 0.330827068)

# Create a data frame
df <- data.frame(
  XIST_negative = XIST_negative,
  XIST_positive = XIST_positive
)

# Reshape the data into long format
df_melt <- melt(df, variable.name = "Condition", value.name = "POLA1_Volume")

# Create the violin plot with the box plot overlay and individual data points
p1<- ggplot(df_melt, aes(x = Condition, y = POLA1_Volume, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5) +   # Violin plot (trim = FALSE shows full distribution)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +  # Overlay box plot
  geom_point(position = position_dodge(width = 0.9), shape = 21, color = "black", fill = "black", size = 2) +  # Add individual data points
  labs(title = "POLA1 Volume: XIST-negative vs XIST-positive", 
       x = "Condition", y = "Proportion of POLA1 Volume") +
  theme_minimal() +
  scale_fill_manual(values = c("XIST_negative" = "blue", "XIST_positive" = "red"))
p1

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "proportion_of_POLA1_volume quantification_naive.pdf"

ggsave(file.path(output_path, filename), plot = p1, height = 8, width = 8)

# Load the required libraries
library(ggplot2)
library(reshape2)

# Sample Data (Add your actual cell identifiers if available)
cell_id <- 1:30  # Example cell identifiers
XIST_negative <- c(0.768566494, 0.647619048, 0.684466019, 0.655882353, 0.611111111, 0.704453441, 0.686746988, 0.594405594, 0.90821256, 0.833333333, 0.909090909, 0.93697479, 0.738035264, 0.715133531, 0.657894737, 0.907407407, 0.620535714, 0.616666667, 0.577922078, 0.646258503, 0.715811966, 0.674698795, 0.767045455, 0.783898305, 0.639175258, 0.618181818, 0.612440191, 0.795081967, 0.733333333, 0.669172932)
XIST_positive <- c(0.231433506, 0.352380952, 0.315533981, 0.344117647, 0.388888889, 0.295546559, 0.313253012, 0.405594406, 0.09178744, 0.166666667, 0.090909091, 0.06302521, 0.261964736, 0.284866469, 0.342105263, 0.092592593, 0.379464286, 0.383333333, 0.422077922, 0.353741497, 0.284188034, 0.325301205, 0.232954545, 0.216101695, 0.360824742, 0.381818182, 0.387559809, 0.204918033, 0.266666667, 0.330827068)
# Create a data frame with cell IDs
df <- data.frame(
  Cell_ID = cell_id,
  XIST_negative = XIST_negative,
  XIST_positive = XIST_positive
)

# Reshape the data into long format for plotting
df_melt <- melt(df, id.vars = "Cell_ID", variable.name = "Condition", value.name = "POLA1_Volume")

# Create the bar plot
p2<- ggplot(df_melt, aes(x = Cell_ID, y = POLA1_Volume, fill = Condition)) +
  geom_col(position = "dodge") +  # Use geom_col() for a bar plot
  labs(title = "POLA1 Volume: XIST-negative vs XIST-positive by Cell", 
       x = "Cell ID", y = "POLA1 Volume") +
  theme_minimal() +
  scale_fill_manual(values = c("XIST_negative" = "blue", "XIST_positive" = "red"))
p2

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "proportion_of_POLA1_volume quantification_naive_bar_plot.pdf"

ggsave(file.path(output_path, filename), plot = p2, height = 8, width = 8)

# Load the required libraries
library(ggplot2)
library(reshape2)

# Data (replace with your actual data)
XIST_negative <- c(2.656734645, 1.463255941, 2.211724817, 1.759499237, 1.736872669, 2.470629942, 2.891704604, 2.656734645, 1.948810124, 0.402375425, 1.642707154, 1.773734291, 2.035842452, 2.091913399, 1.415068495, 1.959837148, 1.808236784, 2.105565991, 2.095248751, 1.653370066, 2.901213663, 1.858108535, 2.473708781, 2.183585871, 1.872816405, 1.690363592, 1.57057444, 1.724571066, 1.762339585, 1.333868086)
XIST_positive <- c(1.360811881, 1.221402987, 1.534118431, 1.256430638, 1.147830731, 1.830794629, 1.834093083, 1.582531727, 0.625880713, 0.1, 0.455142272, 0.445518569, 1.360811881, 1.28129882, 1.412766132, 0.750620135, 1.774121929, 1.825945791, 1.904187645, 1.34884529, 2.02835221, 1.651455912, 1.426915446, 1.49639991, 1.460288084, 1.514378029, 1.406167354, 1.199578367, 1.333868086, 1.075031316)

# Create a data frame
df <- data.frame(
  XIST_negative = XIST_negative,
  XIST_positive = XIST_positive
)

# Reshape the data into long format
df_melt <- melt(df, variable.name = "Condition", value.name = "POLA1_maxferet_diameter")

# Create the violin plot with the box plot overlay and individual data points
p1<- ggplot(df_melt, aes(x = Condition, y = POLA1_maxferet_diameter, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5) +   # Violin plot (trim = FALSE shows full distribution)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +  # Overlay box plot
  geom_point(position = position_dodge(width = 0.9), shape = 21, color = "black", fill = "black", size = 2) +  # Add individual data points
  labs(title = "POLA1 max_feret_diameter: XIST-negative vs XIST-positive", 
       x = "Condition", y = "POLA1 max feret diameter") +
  theme_minimal() +
  scale_fill_manual(values = c("XIST_negative" = "blue", "XIST_positive" = "red"))
p1

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "POLA1_max_feret_diameter_quantification_naive.pdf"

ggsave(file.path(output_path, filename), plot = p1, height = 8, width = 8)

# Load the required libraries
library(ggplot2)
library(reshape2)

# Sample Data (Add your actual cell identifiers if available)
cell_id <- 1:30  # Example cell identifiers
XIST_negative <- c(2.656734645, 1.463255941, 2.211724817, 1.759499237, 1.736872669, 2.470629942, 2.891704604, 2.656734645, 1.948810124, 0.402375425, 1.642707154, 1.773734291, 2.035842452, 2.091913399, 1.415068495, 1.959837148, 1.808236784, 2.105565991, 2.095248751, 1.653370066, 2.901213663, 1.858108535, 2.473708781, 2.183585871, 1.872816405, 1.690363592, 1.57057444, 1.724571066, 1.762339585, 1.333868086)
XIST_positive <- c(1.360811881, 1.221402987, 1.534118431, 1.256430638, 1.147830731, 1.830794629, 1.834093083, 1.582531727, 0.625880713, 0.1, 0.455142272, 0.445518569, 1.360811881, 1.28129882, 1.412766132, 0.750620135, 1.774121929, 1.825945791, 1.904187645, 1.34884529, 2.02835221, 1.651455912, 1.426915446, 1.49639991, 1.460288084, 1.514378029, 1.406167354, 1.199578367, 1.333868086, 1.075031316)

# Create a data frame with cell IDs
df <- data.frame(
  Cell_ID = cell_id,
  XIST_negative = XIST_negative,
  XIST_positive = XIST_positive
)

# Reshape the data into long format for plotting
df_melt <- melt(df, id.vars = "Cell_ID", variable.name = "Condition", value.name = "POLA1_max_feret_diameter")

# Create the bar plot
p4<- ggplot(df_melt, aes(x = Cell_ID, y = POLA1_max_feret_diameter, fill = Condition)) +
  geom_col(position = "dodge") +  # Use geom_col() for a bar plot
  labs(title = "POLA1 max_feret_diameter: XIST-negative vs XIST-positive by Cell", 
       x = "Cell ID", y = "POLA1_max_feret_diameter") +
  theme_minimal() +
  scale_fill_manual(values = c("XIST_negative" = "blue", "XIST_positive" = "red"))
p4

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "POLA1_max_feret_diameter_quantification_naive_bar_plot.pdf"

ggsave(file.path(output_path, filename), plot = p4, height = 8, width = 8)

# Load the required libraries
library(ggplot2)
library(reshape2)

# Data (replace with your actual data)
XIST_negative <- c(0.661282857, 0.545043516, 0.590447777, 0.583401906, 0.602097487, 0.574374815, 0.611897672, 0.626696794, 0.756910343, 0.800945677, 0.783043403, 0.799248397, 0.599366981, 0.620154696, 0.500407089, 0.72306513, 0.504761507, 0.535561409, 0.523886004, 0.550716677, 0.588533298, 0.529441349, 0.634182797, 0.593368019, 0.561883497, 0.52745706, 0.527615275, 0.589768446, 0.569192952, 0.55372511)
XIST_positive <- c(0.338717143, 0.454956484, 0.409552223, 0.416598094, 0.397902513, 0.425625185, 0.388102328, 0.373303206, 0.243089657, 0.199054323, 0.216956597, 0.200751603, 0.400633019, 0.379845304, 0.499592911, 0.27693487, 0.495238493, 0.464438591, 0.476113996, 0.449283323, 0.411466702, 0.470558651, 0.365817203, 0.406631981, 0.438116503, 0.47254294, 0.472384725, 0.410231554, 0.430807048, 0.44627489)

# Create a data frame
df <- data.frame(
  XIST_negative = XIST_negative,
  XIST_positive = XIST_positive
)

# Reshape the data into long format
df_melt <- melt(df, variable.name = "Condition", value.name = "POLA1_maxferet_diameter")

# Create the violin plot with the box plot overlay and individual data points
p2<- ggplot(df_melt, aes(x = Condition, y = POLA1_maxferet_diameter, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5) +   # Violin plot (trim = FALSE shows full distribution)
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +  # Overlay box plot
  geom_point(position = position_dodge(width = 0.9), shape = 21, color = "black", fill = "black", size = 2) +  # Add individual data points
  labs(title = "POLA1 Volume: XIST-negative vs XIST-positive", 
       x = "Condition", y = "Proportion of POLA1 max feret diameter") +
  theme_minimal() +
  scale_fill_manual(values = c("XIST_negative" = "blue", "XIST_positive" = "red"))
p2

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "proportion_of_POLA1_max_feret_diameter_quantification_naive.pdf"

ggsave(file.path(output_path, filename), plot = p2, height = 8, width = 8)

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "Proportion_of_POLA1_max_feret_diameter_quantification_naive_bar_plot.pdf"

ggsave(file.path(output_path, filename), plot = p6, height = 8, width = 8)

library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
library(broom) 

# Example data
data_THOC2 <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP3', 'EXP3', 'EXP3'),
    Cell_type = c('NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC'),
  Monoallelic = c(13, 108, 100, 23, 90, 81, 6, 110, 64),
  Biallelic = c(65, 29, 35, 92, 25, 17, 75, 36, 18)
)

print(data_THOC2)

# Load necessary libraries
library(ggplot2)

# Example data
df <- data.frame(
  Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP3', 'EXP3', 'EXP3'),
  Cell_type = c('NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC'),
  Monoallelic = c(13, 108, 100, 23, 90, 81, 6, 110, 64),
  Biallelic = c(65, 29, 35, 92, 25, 17, 75, 36, 18)
)

# Set 'NAIVE' as the reference level for Cell_type
df$Cell_type <- factor(df$Cell_type, levels = c('NAIVE', 'EXMC', 'TSC'))

# Create a binary outcome variable (using the total counts)
df$total_counts <- df$Monoallelic + df$Biallelic
df$proportion_monoallelic <- df$Monoallelic / df$total_counts

# Logistic regression (using the binomial family), with NAIVE as the reference level
logistic_model <- glm(cbind(Monoallelic, Biallelic) ~ Cell_type, data = df, family = binomial)

# View summary of the model
summary(logistic_model)

# Get the odds ratios (exponentiated coefficients)
exp(coef(logistic_model))

# Optional: Predict probabilities based on the model
df$predicted_prob <- predict(logistic_model, type = "response")

# Plot predicted probabilities for monoallelic expression by Cell_type
ggplot(df, aes(x = Cell_type, y = predicted_prob, fill = Cell_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Predicted Probability of Monoallelic Expression by Cell Type",
       x = "Cell Type", y = "Predicted Probability") +
  theme_minimal()

# Fit the logistic regression model (as done previously)
logistic_model <- glm(cbind(Monoallelic, Biallelic) ~ Cell_type, data = df, family = binomial)

# Extract the model summary
model_summary <- summary(logistic_model)

# Extract the p-values from the model summary
p_values <- coef(model_summary)[, "Pr(>|z|)"]

# Display the exact p-values
p_values




# Load necessary libraries
library(tidyverse)
library(ggpubr)

# Your original data
data_XIST_qPCR <- data.frame(
  Experiment = rep(c('EXP1', 'EXP2', 'EXP3'), each = 4),
  Cell_type = rep(c('WT-nodox', 'WT-DOX', 'XIST-IKD-NODOX', 'XIST-IKD-DOX'), times = 3),
  expression = c(0.030014314, 0.012718136, 0.012383474, 0.000382648, 0.013696084, 0.01595841, 0.012135414, 0.001382585, 0.023215975, 0.014499767, 0.009308463, 0.000202273)
)

# Long format (though not strictly needed here, kept for completeness)
data_long_XIST <- data_XIST_qPCR %>%
  pivot_longer(cols = c("expression"), names_to = "ExpressionType", values_to = "Proportion")

# Define comparisons
my_comparisons <- list(
  c("WT-nodox", "WT-DOX"),
  c("XIST-IKD-NODOX", "XIST-IKD-DOX")
)

# Plot
p_XIST <- ggplot(data_long_XIST, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +
  geom_point(aes(color = ExpressionType), position = position_jitter(width = 0.2), 
             size = 3.5, shape = 21, stroke = 1.8) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  labs(title = "XIST expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("expression" = "#ce4751")) +
  scale_color_manual(values = c("expression" = "#b03038")) +
  theme_minimal() +
  theme(legend.position = "top")

# Display plot
p_XIST

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "XIST_qPCR.pdf"

ggsave(file.path(output_path, filename), plot = p_XIST, height = 8, width = 8)

# Load necessary libraries
library(tidyverse)
library(ggpubr)

# Your original data
data_NANOG_qPCR <- data.frame(
  Experiment = rep(c('EXP1', 'EXP2', 'EXP3'), each = 4),
  Cell_type = rep(c('WT-nodox', 'WT-DOX', 'XIST-IKD-NODOX', 'XIST-IKD-DOX'), times = 3),
  expression = c(0.139499598, 0.124999627, 0.108453415, 0.123874351, 0.131822081, 0.135387757, 0.107319926, 0.132274791, 0.139642099, 0.11445616, 0.117093842, 0.139779703)
)

# Long format (though not strictly needed here, kept for completeness)
data_long_NANOG <- data_NANOG_qPCR %>%
  pivot_longer(cols = c("expression"), names_to = "ExpressionType", values_to = "Proportion")

# Define comparisons
my_comparisons <- list(
  c("WT-nodox", "WT-DOX"),
  c("XIST-IKD-NODOX", "XIST-IKD-DOX")
)

# Plot
p_NANOG <- ggplot(data_long_NANOG, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +
  geom_point(aes(color = ExpressionType), position = position_jitter(width = 0.2), 
             size = 3.5, shape = 21, stroke = 1.8) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  labs(title = "NANOG expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("expression" = "#ce4751")) +
  scale_color_manual(values = c("expression" = "#b03038")) +
  theme_minimal() +
  theme(legend.position = "top")

# Display plot
p_NANOG

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "NANOG_qPCR.pdf"

ggsave(file.path(output_path, filename), plot = p_NANOG, height = 8, width = 8)

# Load necessary libraries
library(tidyverse)
library(ggpubr)

# Your original data
data_OCT4_qPCR <- data.frame(
  Experiment = rep(c('EXP1', 'EXP2', 'EXP3'), each = 4),
  Cell_type = rep(c('WT-nodox', 'WT-DOX', 'XIST-IKD-NODOX', 'XIST-IKD-DOX'), times = 3),
  expression = c(0.61624113, 0.654564598, 0.573677141, 0.559986025, 0.63218087, 0.64889699, 0.52619785, 0.57157933, 0.663074905, 0.625763963, 0.570474509, 0.619766427)
)

# Long format (though not strictly needed here, kept for completeness)
data_long_OCT4 <- data_OCT4_qPCR %>%
  pivot_longer(cols = c("expression"), names_to = "ExpressionType", values_to = "Proportion")

# Define comparisons
my_comparisons <- list(
  c("WT-nodox", "WT-DOX"),
  c("XIST-IKD-NODOX", "XIST-IKD-DOX")
)

# Plot
p_OCT4 <- ggplot(data_long_OCT4, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +
  geom_point(aes(color = ExpressionType), position = position_jitter(width = 0.2), 
             size = 3.5, shape = 21, stroke = 1.8) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  labs(title = "OCT4 expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("expression" = "#ce4751")) +
  scale_color_manual(values = c("expression" = "#b03038")) +
  theme_minimal() +
  theme(legend.position = "top")

# Display plot
p_OCT4

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "OCT4_qPCR.pdf"

ggsave(file.path(output_path, filename), plot = p_OCT4, height = 8, width = 8)

# Load necessary libraries
library(tidyverse)
library(ggpubr)

# Your original data
data_KLF4_qPCR <- data.frame(
  Experiment = rep(c('EXP1', 'EXP2', 'EXP3'), each = 4),
  Cell_type = rep(c('WT-nodox', 'WT-DOX', 'XIST-IKD-NODOX', 'XIST-IKD-DOX'), times = 3),
  expression = c(0.03499367, 0.03161026, 0.022204256, 0.039993015, 0.032579541, 0.036579794, 0.026674901, 0.032199797, 0.035408137, 0.027786742, 0.028632478, 0.030339067)
)

# Long format (though not strictly needed here, kept for completeness)
data_long_KLF4 <- data_KLF4_qPCR %>%
  pivot_longer(cols = c("expression"), names_to = "ExpressionType", values_to = "Proportion")

# Define comparisons
my_comparisons <- list(
  c("WT-nodox", "WT-DOX"),
  c("XIST-IKD-NODOX", "XIST-IKD-DOX")
)

# Plot
p_KLF4 <- ggplot(data_long_KLF4, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +
  geom_point(aes(color = ExpressionType), position = position_jitter(width = 0.2), 
             size = 3.5, shape = 21, stroke = 1.8) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  labs(title = "KLF4 expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("expression" = "#ce4751")) +
  scale_color_manual(values = c("expression" = "#b03038")) +
  theme_minimal() +
  theme(legend.position = "top")

# Display plot
p_KLF4

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "KLF4_qPCR.pdf"

ggsave(file.path(output_path, filename), plot = p_KLF4, height = 8, width = 8)

# Load necessary libraries
library(tidyverse)
library(ggpubr)

# Your original data
data_LBP9_qPCR <- data.frame(
  Experiment = rep(c('EXP1', 'EXP2', 'EXP3'), each = 4),
  Cell_type = rep(c('WT-nodox', 'WT-DOX', 'XIST-IKD-NODOX', 'XIST-IKD-DOX'), times = 3),
  expression = c(0.037313403, 0.033966206, 0.02501069, 0.027077715,
                 0.033582705, 0.034315127, 0.02766662, 0.031890714,
                 0.036288658, 0.03439112, 0.028966591, 0.028703949)
)

# Long format (though not strictly needed here, kept for completeness)
data_long_LBP9 <- data_LBP9_qPCR %>%
  pivot_longer(cols = c("expression"), names_to = "ExpressionType", values_to = "Proportion")

# Define comparisons
my_comparisons <- list(
  c("WT-nodox", "WT-DOX"),
  c("XIST-IKD-NODOX", "XIST-IKD-DOX")
)

# Plot
p_LBP9 <- ggplot(data_long_LBP9, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +
  geom_point(aes(color = ExpressionType), position = position_jitter(width = 0.2), 
             size = 3.5, shape = 21, stroke = 1.8) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.format") +
  labs(title = "LBP9 expression",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("expression" = "#ce4751")) +
  scale_color_manual(values = c("expression" = "#b03038")) +
  theme_minimal() +
  theme(legend.position = "top")

# Display plot
p_LBP9

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "LBP9_qPCR.pdf"

ggsave(file.path(output_path, filename), plot = p_LBP9, height = 8, width = 8)

library(tidyverse)

data_path <- "/lustre1/project/stg_00041/Amitesh/DAPI_COMPACTION/"  # <-- change this to your actual path

all_data <- list.files(path = data_path, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE) %>%
  set_names() %>%
  map_dfr(~ read_table(.x, col_names = c("Distance", "DAPI", "XIST")) %>%
            mutate(
              Nucleus = tools::file_path_sans_ext(basename(.x)),
              Condition = basename(dirname(.x))  # gets the folder name
            ),
          .id = "File")

all_data

df_clean <- all_data %>% filter(Distance != "Distance_(microns)")


df_clean

df_clean <- df_clean %>%
  mutate(
    Distance = as.numeric(Distance),
    DAPI = as.numeric(DAPI),
    XIST = as.numeric(XIST)
  )

# Define the XIST cloud as top N% of XIST intensity within each nucleus
library(dplyr)

df_cloud <- df_clean %>%
  group_by(Nucleus) %>%
  mutate(
    XIST_threshold = quantile(XIST, 0.92, na.rm = TRUE),  # top 10% as cloud
    in_cloud = XIST >= XIST_threshold
  )

df_cloud

df_dapi_summary <- df_cloud %>%
  group_by(Nucleus, Condition, in_cloud) %>%
  summarise(mean_DAPI = mean(DAPI, na.rm = TRUE), .groups = "drop")

ggplot(df_dapi_summary, aes(x = in_cloud, y = mean_DAPI, fill = in_cloud)) +
  geom_boxplot() +
  facet_wrap(~Condition) +
  labs(x = "Region", y = "Mean DAPI intensity", fill = "In XIST cloud") +
  scale_x_discrete(labels = c("Outside", "Inside")) +
  theme_minimal()


library(tidyr)
df_wide <- df_dapi_summary %>%
  pivot_wider(names_from = in_cloud, values_from = mean_DAPI, 
              names_prefix = "DAPI_")


df_wide <- df_wide %>%
  mutate(DAPI_ratio = DAPI_TRUE / DAPI_FALSE)


library(ggplot2)
library(ggpubr)

ggplot(df_wide, aes(x = Condition, y = DAPI_ratio, fill = Condition)) +
  geom_boxplot() +
  scale_fill_manual(values = c("naive" = "#66c2a5", 
                               "TSC" = "#fc8d62", 
                               "exmc" = "#8da0cb")) +
  stat_compare_means(comparisons = list(c("naive", "TSC"), 
                                        c("naive", "exmc")),
                     method = "wilcox.test", 
                     label = "p.signif") +
  labs(y = "Inside/Outside DAPI intensity ratio",
       x = "Cell type (Condition)") +
  theme_minimal()


# Load the required libraries
library(ggplot2)
library(reshape2)

# Data (replace with your actual data)
Naive <- c(364.301,236.478,203.975,269.958,249.514,326.386,260.101,210.923,204.364,220.714,221.624,386.758,232.999,174.716,191.737,224.979,211.239,274.56,289.064,221.542,301.187,382.953,233.862,201.7,233.896,495.489,455.798,441.741,354.3,215.957)
EXMC <- c(434.082,365.727,504.937,377.553,297.825,414.184,494.516,585.145,610.914,122.785,290.575,395.014,375.383,291.352,375.872,568.923,314.993,306.076,446.523,407.905,775.817,311.366,418.32,406.051,307.11,300.916,458.304,495.497,477.594,397.441)
TSC <- c(714.122,619.363,583.782,617.915,574.45,464.666,554.5,736.261,850.423,1033.823,839.665,674.791,675.285,542.392,635.396,683.933,475.829,371.677,626.702,775.681,646.082,504.641,599.586,526.511,311.089,291.384,201.405,382.347,304.988,399.949)
max_length <- max(length(Naive), length(EXMC), length(TSC))

Naive <- c(Naive, rep(NA, max_length - length(Naive)))
EXMC <- c(EXMC, rep(NA, max_length - length(EXMC)))
TSC <- c(TSC, rep(NA, max_length - length(TSC)))

# Create a data frame
df_1 <- data.frame(
  Naive = Naive,
  EXMC = EXMC,
    TSC = TSC
)

df_1

df_melt_1 <- melt(df_1, variable.name = "Condition", value.name = "H3K27ME3_gray_intensity_value")

df_melt_1


library(ggpubr)

# Define comparisons
my_comparisons <- list(c("Naive", "EXMC"), c("Naive", "TSC"), c("EXMC", "TSC"))

# Create plot
p3 <- ggplot(df_melt_1, aes(x = Condition, y = H3K27ME3_gray_intensity_value, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.9), shape = 21, color = "black", fill = "black", size = 2) +
  labs(title = "H3K27ME3_gray_intensity_value (Naive, EXMC, TSC)", 
       x = "Condition", y = "H3K27ME3_gray_intensity_value") +
  theme_minimal() +
  scale_fill_manual(values = c("Naive" = "blue", "EXMC" = "green", "TSC" = "red")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format")
p3

output_path <- "/lustre1/project/stg_00041/Amitesh/H3K27me3_gray_value_quantificaion/"
filename <- "H3K27ME3_gray_intensity_value_naive_EXMC_TSC.pdf"

ggsave(file.path(output_path, filename), plot = p3, height = 8, width = 8)



# Load the required libraries
library(ggplot2)
library(reshape2)

# Data (replace with your actual data)
Naive <- c(705.72,551.344,454.948,449.968,381.734,821.647,678.31,668.529,529.853,558.285,874.441,768.813,514.835,555.5,481.341,616.167,384.022,517.63,511.282,617.648,657.491,659,460.508,378.966,524.252,785.605,943.75,941.741,848.875,663.819)
EXMC <- c(1416.2,1146.569,1209.77,1311.933,927.148,1150.396,1375.516,1530.055,1512.861,991.686,877.618,787.31,1082.238,704.771,1164.226,937.978,713.927,720.841,1230.264,1190.25,1605.018,704.704,1092.957,943.464,1117.345,683.165,971.214,1387.361,1079.163,1128.094)
TSC <- c(1402.512,1183.578,1664.908,1183.092,1754.532,1227.295,1486.684,1410.459,1580.27,1634.736,1785.114,1243.853,1077.824,1152,1142.384,1428.145,974.568,1057.507,1365.966,1554.172,1510.718,1299.842,1341.282,1365.183,818.605,709.89,874.74,698.524,664.112,707.025)
max_length <- max(length(Naive), length(EXMC), length(TSC))

Naive <- c(Naive, rep(NA, max_length - length(Naive)))
EXMC <- c(EXMC, rep(NA, max_length - length(EXMC)))
TSC <- c(TSC, rep(NA, max_length - length(TSC)))

# Create a data frame
df_2 <- data.frame(
  Naive = Naive,
  EXMC = EXMC,
    TSC = TSC
)

df_2

df_melt_2 <- melt(df_2, variable.name = "Condition", value.name = "H3K27ME3_gray_intensity_value_inside_Xi")

df_melt_2


library(ggpubr)

# Define comparisons
my_comparisons <- list(c("Naive", "EXMC"), c("Naive", "TSC"), c("EXMC", "TSC"))

# Create plot
p4 <- ggplot(df_melt_2, aes(x = Condition, y = H3K27ME3_gray_intensity_value_inside_Xi, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA) +
  geom_point(position = position_dodge(width = 0.9), shape = 21, color = "black", fill = "black", size = 2) +
  labs(title = "H3K27ME3_gray_intensity_value_inside_Xi (Naive, EXMC, TSC)", 
       x = "Condition", y = "H3K27ME3_gray_intensity_value_inside_Xi") +
  theme_minimal() +
  scale_fill_manual(values = c("Naive" = "blue", "EXMC" = "green", "TSC" = "red")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format")
p4

output_path <- "/lustre1/project/stg_00041/Amitesh/H3K27me3_gray_value_quantificaion/"
filename <- "H3K27ME3_gray_intensity_value_inside_Xi_naive_EXMC_TSC.pdf"

ggsave(file.path(output_path, filename), plot = p4, height = 8, width = 8)



# Example data
data_THOC2 <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP3', 'EXP3', 'EXP3'),
    Cell_type = c('NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC', 'NAIVE', 'EXMC', 'TSC'),
  Monoallelic = c(71.75572519,95.39170507,96.36363636,73.29192547,97.20670391,98.16513761,75.75757576,97.34513274,96.17486339),
  Biallelic = c(1.526717557,0,0,0,0,0,0.432900433,0,0),
    negative = c(26.71755725,4.608294931,3.636363636,26.70807453,2.793296089,1.834862385,23.80952381,2.654867257,3.825136612)
)

# Reshape data to long format
library(tidyr)
data_long <- data_THOC2 %>%
  pivot_longer(cols = c("Biallelic", "Monoallelic", "negative"), names_to = "ExpressionType", values_to = "Proportion")

print(data_long)

# Set desired stacking order BEFORE plotting
data_long$ExpressionType <- factor(data_long$ExpressionType,
                                   levels = c("Monoallelic", "Biallelic", "negative"))

# Then replot
p_THOC2 <- ggplot(data_long, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +
  geom_point(aes(color = ExpressionType),
             position = position_jitter(width = 0.2), size = 2, shape = 21, stroke = 1.5) +
  labs(title = "Proportion of Monoallelic, Biallelic, and Negative Expression",
       x = "Cell Type", y = "Proportion (%)") +
  scale_fill_manual(values = c("Monoallelic" = mono_color,
                               "Biallelic" = bi_color,
                               "negative" = neg_color)) +
  scale_color_manual(values = c("Monoallelic" = mono_color,
                                "Biallelic" = bi_color,
                                "negative" = neg_color)) +
  theme_minimal() +
  theme(legend.position = "top") +
  stat_compare_means(aes(group = ExpressionBinary), 
                     method = "wilcox.test", 
                     label = "p.signif", 
                     label.y = 120, 
                     hide.ns = TRUE)

p_THOC2


# Define comparisons
my_comparisons <- list(c("NAIVE", "EXMC"), c("NAIVE", "TSC"))

library(ggpubr)

# Ensure correct stacking order
data_long$ExpressionType <- factor(data_long$ExpressionType,
                                   levels = c("Monoallelic", "Biallelic", "negative"))

# Plot
p_THOC2 <- ggplot(data_long, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +
  geom_point(aes(color = ExpressionType),
             position = position_jitter(width = 0.2), size = 2, shape = 21, stroke = 1.5) +
  scale_fill_manual(values = c("Monoallelic" = mono_color,
                               "Biallelic" = bi_color,
                               "negative" = neg_color)) +
  scale_color_manual(values = c("Monoallelic" = mono_color,
                                "Biallelic" = bi_color,
                                "negative" = neg_color)) +
  labs(title = "Proportion of Monoallelic, Biallelic, and Negative Expression",
       x = "Cell Type", y = "Proportion (%)") +
  theme_minimal() +
  theme(legend.position = "top") +
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     aes(group = Cell_type),
                     label = "p.signif",
                     label.y = c(110, 120)) +  # Adjust for visibility
  ylim(0, 130)

p_THOC2

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "XIST_with_negative cells_with_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p_THOC2, height = 8, width = 8)



# Example data
data_pol2 <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2'),
    Cell_type = c('D0', 'D4', 'D8', 'D10', 'D12', 'D14', 'D16', 'D0', 'D4', 'D8', 'D10', 'D12', 'D14', 'D16'),
  exclusion = c(3, 16, NA, 42, NA, 45, 80, 19.5122, 31.3752, 67.3913, 66.9014, 61.6, 75, 80),
  noexclusion = c(97,84, NA, 58, NA, 55, 20, 80.4878, 68.6275, 32.6087, 33.098, 38.4, 25, 20)
    )
# Reshape data to long format
library(tidyr)
data_long_pol2 <- data_pol2 %>%
  pivot_longer(cols = c("exclusion", "noexclusion"), names_to = "ExpressionType", values_to = "Proportion")

#view(data_long_GAT3_POLA1)
#tail(data_long_GAT3_POLA1)
dim(data_long_pol2)
head(data_long_pol2, 50)

# Set factor levels for Cell_type in the desired order
data_long_pol2$Cell_type <- factor(data_long_pol2$Cell_type, levels = c("D0", "D4", "D8", "D10", "D12", "D14", "D16"))

# Load necessary libraries
library(ggplot2)

# Custom color palette for bars
exclusion <- "#94a096"
no_exclusion <- "#ce4751"

# New color palette for dots (distinct but corresponding)
exclusion_dot_color <- "#5b6d7c"  # Darker shade of bi_color
noexclusion_dot_color <- "#b03038"  # Darker shade of mono_color

# Plot with stacked bars and bold dots for individual experiments
p_pol2 <- ggplot(data_long_pol2, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +  # Bars remain the same
  geom_point(aes(x = Cell_type, y = Proportion, color = ExpressionType),  # Different colors for dots
             position = position_jitter(width = 0.2, height = 0), 
             size = 3.5, shape = 21, stroke = 1.8) +  # Bold dots with larger size and stroke
  labs(title = "Proportion of exclusion vs. noexclusion Expression IN pol2 TIMECOURSE",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("exclusion" = exclusion, "noexclusion" = no_exclusion)) +  # Keep bar colors
  scale_color_manual(values = c("exclusion" = exclusion_dot_color, "noexclusion" = noexclusion_dot_color)) +  # Different dot colors
  theme_minimal() +
  theme(legend.position = "top")

p_pol2

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "pol2_TIMECOURSE_with_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p_pol2, height = 8, width = 8)



# Example data
data_pol2_GATA3 <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2'),
    Cell_type = c('D0', 'D4', 'D8', 'D10', 'D12', 'D14', 'D16', 'D0', 'D4', 'D8', 'D10', 'D12', 'D14', 'D16'),
  exclusion = c(NA, 24, NA, 57, 67, 82, 85, NA, 33.333, 67.3913, 68.80734, 65.78947, 81.15942, 78.04878),
  noexclusion = c(NA,76, NA, 43, 33, 18, 15, NA, 67.666, 32.6087, 31.19266, 34.21053, 18.84058, 21.95122)
    )
# Reshape data to long format
library(tidyr)
data_long_pol2_GATA3 <- data_pol2_GATA3 %>%
  pivot_longer(cols = c("exclusion", "noexclusion"), names_to = "ExpressionType", values_to = "Proportion")

#view(data_long_GAT3_POLA1)
#tail(data_long_GAT3_POLA1)
dim(data_long_pol2_GATA3)
head(data_long_pol2_GATA3, 50)

# Set factor levels for Cell_type in the desired order
data_long_pol2_GATA3$Cell_type <- factor(data_long_pol2_GATA3$Cell_type, levels = c("D0", "D4", "D8", "D10", "D12", "D14", "D16"))

# Load necessary libraries
library(ggplot2)

# Custom color palette for bars
exclusion <- "#94a096"
no_exclusion <- "#ce4751"

# New color palette for dots (distinct but corresponding)
exclusion_dot_color <- "#5b6d7c"  # Darker shade of bi_color
noexclusion_dot_color <- "#b03038"  # Darker shade of mono_color

# Plot with stacked bars and bold dots for individual experiments
p_pol2_gata3 <- ggplot(data_long_pol2_GATA3, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +  # Bars remain the same
  geom_point(aes(x = Cell_type, y = Proportion, color = ExpressionType),  # Different colors for dots
             position = position_jitter(width = 0.2, height = 0), 
             size = 3.5, shape = 21, stroke = 1.8) +  # Bold dots with larger size and stroke
  labs(title = "Proportion of exclusion vs. noexclusion Expression IN pol2 TIMECOURSE with GATA3 positive cells",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("exclusion" = exclusion, "noexclusion" = no_exclusion)) +  # Keep bar colors
  scale_color_manual(values = c("exclusion" = exclusion_dot_color, "noexclusion" = noexclusion_dot_color)) +  # Different dot colors
  theme_minimal() +
  theme(legend.position = "top")

p_pol2_gata3

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "pol2_TIMECOURSE_with_GATA3_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p_pol2_gata3, height = 8, width = 8)



# Example data
data_pol2_GATA6 <- data.frame(
     Experiment = c('EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP1', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2', 'EXP2'),
    Cell_type = c('D0', 'D4', 'D8', 'D10', 'D12', 'D14', 'D16', 'D0', 'D4', 'D8', 'D10', 'D12', 'D14', 'D16'),
  exclusion = c(NA, 31, NA, 71, 74, 82, 63, NA, 31.3725, 59.2593, 63.2653, 55.102, 66.6667, 80),
  noexclusion = c(NA,69, NA, 29, 26, 18, 37, NA, 68.6275, 40.7404, 36.7347, 44.898, 33.3333, 20)
    )
# Reshape data to long format
library(tidyr)
data_long_pol2_GATA6 <- data_pol2_GATA6 %>%
  pivot_longer(cols = c("exclusion", "noexclusion"), names_to = "ExpressionType", values_to = "Proportion")

#view(data_long_GAT3_POLA1)
#tail(data_long_GAT3_POLA1)
dim(data_long_pol2_GATA6)
head(data_long_pol2_GATA6, 50)

# Set factor levels for Cell_type in the desired order
data_long_pol2_GATA6$Cell_type <- factor(data_long_pol2_GATA6$Cell_type, levels = c("D0", "D4", "D8", "D10", "D12", "D14", "D16"))

# Load necessary libraries
library(ggplot2)

# Custom color palette for bars
exclusion <- "#94a096"
no_exclusion <- "#ce4751"

# New color palette for dots (distinct but corresponding)
exclusion_dot_color <- "#5b6d7c"  # Darker shade of bi_color
noexclusion_dot_color <- "#b03038"  # Darker shade of mono_color

# Plot with stacked bars and bold dots for individual experiments
p_pol2_gata6 <- ggplot(data_long_pol2_GATA6, aes(x = Cell_type, y = Proportion, fill = ExpressionType)) +
  geom_bar(stat = "summary", fun = "mean", position = "stack", width = 0.5) +  # Bars remain the same
  geom_point(aes(x = Cell_type, y = Proportion, color = ExpressionType),  # Different colors for dots
             position = position_jitter(width = 0.2, height = 0), 
             size = 3.5, shape = 21, stroke = 1.8) +  # Bold dots with larger size and stroke
  labs(title = "Proportion of exclusion vs. noexclusion Expression IN pol2 TIMECOURSE with GATA6 positive cells",
       x = "Cell Type",
       y = "Proportion (%)") +
  scale_fill_manual(values = c("exclusion" = exclusion, "noexclusion" = no_exclusion)) +  # Keep bar colors
  scale_color_manual(values = c("exclusion" = exclusion_dot_color, "noexclusion" = noexclusion_dot_color)) +  # Different dot colors
  theme_minimal() +
  theme(legend.position = "top")

p_pol2_gata6

output_path <- "/lustre1/project/stg_00041/Amitesh/FISH_quantification_random/"
filename <- "pol2_TIMECOURSE_with_GATA6_individual_experiments.pdf"

ggsave(file.path(output_path, filename), plot = p_pol2_gata6, height = 8, width = 8)

library(ggplot2)
library(dplyr)
library(tidyr)

library(readxl)
library(broom) 

data_pol2_mixed <- read.xlsx("/lustre1/project/stg_00041/Amitesh/STST_POL2_DYNAMICS/MIXED_POL2.xlsx")

data_pol2_mixed

alpha <- 0.01
data_long <- data_pol2_mixed %>%
  select(EXP, DAY, Exclusion, Noexclusion) %>%
  pivot_longer(cols = c(Exclusion, Noexclusion), 
               names_to = "Type", 
               values_to = "Count") %>%
  mutate(Expression_type = ifelse(Type == "Noexclusion", 1, 0))

head(data_long)

#Take D4, D10,D14
data_long$DAY <- factor(data_long$DAY)#relevel will only work on DAY if it is a factor
data_long$DAY <- relevel(data_long$DAY, ref = "D16")
Tc_model <- glm(Expression_type ~ DAY + factor(EXP),
                  data = data_long,
                  weights = Count,  # Using Count as weights
                  family = binomial)

summary(Tc_model) #returns results at alpha 0.05

Tc_model_summary <- summary(Tc_model) #and save to variable
p_values <- coef(Tc_model_summary)[, "Pr(>|z|)"]
print(p_values)

Tc_model_results <- broom::tidy(Tc_model) %>% 
  mutate(Significant = p.value < alpha)

print(Tc_model_results) #returns results at alpha 0.01

install.packages("emmeans")
library(emmeans)
pairwise_comparisons <- emmeans(Tc_model, pairwise ~ DAY, type = "response")
summary(pairwise_comparisons$contrasts)
summary(pairwise_comparisons$contrasts, adjust = "bonferroni")


output_dir <- "/lustre1/project/stg_00041/Amitesh/STST_POL2_DYNAMICS/"

write_xlsx(Tc_model_results, file.path(output_dir, "Tc_model_results_mixed_pol2_dynamics.xlsx"))

pairwise_df <- summary(pairwise_comparisons$contrasts, adjust = "bonferroni") %>% as.data.frame()
write_xlsx(pairwise_df, file.path(output_dir, "pairwise_comparisons_bonferroni_mixed_pol2_dynamics.xlsx"))

install.packages("writexl")


library(writexl)

data_pol2_GATA3 <- read.xlsx("/lustre1/project/stg_00041/Amitesh/STST_POL2_DYNAMICS/GATA3_POL2.xlsx")

data_pol2_GATA3

alpha <- 0.01
data_long_GATA3 <- data_pol2_GATA3 %>%
  select(EXP, DAY, Exclusion, Noexclusion) %>%
  pivot_longer(cols = c(Exclusion, Noexclusion), 
               names_to = "Type", 
               values_to = "Count") %>%
  mutate(Expression_type = ifelse(Type == "Noexclusion", 1, 0))

head(data_long_GATA3)

#Take D4, D10,D14
data_long_GATA3$DAY <- factor(data_long_GATA3$DAY)#relevel will only work on DAY if it is a factor
data_long_GATA3$DAY <- relevel(data_long_GATA3$DAY, ref = "D16")
Tc_model_GATA3 <- glm(Expression_type ~ DAY + factor(EXP),
                  data = data_long_GATA3,
                  weights = Count,  # Using Count as weights
                  family = binomial)

summary(Tc_model_GATA3) #returns results at alpha 0.05

Tc_model_summary_GATA3 <- summary(Tc_model_GATA3) #and save to variable
p_values_GATA3 <- coef(Tc_model_summary_GATA3)[, "Pr(>|z|)"]
print(p_values_GATA3)

Tc_model_results_GATA3 <- broom::tidy(Tc_model_GATA3) %>% 
  mutate(Significant = p.value < alpha)

print(Tc_model_results_GATA3) #returns results at alpha 0.01

library(emmeans)
pairwise_comparisons_GATA3 <- emmeans(Tc_model_GATA3, pairwise ~ DAY, type = "response")
summary(pairwise_comparisons_GATA3$contrasts)
summary(pairwise_comparisons_GATA3$contrasts, adjust = "bonferroni")

output_dir <- "/lustre1/project/stg_00041/Amitesh/STST_POL2_DYNAMICS/"

write_xlsx(Tc_model_results_GATA3, file.path(output_dir, "Tc_model_results_GATA3_pol2_dynamics.xlsx"))

pairwise_df <- summary(pairwise_comparisons_GATA3$contrasts, adjust = "bonferroni") %>% as.data.frame()
write_xlsx(pairwise_df, file.path(output_dir, "pairwise_comparisons_bonferroni_GATA3_pol2_dynamics.xlsx"))



data_pol2_GATA6 <- read.xlsx("/lustre1/project/stg_00041/Amitesh/STST_POL2_DYNAMICS/GATA6_POL2.xlsx")

data_pol2_GATA6

alpha <- 0.01
data_long_GATA6 <- data_pol2_GATA6 %>%
  select(EXP, DAY, Exclusion, Noexclusion) %>%
  pivot_longer(cols = c(Exclusion, Noexclusion), 
               names_to = "Type", 
               values_to = "Count") %>%
  mutate(Expression_type = ifelse(Type == "Noexclusion", 1, 0))

head(data_long_GATA6)

#Take D4, D10,D14
data_long_GATA6$DAY <- factor(data_long_GATA6$DAY)#relevel will only work on DAY if it is a factor
data_long_GATA6$DAY <- relevel(data_long_GATA6$DAY, ref = "D16")
Tc_model_GATA6 <- glm(Expression_type ~ DAY + factor(EXP),
                  data = data_long_GATA6,
                  weights = Count,  # Using Count as weights
                  family = binomial)

summary(Tc_model_GATA6) #returns results at alpha 0.05

Tc_model_summary_GATA6 <- summary(Tc_model_GATA6) #and save to variable
p_values_GATA6 <- coef(Tc_model_summary_GATA6)[, "Pr(>|z|)"]
print(p_values_GATA6)

Tc_model_results_GATA6 <- broom::tidy(Tc_model_GATA6) %>% 
  mutate(Significant = p.value < alpha)

print(Tc_model_results_GATA6) #returns results at alpha 0.01

library(emmeans)
pairwise_comparisons_GATA6 <- emmeans(Tc_model_GATA6, pairwise ~ DAY, type = "response")
summary(pairwise_comparisons_GATA6$contrasts)
summary(pairwise_comparisons_GATA6$contrasts, adjust = "bonferroni")

output_dir <- "/lustre1/project/stg_00041/Amitesh/STST_POL2_DYNAMICS/"

write_xlsx(Tc_model_results_GATA3, file.path(output_dir, "Tc_model_results_GATA6_pol2_dynamics.xlsx"))

pairwise_df <- summary(pairwise_comparisons_GATA6$contrasts, adjust = "bonferroni") %>% as.data.frame()
write_xlsx(pairwise_df, file.path(output_dir, "pairwise_comparisons_bonferroni_GATA6_pol2_dynamics.xlsx"))


