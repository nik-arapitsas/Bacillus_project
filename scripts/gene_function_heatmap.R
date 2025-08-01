######################################################################################################
# script name: gene_function_heatmap.R
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
######################################################################################################
# GOAL:
# Aim of this script is to create a heatmap illustrating the distribution of gene functions across the
# 25 isolates.
######################################################################################################
# usage:./gene_function_heatmap.R
# complete path: /home/nik_arapitsas/Documents/Bacillus_project/scripts/gene_function_heatmap.R
######################################################################################################

# Load required packages
install.packages("pheatmap")
library(pheatmap)
library(readxl)

setwd("/media/sarlab/DATA/Bacillus_project/Bacillus_project_gene_mining")

# Load your data
gene_function_distribution <- read_excel("Supplementary_Material_Functions_Distribution_per_Isolate.xlsx")  

# Set row names 
gene_function_distribution <- as.data.frame(gene_function_distribution)
rownames(gene_function_distribution) <- gene_function_distribution$Isolate_ID  
gene_function_distribution <- gene_function_distribution[ , -1]

min_val <- min(gene_function_distribution)
max_val <- max(gene_function_distribution)

# Set two midpoints

# Ensure increasing, valid breaks
b1 <- max(min_val, 20)
b2 <- max(b1 + 0.01, 40)
b3 <- max(b2 + 0.01, 80)
b4 <- max(b3 + 0.01, max_val + 1)

# Create 4 color segments (3 breakpoints = 4 color bins)
n_segments <- 40  
breaks <- c(
  seq(min_val, b1, length.out = n_segments),
  seq(b1 + 0.01, b2, length.out = n_segments),
  seq(b2 + 0.01, b3, length.out = n_segments),
  seq(b3 + 0.01, b4, length.out = n_segments)
)

# Make sure breaks are unique
breaks <- unique(breaks)

# Set the palette
colors <- c(
  colorRampPalette(c("#3B84CE", "#67a9cf"))(n_segments - 1),         
  colorRampPalette(c("#67a9cf", "#fddbc7"))(n_segments - 1),         
  colorRampPalette(c("#fddbc7", "#f4a582"))(n_segments - 1),         
  colorRampPalette(c("#f4a582", "#b2182b"))(length(breaks) - 3 * (n_segments - 1))  
)

colnames(gene_function_distribution) <- gsub("_", " ", colnames(gene_function_distribution))
colnames(gene_function_distribution) <- gsub(
  "Plant Hormone & VOCs production",
  "Plant Hormone\n& VOCs production",
  colnames(gene_function_distribution)
)

png("gene_function_heatmap.png", width = 3000, height = 3000, res = 300)

pheatmap(
  gene_function_distribution,
  color = colors,
  breaks = breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.0f",
  border_color = "black",   
  number_color = "black",
  fontsize_number = 12,
  fontsize_row = 12, 
  fontsize_col = 12,
  angle_col = 45,  
  main = "Gene Counts per Isolate"
)

dev.off()
