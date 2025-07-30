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
install.packages("viridis")
install.packages("viridisLite")
library(pheatmap)
library(viridis)
library(readxl)

# Load your data
gene_function_distribution <- read_excel("/media/sarlab/DATA/Bacillus_project/Bacillus_project_gene_mining/Supplementary_Material_Functions_Distribution_per_Isolate.xlsx")  # Replace with your filename

# Set row names 
gene_function_distribution <- as.data.frame(gene_function_distribution)
rownames(gene_function_distribution) <- gene_function_distribution$Isolate_ID  
gene_function_distribution <- gene_function_distribution[ , -1]

range_val <- max(abs(gene_function_distribution - 80))  # Calculate max deviation from 80
breaks <- seq(80 - range_val, 80 + range_val, length.out = 100)

# Create blue-white-red diverging color palette
palette <- colorRampPalette(c("blue", "white", "red"))(length(breaks))

# Plot heatmap
pheatmap(
  gene_function_distribution,
  color = palette,
  breaks = breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.0f",
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Functional Gene Counts per Isolate",
  border_color = "grey90",
  angle_col = 45
)
pheatmap(
  gene_function_distribution,
  color = colorRampPalette(viridis(100))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.0f",
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Functional Gene Counts per Isolate",
  border_color = "grey90",
  angle_col = 45
)
