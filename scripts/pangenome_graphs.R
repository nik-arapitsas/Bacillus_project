######################################################################################################
# script name: pangenome_graphs.R
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
######################################################################################################
# GOAL:
# Aim of this script is to create the collective barplot and upside plots of the pangenome analysis
######################################################################################################
# usage:./pangenome_graphs.R
# complete path: /home/nik_arapitsas/Documents/Bacillus_project/scripts/pangenome_graphs.R
######################################################################################################


# 1) Collective Barplot 

library(tidyverse)

# List of your actual .txt files with labels
setwd("/media/sarlab/DATA/Bacillus_project/Bacillus_project_anvio/Bacillus_project_anvio_Graphs")

pangenome_tables <- list(
  "SRL179 and relatives" = "SRL179_gene_cluster_output.txt",
  "SRL337 and relatives" = "SRL337_gene_cluster_output.txt",
  "SRL543 and relatives" = "SRL543_gene_cluster_output.txt",
  "SRL368 and relatives" = "SRL368_gene_cluster_output.txt"
)

# Function to process each table
process_table <- function(file, label) {
  df <- read_delim(file, delim = "\t", col_types = cols())

  n_genomes <- length(unique(df$genome_name))

  df_summary <- df %>%
    distinct(gene_cluster_id, genome_name) %>%
    count(gene_cluster_id) %>%
    mutate(category = case_when(
      n == n_genomes ~ "core",
      n == 1 ~ "species_specific",
      TRUE ~ "accessory"
    )) %>%
    count(category) %>%
    mutate(
      total = sum(n),
      percent = round(n / total * 100, 1),
      group = label
    )

  return(df_summary %>% select(group, category, n, percent))
}

# Process and bind all tables
all_pangenome_data <- bind_rows(
  lapply(names(pangenome_tables), function(name) {
    process_table(pangenome_tables[[name]], name)
  })
)

# Set factor levels for plotting order
all_pangenome_data$category <- recode(all_pangenome_data$category,
                            species_specific = "species-specific")
all_pangenome_data$category <- factor(all_pangenome_data$category, levels = c("core", "accessory", "species-specific"))
all_pangenome_data$group <- factor(all_pangenome_data$group, levels = names(pangenome_tables))

# Plot
pangenome_barplot <- ggplot(all_pangenome_data, aes(x = percent, y = group, fill = category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  geom_text(aes(label = paste0(percent, "%")),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4.5, fontface = "bold") +
  scale_fill_manual(
    values = c(
      "core" = "#56B4E9",           # sky blue (Okabe-Ito)
  "accessory" = "#CC79A7",      # reddish purple (Okabe-Ito)
  "species-specific" = "#E69F00" # orange (Okabe-Ito)
    )) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(
    x = "Pangenome Composition (%)",
    y = "Pangenomes",
    fill = "Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black", size = 0.5),  # Add black axis lines
    axis.ticks.length = unit(0.1, 'cm'),  # Set the tick length to be smaller
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center align title
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = -1, 
                               size = 12, color = "black", family = "sans", 
                               margin = margin(t = 2.5)),  # Darker and more spaced letters
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 12, color = "black", family = "sans", 
                               margin = margin(t = 2.0)),
    plot.margin = margin(10, 10, 10, 10)  # Add extra space around the plot
  )

# Save the graph 

ggsave(paste0("/media/sarlab/DATA/Bacillus_project/Bacillus_project_anvio/Bacillus_project_anvio_Graphs/pangenome_barplot",".png"),
       plot= pangenome_barplot, 
       height = 20, 
       width = 25,
       dpi = 300, 
       units="cm",
       device="png")


# 2) Upset Plots 

install.packages("UpSetR")
library(UpSetR)

# i) SRL179 

# Load the table
srl179 <- read_delim("SRL179_gene_cluster_output.txt", delim = "\t", col_types = cols())

# Create presence/absence matrix
srl179_matrix <- srl179 %>%
  distinct(gene_cluster_id, genome_name) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = genome_name, values_from = present, values_fill = 0) %>%
  column_to_rownames("gene_cluster_id")

# UpSet plot
png("SRL179_upset.png", width = 4300, height = 3440, res = 300)
upset(srl179_matrix,
      sets = colnames(srl179_matrix),
      order.by = "freq",
      mainbar.y.label = "\n\n\n\n\n\n\n\nGene Cluster Intersections - SRL179 and relatives",
      sets.x.label = "Gene Clusters per Genome",
      text.scale = c(1.3, 1.5, 1.5, 1.5, 1.5, 1.3))
dev.off()

# ii) SRL337 

# Load the table
srl337 <- read_delim("SRL337_gene_cluster_output.txt", delim = "\t", col_types = cols())

# Create presence/absence matrix
srl337_matrix <- srl337 %>%
  distinct(gene_cluster_id, genome_name) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = genome_name, values_from = present, values_fill = 0) %>%
  column_to_rownames("gene_cluster_id")

# UpSet plot
png("SRL337_upset.png", width = 4300, height = 3440, res = 300)
upset(srl337_matrix,
      sets = colnames(srl337_matrix),
      order.by = "freq",
      mainbar.y.label = "\n\n\n\n\n\n\n\nGene Cluster Intersections - SRL337 and relatives",
      sets.x.label = "Gene Clusters per Genome",
      text.scale = c(1.3, 1.5, 1.5, 1.5, 1.5, 1.3))
dev.off()

# iii) SRL543

# Load the table
srl543 <- read_delim("SRL543_gene_cluster_output.txt", delim = "\t", col_types = cols())

# Create presence/absence matrix
srl543_matrix <- srl543 %>%
  distinct(gene_cluster_id, genome_name) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = genome_name, values_from = present, values_fill = 0) %>%
  column_to_rownames("gene_cluster_id")

# UpSet plot
png("SRL543_upset.png", width = 4300, height = 3440, res = 300)
upset(srl543_matrix,
      sets = colnames(srl543_matrix),
      order.by = "freq",
      mainbar.y.label = "\n\n\n\n\n\n\n\nGene Cluster Intersections - SRL543 and relatives",
      sets.x.label = "Gene Clusters per Genome",
      text.scale = c(1.3, 1.5, 1.5, 1.5, 1.5, 1.3))
dev.off()

# iv) SRL368

# Load the table
srl368 <- read_delim("SRL368_gene_cluster_output.txt", delim = "\t", col_types = cols())

# Create presence/absence matrix
srl368_matrix <- srl368 %>%
  distinct(gene_cluster_id, genome_name) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = genome_name, values_from = present, values_fill = 0) %>%
  column_to_rownames("gene_cluster_id")

# UpSet plot
png("SRL368_upset.png", width = 4300, height = 3440, res = 300)
upset(srl368_matrix,
      sets = colnames(srl368_matrix),
      order.by = "freq",
      mainbar.y.label = "\n\n\n\n\n\n\n\nGene Cluster Intersections - SRL368 and relatives",
      sets.x.label = "Gene Clusters per Genome",
      text.scale = c(1.3, 1.5, 1.5, 1.5, 1.5, 1.3))
dev.off() 