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
ggplot(all_pangenome_data, aes(x = group, y = percent, fill = category)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  geom_text(aes(label = paste0(percent, "%")),
            position = position_stack(vjust = 0.5),
            color = "white", size = 4.5, fontface = "bold") +
  scale_fill_manual(
    values = c(
      "core" = "#7171be",
      "accessory" = "#ff7f0e",
      "species-specific" = "#2ca02c"
    )
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
  labs(
    x = "Pangenomes",
    y = "Gene Percentage (%)",
    fill = "Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
    axis.ticks.length = unit(0.1, 'cm'),  # Set the tick length to be smaller
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center align title
    axis.text.x = element_text(angle = 30, hjust = 1, 
                               size = 10, color = "black", family = "sans", 
                               margin = margin(t = 2.5)),  # Darker and more spaced letters
    axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black", family = "sans", 
                               margin = margin(t = 2.0)),
    plot.margin = margin(10, 10, 10, 10)  # Add extra space around the plot
  )
