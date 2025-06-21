######################################################################################################
# script name: bgcs_boxplot_pergenus.R
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
######################################################################################################
# GOAL:
# Aim of this script is to generate a boxplot of BGC counts per genus. 
######################################################################################################
# usage:./bgcs_boxplot_pergenus.R
# complete path: /home/nik_arapitsas/Documents/Bacillus_project/scripts/bgcs_boxplot_pergenus.R
######################################################################################################

# Load the necessary libraries 
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)

# Read the CSV file that contains the Isolate ID, BGC type, BGC count and Similarity Confidence without the species name 
bgcs_perisolate <- read_csv("/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/bgc_similarity_with_species.csv")

# Keep only the Genus name for every isolate 
bgcs_perisolate_genera <- bgcs_perisolate %>%
  mutate(genus = word(Species, 1))  # Extracts the first word (genus)

# Count BGCs per isolate, group them, and dsplay just the genus name of every isolate
bgc_perisolate_genera_grouped <- bgcs_perisolate_genera %>%
  group_by(IsolateID, genus) %>%
  summarise(bgc_count = sum(BGC_Count), .groups = "drop")

# Count the number of isolates in each genus
genus_sample_sizes <- bgc_perisolate_genera_grouped %>%
  count(genus) %>%
  rename(n = n)


# Define custom colors

genus_colors <- c(
  "Bacillus" = "#7171be",         
  "Paenibacillus" = "#d4c63a",   
  "Neobacillus" = "#ff94b4",      
  "Peribacillus" = "#a26324",     
  "Cytobacillus" = "#8f9ed7",
  "Rossellomorea" = "#f8d48c"     
)

# Create the graph

bgc_perisolate_genera_grouped_graph <- ggplot(bgc_perisolate_genera_grouped, aes(x = bgc_count, y = reorder(genus, bgc_count, median), fill = genus)) +
  geom_boxplot(outlier.shape = 18) +
  scale_fill_manual(values = genus_colors) +
  geom_text(data = genus_sample_sizes,
          aes(x = 2, y = genus, label = paste0("n: ", n)),
          hjust = 0, vjust = 0.5,
          inherit.aes = FALSE, size = 3.5) +
  theme_minimal() +
  scale_y_discrete(labels = function(x) parse(text = paste0("italic('", x, "')"))) +
  labs(
    title = "Distribution of antiSMASH regions by genus",
    x = "antiSMASH regions (count)",
    y = NULL
  ) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
    axis.ticks.length = unit(0.1, 'cm'),  # Set the tick length to be smaller
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center align title
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10, color = "black", family = "sans", margin = margin(t = 2.5)),
    axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black", family = "sans", margin = margin(t = 2.0)),
    plot.margin = margin(10, 10, 10, 10)
) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  coord_cartesian(xlim = c(0, 20))

# Save the graph 

ggsave(paste0("/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/bgc_perisolate_genera_grouped_graph",".png"),
       plot= bgc_perisolate_genera_grouped_graph, 
       height = 20, 
       width = 50,
       dpi = 300, 
       units="cm",
       device="png")