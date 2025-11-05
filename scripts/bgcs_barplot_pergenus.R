######################################################################################################
# script name: bgcs_barplot_pergenus.R
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
######################################################################################################
# GOAL:
# Aim of this script is to generate a barplot of BGC counts per genus.
######################################################################################################
# usage:./bgcs_barplot_pergenus.R
# complete path: /home/nik_arapitsas/Documents/Bacillus_project/scripts/bgcs_barplot_pergenus.R
######################################################################################################

# Load the necessary libraries
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(forcats)

# Read the CSV file that contains the Isolate ID, BGC type, BGC count and Similarity Confidence without the species name 
bgcs_perisolate <- read_csv("bgc_similarity_with_species.csv")

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
  "Bacillus" = "#009E73",        
  "Paenibacillus" = "#D55E00",   
  "Neobacillus" = "#F0E442",     
  "Peribacillus" = "#0072B2",    
  "Cytobacillus" = "#56B4E9",    
  "Rossellomorea" = "#E69F00"    
)

# Create a barplot

# Summarize BGCs per genus
bgc_per_genus <- bgcs_perisolate_genera %>%
  group_by(genus) %>%
  summarise(total_bgcs = sum(BGC_Count), .groups = "drop") %>%
  mutate(genus = fct_reorder(genus, total_bgcs))  # order by BGCs

# Summarize mean and SE per genus
bgc_stats <- bgc_perisolate_genera_grouped %>%
  group_by(genus) %>%
  summarise(
    mean_bgc = mean(bgc_count),
    se_bgc = sd(bgc_count) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(genus = fct_reorder(genus, mean_bgc))  # order for plotting

# Create the barplot
bgc_perisolate_genera_grouped_barplot <- ggplot(bgc_stats, aes(x = genus, y = mean_bgc, fill = genus)) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_bgc - se_bgc, ymax = mean_bgc + se_bgc),
    width = 0.2 
  ) +
  geom_text(
  aes(
    label = paste0("n = ", n),
    y = ifelse(is.na(se_bgc) | is.nan(se_bgc), mean_bgc + 0.5, mean_bgc + se_bgc + 0.5)
  ),
  size = 3,
  color = "black",
  fontface = "bold"
) +
  scale_fill_manual(values = genus_colors) +  # your custom palette
  labs(
    title = "Average number of antiSMASH regions per genus",
    x = "Genus",
    y = "Mean antiSMASH region count ± SE"
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
        axis.ticks.length = unit(0.1, 'cm'),  # Set the tick length to be smaller
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center align title
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, 
                                   size = 10, color = "black", family = "sans", 
                                   margin = margin(t = 2.5), face = "italic"),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black", family = "sans", 
                                   margin = margin(t = 2.0)),
        legend.text = element_text(size = 10),  # Reduce font size of legend
        legend.spacing.y = unit(0.2, 'cm'),
        legend.key.size = unit(0.5, "cm"),  # Adjust size of legend keys (squares)
        plot.margin = margin(10, 10, 10, 10),
        legend.position = "none"
  ) +  
  scale_y_continuous(breaks = 0:16, expand = expansion(mult = c(0, 0.01))) +
  coord_cartesian(ylim = c(0, 16))

# Save the graph 

ggsave(paste0("bgc_perisolate_genera_grouped_barplot",".png"),
       plot= bgc_perisolate_genera_grouped_barplot, 
       height = 15, 
       width = 20,
       dpi = 300, 
       units="cm",
       device="png")
