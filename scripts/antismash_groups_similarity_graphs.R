###########################################################################################################
# script name: antismash_groups_similarity_graphs.R
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
###########################################################################################################
# GOAL:
# Aim of this script is to plot the BGC count per isolate based on Similarity Confidence and a collective 
# plot that depicts the distribution of BGCs across the isolates based both on their category and their 
# similarity confidence to known BGCs
###########################################################################################################
# usage:./antismash_groups_similarity_graphs.R
# complete path: /home/nik_arapitsas/Documents/Bacillus_project/scripts/antismash_groups_similarity_graphs.R
###########################################################################################################

# Load the necessary libraries 
library(readr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)

# 1). Create the Plot that depicts the BGC count based on Similarity Confidence per isolate (Supplementary Fig. S10)

# Read the CSV file that contains the Isolate ID, and BGC count based on Similarity Confidence
antiSMASH_regions_similarity_count <- read_csv("antismash_table.csv")

# Rename columns to make them easier to work with
antiSMASH_regions_similarity_count <- antiSMASH_regions_similarity_count %>%
  rename(IsolateID = `Isolate ID`,
         Low = `Low Similarity`,
         Medium = `Medium Similarity`,
         High = `High Similarity`,
         Undefined = `Undefined Similarity`)

# Pivot to long format for ggplot
antiSMASH_regions_similarity_count_long <- antiSMASH_regions_similarity_count %>%
  select(IsolateID, High, Medium, Low, Undefined) %>%
  pivot_longer(cols = High:Undefined, names_to = "Similarity", values_to = "Count")

# Set stacking order
antiSMASH_regions_similarity_count_long$Similarity <- factor(antiSMASH_regions_similarity_count_long$Similarity, levels = c("Undefined", "Low", "Medium", "High"))

# Define custom colors
similarity_colors <- c(
  "High" = "#009E73",
  "Medium" = "#56B4E9",
  "Low" = "#ff94b4",
  "Undefined" = "#7171be"
)

# Plot
antiSMASH_regions_similarity_count_plot <- ggplot(antiSMASH_regions_similarity_count_long, aes(x = IsolateID, y = Count, fill = Similarity)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = similarity_colors, name = "Similarity Confidence") +
  scale_x_discrete(expand=expansion(add=c(.8, .8))) +
  scale_y_continuous(expand=expansion(add=c(0, 0))) +
  labs(title = "antiSMASH region count per isolate",
       x = NULL,
       y = "Number of antiSMASH Regions") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
        panel.border = element_blank(),
        axis.ticks.length = unit(0.1, 'cm'),  # Set the tick length to be smaller
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center align title
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, 
                                   size = 10, color = "black", family = "sans", 
                                   margin = margin(t = 2.5)),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black", family = "sans", 
                                   margin = margin(t = 2.0)),
        legend.text = element_text(size = 10),  # Reduce font size of legend
        legend.spacing.y = unit(0.2, 'cm'),
        legend.key.size = unit(0.5, "cm"),  # Adjust size of legend keys (squares)
        plot.margin = margin(10, 10, 10, 10)) +  # Add extra space around the plot
          coord_cartesian(ylim = c(0,20))

ggsave(paste0("antiSMASH_regions_similarity_count_plot",".png"),
       plot=antiSMASH_regions_similarity_count_plot, 
       height = 22, 
       width = 55,
       dpi = 300, 
       units="cm",
       device="png")

# 2). Create the Collective plot that depicts the distribution of BGCs across the isolates based both on their category and their similarity confidence to known BGCs (Figure 6)

# Read the CSV file that contains the Isolate ID, BGC type, BGC count and Similarity Confidence 
bgcs_perisolate <- read_csv("bgc_similarity_with_species.csv")

bgcs_perisolate_bgc <- bgcs_perisolate |>
    mutate(Category = case_when(
        BGC_Type %in% c("NRPS", "NRPS-like") ~ "NRPS",
        BGC_Type %in% c("terpene", "terpene-precursor") ~ "terpene",
        BGC_Type %in% c("T3PKS", "transAT-PKS", "PKS-like", "HR-T2PKS") ~ "PKS",
        BGC_Type == "NI-siderophore" ~ "NI-siderophore",
        grepl("NRPS.*PKS|PKS.*NRPS", BGC_Type) ~ "NRPS-PKS hybrids",
        BGC_Type %in% c("RiPP-like", "azole-containing-RiPP") ~ "RiPPs",
        BGC_Type %in% c("NRP-metallophore.NRPS.RiPP-like.terpene-precursor",
                    "NRP-metallophore.NRPS", 
                    "NRPS.RRE-containing",
                    "NRPS.terpene",
                    "NRPS.betalactone",
                    "CDPS.NRPS") ~ "NRPS-other hybrids",
        TRUE ~ "others"
        ))

# Define custom colors

bgc_group_colors <- c(
  "NRPS" = "#0072B2",
  "terpene" = "#E69F00",
  "PKS" = "#D55E00",
  "NI-siderophore" = "#009E73",
  "NRPS-PKS hybrids" = "#F0E442",
  "RiPPs" = "#CC79A7",
  "NRPS-other hybrids" = "#56B4E9",
  "others" = "#999999"
)

# Connect all the similarity categories in a common graph

bgcs_perisolate_summary <- bgcs_perisolate_bgc |>
    group_by(IsolateID,Similarity,Category) |>
    summarise(Count=n(), .groups = "keep")

bgcs_complete <- bgcs_perisolate_bgc |>
    group_by(IsolateID,Similarity,Category) |>
    summarise(Count=n(), .groups = "keep") %>%
    ungroup() %>%
    complete(IsolateID, Similarity, Category, fill = list(Count = 0)) %>%
    mutate(Similarity = factor(Similarity,
                               levels = c("High", "Medium", "Low", "Undefined"),
                               labels = c("High similarity confidence",
                                          "Medium similarity confidence",
                                          "Low similarity confidence",
                                          "Undefined similarity confidence")))

bgc_types_plot <- ggplot(bgcs_complete, aes(x = IsolateID, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = bgc_group_colors) +
  theme_bw() +
  labs(title = "BGCs per Isolate by Category",
       x = "Isolate ID",
       y = "Count of BGCs",
       fill = "BGC Category") +
  theme(panel.grid = element_blank(), #remove grid lines
        axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
        axis.ticks.length = unit(0.1, 'cm'),  # Set the tick length to be smaller
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center align title
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, 
                                   size = 10, color = "black", family = "sans", 
                                   margin = margin(t = 2.5)),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black", family = "sans", 
                                   margin = margin(t = 2.0)),
        legend.text = element_text(size = 10),  # Reduce font size of legend
        legend.spacing.y = unit(0.2, 'cm'),
        strip.placement = "outside",
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        legend.key.size = unit(0.5, "cm"),  # Adjust size of legend keys (squares)
        plot.margin = margin(10, 10, 10, 10)) +  # Add extra space around the plot
        scale_y_continuous(breaks = 0:12, expand = expansion(mult = c(0, 0.01))) +
        coord_cartesian(ylim = c(0, 12)) + 
        facet_wrap(~Similarity, ncol = 1, strip.position = "top", scales = "free")

# Save the graph 

ggsave(paste0("bgc_types_plot_colorblind",".png"),
       plot= bgc_types_plot, 
       height = 40, 
       width = 40,
       dpi = 300, 
       units="cm",
       device="png")
