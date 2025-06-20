######################################################################################################
# script name: antismash_groups_graph.R
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
######################################################################################################
# GOAL:
# Aim of this script is to create graphs plotting the types of BGC per isolate containing also info 
# about their similarity with kown clusters
######################################################################################################
# usage:./antismash_groups_graph.R
# complete path: /home/nik_arapitsas/Documents/Bacillus_project/scripts/antismash_groups_graph.R
######################################################################################################

# Load the necessary libraries 
library(readr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the CSV file that contains the Isolate ID, BGC type, BGC count and Similarity Confidence without the species name 
bgcs_perisolate <- read_csv("/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/bgc_similarity_with_species.csv")

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
                    "NI-siderophore.terpene",
                    "NRPS.betalactone") ~ "complex",
        TRUE ~ "others"
        )) 

# A) Depict only the BGCs with Undefined similarity

# Select the regions with Undefined similarity

bgcs_undefined <- bgcs_perisolate_bgc %>%
  filter(Similarity == "Undefined")

# Define the BGC groups

# Summarize Undefined BGCs per isolate and category

bgcs_summary <- bgcs_undefined %>%
  group_by(IsolateID, Category) %>%
  summarise(Count = n(), .groups = 'drop')

# Define custom colors

bgc_group_colors <- c(
  "NRPS" = "#7171be",
  "terpene" = "#d4c63a",
  "PKS" = "#ff94b4",
  "NI-siderophore" = "#a26324",
  "NRPS-PKS hybrids" = "#8f9ed7",
  "RiPPs" = "#c18563",
  "complex" = "#c47ece",
  "others" = "#f8d48c"
)

# Plot stacked bar graph

undefined_bgc_types_plot <- ggplot(bgcs_summary, aes(x = IsolateID, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = bgc_group_colors) +
  theme_minimal() +
  labs(title = "Undefined BGCs per Isolate by Category",
       x = "Isolate ID",
       y = "Count of Undefined BGCs",
       fill = "BGC Category") +
  theme(axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
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
        scale_y_continuous(breaks = 0:11, expand = expansion(mult = c(0, 0.01))) +
        coord_cartesian(ylim = c(0, 11))

# Save the graph 

ggsave(paste0("/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/undefined_bgc_types_plot_colorblind",".png"),
       plot= undefined_bgc_types_plot, 
       height = 20, 
       width = 50,
       dpi = 300, 
       units="cm",
       device="png")

# B) Depict only the BGCs with Low similarity

# Select the regions with Low similarity

bgcs_low <- bgcs_perisolate_bgc %>%
  filter(Similarity == "Low")

# Define the BGC groups


# Summarize Undefined BGCs per isolate and category

bgcs_low_summary <- bgcs_low %>%
  group_by(IsolateID, Category) %>%
  summarise(Count = n(), .groups = 'drop')

# Plot stacked bar graph

low_bgc_types_plot <- ggplot(bgcs_low_summary, aes(x = IsolateID, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = bgc_group_colors) +
  theme_minimal() +
  labs(title = "Low BGCs per Isolate by Category",
       x = "Isolate ID",
       y = "Count of Low BGCs",
       fill = "BGC Category") +
  theme(axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
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
        scale_y_continuous(breaks = 0:11, expand = expansion(mult = c(0, 0.01))) +
        coord_cartesian(ylim = c(0, 5))

# Save the graph 

ggsave(paste0("/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/low_bgc_types_plot_colorblind",".png"),
       plot= low_bgc_types_plot, 
       height = 20, 
       width = 50,
       dpi = 300, 
       units="cm",
       device="png")

# C) Depict only the BGCs with Low similarity

# Select the regions with Medium similarity

bgcs_medium <- bgcs_perisolate_bgc %>%
  filter(Similarity == "Medium")

# Define the BGC groups

# Summarize Undefined BGCs per isolate and category

bgcs_medium_summary <- bgcs_medium %>%
  group_by(IsolateID, Category) %>%
  summarise(Count = n(), .groups = 'drop')

# Plot stacked bar graph

medium_bgc_types_plot <- ggplot(bgcs_medium_summary, aes(x = IsolateID, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = bgc_group_colors) +
  theme_minimal() +
  labs(title = "Medium BGCs per Isolate by Category",
       x = "Isolate ID",
       y = "Count of Medium BGCs",
       fill = "BGC Category") +
  theme(axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
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
        scale_y_continuous(breaks = 0:11, expand = expansion(mult = c(0, 0.01))) +
        coord_cartesian(ylim = c(0, 5))

# Save the graph 

ggsave(paste0("/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/medium_bgc_types_plot_colorblind",".png"),
       plot= medium_bgc_types_plot, 
       height = 20, 
       width = 50,
       dpi = 300, 
       units="cm",
       device="png")

# D) Depict only the BGCs with High similarity

# Select the regions with High similarity

bgcs_high <- bgcs_perisolate_bgc %>%
  filter(Similarity == "High")

# Define the BGC groups


# Summarize Undefined BGCs per isolate and category

bgcs_high_summary <- bgcs_high %>%
  group_by(IsolateID, Category) %>%
  summarise(Count = n(), .groups = 'drop')

# Plot stacked bar graph

high_bgc_types_plot <- ggplot(bgcs_high_summary, aes(x = IsolateID, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = bgc_group_colors) +
  theme_minimal() +
  labs(title = "Heigh BGCs per Isolate by Category",
       x = "Isolate ID",
       y = "Count of High BGCs",
       fill = "BGC Category") +
  theme(axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
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
        scale_y_continuous(breaks = 0:11, expand = expansion(mult = c(0, 0.01))) +
        coord_cartesian(ylim = c(0, 8))

# Save the graph 

ggsave(paste0("/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/high_bgc_types_plot_colorblind",".png"),
       plot= high_bgc_types_plot, 
       height = 20, 
       width = 50,
       dpi = 300, 
       units="cm",
       device="png")

# D) Connect all the similarity categories in a common graph

bgcs_perisolate_summary <- bgcs_perisolate_bgc |>
    group_by(IsolateID,Similarity,Category) |>
    summarise(Count=n(), .groups = "keep")

bgcs_complete <- bgcs_perisolate_summary %>%
  ungroup() %>%
  complete(IsolateID, Similarity, Category, fill = list(Count = 0))

bgc_types_plot <- ggplot(bgcs_complete, aes(x = IsolateID, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = bgc_group_colors) +
  theme_minimal() +
  labs(title = "BGCs per Isolate by Category",
       x = "Isolate ID",
       y = "Count of BGCs",
       fill = "BGC Category") +
  theme(axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
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

ggsave(paste0("/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/bgc_types_plot_colorblind",".png"),
       plot= bgc_types_plot, 
       height = 40, 
       width = 40,
       dpi = 300, 
       units="cm",
       device="png")
