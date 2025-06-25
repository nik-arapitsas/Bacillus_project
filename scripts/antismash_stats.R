######################################################################################################
# script name: antismash_stats.R
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
######################################################################################################
# GOAL:
# Aim of this script is to help providing some stats for the antiSMASH output
######################################################################################################
# usage:./antismash_stats.R
# complete path: /home/nik_arapitsas/Documents/Bacillus_project/scripts/antismash_stats.R
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
                    "NRPS.betalactone",
                    "CDPS.NRPS") ~ "NRPS-other hybrids",
        TRUE ~ "others"
        )) 

# Select the regions with Undefined similarity

bgcs_undefined <- bgcs_perisolate_bgc %>%
  filter(Similarity == "Undefined")

# Summarize Undefined BGCs and categories

bgcs_summary_undefined <- bgcs_undefined %>%
  group_by(BGC_Type) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  arrange(desc(Count))

bgcs_perisolate_bgc_summary_undefined <- bgcs_undefined %>%
  group_by(Category) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  arrange(desc(Count)) 

# Select the regions with Undefined similarity

bgcs_low <- bgcs_perisolate_bgc %>%
  filter(Similarity == "Low")

# Summarize Low BGCs and categories

bgcs_summary_low <- bgcs_low %>%
  group_by(BGC_Type) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  arrange(desc(Count))

bgcs_perisolate_bgc_summary_low <- bgcs_low %>%
  group_by(Category) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  arrange(desc(Count)) 

# Total BGCs

bgcs_undefined <- bgcs_perisolate_bgc %>%
  filter(Similarity == "Undefined")

# Summarize Undefined BGCs and categories

bgcs_summary <- bgcs_perisolate %>%
  group_by(BGC_Type) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  arrange(desc(Count))

bgcs_perisolate_bgc_summary <- bgcs_perisolate_bgc %>%
  group_by(Category) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Percentage = 100 * Count / sum(Count)) %>%
  arrange(desc(Count)) 

# Summarize High BGCs and categories

bgcs_high <- bgcs_perisolate_bgc %>%
  filter(Similarity == "High")

bgcs_summary_high <- bgcs_high %>%
  group_by(BGC_Type) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  arrange(desc(Count))

bgcs_perisolate_bgc_summary_high <- bgcs_high %>%
  group_by(Category) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Percentage = 100 * Count / sum(Count)) %>%
  arrange(desc(Count)) 

# Summarize Medium BGCs and categories

bgcs_medium <- bgcs_perisolate_bgc %>%
  filter(Similarity == "Medium")

bgcs_summary_medium <- bgcs_medium %>%
  group_by(BGC_Type) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  arrange(desc(Count))

bgcs_perisolate_bgc_summary_medium <- bgcs_medium %>%
  group_by(Category) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Percentage = 100 * Count / sum(Count)) %>%
  arrange(desc(Count)) 

# Create a wide table

bgc_wide <- bgcs_perisolate_bgc %>%
  group_by(IsolateID, Category) %>%
  summarise(BGC_Count = sum(BGC_Count), .groups = "drop") %>%
  pivot_wider(
    names_from = Category, 
    values_from = BGC_Count,
    values_fill = 0
  )

# Save the wide table 
write.csv(bgc_wide, "/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/bgc_wide_table.csv", row.names = FALSE)

