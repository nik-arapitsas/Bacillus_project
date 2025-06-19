# Load necessary libraries
library(tidyverse)

# Read the CSV file (adjust the path if needed)
df <- read_csv2("/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/antismash_table_bgc_groups.csv")

# Rename columns to make them easier to work with
df <- df %>%
  rename(IsolateID = `Isolate ID`)

similarity_cols <- c("Low", "Medium", "High", "Undefined")

bgc_cols <- c(
  "azole-containing-RiPP", "NI-siderophore", "lassopeptide", "terpene", "RiPP-like", "betalactone", 
  "CDPS.NRPS", "HR-T2PKS", "NRPS.terpene", "NRP-metallophore.NRPS", "NRPS.RRE-containing", 
  "NRPS.NRPS-like.T3PKS.transAT-PKS", "NRP-metallophore.NRPS.RiPP-like.terpene-precursor", 
  "NRPS", "T3PKS", "opine-like-metallophore", "lanthipeptide-class-i", 
  "cyclic-lactone-autoinducer.lanthipeptide-class-ii", "NRPS.T1PKS.transAT-PKS", 
  "NRPS.transAT-PKS", "NRPS-like", "proteusin", "lanthipeptide-class-ii", "RRE-containing", 
  "NI-siderophore.terpene", "NRPS.T1PKS.lanthipeptide-class-ii", "NRPS.betalactone", "CDPS", 
  "sactipeptide", "other", "epipeptide", "transAT-PKS", "PKS-like", "ranthipeptide", 
  "phosphonate", "lanthipeptide-class-iii"
)

# Pivot to long format

df_bgc_long <- df %>%
  select(`IsolateID`, all_of(bgc_cols)) %>%
  pivot_longer(cols = all_of(bgc_cols),
               names_to = "BGC_Type",
               values_to = "BGC_Count") %>%
  filter(BGC_Count > 0)  # optional: keep only present clusters

df_bgc_expanded <- df_bgc_long %>%
  uncount(weights = BGC_Count, .remove = FALSE) %>%
  mutate(BGC_Count = 1)

write_csv(df_bgc_expanded, "/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/antismash_bgc_groups_per_isolate.csv")

