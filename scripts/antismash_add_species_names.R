######################################################################################################
# script name: antismash_add_species_names.R
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
######################################################################################################
# GOAL:
# Aim of this script is to add the species names for every isolate in the 
# antismash_bgc_groups_per_isolate_similarity_added table 
######################################################################################################
# usage:./antismash_add_species_names.R
# complete path: /home/nik_arapitsas/Documents/Bacillus_project/scripts/antismash_add_species_names.R
######################################################################################################


# Read the CSV file that contains the Isolate ID, BGC type, BGC count and Similarity Confidence without the species name 
bgc_typecount_similarity_perisolate <- read_csv("/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/antismash_bgc_groups_per_isolate_similarity_added.csv")

# Read the CSV file that contains the Isolate ID and the species names 
table_with_species_names <- read_csv2("/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/antismash_table_bgc_groups.csv")

table_with_species_names <- table_with_species_names %>%
  rename(IsolateID = `Isolate ID`,
         Species = `Species name`)

# Add the species names to every isolate on the bgc_typecount_similarity_perisolate table
bgc_typecount_similarity_perisolate <- merge(bgc_typecount_similarity_perisolate, table_with_species_names[, c("IsolateID", "Species")], by = "IsolateID", all.x = TRUE)

#Reorder columns so Species comes right after IsolateID
column_names <- colnames(bgc_typecount_similarity_perisolate)
new_order <- c("IsolateID", "Species", setdiff(column_names, c("IsolateID", "Species")))
bgc_typecount_similarity_perisolate <- bgc_typecount_similarity_perisolate[, new_order]

# Remove underscores and trailing single-letter suffixes from the end
bgc_typecount_similarity_perisolate$Species <- gsub("_[A-Z]$", "", bgc_typecount_similarity_perisolate$Species)  # Remove "_A", "_S", etc. at end

# Remove internal GTDB-style suffixes like "_A", "_AB", etc.
bgc_typecount_similarity_perisolate$Species <- gsub("_[A-Z]+ ", " ", bgc_typecount_similarity_perisolate$Species)

# Replace any remaining underscores with spaces (just in case)
bgc_typecount_similarity_perisolate$Species <- gsub("_", " ", bgc_typecount_similarity_perisolate$Species)
