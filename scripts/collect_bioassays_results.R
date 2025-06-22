######################################################################################################
# script name: collect_bioassays_results.R
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
######################################################################################################
# GOAL:
# Aim of this script is to provide a table that includes th bioassay result for everyone of our 25 isolates
######################################################################################################
# usage:./collect_bioassays_results.R
# complete path: /home/nik_arapitsas/Documents/Bacillus_project/scripts/collect_bioassays_results.R
######################################################################################################

# Install and load the appropriate packages
install.packages("readxl")
install.packages("writexl")
library("readxl")
library("writexl")

# Open the .xlsx file
collective_table_with_bioassays_results <- read_excel("/Users/user/Desktop/All_strains_Bioassay_results_inhibition intensity added.xlsx")

# Define the numbers of the isolates of interest
target_ids <- c(152, 163, 179, 215, 218, 221, 224, 244, 266,
                335, 337, 340, 342, 368, 369, 374, 379, 389,
                398, 543, 544, 571, 656, 658, 662)

# Filter the rows where the first column matches one of the isolate IDs
isolates_collective_table_bioassays <- collective_table_with_bioassays_results[collective_table_with_bioassays_results$`Strain No` %in% target_ids, ]

# Add "SRL" prefix to "Strain No"
isolates_collective_table_bioassays$`Strain No` <- paste0("SRL", isolates_collective_table_bioassays$`Strain No`)

# Save the final table to a new Excel file
write_xlsx(isolates_collective_table_bioassays, "/Users/user/Desktop/isolates_collective_table_bioassays.xlsx")
