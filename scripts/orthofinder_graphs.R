# Load necessary library
library(ggplot2)

# 1 Number of Isolate-Specific Orthogroups per Isolate 

# Read the sorted data into R
data <- read.table("/home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Graphs/species_specific_orthogroups.txt", header=FALSE, col.names=c("Isolates", "IsolatesSpecificOrthogroups"))

# Convert Species to factor to maintain sorting order
  data$Isolates <- factor(data$Isolates, levels = data$Isolates[order(data$IsolatesSpecificOrthogroups)])

#Check for duplicates
duplicates <- data[duplicated(data$Species), ]
print(duplicates)

# Create the bar plot

#With legend in y axis

ggplot(data, aes(x = Isolates, y = IsolatesSpecificOrthogroups)) +
  geom_bar(stat = "identity", fill = "black", position = "identity") +  # Set bars to black
  theme_minimal() +
  labs(title = "Number of Isolate-Specific Orthogroups",
       x = "", y = "Isolate-Specific Orthogroups") +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
    axis.ticks.length = unit(0.1, 'cm'),  # Set the tick length to be smaller
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center align title
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, 
                               size = 10, color = "black", family = "sans", 
                               margin = margin(t = 2.5)),  # Darker and more spaced letters
    plot.margin = margin(10, 10, 10, 10),  # Add extra space around the plot
    legend.position = "none"  # Remove the legend
  ) +
  scale_y_continuous(breaks = c(5, 10, 25, 50, 110)) +
  coord_cartesian(expand = FALSE, ylim = c(0, max(data$IsolatesSpecificOrthogroups) * 1.04))  # Extend the x-axis


#Without legend in y axis

ggplot(data, aes(x = Isolates, y = IsolatesSpecificOrthogroups)) +
  geom_bar(stat = "identity", fill = "black", position = "identity") +  # Set bars to black
  theme_minimal() +
  labs(title = "Number of Isolate-Specific Orthogroups",
       x = "", y = "") +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
    axis.ticks.length = unit(0.1, 'cm'),  # Set the tick length to be smaller
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center align title
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, 
                               size = 10, color = "black", family = "sans", 
                               margin = margin(t = 2.5)),  # Darker and more spaced letters
    plot.margin = margin(10, 10, 1, 0.5),  # Add extra space around the plot
    legend.position = "none"  # Remove the legend
  ) +
  scale_y_continuous(breaks = c(5, 10, 25, 50, 110)) +
  coord_cartesian(expand = FALSE, ylim = c(0, max(data$IsolatesSpecificOrthogroups) * 1.04)) # Extend the x-axis

#Drafts for 1

library(grid)
# Manually adjust spacing of x-axis labels
annotation_custom(grob = textGrob("Isolates, test. samples", gp = gpar(fontsize = 10, col = "black", cex = 1.2)))

  scale_y_continuous(breaks = seq(0, 110, 10))
ylim = c(0, max(data$IsolatesSpecificOrthogroups) * 1.1)

# 2 Prepare graph for percentage of genes from each isolate assigned to orthogroups

# Read the sorted data into R
data_gene_perc <- read.table("/home/nik_arapitsas/Documents/Bacillus_project/Results/orthofinder/Results_Feb03/Graphs/percofgenes_inogs_per_isolate.txt", header=FALSE, col.names=c("Isolates", "PercentageinOgs"))

# Convert Species to factor to maintain sorting order
data_gene_perc$Isolates <- factor(data_gene_perc$Isolates, levels = data_gene_perc$Isolates)

# Create the bar plot

#With legend in y axis

ggplot(data_gene_perc, aes(x = Isolates, y = PercentageinOgs)) +
  geom_bar(stat = "identity", fill = "black", position = "identity") +  # Set bars to black
  geom_hline(yintercept = 100, linetype = "dashed", color = "black", size = 0.2) +  # Dashed line at 100%
  theme_minimal() +
  labs(title = "Percentage of Genes in Orthogroups",
       x = "", y = "Percentage of genes (%)") +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
    axis.ticks.length = unit(0.1, 'cm'),  # Set the tick length to be smaller
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center align title
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, 
                               size = 10, color = "black", family = "sans", 
                               margin = margin(t = 2.5)),  # Darker and more spaced letters
    plot.margin = margin(10, 10, 10, 10),  # Add extra space around the plot
    legend.position = "none"  # Remove the legend
  ) +
  coord_cartesian(expand = FALSE)  # Extend the x-axis
