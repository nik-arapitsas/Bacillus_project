library(ggplot2)
library(dplyr)

# Read your data
region_count_per_isolate <- read.table("/home/nik_arapitsas/Documents/Bacillus_project/Results/Antismash/Antismash_Output/region_counts.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Create the graph using ggplot
ggplot(region_count_per_isolate, aes(x = reorder(isolate, -region_count), y = region_count, fill = isolate)) +
  geom_bar(stat = "identity", width = 0.6, fill = "black") +
  coord_flip() +
  labs(
    title = "antiSMASH regions per isolate",
    x = "Isolate",
    y = "Region count"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),  # Remove grid lines 
        legend.position = "none",
        axis.line = element_line(color = "black", size = 0.2),  # Add black axis lines
        axis.ticks.length = unit(0.1, 'cm'),  # Set the tick length to be smaller
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center align title
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, 
                                   size = 10, color = "black", family = "sans", 
                                   margin = margin(t = 2.5)),  # Darker and more spaced letters
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black", family = "sans", 
                                   margin = margin(t = 2.0)),
        plot.margin = margin(10, 10, 1, 0.5))  # Add extra space around the plot)

