# Load necessary libraries
library(tidyverse)

# Read the CSV file (adjust the path if needed)
df <- read_csv("/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/antismash_table.csv")

# Rename columns to make them easier to work with
df <- df %>%
  rename(IsolateID = `Isolate ID`,
         Low = `Low Similarity`,
         Medium = `Medium Similarity`,
         High = `High Similarity`,
         Undefined = `Undefined Similarity`)

# Pivot to long format for ggplot
df_long <- df %>%
  select(IsolateID, High, Medium, Low, Undefined) %>%
  pivot_longer(cols = High:Undefined, names_to = "Similarity", values_to = "Count")

# Set stacking order
df_long$Similarity <- factor(df_long$Similarity, levels = c("Undefined", "Low", "Medium", "High"))

# Define custom colors
similarity_colors <- c(
  "High" = "#f8d48c",
  "Medium" = "#d4c63a",
  "Low" = "#ff94b4",
  "Undefined" = "#7171be"
)

# Plot
antiSMASH_regions_similarity_count_plot <- ggplot(df_long, aes(x = IsolateID, y = Count, fill = Similarity)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = similarity_colors, name = "Similarity Confidence") +
  scale_x_discrete(expand=expansion(add=c(.8, .8))) +
  scale_y_continuous(expand=expansion(add=c(0, 0))) +
  labs(title = "antiSMASH region count per isolate",
       x = NULL,
       y = "Number of antiSMASH Regions") +
  theme_minimal() +
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
          coord_cartesian(ylim = c(0,20))

ggsave(paste0("/media/sarlab/DATA/Bacillus_project/Antismash_Graphs/antiSMASH_regions_similarity_count_plot",".png"),
       plot=antiSMASH_regions_similarity_count_plot, 
       height = 20, 
       width = 50,
       dpi = 300, 
       units="cm",
       device="png")
