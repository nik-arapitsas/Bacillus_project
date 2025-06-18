# Load necessary library
library(ggplot2)
install.packages("ggpubr", dependencies = TRUE)
install.packages("RcppEigen")
install.packages("ggpubr")
library(ggpubr)

# 1 Number of Isolate-Specific Orthogroups per Isolate 

# Read the sorted data into R
data <- read.table("/media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/species_specific_orthogroups.txt", header=FALSE, col.names=c("Isolates", "IsolatesSpecificOrthogroups"))

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

IsolateSpecificOgs_plot <- ggplot(data, aes(x = Isolates, y = IsolatesSpecificOrthogroups)) +
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
    axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black", family = "sans", 
                               margin = margin(t = 2.0)),
    plot.margin = margin(10, 10, 1, 0.5),  # Add extra space around the plot
    legend.position = "none"  # Remove the legend
  ) +
  scale_y_continuous(breaks = c(10, 25, 50, 110)) +
  coord_cartesian(expand = FALSE, ylim = c(0, max(data$IsolatesSpecificOrthogroups) * 1.04)) # Extend the x-axis

ggsave(paste0("/media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/no_of_isolate_specific_ogs",".png"),
       plot=IsolateSpecificOgs_plot, 
       height = 20, 
       width = 50,
       dpi = 300, 
       units="cm",
       device="png")

#Drafts for 1

library(grid)
# Manually adjust spacing of x-axis labels
annotation_custom(grob = textGrob("Isolates, test. samples", gp = gpar(fontsize = 10, col = "black", cex = 1.2)))

  scale_y_continuous(breaks = seq(0, 110, 10))
ylim = c(0, max(data$IsolatesSpecificOrthogroups) * 1.1)

# 2 Prepare graph for percentage of genes from each isolate assigned to orthogroups

# Read the sorted data into R
data_gene_perc <- read.table("/media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/percofgenes_inogs_per_isolate.txt", header=FALSE, col.names=c("Isolates", "PercentageinOgs"))

# Convert Species to factor to maintain sorting order
data_gene_perc$Isolates <- factor(data_gene_perc$Isolates, levels = data_gene_perc$Isolates)

# Create the bar plot

#With legend in y axis

PercOfOgs_plot <- ggplot(data_gene_perc, aes(x = Isolates, y = PercentageinOgs)) +
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
    axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black", family = "sans", 
                               margin = margin(t = 2.0)),
    plot.margin = margin(10, 10, 10, 10),  # Add extra space around the plot
    legend.position = "none"  # Remove the legend
  ) +
  coord_cartesian(expand = FALSE)  # Extend the x-axis

ggsave(paste0("/media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/perc_of_genes_in_ogs_plot",".png"),
       plot=PercOfOgs_plot, 
       height = 20, 
       width = 50,
       dpi = 300, 
       units="cm",
       device="png")


# 3 Prepare graph for orthogroups in every isolate with orthogroups in all or any isolates

library(reshape2)

# Read the sorted data into R
data_gene_all_any <- read.table("/media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/orthogroupcount_in_isolates.txt", header=TRUE, sep = "\t")

# Convert Species to factor to maintain sorting order
data_gene_all_any$Isolates <- factor(data_gene_perc$Isolates, levels = data$Isolates)

data_long <- melt(data_gene_all_any, id.vars = "Isolates", variable.name = "Category", value.name = "Count")

data_any <- subset(data_long, Category == "Partially.Shared.Orthogroups")   # Blue bars (background)
data_all <- subset(data_long, Category == "Core.Orthogroups")  # Green bars (foreground)

# Plot without legend

ggplot() +
  geom_bar(data=data_any, aes(x=Isolates, y=Count), stat="identity", fill="blue") +   # Blue first
  geom_bar(data=data_all, aes(x=Isolates, y=Count), stat="identity", fill="green") +  # Green on top
  labs(y="Number of Genes", x="", fill="Category") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.2),
        axis.ticks.length = unit(0.1, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center align title
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10, color = "black", family = "sans", margin = margin(t = 2.5)),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black", family = "sans", 
          margin = margin(t = 2.0)),
        plot.margin = margin(10, 10, 10, 10)) +
  scale_y_continuous(expand = c(0, 0))   # Rotate x-axis labels 

# Plot with legend

data_combined <- rbind(data_any, data_all)
data_combined$Category <- factor(data_combined$Category, 
                                 levels = c("Partially.Shared.Orthogroups", "Core.Orthogroups"))

SharedOgs_plot <- ggplot(data_combined, aes(x = Isolates, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "identity") +  # Maintain layering
  scale_fill_manual(values = c("Partially.Shared.Orthogroups" = "blue", 
                               "Core.Orthogroups" = "green"),
                    labels = c("Partially.Shared.Orthogroups" = "in any isolate", 
                               "Core.Orthogroups" = "in all isolates")) +
  labs(title = "Number of Genes with Orthogroups", y = "", x = "", fill = "Category") +
  guides(fill = guide_legend(title = NULL, byrow = TRUE)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black", size = 0.2),
        axis.ticks.length = unit(0.1, "cm"),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        legend.position = c(0.89, 0.89),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.text = element_text(size = 10),  # Reduce font size of legend
        legend.spacing.y = unit(0.2, 'cm'),
        legend.key.size = unit(0.5, "cm"),  # Adjust size of legend keys (squares)
        legend.margin = margin(-2, 15, 4, 1.2),  # Reduce legend margins
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Center align title
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 10, color = "black", family = "sans", 
                                   margin = margin(t = 2.5)),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5, size = 10, color = "black", family = "sans", 
                                   margin = margin(t = 2.0)),
        plot.margin = margin(10, 1, 1, 0.1)) + 
  coord_cartesian(expand = FALSE, ylim = c(0,6000))

ggsave(paste0("/media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/no_of_genes_with_ogs_inanyorall",".png"),
       plot=SharedOgs_plot, 
       height = 20, 
       width = 50,
       dpi = 300, 
       units="cm",
       device="png")

SharedOgs_plot <- SharedOgs_plot + theme(legend.position = "bottom", legend.margin = margin(t = -18, unit = "pt"))
  
collective_orthofinder_graph <- ggarrange(IsolateSpecificOgs_plot, PercOfOgs_plot, SharedOgs_plot,
          labels = c("A", "B", "C"),                                
          ncol = 1, 
          nrow = 3, 
          align = "v")  # Align vertically by x-axis

ggsave(paste0("/media/sarlab/DATA/Bacillus_project/Bacillus_project_orthofinder/Results_Jun11/Graphs/collective_orthofinder_graph",".png"),
       plot=collective_orthofinder_graph, 
       height = 30, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")



# Draft 3

data_long$Category <- factor(data_long$Category, levels = c("Core.Orthogroups", "Partially.Shared.Orthogroups"))
# Create the bar plot

ggplot(data_long, aes(x=Isolates, y=Count, fill=Category)) +
  geom_bar(stat="identity", position="identity", alpha=0.8) +  # Stacked bars with transparency
  scale_fill_manual(values=c("Partially.Shared.Orthogroups"="blue", "Core.Orthogroups"="green")) +  # Set custom colors
  labs(y="Number of Genes", x="", fill="Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1))  # Rotate x-axis labels for better readability
