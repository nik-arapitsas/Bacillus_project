######################################################################################################
# script name: antismash_bgcs_pie_charts.R
# developed by: Nikolaos P. Arapitsas
# framework: SarrisLab
######################################################################################################
# GOAL:
# Aim of this script is to create pie charts plotting the types of BGCs per isolate 
######################################################################################################
# usage:./antismash_bgcs_pie_charts.R
# complete path: /home/nik_arapitsas/Documents/Bacillus_project/scripts/antismash_bgcs_pie_charts.R
######################################################################################################

# Load the necessary libraries 
library(readr)
library(tidyverse)
install.packages("patchwork")
library(patchwork)
library(cowplot)

# Read the CSV file that contains the Isolate ID, BGC type, BGC count and Similarity Confidence without the species name 
bgcs_perisolate <- read_csv("bgc_similarity_with_species.csv")

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

# Express every BGC group number as a percent of the total BGCs for every isolate
bgcs_percent <- bgcs_perisolate %>%
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
  )) %>%
  group_by(IsolateID, Category) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(IsolateID) %>%
  mutate(Percent = Count / sum(Count) * 100)

# Create pie charts
make_pie <- function(data, isolate_id, colors, show_legend = FALSE) {
  ggplot(data, aes(x = "", y = Percent, fill = Category)) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = colors) +
    theme_void() +
    ggtitle(isolate_id) +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      legend.position = if (show_legend) "right" else "none"
    )
}

# List of isolates
isolates <- unique(bgcs_percent$IsolateID)

# Create one pie chart per isolate
plots <- lapply(seq_along(isolates), function(i) {
  iso <- isolates[i]
  make_pie(
    bgcs_percent %>% filter(IsolateID == iso),
    isolate_id = iso,
    colors = bgc_group_colors,
    show_legend = FALSE
  )
})

# Prepare the legend for the collective graph
dummy_legend_data <- data.frame(
  Category = names(bgc_group_colors),
  Percent = rep(1, length(bgc_group_colors)),  # equal slices
  IsolateID = "legend"
)

legend_plot <- make_pie(
  data = dummy_legend_data,
  isolate_id = "",
  colors = bgc_group_colors,
  show_legend = TRUE
) + theme(
  plot.title = element_blank(),
  legend.position = "bottom",
  legend.background = element_rect(fill = "white", color = NA),       
  legend.box.background = element_rect(fill = "white", color = NA)    
)

legend <- get_legend(legend_plot)

# Combine all graphs into one figure
pie_grid <- wrap_plots(plots, ncol = 5) +
  plot_layout(guides = "collect") & 
  theme(plot.margin = margin(0, 0, 0, 0))  # the 25 pies

bgcs_perisolate_pies <- plot_grid(
  pie_grid,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.2)
)

ggsave(paste0("bgcs_perisolate_pies",".png"),
       plot= bgcs_perisolate_pies, 
       height = 30, 
       width = 25,
       dpi = 300,
       bg = "white", 
       units="cm",
       device="png")
