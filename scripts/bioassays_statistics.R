#!/usr/bin/env Rscript
## Script name: bioassys_statistics.R
##
## Purpose of script: 
## Data summary and statistics of experiments.
##
## Author: Savvas Paragkamian
##
## Date Created: 2025-06-20
##
library(tidyverse)
library(readxl)
#library(vegan)

# load data

all_experiments_l <- read_excel("../data/in_vitro_pathogen_inhibition.xlsx",sheet="In vitro data") |>
    mutate(variable = paste0("var_",variable,sep="")) |>
    mutate(variable = gsub(" ","_",variable))

# check uniqueness of the line

all_experiments_distinct_v <- all_experiments_l |>
    group_by(microbe_id,pathogen,batch_id,plate_number,type_bioassay, variable, date,DPI) |>
    mutate(n=n()) |>
    ungroup()

# hyphea are counted 30 times in each plate
table(all_experiments_distinct_v$n,all_experiments_distinct_v$variable)

# bioassays summary
microbes_bioassay_summary <- all_experiments |>
    distinct(microbe_id,type_bioassay) |> 
    group_by(microbe_id) |> 
    summarise(n_type_bioassay=n(),
              conditions_tested=str_c(type_bioassay, collapse = ","))

# the control is named PDA, each batch id has it's control.
unique(microbes_bioassay_summary$microbe_id)

# summary of batches without the control-PDA
microbes_batches_summary <- all_experiments_l |>
    distinct(microbe_id,type_bioassay,pathogen, batch_id) |> 
    filter(microbe_id!="PDA") |>
    group_by(batch_id,pathogen,type_bioassay) |> 
    summarise(n_microbes=n(),
              microbe_id=str_c(microbe_id, collapse = ","),
              .groups="keep") 

# is a microbe tested in more than one batch? 
microbes_n_batches <- all_experiments_l |>
    distinct(microbe_id,type_bioassay,pathogen, batch_id) |> 
    filter(microbe_id!="PDA") |>
    group_by(microbe_id,type_bioassay,pathogen) |> 
    summarise(n_batch=n(), .groups="keep")

print("how many times appears the microbe_id,type_bioassay,pathogen?")
unique(microbes_n_batches$n_batch)

pathogen_summary <- microbes_batches_summary |>
    group_by(pathogen,type_bioassay) |>
    summarise(n_microbes=sum(n_microbes),
              microbe_id=str_c(microbe_id, collapse = ","),
              .groups="keep") |>
    mutate(n_unique = length(unique(str_split(microbe_id, ",")[[1]])))

write_delim(pathogen_summary, "../results/pathogen_summary.tsv",delim="\t")

### plates summary

plates_summary <- all_experiments_l |>
    group_by(microbe_id,pathogen,batch_id,type_bioassay, variable, DPI, date) |>
    summarise(plates=n(), .groups="keep")

plates_summary_unique <- all_experiments_l |>
    distinct(microbe_id,pathogen,batch_id,type_bioassay, plate_number) 

print("number of plates:")
nrow(plates_summary_unique)

print("number of plates without controls:")
plates_summary_unique |> filter(microbe_id!="PDA") |> nrow()

### dpi summary
### dpi is needed only for the var_Radial_growth_mm in order to calculate 
### the growth rate

dpi_summary <- all_experiments_l |>
    distinct(pathogen,batch_id,type_bioassay, DPI)

## check that date and DPI are the same in uniqueness

###################### radial growth ######################

radial_growth_rate_data <- all_experiments_l |> 
    filter(variable=="var_Radial_growth_mm") |>
    arrange(microbe_id, pathogen, batch_id, plate_number, type_bioassay, variable, date, DPI) |>
    group_by(microbe_id, pathogen, batch_id, plate_number, type_bioassay, variable) |>
    mutate(
           prev_value = lag(value, default = 5),
           prev_DPI   = lag(DPI, default = 0),
           growth_step = (value - prev_value) / (DPI - prev_DPI)
           ) |>
    summarise(
              growth_rate = sum(growth_step, na.rm = TRUE)/n(),
              .groups = "drop") |>
    rename("value"="growth_rate") |>
    mutate(variable="var_Radial_growth_rate")

all_experiments_final <- all_experiments_l |>
    filter(variable!="var_Radial_growth_mm") |>
    bind_rows(radial_growth_rate_data)


############################ descriptives ###########################
#### first seperate the different variables in order to 
### performe calculations for growth rate, mean, sd.

experiments_list <- split(all_experiments_final, all_experiments_final$variable)
variables <- names(experiments_list)

experiments_stats_l <- list()

for (var in seq_along(variables)){

    print(var)
    print(variables[var])
    variable_data <- experiments_list[[variables[var]]] 

    variable_summary <- variable_data |>
        group_by(microbe_id,pathogen,batch_id,type_bioassay,variable) |>
        summarise(mean=mean(value),
                  sderr=sd(value,na.rm=T) / sqrt(n())
                  ) |>
        ungroup() |>
        mutate(batch_bioassay=paste0(type_bioassay," batch ",batch_id,sep=""))

    experiments_stats_l[[var]] <- variable_summary
    
    fig_barplots <- ggplot(variable_summary,
                      aes(x = microbe_id,
                          y = mean)) + 
                geom_col(width = 0.6,
                         fill="lightseagreen")+
                geom_errorbar(aes(ymin = mean - sderr,
                                  ymax = mean + sderr),
                              width = 0.1) +  # Error bars
                labs(
                     title = "Mean Values with Standard Error",
                     x = "Microbe ID",
                     y = paste0(variables[var])) +
                theme_bw() +
                facet_grid(batch_bioassay ~ pathogen) +
                theme(
                      legend.position = "top",
                      legend.direction = "horizontal",
                      legend.box = "horizontal",
                      axis.text.y = element_text(size = 15),
                      axis.text.x = element_text(angle = 90,
                                                 hjust = 1,
                                                 size = 15),
                      axis.title.y = element_text(size = 16),
                      strip.text.x = element_text(face = "italic")
                ) 
    
    ggsave(paste0("../figures/phytopathogen_invitro_",variables[var],"_bar.png"),
           plot=fig_barplots, 
           height = 30, 
           width = 50,
           dpi = 300, 
           units="cm",
           device="png")
    
}

experiments_stats <- bind_rows(experiments_stats_l)

############################# Normalisation ##############################
## how to compare with the control values?
# Calculate percentage change
# [(Vc-Vt)/Vc] ×100 
all_experiments_norm <- experiments_stats |>
    group_by(pathogen, batch_id, type_bioassay, variable) |>
    mutate(control_value = mean[microbe_id == "PDA"], #Extract control value
         percent_change = ((mean - control_value) / control_value) * 100,
         difference_value = mean - control_value) |>
  ungroup() |>
  filter(microbe_id!="PDA") |>
  select(-control_value) # Optional: Remove intermediate control_value column

all_experiments_percent_change <- all_experiments_norm |>
    dplyr::select(batch_id,microbe_id,pathogen,type_bioassay,variable,percent_change,mean,sderr) |>  # Remove the original value column
    mutate(percent_change=round(percent_change,2)) 


################################# Statistics ##############################
all_experiments_nest <- all_experiments_final |>
    group_by(batch_id, type_bioassay, pathogen, variable) |>
    nest() 

## anova and tukey
statistics <- all_experiments_nest |>
    mutate(
           model = map(data, ~ aov(value ~ microbe_id, data = .x)),
           anova = map(model,~ broom::tidy(.x)),
           tukey = map(model, ~ TukeyHSD(.x)$microbe_id %>% as.data.frame() %>%
                       rownames_to_column("comparison") %>% rename("p.value"="p adj"))
    )

## extract anova values
anova_pvalues <- statistics |>
    select(batch_id, type_bioassay, pathogen, variable, anova) |>
    unnest(anova) |>
    select(batch_id, type_bioassay, pathogen, variable, term, p.value) |>
    rename(anova_term = term, anova_p = p.value) |>
    filter(anova_term=="microbe_id") |>
    ungroup()

## extract tukey values and filter the comparisons with PDA (control)
tukey_pvalues <- statistics |>
    select(batch_id, type_bioassay, pathogen, variable, tukey) |>
    unnest(tukey) |>
    select(batch_id, type_bioassay, pathogen, variable, comparison, p.value) |>
    rename(tukey_comparison = comparison, tukey_p = p.value) |>
    separate(tukey_comparison, into=c("microbe_id", "comparison"),sep="-") |>
    filter(comparison=="PDA") |>
    ungroup()


all_experiments_statistics <- all_experiments_percent_change |>
    left_join(anova_pvalues,by=c("batch_id", "type_bioassay", "pathogen", "variable")) |>
    left_join(tukey_pvalues,by=c("batch_id", "type_bioassay", "pathogen", "variable","microbe_id"))

write_delim(all_experiments_statistics,"../results/phytopathogen_invitro_statistics.tsv",delim="\t")


### significant pairwise microbes with values higher than control

control_pairwise_sig_sum <- all_experiments_statistics |>
    filter(percent_change < 0) |>
    filter(anova_p < 0.05) |>
    filter(tukey_p < 0.05) |>
    mutate(
           significance = case_when(
                                    tukey_p < 0.001 ~ "***",
                                    tukey_p < 0.01 ~ "**",
                                    tukey_p < 0.05 ~ "*",
                                    TRUE ~ ""  # For non-significant p-values
           )
    ) |>
    mutate(variable=gsub("var_","",variable)) |>
    mutate(variable=gsub("_"," ",variable)) |>
    mutate(facet_vars=paste0(pathogen," - ",type_bioassay,sep=""))


microbe_heatmap <- ggplot() +
    geom_tile(data = control_pairwise_sig_sum,
              aes(x = microbe_id, y = variable, fill = percent_change),
              color = "white",       # <- tile borders
              linewidth = 0.5,       # <- optional: control border thickness
              alpha = 1,
              show.legend = TRUE) +
    geom_text(data = control_pairwise_sig_sum,
              aes(x = microbe_id, y = variable, label = significance,color=significance),
              show.legend = TRUE,
              size = 4) +
    scale_fill_gradient(
      low = "#0072B2",
      high = "gray87",
      breaks = waiver(),
      n.breaks = 4,
      limits = c(
        min(control_pairwise_sig_sum$percent_change),
        max(control_pairwise_sig_sum$percent_change)
      ),
      na.value = "white",
      guide = "colorbar"
    ) +
    scale_color_manual(
                       name = "Significance\n(p-value)",
                       values = c("***" = "black", "**" = "black", "*" = "black"),
                       breaks = c("***", "**", "*"),
                       labels = c("(<0.001)", "(<0.01)", "(<0.05)"))+
    guides(fill = guide_colorbar(
                                 ticks = FALSE,
                                 title = "Percent change",
                                 label.vjust = 0.8,
                                 title.vjust = 0.8
                                 ),
           color = guide_legend(
                                override.aes = list(fill = "white",
                                                    label = c("***", "**", "*"),
                                                    size = 5))
           ) +

   # coord_fixed() +  # square tiles
    theme_bw() +
    facet_wrap(~ facet_vars, nrow = 4,scales = "free_y") +
    theme(
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
          panel.grid = element_blank(),  # no gridlines
          axis.text.x = element_text(face = "bold", angle = 90, hjust = 0),
      axis.text.y = element_text(face = "bold"),
      axis.text = element_text(size = 13),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.title = element_text(size = 9),
      strip.background = element_blank(),
      strip.text.x = element_text(face = "bold.italic", size = 12)
    ) 


ggsave("../figures/significant_microbes_control.png",
       plot = microbe_heatmap,
       width = 40,
       height = 20,
       units='cm', 
       device = "png",
       dpi = 300)

