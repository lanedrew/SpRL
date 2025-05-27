########################################################
#### This script generates Figure 8 from the paper. ####
########################################################

## Load the relevant packages
library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(gtools)
library(latex2exp)
library(ggpubr)
library(purrr)
library(stringr)
library(kableExtra)
library(tidyr)
library(patchwork)
library(data.table)


## Set a colorblind friendly palate to use for visualizations
cbbPalette <- c("#56B4E9", "#009E73", "#D55E00")


## Define a function to add the filename to each line of the imported dataset
read_plus <- function(flnm) {
  fread(flnm) %>%
    .[, filename := flnm] # Add the filename column directly in the data.table
}


## Read in all of the results .csv's into one tibble with the filename attached to each line
sprl_sim_table <-
  list.files(path = "./1_simulation_study/simulation_results/NDM/linkage_processed",
             pattern = "1_alpha", 
             full.names = T) |>  
  mixedsort() |>
  lapply(read_plus) |>
  rbindlist(use.names = TRUE, fill = TRUE)


## Modify the names from the files
myNames <- gsub("./1_simulation_study/simulation_results/NDM/linkage_processed/", "", sprl_sim_table$filename, fixed = TRUE)
myNames <- gsub("_NDM_linkage_metrics.csv", "", myNames, fixed = TRUE)
variables <- do.call(rbind, str_split(myNames, "_"))

sprl_sim_table <- sprl_sim_table %>% 
  mutate(filename = myNames,
         density = variables[,1],
         noise = variables[,3],
         alpha = as.numeric(variables[,5]),
         index = as.numeric(variables[,7]))

ndm_sim_table <- sprl_sim_table


## Read in all of the results .csv's into one tibble with the filename attached to each line
sprl_sim_table <-
  list.files(path = "./1_simulation_study/simulation_results/two_stage/linkage_processed/",
             pattern = "1_alpha", 
             full.names = T) |>  
  mixedsort() |>
  lapply(read_plus) |>
  rbindlist(use.names = TRUE, fill = TRUE)


## Modify the names from the files
myNames <- gsub("./1_simulation_study/simulation_results/two_stage/linkage_processed/", "", sprl_sim_table$filename, fixed = TRUE)
myNames <- gsub("_N_thresh_two_stage_linkage_metrics.csv", "", myNames, fixed = TRUE)

variables <- do.call(rbind, str_split(myNames, "_"))

sprl_sim_table <- sprl_sim_table %>% 
  mutate(filename = myNames,
         density = variables[,1],
         alpha = as.numeric(variables[,5]),
         noise = variables[,3],
         index = as.numeric(variables[,7]))

rl_sim_table <- sprl_sim_table

linkage_summary_table_la <- rl_sim_table %>%
  group_by(noise, density, alpha, index) %>%
  summarise(precision_mean = mean(precision), precision_q5 = quantile(precision, .05), precision_q95 = quantile(precision, .95),
            recall_mean = mean(recall), recall_q5 = quantile(recall, .05), recall_q95 = quantile(recall, .95))


linkage_join <- full_join(linkage_summary_table_la, ndm_sim_table, by = c("noise", "density", "index", "alpha"))


linkage_colors <- c("Linkage-Averaging" = "#56B4E9", "Nearest Distance Matching" = "#009E73")

prec_plot <- linkage_join |>  mutate(noise = factor(noise,
                                                    levels = c("small", "medium", "large"),
                                                    labels = c("Small Noise", "Medium Noise", "Large Noise")),
                                     density = factor(density,
                                                      levels = c("low", "med", "high"),
                                                      labels = c("Low Density", "Medium Density", "High Density"))) |>
  ggplot(aes(x = precision_mean, y = index)) +
  geom_point(alpha = .75, aes(color = "Linkage-Averaging")) +
  geom_point(aes(x = precision, y = index, color = "Nearest Distance Matching"), alpha = .75) +
  geom_linerange(aes(xmin = precision_q5, xmax = precision_q95, color = "Linkage-Averaging"), alpha = .5) +
  facet_grid(noise ~ density) +
  labs(x = "Precision", y = "Dataset Index", color = "Linkage Approach") +
  scale_color_manual(values = linkage_colors) +
  theme_bw() + theme(text = element_text(family = "serif", size = 14), legend.position = "bottom")


rec_plot <- linkage_join |>  mutate(noise = factor(noise,
                                                   levels = c("small", "medium", "large"),
                                                   labels = c("Small Noise", "Medium Noise", "Large Noise")),
                                    density = factor(density,
                                                     levels = c("low", "med", "high"),
                                                     labels = c("Low Density", "Medium Density", "High Density"))) |>
  ggplot(aes(x = recall_mean, y = index)) +
  geom_point(alpha = .75, aes(color = "Linkage-Averaging")) +
  geom_point(aes(x = recall, y = index, color = "Nearest Distance Matching"), alpha = .75) +
  geom_linerange(aes(xmin = recall_q5, xmax = recall_q95, color = "Linkage-Averaging"), alpha = .5) +
  facet_grid(noise ~ density) +
  labs(x = "Recall", y = "Dataset Index", color = "Linkage Approach") +
  scale_color_manual(values = linkage_colors) +
  theme_bw() + theme(text = element_text(family = "serif", size = 14), legend.position = "bottom")

## Create and save the two panel plot (Figure 8)
rl_metrics <- prec_plot + rec_plot + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") &
  theme(legend.position='bottom')

ggsave(filename = "8_precision_recall_ndm_la_comparison_plot.png", plot = rl_metrics, path = "./3_figures_and_tables/",
       width = 30, height = 25, units = "cm", dpi = "retina")
