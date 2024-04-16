###########################################################################################
#### This script processes all of the results from the joint and two-stage SpRL models ####
###########################################################################################

## Load the relevant packages ##
# library(tidyverse)
library(stringr)
library(dplyr)
library(gridExtra)
library(gtools)
library(latex2exp)
library(ggpubr)
library(ggplot2)
library(readr) ## load and save results
library(data.table)
library(patchwork)

## Set seed for reproducibility ----
set.seed(90210)
cbbPalette <- c("#56B4E9", "#009E73", "#D55E00")

## Read in and process the specified dataset
filename <- paste0("./code/growth_sim_data_F23/med_dens_medium_noise_1_alpha_sim_8.csv")
scan_data <- read_csv(filename, show_col_types = FALSE)

scan_data |> filter(!file == 0) |> ggplot(aes(x = x, y = y,
                                              color = factor(file,
                                                             levels = c("1", "2"),
                                                             labels = c("Observed Data 1", "Observed Data 2")))) +
  geom_point(aes(size = size), alpha = .5, show.legend = FALSE) +
  scale_color_manual(values = cbbPalette[1:2]) +
  labs(x = "", y = "", size = TeX('Size in $m^3$'), color = "Dataset") +
  theme_bw() + theme(axis.text = element_blank(), legend.position = "bottom", text = element_text(family = "serif", size = 14)) -> sim_med_plot


## For medium density area of size 100m^2
a_x <- 326996
a_y <- 4311239
b_x <- 327096
b_y <- 4311339

## Read in and process the specified dataset
file2015 <- read_csv("./data/UER_lidar_canopy_segmentation/crown_attributes_2015.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(XTOP, YTOP, CANVOL2015) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2015) %>%
  mutate(file = 1)
file2019 <- read_csv("./data/UER_lidar_canopy_segmentation/crown_attributes_2019.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(XTOP, YTOP, CANVOL2019) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2019) %>%
  mutate(file = 2)


scan_data_emp <- rbind(file2015, file2019)
scan_data_emp |> ggplot(aes(x = x, y = y,
                            color = factor(file,
                            levels = c("1", "2"),
                            labels = c("Observed Data 1", "Observed Data 2")))) +
  geom_point(aes(size = size), alpha = .5) +
  scale_color_manual(values = cbbPalette[1:2]) +
  labs(x = "", y = "", size = TeX('Size in $m^3$'), color = "Dataset") +
  theme_bw() + theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14), legend.position = "bottom") -> emp_med_plot

data_comp_plot <- sim_med_plot + emp_med_plot + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") &
  theme(legend.position='bottom')

ggsave(filename = "med_dens_data_comp_plot.png", plot = data_comp_plot, path = "./plots/F23/",
       width = 20, height = 12, units = "cm", dpi = "retina")

