########################################################
#### This script generates Figure 7 from the paper. ####
########################################################

## Load the relevant packages
library(stringr)
library(dplyr)
library(gtools)
library(latex2exp)
library(ggpubr)
library(ggplot2)
library(readr)
library(data.table)
library(patchwork)

## Set seed for reproducibility ----
set.seed(90210)

## Set a colorblind friendly palate to use for visualizations
cbbPalette <- c("#56B4E9", "#009E73", "#D55E00")

## Read in and process the specified dataset
filename <- paste0("./1_simulation_study/simulated_data/med_dens_medium_noise_1_alpha_sim_8.csv")
scan_data <- read_csv(filename, show_col_types = FALSE)

## Plot a realization of a simulated medium density dataset
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
file2015 <- read_csv("./resources/empirical_data/UER_lidar_canopy_segmentation/crown_attributes_2015.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(XTOP, YTOP, CANVOL2015) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2015) %>%
  mutate(file = 1)
file2019 <- read_csv("./resources/empirical_data/UER_lidar_canopy_segmentation/crown_attributes_2019.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(XTOP, YTOP, CANVOL2019) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2019) %>%
  mutate(file = 2)

## Plot the medium density subset of the empirical dataset
scan_data_emp <- rbind(file2015, file2019)
scan_data_emp |> ggplot(aes(x = x, y = y,
                            color = factor(file,
                            levels = c("1", "2"),
                            labels = c("Observed Data 1", "Observed Data 2")))) +
  geom_point(aes(size = size), alpha = .5) +
  scale_color_manual(values = cbbPalette[1:2]) +
  labs(x = "", y = "", size = TeX('Size in $m^3$'), color = "Dataset") +
  theme_bw() + theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14), legend.position = "bottom") -> emp_med_plot


## Create and save the two panel data comparison plot (Figure 7)
data_comp_plot <- sim_med_plot + emp_med_plot + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") &
  theme(legend.position='bottom')

ggsave(filename = "7_medium_density_data_comparison_plot.png", plot = data_comp_plot, path = "./3_figures_and_tables/",
       width = 20, height = 12, units = "cm", dpi = "retina")
