###############################################################
#### This script generates Figures 3 and 4 from the paper. ####
###############################################################

## Load libraries and sampler functions ----
library(readr) 
library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(dplyr)
library(purrr)
library(stringr)
library(gtools)
library(ggplot2)
library(mcclust)
library(latex2exp)
library(reshape2)

## Load C++ helper functions
sourceCpp('./resources/code/cpp_code/two_stage_func.cpp')

## Set a colorblind friendly palate to use for visualizations
cbbPalette <- c("#56B4E9", "#009E73", "#D55E00", "#F0E442")

## Define a function to add the filename to each line of the imported dataset
read_plus <- function(flnm) {
  read_csv(flnm, show_col_types = FALSE, skip = 2501, col_names = FALSE) %>% 
    mutate(filename = flnm)
}

## Read in all of the timing results .csv's into one tibble with the filename attached to each line
emp_timing_table <-
  list.files(path = "./2_empirical_analysis/timing_results",
             pattern = "timing_results",
             full.names = T) |>
  mixedsort() |>
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./2_empirical_analysis/timing_results/", "", emp_timing_table$filename, fixed = TRUE)
myNames <- gsub(".csv", "", myNames, fixed = TRUE)
variables <- do.call(rbind, str_split(myNames, "_"))


## Subset the data and generate summary statistics datafrane
emp_timing_table <- emp_timing_table %>%
  mutate(filename = myNames,
         box_size = as.numeric(variables[,6]),
         area = variables[,8]) |> 
  rename(time = X1)
emp_timing_summary <- emp_timing_table |> group_by(box_size, area) |> summarise(mean_time = mean(time))


## Plot the full timing data
emp_timing_summary |>
  ggplot(aes(x = box_size, y = mean_time, group = area, color = as.factor(area), shape = as.factor(area))) +
  geom_point() +
  geom_line() +
  labs(x = "Bounding Box Margin",
       y = "Average Time Per Iteration (Seconds)",
       title = "",
       color = TeX("Domain Size $(m^2)$"),
       shape = TeX("Domain Size $(m^2)$")) +
  scale_colour_manual(values = cbbPalette) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom",
    text = element_text(family = "serif", size = 12)) -> emp_timing_plot


## Plot the subsetted timing data
emp_timing_summary |> filter(box_size <= 10) |> 
  ggplot(aes(x = as.numeric(box_size), y = mean_time, group = area, color = as.factor(area), shape = as.factor(area))) +
  geom_point() +
  geom_line() +
  labs(x = "Bounding Box Margin",
       y = "Average Time Per Iteration (Seconds)",
       title = "",
       color = TeX("Domain Size $(m^2)$"),
       shape = TeX("Domain Size $(m^2)$")) +
  scale_colour_manual(values = cbbPalette) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    text = element_text(family = "serif", size = 12)) -> emp_timing_10_inset


## Create and save the two panel timing results plot (Figure 4 in paper)
timing_plot_all <- (emp_timing_plot + plot_layout(guides = "collect") & theme(legend.position = "bottom")) +
  (emp_timing_10_inset + theme(legend.position = "none")) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")")


ggsave(filename = "4_empirical_timing_comparison_plot.png", plot = timing_plot_all, path = "./3_figures_and_tables/",
       width = 20, height = 10, units = "cm", dpi = "retina")


## Code for generating Figure 3 from the paper
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}


PSM_list <- list()
counter <- 1
for(i in c("1", "2", "3", "5", "10", "25", "50", "100", "150", "200")){
  linkage_file <- paste0("./2_empirical_analysis/timing_results/empirical_linkage_lambda_results_box_", i, "_area_200.csv")
  linkage_results <- fread(file = linkage_file, skip = 2501, header = FALSE) %>% as.matrix()
  linkage_results <- linkage_results + 1
  PSM_list[[counter]] <- comp.psm(linkage_results)
  counter <- counter + 1
}

area_200_cormat <- round(cor(do.call(cbind, lapply(PSM_list, c))), 4)
rownames(area_200_cormat) <- colnames(area_200_cormat) <- c("1", "2", "3", "5", "10", "25", "50", "100", "150", "200")
lower_tri_200 <- get_lower_tri(area_200_cormat)
area_200_cormat_melt <- melt(lower_tri_200, na.rm = TRUE)


## Create and save the correlation heatmap plot (Figure 3)
area_200_cormat_melt |> 
  ggplot(aes(x = as.factor(Var1), y = as.factor(Var2), fill = value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = value), size = 4) +
  scale_fill_gradient(low = cbbPalette[1], high = cbbPalette[2], 
                      limit = c(.95,1), space = "Lab", 
                      name="Correlation") +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.3, 0.6),
    legend.direction = "vertical",
    text = element_text(family = "serif", size = 14)) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5,
                               title.position = "top")) +
  labs(x = "Bounding Box Margin", y = "Bounding Box Margin",
       title = "",
       fill = "Correlation") -> cor_plot_200

ggsave(filename = "3_correlation_plot_200.png", plot = cor_plot_200, path = "./3_figures_and_tables/",
       width = 15, height = 15, units = "cm", dpi = "retina")
