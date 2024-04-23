## Load libraries and sampler functions ----
library(readr) ## load and save results
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

sourceCpp('./code/cpp_code/two_stage_func.cpp')
# cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#56B4E9", "#009E73", "#D55E00", "#F0E442")

## Define a function to add the filename to each line of the imported dataset
read_plus <- function(flnm) {
  read_csv(flnm, show_col_types = FALSE, skip = 2501, col_names = FALSE) %>% 
    mutate(filename = flnm)
}


## Two-stage growth model processing
## Two-stage linkage averaging

n_latents <- read_csv("./code/empirical_data/timing_results/N_latents_per_area.csv",
                      show_col_types = FALSE,
                      col_names = TRUE)

## Read in all of the results .csv's into one tibble with the filename attached to each line
emp_timing_table <-
  list.files(path = "./code/empirical_data/timing_results",
             pattern = "timing_results",
             full.names = T) |>
  mixedsort() |>
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./code/empirical_data/timing_results/", "", emp_timing_table$filename, fixed = TRUE)
myNames <- gsub(".csv", "", myNames, fixed = TRUE)

variables <- do.call(rbind, str_split(myNames, "_"))

# emp_timing_table <- emp_timing_table %>%
#   mutate(filename = myNames,
#          box_size = factor(variables[,6],
#                            levels = c("1", "2", "3", "5", "10", "25", "50", "100", "150", "200", "250", "300")),
#          area = variables[,8]) |> 
#   rename(time = X1) %>%
#   filter(!box_size == "all")
emp_timing_table <- emp_timing_table %>%
  mutate(filename = myNames,
         box_size = as.numeric(variables[,6]),
         area = variables[,8]) |> 
  rename(time = X1) %>%
  filter(!box_size == "all")

emp_timing_summary <- emp_timing_table |> group_by(box_size, area) |> summarise(mean_time = mean(time))

# emp_timing_summary |> 
#   ggplot(aes(x = box_size, y = mean_time, group = area, color = as.factor(area), shape = as.factor(area))) +
#   geom_point() +
#   geom_line() +
#   labs(x = "Bounding Box Margin",
#        y = "Average Time Per Iteration (Seconds)",
#        title = "",
#        color = TeX("Domain Size $(m^2)$")) +
#   scale_x_discrete(labels = c(TeX("$\\pm 2$"), TeX("$\\pm 5$"),
#                               TeX("$\\pm 10$"), TeX("$\\pm 25"),
#                               TeX("$\\pm 50$"), TeX("$\\pm 100"),
#                               TeX("$\\pm 150$"), TeX("$\\pm 200"),
#                               TeX("$\\pm 250$"), TeX("$\\pm 300"))) +
#   scale_colour_manual(values = cbbPalette) +
#   theme_bw() +
#   theme(
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.ticks = element_blank(),
#     legend.justification = c(1, 0),
#     legend.position = c(0.4, 0.7),
#     legend.direction = "vertical")


# emp_timing_table <- emp_timing_table %>%
#   mutate(filename = myNames,
#          box_size = as.numeric(variables[,6]),
#          area = variables[,8]) |> 
#   rename(time = X1) %>%
#   filter(!box_size == "all")
emp_timing_summary <- emp_timing_table |> group_by(box_size, area) |> summarise(mean_time = mean(time))
emp_timing_summary$n_latents <- NA
for(i in 1:nrow(emp_timing_summary)){
  if(emp_timing_summary$area[i] == "100"){
    emp_timing_summary$n_latents[i] <- n_latents$N[1]
  }else if(emp_timing_summary$area[i] == "150"){
    emp_timing_summary$n_latents[i] <- n_latents$N[2]
  }else if(emp_timing_summary$area[i] == "200"){
    emp_timing_summary$n_latents[i] <- n_latents$N[3]
  }else if(emp_timing_summary$area[i] == "300"){
    emp_timing_summary$n_latents[i] <- n_latents$N[4]
  }
}

emp_timing_summary <- emp_timing_summary |> mutate(scaled_mean_time = mean_time / n_latents)

emp_timing_summary |>
  ggplot(aes(x = box_size, y = mean_time, group = area, color = as.factor(area), shape = as.factor(area))) +
  geom_point() +
  geom_line() +
  labs(x = "Bounding Box Margin",
       y = "Average Time Per Iteration (Seconds)",
       title = "",
       color = TeX("Domain Size $(m^2)$"),
       shape = TeX("Domain Size $(m^2)$")) +
  # scale_x_discrete(labels = c(TeX("$\\pm 2$"), TeX("$\\pm 5$"),
  #                             TeX("$\\pm 10$"), TeX("$\\pm 25"),
  #                             TeX("$\\pm 50$"), TeX("$\\pm 100"),
  #                             TeX("$\\pm 150$"), TeX("$\\pm 200"))) +
  scale_colour_manual(values = cbbPalette) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    # legend.justification = c(1, 0),
    # legend.position = c(0.2, 0.7),
    # legend.direction = "vertical",
    legend.position = "bottom",
    text = element_text(family = "serif", size = 12)) -> emp_timing_plot

ggsave(filename = "emp_timing_plot.png", plot = emp_timing_plot, path = "./plots/F23/",
       width = 15, height = 15, units = "cm", dpi = "retina")


emp_timing_summary |> filter(box_size <= 25) |> 
  ggplot(aes(x = as.numeric(box_size), y = mean_time, group = area, color = as.factor(area), shape = as.factor(area))) +
  geom_point() +
  geom_line() +
  labs(x = "Bounding Box Margin",
       y = "Average Time Per Iteration (Seconds)",
       title = "",
       color = TeX("Domain Size $(m^2)$"),
       shape = TeX("Domain Size $(m^2)$")) +
  # scale_x_discrete(labels = c(TeX("$\\pm 2$"), TeX("$\\pm 5$"),
  #                             TeX("$\\pm 10$"), TeX("$\\pm 25"),
  #                             TeX("$\\pm 50$"), TeX("$\\pm 100"),
  #                             TeX("$\\pm 150$"), TeX("$\\pm 200"))) +
  scale_colour_manual(values = cbbPalette) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.2, 0.7),
    legend.direction = "vertical",
    text = element_text(family = "serif")) -> emp_timing_25_inset



emp_timing_summary |> filter(box_size <= 10) |> 
  ggplot(aes(x = as.numeric(box_size), y = mean_time, group = area, color = as.factor(area), shape = as.factor(area))) +
  geom_point() +
  geom_line() +
  labs(x = "Bounding Box Margin",
       y = "Average Time Per Iteration (Seconds)",
       title = "",
       color = TeX("Domain Size $(m^2)$"),
       shape = TeX("Domain Size $(m^2)$")) +
  # scale_x_discrete(labels = c(TeX("$\\pm 2$"), TeX("$\\pm 5$"),
  #                             TeX("$\\pm 10$"), TeX("$\\pm 25"),
  #                             TeX("$\\pm 50$"), TeX("$\\pm 100"),
  #                             TeX("$\\pm 150$"), TeX("$\\pm 200"))) +
  scale_colour_manual(values = cbbPalette) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    # legend.justification = c(1, 0),
    # legend.position = c(0.2, 0.7),
    # legend.direction = "vertical",
    text = element_text(family = "serif", size = 12)) -> emp_timing_10_inset



# emp_timing_summary |> 
#   ggplot(aes(x = box_size, y = scaled_mean_time, group = area, color = as.factor(area), shape = as.factor(area))) +
#   geom_point() +
#   geom_line() +
#   labs(x = "Bounding Box Margin",
#        y = "Average Time Per Iteration (Seconds)",
#        title = "",
#        color = TeX("Domain Size $(m^2)$"),
#        shape = TeX("Domain Size $(m^2)$")) +
#   # scale_x_discrete(labels = c(TeX("$\\pm 2$"), TeX("$\\pm 5$"),
#   #                             TeX("$\\pm 10$"), TeX("$\\pm 25"),
#   #                             TeX("$\\pm 50$"), TeX("$\\pm 100"),
#   #                             TeX("$\\pm 150$"), TeX("$\\pm 200"))) +
#   scale_colour_manual(values = cbbPalette) +
#   ylim(0, 5) +
#   theme_bw() +
#   theme(
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.ticks = element_blank(),
#     legend.justification = c(1, 0),
#     legend.position = c(0.4, 0.7),
#     legend.direction = "vertical",
#     text = element_text(family = "serif", size = 14))


timing_plot_all <- (emp_timing_plot + plot_layout(guides = "collect") & theme(legend.position = "bottom")) +
  (emp_timing_10_inset + theme(legend.position = "none")) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")")


ggsave(filename = "emp_timing_plot_all.png", plot = timing_plot_all, path = "./plots/F23/",
       width = 20, height = 10, units = "cm", dpi = "retina")


library(reshape2)
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

PSM_list <- list()
counter <- 1
for(i in c("1", "2", "3", "5", "10", "25", "50", "100")){
  linkage_file <- paste0("./code/empirical_data/timing_results/empirical_linkage_lambda_results_box_", i, "_area_100.csv")
  # linkage_results <- fread(file = linkage_file, header = TRUE) %>% as.matrix()
  linkage_results <- fread(file = linkage_file, skip = 2501, header = FALSE) %>% as.matrix()
  linkage_results <- linkage_results + 1
  PSM_list[[counter]] <- comp.psm(linkage_results)
  counter <- counter + 1
}

area_100_cormat <- round(cor(do.call(cbind, lapply(PSM_list, c))), 4)
rownames(area_100_cormat) <- colnames(area_100_cormat) <- c("1", "2", "3", "5", "10", "25", "50", "100")
lower_tri_100 <- get_lower_tri(area_100_cormat)
area_100_cormat_melt <- melt(lower_tri_100, na.rm = TRUE)

area_100_cormat_melt |> 
  ggplot(aes(x = as.factor(Var1), y = as.factor(Var2), fill = value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = value), size = 4) +
  scale_fill_gradient(low = cbbPalette[1], high = cbbPalette[2], 
                      limit = c(.95,1), space = "Lab", 
                      name="Correlation") +
  theme_bw() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.3, 0.6),
    legend.direction = "vertical",
    text = element_text(family = "serif")) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5,
                               title.position = "top")) +
  labs(x = "", y = "",
       title = "",
       fill = "Correlation") -> cor_plot_100


PSM_list <- list()
counter <- 1
for(i in c("1", "2", "3", "5", "10", "25", "50", "100", "150")){
  linkage_file <- paste0("./code/empirical_data/timing_results/empirical_linkage_lambda_results_box_", i, "_area_150.csv")
  # linkage_results <- fread(file = linkage_file, header = TRUE) %>% as.matrix()
  linkage_results <- fread(file = linkage_file, skip = 2501, header = FALSE) %>% as.matrix()
  linkage_results <- linkage_results + 1
  PSM_list[[counter]] <- comp.psm(linkage_results)
  counter <- counter + 1
}

area_150_cormat <- round(cor(do.call(cbind, lapply(PSM_list, c))), 4)
rownames(area_150_cormat) <- colnames(area_150_cormat) <- c("1", "2", "3", "5", "10", "25", "50", "100", "150")
lower_tri_150 <- get_lower_tri(area_150_cormat)
area_150_cormat_melt <- melt(lower_tri_150, na.rm = TRUE)

area_150_cormat_melt |> 
  ggplot(aes(x = as.factor(Var1), y = as.factor(Var2), fill = value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = value), size = 4) +
  scale_fill_gradient(low = cbbPalette[1], high = cbbPalette[2], 
                      limit = c(.95,1), space = "Lab", 
                      name="Correlation") +
  theme_bw() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.3, 0.6),
    # legend.position = "none",
    legend.direction = "vertical",
    text = element_text(family = "serif")) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5,
                               title.position = "top")) +
  labs(x = "", y = "",
       title = "",
       fill = "Correlation") -> cor_plot_150





PSM_list <- list()
counter <- 1
for(i in c("1", "2", "3", "5", "10", "25", "50", "100", "150", "200")){
  linkage_file <- paste0("./code/empirical_data/timing_results/empirical_linkage_lambda_results_box_", i, "_area_200.csv")
  # linkage_results <- fread(file = linkage_file, header = TRUE) %>% as.matrix()
  linkage_results <- fread(file = linkage_file, skip = 2501, header = FALSE) %>% as.matrix()
  linkage_results <- linkage_results + 1
  PSM_list[[counter]] <- comp.psm(linkage_results)
  counter <- counter + 1
}

area_200_cormat <- round(cor(do.call(cbind, lapply(PSM_list, c))), 4)
rownames(area_200_cormat) <- colnames(area_200_cormat) <- c("1", "2", "3", "5", "10", "25", "50", "100", "150", "200")
lower_tri_200 <- get_lower_tri(area_200_cormat)
area_200_cormat_melt <- melt(lower_tri_200, na.rm = TRUE)

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
    # legend.position = "none",
    legend.direction = "vertical",
    text = element_text(family = "serif", size = 14)) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5,
                               title.position = "top")) +
  labs(x = "Bounding Box Margin", y = "Bounding Box Margin",
       title = "",
       fill = "Correlation") -> cor_plot_200

ggsave(filename = "cor_plot_200.png", plot = cor_plot_200, path = "./plots/F23/",
       width = 15, height = 15, units = "cm", dpi = "retina")

# ((cor_plot_100 + cor_plot_150) / (cor_plot_200 + plot_layout(widths = c(2, 2))) ) + plot_annotation(tag_levels = 'a', tag_suffix = ")")
