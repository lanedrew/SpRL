
###########################################################################################
#### This script processes all of the results from the joint and two-stage SpRL models ####
###########################################################################################

## Load the relevant packages ##
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
library(ggstance)

cbbPalette <- c("#56B4E9", "#009E73", "#D55E00")


## Define a function to add the filename to each line of the imported dataset
read_plus <- function(flnm) {
  read_csv(flnm, show_col_types = FALSE) %>% 
    mutate(filename = flnm)
}


## Two-stage linkage averaging growth results

## Read in all of the results .csv's into one tibble with the filename attached to each line
emp_res_table <-
  list.files(path = "./code/empirical_data",
             pattern = "emp_pooled_LA_covars_all_growth_cutoff_100_index_",
             full.names = T) |>
  mixedsort() |>
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./code/empirical_data/", "", emp_res_table$filename, fixed = TRUE)
myNames <- gsub(".csv", "", myNames, fixed = TRUE)

variables <- do.call(rbind, str_split(myNames, "_"))

emp_res_table <- emp_res_table %>%
  mutate(filename = myNames,
         index = as.numeric(variables[,10])) %>%
  select(-"lp__")

emp_mat <- emp_res_table |> select(-c("filename", "index")) |> as.matrix()

emp_means <- apply(emp_mat, 2, mean)
emp_quants <- apply(emp_mat, 2, quantile, probs = c(.05, .95))

emp_all <- rbind(emp_means, emp_quants)

colnames(emp_all) <- c("b", "tau", "beta_0",
                       "TWI",
                       "southness",
                       "SPP",
                       "GDD",
                       "TWI_SPP",
                       "TWI_GDD",
                       "LNV",
                       "RSI",
                       "ND",
                       "alpha")

emp_proc <- data.frame(t(emp_all)) |>
  mutate(variable = colnames(emp_all)) |>
  mutate(variable = factor(variable, levels = rev(c("b", "tau", "beta_0", "TWI",
                                                    "southness",
                                                    "SPP",
                                                    "GDD",
                                                    "TWI_SPP",
                                                    "TWI_GDD",
                                                    "LNV",
                                                    "RSI",
                                                    "ND",
                                                    "alpha"))),
         coverage = NA)


for(i in 1:nrow(emp_proc)){
  if(0 > emp_proc$X5.[i] & 0 < emp_proc$X95.[i]){
    emp_proc$coverage[i] <- 0
  }else{
    emp_proc$coverage[i] <- 1
  }
}




# emp_proc |> 
#   filter(variable != "b") |>
#   ggplot(aes(color = as.factor(coverage))) +
#   geom_segment(aes(x = X5., xend = X95., y = variable, yend = variable)) +
#   geom_point(aes(x = emp_means, y = variable, shape = as.factor(coverage))) +
#   geom_vline(xintercept = 0, linetype = "dashed", alpha = .5) +
#   scale_color_manual(values = cbbPalette) +
#   labs(x = "", y = "Variable", title = "90% CIs for Growth Parameters with 100% Mortality Theshold") +
#   theme_classic() +
#   theme(legend.position = "none")




## Near distance matching growth results

emp_res_table_ndm <-
  list.files(path = "./code/naive_growth_models",
             pattern = "near_distance_matching_covars_all_growth_cutoff_100.csv",
             full.names = T) |>
  mixedsort() |>
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./code/naive_growth_models/", "", emp_res_table_ndm$filename, fixed = TRUE)
myNames <- gsub(".csv", "", myNames, fixed = TRUE)

variables <- do.call(rbind, str_split(myNames, "_"))

emp_res_table_ndm <- emp_res_table_ndm %>%
  mutate(filename = myNames) %>%
  select(-"lp__")

emp_mat_ndm <- emp_res_table_ndm |> select(-c("filename")) |> as.matrix()

emp_means_ndm <- apply(emp_mat_ndm, 2, mean)
emp_quants_ndm <- apply(emp_mat_ndm, 2, quantile, probs = c(.05, .95))

emp_all_ndm <- rbind(emp_means_ndm, emp_quants_ndm)
colnames(emp_all_ndm) <- c("b", "tau", "beta_0",
                           "TWI",
                           "southness",
                           "SPP",
                           "GDD",
                           "TWI_SPP",
                           "TWI_GDD",
                           "LNV",
                           "RSI",
                           "ND",
                           "alpha")

emp_ndm_proc <- data.frame(t(emp_all_ndm)) |>
  mutate(variable = colnames(emp_all_ndm)) |>
  mutate(variable = factor(variable, levels = rev(c("b", "tau", "beta_0", "TWI",
                                                    "southness",
                                                    "SPP",
                                                    "GDD",
                                                    "TWI_SPP",
                                                    "TWI_GDD",
                                                    "LNV",
                                                    "RSI",
                                                    "ND",
                                                    "alpha"))),
         coverage = NA) |> 
  rename(emp_means = emp_means_ndm)


for(i in 1:nrow(emp_ndm_proc)){
  if(0 > emp_ndm_proc$X5.[i] & 0 < emp_ndm_proc$X95.[i]){
    emp_ndm_proc$coverage[i] <- 0
  }else{
    emp_ndm_proc$coverage[i] <- 1
  }
}




# emp_ndm_proc |> 
#   filter(variable != "b") |>
#   ggplot(aes(color = as.factor(coverage))) +
#   geom_segment(aes(x = X5., xend = X95., y = variable, yend = variable)) +
#   geom_point(aes(x = emp_means, y = variable, shape = as.factor(coverage))) +
#   geom_vline(xintercept = 0, linetype = "dashed", alpha = .5) +
#   scale_color_manual(values = cbbPalette) +
#   labs(x = "", y = "Variable", title = "90% CIs for Growth Parameters with 100% Mortality Theshold") +
#   theme_classic() +
#   theme(legend.position = "none")


  

## Polygon overlap matching growth results

emp_res_table_pom <-
  list.files(path = "./code/naive_growth_models",
             pattern = "polygon_overlap_matching_2015_covars_all_growth_cutoff_100.csv",
             full.names = T) |>
  mixedsort() |>
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./code/naive_growth_models/", "", emp_res_table_pom$filename, fixed = TRUE)
myNames <- gsub(".csv", "", myNames, fixed = TRUE)

variables <- do.call(rbind, str_split(myNames, "_"))

emp_res_table_pom <- emp_res_table_pom %>%
  mutate(filename = myNames) %>%
  select(-"lp__")

emp_mat_pom <- emp_res_table_pom |> select(-c("filename")) |> as.matrix()

emp_means_pom<- apply(emp_mat_pom, 2, mean)
emp_quants_pom <- apply(emp_mat_pom, 2, quantile, probs = c(.05, .95))

emp_all_pom <- rbind(emp_means_pom, emp_quants_pom)
colnames(emp_all_pom) <- c("b", "tau", "beta_0",
                           "TWI",
                           "southness",
                           "SPP",
                           "GDD",
                           "TWI_SPP",
                           "TWI_GDD",
                           "LNV",
                           "RSI",
                           "ND",
                           "alpha")

emp_pom_proc <- data.frame(t(emp_all_pom)) |>
  mutate(variable = colnames(emp_all_pom)) |>
  mutate(variable = factor(variable, levels = rev(c("b", "tau", "beta_0", "TWI",
                                                    "southness",
                                                    "SPP",
                                                    "GDD",
                                                    "TWI_SPP",
                                                    "TWI_GDD",
                                                    "LNV",
                                                    "RSI",
                                                    "ND",
                                                    "alpha"))),
         coverage = NA) |> 
  rename(emp_means = emp_means_pom)


for(i in 1:nrow(emp_pom_proc)){
  if(0 > emp_pom_proc$X5.[i] & 0 < emp_pom_proc$X95.[i]){
    emp_pom_proc$coverage[i] <- 0
  }else{
    emp_pom_proc$coverage[i] <- 1
  }
}




# emp_pom_proc |> 
#   filter(variable != "b") |>
#   ggplot(aes(color = as.factor(coverage))) +
#   geom_segment(aes(x = X5., xend = X95., y = variable, yend = variable)) +
#   geom_point(aes(x = emp_means, y = variable, shape = as.factor(coverage))) +
#   geom_vline(xintercept = 0, linetype = "dashed", alpha = .5) +
#   scale_color_manual(values = cbbPalette) +
#   labs(x = "", y = "Variable", title = "90% CIs for Growth Parameters with 100% Mortality Theshold") +
#   theme_classic() +
#   theme(legend.position = "none")





# Combined results

emp_full <- data.frame(rbind(emp_proc, emp_ndm_proc, emp_pom_proc), linkage = rep(c("LA", "NDM", "POM"), each = nrow(emp_proc))) |> 
  mutate(linkage = factor(linkage, levels = c("LA", "NDM", "POM"),
                          labels = c("Linkage Averaging",
                                     "Near Distance Matching",
                                     "Polygon Overlap Matching")))

emp_full |>
  filter(variable != "b") |>
  ggplot(aes(y = variable, color = linkage, xmin = X5., xmax = X95.)) +
  geom_linerangeh(position = position_dodge(.5)) +
  geom_point(aes(x = emp_means, y = variable, shape = as.factor(linkage)), position = ggstance::position_dodgev(height = .5, preserve = "total")) +
  # geom_vline(xintercept = 0, linetype = "dashed", alpha = .5) +
  scale_color_manual(values = cbbPalette) +
  scale_y_discrete(labels = rev(c(TeX("$\\tau$"), TeX("$\\beta_0$"),
                              TeX("TWI"), TeX("Southness"),
                              TeX("SPP"), TeX("GDD"),
                              TeX("SPP$\\times$TWI"), TeX("GDD$\\times$TWI"),
                              TeX("LNV"), TeX("RSI"),
                              TeX("ND"), TeX("\\alpha")))) +
  labs(x = "", y = "", title = "", color = "Linkage Procedure", shape = "Linkage Procedure") +
  theme_bw() +
  theme(
    legend.justification = c(1, 0),
    legend.position = c(0.8, 0.3),
    legend.direction = "vertical",
    text = element_text(family = "serif")) -> coverage_comp_plot

ggsave(filename = "coverage_comp_plot.png", plot = coverage_comp_plot, path = "./plots/F23/",
       width = 15, height = 15, units = "cm", dpi = "retina")




emp_full |>
  filter(!variable %in% c("b", "tau", "beta_0", "alpha")) |>
  ggplot(aes(x = variable, color = linkage, ymin = X5., ymax = X95.)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .5, color = "black") +
  geom_linerange(position = position_dodge(.5)) +
  geom_point(aes(x = variable, y = emp_means, shape = as.factor(linkage)), position = position_dodge(width = .5, preserve = "total")) +
  scale_color_manual(values = cbbPalette) +
  scale_x_discrete(limits = rev,
                   labels = (c(TeX("TWI"), TeX("Southness"),
                                  TeX("SPP"), TeX("GDD"),
                                  TeX("SPP$\\times$TWI"), TeX("GDD$\\times$TWI"),
                                  TeX("LNV"), TeX("RSI"),
                                  TeX("ND")))) +
  labs(x = "", y = "", title = "", color = "Linkage Procedure", shape = "Linkage Procedure") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "serif", size = 14)) -> coverage_comp_plot_covars


ggsave(filename = "coverage_comp_plot_covars.png", plot = coverage_comp_plot_covars, path = "./plots/F23/",
       width = 20, height = 15, units = "cm", dpi = "retina")

