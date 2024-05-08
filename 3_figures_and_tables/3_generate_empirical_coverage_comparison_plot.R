########################################################
#### This script generates Figure 5 from the paper. ####
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
library(ggstance)

## Set a colorblind friendly palate to use for visualizations
cbbPalette <- c("#56B4E9", "#009E73", "#D55E00")

## Define a function to add the filename to each line of the imported dataset
read_plus <- function(flnm) {
  read_csv(flnm, show_col_types = FALSE) %>% 
    mutate(filename = flnm)
}


## Two-stage linkage averaging growth results
## Read in all of the results .csv's into one tibble with the filename attached to each line
emp_res_table <-
  list.files(path = "./2_empirical_analysis/model_results/growth_model/LA",
             pattern = "model_skew_t",
             full.names = T) |>
  mixedsort() |>
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./2_empirical_analysis/model_results/growth_model/LA/", "", emp_res_table$filename, fixed = TRUE)
myNames <- gsub(".csv", "", myNames, fixed = TRUE)
variables <- do.call(rbind, str_split(myNames, "_"))

emp_res_table <- emp_res_table %>%
  mutate(filename = myNames,
         index = as.numeric(variables[,15])) %>%
  select(-"lp__") %>%
  rename(delta = lambda,
         omega = q)

emp_mat <- emp_res_table |> select(-c("filename", "index")) |> as.matrix()

emp_means <- apply(emp_mat, 2, mean)
emp_quants <- apply(emp_mat, 2, quantile, probs = c(.05, .95))

emp_all <- rbind(emp_means, emp_quants)

colnames(emp_all) <- c("gamma", "tau", "beta_0",
                       "TWI",
                       "FA",
                       "SP",
                       "GDD",
                       "TWI_SP",
                       "TWI_GDD",
                       "LNV",
                       "RSI",
                       "ND",
                       "alpha",
                       "delta",
                       "omega")

emp_proc <- data.frame(t(emp_all)) |>
  mutate(variable = colnames(emp_all)) |>
  mutate(variable = factor(variable, levels = rev(c("gamma", "tau", "beta_0", "TWI",
                                                    "FA",
                                                    "SP",
                                                    "GDD",
                                                    "TWI_SP",
                                                    "TWI_GDD",
                                                    "LNV",
                                                    "RSI",
                                                    "ND",
                                                    "alpha",
                                                    "delta",
                                                    "omega"))))


## Near distance matching growth results
emp_res_table_ndm <-
  list.files(path = "./2_empirical_analysis/model_results/growth_model/NDM",
             pattern = "model_skew_t",
             full.names = T) |>
  mixedsort() |>
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./2_empirical_analysis/model_results/growth_model/NDM/", "", emp_res_table_ndm$filename, fixed = TRUE)
myNames <- gsub(".csv", "", myNames, fixed = TRUE)
variables <- do.call(rbind, str_split(myNames, "_"))

emp_res_table_ndm <- emp_res_table_ndm %>%
  mutate(filename = myNames)

emp_mat_ndm <- emp_res_table_ndm |> select(-c("filename")) |> as.matrix()

emp_means_ndm <- apply(emp_mat_ndm, 2, mean)
emp_quants_ndm <- apply(emp_mat_ndm, 2, quantile, probs = c(.05, .95))

emp_all_ndm <- rbind(emp_means_ndm, emp_quants_ndm)
colnames(emp_all_ndm) <- c("gamma", "tau", "beta_0",
                           "TWI",
                           "FA",
                           "SP",
                           "GDD",
                           "TWI_SP",
                           "TWI_GDD",
                           "LNV",
                           "RSI",
                           "ND",
                           "alpha",
                           "delta",
                           "omega")

emp_ndm_proc <- data.frame(t(emp_all_ndm)) |>
  mutate(variable = colnames(emp_all_ndm)) |>
  mutate(variable = factor(variable, levels = rev(c("gamma", "tau", "beta_0", "TWI",
                                                    "FA",
                                                    "SP",
                                                    "GDD",
                                                    "TWI_SP",
                                                    "TWI_GDD",
                                                    "LNV",
                                                    "RSI",
                                                    "ND",
                                                    "alpha",
                                                    "delta",
                                                    "omega")))) |> 
  rename(emp_means = emp_means_ndm)


## Polygon overlap matching growth results
emp_res_table_pom <-
  list.files(path = "./2_empirical_analysis/model_results/growth_model/POM",
             pattern = "model_skew_t",
             full.names = T) |>
  mixedsort() |>
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./2_empirical_analysis/model_results/growth_model/POM/", "", emp_res_table_pom$filename, fixed = TRUE)
myNames <- gsub(".csv", "", myNames, fixed = TRUE)
variables <- do.call(rbind, str_split(myNames, "_"))

emp_res_table_pom <- emp_res_table_pom %>%
  mutate(filename = myNames) 

emp_mat_pom <- emp_res_table_pom |> select(-c("filename")) |> as.matrix()

emp_means_pom<- apply(emp_mat_pom, 2, mean)
emp_quants_pom <- apply(emp_mat_pom, 2, quantile, probs = c(.05, .95))

emp_all_pom <- rbind(emp_means_pom, emp_quants_pom)
colnames(emp_all_pom) <- c("gamma", "tau", "beta_0",
                           "TWI",
                           "FA",
                           "SP",
                           "GDD",
                           "TWI_SP",
                           "TWI_GDD",
                           "LNV",
                           "RSI",
                           "ND",
                           "alpha",
                           "delta",
                           "omega")

emp_pom_proc <- data.frame(t(emp_all_pom)) |>
  mutate(variable = colnames(emp_all_pom)) |>
  mutate(variable = factor(variable, levels = rev(c("gamma", "tau", "beta_0", "TWI",
                                                    "FA",
                                                    "SP",
                                                    "GDD",
                                                    "TWI_SP",
                                                    "TWI_GDD",
                                                    "LNV",
                                                    "RSI",
                                                    "ND",
                                                    "alpha",
                                                    "delta",
                                                    "omega")))) |> 
  rename(emp_means = emp_means_pom)


# Combined results
emp_full <- data.frame(rbind(emp_proc, emp_ndm_proc, emp_pom_proc), linkage = rep(c("LA", "NDM", "POM"), each = nrow(emp_proc))) |> 
  mutate(linkage = factor(linkage, levels = c("LA", "NDM", "POM"),
                          labels = c("Linkage-Averaging",
                                     "Nearest Distance Matching",
                                     "Polygon Overlap Matching")))

## Plot and save the empirical coverage comparison plot (Figure 5)
emp_full |>
  filter(!variable %in% c("gamma", "tau", "alpha", "delta", "omega")) |>
  ggplot(aes(x = variable, color = linkage, ymin = X5., ymax = X95.)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .5, color = "black") +
  geom_linerange(position = position_dodge(.5)) +
  geom_point(aes(x = variable, y = emp_means, shape = as.factor(linkage)), position = position_dodge(width = .5, preserve = "total")) +
  scale_color_manual(values = cbbPalette) +
  scale_x_discrete(limits = rev,
                   labels = (c(TeX("$\\beta_0$"), TeX("TWI"), TeX("FA"),
                               TeX("SP"), TeX("GDD"),
                               TeX("SP$\\times$TWI"), TeX("GDD$\\times$TWI"),
                               TeX("LNV"), TeX("RSI"),
                               TeX("ND")))) +
  labs(x = "", y = "", title = "", color = "Linkage Approach", shape = "Linkage Approach") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "serif", size = 14)) -> coverage_comp_plot_covars


ggsave(filename = "5_coverage_comparison_plot.png", plot = coverage_comp_plot_covars, path = "./3_figures_and_tables/",
       width = 20, height = 12, units = "cm", dpi = "retina")