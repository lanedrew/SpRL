######################################################################
#### This script generates the basis for Table 3 from the paper. #####
######################################################################

## Load the relevant packages
library(readr)
library(dplyr)
library(gtools)
library(purrr)
library(stringr)
library(tidyr)


## Define a function to read in the data
read_plus <- function(flnm) {
  read_csv(flnm, show_col_types = FALSE) %>%
    mutate(filename = flnm)
}


## True Linkage Processing
## Read in the data
sprl_growth_sim_table <-
  list.files(path = "./1_simulation_study/simulation_results/true_linkage",
             pattern = "1_alpha",
             full.names = TRUE) %>%
  mixedsort() %>%
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./1_simulation_study/simulation_results/true_linkage/", "", sprl_growth_sim_table$filename, fixed = TRUE)
myNames <- gsub("_growth_results_TL.csv", "", myNames, fixed = TRUE)
variables <- do.call(rbind, str_split(myNames, "_"))


sprl_growth_sim_table <- sprl_growth_sim_table %>%
  mutate(filename = myNames,
         density = variables[,1],
         noise = variables[,3],
         true_alpha = variables[,5],
         index = as.numeric(variables[,7])) %>%
  rename(beta_1 = V4,
         beta_2 = V5,
         beta_3 = V6,
         beta_4 = V7)


growth_summary_table <- sprl_growth_sim_table %>%
  group_by(density, noise, true_alpha, index) %>%
  summarise(gamma_est = mean(gamma), alpha_est = mean(alpha), tau_est = mean(tau), beta_0_est = mean(beta_0),
            beta_1_est = mean(beta_1), beta_2_est = mean(beta_2), beta_3_est = mean(beta_3), beta_4_est = mean(beta_4),
            gamma_q5 = quantile(gamma, c(.05)), alpha_q5 = quantile(alpha, c(.05)), tau_q5 = quantile(tau, c(.05)), beta_0_q5 = quantile(beta_0, c(.05)),
            beta_1_q5 = quantile(beta_1, c(.05)), beta_2_q5 = quantile(beta_2, c(.05)), beta_3_q5 = quantile(beta_3, c(.05)), beta_4_q5 = quantile(beta_4, c(.05)),
            gamma_q95 = quantile(gamma, c(.95)), alpha_q95 = quantile(alpha, c(.95)), tau_q95 = quantile(tau, c(.95)), beta_0_q95 = quantile(beta_0, c(.95)),
            beta_1_q95 = quantile(beta_1, c(.95)), beta_2_q95 = quantile(beta_2, c(.95)), beta_3_q95 = quantile(beta_3, c(.95)), beta_4_q95 = quantile(beta_4, c(.95)))


test_table <- growth_summary_table %>%
  pivot_longer(cols = gamma_est:beta_4_q95,
               names_to = c(".value", "estimate_type"),
               names_pattern = "(.*)_(.*)") %>%
  pivot_longer(cols = gamma:beta_4,
               names_to = "parameter") %>%
  pivot_wider(names_from = estimate_type,
              values_from = value) %>%
  mutate(coverage = NA)


coverage_func <- function(parameter, true_value, lower_bound, upper_bound){
  coverage_vector <- rep(NA, length(parameter))
  for(i in 1:length(parameter)){
    if(parameter[i] == "gamma"){
      coverage_vector[i] <- ifelse(lower_bound[i] < true_value[1] & true_value[1] < upper_bound[i], 1, 0)
    }else if(parameter[i] == "alpha"){
      coverage_vector[i] <- ifelse(lower_bound[i] < true_value[2] & true_value[2] < upper_bound[i], 1, 0)
    }else if(parameter[i] == "tau"){
      coverage_vector[i] <- ifelse(lower_bound[i] < true_value[3] & true_value[3] < upper_bound[i], 1, 0)
    }else if(parameter[i] == "beta_0"){
      coverage_vector[i] <- ifelse(lower_bound[i] < true_value[4] & true_value[4] < upper_bound[i], 1, 0)
    }else if(parameter[i] == "beta_1"){
      coverage_vector[i] <- ifelse(lower_bound[i] < true_value[5] & true_value[5] < upper_bound[i], 1, 0)
    }else if(parameter[i] == "beta_2"){
      coverage_vector[i] <- ifelse(lower_bound[i] < true_value[6] & true_value[6] < upper_bound[i], 1, 0)
    }else if(parameter[i] == "beta_3"){
      coverage_vector[i] <- ifelse(lower_bound[i] < true_value[7] & true_value[7] < upper_bound[i], 1, 0)
    }else if(parameter[i] == "beta_4"){
      coverage_vector[i] <- ifelse(lower_bound[i] < true_value[8] & true_value[8] < upper_bound[i], 1, 0)
    }
  }
  return(coverage_vector)
}


test_table$coverage <- coverage_func(test_table$parameter, true_value = c(12, 1, .5, 3, .5, -.5, .5, -.5),
                                     lower_bound = test_table$q5, upper_bound = test_table$q95)


TL_cov_table <- test_table %>%
  group_by(density, noise, true_alpha, parameter) %>%
  summarise(cov_prob = mean(coverage)) %>%
  mutate(linkage = "TL")

print(paste0("TL processing complete"))

## NDM Processing
## Read in the data
sprl_growth_sim_table <-
  list.files(path = "./1_simulation_study/simulation_results/NDM/growth_model_estimates",
             pattern = "1_alpha",
             full.names = TRUE) %>%
  mixedsort() %>%
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./1_simulation_study/simulation_results/NDM/growth_model_estimates/", "", sprl_growth_sim_table$filename, fixed = TRUE)
myNames <- gsub("_growth_results_NDM.csv", "", myNames, fixed = TRUE)
variables <- do.call(rbind, str_split(myNames, "_"))

sprl_growth_sim_table <- sprl_growth_sim_table %>%
  mutate(filename = myNames,
         density = variables[,1],
         noise = variables[,3],
         true_alpha = variables[,5],
         index = as.numeric(variables[,8])) %>%
  rename(beta_1 = beta.1,
         beta_2 = beta.2,
         beta_3 = beta.3,
         beta_4 = beta.4)


growth_summary_table <- sprl_growth_sim_table %>%
  group_by(density, noise, true_alpha, index) %>%
  summarise(gamma_est = mean(gamma), alpha_est = mean(alpha), tau_est = mean(tau), beta_0_est = mean(beta_0),
            beta_1_est = mean(beta_1), beta_2_est = mean(beta_2), beta_3_est = mean(beta_3), beta_4_est = mean(beta_4),
            gamma_q5 = quantile(gamma, c(.05)), alpha_q5 = quantile(alpha, c(.05)), tau_q5 = quantile(tau, c(.05)), beta_0_q5 = quantile(beta_0, c(.05)),
            beta_1_q5 = quantile(beta_1, c(.05)), beta_2_q5 = quantile(beta_2, c(.05)), beta_3_q5 = quantile(beta_3, c(.05)), beta_4_q5 = quantile(beta_4, c(.05)),
            gamma_q95 = quantile(gamma, c(.95)), alpha_q95 = quantile(alpha, c(.95)), tau_q95 = quantile(tau, c(.95)), beta_0_q95 = quantile(beta_0, c(.95)),
            beta_1_q95 = quantile(beta_1, c(.95)), beta_2_q95 = quantile(beta_2, c(.95)), beta_3_q95 = quantile(beta_3, c(.95)), beta_4_q95 = quantile(beta_4, c(.95)))


test_table <- growth_summary_table %>%
  pivot_longer(cols = gamma_est:beta_4_q95,
               names_to = c(".value", "estimate_type"),
               names_pattern = "(.*)_(.*)") %>%
  pivot_longer(cols = gamma:beta_4,
               names_to = "parameter") %>%
  pivot_wider(names_from = estimate_type,
              values_from = value) %>%
  mutate(coverage = NA)


test_table$coverage <- coverage_func(test_table$parameter, true_value = c(12, 1, .5, 3, .5, -.5, .5, -.5),
                                     lower_bound = test_table$q5, upper_bound = test_table$q95)


NDM_cov_table <- test_table %>%
  group_by(density, noise, true_alpha, parameter) %>%
  summarise(cov_prob = mean(coverage)) %>%
  mutate(linkage = "NDM")

print(paste0("NDM processing complete"))

## Two-stage growth model processing
## Define a function to read in the data
read_plus <- function(flnm) {
  read_csv(flnm, show_col_types = FALSE) %>%
    mutate(filename = flnm)
}


## Read in the data
sprl_growth_sim_table <-
  list.files(path = "./1_simulation_study/simulation_results/two_stage/linkage_avg",
             pattern = "1_alpha",
             full.names = TRUE) %>%
  mixedsort() %>%
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./1_simulation_study/simulation_results/two_stage/linkage_avg/", "", sprl_growth_sim_table$filename, fixed = TRUE)
myNames <- gsub("_growth_results_linkage_avg.csv", "", myNames, fixed = TRUE)
variables <- do.call(rbind, str_split(myNames, "_"))

sprl_growth_sim_table <- sprl_growth_sim_table %>%
  mutate(filename = myNames,
         density = variables[,1],
         noise = variables[,3],
         true_alpha = variables[,5],
         index = as.numeric(variables[,8])) %>%
  rename(beta_1 = V4,
         beta_2 = V5,
         beta_3 = V6,
         beta_4 = V7)


growth_summary_table <- sprl_growth_sim_table %>%
  group_by(density, noise, true_alpha, index) %>%
  summarise(gamma_est = mean(gamma), alpha_est = mean(alpha), tau_est = mean(tau), beta_0_est = mean(beta_0),
            beta_1_est = mean(beta_1), beta_2_est = mean(beta_2), beta_3_est = mean(beta_3), beta_4_est = mean(beta_4),
            gamma_q5 = quantile(gamma, c(.05)), alpha_q5 = quantile(alpha, c(.05)), tau_q5 = quantile(tau, c(.05)), beta_0_q5 = quantile(beta_0, c(.05)),
            beta_1_q5 = quantile(beta_1, c(.05)), beta_2_q5 = quantile(beta_2, c(.05)), beta_3_q5 = quantile(beta_3, c(.05)), beta_4_q5 = quantile(beta_4, c(.05)),
            gamma_q95 = quantile(gamma, c(.95)), alpha_q95 = quantile(alpha, c(.95)), tau_q95 = quantile(tau, c(.95)), beta_0_q95 = quantile(beta_0, c(.95)),
            beta_1_q95 = quantile(beta_1, c(.95)), beta_2_q95 = quantile(beta_2, c(.95)), beta_3_q95 = quantile(beta_3, c(.95)), beta_4_q95 = quantile(beta_4, c(.95)))


test_table <- growth_summary_table %>%
  pivot_longer(cols = gamma_est:beta_4_q95,
               names_to = c(".value", "estimate_type"),
               names_pattern = "(.*)_(.*)") %>%
  pivot_longer(cols = gamma:beta_4,
               names_to = "parameter") %>%
  pivot_wider(names_from = estimate_type,
              values_from = value) %>%
  mutate(coverage = NA)


test_table$coverage <- coverage_func(test_table$parameter, true_value = c(12, 1, .5, 3, .5, -.5, .5, -.5),
                                     lower_bound = test_table$q5, upper_bound = test_table$q95)


LA_cov_table <- test_table %>%
  group_by(density, noise, true_alpha, parameter) %>%
  summarise(cov_prob = mean(coverage)) %>%
  mutate(linkage = "LA")

print(paste0("LA processing complete"))

## Create a table of all results
growth_cov_table <- as.data.frame(rbind(TL_cov_table, NDM_cov_table, LA_cov_table))
growth_cov_table_mod <- growth_cov_table |>
  pivot_wider(names_from = c("parameter"),
              values_from = c("cov_prob")) |>
  mutate(noise = factor(noise,
                        levels = c("small", "medium", "large"),
                        labels = c("Small", "Medium", "Large")),
         density = factor(density,
                          levels = c("low", "med", "high"),
                          labels = c("Low", "Medium", "High")),
         linkage = factor(linkage,
                          levels = c("TL", "LA", "NDM")))


## Generate the LaTeX for Table 3
growth_cov_table_mod |>
  arrange(density, noise, true_alpha, linkage) |>
  select(-c(true_alpha)) |>
  kbl(format = "latex", booktabs = TRUE, align = "lllcccccccc", linesep = "", escape = FALSE, digits = 4,
      col.names = c("Density", "Noise", "Linkage Approach", "$\\alpha$", "$\\beta_0$", "$\\beta_1$", "$\\beta_2$", "$\\beta_3$", "$\\beta_4$", "$\\gamma$", "$\\tau^2$"),
      label = "growth_cov_table", caption = c("")) |>
  kable_classic() |>
  collapse_rows(columns = 1:3) |>
  add_header_above(header = c(" " = 3, "Empirical Coverage by Parameter" = 8), align = "c",
                   bold = TRUE)
