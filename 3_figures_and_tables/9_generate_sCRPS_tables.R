####################################################################################
#### This script pools the sCRPS scores from the LA model runs into a single   #####
#### pooled value and generates Tables 1 and 2 from the paper supplement.      #####
####################################################################################

## Load the necessary libraries
library(readr)
library(dplyr)
library(purrr)
library(loo)
library(tidyr)
library(kableExtra)

## Read in and collate the sCRPS estimates for each model
scrps_list_skew_t <- list()
for(i in 1:100){
  
  scrps_filename <- paste0("./2_empirical_analysis/model_results/growth_model/LA/emp_pooled_LA_N_25_model_skew_t_growth_cutoff_90_index_", i, "_scrps_ests.csv")
  scrps_list_skew_t[[i]] <- read_csv(scrps_filename, show_col_types = FALSE)$calc_scrps
  
}

full_scrps_skew_t <- do.call(c, scrps_list_skew_t)


scrps_list_skew_n <- list()
for(i in 1:100){
  
  scrps_filename <- paste0("./2_empirical_analysis/model_results/growth_model/LA/emp_pooled_LA_N_25_model_skew_normal_growth_cutoff_90_index_", i, "_scrps_ests.csv")
  scrps_list_skew_n[[i]] <- read_csv(scrps_filename, show_col_types = FALSE)$calc_scrps
  
}

full_scrps_skew_n <- do.call(c, scrps_list_skew_n)


scrps_list_normal <- list()
for(i in 1:100){
  
  scrps_filename <- paste0("./2_empirical_analysis/model_results/growth_model/LA/emp_pooled_LA_N_25_model_normal_growth_cutoff_90_index_", i, "_scrps_ests.csv")
  scrps_list_normal[[i]] <- read_csv(scrps_filename, show_col_types = FALSE)$calc_scrps
  
}

full_scrps_normal <- do.call(c, scrps_list_normal)


scrps_list_MLR <- list()
for(i in 1:100){
  
  scrps_filename <- paste0("./2_empirical_analysis/model_results/growth_model/LA/emp_pooled_LA_N_25_model_MLR_growth_cutoff_90_index_", i, "_scrps_ests.csv")
  scrps_list_MLR[[i]] <- read_csv(scrps_filename, show_col_types = FALSE)$calc_scrps
  
}

full_scrps_MLR <- do.call(c, scrps_list_MLR)


## Combine and save the sCRPS results
all_scrps <- cbind(full_scrps_skew_t, full_scrps_skew_n, full_scrps_normal, full_scrps_MLR)
scrps_comp <- apply(all_scrps, 2, function(x) c(mean(x), sd(x)/sqrt(length(x))))

all_scrps_filename <- paste0("./2_empirical_analysis/model_results/growth_model/LA/emp_pooled_LA_N_25_growth_cutoff_90_scrps_ests_combined.csv")
write_csv(as.data.frame(scrps_comp), all_scrps_filename)

## Read in the saved values
scrps_vals <- read_csv(all_scrps_filename, show_col_types = FALSE)

## Tidy the dataframe
scrps_vals <- scrps_vals |> 
  pivot_longer(everything(),
               names_to = "Model",
               values_to = "Estimates") |> 
  mutate(Type = rep(c("Mean", "SE"), each = 4)) |> 
  pivot_wider(names_from = Type,
              values_from = Estimates) |> 
  mutate(Model = str_replace_all(Model, "full_scrps_", ""))

## Generate Table 1
scrps_vals |> kbl(format = "latex", booktabs = TRUE, align = "lrr", linesep = "", escape = FALSE, digits = 6,
                  label = "scrps_table", caption = c("")) |> 
  kable_classic() |> 
  add_header_above(header = c(" " = 1, "sCRPS Value" = 2), align = "c",
                   bold = TRUE)


## Read in and collate the results for the POM and NDM models
model_types <- c("skew_t", "skew_normal", "normal", "MLR")
POM_results <- list()
NDM_results <- list()

for(i in model_types){
  
  NDM_scrps_filename <- paste0("./2_empirical_analysis/model_results/growth_model/NDM/nearest_distance_matching_model_", model_type, "_growth_cutoff_", mort_threshold, "_scrps_ests.csv")
  POM_scrps_filename <- paste0("./2_empirical_analysis/model_results/growth_model/POM/polygon_overlap_matching_model_", model_type, "_growth_cutoff_", mort_threshold, "_scrps_ests.csv")
  
  NDM_results[[i]] <- read_csv(NDM_scrps_filename, show_col_types = FALSE)
  POM_results[[i]] <- read_csv(POM_scrps_filename, show_col_types = FALSE)
  
}

scrps_vals2 <- rbind(do.call(rbind,
                             lapply(POM_results,
                                    function(x) c(mean(x), sd(x)/sqrt(length(x))))),
                     do.call(rbind,
                             lapply(NDM_results,
                                    function(x) c(mean(x), sd(x)/sqrt(length(x)))))) |> 
  as.data.frame()
col.names(scrps_vals2) <- c("Estimate", "SE")

## Tidy the dataframe
scrps_vals2 <- scrps_vals2 |> 
  mutate(Linkage = rep(c("POM", "NDM"), each = 4),
         Model = rep(c("Skewed t", "Skew Normal", "Normal", "MLR"), times = 2)) |> 
  mutate(Linkage = factor(Linkage, levels = c("POM", "NDM")),
         Model = factor(Model, levels = c("Skewed t", "Skew Normal", "Normal", "MLR"))) |> 
  relocate(c("Linkage", "Model"), .before = "Estimate" )

## Generate Table 2
scrps_vals2 |> kbl(format = "latex", booktabs = TRUE, align = "llrr", linesep = "", escape = FALSE, digits = 4,
                   label = "scrps_table", caption = c("")) |> 
  kable_classic() |> 
  collapse_rows(1) |> 
  add_header_above(header = c(" " = 2, "sCRPS Value" = 2), align = "c",
                   bold = TRUE)
