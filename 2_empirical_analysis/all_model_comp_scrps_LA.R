####################################################################################
#### This script pools the sCRPS scores from the LA model runs into a single   #####
#### pooled value.                                                             #####
####################################################################################

## Load the necessary libraries
library(readr)
library(dplyr)
library(purrr)
library(loo)


scrps_list_skew_t <- list()

for(i in 1:100){
  
  scrps_filename <- paste0("./code/empirical_data/final_run/emp_pooled_LA_N_25_model_skew_t_covars_all_growth_cutoff_90_index_", i, "_scrps_ests.csv")
  scrps_list_skew_t[[i]] <- read_csv(scrps_filename, show_col_types = FALSE)$calc_scrps
  
}

full_scrps_skew_t <- do.call(c, scrps_list_skew_t)


scrps_list_skew_n <- list()

for(i in 1:100){
  
  scrps_filename <- paste0("./code/empirical_data/final_run/emp_pooled_LA_N_25_model_skew_normal_covars_all_growth_cutoff_90_index_", i, "_scrps_ests.csv")
  scrps_list_skew_n[[i]] <- read_csv(scrps_filename, show_col_types = FALSE)$calc_scrps
  
}

full_scrps_skew_n <- do.call(c, scrps_list_skew_n)


scrps_list_normal <- list()

for(i in 1:100){
  
  scrps_filename <- paste0("./code/empirical_data/final_run/emp_pooled_LA_N_25_model_normal_covars_all_growth_cutoff_90_index_", i, "_scrps_ests.csv")
  scrps_list_normal[[i]] <- read_csv(scrps_filename, show_col_types = FALSE)$calc_scrps
  
}

full_scrps_normal <- do.call(c, scrps_list_normal)


scrps_list_MLR <- list()

for(i in 1:100){
  
  scrps_filename <- paste0("./code/empirical_data/final_run/emp_pooled_LA_N_25_model_MLR_covars_all_growth_cutoff_90_index_", i, "_scrps_ests.csv")
  scrps_list_MLR[[i]] <- read_csv(scrps_filename, show_col_types = FALSE)$calc_scrps
  
}

full_scrps_MLR <- do.call(c, scrps_list_MLR)


all_scrps <- cbind(full_scrps_skew_t, full_scrps_skew_n, full_scrps_normal, full_scrps_MLR)
scrps_comp <- apply(all_scrps, 2, function(x) c(mean(x), sd(x)/sqrt(length(x))))

all_scrps_filename <- paste0("./code/empirical_data/final_run/emp_pooled_LA_N_25_covars_all_growth_cutoff_90_scrps_ests_combined.csv")
write_csv(as.data.frame(scrps_comp), all_scrps_filename)



