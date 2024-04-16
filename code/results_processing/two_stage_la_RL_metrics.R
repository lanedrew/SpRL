# pass from command line ----
# Rscript code/scripts/r/SpRL_sim_study.R
args <- commandArgs(trailingOnly=TRUE)

## Load libraries ----
library(readr) ## load and save results
library(codetools)
library(Rcpp)
library(RcppDist)
library(RcppArmadillo)
library(terra)
library(dplyr)
library(mvtnorm)

sourceCpp('./code/cpp_code/two_stage_func.cpp')
sourceCpp('./code/cpp_code/eval_links.cpp')


density_levels <- args[1]
noise_levels <- args[2]
alpha_vals <- args[3]
index <- 1:100

parameter_combos <- expand.grid(density_levels, noise_levels, alpha_vals, index)
names(parameter_combos) <- c("density", "noise", "alpha", "index")

pe_prec_rec_F1 <- list()

for(i in 1:nrow(parameter_combos)){
  print(paste0("Density: ", parameter_combos[i, "density"], ",
               Noise: ", parameter_combos[i, "noise"], ",
               Alpha: ", parameter_combos[i, "alpha"], ",
               Index: ", parameter_combos[i, "index"]))
  linkage_file <- paste0("./code/growth_sim_results/two_stage/linkage/", parameter_combos[i, "density"], "_density_",
                         parameter_combos[i, "noise"], "_noise_",
                         parameter_combos[i, "alpha"], "_alpha_",
                         parameter_combos[i, "index"], "_index_two_stage_linkage_results.csv")
  
  filename <- paste0("./code/growth_sim_data_F23/", parameter_combos[i, "density"], "_dens_",
                     parameter_combos[i, "noise"], "_noise_",
                     parameter_combos[i, "alpha"], "_alpha_sim_",
                     parameter_combos[i, "index"], ".csv")
  scan_data <- read_csv(filename, show_col_types = FALSE)
  scan_data <- scan_data %>% filter(!file == 0)
  
  
  linkage_sample <- read_csv(linkage_file, show_col_types = FALSE) %>% as.matrix()
  # lambda <- linkage_sample[-c(1:1000),] + 1
  lambda <- linkage_sample + 1
  
  prec_rec_F1_eval <- eval_links(z = lambda, true_id = scan_data$id)
  prec_rec_F1_dist <- cbind(apply(prec_rec_F1_eval, 2, mean), apply(prec_rec_F1_eval, 2, sd))
  colnames(prec_rec_F1_dist) <- c("precision_mean", "recall_mean", "F1_score_mean", "precision_sd", "recall_sd", "F1_score_sd")
  pe_prec_rec_F1[[i]] <- prec_rec_F1_dist
  
}

pe_prec_rec_F1_df <- do.call(rbind, pe_prec_rec_F1)
parameter_combos_all <- cbind(parameter_combos, pe_prec_rec_F1_df)

write_csv(as.data.frame(parameter_combos_all),
          paste0("./code/growth_sim_results/two_stage/linkage_processed/la_metrics/",
                 density_levels, "_dens_",
                 noise_levels, "_noise_",
                 alpha_vals, "_alpha_la_two_stage_linkage_metrics.csv"))
