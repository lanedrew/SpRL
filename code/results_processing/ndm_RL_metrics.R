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

ndm_prec_rec_F1 <- list()

for(i in 1:nrow(parameter_combos)){
  print(paste0("Density: ", parameter_combos[i, "density"], ",
               Noise: ", parameter_combos[i, "noise"], ",
               Alpha: ", parameter_combos[i, "alpha"], ",
               Index: ", parameter_combos[i, "index"]))
  linkage_file <- paste0("./code/growth_sim_results/NDM/linkage_processed/", parameter_combos[i, "density"], "_density_",
                         parameter_combos[i, "noise"], "_noise_",
                         parameter_combos[i, "alpha"], "_alpha_",
                         parameter_combos[i, "index"], "_index_NDM_linkage_metrics.csv")
  
  
  linkage_metrics <- read_csv(linkage_file, show_col_types = FALSE) %>% as.matrix()
  ndm_prec_rec_F1[[i]] <- linkage_metrics
  
}

ndm_prec_rec_F1_df <- do.call(rbind, ndm_prec_rec_F1)
parameter_combos_all <- cbind(parameter_combos, ndm_prec_rec_F1_df)

write_csv(as.data.frame(parameter_combos_all),
          paste0("./code/growth_sim_results/NDM/linkage_processed/linkage_metrics/",
                 density_levels, "_dens_",
                 noise_levels, "_noise_",
                 alpha_vals, "_alpha_NDM_linkage_metrics.csv"))
