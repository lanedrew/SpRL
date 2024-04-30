############################################################################################################
#### This script generates a pooled sample from the record linkage model chains for the empirical data. ####
############################################################################################################

## Load libraries and sampler functions ----
library(readr) ## load and save results
library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(dplyr)

## Source necessary C++ helper functions
sourceCpp('./resources/code/cpp_code/two_stage_func.cpp')

## Set seed for reproducibility ----
set.seed(90210)

linkage_list <- list()
latent_list <- list()

## Load and collate the chains
for(i in 1:4){
  linkage_file <- paste0("./2_empirical_analysis/model_results/record_linkage_model/empirical_linkage_lambda_results_N_25_sigma2_prior_weak_chain_", i, ".csv")
  linkage_sample <- fread(file = linkage_file, header = FALSE, skip = 2500) %>% as.matrix()
  thin_n <- floor(nrow(linkage_sample)/5)
  thin_seq <- (1:thin_n)*5
  linkage_list[[i]] <- linkage_sample[thin_seq,]
  
  latent_file <- paste0("./2_empirical_analysis/model_results/record_linkage_model/empirical_linkage_s_results_N_25_sigma2_prior_weak_chain_", i, ".csv")
  latent_sample <- fread(file = latent_file, header = FALSE, skip = 2500) %>% as.matrix()
  latent_list[[i]] <- latent_sample[thin_seq,]
}

pooled_linkage <- do.call(rbind, linkage_list) %>% as.data.frame()
pooled_latents <- do.call(rbind, latent_list) %>% as.data.frame()

## Save the pooled chains for lambda and s
write_csv(pooled_linkage, file = "./2_empirical_analysis/model_results/record_linkage_model/empirical_linkage_lambda_pooled_N_25.csv", col_names = FALSE)
write_csv(pooled_latents, file = "./2_empirical_analysis/model_results/record_linkage_model/empirical_linkage_s_pooled_N_25.csv", col_names = FALSE)
