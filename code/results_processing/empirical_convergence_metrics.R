####################################################################################
#### This script generates a pooled sample from the record linkage model chains ####
####################################################################################

## Load libraries and sampler functions ----
library(readr) ## load and save results
library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(dplyr)

sourceCpp('./code/cpp_code/two_stage_func.cpp')

## Set seed for reproducibility ----
set.seed(90210)

linkage_list <- list()
latent_list <- list()
n_list <- list()

for(i in 1:4){
  linkage_file <- paste0("./code/empirical_data/empirical_linkage_lambda_results_N_25_sigma2_prior_weak_chain_", i, ".csv")
  linkage_sample <- fread(file = linkage_file, header = FALSE, skip = 2500) %>% as.matrix()
  # thin_n <- floor(nrow(linkage_sample)/5)
  # thin_seq <- (1:thin_n)*5
  # linkage_list[[i]] <- linkage_sample[thin_seq,]
  n_list[[i]] <- apply(linkage_sample, 1, function(x) length(unique(x)))
  
  # latent_file <- paste0("./code/empirical_data/empirical_linkage_s_results_sigma2_prior_weak_chain_", i, ".csv")
  # latent_sample <- fread(file = latent_file, header = FALSE, skip = 2500) %>% as.matrix()
  # latent_list[[i]] <- latent_sample[thin_seq,]
}
n_length <- lapply(n_list, length)

n_latents_df <- data.frame(N = do.call(c, n_list),
                           seq_val = c(seq_len(n_length[[1]]), seq_len(n_length[[2]]),
                                       seq_len(n_length[[3]]), seq_len(n_length[[4]])),
                           chain = rep(c(1, 2, 3, 4), times = do.call(c, n_length)))

n_latents_df |> ggplot(aes(x = seq_val, y = N, group = as.factor(chain), color = as.factor(chain))) +
  geom_line()

Rhat(as.matrix(cbind(n_list[[1]][1:3900], n_list[[2]][1:3900], n_list[[3]][1:3900], n_list[[4]][1:3900])))

# pooled_linkage <- do.call(rbind, linkage_list) %>% as.data.frame()
# pooled_latents <- do.call(rbind, latent_list) %>% as.data.frame()
# 
# write_csv(pooled_linkage, file = "./code/empirical_data/empirical_linkage_lambda_pooled.csv", col_names = FALSE)
# write_csv(pooled_latents, file = "./code/empirical_data/empirical_linkage_s_pooled.csv", col_names = FALSE)
