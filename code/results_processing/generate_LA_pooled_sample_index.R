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


# linkage_file <- paste0("./code/empirical_data/empirical_linkage_lambda_pooled.csv")
#   # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE, skip = 500) %>% as.matrix()
# linkage_sample <- fread(file = linkage_file, header = TRUE) %>% as.matrix()
# n_latents <- apply(linkage_sample, 1, function(x) length(unique(x)))
# plot(n_latents, type = 'l', main = "Pooled Chain")
# 
# sample_index <- sample(1:nrow(linkage_sample), 100, replace = FALSE)
# write_csv(as.data.frame(sample_index), "./code/empirical_data/LA_sample_index_pooled.csv", col_names = FALSE)

linkage_file <- paste0("./code/empirical_data/empirical_linkage_lambda_pooled_N_25.csv")
# linkage_list[[i]] <- fread(file = linkage_file, header = FALSE, skip = 500) %>% as.matrix()
linkage_sample <- fread(file = linkage_file, header = TRUE) %>% as.matrix()
n_latents <- apply(linkage_sample, 1, function(x) length(unique(x)))
plot(n_latents, type = 'l', main = "Pooled Chain")

sample_index <- sample(1:nrow(linkage_sample), 100, replace = FALSE)
write_csv(as.data.frame(sample_index), "./code/empirical_data/LA_sample_index_pooled_N_25.csv", col_names = FALSE)


# latent_file <- paste0("./code/empirical_data/empirical_linkage_s_pooled.csv")
# # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE, skip = 500) %>% as.matrix()
# latent_sample <- fread(file = latent_file, header = TRUE) %>% as.matrix()
