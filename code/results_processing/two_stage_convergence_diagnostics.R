#########################################################################
#### This script generates a point estimate for the linkage from the ####
#### two-stage model run on the empirical data                       ####
#########################################################################


## Load libraries and sampler functions ----
library(readr) ## load and save results
library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(dplyr)
library(rstan)

sourceCpp('./code/cpp_code/two_stage_func.cpp')

## Set seed for reproducibility ----
set.seed(90210)

burn_in <- 1:1000
linkage_list <- list()
n_latents_list <- list()

par(mfrow = c(2,2))
for(i in 1:4){
  linkage_file <- paste0("./code/empirical_data/empirical_linkage_lambda_results_sigma2_prior_weak_chain_", i, ".csv")
  linkage_list[[i]] <- fread(file = linkage_file, header = FALSE, skip = 2500) %>% as.matrix()
  # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE) %>% as.matrix()
  n_latents_list[[i]] <- apply(linkage_list[[i]], 1, function(x) length(unique(x)))
  plot(n_latents_list[[i]], type = 'l', main = paste0("Chain ", i))
}


par(mfrow = c(2,2))
for(i in 1:4){
  # linkage_file <- paste0("./code/empirical_data/empirical_linkage_lambda_results_sigma2_prior_strong_chain_", i, ".csv")
  # # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE, skip = 500) %>% as.matrix()
  # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE) %>% as.matrix()
  # n_latents_list[[i]] <- apply(linkage_list[[i]], 1, function(x) length(unique(x)))
  plot(n_latents_list[[i]][-burn_in], type = 'l', main = paste0("Chain ", i))
}

plot(n_latents_list[[1]][1000:2700], type = 'l', col = "blue")
plot(n_latents_list[[2]][1000:2700], type = 'l', col = "red")
plot(n_latents_list[[3]][1000:2700], type = 'l', col = "green", add = TRUE)
plot(n_latents_list[[4]][1000:2700], type = 'l', col = "yellow", add = TRUE)



res_dat <- data.frame(n_lats = c(n_latents_list[[1]][1001:2700],
                                 n_latents_list[[2]][1001:2700],
                                 n_latents_list[[3]][1001:2700],
                                 n_latents_list[[4]][1001:2700]),
                      index = rep(1:1700, 4),
                      dataset = rep(c(1,2,3,4), each = 1700))
ggplot(data= res_dat, aes(x = index, y = n_lats, group = dataset, color = as.factor(dataset))) +
  geom_line()


test_df <- data.frame(n_latents_list[[1]][1000:2700], n_latents_list[[2]][1000:2700],
                      n_latents_list[[3]][1000:2700], n_latents_list[[4]][1000:2700])
Rhat(as.matrix(test_df))


sigma2_list <- list()

par(mfrow = c(2,2))
for(i in 1:4){
  sigma2_file <- paste0("./code/empirical_data/empirical_linkage_sigma2_results_sigma2_prior_strong_chain_", i, ".csv")
  # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE, skip = 500) %>% as.matrix()
  sigma2_list[[i]] <- fread(file = sigma2_file, header = FALSE) %>% as.matrix()
  plot(sigma2_list[[i]], type = 'l', main = paste0("Chain ", i))
}


par(mfrow = c(2,2))
for(i in 1:4){
  # linkage_file <- paste0("./code/empirical_data/empirical_linkage_lambda_results_sigma2_prior_strong_chain_", i, ".csv")
  # # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE, skip = 500) %>% as.matrix()
  # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE) %>% as.matrix()
  # n_latents_list[[i]] <- apply(linkage_list[[i]], 1, function(x) length(unique(x)))
  plot(sigma2_list[[i]][-burn_in], type = 'l', main = paste0("Chain ", i))
}

test_df2 <- data.frame(sigma2_list[[1]][1000:2100], sigma2_list[[2]][1000:2100],
                      sigma2_list[[3]][1000:2100], sigma2_list[[4]][1000:2100])
Rhat(as.matrix(test_df2))








linkage_list <- list()
n_latents_list <- list()

par(mfrow = c(2,2))
for(i in 1:4){
  linkage_file <- paste0("./code/empirical_data/empirical_linkage_lambda_results_sigma2_prior_weak_chain_", i, ".csv")
  # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE, skip = 500) %>% as.matrix()
  linkage_list[[i]] <- fread(file = linkage_file, header = FALSE) %>% as.matrix()
  n_latents_list[[i]] <- apply(linkage_list[[i]], 1, function(x) length(unique(x)))
  plot(n_latents_list[[i]], type = 'l', main = paste0("Chain ", i))
}


par(mfrow = c(2,2))
for(i in 1:4){
  # linkage_file <- paste0("./code/empirical_data/empirical_linkage_lambda_results_sigma2_prior_strong_chain_", i, ".csv")
  # # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE, skip = 500) %>% as.matrix()
  # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE) %>% as.matrix()
  # n_latents_list[[i]] <- apply(linkage_list[[i]], 1, function(x) length(unique(x)))
  plot(n_latents_list[[i]][-burn_in], type = 'l', main = paste0("Chain ", i))
}

test_df <- data.frame(n_latents_list[[1]][1500:2100], n_latents_list[[2]][1500:2100],
                      n_latents_list[[3]][1500:2100], n_latents_list[[4]][1500:2100])
Rhat(as.matrix(test_df))


sigma2_list <- list()

par(mfrow = c(2,2))
for(i in 1:4){
  sigma2_file <- paste0("./code/empirical_data/empirical_linkage_sigma2_results_sigma2_prior_weak_chain_", i, ".csv")
  # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE, skip = 500) %>% as.matrix()
  sigma2_list[[i]] <- fread(file = sigma2_file, header = FALSE) %>% as.matrix()
  plot(sigma2_list[[i]], type = 'l', main = paste0("Chain ", i))
}


par(mfrow = c(2,2))
for(i in 1:4){
  # linkage_file <- paste0("./code/empirical_data/empirical_linkage_lambda_results_sigma2_prior_strong_chain_", i, ".csv")
  # # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE, skip = 500) %>% as.matrix()
  # linkage_list[[i]] <- fread(file = linkage_file, header = FALSE) %>% as.matrix()
  # n_latents_list[[i]] <- apply(linkage_list[[i]], 1, function(x) length(unique(x)))
  plot(sigma2_list[[i]][1000:2100][(1:110)*10], type = 'l', main = paste0("Chain ", i))
}

test_df2 <- data.frame(sigma2_list[[1]][1000:2100][(1:110)*10], sigma2_list[[2]][1000:2100][(1:110)*10],
                       sigma2_list[[3]][1000:2100][(1:110)*10], sigma2_list[[4]][1000:2100][(1:110)*10])
Rhat(as.matrix(test_df2))
