#### This script runs the extension of the two stage linkage model sampler for the empirical data. ####
#######################################################################################################

# Arguments from command line ----
args <- commandArgs(trailingOnly=TRUE)
chain <- args[1]
sigma2_prior <- args[2]
N_threshold <- as.numeric(args[3])

## Load libraries and sampler functions ----
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(readr) ## load and save results
library(spatstat)
library(terra)
library(mvtnorm)
library(data.table)

## Load C++ helper functions
sourceCpp('./resources/code/cpp_code/two_stage_func.cpp')

## Set seed for reproducibility ----
seed_num <- 90210 + chain - 1
set.seed(seed_num)

## Specify the relevant indexes
a_x <- 326096
a_y <- 4309939
b_x <- 328096
b_y <- 4311939

## Read in and process the specified dataset
file2015 <- read_csv("./resources/empirical_data/UER_lidar_canopy_segmentation/crown_attributes_2015.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(XTOP, YTOP, CANVOL2015) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2015) %>%
  mutate(file = 1)
file2019 <- read_csv("./resources/empirical_data/UER_lidar_canopy_segmentation/crown_attributes_2019.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(XTOP, YTOP, CANVOL2019) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2019) %>%
  mutate(file = 2)


scan_data <- rbind(file2015, file2019)

if(nrow(scan_data) == 0) stop("Subset is empty.", call.=FALSE)


## Obtain the file sizes and cumulative index m
file_size <- scan_data %>% 
  dplyr::select(file) %>%
  dplyr::group_by(file) %>%
  dplyr::summarise(n_i = n()) %>%
  dplyr::select(n_i)
file_size <- unlist(file_size)

N <- floor(((100 + N_threshold)/100)*max(file_size))

## Update the spatial domain for the latents from D to D*
epsilon <- 1
a_x_exp <- a_x - epsilon
a_y_exp <- a_y - epsilon
b_x_exp <- b_x + epsilon
b_y_exp <- b_y + epsilon


Y_mat <- scan_data %>% as.matrix()
D_bounds <- c(a_x_exp, b_x_exp, a_y_exp, b_y_exp)

## Specify results file names
lambda_file <- paste0("./2_empirical_analysis/model_results/record_linkage_model/empirical_linkage_lambda_results_N_",
                      N_threshold, "_sigma2_prior_", sigma2_prior, "_chain_", chain,".csv")
s_file <- paste0("./2_empirical_analysis/model_results/record_linkage_model/empirical_linkage_s_results_N_",
                 N_threshold, "_sigma2_prior_", sigma2_prior, "_chain_", chain,".csv")
t_file <- paste0("./2_empirical_analysis/model_results/record_linkage_model/empirical_linkage_t_results_N_",
                 N_threshold, "_sigma2_prior_", sigma2_prior, "_chain_", chain,".csv")
sigma2_file <- paste0("./2_empirical_analysis/model_results/record_linkage_model/empirical_linkage_sigma2_results_N_",
                      N_threshold, "_sigma2_prior_", sigma2_prior, "_chain_", chain,".csv")

## Load the value of the previous chain to continue tht sampler
lambda_chain <- fread(file = lambda_file, header = FALSE, skip = 2499) %>% as.matrix()
s_chain <- fread(file = s_file, header = FALSE, skip = 2499) %>% as.matrix()
t_chain <- fread(file = t_file, header = FALSE, skip = 2499) %>% as.matrix()
sigma2_chain <- fread(file = sigma2_file, header = FALSE, skip = 2499) %>% as.matrix()

## Specify initial values and hyperparameters
lambda <- c(lambda_chain[nrow(lambda_chain),])
s <- vector2matrix(c(s_chain[nrow(s_chain),]), dim = c(ncol(s_chain)/2, 2))
sigma2 <- c(sigma2_chain)[nrow(sigma2_chain)]
t <- c(t_chain[nrow(t_chain),])

if(sigma2_prior == "weak"){
  c_sigma = .0001
  d_sigma = .0001
} else if(sigma2_prior == "strong"){
  c_sigma = 100
  d_sigma = 10
}

hypers <- list("c_sigma" = c_sigma, "d_sigma" = d_sigma,
               "kappa" = 100, "nu" = 0, "sigma2_theta_min" = .0000005, "sigma2_theta_rate" = 1/2,
               "sigma2_t" = .0001^2)

inits <- list("lambda" = lambda,
              "s" = s,
              "sigma2" = sigma2,
              "theta" = 0,
              "sigma2_theta" = .001^2,
              "t" = t)

## Run the model
linkage_post <- run_mcmc_SpRL_linkage_arma_ft(n_iter = 5000, Y_mat = Y_mat, N = N, m = c(file_size),
                                              D_bounds = D_bounds, dist = 3.0, init_vals = inits,
                                              hyperparameters = hypers, sigma2_range = c(0, 3.175),
                                              file_name_lambda = lambda_file, file_name_s = s_file,
                                              file_name_t = t_file, file_name_theta = "",
                                              file_name_sigma2 = sigma2_file, verbose = TRUE)

## Save the model output for the second stage
res_rds <- paste0("./2_empirical_analysis/model_results/record_linkage_model/empirical_linkage_full_results_N_",
                  N_threshold, "_sigma2_prior_", sigma2_prior, "_chain_", chain,"_part2.rds")
saveRDS(linkage_post, res_rds)