###############################################################################################################
#### This script defines and runs the two stage linkage model sampler for the spatial record linkage model ####
###############################################################################################################

# pass from command line ----
# Rscript code/two_stage_model/two_stage_emp_full_ext.R
args <- commandArgs(trailingOnly=TRUE)
chain <- args[1]

## Load libraries and sampler functions ----
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(readr) ## load and save results
library(spatstat)
library(terra)
library(mvtnorm)
library(data.table)

## Load Cpp functions
sourceCpp("./code/cpp_code/two_stage_func.cpp")

## Set seed for reproducibility ----
# seed_num <- 90210 + chain - 1
# set.seed(seed_num)
set.seed(90210)

## Specify the relevant indexes
a_x <- 326096
a_y <- 4309939
b_x <- 328096
b_y <- 4311939

## Read in and process the specified dataset
file2015 <- read_csv("./data/UER_lidar_canopy_segmentation/crown_attributes_2015.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(XTOP, YTOP, CANVOL2015) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2015) %>%
  mutate(file = 1)
file2019 <- read_csv("./data/UER_lidar_canopy_segmentation/crown_attributes_2019.csv", show_col_types = FALSE) %>%
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

N <- floor(1.1*max(file_size))


## Update the spatial domain for the latents from D to D*
epsilon <- 1
a_x_exp <- a_x - epsilon
a_y_exp <- a_y - epsilon
b_x_exp <- b_x + epsilon
b_y_exp <- b_y + epsilon


Y_mat <- scan_data %>% as.matrix()
D_bounds <- c(a_x_exp, b_x_exp, a_y_exp, b_y_exp)


lambda_file <- paste0("./code/empirical_data/empirical_linkage_lambda_results_chain_", chain,".csv")
s_file <- paste0("./code/empirical_data/empirical_linkage_s_results_chain_", chain,".csv")
t_file <- paste0("./code/empirical_data/empirical_linkage_t_results_chain_", chain,".csv")
sigma2_file <- paste0("./code/empirical_data/empirical_linkage_sigma2_results_chain_", chain,".csv")

lambda_chain <- fread(file = lambda_file, header = FALSE, skip = 2499) %>% as.matrix()
s_chain <- fread(file = s_file, header = FALSE, skip = 2499) %>% as.matrix()
t_chain <- fread(file = t_file, header = FALSE, skip = 2499) %>% as.matrix()
sigma2_chain <- fread(file = sigma2_file, header = FALSE, skip = 2499) %>% as.matrix()


lambda <- c(lambda_chain)
s <- vector2matrix(c(s_chain), dim = c(ncol(s_chain)/2, 2))
sigma2 <- c(sigma2_chain)
t <- c(t_chain)


hypers <- list("c_sigma" = .0001, "d_sigma" = .0001,
               "kappa" = 100, "nu" = 0, "sigma2_theta_min" = .0000005, "sigma2_theta_rate" = 1/2,
               "sigma2_t" = .0001^2)

inits <- list("lambda" = lambda,
              "s" = s,
              "sigma2" = sigma2,
              "theta" = 0,
              "sigma2_theta" = .001^2,
              "t" = t)

linkage_post <- run_mcmc_SpRL_linkage_arma_ft(n_iter = 5000, Y_mat = Y_mat, N = N, m = c(file_size),
                                              D_bounds = D_bounds, dist = 2.0, init_vals = inits,
                                              hyperparameters = hypers, sigma2_range = c(0, 3.175),
                                              file_name_lambda = lambda_file, file_name_s = s_file,
                                              file_name_t = t_file, file_name_theta = "",
                                              file_name_sigma2 = sigma2_file, verbose = TRUE)


res_rds <- paste0("./code/empirical_data/empirical_linkage_full_results_chain_", chain,"_part2.rds")
saveRDS(linkage_post, res_rds)