# pass from command line ----
# Rscript code/scripts/r/SpRL_sim_study.R
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) stop("Pass in density (low or med or high), noise (small or medium or large),
                            alpha (1 or 2 or 3), index (integer from 1 to 100), and N threshold value.", call.=FALSE)
if (!(args[1] %in% c("low", "med", "high"))) stop("Pass in the density (low or med or high)", call.=FALSE)
if (!(args[2] %in% c("small", "medium", "large"))) stop("Pass in the noise level (small or medium or large)", call.=FALSE)

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

## Read in the raster data for the covariates of interest
slope.rast <- rast('./data/Snodgrass_slope_1m.tif')
southness.rast <- rast('./data/Snodgrass_aspect_southness_1m.tif')
wetness.rast <- rast('./data/Snodgrass_wetness_index_1m.tif')
DEM.rast <- rast('./data/Snodgrass_DEM_1m.tif')

## Set seed for reproducibility ----
set.seed(90210)

## Set the priors and relevant indexes
sigma2.prior <- c("uninformative")
density <- args[1]
noise <- args[2]
alpha <- args[3]
index <- args[4]
N_threshold <- as.numeric(args[5])

print(paste0("Density: ", density, ", Noise: ", noise, ", Alpha: ", alpha, ", Index: ", index))

## Read in the specified dataset
# filename <- paste0("./code/growth_sim_data_3/", density, "_dens_", noise, "_noise_sim_", index, ".csv")
filename <- paste0("./code/growth_sim_data_F23/", density, "_dens_", noise, "_noise_", alpha, "_alpha_sim_", index, ".csv")
scan_data <- read_csv(filename)
scan_data <- scan_data %>% filter(!file == 0)

## Specify file paths for saving the results
# linkage_file <- paste0("./code/growth_sim_results/two_stage/linkage/", density, "_density_", noise, "_noise_", index, "_index_two_stage_linkage_results.csv")
# params_file <- paste0("./code/growth_sim_results/two_stage/params/", density, "_density_", noise, "_noise_", index, "_index_two_stage_parameter_results.csv")
# raw_data_file <- paste0("./code/growth_sim_results/two_stage/raw_data/", density, "_density_", noise, "_noise_", index, "_index_two_stage_raw_data.RDS")
linkage_file <- paste0("./code/growth_sim_results/two_stage/linkage/", density, "_density_", noise, "_noise_", alpha, "_alpha_", index, "_index_", N_threshold, "_N_thresh_two_stage_linkage_results.csv")
params_file <- paste0("./code/growth_sim_results/two_stage/params/", density, "_density_", noise, "_noise_", alpha, "_alpha_", index, "_index_", N_threshold, "_N_thresh_two_stage_parameter_results.csv")
raw_data_file <- paste0("./code/growth_sim_results/two_stage/raw_data/", density, "_density_", noise, "_noise_", alpha, "_alpha_", index, "_index_", N_threshold, "_N_thresh_two_stage_raw_data.RDS")

## Obtain the file sizes and cumulative index m
file_size <- scan_data %>% 
  dplyr::select(file) %>%
  dplyr::group_by(file) %>%
  dplyr::summarise(n_i = n()) %>%
  dplyr::select(n_i)
file_size <- unlist(file_size)

N <- floor(((100 + N_threshold)/100)*max(file_size))

# N <- floor(1.1*max(file_size))
# if(N_threshold == 110){
#   N <- floor(1.1*max(file_size))
# }else if(N_threshold == 125){
#   N <- floor(1.25*max(file_size))
# }
# N <- floor(N_threshold*max(file_size))


## Specify the limits of the spatial domains given the specified density
if(density == "low"){
  
  a_x <- 326496
  a_y <- 4311439
  b_x <- 326596
  b_y <- 4311539
  
}else if(density == "med"){
  
  a_x <- 326996
  a_y <- 4311239
  b_x <- 327096
  b_y <- 4311339
  
}else if(density == "high"){
  
  a_x <- 327096
  a_y <- 4311239
  b_x <- 327196
  b_y <- 4311339
  
}

## Update the spatial domain for the latents from D to D*
epsilon <- 1
a_x_exp <- a_x - epsilon
a_y_exp <- a_y - epsilon
b_x_exp <- b_x + epsilon
b_y_exp <- b_y + epsilon

## Crop the rasters to D* and discard the originals
south <- scale(crop(southness.rast, ext(a_x - 15, b_x + 15, a_y - 15, b_y + 15)))
slope <- scale(crop(slope.rast, ext(a_x - 15, b_x + 15, a_y - 15, b_y + 15)))
wetness <- scale(crop(wetness.rast, ext(a_x - 15, b_x + 15, a_y - 15, b_y + 15)))
DEM <- scale(crop(DEM.rast, ext(a_x - 15, b_x + 15, a_y - 15, b_y + 15)))
rm(southness.rast, slope.rast, wetness.rast, DEM.rast)

Y_mat <- scan_data[,-c(4,6)] %>% as.matrix()
D_bounds <- c(a_x_exp, b_x_exp, a_y_exp, b_y_exp)

raster_list <- c(list(south), list(slope), list(wetness), list(DEM))

# lambda <- sample(0:(N-1), nrow(Y_mat), replace = TRUE)
# sample_index <- sample(1:nrow(Y_mat), N)
# s <- t(apply(Y_mat[sample_index,1:2], 1, function(x) rmvnorm(1, x, diag(.1, 2))))
s <- t(apply(Y_mat[1:N,1:2], 1, function(x) rmvnorm(1, x, diag(.1, 2))))
lambda <- c(0:(N-1), sample(0:(N-1), nrow(Y_mat) - N, replace = FALSE))
t <- c(0, 0)

hypers <- list("c_sigma" = 1.001, "d_sigma" = .0001,
               "kappa" = 100, "nu" = 0, "sigma2_theta_min" = .000000001, "sigma2_theta_rate" = 1/2,
               "sigma2_t" = .005^2)


inits <- list("lambda" = lambda,
              "s" = s,
              "sigma2" = .1,
              "theta" = 0,
              "sigma2_theta" = .0001^2,
              "t" = t)
res_file <- paste0("")


posterior_draws <- run_mcmc_SpRL_linkage_arma(n_iter = 5000, Y_mat = Y_mat, N = N, m = c(file_size),
                                              D_bounds = D_bounds, dist = 3.0, init_vals = inits,
                                              hyperparameters = hypers, sigma2_range = c(0, 3.175),
                                              file_name = "res_file", verbose = TRUE)


burn_in <- 1:2501
linkage_sample <- posterior_draws$lambda[-burn_in,] + 1
linkage_res <- do.call(cbind, (eval_links(z = linkage_sample, true_id = scan_data$id)))
results_df <- data.frame(sigma2 = posterior_draws$sigma2[-burn_in],
                         theta = posterior_draws$theta[-burn_in],
                         t1 = posterior_draws$t[-burn_in,1],
                         t2 = posterior_draws$t[-burn_in,2])

write_csv(as.data.frame(posterior_draws$lambda[-burn_in,]), linkage_file)
write_csv(results_df, params_file)
write_rds(posterior_draws$s[,,-burn_in], raw_data_file)
# write_csv(as.data.frame(linkage_res),
#           paste0("./code/growth_sim_results/two_stage/linkage_processed/",
#                  density, "_density_", noise, "_noise_",
#                  index, "_index_two_stage_linkage_metrics.csv"))
write_csv(as.data.frame(linkage_res),
          paste0("./code/growth_sim_results/two_stage/linkage_processed/",
                 density, "_density_", noise, "_noise_", alpha, "_alpha_", index,
                 "_index_", N_threshold, "_N_thresh_two_stage_linkage_metrics.csv"))
