###############################################################################################################
#### This script runs the two stage linkage model sampler for the spatial record linkage model for varying ####
#### sizes of spatial domains and bounding box margins to obtain timing and linkage results from the       ####
#### different combinations to assess the performance gains and accuracy of the approximation.             ####
###############################################################################################################

# Arguments from command line ----
args <- commandArgs(trailingOnly=TRUE)
box_margin <- as.numeric(args[1])
area_size <- args[2]

## Load libraries and sampler functions ----
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(readr) 
library(spatstat)
library(terra)
library(mvtnorm)
library(data.table)

## Load C++ helper functions
sourceCpp('./resources/code/cpp_code/two_stage_func.cpp')

## Set seed for reproducibility ----
set.seed(90219)

## Specify the initial spatial domain
a_x <- 326996
a_y <- 4311239
b_x <- 327096
b_y <- 4311339


## Specify the area given the area_size argument
if(area_size == "100"){
  
  ## For area of size 100m^2
  a_x <- a_x - 0
  a_y <- a_y - 0
  b_x <- b_x + 0
  b_y <- b_y + 0
  
}else if(area_size == "150"){
  
  ## For area of size 150m^2
  a_x <- a_x - 25
  a_y <- a_y - 25
  b_x <- b_x + 25
  b_y <- b_y + 25
  
}else if(area_size == "200"){
  
  ## For area of size 200m^2
  a_x <- a_x - 50
  a_y <- a_y - 50
  b_x <- b_x + 50
  b_y <- b_y + 50
  
}else if(area_size == "300"){
  
  ## For area of size 300m^2
  a_x <- a_x - 100
  a_y <- a_y - 100
  b_x <- b_x + 100
  b_y <- b_y + 100
  
}



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

N <- floor(1.1*max(file_size))

## Update the spatial domain for the latents from D to D*
epsilon <- 1
a_x_exp <- a_x - epsilon
a_y_exp <- a_y - epsilon
b_x_exp <- b_x + epsilon
b_y_exp <- b_y + epsilon


Y_mat <- scan_data %>% as.matrix()
D_bounds <- c(a_x_exp, b_x_exp, a_y_exp, b_y_exp)

## Specify the initial values and the hyperparameters
s <- t(apply(Y_mat[1:N,1:2], 1, function(x) rmvnorm(1, x, diag(.1, 2))))
lambda <- c(0:(N-1), sample(0:(N-1), nrow(Y_mat) - N, replace = FALSE))
t <- c(0, 0)
c_sigma <- .0001
d_sigma <- .0001
sigma2 <- .1^2


hypers <- list("c_sigma" = c_sigma, "d_sigma" = d_sigma,
               "kappa" = 100, "nu" = 0, "sigma2_theta_min" = .0000005, "sigma2_theta_rate" = 1/2,
               "sigma2_t" = .0001^2)

inits <- list("lambda" = lambda,
              "s" = s,
              "sigma2" = sigma2,
              "theta" = 0,
              "sigma2_theta" = .001^2,
              "t" = t)


## Specify the results filenames
lambda_file <- paste0("./2_empirical_analysis/timing_results/empirical_linkage_lambda_results_box_", box_margin, "_area_", area_size, ".csv")
s_file <- paste0("./2_empirical_analysis/timing_results/empirical_linkage_s_results_box_", box_margin, "_area_", area_size, ".csv")
t_file <- paste0("./2_empirical_analysis/timing_results/empirical_linkage_t_results_box_", box_margin, "_area_", area_size, ".csv")
sigma2_file <- paste0("./2_empirical_analysis/timing_results/empirical_linkage_sigma2_results_box_", box_margin, "_area_", area_size, ".csv")

if(!file.exists(lambda_file)){

  linkage_post <- run_mcmc_SpRL_linkage_arma_ft_timing(n_iter = 5000, Y_mat = Y_mat, N = N, m = c(file_size),
                                                       D_bounds = D_bounds, dist = box_margin, init_vals = inits,
                                                       hyperparameters = hypers, sigma2_range = c(0, 3.175),
                                                       file_name_lambda = lambda_file, file_name_s = s_file,
                                                       file_name_t = t_file, file_name_theta = "",
                                                       file_name_sigma2 = sigma2_file, verbose = TRUE)
  
  ## Save the timing results
  timing <- as.data.frame(linkage_post$iter_timing)
  write_csv(timing, paste0("./2_empirical_analysis/timing_results/empirical_linkage_timing_results_box_", box_margin, "_area_", area_size, ".csv"))

}

