###########################################################################################
#### This script runs the growth model using the POM linkage on the empirical dataset. ####
###########################################################################################

# Arguments from command line ----
args <- commandArgs(trailingOnly=TRUE)
mort_threshold <- args[1]
model_type <- args[2]

## Load libraries and sampler functions ----
library(rstan) ## for growth model fit
library(readr) ## load and save results
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(terra)
library(dplyr)
library(data.table)

## Optional parallelization for rstan
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())

## Source necessary C++ helper functions
sourceCpp('./resources/code/cpp_code/two_stage_func.cpp')

## Set seed for reproducibility ----
set.seed(90210)

# Define the spatial domain
a_x <- 326096
a_y <- 4309939
b_x <- 328096
b_y <- 4311939
a_x2 <- a_x + 15
a_y2 <- a_y + 15
b_x2 <- b_x - 15
b_y2 <- b_y - 15


## Read in and process the specified dataset
file2015 <- read_csv("./resources/empirical_data/UER_lidar_canopy_segmentation/crown_attributes_2015.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(XTOP, YTOP, CANVOL2015, deltaCANVOL, propCANVOL) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2015, growth = deltaCANVOL, prop_growth = propCANVOL)

in_bounds <- which(file2015$x > a_x2 & file2015$x < b_x2 & file2015$y > a_y2 & file2015$y < b_y2)
scan_data <- file2015[in_bounds,]

CM <- read_csv("./2_empirical_analysis/competition_metrics_2015.csv", show_col_types = FALSE) %>%
  filter(XTOP > a_x2 & XTOP < b_x2 & YTOP > a_y2 & YTOP < b_y2) %>%
  select(LNV_norm, RSI_norm, ND_norm) %>%
  as.matrix()


## Read in the raster data for the covariates of interest
southness.rast <- scale(rast('./resources/empirical_data/Snodgrass_aspect_southness_1m.tif'))
wetness.rast <- scale(rast('./resources/empirical_data/Snodgrass_wetness_index_1m.tif'))
GDD.rast <- scale(rast('./resources/empirical_data/Snodgrass_Degree_Days_2013_2019.tif'))
SPP.rast <- scale(rast('./resources/empirical_data/Snodgrass_Snowpack_Persistence_DOY_2013_2019.tif'))
  
  
## Crop the rasters to D* and discard the originals
southness <- crop(southness.rast, ext(a_x, b_x, a_y, b_y))
wetness <- crop(wetness.rast, ext(a_x, b_x, a_y, b_y))
spp <- crop(SPP.rast, ext(a_x, b_x, a_y, b_y))
gdd <- crop(GDD.rast, ext(a_x, b_x, a_y, b_y))
rm(southness.rast, wetness.rast, SPP.rast, GDD.rast)
  
raster_list <- c(list(wetness), list(southness),
                 lapply(2:6, function(x) spp[[x]]), lapply(2:6, function(x) gdd[[x]]),
                 lapply(2:6, function(x) spp[[x]]*wetness), lapply(2:6, function(x) gdd[[x]]*wetness))

linked_data <- cbind(scan_data, CM)
linked_data$est_mort <- NA
if(mort_threshold == "80"){
  for(j in 1:nrow(linked_data)){
    if(linked_data$prop_growth[[j]] < -.2 | linked_data$prop_growth[[j]] > .6){
      linked_data$est_mort[[j]] <- 1
    } else{
      linked_data$est_mort[[j]] <- 0
    }
  }
}else if(mort_threshold == "90"){
  for(j in 1:nrow(linked_data)){
    if(linked_data$prop_growth[[j]] < -.1 | linked_data$prop_growth[[j]] > .6){
      linked_data$est_mort[[j]] <- 1
    } else{
      linked_data$est_mort[[j]] <- 0
    }
  }
}else if(mort_threshold == "100"){
  for(j in 1:nrow(linked_data)){
    if(linked_data$prop_growth[[j]] < 0 | linked_data$prop_growth[[j]] > .6){
      linked_data$est_mort[[j]] <- 1
    } else{
      linked_data$est_mort[[j]] <- 0
    }
  }
}


# Obtain the quantities necessary for the growth model
G <- linked_data$est_growth
M <- as.numeric(linked_data$est_mort)
S_1 <- as.numeric(linked_data$size.x)
S_2 <- as.numeric(linked_data$size.y)
s <- cbind(linked_data$x.x, linked_data$y.x)
X <- update_covars_arma(s, raster_list)
gc_index <- which(M == 0)
N <- length(gc_index)

X_new <- cbind(X[,2:3],
               apply(X[,4:8], 1, median),
               apply(X[,9:13], 1, median),
               apply(X[,14:18], 1, median),
               apply(X[,19:23], 1, median),
               linked_data$lnv_norm,
               linked_data$rsi_norm,
               linked_data$nd_norm)


# Run the growth model
stan_growth_data_2stage = list(N = N, # Number of Obs
                               K = ncol(X_new),
                               G = G[gc_index],
                               S = S_1[gc_index],
                               X = X_new[gc_index,],
                               mu_0 = rep(0, ncol(X_new)),
                               sig20 = diag(rep(2.5, ncol(X_new))),
                               c = .001,
                               d = .001,
                               a_alpha = 1,
                               b_alpha = 1,
                               nu_beta = 2,
                               nu_delta = 2,
                               mu_delta = 0,
                               sigma_delta = .1,
                               h_1 = .1,
                               gamma_max = max(file2015$size)) 


if(model_type == "skew_t"){
  
  stanfit_2stage <- stan(file = "./resources/code/STAN_code/STAN_growth_mod_skew_t.stan", # Stan file
                         data = stan_growth_data_2stage, # Data
                         warmup = 10000, # Number of iteration to burn-in
                         iter = 20000, # Total number of iterations
                         chains = 4, # Number of chains to run
                         thin = 10) # Number of chains to run
  
}else if(model_type == "skew_normal"){
  
  stanfit_2stage <- stan(file = "./resources/code/STAN_code/STAN_growth_mod_skew_normal.stan", # Stan file
                         data = stan_growth_data_2stage, # Data
                         warmup = 10000, # Number of iteration to burn-in
                         iter = 20000, # Total number of iterations
                         chains = 4, # Number of chains to run
                         thin = 10) # Number of chains to run
  
}else if(model_type == "normal"){
  
  stanfit_2stage <- stan(file = "./resources/code/STAN_code/STAN_growth_mod_alpha.stan", # Stan file
                         data = stan_growth_data_2stage, # Data
                         warmup = 10000, # Number of iteration to burn-in
                         iter = 20000, # Total number of iterations
                         chains = 4, # Number of chains to run
                         thin = 10) # Number of chains to run
  
}else if(model_type == "MLR"){
  
  X_MLR <- cbind(scale(S_1[gc_index]), X_new[gc_index,])
  
  stan_growth_data_2stage_MLR = list(N = N, # Number of Obs
                                     K = ncol(X_MLR),
                                     G = G[gc_index],
                                     X = X_MLR,
                                     mu_0 = rep(0, ncol(X_MLR)),
                                     sig20 = diag(rep(2.5, ncol(X_MLR))))
  
  stanfit_2stage <- stan(file = "./resources/code/STAN_code/STAN_growth_mod_MLR.stan", # Stan file
                         data = stan_growth_data_2stage_MLR, # Data
                         warmup = 10000, # Number of iteration to burn-in
                         iter = 20000, # Total number of iterations
                         chains = 4, # Number of chains to run
                         thin = 10)
  
}

growth_results <- rstan::extract(stanfit_2stage, permuted = TRUE)
growth_results_df <- as.data.frame(growth_results[!names(growth_results) %in% c("y_rep1", "y_rep2", "mu")])

y_rep1 <- as.data.frame(growth_results["y_rep1"])
y_rep2 <- as.data.frame(growth_results["y_rep2"])
calc_scrps <- scrps(as.matrix(y_rep1), as.matrix(y_rep2), G[gc_index])$pointwise

## Save the results
pom_growths_file <- paste0("./2_empirical_analysis/model_results/growth_model/POM/polygon_overlap_matching_growth_cutoff_", mort_threshold, "_growths.csv")
pom_results_file <- paste0("./2_empirical_analysis/model_results/growth_model/POM/polygon_overlap_matching_model_", model_type, "_growth_cutoff_", mort_threshold, ".csv")
rep1_results_file <- paste0("./2_empirical_analysis/model_results/growth_model/POM/polygon_overlap_matching_model_", model_type, "_growth_cutoff_", mort_threshold, "_rep1.csv")
rep2_results_file <- paste0("./2_empirical_analysis/model_results/growth_model/POM/polygon_overlap_matching_model_", model_type, "_growth_cutoff_", mort_threshold, "_rep2.csv")
scrps_results_file <- paste0("./2_empirical_analysis/model_results/growth_model/POM/polygon_overlap_matching_model_", model_type, "_growth_cutoff_", mort_threshold, "_scrps_ests.csv")

if(!file.exists(pom_growths_file)){
  write_csv(data.frame(growth = G[gc_index]), pom_growths_file)
}
write_csv(growth_results_df, file = pom_results_file)
write_csv(y_rep1, file = rep1_results_file)
write_csv(y_rep2, file = rep2_results_file)
write_csv(as.data.frame(calc_scrps), file = scrps_results_file)
