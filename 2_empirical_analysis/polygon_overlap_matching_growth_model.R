#########################################################################################################
#### This script runs the growth portion of the two-stage model for the spatial record linkage model ####
#########################################################################################################

# pass from command line ----
args <- commandArgs(trailingOnly=TRUE)
covars <- args[1]
# index <- as.numeric(args[2])
mort_threshold <- args[3]


## Load libraries and sampler functions ----
library(rstan) ## for growth model fit
library(readr) ## load and save results
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(terra)
library(dplyr)
library(data.table)


sourceCpp('./code/cpp_code/two_stage_func.cpp')
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())

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
file2015 <- read_csv("./data/UER_lidar_canopy_segmentation/crown_attributes_2015.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(XTOP, YTOP, CANVOL2015, deltaCANVOL, propCANVOL) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2015, growth = deltaCANVOL, prop_growth = propCANVOL)

in_bounds <- which(file2015$x > a_x2 & file2015$x < b_x2 & file2015$y > a_y2 & file2015$y < b_y2)
scan_data <- file2015[in_bounds,]

CM <- read_csv("./data/comp_metrics_2015_RSI.csv", show_col_types = FALSE) %>%
  filter(XTOP > a_x2 & XTOP < b_x2 & YTOP > a_y2 & YTOP < b_y2) %>%
  select(LNV_norm, RSI_norm, ND_norm) %>%
  as.matrix()

if(covars == "all"){
  
  ## Read in the raster data for the covariates of interest
  southness.rast <- scale(rast('./data/Snodgrass_aspect_southness_1m.tif'))
  wetness.rast <- scale(rast('./data/Snodgrass_wetness_index_1m.tif'))
  GDD.rast <- scale(rast('./data/Snodgrass_Degree_Days_2013_2019.tif'))
  SPP.rast <- scale(rast('./data/Snodgrass_Snowpack_Persistence_DOY_2013_2019.tif'))
  
  
  ## Crop the rasters to D* and discard the originals
  southness <- crop(southness.rast, ext(a_x, b_x, a_y, b_y))
  wetness <- crop(wetness.rast, ext(a_x, b_x, a_y, b_y))
  spp <- crop(SPP.rast, ext(a_x, b_x, a_y, b_y))
  gdd <- crop(GDD.rast, ext(a_x, b_x, a_y, b_y))
  rm(southness.rast, wetness.rast, SPP.rast, GDD.rast)
  
  raster_list <- c(list(wetness), list(southness),
                   lapply(2:6, function(x) spp[[x]]), lapply(2:6, function(x) gdd[[x]]),
                   lapply(2:6, function(x) spp[[x]]*wetness), lapply(2:6, function(x) gdd[[x]]*wetness))
  
  
}else if(covars == "subset"){
  
  ## Read in the raster data for the covariates of interest
  southness.rast <- scale(rast('./data/Snodgrass_aspect_southness_1m.tif'))
  wetness.rast <- scale(rast('./data/Snodgrass_wetness_index_1m.tif'))
  GDD.rast <- scale(rast('./data/Snodgrass_Degree_Days_2013_2019.tif'))
  SPP.rast <- scale(rast('./data/Snodgrass_Snowpack_Persistence_DOY_2013_2019.tif'))
  
  
  ## Crop the rasters to D* and discard the originals
  southness <- crop(southness.rast, ext(a_x, b_x, a_y, b_y))
  wetness <- crop(wetness.rast, ext(a_x, b_x, a_y, b_y))
  spp <- crop(SPP.rast, ext(a_x, b_x, a_y, b_y))
  gdd <- crop(GDD.rast, ext(a_x, b_x, a_y, b_y))
  rm(southness.rast, wetness.rast, SPP.rast, GDD.rast)
  
  raster_list <- c(list(wetness), list(southness),
                   lapply(2:6, function(x) spp[[x]]), lapply(2:6, function(x) gdd[[x]]),
                   lapply(2:6, function(x) spp[[x]]*wetness), lapply(2:6, function(x) gdd[[x]]*wetness))
  
}



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
G <- linked_data$growth
M <- as.numeric(linked_data$est_mort)
S_1 <- as.numeric(linked_data$size)
s <- cbind(linked_data$x, linked_data$y)
X <- update_covars_arma(s, raster_list)
if(covars == "all"){
  
  X_new <- cbind(X[,2:3],
                 apply(X[,4:8], 1, median),
                 apply(X[,9:13], 1, median),
                 apply(X[,14:18], 1, median),
                 apply(X[,19:23], 1, median),
                 linked_data$LNV_norm,
                 linked_data$RSI_norm,
                 linked_data$ND_norm)
  
}else if(covars == "subset"){
  
  X_new <- cbind(X[,2:3],
                 apply(X[,4:8], 1, median),
                 apply(X[,9:13], 1, median),
                 apply(X[,14:18], 1, median),
                 apply(X[,19:23], 1, median))
  
}

# Run the growth model
stan_growth_data_2stage = list(N = length(G[which(M == 0)]), # Number of Obs
                               G = G[which(M == 0)],
                               mu_0 = rep(0, ncol(X_new)),
                               S = S_1[which(M == 0)],
                               sig20 = diag(rep(2.5, ncol(X_new))),
                               K = ncol(X_new),
                               X = X_new[which(M == 0),],
                               c = .0001,
                               d = .0001,
                               max_tau2 = 75,
                               a_alpha = 1,
                               b_alpha = 1)  


stanfit_2stage <- stan(file = "./code/STAN_models/STAN_growth_mod_alpha.stan", # Stan file
                       data = stan_growth_data_2stage, # Data
                       warmup = 10000, # Number of iteration to burn-in
                       iter = 20000, # Total number of iterations
                       chains = 4, # Number of chains to run
                       thin = 10)

growth_results <- rstan::extract(stanfit_2stage, permuted = TRUE)
growth_results <- as.data.frame(growth_results)


## Save the results
growth_results_file <- paste0("./code/naive_growth_models/polygon_overlap_matching_2015_covars_", covars, "_growth_cutoff_", mort_threshold, ".csv")

write_csv(growth_results, file = growth_results_file)
