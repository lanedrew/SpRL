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

in_bounds <- which(file2015$x > a_x2 & file2015$x < b_x2 & file2015$y > a_y2 & file2015$y < b_y2)
scan_data <- rbind(file2015, file2019)

in_bounds_full <- which((scan_data$x > a_x2 & scan_data$x < b_x2 & scan_data$y > a_y2 & scan_data$y < b_y2 & scan_data$file == 1) | (scan_data$file == 2))

# CM <- read_csv("./data/comp_metrics_2015.csv", show_col_types = FALSE) %>%
#   filter(XTOP > a_x2 & XTOP < b_x2 & YTOP > a_y2 & YTOP < b_y2) %>%
#   select(LNV_norm, AND_norm, ND_norm) %>%
#   as.matrix()
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


linkage_file <- paste0("./code/naive_growth_models/near_distance_matching_lambda.csv")
linkage_sample <- read_csv(file = linkage_file, col_names = TRUE) %>% as.matrix()
linkage_sample <- linkage_sample + 1


## Run the growth model for the thinned lambda configurations with appropriate latents
data_per_it <- scan_data[in_bounds_full,]

data_per_it$id <- linkage_sample[in_bounds_full]
lnv_norm <- unsplit(lapply(split(CM[,1], data_per_it$id[1:length(in_bounds)]), mean), data_per_it$id[1:length(in_bounds)])
rsi_norm <- unsplit(lapply(split(CM[,2], data_per_it$id[1:length(in_bounds)]), mean), data_per_it$id[1:length(in_bounds)])
nd_norm <- unsplit(lapply(split(CM[,3], data_per_it$id[1:length(in_bounds)]), mean), data_per_it$id[1:length(in_bounds)])
data_per_it$x <- unsplit(lapply(split(data_per_it$x, data_per_it$id), mean), data_per_it$id)
data_per_it$y <- unsplit(lapply(split(data_per_it$y, data_per_it$id), mean), data_per_it$id)

data_file_1 <- data_per_it %>%
  filter(file == 1) %>%
  mutate(lnv_norm = lnv_norm,
         rsi_norm = rsi_norm,
         nd_norm = nd_norm)
data_file_2 <- data_per_it %>%
  filter(file == 2)

if(length(data_file_1$id) != length(unique(data_file_1$id))){
  
  N_1 <- unique(data_file_1$id)
  links_within <- list()
  
  for(j in 1:length(N_1)){
    
    links_within[[j]] <- which(data_file_1$id == N_1[j])
    
  }
  
  if(any(lapply(links_within, length) > 1)){
    
    linked_records <- which(lapply(links_within, length) > 1)
    
    for(k in 1:length(linked_records)){
      
      data_file_1$size[links_within[[linked_records[k]]]] <- sum(data_file_1$size[links_within[[linked_records[k]]]])
      # data_file_1$lnv_norm[links_within[[linked_records[k]]]] <- mean(data_file_1$lnv_norm[links_within[[linked_records[k]]]])
      # data_file_1$and_norm[links_within[[linked_records[k]]]] <- mean(data_file_1$and_norm[links_within[[linked_records[k]]]])
      # data_file_1$nd_norm[links_within[[linked_records[k]]]] <- mean(data_file_1$nd_norm[links_within[[linked_records[k]]]])
      
    }
  }
  
  data_file_1 <- data_file_1[!duplicated(data_file_1$id),]
  
}

if(length(data_file_2$id) != length(unique(data_file_2$id))){
  
  N_2 <- unique(data_file_2$id)
  links_within <- list()
  
  for(j in 1:length(N_2)){
    
    links_within[[j]] <- which(data_file_2$id == N_2[j])
    
  }
  
  if(any(lapply(links_within, length) > 1)){
    
    linked_records <- which(lapply(links_within, length) > 1)
    
    for(k in 1:length(linked_records)){
      
      data_file_2$size[links_within[[linked_records[k]]]] <- sum(data_file_2$size[links_within[[linked_records[k]]]])
      
    }
  }
  
  data_file_2 <- data_file_2[!duplicated(data_file_2$id),]
  
}

## Merge the two files after merging records within the same files
linked_data <- inner_join(data_file_1, data_file_2, by = "id")
linked_data <- linked_data %>% mutate(est_growth = (size.y - size.x)/4,
                                      delta_can = (size.y - size.x)/size.x,
                                      log_est_growth = log(est_growth))


linked_data$est_mort <- NA
if(mort_threshold == "80"){
  for(j in 1:nrow(linked_data)){
    if(linked_data$delta_can[[j]] < -.2 | linked_data$delta_can[[j]] > .6){
      linked_data$est_mort[[j]] <- 1
    } else{
      linked_data$est_mort[[j]] <- 0
    }
  }
}else if(mort_threshold == "90"){
  for(j in 1:nrow(linked_data)){
    if(linked_data$delta_can[[j]] < -.1 | linked_data$delta_can[[j]] > .6){
      linked_data$est_mort[[j]] <- 1
    } else{
      linked_data$est_mort[[j]] <- 0
    }
  }
}else if(mort_threshold == "100"){
  for(j in 1:nrow(linked_data)){
    if(linked_data$delta_can[[j]] < 0 | linked_data$delta_can[[j]] > .6){
      linked_data$est_mort[[j]] <- 1
    } else{
      linked_data$est_mort[[j]] <- 0
    }
  }
}


# Obtain the quantities necessary for the growth model
G <- linked_data$est_growth
# G <- linked_data$log_est_growth
M <- as.numeric(linked_data$est_mort)
S_1 <- as.numeric(linked_data$size.x)
S_2 <- as.numeric(linked_data$size.y)
s <- cbind(linked_data$x.x, linked_data$y.x)
X <- update_covars_arma(s, raster_list)
if(covars == "all"){
  
  X_new <- cbind(X[,2:3],
                 apply(X[,4:8], 1, median),
                 apply(X[,9:13], 1, median),
                 apply(X[,14:18], 1, median),
                 apply(X[,19:23], 1, median),
                 linked_data$lnv_norm,
                 linked_data$rsi_norm,
                 linked_data$nd_norm)
  
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


ndm_results_file <- paste0("./code/naive_growth_models/near_distance_matching_covars_", covars, "_growth_cutoff_", mort_threshold, ".csv")

write_csv(growth_results, file = ndm_results_file)
