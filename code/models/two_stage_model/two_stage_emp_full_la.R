#########################################################################################################
#### This script runs the growth portion of the two-stage model for the spatial record linkage model ####
#########################################################################################################

# pass from command line ----
args <- commandArgs(trailingOnly=TRUE)
covars <- args[1]
index <- as.numeric(args[2])
mort_threshold <- args[3]
model_type <- args[4]


## Load libraries and sampler functions ----
library(rstan) ## for growth model fit
library(readr) ## load and save results
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(terra)
library(dplyr)
library(data.table)
library(loo)


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


# linkage_file <- paste0("./code/empirical_data/empirical_linkage_lambda_pooled.csv")
# latent_file <- paste0("./code/empirical_data/empirical_linkage_s_pooled.csv")
# sample_index <- fread(file = "./code/empirical_data/LA_sample_index_pooled.csv", header = FALSE) %>% as.matrix()
linkage_file <- paste0("./code/empirical_data/empirical_linkage_lambda_pooled_N_25.csv")
latent_file <- paste0("./code/empirical_data/empirical_linkage_s_pooled_N_25.csv")
sample_index <- fread(file = "./code/empirical_data/LA_sample_index_pooled_N_25.csv", header = FALSE) %>% as.matrix()
current_index <- sample_index[index]


linkage_sample <- fread(file = linkage_file, skip = current_index, header = FALSE, nrows = 1) %>% as.matrix()
latent_sample <- fread(file = latent_file, skip = current_index, header = FALSE, nrows = 1) %>% as.matrix()

linkage_sample <- linkage_sample + 1
s_configs <- vector2matrix(latent_sample, dim = c(ncol(latent_sample)/2, 2))
rm(latent_sample)

## Obtain the thinned linkage and corresponding latents for the two-stage model

N <- dim(s_configs)[2]


## Run the growth model for the thinned lambda configurations with appropriate latents
data_per_it <- scan_data[in_bounds_full,]

data_per_it$id <- linkage_sample[in_bounds_full]
lnv_norm <- unsplit(lapply(split(CM[,1], data_per_it$id[1:length(in_bounds)]), mean), data_per_it$id[1:length(in_bounds)])
rsi_norm <- unsplit(lapply(split(CM[,2], data_per_it$id[1:length(in_bounds)]), mean), data_per_it$id[1:length(in_bounds)])
nd_norm <- unsplit(lapply(split(CM[,3], data_per_it$id[1:length(in_bounds)]), mean), data_per_it$id[1:length(in_bounds)])
data_per_it$x <- s_configs[c(data_per_it$id), 1]
data_per_it$y <- s_configs[c(data_per_it$id), 2]


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
                                      delta_can = (size.y - size.x)/size.x)


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
M <- as.numeric(linked_data$est_mort)
S_1 <- as.numeric(linked_data$size.x)
s <- cbind(linked_data$x.x, linked_data$y.x)
X <- update_covars_arma(s, raster_list)
gc_index <- which(M == 0)
N <- length(gc_index)

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
  
  stanfit_2stage <- stan(file = "./code/STAN_models/STAN_growth_mod_skew_t_test2.stan", # Stan file
                         data = stan_growth_data_2stage, # Data
                         iter = 5000,
                         chains = 4,
                         thin = 5) # Number of chains to run
  
}else if(model_type == "skew_normal"){
  
  stanfit_2stage <- stan(file = "./code/STAN_models/STAN_growth_mod_skew_normal2.stan", # Stan file
                         data = stan_growth_data_2stage, # Data
                         iter = 5000,
                         chains = 4,
                         thin = 5) # Number of chains to run
  
}else if(model_type == "normal"){
  
  stanfit_2stage <- stan(file = "./code/STAN_models/STAN_growth_mod_alpha.stan", # Stan file
                         data = stan_growth_data_2stage, # Data
                         iter = 5000,
                         chains = 4,
                         thin = 5) # Number of chains to run
  
}else if(model_type == "MLR"){
  
  X_MLR <- cbind(scale(S_1[gc_index]), X_new[gc_index,])
  
  stan_growth_data_2stage_MLR = list(N = N, # Number of Obs
                                     K = ncol(X_MLR),
                                     G = G[gc_index],
                                     X = X_MLR,
                                     mu_0 = rep(0, ncol(X_MLR)),
                                     sig20 = diag(rep(2.5, ncol(X_MLR))))
  
  stanfit_2stage <- stan(file = "./code/STAN_models/STAN_growth_mod_MLR.stan", # Stan file
                         data = stan_growth_data_2stage_MLR, # Data
                         iter = 5000,
                         chains = 4,
                         thin = 5)
  
}


growth_results <- rstan::extract(stanfit_2stage, permuted = TRUE)
growth_results_df <- as.data.frame(growth_results[!names(growth_results) %in% c("y_rep1", "y_rep2", "mu")])

y_rep1 <- as.data.frame(growth_results["y_rep1"])
y_rep2 <- as.data.frame(growth_results["y_rep2"])

calc_scrps <- scrps(as.matrix(y_rep1), as.matrix(y_rep2), G[gc_index])$pointwise


## Save the results
# la_results_file <- paste0("./code/empirical_data/emp_LA_covars_", covars, "_growth_cutoff_", mort_threshold, "_index_", index, ".csv")
la_results_file <- paste0("./code/empirical_data/final_run/emp_pooled_LA_N_25_model_", model_type, "_covars_", covars, "_growth_cutoff_", mort_threshold, "_index_", index, ".csv")
rep1_results_file <- paste0("./code/empirical_data/final_run/emp_pooled_LA_N_25_model_", model_type, "_covars_", covars, "_growth_cutoff_", mort_threshold, "_index_", index, "_rep1.csv")
rep2_results_file <- paste0("./code/empirical_data/final_run/emp_pooled_LA_N_25_model_", model_type, "_covars_", covars, "_growth_cutoff_", mort_threshold, "_index_", index, "_rep2.csv")
scrps_results_file <- paste0("./code/empirical_data/final_run/emp_pooled_LA_N_25_model_", model_type, "_covars_", covars, "_growth_cutoff_", mort_threshold, "_index_", index, "_scrps_ests.csv")


write_csv(growth_results_df, file = la_results_file)
write_csv(y_rep1, file = rep1_results_file)
write_csv(y_rep2, file = rep2_results_file)
write_csv(as.data.frame(calc_scrps), file = scrps_results_file)
