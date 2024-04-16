#########################################################################################################
#### This script runs the growth portion of the two-stage model for the spatial record linkage model ####
#########################################################################################################

# pass from command line ----
# Rscript code/joint_growth_sim_study.R
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) stop("Pass in density (low or med or high), 
                            noise level (small or medium or large),
                            alpha level (1 or 2 or 3),
                            the index (1 to 100),
                            the N threshold,
                            and the linkage set", call.=FALSE)
if (!(args[1] %in% c("low", "med", "high"))) stop("Pass in the density (low or medium or high)", call.=FALSE)
if (!(args[2] %in% c("small", "medium", "large"))) stop("Pass in the noise level (small or medium or large)", call.=FALSE)


## Load libraries and sampler functions ----
library(rstan) ## for growth model fit
library(readr) ## load and save results
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(terra)
library(GreedyEPL) ## for linkage point estimate
library(dplyr)


sourceCpp('./code/cpp_code/two_stage_func.cpp')
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())


## Read in the raster data for the covariates of interest
slope.rast <- rast('./data/Snodgrass_slope_1m.tif')
southness.rast <- rast('./data/Snodgrass_aspect_southness_1m.tif')
wetness.rast <- rast('./data/Snodgrass_wetness_index_1m.tif')
DEM.rast <- rast('./data/Snodgrass_DEM_1m.tif')

## Set seed for reproducibility ----
set.seed(90210)

## Set the priors and relevant indexes
density <- args[1]
noise <- args[2]
alpha <- args[3]
index <- args[4]
N_threshold <- args[5]
linkage_set <- args[6]

print(paste0("LA Density: ", density, ", Noise: ", noise, ", Alpha: ", alpha, ", Index: ", index, ", N Threshold: ", N_threshold))

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


## Crop the rasters to D* and discard the originals
southness <- scale(crop(southness.rast, ext(a_x - 15, b_x + 15, a_y - 15, b_y + 15)))
slope <- scale(crop(slope.rast, ext(a_x - 15, b_x + 15, a_y - 15, b_y + 15)))
wetness <- scale(crop(wetness.rast, ext(a_x - 15, b_x + 15, a_y - 15, b_y + 15)))
DEM <- scale(crop(DEM.rast, ext(a_x - 15, b_x + 15, a_y - 15, b_y + 15)))
rm(southness.rast, slope.rast, wetness.rast, DEM.rast)


if(N_threshold == "10"){
  linkage_file <- paste0("./code/growth_sim_results/two_stage/linkage/", density, "_density_", noise, "_noise_",
                         alpha, "_alpha_", index, "_index_two_stage_linkage_results.csv")
  latent_file <- paste0("./code/growth_sim_results/two_stage/raw_data/", density, "_density_", noise, "_noise_",
                        alpha, "_alpha_", index, "_index_two_stage_raw_data.RDS")
}else if(N_threshold == "25"){
  linkage_file <- paste0("./code/growth_sim_results/two_stage/linkage/", density, "_density_", noise, "_noise_", alpha, "_alpha_", index, "_index_", N_threshold, "_N_thresh_two_stage_linkage_results.csv")
  latent_file <- paste0("./code/growth_sim_results/two_stage/raw_data/", density, "_density_", noise, "_noise_", alpha, "_alpha_", index, "_index_", N_threshold, "_N_thresh_two_stage_raw_data.RDS")
  
}
# linkage_file <- paste0("./code/growth_sim_results/two_stage/linkage/", density, "_density_", noise, "_noise_",
#                        index, "_index_two_stage_linkage_results.csv")
# latent_file <- paste0("./code/growth_sim_results/two_stage/raw_data/", density, "_density_", noise, "_noise_",
#                       index, "_index_two_stage_raw_data.RDS")
# linkage_file <- paste0("./code/growth_sim_results/two_stage/linkage/", density, "_density_", noise, "_noise_", alpha, "_alpha_", index, "_index_", N_threshold, "_N_thresh_two_stage_linkage_results.csv")
# latent_file <- paste0("./code/growth_sim_results/two_stage/raw_data/", density, "_density_", noise, "_noise_", alpha, "_alpha_", index, "_index_", N_threshold, "_N_thresh_two_stage_raw_data.RDS")

latent_index <- read_csv("./code/growth_sim_results/two_stage/LA_sim_sample_index.csv",
                         col_names = FALSE,
                         show_col_types = FALSE)

linkage_sample <- read_csv(linkage_file, show_col_types = FALSE) %>% as.matrix()
latent_sample <- read_rds(latent_file)


# data_file <- paste0("./code/growth_sim_data_3/", density, "_dens_medium_noise_alpha_3_sim_", index, ".csv")
data_file <- paste0("./code/growth_sim_data_F23/", density, "_dens_", noise,
                    "_noise_", alpha, "_alpha_sim_", index, ".csv")
scan_data <- read_csv(data_file, show_col_types = FALSE)
scan_data <- scan_data %>% filter(!file == 0)


## Obtain the thinned linkage and corresponding latents for the two-stage model
# latent_index <- sample(1:nrow(linkage_sample), 100, replace = FALSE)
if(linkage_set == "1"){
  lambda_configs <- linkage_sample[latent_index$X1[1:50],] + 1
  s_configs <- latent_sample[,,latent_index$X1[1:50]]
}else if(linkage_set == "2"){
  lambda_configs <- linkage_sample[latent_index$X1[51:100],] + 1
  s_configs <- latent_sample[,,latent_index$X1[51:100]]
}
# lambda_configs <- linkage_sample[latent_index,] + 1
# s_configs <- latent_sample[,,latent_index]
N <- dim(s_configs)[2]


## Create storage for two-stage growth model results
growth_results <- list()

## Run the growth model for the thinned lambda configurations with appropriate latents
for(i in 1:nrow(lambda_configs)){
  
  dat <- scan_data
  
  dat$id <- lambda_configs[i,]
  dat$x <- s_configs[c(dat$id), 1, i]
  dat$y <- s_configs[c(dat$id), 2, i]
  
  dat_file_1 <- dat %>% filter(file == 1)
  dat_file_2 <- dat %>% filter(file == 2)
  
  if(length(dat_file_1$id) != length(unique(dat_file_1$id))){
    
    links_within <- list()
    
    for(j in 1:N){
      
      links_within[[j]] <- which(dat_file_1$id == j)
      
    }
    
    if(any(lapply(links_within, length) > 1)){
      
      linked_records <- which(lapply(links_within, length) > 1)
      
      for(k in 1:length(linked_records)){
        
        dat_file_1$size[links_within[[linked_records[k]]]] <- sum(dat_file_1$size[links_within[[linked_records[k]]]])
        
      }
    }
    
    dat_file_1 <- dat_file_1[!duplicated(dat_file_1$size),]
    
  }else if(length(dat_file_2$id) != length(unique(dat_file_2$id))){
    
    links_within <- list()
    
    for(j in 1:N){
      
      links_within[[j]] <- which(dat_file_2$id == j)
      
    }
    
    if(any(lapply(links_within, length) > 1)){
      
      linked_records <- which(lapply(links_within, length) > 1)
      
      for(k in 1:length(linked_records)){
        
        dat_file_2$size[links_within[[linked_records[k]]]] <- sum(dat_file_2$size[links_within[[linked_records[k]]]])
        
      }
    }
    
    dat_file_2 <- dat_file_2[!duplicated(dat_file_2$size),]
    
  }
  
  ## Merge the two files after merging records within the same files
  linked_data <- left_join(dat_file_1, dat_file_2, by = "id")
  linked_data <- linked_data %>% mutate(est_growth = (size.y - size.x)/1,
                                        delta_can = (size.y - size.x)/size.x)
  if(any(is.na(linked_data$file.y)) == TRUE){
    linked_data <- linked_data[-which(is.na(linked_data$file.y)),]
  }
  linked_data$est_mort <- NA
  for(j in 1:nrow(linked_data)){
    if(linked_data$delta_can[[j]] < -.1 | linked_data$delta_can[[j]] > .6){
      linked_data$est_mort[[j]] <- 1
    } else{
      linked_data$est_mort[[j]] <- 0
    }
    
  }
  
  # Obtain the quantities necessary for the growth model
  G <- linked_data$est_growth
  M <- as.numeric(linked_data$est_mort)
  S_1 <- as.numeric(linked_data$size.x)
  S_2 <- as.numeric(linked_data$size.y)
  s <- cbind(linked_data$x.x, linked_data$y.x)
  X <- cbind(terra::extract(southness, s, method = "bilinear"), terra::extract(slope, s, method = "bilinear"),
             terra::extract(wetness, s, method = "bilinear"), terra::extract(DEM, s, method = "bilinear"))
  gc_index <- which(M == 0)
  
  # Run the growth model
  stan_growth_data_2stage = list(N = length(G[gc_index]), # Number of Obs
                                 G = G[gc_index],
                                 mu_0 = rep(0, ncol(X)),
                                 S = S_1[gc_index],
                                 sig20 = diag(rep(2.5, ncol(X))),
                                 K = ncol(X),
                                 X = X[gc_index,],
                                 c = .0001,
                                 d = .0001,
                                 a_gamma = 0,
                                 gamma_max = 1500,
                                 max_tau2 = 75,
                                 a_alpha = 0,
                                 b_alpha = 5) 
  
  # stanfit_2stage <- stan(file = "./code/STAN_models/STAN_growth_mod.stan", # Stan file
  #                        data = stan_growth_data_2stage, # Data
  #                        warmup = 10000, # Number of iteration to burn-in
  #                        iter = 15000, # Total number of iterations
  #                        chains = 4, # Number of chains to run
  #                        thin = 10)
  stanfit_2stage <- stan(file = "./code/STAN_models/STAN_growth_mod_alpha.stan", # Stan file
                         data = stan_growth_data_2stage, # Data
                         warmup = 10000, # Number of iteration to burn-in
                         iter = 15000, # Total number of iterations
                         chains = 4, # Number of chains to run
                         thin = 10)
  
  growth_results[[i]] <- do.call(cbind, rstan::extract(stanfit_2stage, permuted = TRUE))
  print(paste0("Growth Model Iteration ", i, "/", nrow(lambda_configs)))
  
}

growth_results <- do.call(rbind, growth_results)
growth_results <- as.data.frame(growth_results)


## Save the results
# la_results_file <- paste0("./code/growth_sim_results/two_stage/linkage_avg/", density, "_density_", noise, "_noise_", index, "_growth_results_linkage_avg.csv")
# la_results_file <- paste0("./code/growth_sim_results/two_stage/linkage_avg/", density, "_density_", noise, "_noise_", alpha, "_alpha_sim_", index, "_growth_results_linkage_avg.csv")
la_results_file <- paste0("./code/growth_sim_results/two_stage/linkage_avg/", density, "_density_", noise, "_noise_", alpha, "_alpha_sim_", index, "_index_", linkage_set, "_set_", N_threshold, "_N_thresh_growth_results_linkage_avg.csv")


write_csv(as.data.frame(growth_results), file = la_results_file)
