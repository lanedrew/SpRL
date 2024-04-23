# pass from command line ----
# Rscript code/scripts/r/SpRL_sim_study.R
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) stop("Pass in density (low or med or high), noise (small or medium or large),
                            alpha (1 or 2 or 3), and index (integer from 1 to 100).", call.=FALSE)
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
library(data.table)
library(rstan)

rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())

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
mort_threshold <- "90"
sigma2.prior <- c("uninformative")
density <- args[1]
noise <- args[2]
alpha <- args[3]
index <- args[4]

print(paste0("NDM Density: ", density, ", Noise: ", noise, ", Alpha: ", alpha, ", Index: ", index))

## Read in the specified dataset
filename <- paste0("./code/growth_sim_data_F23/", density, "_dens_", noise, "_noise_", alpha, "_alpha_sim_", index, ".csv")
scan_data <- read_csv(filename)
scan_data <- scan_data %>% filter(!file == 0)

file1_mat <- scan_data %>% filter(file == 1) %>% select(x, y) %>% as.matrix()
file2_mat <- scan_data %>% filter(file == 2) %>% select(x, y) %>% as.matrix()

lambda <- nearest_distance_matching(file_1 = file1_mat, file_2 = file2_mat, dist = 2)

linkage_file <- paste0("./code/growth_sim_results/NDM/linkage/", density, "_dens_", noise, "_noise_", alpha, "_alpha_sim_", index, "_NDM_linkage_results.csv")
write_csv(as.data.frame(lambda), linkage_file)

linkage_res <- do.call(cbind, eval_links(z = t(lambda), true_id = scan_data$id))
write_csv(as.data.frame(linkage_res),
          paste0("./code/growth_sim_results/NDM/linkage_processed/",
                 density, "_density_", noise, "_noise_", alpha, "_alpha_", index,
                 "_index_NDM_linkage_metrics.csv"))


## Obtain the file sizes and cumulative index m
# file_size <- scan_data %>% 
#   dplyr::select(file) %>%
#   dplyr::group_by(file) %>%
#   dplyr::summarise(n_i = n()) %>%
#   dplyr::select(n_i)
# file_size <- unlist(file_size)
# 
# N <- floor(1.1*max(file_size))


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
southness <- scale(crop(southness.rast, ext(a_x - 15, b_x + 15, a_y - 15, b_y + 15)))
slope <- scale(crop(slope.rast, ext(a_x - 15, b_x + 15, a_y - 15, b_y + 15)))
wetness <- scale(crop(wetness.rast, ext(a_x - 15, b_x + 15, a_y - 15, b_y + 15)))
DEM <- scale(crop(DEM.rast, ext(a_x - 15, b_x + 15, a_y - 15, b_y + 15)))
rm(southness.rast, slope.rast, wetness.rast, DEM.rast)


data_per_it <- scan_data %>% select(x, y, size, file) %>% mutate(id = lambda)

data_per_it$x <- unsplit(lapply(split(data_per_it$x, data_per_it$id), mean), data_per_it$id)
data_per_it$y <- unsplit(lapply(split(data_per_it$y, data_per_it$id), mean), data_per_it$id)

data_file_1 <- data_per_it %>%
  filter(file == 1)
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
linked_data <- linked_data %>% mutate(est_growth = (size.y - size.x)/1,
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


stanfit_2stage <- stan(file = "./code/STAN_models/STAN_growth_mod_alpha.stan", # Stan file
                       data = stan_growth_data_2stage, # Data
                       warmup = 10000, # Number of iteration to burn-in
                       iter = 25000, # Total number of iterations
                       chains = 4, # Number of chains to run
                       thin = 10)

growth_results <- rstan::extract(stanfit_2stage, permuted = TRUE)
growth_results <- as.data.frame(growth_results)

ndm_results_file <- paste0("./code/growth_sim_results/NDM/growth_model_ests/", density, "_density_", noise, "_noise_", alpha, "_alpha_sim_", index, "_growth_results_NDM.csv")

write_csv(growth_results, file = ndm_results_file)
