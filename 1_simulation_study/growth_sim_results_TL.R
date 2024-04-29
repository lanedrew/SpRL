########################################################
#### This script runs the true linkage growth model ####
########################################################

# pass from command line ----
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) stop("Pass in density (low or med or high), noise (low or medium or high), alpha level (1 or 2 or 3), and index (integer from 1 to 100)", call.=FALSE)
if (!(args[1] %in% c("low", "med", "high"))) stop("Pass in the density (low or med or high)", call.=FALSE)
if (!(args[2] %in% c("small", "medium", "large"))) stop("Pass in the noise level (small or medium or large)", call.=FALSE)


## Load libraries ----
library(tidyr) ## data manip
library(readr) ## load and save results
library(dplyr) ## data manip
library(codetools)
library(rstan) ## for growth model fit
library(terra) ## for manipulating raster data

rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())

## Read in the raster data for the covariates of interest
slope_rast <- rast('./data/Snodgrass_slope_1m.tif')
southness_rast <- rast('./data/Snodgrass_aspect_southness_1m.tif')
wetness_rast <- rast('./data/Snodgrass_wetness_index_1m.tif')
DEM_rast <- rast('./data/Snodgrass_DEM_1m.tif')

## Set seed for reproducibility ----
set.seed(90210)

density <- args[1]
noise <- args[2]
alpha <- args[3]
index <- args[4]

print(paste0("TL Density: ", density, ", Noise: ", noise, ", Alpha: ", alpha, ", Index: ", index))

# filename <- paste0("./code/growth_sim_data_3/", density, "_dens_sim_", index, ".csv")
# filename <- paste0("./code/growth_sim_data_3/", density, "_dens_", tau2, "_tau2_sim_", index, ".csv")
filename <- paste0("./code/growth_sim_data_F23/", density, "_dens_", noise, "_noise_", alpha, "_alpha_sim_", index, ".csv")

og_data <- read_csv(filename, show_col_types = FALSE)

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
epsilon <- 15
a_x_star <- a_x - epsilon
a_y_star <- a_y - epsilon
b_x_star <- b_x + epsilon
b_y_star <- b_y + epsilon

## Crop the rasters to D* and discard the originals
southness <- scale(crop(southness_rast, ext(a_x_star, b_x_star, a_y_star, b_y_star)))
slope <- scale(crop(slope_rast, ext(a_x_star, b_x_star, a_y_star, b_y_star)))
wetness <- scale(crop(wetness_rast, ext(a_x_star, b_x_star, a_y_star, b_y_star)))
DEM <- scale(crop(DEM_rast, ext(a_x_star, b_x_star, a_y_star, b_y_star)))
rm(southness_rast, slope_rast, wetness_rast, DEM_rast)

## Run the growth model for the thinned lambda configurations with appropriate latents

dat <- og_data

linked_data <-
  left_join(dat[which(dat$file == 1), ], dat[which(dat$file == 2), ], by = "id")
linked_data <- linked_data %>% 
                  mutate(est_growth = (size.y - size.x) / 1,
                  delta_can = (size.y - size.x) / size.x)

if (any(is.na(linked_data$file.y)) == TRUE) {
  linked_data <- linked_data[-which(is.na(linked_data$file.y)), ]
}

linked_data$est_mort <- NA
for (j in 1:nrow(linked_data)) {
  if (linked_data$delta_can[[j]] < -.1 | linked_data$delta_can[[j]] > .6) {
    linked_data$est_mort[[j]] <- 1
  } else{
    linked_data$est_mort[[j]] <- 0
  }
}

true_loc <- dat[, c(1, 2, 4)][which(dat$file == 0 & dat$id %in% linked_data$id),]


G <- linked_data$est_growth
M <- as.numeric(linked_data$est_mort)
S_1 <- as.numeric(linked_data$size.x)
s <- as.matrix(true_loc[, 1:2])

covars_2015 <- cbind(terra::extract(southness, s, method = "bilinear"),
                     terra::extract(slope, s, method = "bilinear"),
                     terra::extract(wetness, s, method = "bilinear"),
                      terra::extract(DEM, s, method = "bilinear"))
X <- as.matrix(covars_2015)
gc_index <- which(M == 0)

stan_growth_data_true_linkage = list(N = length(G[gc_index]), # Number of Obs
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



stanfit_true_linkage <- stan(file = "./code/STAN_models/STAN_growth_mod_alpha.stan", # Stan file
                             data = stan_growth_data_true_linkage, # Data
                             warmup = 10000, # Number of iteration to burn-in
                             iter = 25000, # Total number of iterations
                             chains = 4, # Number of chains to run
                             thin = 10)

growth_results <- as.data.frame(do.call(cbind, rstan::extract(stanfit_true_linkage, permuted = TRUE)))


## Save the results
results_file <- paste0("./code/growth_sim_results/true_linkage/", density, "_density_", noise, "_noise_", alpha, "_alpha_", index, "_growth_results_TL.csv")
write_csv(growth_results, file = results_file)

