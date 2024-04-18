#########################################################################################################
#### This script runs the growth portion of the two-stage model for the spatial record linkage model ####
#########################################################################################################

# pass from command line ----
# args <- commandArgs(trailingOnly=TRUE)

covars <- "all"
index <- as.numeric(29)
mort_threshold <- "100"


## Load libraries and sampler functions ----
library(rstan) ## for growth model fit
library(readr) ## load and save results
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(terra)
library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)
library(GGally)
library(purrr)
library(gtools)


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
  
  # ## Read in the raster data for the covariates of interest
  # slope.rast <- scale(rast('./data/Snodgrass_slope_1m.tif'))
  # southness.rast <- scale(rast('./data/Snodgrass_aspect_southness_1m.tif'))
  # wetness.rast <- scale(rast('./data/Snodgrass_wetness_index_1m.tif'))
  # DEM.rast <- scale(rast('./data/Snodgrass_DEM_1m.tif'))
  # GDD.rast <- scale(rast('./data/Snodgrass_Degree_Days_2013_2019.tif'))
  # SPP.rast <- scale(rast('./data/Snodgrass_Snowpack_Persistence_DOY_2013_2019.tif'))
  # 
  # 
  # ## Crop the rasters to D* and discard the originals
  # southness <- crop(southness.rast, ext(a_x, b_x, a_y, b_y))
  # slope <- crop(slope.rast, ext(a_x, b_x, a_y, b_y))
  # wetness <- crop(wetness.rast, ext(a_x, b_x, a_y, b_y))
  # dem <- crop(DEM.rast, ext(a_x, b_x, a_y, b_y))
  # spp <- crop(SPP.rast, ext(a_x, b_x, a_y, b_y))
  # gdd <- crop(GDD.rast, ext(a_x, b_x, a_y, b_y))
  # rm(southness.rast, slope.rast, DEM.rast, wetness.rast, SPP.rast, GDD.rast)
  # 
  # raster_list <- c(list(wetness), list(southness), list(slope), list(dem),
  #                  lapply(1:5, function(x) spp[[x]]), lapply(1:5, function(x) gdd[[x]]),
  #                  lapply(1:5, function(x) spp[[x]]*wetness), lapply(1:5, function(x) gdd[[x]]*wetness))
  
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


linkage_file <- paste0("./code/empirical_data/empirical_linkage_lambda_pooled.csv")
latent_file <- paste0("./code/empirical_data/empirical_linkage_s_pooled.csv")

sample_index <- fread(file = "./code/empirical_data/LA_sample_index_pooled.csv", header = FALSE) %>% as.matrix()
current_index <- sample_index[index]


linkage_sample <- fread(file = linkage_file, skip = current_index, header = FALSE, nrows = 1) %>% as.matrix()
latent_sample <- fread(file = latent_file, skip = current_index, header = FALSE, nrows = 1) %>% as.matrix()

# linkage_sample <- fread(file = linkage_file, header = TRUE) %>% as.matrix()
# latent_sample <- fread(file = latent_file, header = TRUE) %>% as.matrix()



# latent_array <- apply(latent_sample, 1, vector2matrix, dim = c(ncol(latent_sample)/2, 2), simplify = FALSE)


# linkage_sample <- linkage_sample[current_index,] + 1
# s_configs <- vector2matrix(latent_sample[current_index,], dim = c(ncol(latent_sample)/2, 2))
# rm(latent_sample)
linkage_sample <- linkage_sample + 1
s_configs <- vector2matrix(latent_sample, dim = c(ncol(latent_sample)/2, 2))
rm(latent_sample)


## Obtain the thinned linkage and corresponding latents for the two-stage model

N <- dim(s_configs)[2]


## Run the growth model for the thinned lambda configurations with appropriate latents
data_per_it <- scan_data[in_bounds_full,]
data_per_it$id <- linkage_sample[in_bounds_full]
lnv_norm <- unsplit(lapply(split(CM[,1], data_per_it$id[1:length(in_bounds)]), mean), data_per_it$id[1:length(in_bounds)])
# and_norm <- unsplit(lapply(split(CM[,2], data_per_it$id[1:length(in_bounds)]), mean), data_per_it$id[1:length(in_bounds)])
rsi_norm <- unsplit(lapply(split(CM[,2], data_per_it$id[1:length(in_bounds)]), mean), data_per_it$id[1:length(in_bounds)])
nd_norm <- unsplit(lapply(split(CM[,3], data_per_it$id[1:length(in_bounds)]), mean), data_per_it$id[1:length(in_bounds)])
data_per_it$x <- s_configs[c(data_per_it$id), 1]
data_per_it$y <- s_configs[c(data_per_it$id), 2]

# data_file_1 <- data_per_it %>%
#   filter(file == 1) %>%
#   mutate(lnv_norm = lnv_norm,
#          and_norm = and_norm,
#          nd_norm = nd_norm)
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



X_new <- cbind(X[,2:3],
               apply(X[,4:8], 1, median),
               apply(X[,9:13], 1, median),
               apply(X[,14:18], 1, median),
               apply(X[,19:23], 1, median),
               linked_data$lnv_norm,
               linked_data$rsi_norm,
               linked_data$nd_norm)
colnames(X_new) <- c("TWI", "Southness", "SPP", "GDD", "SPP_TWI", "GDD_TWI", "LNV", "RSI", "ND")
X_new_df <- as.data.frame(X_new)



model_data <- data.frame(growth = G[which(M == 0)], 
                         size = S_1[which(M == 0)],
                         X_new_df[which(M == 0),])

model_data |> ggpairs(lower = list(continuous = wrap("points", alpha = .1)))

N <- length(which(M == 0))

post2 <- stan_glm(log(growth) ~ log(size) + TWI + Southness + SPP + GDD + SPP_TWI + GDD_TWI + LNV + RSI + ND,
                  data = model_data,
                  family = gaussian(link = "identity"),
                  prior = normal(0, 2.5, autoscale = FALSE),
                  prior_intercept = normal(0, 10),
                  warmup = 10000,
                  iter = 20000,
                  seed = 90210)

cbbPalette <- c("#56B4E9", "#009E73", "#D55E00")

post2_df <- data.frame(growths = c(model_data$growth, exp(predict(post2, model_data, type = "response"))),
                       type = rep(c("observed", "predicted"), each = N))

post2_resids <- data.frame(fitted = exp(predict(post2, model_data, type = "response")),
                           resids = model_data$growth - exp(predict(post2, model_data, type = "response")))

post2_df |> ggplot(aes(x = growths, color = type, group = type)) +
  geom_density() +
  scale_color_manual(values = c(cbbPalette[1], cbbPalette[2])) +
  labs(x = "Annual Growth", y = "Density", title = "Log-Log Linear Growth Model", color = "Density Type") +
  theme_bw() +
  ylim(c(0, .25)) -> log_log_mod_pred

post2_resids |> ggplot(aes(x = fitted, y = resids)) +
  # geom_point(alpha = .01) +
  geom_hex() +
  scale_fill_gradient(limits = c(0, 5000),
                      low = cbbPalette[1], high = cbbPalette[2]) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = "Predicted Annual Growth", y = "Residual", title = "Residuals vs Predicted Values Log-Log Linear Model") +
  xlim(c(0, 30)) + ylim(c(-30, 130)) -> log_log_mod_res


## Define a function to add the filename to each line of the imported dataset
read_plus <- function(flnm) {
  read_csv(flnm, show_col_types = FALSE) %>% 
    mutate(filename = flnm)
}


## Two-stage growth model processing
## Two-stage linkage averaging

## Read in all of the results .csv's into one tibble with the filename attached to each line
emp_res_table <-
  list.files(path = "./code/empirical_data/log_data",
             pattern = "emp_pooled_LA_covars_all_growth_cutoff_100_index_",
             full.names = T) |>
  mixedsort() |>
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./code/empirical_data/log_data/", "", emp_res_table$filename, fixed = TRUE)
myNames <- gsub(".csv", "", myNames, fixed = TRUE)

variables <- do.call(rbind, str_split(myNames, "_"))

emp_res_table <- emp_res_table %>%
  mutate(filename = myNames,
         index = as.numeric(variables[,11])) %>%
  select(-"lp__")



emp_mat <- emp_res_table |> select(-c("filename", "index")) |> as.matrix()

emp_means <- apply(emp_mat, 2, mean)
# emp_means <- apply(emp_mat, 2, median)
emp_quants <- apply(emp_mat, 2, quantile, probs = c(.05, .95))
emp_all <- rbind(emp_means, emp_quants)




N = length(which(M == 0))

X_new_all <- cbind(X[,1:3],
                   apply(X[,4:8], 1, median),
                   apply(X[,9:13], 1, median),
                   apply(X[,14:18], 1, median),
                   apply(X[,19:23], 1, median),
                   linked_data$lnv_norm,
                   linked_data$rsi_norm,
                   linked_data$nd_norm)

pred_means_all <- pred_means(X = X_new_all[which(M == 0),], beta = emp_means[3:12], gamma = emp_means[1],
                             alpha = emp_means[13], V = S_1[which(M == 0)], N = length(which(M == 0)))
# all_growths_df <- data.frame(growths = c(G[which(M == 0)], exp(pred_means_all)), type = rep(c("observed", "predicted"), each = N))
all_growths_df <- data.frame(growths = c(G[which(M == 0)], exp(pred_means_all)),
                             type = rep(c("observed", "predicted"), each = N))
all_growths_df |> ggplot(aes(x = growths, color = type, group = type)) +
  scale_color_manual(values = c(cbbPalette[1], cbbPalette[2])) +
  geom_density() +
  # labs(x = "Annual Growth", y = "Density", title = "Density Comparison w/ Competition Metrics", color = "Density Type") +
  labs(x = "Annual Growth", y = "Density", title = "Log Growth MM Model", color = "Density Type") +
  theme_bw() +
  ylim(c(0, .25)) -> all_dens

all_df <- data.frame(obs_growth = G[which(M == 0)], pred_growth = exp(pred_means_all), residuals = G[which(M == 0)] - exp(pred_means_all))
all_df |> ggplot(aes(x = pred_growth, y = residuals)) + 
  # geom_point(alpha = .5) +
  geom_hex() +
  scale_fill_gradient(limits = c(0, 5000),
                      low = cbbPalette[1], high = cbbPalette[2]) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = "Predicted Annual Growth", y = "Residual", title = "Residuals vs Predicted Values Log MM Model") +
  xlim(c(0, 30)) + ylim(c(-30, 130)) -> all_res


emp_res_table <-
  list.files(path = "./code/empirical_data/",
             pattern = "emp_pooled_LA_covars_all_growth_cutoff_100_index_",
             full.names = T) |>
  mixedsort() |>
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./code/empirical_data/", "", emp_res_table$filename, fixed = TRUE)
myNames <- gsub(".csv", "", myNames, fixed = TRUE)

variables <- do.call(rbind, str_split(myNames, "_"))

emp_res_table <- emp_res_table %>%
  mutate(filename = myNames,
         index = as.numeric(variables[,10])) %>%
  select(-"lp__")



emp_mat <- emp_res_table |> select(-c("filename", "index")) |> as.matrix()

emp_means <- apply(emp_mat, 2, mean)
# emp_means <- apply(emp_mat, 2, median)
emp_quants <- apply(emp_mat, 2, quantile, probs = c(.05, .95))
emp_all <- rbind(emp_means, emp_quants)




N = length(which(M == 0))

X_new_all <- cbind(X[,1:3],
                   apply(X[,4:8], 1, median),
                   apply(X[,9:13], 1, median),
                   apply(X[,14:18], 1, median),
                   apply(X[,19:23], 1, median),
                   linked_data$lnv_norm,
                   linked_data$rsi_norm,
                   linked_data$nd_norm)

pred_means_all <- pred_means(X = X_new_all[which(M == 0),], beta = emp_means[3:12], gamma = emp_means[1],
                             alpha = emp_means[13], V = S_1[which(M == 0)], N = length(which(M == 0)))
all_growths_df <- data.frame(growths = c(G[which(M == 0)], pred_means_all), type = rep(c("observed", "predicted"), each = N))
all_growths_df |> ggplot(aes(x = growths, color = type, group = type)) +
  geom_density() +
  scale_color_manual(values = c(cbbPalette[1], cbbPalette[2])) +
  # labs(x = "Annual Growth", y = "Density", title = "Density Comparison w/ Competition Metrics", color = "Density Type") +
  labs(x = "Annual Growth", y = "Density", title = "Growth MM Model", color = "Density Type") +
  theme_bw() +
  ylim(c(0, .25)) -> all_dens2

all_df <- data.frame(obs_growth = G[which(M == 0)], pred_growth = pred_means_all, residuals = G[which(M == 0)] - pred_means_all)
all_df |> ggplot(aes(x = pred_growth, y = residuals)) + 
  # geom_point(alpha = .5) +
  geom_hex() +
  scale_fill_gradient(limits = c(0, 5000),
                      low = cbbPalette[1], high = cbbPalette[2]) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = "Predicted Annual Growth", y = "Residual", title = "Residuals vs Predicted Values MM Model") +
  xlim(c(0, 30)) + ylim(c(-30, 130)) -> all_res2




results_file <- paste0("./code/empirical_data/gamma_MM_model/gamma_emp_pooled_LA_covars_all_growth_cutoff_100_index_29.csv")

growth_res <- fread(results_file, header = TRUE)

growth_summary <- growth_res |>
  select(c("beta_0", "betas.1", "betas.2", "betas.3", "betas.4", "betas.5",
           "betas.6", "betas.7", "betas.8", "betas.9", "gamma", "alpha", "phi")) |> 
  as.matrix()

means_growth <- apply(growth_summary, 2, mean)
# means_growth <- apply(growth_summary, 2, median)
q_5 <- apply(growth_summary, 2, quantile, probs = c(.05, .95))
q_95 <- apply(growth_summary, 2, quantile, probs = .95)

emp_all <- rbind(means_growth, q_5)

colnames(emp_all) <- c("beta_0",
                       "TWI",
                       "southness",
                       "SPP",
                       "GDD",
                       "TWI_SPP",
                       "TWI_GDD",
                       "LNV",
                       "RSI",
                       "ND",
                       "gamma",
                       "alpha",
                       "phi")


emp_proc <- data.frame(t(emp_all)) |>
  mutate(variable = colnames(emp_all)) |>
  mutate(variable = factor(variable, levels = rev(c("beta_0",
                                                    "TWI",
                                                    "southness",
                                                    "SPP",
                                                    "GDD",
                                                    "TWI_SPP",
                                                    "TWI_GDD",
                                                    "LNV",
                                                    "RSI",
                                                    "ND",
                                                    "gamma",
                                                    "alpha",
                                                    "phi"))),
         coverage = NA)


for(i in 1:nrow(emp_proc)){
  if(0 > emp_proc$X5.[i] & 0 < emp_proc$X95.[i]){
    emp_proc$coverage[i] <- 0
  }else{
    emp_proc$coverage[i] <- 1
  }
}


N = length(which(M == 0))

pred_means_all <- pred_means(X = X_new_all[which(M == 0),], beta = means_growth[1:10], gamma = means_growth[11],
                             alpha = means_growth[12], V = S_1[which(M == 0)], N = length(which(M == 0)))
all_growths_df <- data.frame(growths = c(G[which(M == 0)], exp(pred_means_all)), type = rep(c("observed", "predicted"), each = N))
all_growths_df |> ggplot(aes(x = growths, color = type, group = type)) +
  geom_density() +
  scale_color_manual(values = c(cbbPalette[1], cbbPalette[2])) +
  labs(x = "Annual Growth", y = "Density", title = "Gamma MM Model", color = "Density Type") +
  theme_bw() +
  ylim(c(0, .25)) -> p1

all_df <- data.frame(obs_growth = G[which(M == 0)], pred_growth = exp(pred_means_all), residuals = G[which(M == 0)] - exp(pred_means_all))
all_df |> ggplot(aes(x = pred_growth, y = residuals)) + 
  geom_hex() +
  scale_fill_gradient(limits = c(0, 5000),
                      low = cbbPalette[1], high = cbbPalette[2]) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  labs(x = "Predicted Annual Growth", y = "Residuals",
       title = "Residuals vs Predicted Values Gamma MM Model") +
  xlim(c(0, 30)) + ylim(c(-30, 130)) -> p2



(all_dens2 + all_dens + p1 + log_log_mod_pred + plot_layout(guides = "collect", ncol = 4) &
    theme(legend.position = "right")) / 
  (all_res2 + all_res + p2 + log_log_mod_res + plot_layout(guides = "collect", ncol = 4) &
     theme(legend.position = "right",
           plot.title = element_blank()))







