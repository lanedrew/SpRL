########################################################
#### This script generates Figure 6 from the paper. ####
########################################################

## Load the relevant packages
library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(gtools)
library(latex2exp)
library(purrr)
library(stringr)
library(tidyr)
library(patchwork)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(terra)
library(data.table)
library(GGally)
library(latex2exp)

## Load C++ helper functions
sourceCpp('./resources/code/cpp_code/two_stage_func.cpp')

## Set a colorblind friendly palate to use for visualizations
cbbPalette <- c("#56B4E9", "#009E73", "#D55E00")

## Define a function to add the filename to each line of the imported dataset
read_plus <- function(flnm) {
  read_csv(flnm, show_col_types = FALSE) %>% 
    mutate(filename = flnm)
}


## Set seed for reproducibility ----
set.seed(90210)

## Two-stage linkage averaging growth results
## Read in all of the results .csv's into one tibble with the filename attached to each line
emp_res_table <-
  list.files(path = "./2_empirical_analysis/model_results/growth_model/LA",
             pattern = "model_skew_t",
             full.names = T) |>
  mixedsort() |>
  map_df(~read_plus(.))


## Modify the names from the files
myNames <- gsub("./2_empirical_analysis/model_results/growth_model/LA/", "", emp_res_table$filename, fixed = TRUE)
myNames <- gsub(".csv", "", myNames, fixed = TRUE)
variables <- do.call(rbind, str_split(myNames, "_"))

emp_res_table <- emp_res_table %>%
  mutate(filename = myNames,
         index = as.numeric(variables[,15])) %>%
  select(-"lp__") %>%
  rename(delta = lambda,
         omega = q)

emp_mat <- emp_res_table |> select(-c("filename", "index")) |> as.matrix()

emp_means <- apply(emp_mat, 2, mean)
emp_quants <- apply(emp_mat, 2, quantile, probs = c(.05, .95))

emp_all <- rbind(emp_means, emp_quants)

colnames(emp_all) <- c("gamma", "tau", "beta_0",
                       "TWI",
                       "southness",
                       "SP",
                       "GDD",
                       "TWI_SP",
                       "TWI_GDD",
                       "LNV",
                       "RSI",
                       "ND",
                       "alpha",
                       "delta",
                       "omega")

emp_proc <- data.frame(t(emp_all)) |>
  mutate(variable = colnames(emp_all)) |>
  mutate(variable = factor(variable, levels = rev(c("gamma", "tau", "beta_0", "TWI",
                                                    "southness",
                                                    "SP",
                                                    "GDD",
                                                    "TWI_SP",
                                                    "TWI_GDD",
                                                    "LNV",
                                                    "RSI",
                                                    "ND",
                                                    "alpha",
                                                    "delta",
                                                    "omega"))))


## Specify necessary arguments
index <- as.numeric(29)
mort_threshold <- "90"

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
  select(XTOP, YTOP, CANVOL2015, ZTOP) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2015, height = ZTOP) %>%
  mutate(file = 1)
file2019 <- read_csv("./resources/empirical_data/UER_lidar_canopy_segmentation/crown_attributes_2019.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(XTOP, YTOP, CANVOL2019, ZTOP) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2019, height = ZTOP) %>%
  mutate(file = 2)

in_bounds <- which(file2015$x > a_x2 & file2015$x < b_x2 & file2015$y > a_y2 & file2015$y < b_y2)
scan_data <- rbind(file2015, file2019)

in_bounds_full <- which((scan_data$x > a_x2 & scan_data$x < b_x2 & scan_data$y > a_y2 & scan_data$y < b_y2 & scan_data$file == 1) | (scan_data$file == 2))

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
  
  
linkage_file <- paste0("./2_empirical_analysis/model_results/record_linkage_model/empirical_linkage_lambda_pooled_N_25.csv")
latent_file <- paste0("./2_empirical_analysis/model_results/record_linkage_model/empirical_linkage_s_pooled_N_25.csv")

sample_index <- fread(file = "./2_empirical_analysis/empirical_data/LA_sample_index_pooled_N_25.csv", header = FALSE) %>% as.matrix()
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
      data_file_1$height[links_within[[linked_records[k]]]] <- max(data_file_1$height[links_within[[linked_records[k]]]])
      
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
      data_file_2$height[links_within[[linked_records[k]]]] <- max(data_file_2$height[links_within[[linked_records[k]]]])
      
    }
  }
  
  data_file_2 <- data_file_2[!duplicated(data_file_2$id),]
  
}

## Merge the two files after merging records within the same files
linked_data <- inner_join(data_file_1, data_file_2, by = "id")
linked_data <- linked_data %>% mutate(est_growth = (size.y - size.x)/4,
                                      delta_can = (size.y - size.x)/size.x,
                                      log_est_growth = log(est_growth),
                                      est_growth_height = (height.y - height.x)/4,
                                      delta_height = (height.y - height.x)/height.x)


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

X_new <- cbind(rep(1, length(G)),
               X[,2:3],
               apply(X[,4:8], 1, median),
               apply(X[,9:13], 1, median),
               apply(X[,14:18], 1, median),
               apply(X[,19:23], 1, median),
               linked_data$lnv_norm,
               linked_data$rsi_norm,
               linked_data$nd_norm)[gc_index,]


## Generate predictions and plot for SP
low_spp_obs <- c(apply(X_new[,1:3], 2, median), quantile(X_new[,4], .2), apply(X_new[,-c(1:4)], 2, median))
high_spp_obs <- c(apply(X_new[,1:3], 2, median), quantile(X_new[,4], .8), apply(X_new[,-c(1:4)], 2, median))
low_spp_obs[6] <- low_spp_obs[2]*low_spp_obs[4]
high_spp_obs[6] <- high_spp_obs[2]*high_spp_obs[4]


## Create a function to predict growth given a set of parameters and covariates and size
growth_func <- function(post_summary, obs_vals, size){
  
  growth <- (post_summary[,3:12] %*% as.matrix(obs_vals) * size^post_summary[1,13]) / (post_summary[1,1]^post_summary[1,13] + size^post_summary[1,13])
  return(t(growth))
  
}


size_vec <- seq(min(S_1[gc_index]), max(S_1[gc_index]), by = .1)
growth_pred_low <- data.frame(do.call(rbind, lapply(size_vec, function(x) growth_func(emp_all, low_spp_obs, x))), sizes = size_vec)
growth_pred_high <- data.frame(do.call(rbind, lapply(size_vec, function(x) growth_func(emp_all, high_spp_obs, x))), sizes = size_vec)

growth_pred_all <- data.frame(rbind(growth_pred_low,
                                    growth_pred_high),
                              dataset = rep(c("low_spp", "high_spp"),
                                            each = length(size_vec))) |> 
  mutate(dataset = factor(dataset, 
                          levels = c("low_spp", "high_spp"),
                          labels = c("20th Quantile", "80th Quantile")))

growth_pred_all |> ggplot(aes(x = sizes, y = emp_means, group = dataset, color = dataset)) +
  geom_line() +
  geom_line(aes(y = X5.), linetype = "dashed", alpha = .5) +
  geom_line(aes(y = X95.), linetype = "dashed", alpha = .5) +
  scale_color_manual(values = cbbPalette[-3]) +
  ylim(0, 7) +
  labs(y = TeX('Annual Canopy Volume Growth in $m^3'), x = TeX('Canopy Volume in $m^3'),
       title = "",
       color = "Covariate Quantile") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "serif", size = 14)) -> spp_plot


## Generate predictions and plot for GDD
low_gdd_obs <- c(apply(X_new[,1:4], 2, median), quantile(X_new[,5], .2), apply(X_new[,-c(1:5)], 2, median))
high_gdd_obs <- c(apply(X_new[,1:4], 2, median), quantile(X_new[,5], .8), apply(X_new[,-c(1:5)], 2, median))
low_gdd_obs[7] <- low_gdd_obs[2]*low_gdd_obs[5]
high_gdd_obs[7] <- high_gdd_obs[2]*high_gdd_obs[5]


size_vec <- seq(min(S_1[gc_index]), max(S_1[gc_index]), by = .1)
growth_pred_low <- data.frame(do.call(rbind, lapply(size_vec, function(x) growth_func(emp_all, low_gdd_obs, x))), sizes = size_vec)
growth_pred_high <- data.frame(do.call(rbind, lapply(size_vec, function(x) growth_func(emp_all, high_gdd_obs, x))), sizes = size_vec)


growth_pred_all <- data.frame(rbind(growth_pred_low,
                                    growth_pred_high),
                              dataset = rep(c("low_gdd", "high_gdd"),
                                            each = length(size_vec))) |> 
  mutate(dataset = factor(dataset, 
                          levels = c("low_gdd", "high_gdd"),
                          labels = c("20th Quantile", "80th Quantile")))


growth_pred_all |> ggplot(aes(x = sizes, y = emp_means, group = dataset, color = dataset)) +
  geom_line() +
  geom_line(aes(y = X5.), linetype = "dashed", alpha = .5) +
  geom_line(aes(y = X95.), linetype = "dashed", alpha = .5) +
  scale_color_manual(values = cbbPalette[-3]) +
  ylim(0, 7) +
  labs(y = TeX('Annual Canopy Volume Growth in $m^3'), x = TeX('Canopy Volume in $m^3'),
       title = "",
       color = "Covariate Quantile") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "serif", size = 14)) -> gdd_plot

## Generate predictions and plot for GDD/SPP interactions
low_gdd_spp_obs <- c(apply(X_new[,1:3], 2, median), quantile(X_new[,4], .2), quantile(X_new[,5], .2), apply(X_new[,-c(1:5)], 2, median))
high_gdd_spp_obs <- c(apply(X_new[,1:3], 2, median), quantile(X_new[,4], .8), quantile(X_new[,5], .8), apply(X_new[,-c(1:5)], 2, median))
low_gdd_spp_obs[6] <- low_gdd_spp_obs[2]*low_gdd_spp_obs[4]
high_gdd_spp_obs[6] <- high_gdd_spp_obs[2]*high_gdd_spp_obs[4]
low_gdd_spp_obs[7] <- low_gdd_spp_obs[2]*low_gdd_spp_obs[5]
high_gdd_spp_obs[7] <- high_gdd_spp_obs[2]*high_gdd_spp_obs[5]


size_vec <- seq(min(S_1[gc_index]), max(S_1[gc_index]), by = .1)
growth_pred_low <- data.frame(do.call(rbind, lapply(size_vec, function(x) growth_func(emp_all, low_gdd_spp_obs, x))), sizes = size_vec)
growth_pred_high <- data.frame(do.call(rbind, lapply(size_vec, function(x) growth_func(emp_all, high_gdd_spp_obs, x))), sizes = size_vec)


growth_pred_all <- data.frame(rbind(growth_pred_low,
                                    growth_pred_high),
                              dataset = rep(c("low_gdd_spp", "high_gdd_spp"),
                                            each = length(size_vec))) |> 
  mutate(dataset = factor(dataset, 
                          levels = c("low_gdd_spp", "high_gdd_spp"),
                          labels = c("Dry/Cold", "Wet/Hot")))

growth_pred_all |> ggplot(aes(x = sizes, y = emp_means, group = dataset, color = dataset)) +
  geom_line() +
  geom_line(aes(y = X5.), linetype = "dashed", alpha = .5) +
  geom_line(aes(y = X95.), linetype = "dashed", alpha = .5) +
  scale_color_manual(values = cbbPalette[-3]) +
  ylim(0, 8) +
  labs(y = TeX('Annual Canopy Volume Growth in $m^3'), x = TeX('Canopy Volume in $m^3'),
       title = "Combined Growth Impact of SP and GDD",
       color = "Combined Conditions") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "serif", size = 14)) -> gdd_spp_plot

## Create the three panel growth curve plot and save it (Figure 6)
growth_comp_plot <- spp_plot + gdd_plot + gdd_spp_plot +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") &
  theme(legend.position='bottom')

ggsave(filename = "6_predicted_growth_curves_plot.png", plot = growth_comp_plot, path = "./3_figures_and_tables/",
       width = 36, height = 14, units = "cm", dpi = "retina")
