####################################################################
#### This script generates Figure 2 from the paper supplement. #####
####################################################################

## Load libraries
library(readr)
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
library(latex2exp)

## Set a colorblind friendly palate to use for visualizations
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Set necessary arguments
index <- as.numeric(29)
mort_threshold <- "90"
model_type <- c("skew_t", "skew_normal", "normal", "MLR")

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

in_bounds <- which(file2015$x > a_x2 & file2015$x < b_x2 & file2015$y > a_y2 & file2015$y < b_y2)
scan_data <- rbind(file2015, file2019)

in_bounds_full <- which((scan_data$x > a_x2 & scan_data$x < b_x2 & scan_data$y > a_y2 & scan_data$y < b_y2 & scan_data$file == 1) | (scan_data$file == 2))


#### RL ####
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
gc_index <- which(M == 0)

RL_dens_list <- list()
RL_dens_list[["observed"]] <- data.frame(growth = G[gc_index], replicate = rep("observed", sum(gc_index)))

for(i in model_type){
  
  replicate_data <- read_csv(paste0("./2_empirical_analysis/model_results/growth_model/LA/emp_pooled_LA_N_25_model_", i,
                                    "_growth_cutoff_", mort_threshold, "_index_", index, "_rep1.csv"))
  RL_dens_list[[i]] <- data.frame(growth = replicate_data[1000,], replicate = rep(i, ncol(replicate_data)))
  
}

dens_comp_rl <- do.call(cbind, RL_dens_list) |> 
  mutate(replicate = factor(type,
                            levels = c("observed", "skew_t", "skew_normal", "normal", "MLR"),
                            labels = c("Observed", "Skewed t", "Skew Normal", "Normal", "MLR")))
dens_comp_rl |> ggplot(aes(x = growth, group = as.factor(replicate), color = replicate)) +
  geom_density() +
  theme_bw() +
  scale_color_manual(values = cbbPalette) +
  labs(x = TeX("Annual Growth $(m^3)$"), y = "Density",
       title = "Record Linkage",
       color = "Replicate Type") +
  theme(text = element_text(family = "serif", size = 14)) -> rl_dens_plot


#### POM ####
POM_growth <- read_csv(paste0("./2_empirical_analysis/model_results/growth_model/POM/polygon_overlap_matching_growth_cutoff_", mort_threshold, "_growths.csv"))

POM_dens_list <- list()
POM_dens_list[["observed"]] <- data.frame(growth = POM_growth$growth, replicate = rep("observed", length(POM_growth$growth)))

for(i in model_type){
  
  replicate_data <- read_csv(paste0("./2_empirical_analysis/model_results/growth_model/POM/polygon_overlap_matching_model_",
                                    model_type, "_growth_cutoff_", mort_threshold, "_rep1.csv"))
  POM_dens_list[[i]] <- data.frame(growth = replicate_data[1000,], replicate = rep(i, ncol(replicate_data)))
  
}

dens_comp_pom <- do.call(cbind, POM_dens_list) |> 
  mutate(replicate = factor(type,
                            levels = c("observed", "skew_t", "skew_normal", "normal", "MLR"),
                            labels = c("Observed", "Skewed t", "Skew Normal", "Normal", "MLR")))
dens_comp_pom |> ggplot(aes(x = growth, group = as.factor(replicate), color = replicate)) +
  geom_density() +
  theme_bw() +
  scale_color_manual(values = cbbPalette) +
  labs(x = TeX("Annual Growth $(m^3)$"), y = "Density",
       title = "Polygon Overlap Matching",
       color = "Replicate Type") +
  theme(text = element_text(family = "serif", size = 14)) -> pom_dens_plot


#### NDM ####
NDM_growth <- read_csv(paste0("./2_empirical_analysis/model_results/growth_model/NDM/nearest_distance_matching_growth_cutoff_", mort_threshold, "_growths.csv"))

NDM_dens_list <- list()
NDM_dens_list[["observed"]] <- data.frame(growth = NDM_growth$growth, replicate = rep("observed", length(NDM_growth$growth)))

for(i in model_type){
  
  replicate_data <- read_csv(paste0("./2_empirical_analysis/model_results/growth_model/NDM/polygon_overlap_matching_model_",
                                    model_type, "_growth_cutoff_", mort_threshold, "_rep1.csv"))
  NDM_dens_list[[i]] <- data.frame(growth = replicate_data[1000,], replicate = rep(i, ncol(replicate_data)))
  
}

dens_comp_ndm <- do.call(cbind, NDM_dens_list) |> 
  mutate(replicate = factor(type,
                            levels = c("observed", "skew_t", "skew_normal", "normal", "MLR"),
                            labels = c("Observed", "Skewed t", "Skew Normal", "Normal", "MLR")))
dens_comp_ndm |> ggplot(aes(x = growth, group = as.factor(replicate), color = replicate)) +
  geom_density() +
  theme_bw() +
  scale_color_manual(values = cbbPalette) +
  labs(x = TeX("Annual Growth $(m^3)$"), y = "Density",
       title = "Nearest Distance Matching",
       color = "Replicate Type") +
  theme(text = element_text(family = "serif", size = 14)) -> pom_dens_plot


## Generate and save the four panel density comparison plot (Figure 2)
all_dens_plot <- rl_dens_plot + pom_dens_plot + ndm_dens_plot + guide_area() + 
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a', tag_suffix = ")")

ggsave(filename = "10_replicated_density_comparison_plot.png", plot = all_dens_plot, path = "./3_figures_and_tables/",
       width = 20, height = 20, units = "cm", dpi = "retina")
