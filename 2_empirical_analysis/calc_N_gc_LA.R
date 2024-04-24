## Load libraries ----
library(tidyr) ## data manip
library(readr) ## load and save results
library(dplyr) ## data manip
library(codetools)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)

sourceCpp('./code/cpp_code/two_stage_func.cpp')

## Set seed for reproducibility ----
set.seed(90210)

settings_grid <- expand.grid(c("low", "med", "high"), c("small", "medium", "large"), 1, 1:100)
colnames(settings_grid) <- c("density", "noise", "alpha", "index")


N_gc_LA <- matrix(NA, nrow = nrow(settings_grid), ncol = 100)

for(m in 1:nrow(settings_grid)){
  
  density <- settings_grid$density[m]
  noise <- settings_grid$noise[m]
  alpha <- settings_grid$alpha[m]
  index <- settings_grid$index[m]
  
  print(paste0("LA Density: ", density, ", Noise: ", noise, ", Alpha: ", alpha, ", Index: ", index))
  
  linkage_file <- paste0("./code/growth_sim_results/two_stage/linkage/", density, "_density_", noise, "_noise_", alpha, "_alpha_", index, "_index_25_N_thresh_two_stage_linkage_results.csv")
  latent_file <- paste0("./code/growth_sim_results/two_stage/raw_data/", density, "_density_", noise, "_noise_", alpha, "_alpha_", index, "_index_25_N_thresh_two_stage_raw_data.RDS")
  
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
  lambda_configs <- linkage_sample[latent_index$X1,] + 1
  s_configs <- latent_sample[,,latent_index$X1]
  
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
    gc_index <- which(M == 0)
    
    N_gc_LA[m,i] <- length(gc_index)
    
  }
  
}

write_csv(as.data.frame(N_gc_LA), file = "./code/growth_sim_results/N_gc_LA.csv")