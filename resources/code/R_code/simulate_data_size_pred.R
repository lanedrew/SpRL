## Load the necessary data simulation and naive sampler functions
source('./code/data_simulation/data_sim_functions.R')
library(readr) ## load and save results
library(terra)
library(xgboost)
library(parsnip)
library(tmvtnorm) ## for generating tmvn random variables
library(spatstat)
library(mvtnorm)

## Read in the raster data for the covariates of interest
slope_rast <- rast('./data/Snodgrass_slope_1m.tif')
southness_rast <- rast('./data/Snodgrass_aspect_southness_1m.tif')
wetness_rast <- rast('./data/Snodgrass_wetness_index_1m.tif')
DEM_rast <- rast('./data/Snodgrass_DEM_1m.tif')

## Define a function to simulate the data
sim_data <- function(sigma2_t, theta, density = c("low", "med", "high"), noise = c("small", "medium", "large"),
                     n.cycles, tau2, betas, b, phi = NULL, eta = NULL, gammas = NULL, alpha){
  
  ## Define the spatial domain based on the density and relevant interaction radius and rate of decay for the marked point process
  if(density == "low"){
    
    a_x <- 326496
    a_y <- 4311439
    b_x <- 326596
    b_y <- 4311539
    
    sizes_df <- as.data.frame(read_csv("./data/Expanded_Area_Data_SpRL/low_dens_data_exp.csv", show_col_types = FALSE))
    sizes <- sizes_df$size
    
    load('./code/data_simulation/size_models/low_dens_size_mod.RData')
    load('./code/data_simulation/size_models/low_dens_int_rad.RData')
    
  }else if(density == "med"){
    
    a_x <- 326996
    a_y <- 4311239
    b_x <- 327096
    b_y <- 4311339
    
    sizes_df <- as.data.frame(read_csv("./data/Expanded_Area_Data_SpRL/med_dens_data_exp_2.csv", show_col_types = FALSE))
    sizes <- sizes_df$size
    
    load('./code/data_simulation/size_models/med_dens_size_mod.RData')
    load('./code/data_simulation/size_models/med_dens_int_rad.RData')
    
  }else if(density == "high"){
    
    a_x <- 327096
    a_y <- 4311239
    b_x <- 327196
    b_y <- 4311339
    
    sizes_df <- as.data.frame(read_csv("./data/Expanded_Area_Data_SpRL/high_dens_data_exp_2.csv", show_col_types = FALSE))
    sizes <- sizes_df$size
    
    load('./code/data_simulation/size_models/high_dens_size_mod.RData')
    load('./code/data_simulation/size_models/high_dens_int_rad.RData')
    
  }
  
  N_exp <- length(sizes)
  a_x_exp <- a_x - 15
  a_y_exp <- a_y - 15
  b_x_exp <- b_x + 15
  b_y_exp <- b_y + 15
  
  xy <- as.matrix(cbind(sizes_df$x, sizes_df$y))
  colnames(xy) <- c("x", "y")
  distance_matrix <- as.matrix(dist(xy, method = "euclidean"))
  
  
  ## Obtain the relevant subsets of the covariate raster data and midpoint
  mid <- matrix(c((a_x_exp+b_x_exp)/2, (a_y_exp+b_y_exp)/2), nrow = 1)
  south <- scale(crop(southness_rast, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp)))
  slope <- scale(crop(slope_rast, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp)))
  wetness <- scale(crop(wetness_rast, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp)))
  DEM <- scale(crop(DEM_rast, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp)))
  
  
  ## Set the measurement error according to the error specified
  if(noise == "small"){
    sigma2 <- (.25)^2
  }else if(noise == "medium"){
    sigma2 <- (.35)^2
  }else if(noise == "large"){
    sigma2 <- (.45)^2
  }
  
  ## Generate a translation given the variance parameter sigma2_t
  t_list <- list()
  t_list[[1]] <- matrix(data = c(0,0), nrow = 1)
  t_list[[2]] <-  matrix(data = c(.025, .025), nrow = 1)
  # for(i in 2:2){
  #   t_list[[i]] <- matrix(data = rmvnorm(1, mean = c(0,0), sigma = diag(sigma2_t, 2)), nrow = 1)
  # }
  
  ## Generate the rotation matrix given the rotation angle theta
  R_list <- list()
  R_list[[1]] <- diag(2)
  for(i in 2:2){
    R_list[[i]] <- matrix(data = c(cos(theta[i]), -sin(theta[i]), sin(theta[i]), cos(theta[i])), nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  
  ## Generate the latent marked point pattern with size based interaction
  s <- matrix(c(runif(1, a_x_exp, b_x_exp), runif(1, a_y_exp, b_y_exp), rep(NA, (N_exp - 1)*2)), ncol=2, byrow = TRUE)  # x and y corresponding to M_n(t_0)
  s_sizes <- rep(NA, N_exp)
  int_rad <- rep(NA, N_exp)
  # s_one_covars <- cbind(terra::extract(south, vect(s[1,, drop = FALSE]), method = "bilinear", ID = FALSE), 
  #                       terra::extract(slope, vect(s[1,, drop = FALSE]), method = "bilinear", ID = FALSE),
  #                       terra::extract(wetness, vect(s[1,, drop = FALSE]), method = "bilinear", ID = FALSE),
  #                       terra::extract(DEM, vect(s[1,, drop = FALSE]), method = "bilinear", ID = FALSE))
  s_one_covars <- cbind(terra::extract(south, s[1,, drop = FALSE], method = "bilinear"), 
                        terra::extract(slope, s[1,, drop = FALSE], method = "bilinear"),
                        terra::extract(wetness, s[1,, drop = FALSE], method = "bilinear"),
                        terra::extract(DEM, s[1,, drop = FALSE], method = "bilinear"))
  
  s_covars <- data.frame(Snodgrass_aspect_southness_1m = rep(NA, N_exp), Snodgrass_slope_1m = rep(NA, N_exp), Snodgrass_wetness_index_1m = rep(NA, N_exp),
                         Snodgrass_DEM_1m = rep(NA, N_exp), south.nbr = rep(NA, N_exp), slope.nbr = rep(NA, N_exp),
                         wet.nbr = rep(NA, N_exp), DEM.nbr = rep(NA, N_exp), near.nbr.dist = rep(NA, N_exp),
                         near.nbr.num = rep(NA, N_exp), avg.nbr.dist.15 = rep(NA, N_exp), near.nbr.size = rep(NA, N_exp),
                         near.nbr.size.all = rep(NA, N_exp), near.nbr.size.dist.ratio = rep(NA, N_exp),
                         x = rep(NA, N_exp), y = rep(NA, N_exp), age = rep(NA, N_exp))
  s_covars[1,] <- c(s_one_covars,
                    s_one_covars$Snodgrass_aspect_southness_1m,
                    s_one_covars$Snodgrass_slope_1m,
                    s_one_covars$Snodgrass_wetness_index_1m,
                    s_one_covars$Snodgrass_DEM_1m,
                    max(distance_matrix),
                    0,
                    max(distance_matrix),
                    min(sizes),
                    min(sizes),
                    min(sizes)/max(distance_matrix),
                    s[1,1],
                    s[1,2],
                    0)
  
  
  s_sizes[1] <- predict(size.mod, new_data = s_covars[1,])[[1]]
  int_rad [1]<- interaction.radius(int.rad.mod, s_sizes[1])
  
  
  for(i in 2:N_exp){
    
    while(is.na(s[i,1])){
      
      newp <- matrix(c(runif(1, a_x_exp, b_x_exp), runif(1, a_y_exp, b_y_exp)), nrow = 1)
      # newp_covars <- cbind(terra::extract(south, vect(newp[1,, drop = FALSE]), method = "bilinear", ID = FALSE),
      #                      terra::extract(slope, vect(newp[1,, drop = FALSE]), method = "bilinear", ID = FALSE),
      #                      terra::extract(wetness, vect(newp[1,, drop = FALSE]), method = "bilinear", ID = FALSE),
      #                      terra::extract(DEM, vect(newp[1,, drop = FALSE]), method = "bilinear", ID = FALSE))
      newp_covars <- cbind(terra::extract(south, newp[1,, drop = FALSE], method = "bilinear"),
                           terra::extract(slope, newp[1,, drop = FALSE], method = "bilinear"),
                           terra::extract(wetness, newp[1,, drop = FALSE], method = "bilinear"),
                           terra::extract(DEM, newp[1,, drop = FALSE], method = "bilinear"))
      
      dist_list <- apply(s[1:(i - 1),, drop = FALSE], 1, function(x) pq.dist(p = x, q = newp))
      close_points_15 <- unique(which(dist_list < 15))
      close_sizes_15 <- ifelse(length(close_points_15) > 0, s_sizes[close_points_15], 0)
      newp_covars$south.nbr <- sum(s_covars$Snodgrass_aspect_southness_1m[close_points_15])
      newp_covars$slope.nbr <- sum(s_covars$Snodgrass_slope_1m[close_points_15])
      newp_covars$wet.nbr <- sum(s_covars$Snodgrass_wetness_index_1m[close_points_15])
      newp_covars$DEM.nbr <- sum(s_covars$Snodgrass_DEM_1m[close_points_15])
      newp_covars$near.nbr.dist <- min(dist_list)
      newp_covars$near.nbr.num <- length(close_points_15)
      newp_covars$avg.nbr.dist.15 <- mean(dist_list[close_points_15])
      newp_covars$near.nbr.size <- s_sizes[unique(which(dist_list == newp_covars$near.nbr.dist))]
      newp_covars$near.nbr.size.all <- mean(close_sizes_15)
      if(length(close_points_15) == 0){
        newp_covars$avg.nbr.dist.15 <- min(dist_list)
        newp_covars$near.nbr.size.all <- s_sizes[unique(which(dist_list == newp_covars$near.nbr.dist))]
        newp_covars$south.nbr <- newp_covars$Snodgrass_aspect_southness_1m
        newp_covars$slope.nbr <- newp_covars$Snodgrass_slope_1m
        newp_covars$wet.nbr <- newp_covars$Snodgrass_wetness_index_1m
        newp_covars$DEM.nbr <- newp_covars$Snodgrass_DEM_1m
        newp_covars$near.nbr.size.all <- 0
      }
      newp_covars$near.nbr.size.dist.ratio <- newp_covars$near.nbr.size[[1]]/newp_covars$near.nbr.dist[[1]]
      newp_covars$x <- newp[1]
      newp_covars$y <- newp[2]
      newp_covars$age <- sizes[1] - sizes[i]
      
      newp_size <- predict(size.mod, new_data = newp_covars)[[1]]
      if(newp_size > 0){
        newp_int_rad <- interaction.radius(int.rad.mod, newp_size)
        
        which_points <- which(dist_list < int_rad[1:(i - 1)])
        
        if(length(which_points) > 0){
          accept <- as.numeric(runif(length(which_points)) < .01)
          
          if(length(which_points) == sum(accept)){
            s[i,] <- newp
            s_sizes[i] <- newp_size
            s_covars[i,] <- newp_covars
            int_rad[i] <- newp_int_rad
          }
        }else {
          s[i,] <- newp
          s_sizes[i] <- newp_size
          s_covars[i,] <- newp_covars
          int_rad[i] <- newp_int_rad
        }
      }
      
    }
    
  }
  
  ## Obtain the ordered sizes, weights for the recruit mixture distribution, and the number of recruits
  # sizes_order <- rev(sort(s_sizes))
  sizes_order <- s_sizes
  size_weights <- sizes_order/sum(sizes_order)
  n_recruits <- floor(sum(sizes_order))
  
  ## Sample the recruits from a mixture of t(1) distributions centered at the parent points
  parent_sample <- sample(seq_len(N_exp), n_recruits, prob = size_weights, replace = TRUE)
  recruits <- rmvt(n_recruits, sigma = diag(2), df = 1) + s[parent_sample,]
  recruits <- recruits[which(a_x_exp < recruits[,1] & recruits[,1] < b_x_exp & a_y_exp < recruits[,2] & recruits[,2] < b_y_exp),]
  
  ## Calculate the distances of each parent point to the recruits
  recruit_dists <- crossdist(as.ppp(recruits, c(a_x_exp, b_x_exp, a_y_exp, b_y_exp)), as.ppp(s, c(a_x_exp, b_x_exp, a_y_exp, b_y_exp)))
  
  ## Restrict parent points to 100m^2 area
  ## accounts for edge effects and possibility of recruits seeded by parents outside the domain
  selected_s <- which(a_x < s[,1] & s[,1] < b_x & a_y < s[,2] & s[,2] < b_y)
  N <- nrow(s)
  N_obs <- length(selected_s)
  
  ## Thin the recruits based on the spatial domain and proximity to parents points
  ## removing all points that fall within the interaction radii of the parent points
  ## and assign sizes to the recruits according to a Beta(1,5) distribution
  remove_points <- list()
  for(i in seq_len(N)){
    
    if(min(recruit_dists[,i]) < int_rad[i]){
      
      if(runif(1) > .1){
        
        remove_points[[i]] <- which(recruit_dists[,i] < int_rad[i])
        
      }
    }
  }
  
  remove_points <- unique(unlist(remove_points))
  selected_recruits <- recruits[-remove_points,]
  
  
  if(density == "low" & alpha == 1){
    rs_b <- 4.5
  } else if(density == "low" & alpha == 2){
    rs_b <- 1.8
  } else if(density == "low" & alpha == 3){
    rs_b <- 1.1
  } else if(density == "med" & alpha == 1){
    rs_b <- 4.25
  } else if(density == "med" & alpha == 2){
    rs_b <- 1.4
  } else if(density == "med" & alpha == 3){
    rs_b <- .95
  } else if(density == "high" & alpha == 1){
    rs_b <- 4
  } else if(density == "high" & alpha == 2){
    rs_b <- 1.2
  } else if(density == "high" & alpha == 3){
    rs_b <- .75
  }
  
  
  size_threshold <- min(sizes_order)
  # recruit_sizes <- rbeta(nrow(selected_recruits), 1, 4.5)*size_threshold
  recruit_sizes <- rbeta(nrow(selected_recruits), 1, rs_b)*size_threshold
  
  
  ## Obtain the set of observed points and potential recruits
  s <- rbind(s, selected_recruits)
  s_sizes <- cbind(s, c(sizes_order, rev(sort(recruit_sizes))))
  
  
  ## Sample the covariate values for the latents and recruits
  # scaled_covars <- cbind(rep(1, nrow(s)), terra::extract(south, vect(s), method = "bilinear", ID = FALSE),
  #                        terra::extract(slope, vect(s), method = "bilinear", ID = FALSE),
  #                        terra::extract(wetness, vect(s), method = "bilinear", ID = FALSE),
  #                        terra::extract(DEM, vect(s), method = "bilinear", ID = FALSE))
  scaled_covars <- cbind(rep(1, nrow(s)), terra::extract(south, s, method = "bilinear"),
                         terra::extract(slope, s, method = "bilinear"),
                         terra::extract(wetness, s, method = "bilinear"),
                         terra::extract(DEM, s, method = "bilinear"))
  
  
  ## Generate the observed point pattern for the first file:
  ## parent points are considered with noise added on top
  X_list <- list()
  lambda_list <- list()
  size_list <- list()
  size_list[[1]] <- sizes_order
  lambda_list[[1]] <- seq_len(N)
  X <- t(apply(s[lambda_list[[1]],], 1, function(x) rmvnorm(1, x, diag(sigma2, 2))))
  
  X_list[[1]] <- cbind(X, sizes_order)
  
  
  ## Obtain the starting values for the growth/mortality cycle
  size_list[[2]] <- as.numeric(c(sizes_order, rev(sort(recruit_sizes))))
  lambda_list[[2]] <- seq_len(nrow(s))
  
  # ## Run the annual growth/mortality algorithm for selected # of cycles
  # one_cycle <- growth.mortality.cycle(n.cycles = n.cycles, betas = betas, b = b, tau2 = tau2, gammas = gammas, phi = phi, eta = eta,
  #                                     id = lambda_list[[2]], covars = as.matrix(scaled_covars),
  #                                     size = size_list[[2]], size_threshold = size_threshold, N_obs = N_obs)
  # size_list[[3]] <- c(one_cycle[[1]])
  a <- c(as.matrix(scaled_covars)%*%as.matrix(betas))
  mu <- (a*(size_list[[2]]^alpha))/(size_list[[2]]^alpha + b^alpha)
  # W <- diag(sizes_order)
  W <- diag(length(sizes_order))
  
  # size_list[[3]] <- size_list[[2]] + c(rnorm(N, mu[1:N], sqrt(tau2)), rep(0, length(mu[-c(1:N)])))
  # size_list[[3]] <- size_list[[2]] + c(rnorm(N, mu[1:N], sqrt(tau2)), (mu[-c(1:N)]))
  size_list[[3]] <- size_list[[2]] + c(rmvnorm(n = 1, mean = c(mu)[1:N], sigma = tau2*W), (mu[-c(1:N)]))
  
  
  
  
  ## Generate the observed point pattern for the second file:
  ## all points that survive the growth/mortality cycle are considered
  ## if they exist across both files, a growth angle, and growth rate are determined based on neighboring points and size
  ## if they only exist in the second file, noise is added to the location
  ## all points in the second file are rotated and translated
  # lambda_list[[3]] <- which(one_cycle[[1]] > 0)
  lambda_list[[3]] <- which(size_list[[3]] > size_threshold)
  lambda_2_index <- seq_len(length(lambda_list[[3]]))
  s_sizes_2 <- cbind(s[lambda_list[[3]],], size_list[[3]][lambda_list[[3]]])
  X <- matrix(NA, nrow = length(lambda_list[[3]]), ncol = 2)
  # for(j in 1:length(lambda_list[[3]])){
  #   if(any(lambda_list[[1]] == lambda_list[[3]][j])){
  #     close_points <- neighbors(s_sizes_2[lambda_2_index[j],], s_sizes_2[-lambda_2_index[j],])
  #     if(any(close_points > 0)){
  #       if(which.max(close_points) == 1){
  #         angle <- runif(1, pi, (4/3)*pi)
  #       }else if(which.max(close_points) == 2){
  #         angle <- runif(1, (4/3)*pi, (5/3)*pi)
  #       }else if(which.max(close_points) == 3){
  #         angle <- runif(1, (5/3)*pi, 2*pi)
  #       }else if(which.max(close_points) == 4){
  #         angle <- runif(1, 0, pi/3)
  #       }else if(which.max(close_points) == 5){
  #         angle <- runif(1, pi/3, (2/3)*pi)
  #       }else if(which.max(close_points) == 6){
  #         angle <- runif(1, (2/3)*pi, pi)
  #       }
  #     } else{
  #       angle <- runif(1, 0, 2*pi)
  #     }
  #     size_j <- size_list[[3]][lambda_list[[3]][j]]
  #     growth_rate <- (1.25*size_j)/(600 + size_j)
  #     X_grow <- s[lambda_list[[3]][j],] + matrix(c(growth_rate*cos(angle), growth_rate*sin(angle)), nrow = 1)
  #     X_grow_trans <- (X_grow - mid)%*%t(R_list[[2]]) + t_list[[2]] + mid
  #     X[j,] <- rmvnorm(1, c(X_grow_trans), diag(sigma2,2))
  #   }else {
  #     X_trans <- (s[lambda_list[[3]][j],] - mid)%*%t(R_list[[2]]) + t_list[[2]] + mid
  #     X[j,] <- rmvnorm(1, c(X_trans), diag(sigma2,2))
  #   }
  # }
  for(j in 1:length(lambda_list[[3]])){
    X_trans <- (s[lambda_list[[3]][j],] - mid)%*%t(R_list[[2]]) + t_list[[2]] + mid
    X[j,] <- rmvnorm(1, c(X_trans), diag(sigma2,2))
  }
  
  X_2 <- cbind(X, size_list[[3]][lambda_list[[3]]])
  X_list[[2]] <- X_2[which(X_2[,3] > size_threshold),]
  
  
  ## Obtain the full simulated data and the observed data without latents
  # sim_data <- as.data.frame(cbind(rbind(s_sizes, do.call(rbind, X_list)),
  #                                 c(seq(1:(nrow(s))), c(lambda_list[[1]], lambda_list[[3]][which(X_2[,3] > size_threshold)])),
  #                                 c(rep(1, nrow(s)), rep(c(2:3), c(nrow(X_list[[1]]), nrow(X_list[[2]])))),
  #                                 c(rep(0, nrow(s)), rep(0, nrow(X_list[[1]])), one_cycle[[3]][which(X_2[,3] > size_threshold)])))
  sim_data <- as.data.frame(cbind(rbind(s_sizes, do.call(rbind, X_list)),
                                  c(seq(1:(nrow(s))), c(lambda_list[[1]], lambda_list[[3]][which(X_2[,3] > size_threshold)])),
                                  c(rep(0, nrow(s)), rep(c(1:2), c(nrow(X_list[[1]]), nrow(X_list[[2]])))),
                                  c(rep(0, nrow(s)), rep(0, nrow(X_list[[1]])), rep(0, nrow(X_list[[2]])))))
  names(sim_data) <- c("x", "y", "size", "id", "file", "mort")
  
  X <- data.frame(x = sim_data$x[(nrow(s)+1):nrow(sim_data)],
                  y = sim_data$y[(nrow(s)+1):nrow(sim_data)],
                  size = sim_data$size[(nrow(s)+1):nrow(sim_data)],
                  id = sim_data$id[(nrow(s)+1):nrow(sim_data)],
                  file = c(rep(c(1:2), c(nrow(X_list[[1]]), nrow(X_list[[2]])))),
                  mort = sim_data$mort[(nrow(s)+1):nrow(sim_data)])
  
  X_obs <- X[which(a_x < X$x & X$x < b_x & a_y < X$y & X$y < b_y),]
  
  # all_data <- rbind(sim_data[which(sim_data$id <= N_exp & sim_data$file == 0),],
  #                   X_obs)
  all_data <- rbind(sim_data[which(sim_data$id %in% X_obs$id & sim_data$file == 0),],
                    X_obs)
  names(all_data) <- c("x", "y", "size", "id", "file", "mort")
  
  
  out <- list()
  out[[1]] <- sim_data
  out[[2]] <- c(N, sigma2_t, theta, density, noise, N_obs)
  out[[3]] <- X
  out[[4]] <- s_sizes
  out[[5]] <- list(lambda_list[[1]], lambda_list[[3]][which(X_2[,3] > size_threshold)])
  out[[6]] <- X_obs
  out[[7]] <- list(lambda_list[[1]][which(a_x < X_list[[1]][,1] & X_list[[1]][,1] < b_x & a_y < X_list[[1]][,2] & X_list[[1]][,2] < b_y)],
                   lambda_list[[3]][which(X_2[,3] > size_threshold & a_x < X_2[,1] & X_2[,1] < b_x & a_y < X_2[,2] & X_2[,2] < b_y)])
  out[[8]] <- all_data
  
  return(out)
}
