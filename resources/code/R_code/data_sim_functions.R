################################################################################################
#### This script contains the support functions for the function to generate simulated data ####
#### contained in the script simulate_data_size_pred.R.                                      ####
################################################################################################

## Load necessary packages
library(tidyr)
library(Matrix)
library(ggplot2)

### Define all supporting functions used in the data simulation function ####

## Function to find distance between points
pq.dist <- function(p, q){
  sqrt((q[1] - p[1])^2 + (q[2] - p[2])^2)
}


## Function to find the angle between two points
pq.angle <- function(p, q){
  if(p[1] == q[1] & p[2] == q[2]){
    0
  } else{
    asin((q[2] - p[2])/pq.dist(p, q))
  }
}


## Function to determine the nearest neighbors and their size by sextant
## around a given point in a radius of 3 meters
neighbors <- function(p, point.mat){
  sext.list <- list()
  for(i in seq_len(nrow(point.mat))){
    dist <- pq.dist(p[-3], point.mat[i,][-3])
    if(dist < 3 & point.mat[i,][3] > p[3]){
      theta <- pq.angle(p[-3], point.mat[i,][-3])
      if(abs(theta) == theta){
        if(0 < theta & theta < pi/3){
          sext.list[[i]] <- c(1, point.mat[i,][3]) 
        } else if(pi/3 < theta & theta < 2*pi/3){
          sext.list[[i]] <- c(2, point.mat[i,][3]) 
        } else if(2*pi/3 < theta & theta < pi){
          sext.list[[i]] <- c(3, point.mat[i,][3])
        }
      } else{
        if(0 < abs(theta) & abs(theta) < pi/3){
          sext.list[[i]] <- c(6, point.mat[i,][3]) 
        } else if(pi/3 < abs(theta) & abs(theta) < 2*pi/3){
          sext.list[[i]] <- c(5, point.mat[i,][3]) 
        } else if(2*pi/3 < abs(theta) & abs(theta) < pi){
          sext.list[[i]] <- c(4, point.mat[i,][3])
        }
      }
    }
  }
  
  sextant <- do.call(rbind, sext.list)
  sext.sizes <- list()
  for(j in seq_len(6)){
    if (any(sextant[,1] == j)){
      sext.sizes[[j]] <- sum(sextant[,2][which(sextant[,1] == j)])
    } else{
      sext.sizes[[j]] <- 0
    }
  }
  
  return(unlist(sext.sizes))
}


## Function to determine the interaction radius of points as a function
## of base rate, decay rate, and index
interaction.radius <- function(model, size){
  predict(model, newdata = data.frame(CANVOL2015 = as.numeric(size)))
}


## Function to calculate distance from recruits to parent points
recruit.parent.dists <- function(i, j, recruits, parents) {
  sqrt(Matrix::rowSums((recruits[i,, drop = FALSE] - parents[j,, drop = FALSE])^2))
}


growth <- function(betas, b, tau2, covars, size, size_threshold, N_obs){
  
  a <- c(covars%*%as.matrix(betas))
  mu <- (a*size)/(size + b)
  
  parent_sizes <- size[1:N_obs]
  parents_growth <- rnorm(length(parent_sizes), mean = mu[1:N_obs], sd = sqrt(tau2))
  recruits_growth <- mu[-c(1:N_obs)]
  growth_update <- c(parents_growth, recruits_growth)
  
  return(growth_update)
}


## Function to estimate annual mortality for a set of points given size and covariates
mortality <- function(n.cycles, phi, eta, gammas, covars, size.t1, size.t2){
  
  alpha <- exp(covars%*%as.matrix(gammas))
  p <- 1 - exp((-n.cycles)*phi + (alpha/eta)*(exp((-eta)*size.t2) - exp((-eta)*size.t1)))
  
  if(length(which(p < 0)) > 0 | length(which(p > 1)) > 0){
    p[which(p < 0)] <- 0
    p[which(p > 1)] <- 1
  }
  
  mort <- rbinom(nrow(covars), 1, p)
  
  return(mort)
}


## Function to estimate yearly cycle of growth and mortality
growth.mortality.cycle <- function(n.cycles, betas, b, tau2, phi, eta, gammas, id, covars, size, size_threshold, N_obs){
  
  ## Estimate the growth for points that survived
  size_year <- list()
  observed_growth <- list()
  size_year[[1]] <- size
  for(i in seq_len(n.cycles)){
    observed_growth[[i]] <- growth(betas = betas, b = b, tau2 = tau2, covars = covars,
                                   size = size_year[[i]], size_threshold = size_threshold,
                                   N_obs = N_obs)
    size_year[[i + 1]] <- size_year[[i]] + observed_growth[[i]]
  }
  
  size_t1 <- size_year[[1]]
  size_t2 <- size_year[[n.cycles + 1]]
  est_growth <- size_t2 - size_t1
  
  ## Estimate the mortality for a given set of points
  est_mort <- mortality(n.cycles = n.cycles, phi = phi, eta = eta, gammas = gammas, covars = covars, size.t1 = size_t1, size.t2 = size_t2)
  
  return(list(size_t2, est_growth, est_mort))
}


## Define a function to plot the simulated data with latents
plot_data <- function(plot.data, overlap){
  ggplot(data = plot.data[[1]], aes(x = V1, y = V2, group = V4)) +
    geom_point(aes(color = as.factor(V5), size = sizes_order), alpha = 0.5) +
    geom_line() +
    scale_color_discrete(name = "Point Pattern", labels = c("Latents", "X1", "X2")) +
    labs(x = "X", y = "Y", size = "Canopy Volume",
         title = paste0("Simulated Data for N=",plot.data[[2]][1],
                        ", sigma2_t=",plot.data[[2]][2],
                        ", theta=",as.numeric(plot.data[[2]][4]),
                        ", density=",plot.data[[2]][5],
                        ", noise=",plot.data[[2]][6],
                        ", overlap %=", overlap))
}


## Define a function to plot the full simulated data
plot_data_full <- function(plot.data, overlap){
  ggplot(data = plot.data[[3]], aes(x = x, y = y, group = id)) +
    geom_point(aes(color = as.factor(file), size = size), alpha = 0.5) +
    geom_line() +
    scale_color_discrete(name = "File", labels = c("X1", "X2")) +
    labs(x = "X", y = "Y", size = "Canopy Volume",
         title = paste0("Simulated Data for N=",plot.data[[2]][1],
                        ", sigma2_t=",plot.data[[2]][2],
                        ", theta=",as.numeric(plot.data[[2]][4]),
                        ", density=",plot.data[[2]][5],
                        ", noise=",plot.data[[2]][6],
                        ", overlap %=", overlap))
}


## Define a function to plot the observed simulated data
plot_data_obs <- function(plot.data, overlap){
  ggplot(data = plot.data[[6]], aes(x = x, y = y, group = id)) +
    geom_point(aes(color = as.factor(file), size = size), alpha = 0.5) +
    geom_line() +
    scale_color_discrete(name = "File", labels = c("X1", "X2")) +
    labs(x = "X", y = "Y", size = "Canopy Volume",
         title = paste0("Simulated Data for N=",plot.data[[2]][7],
                        ", sigma2_t=",plot.data[[2]][2],
                        ", theta=",as.numeric(plot.data[[2]][4]),
                        ", density=",plot.data[[2]][5],
                        ", noise=",plot.data[[2]][6],
                        ", overlap %=", overlap))
}
