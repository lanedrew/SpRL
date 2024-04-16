###################################################################################
#### This script defines the sampler functions for the two-stage model sampler ####
###################################################################################

## Load relevant libraries
library(dplyr) ## data manip
library(tmvtnorm) ## for generating tmvn random variables
library(spatstat) ## for generating spatial point processes
library(codetools)
library(tidyr) ## data manip
library(terra) ## for manipulating raster data
library(purrr) ## Advanced functionals
library(TruncExpFam)


## Define the sampler function for sigma^2
sample_sigma2 <- function(s, lambda, c_sigma, d_sigma, Y, theta, t, mu_D, m) {
  
  n <- nrow(Y)
  
  ## Generate the rotated/translated/centered s_lambda_ij associated with each record
  s_list <- list()
  s_lambda <- s[lambda,]
  m_new <- c(0, m)
  
  for(i in 1:length(m)){
    s_lambda_m <- s_lambda[c((m_new[i]+1):m_new[i+1]),, drop = FALSE]
    s_list[[i]] <- transform_s(s = s_lambda_m, theta = theta[i], t = t[i,], mu_D = mu_D)
  }
  
  s_j <- do.call(rbind, s_list)
  # sigma2 <- 1/rgamma(1, shape = n + c_sigma, rate = (.5 * sum((Y - s_j)^2) + d_sigma))
  sigma2 <- rtruncinvgamma(1, shape = n + c_sigma, rate = (.5 * sum((Y - s_j)^2) + d_sigma), a = 0, b = 1)
  
  return(sigma2)
  
}


## Define a function to perform transformations on the latents
transform_s <- function(s, theta, t, mu_D){
  
  # Generate the rotation matrix
  R_theta <- matrix(data = c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, ncol = 2, byrow = TRUE)
  
  # Generate the rotated/translated s values
  s_trans <- apply(s, 1, function(x) (matrix(data = c(x), nrow = 1) - mu_D)%*%t(R_theta) + t + mu_D)
  
  # Return the transformed values as a matrix
  return(t(s_trans))
  
}


## Define a sampler for s
sample_s <- function(s, lambda, sigma2, Y, N, theta, t, mu_D, m, a_x, a_y, b_x, b_y, file) {
  
  ## Generate the rotation matrix R for each file
  R_list <- list()
  for(i in 1:length(theta)){
    R_list[[i]] <- matrix(data = c(cos(theta[i]), -sin(theta[i]), sin(theta[i]), cos(theta[i])), nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  ## Generate the rotated/translated/centered Y_ij for each file
  Y_list <- list()
  Y_list[[1]] <- Y[which(file == 1),]
  n <- nrow(Y)
  m_new <- c(0,m,n)
  for(i in 2:length(m)){
    Y_m <- Y[c((m_new[i]+1):m_new[i+1]),, drop = FALSE]
    Y_list[[i]] <- transform_Y(Y = Y_m, theta = theta[i], t = t[i,], mu_D = mu_D)
  }
  Y_t <- do.call(rbind, Y_list)
  
  C_lambda <- vector("list", N)
  for(i in seq_len(N)){
    C_lambda[[i]] <- which(lambda == i)
  }
  
  Y_split <- map(C_lambda, ~ matrix(Y_t[.x,], ncol = 2))
  Y_means <- map(Y_split, colMeans)
  Y_means_size <- map(Y_split, nrow)
  
  
  s_star_list <- list()
  for(i in 1:length(Y_split)){
    s_star_list[[i]] <- propose_s_star(s_j_mean = Y_means[[i]], n_j = Y_means_size[[i]],
                                       sigma2 = sigma2, a_x = a_x, a_y = a_y,
                                       b_x = b_x, b_y = b_y)
  }
  s_star <- do.call(rbind, s_star_list)
  
  return(s_star)
  
}

## Define a function to transform Y
transform_Y <- function(Y, theta, t, mu_D){
  
  # Generate the rotation matrix
  R_theta <- matrix(data = c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, ncol = 2, byrow = TRUE)
  
  # Generate the rotated/translated s values
  Y_trans <- t(apply(Y, 1, function(x) (x - t - mu_D)%*%R_theta + mu_D))
  
  return(Y_trans)
  
}

## Define a function to propose s_j*
propose_s_star <- function(s_j_mean, n_j, sigma2, a_x, a_y, b_x, b_y){
  
  n_j <- as.numeric(n_j)
  if(n_j > 0){
    s_star <- mvtnorm::rmvnorm(n = 1, mean = s_j_mean, sigma = diag(sigma2/n_j, 2))
    
    while(s_star[1] < a_x | s_star[1] > b_x | s_star[2] < a_y | s_star[2] > b_y){
      s_star <- mvtnorm::rmvnorm(n = 1, mean = s_j_mean, sigma = diag(sigma2/n_j, 2))
    }
    
  } else {
    s_star <- cbind(runif(1, min = a_x, max = b_x), runif(1, min = a_y, max = b_y))
  }
  
  s_star
  
}


## Updated Sampler for Lambda
## Needs new inputs t - translations, theta - vector of rotation coefficients, mu_D rotation reference point, m - file sizes
sample_lambda <- function(s, sigma2, theta, t, Y, N, m, mu_D, file) {
  
  ## m is a vector of sizes for each file
  n <- nrow(Y)
  lambda <- rep(NA, n)
  
  ## Generate the rotation matrix R for each file
  R_list <- list()
  for(i in 1:length(theta)){
    R_list[[i]] <- matrix(data = c(cos(theta[i]), -sin(theta[i]), sin(theta[i]), cos(theta[i])), nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  s_list <- list()
  s_list[[1]] <- s
  for(i in 2:length(m)){
    s_list[[i]] <- transform_s(s = s, theta = theta[i], t = t[i,], mu_D = mu_D)
  }
  
  
  for(k in 1:n){
    
    sum_sq_mat <- apply(s_list[[file[k]]], 1, function(x) sum((Y[k,] - x)^2))
    logprobs <- -0.5/sigma2 * sum_sq_mat
    
    ## Apply the Gumbel-max trick
    g <- -log(-log(runif(N, 0, 1)))
    
    lambda[k] <- which.max(logprobs + g)
  
  }
  
  return(lambda)
  
}


## Define the sampler function for theta
sample_theta <- function(lambda, sigma2, Y, theta, s, t, mu_D, m, kappa, nu, sigma2_theta){
  
  n <- nrow(Y)
  
  ## Generate proposals for each theta_i for i >= 2
  theta_star <- map(theta, ~ rnorm(1, mean = .x, sd = sqrt(sigma2_theta)))
  
  ## Obtain the S matrix for each theta
  S_matrix <- S_mat(Y = Y, s = s, t = t, lambda = lambda, sigma2 = sigma2, mu_D = mu_D, m = m)
  
  ## Calculate the acceptance ratio for each proposal
  accept_reject <- c()
  accept_reject[1] <- 0
  
  for(i in 2:length(theta)){
    alpha <- exp(logprob(theta = theta[i], theta_star = theta_star[[i]], S_mat = S_matrix[[i]], kappa = kappa, nu = nu))
    
    if(runif(1) < min(1, alpha)) {
      if(abs(theta_star[[i]]) < .01){
        theta[i] <- theta_star[[i]]
        accept_reject[i] <- 1
      }
    }
    
  }
  
  return(list(theta = theta, accept_reject = accept_reject))
  
}


## Define a function to calculate the S matrix for the full conditional of theta  
S_mat <- function(Y, t, sigma2, lambda, s, mu_D, m){
  
  n <- nrow(Y)
  S_mat_list <- list()
  s_lambda <- s[lambda,]
  m.new <- c(0,m,n)
  
  for(i in 2:length(m)){
    s_lambda_m <- s_lambda[c((m.new[i]+1):m.new[i+1]),, drop = FALSE]
    Y_m <- Y[c((m.new[i]+1):m.new[i+1]),, drop = FALSE]
    S_ij_list <- pmap(list(s_x = s_lambda_m[,1], s_y = s_lambda_m[,2], y_x = Y_m[,1], y_y = Y_m[,2]),
                      function(s_x, s_y, y_x, y_y) t(t(matrix(data = c(s_x, s_y), nrow = 1) - mu_D)%*%(matrix(data = c(y_x, y_y), nrow = 1) - t[i,] - mu_D)))
    S_mat_list[[i]] <- (1/(2*sigma2))*Reduce("+", S_ij_list)
  }
  
  return(S_mat_list)
  
}

## Define the log-probability function for calculating the MH-ratio  
logprob <- function(theta, theta_star, S_mat, kappa, nu) {
  
  ## log full conditional for theta
  log_r <- (kappa*cos(nu) + S_mat[1,1] + S_mat[2,2])*(cos(theta_star) - cos(theta)) + (kappa*sin(nu) - S_mat[1,2] + S_mat[2,1])*(sin(theta_star) - sin(theta))
  
  return(log_r)
  
}


## Define the sampler function for t
sample_t <- function(lambda, sigma2, Y, s, theta, sigma2_t, mu_D, m){
  
  n <- nrow(Y)
  
  ## Generate the rotation matrix R for each file
  R_list <- list()
  
  for(i in 2:length(theta)){
    R_list[[i]] <- matrix(data = c(cos(theta[i]), -sin(theta[i]), sin(theta[i]), cos(theta[i])), nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  ## Generate the rotated/translated/centered s_lambda_ij for each file
  s_list <- list()
  s_lambda <- s[lambda,]
  m_new <- c(0,m,n)
  
  for(i in 2:length(m)){
    s_lambda_m <- s_lambda[c((m_new[i]+1):m_new[i+1]),]
    s_list[[i]] <- t(apply(s_lambda_m, 1, function(x) (matrix(data = c(x), nrow = 1) - mu_D)%*%t(R_list[[i]]) + mu_D))
  }
  
  ## Sample t_i for i = 2,...,m
  t <- matrix(NA, nrow = length(m), ncol = 2)
  t[1,] <- matrix(c(0,0), nrow = 1, ncol = 2)
  
  for(i in 2:length(m)){
    mu_t_i <- (sigma2_t/((m_new[i+1] - m_new[i])*sigma2_t + sigma2))*colSums(Y[c((m_new[i]+1):m_new[i+1]),] - s_list[[i]])
    sigma2_t_i <- (sigma2_t*sigma2)/((m_new[i+1] - m_new[i])*sigma2_t + sigma2)
    t[i,] <- mvtnorm::rmvnorm(n = 1, mean = mu_t_i, sigma = diag(sigma2_t_i,2))
  }
  
  return(t)
  
}


## Define a function to adaptively tune the M-H proposal variance
update_sigma2_tune <- function(sigma2_tune, accept_reject, i, min_val, rate){
  
  rate <- mean(accept_reject[c((i-49):i)])
  delta_n <- min(min_val, 1/(i^rate))
  
  if(rate < .44){
    sigma2_tune <- sigma2_tune - delta_n
    
    if(sigma2_tune <= 0){
      sigma2_tune <- delta_n
    }
    
  }else{
    sigma2_tune <- sigma2_tune + delta_n
  }
  
  return(sigma2_tune)
  
}


## Function to obtain linkage results
get_prec_rec <- function(lambda, true_id) {
  n <- length(lambda)
  seq_n <- seq_len(n)
  rec <- data.frame(rec_id = seq_n, lambda = lambda, true_id = true_id)
  true_links <- est_links <- matrix(NA, nrow = n, ncol = n)
  
  for(i in seq_n) {
    for(j in seq_n) {
      r_i <- rec[rec$rec_id == i, ]
      r_j <- rec[rec$rec_id == j, ]
      true_links[i, j] <- r_i$true_id == r_j$true_id
      est_links[i, j] <- r_i$lambda == r_j$lambda
    }
  }
  
  tp <- sum(true_links == est_links & est_links)
  tp_fn <- sum(true_links)
  tp_fp <- sum(est_links)
  
  c(precision = tp/tp_fp, recall = tp/tp_fn, F1_score = (2 / ( (tp/tp_fp)^(-1) + (tp/tp_fn)^(-1))))
}
