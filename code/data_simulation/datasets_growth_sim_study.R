## Source the data simulation functions and relevant functions that it uses
source('./code/data_simulation/simulate_data_size_pred.R')

## Specify the transformation parameters
sigma2_t <- (.1)^2
theta <- c(0, .005)
tau2 <- .5
b <- 12
betas <- c(3, .5, -.5, .5, -.5)
alpha <- 3
noise <- "large"


for(i in 1:100){
  ## Simulate a low, medium, and high density data set on each iteration
  low_dens_sim <- sim_data(sigma2_t = sigma2_t, theta = theta, density = "low", noise = noise, n.cycles = 1,
                           tau2 = tau2, betas = betas, b = b, alpha = alpha)

  med_dens_sim <- sim_data(sigma2_t = sigma2_t, theta = theta, density = "med", noise = noise, n.cycles = 1,
                           tau2 = tau2, betas = betas, b = b, alpha = alpha)

  high_dens_sim <- sim_data(sigma2_t = sigma2_t, theta = theta, density = "high", noise = noise, n.cycles = 1,
                            tau2 = tau2, betas = betas, b = b, alpha = alpha)

  # ## Give the files sequential names by iteration for easy reading later
  # filename_low_dens <- paste0("./code/growth_sim_data_3/low_dens_large_noise_sim_", i,".csv")
  # filename_med_dens <- paste0("./code/growth_sim_data_3/med_dens_large_noise_sim_", i,".csv")
  # filename_high_dens <- paste0("./code/growth_sim_data_3/high_dens_large_noise_sim_", i,".csv")
  
  ## Give the files sequential names by iteration for easy reading later
  filename_low_dens <- paste0("./code/growth_sim_data_3/low_dens_", noise, "_noise_alpha_", alpha, "_sim_", i,".csv")
  filename_med_dens <- paste0("./code/growth_sim_data_3/med_dens_", noise, "_noise_alpha_", alpha, "_sim_", i,".csv")
  filename_high_dens <- paste0("./code/growth_sim_data_3/high_dens_", noise, "_noise_alpha_", alpha, "_sim_", i,".csv")

  ## Write the files to the growth_sim_data folder in the repository
  write_csv(low_dens_sim[[8]], file = filename_low_dens)
  write_csv(med_dens_sim[[8]], file = filename_med_dens)
  write_csv(high_dens_sim[[8]], file = filename_high_dens)
  
  print(paste0("Iteration ", i, " Completed"))
}

