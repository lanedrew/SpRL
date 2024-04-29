#################################################################################
#### This script generates a simulated dataset using the supplied arguments. ####
#################################################################################

## Arguments passed from command line ----
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) stop("Pass in density (low or med or high), tau2 value, alpha value, noise level, and index (integer from 1 to 100)", call.=FALSE)
if (!(args[1] %in% c("low", "med", "high"))) stop("Pass in the density (low or med or high)", call.=FALSE)
if (!(args[4] %in% c("small", "medium", "large"))) stop("Pass in the noise level (small or medium or large)", call.=FALSE)

density <- args[1]
tau2 <- as.numeric(args[2])
alpha <- as.numeric(args[3])
noise <- args[4]
index <- args[5]

print(paste0("Density: ", density, ", Tau2: ", tau2, ", Alpha: ", alpha, ", Index: ", index))

## Source the data simulation function and supporting code
source('./resources/code/R_code/simulate_data_size_pred.R')

## Specify the transformation parameters
sigma2_t <- (.1)^2
theta <- c(0, .005)

## Specify the fixed growth model parameters
betas <- c(3, .5, -.5, .5, -.5)
b <- c(12)

## Simulate a dataset with the given arguments
data_set <- sim_data(sigma2_t = sigma2_t, theta = theta, density = density, noise = noise,
                     tau2 = tau2, betas = betas, b = b, alpha = alpha)
  
## Specify the file name depending on the given arguments
filename_data_set <- paste0("./1_simulation_study/simulated_data/", density, "_dens_", noise, "_noise_", alpha, "_alpha_sim_", index,".csv")

## Write the file to the simulated_data folder in the 1_simulation_study folder in the repository
write_csv(data_set[[8]], file = filename_data_set)
  