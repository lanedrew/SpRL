###################################################################
#### This script generates Figure 1 from the paper supplement. ####
###################################################################

## Load libraries and sampler functions ----
library(readr)
library(data.table)
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(dplyr)

## Load C++ helper functions
sourceCpp('./resources/code/cpp_code/two_stage_func.cpp')

## Set seed for reproducibility ----
set.seed(90210)

## Set a colorblind friendly palate to use for visualizations
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Create a list for storing the number of unique latents
N_list <- list()

## Read in the empirical linkage files and calculate the number of unique latents at each iteration for each chain
for(i in 1:4){
  
  linkage_file <- paste0("./2_empirical_analysis/model_results/record_linkage_model/empirical_linkage_lambda_results_N_25_sigma2_prior_weak_chain_", i, ".csv")
  linkage_sample <- fread(file = linkage_file, header = FALSE) %>% as.matrix()
  N_list[[i]] <- apply(linkage_sample, 1, function(x) length(unique(x)))
  
}

N_length <- lapply(N_list, length)
min_N <- min(unlist(N_length))

n_latents_df <- data.frame(N = c(head(N_list[[1]], min_N), head(N_list[[2]], min_N),
                                 head(N_list[[3]], min_N), head(N_list[[4]], min_N)),
                           seq_val = rep(seq_len(min_N), 4),
                           chain = rep(c(1, 2, 3, 4), each = min_N))

## Plot and save the traceplots for estimated N
n_latents_df |> ggplot(aes(x = seq_val, y = N, group = as.factor(chain), color = as.factor(chain))) +
  geom_line(alpha = .5) +
  scale_color_manual(values = cbbPalette) +
  labs(x = "", main = "", color = "Chain") +
  theme_bw() -> n_plot

ggsave(filename = "8_N_traceplot.png", plot = n_plot, path = "./3_figures_and_tables/",
       width = 20, height = 12, units = "cm", dpi = "retina")

## Calculate Rhat value
N_mat <- as.matrix(do.call(cbind, split(n_latents_df$N, n_latents_df$chain)))
Rhat(N_mat)
