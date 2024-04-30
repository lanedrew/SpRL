###################################################################################################################
#### This script generates a pooled sampling index for the record linkage model chains for the empirical data. ####
###################################################################################################################

## Load libraries and sampler functions ----
library(readr) ## load and save results
library(data.table)
library(dplyr)


## Set seed for reproducibility ----
set.seed(90210)

linkage_file <- paste0("./2_empirical_analysis/empirical_linkage_lambda_pooled_N_25.csv")
linkage_sample <- fread(file = linkage_file, header = TRUE) %>% as.matrix()

sample_index <- sample(1:nrow(linkage_sample), 100, replace = FALSE)
write_csv(as.data.frame(sample_index), "./2_empirical_analysis/empirical_data/LA_sample_index_pooled_N_25.csv", col_names = FALSE)