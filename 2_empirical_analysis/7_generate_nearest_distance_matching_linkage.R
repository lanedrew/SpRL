################################################################################################
#### This script generates the Nearest Distance Matching linkage for the empirical dataset. ####
################################################################################################

## Load libraries and sampler functions ----
library(readr) 
library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(dplyr)
library(data.table)

## Load C++ helper functions
sourceCpp('./resources/code/cpp_code/two_stage_func.cpp')

## Set seed for reproducibility ----
set.seed(90210)

# Define the spatial domain
a_x <- 326096
a_y <- 4309939
b_x <- 328096
b_y <- 4311939


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

file2015_mat <- file2015 |> select(x, y) |> as.matrix()
file2019_mat <- file2019 |> select(x, y) |> as.matrix()

## Generate and save the NDM linkage
lambda <- nearest_distance_matching(file_1 = file2015_mat, file_2 = file2019_mat, dist = 2)
write_csv(as.data.frame(lambda), "./2_empirical_analysis/model_results/record_linkage_model/near_distance_matching_lambda.csv")
