#######################################################################################
#### This script generates the competition metrics for the 2015 empirical dataset. ####
#######################################################################################

## Load libraries ----
library(readr) ## load and save results
library(Rcpp)
library(RcppDist)
library(RcppArmadillo)
library(codetools)
library(dplyr)

## Source necessary C++ helper functions
sourceCpp('./resources/code/cpp_code/two_stage_func.cpp')

## Read in the 2015 data and specify the restricted bounds
tree_data_2015 <- as_tibble(read.csv("./resources/empirical_data/UER_lidar_canopy_segmentation/crown_attributes_2015.csv")) %>%
  filter(LCmajority == 1)
restricted_bounds <- c(326096, 328096, 4309939, 4311939) + c(15, -15, 15, -15)

## Subset the data and calculate the competition metrics for the remaining points
td_2015_con <- tree_data_2015 %>%
  select(XTOP, YTOP, CANVOL2015) 
td_2015_comp_mets <- td_2015_con %>%
  as.matrix() %>%
  calc_comp_indices_arma2() %>%
  as.data.frame()

## Normalize the subsetted variables across the whole dataset
td_2015_comp_mets_norm <- td_2015_comp_mets[which(tree_data_2015$XTOP > restricted_bounds[1] &
                                                    tree_data_2015$XTOP < restricted_bounds[2] &
                                                    tree_data_2015$YTOP > restricted_bounds[3] & 
                                                    tree_data_2015$YTOP < restricted_bounds[4]),] %>%
  transmute(LNV_norm = (V1 - mean(V1))/sd(V1),
            RSI_norm = (V2 - mean(V2))/sd(V2),
            ND_norm = (V3 - mean(V3))/sd(V3))


td_2015_all_vals <- cbind(td_2015_con[which(tree_data_2015$XTOP > restricted_bounds[1] &
                                            tree_data_2015$XTOP < restricted_bounds[2] &
                                            tree_data_2015$YTOP > restricted_bounds[3] & 
                                            tree_data_2015$YTOP < restricted_bounds[4]),],
                          td_2015_comp_mets_norm)

## Save the generated competition metrics
write_csv(td_2015_all_vals, "./2_empirical_analysis/competition_metrics_2015.csv")


## Read in the 2019 data and specify the restricted bounds
tree_data_2019 <- as_tibble(read.csv("./resources/empirical_data/UER_lidar_canopy_segmentation/crown_attributes_2019.csv")) %>%
  filter(LCmajority == 1)
restricted_bounds <- c(326096, 328096, 4309939, 4311939) + c(15, -15, 15, -15)

## Subset the data and calculate the competition metrics for the remaining points
td_2019_con <- tree_data_2019 %>%
  select(XTOP, YTOP, CANVOL2019) 
td_2019_comp_mets <- td_2019_con %>%
  as.matrix() %>%
  calc_comp_indices_arma2() %>%
  as.data.frame()

## Normalize the subsetted variables across the whole dataset
td_2019_comp_mets_norm <- td_2019_comp_mets[which(tree_data_2019$XTOP > restricted_bounds[1] &
                                                  tree_data_2019$XTOP < restricted_bounds[2] &
                                                  tree_data_2019$YTOP > restricted_bounds[3] & 
                                                  tree_data_2019$YTOP < restricted_bounds[4]),] %>%
  transmute(LNV_norm = (V1 - mean(V1))/sd(V1),
            RSI_norm = (V2 - mean(V2))/sd(V2),
            ND_norm = (V3 - mean(V3))/sd(V3))


td_2019_all_vals <- cbind(td_2019_con[which(tree_data_2019$XTOP > restricted_bounds[1] &
                                            tree_data_2019$XTOP < restricted_bounds[2] &
                                            tree_data_2019$YTOP > restricted_bounds[3] & 
                                            tree_data_2019$YTOP < restricted_bounds[4]),],
                          td_2019_comp_mets_norm)

## Save the generated competition metrics
write_csv(td_2019_all_vals, "./2_empirical_analysis/competition_metrics_2019.csv")
