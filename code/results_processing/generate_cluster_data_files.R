#########################################################################################################
#### This script identifies the outlying observed growths resulting from the record linkage model.   ####
#########################################################################################################

## Specify an index for the linkage sample
index <- as.numeric(29)

## Load libraries and sampler functions ----
library(readr) ## load and save results
library(dplyr)
library(data.table)
library(ggplot2)
library(purrr)
library(sf)

## Set seed for reproducibility ----
set.seed(90210)

# Define the spatial domain
a_x <- 326096
a_y <- 4309939
b_x <- 328096
b_y <- 4311939
a_x2 <- a_x + 15
a_y2 <- a_y + 15
b_x2 <- b_x - 15
b_y2 <- b_y - 15


## Read in and process the specified dataset
file2015 <- read_csv("./data/UER_lidar_canopy_segmentation/crown_attributes_2015.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(XTOP, YTOP, CANVOL2015, ZTOP, fid, WKT) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2015, height = ZTOP, tree_id = fid) %>%
  mutate(file = 1)
file2019 <- read_csv("./data/UER_lidar_canopy_segmentation/crown_attributes_2019.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(XTOP, YTOP, CANVOL2019, ZTOP, fid, WKT) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2019, height = ZTOP, tree_id = fid) %>%
  mutate(file = 2)

## Restrict data to trees from the first file within 15m of the boundary for accurate competition metrics
in_bounds <- which(file2015$x > a_x2 & file2015$x < b_x2 & file2015$y > a_y2 & file2015$y < b_y2)
scan_data <- rbind(file2015, file2019)
in_bounds_full <- which((scan_data$x > a_x2 & scan_data$x < b_x2 & scan_data$y > a_y2 & scan_data$y < b_y2 & scan_data$file == 1) | (scan_data$file == 2))

## Obtain the linkage sample given the specified index
linkage_file <- paste0("./code/empirical_data/empirical_linkage_lambda_pooled_N_25.csv")
sample_index <- fread(file = "./code/empirical_data/LA_sample_index_pooled_N_25.csv", header = FALSE) %>% as.matrix()
current_index <- sample_index[index]
linkage_sample <- fread(file = linkage_file, skip = current_index, header = FALSE, nrows = 1) %>% as.matrix()
linkage_sample <- linkage_sample + 1

## Restrict the data to our bounded domain and add the linkage id
data_per_it <- scan_data[in_bounds_full,]
data_per_it$id <- linkage_sample[in_bounds_full]

## Split the data from each file for merging
data_file_1 <- data_per_it %>%
  filter(file == 1)
data_file_2 <- data_per_it %>%
  filter(file == 2)

## Merge the two files by linkage id and calculate observed growth and relative growth
linked_data <- inner_join(data_file_1, data_file_2, by = "id")
linked_data <- linked_data %>% mutate(est_growth = (size.y - size.x)/4,
                                      delta_can = (size.y - size.x)/size.x)


## Full set of linked clusters
linked_data$gc_index <- linked_data$est_growth > (-Inf)
linked_data2 <- linked_data |> filter(gc_index == "TRUE")

## Recombine the data into one dataframe spanning both files
full_data_trans <- data.frame(x = c(linked_data2$x.x, linked_data2$x.y),
                              y = c(linked_data2$y.x, linked_data2$y.y),
                              size = c(linked_data2$size.x, linked_data2$size.y),
                              WKT = c(linked_data2$WKT.x, linked_data2$WKT.y),
                              file = c(linked_data2$file.x, linked_data2$file.y),
                              id = rep(linked_data2$id, 2))

filename_full <- paste0("./data/linked_data_full.csv")
write_csv(full_data_trans, filename_full)

full_data_trans <- full_data_trans |> mutate(WKT = st_as_sfc(WKT))

## Identify and plot a subset of the linkage
filter_index <- which(full_data_trans$x > 326800 & full_data_trans$x < 327000 & full_data_trans$y > 4311500 & full_data_trans$y < 4311700)                              
full_data_trans[filter_index,] |> ggplot(aes(x = x, y = y, group = id, color = as.factor(file))) +
  geom_sf(aes(geometry = WKT, fill = as.factor(file)), alpha = .5) +
  geom_line(color = "black") +
  labs(x = "", y = "", color = "File", fill = "File", title = "Linkage for Full Data Subset")



## Clusters restricted to positive growth
linked_data$gc_index <- linked_data$est_growth > 0
linked_data2 <- linked_data |> filter(gc_index == "TRUE")

## Recombine the data into one dataframe spanning both files
pos_data_trans <- data.frame(x = c(linked_data2$x.x, linked_data2$x.y),
                             y = c(linked_data2$y.x, linked_data2$y.y),
                             size = c(linked_data2$size.x, linked_data2$size.y),
                             WKT = c(linked_data2$WKT.x, linked_data2$WKT.y),
                             file = c(linked_data2$file.x, linked_data2$file.y),
                             id = rep(linked_data2$id, 2))

filename_pos <- paste0("./data/linked_data_pos_growth.csv")
write_csv(pos_data_trans, filename_pos)

pos_data_trans <- pos_data_trans |> mutate(WKT = st_as_sfc(WKT))

## Identify and plot a subset of the linkage
filter_index <- which(pos_data_trans$x > 326800 & pos_data_trans$x < 327000 & pos_data_trans$y > 4311500 & pos_data_trans$y < 4311700)                              
pos_data_trans[filter_index,] |> ggplot(aes(x = x, y = y, group = id, color = as.factor(file))) +
  geom_sf(aes(geometry = WKT, fill = as.factor(file)), alpha = .5) +
  geom_line(color = "black") +
  labs(x = "", y = "", color = "File", fill = "File", title = "Linkage for Positive Data Subset")



## Clusters restricted to "extreme" positive growth
linked_data$gc_index <- linked_data$est_growth > 30
linked_data2 <- linked_data |> filter(gc_index == "TRUE")

## Recombine the data into one dataframe spanning both files
pos_ex_data_trans <- data.frame(x = c(linked_data2$x.x, linked_data2$x.y),
                                y = c(linked_data2$y.x, linked_data2$y.y),
                                size = c(linked_data2$size.x, linked_data2$size.y),
                                WKT = c(linked_data2$WKT.x, linked_data2$WKT.y),
                                file = c(linked_data2$file.x, linked_data2$file.y),
                                id = rep(linked_data2$id, 2))

filename_pos_ex <- paste0("./data/linked_data_pos_growth_ex.csv")
write_csv(pos_ex_data_trans, filename_pos_ex)

pos_ex_data_trans <- pos_ex_data_trans |> mutate(WKT = st_as_sfc(WKT))

## Identify and plot a subset of the linkage
filter_index <- which(pos_ex_data_trans$x > 326800 & pos_ex_data_trans$x < 327000 & pos_ex_data_trans$y > 4311500 & pos_ex_data_trans$y < 4311700)                              
pos_ex_data_trans[filter_index,] |> ggplot(aes(x = x, y = y, group = id, color = as.factor(file))) +
  geom_sf(aes(geometry = WKT, fill = as.factor(file)), alpha = .5) +
  geom_line(color = "black") +
  labs(x = "", y = "", color = "File", fill = "File", title = "Linkage for Extreme Positive Data Subset")


## Calculate number of clusters according to different growth cutoffs
points_comp <- c(linked_data |> select(id) |> unique() |> unlist() |> length(),
                 linked_data |> filter(est_growth > -10) |> select(id) |> unique() |> unlist() |> length(),
                 linked_data |> filter(est_growth > -5) |> select(id) |> unique() |> unlist() |> length(),
                 linked_data |> filter(est_growth > -1) |> select(id) |> unique() |> unlist() |> length(),
                 linked_data |> filter(est_growth > 0) |> select(id) |> unique() |> unlist() |> length(),
                 linked_data |> filter(est_growth > 20) |> select(id) |> unique() |> unlist() |> length(),
                 linked_data |> filter(est_growth > 30) |> select(id) |> unique() |> unlist() |> length())

## Print the number of points and the relative proportions for the different cut points
points_comp
(points_comp/points_comp[1])


##########################################################
#### If working with generated data files, see below! ####
##########################################################

## Full dataset
filename_full <- paste0("./data/linked_data_full.csv")
full_data_trans <- read_csv(filename_full)
full_data_trans <- full_data_trans |> mutate(WKT = st_as_sfc(WKT))

## Identify and plot a subset of the linkage
filter_index <- which(full_data_trans$x > 326800 & full_data_trans$x < 327000 & full_data_trans$y > 4311500 & full_data_trans$y < 4311700)                              
full_data_trans[filter_index,] |> ggplot(aes(x = x, y = y, group = id, color = as.factor(file))) +
  geom_sf(aes(geometry = WKT, fill = as.factor(file)), alpha = .5) +
  geom_line(color = "black") +
  labs(x = "", y = "", color = "File", fill = "File", title = "Linkage for Full Data Subset")


## Clusters restricted to positive growth
filename_pos <- paste0("./data/linked_data_pos_growth.csv")
pos_data_trans <- read_csv(filename_pos)
pos_data_trans <- pos_data_trans |> mutate(WKT = st_as_sfc(WKT))

## Identify and plot a subset of the linkage
filter_index <- which(pos_data_trans$x > 326800 & pos_data_trans$x < 327000 & pos_data_trans$y > 4311500 & pos_data_trans$y < 4311700)                              
pos_data_trans[filter_index,] |> ggplot(aes(x = x, y = y, group = id, color = as.factor(file))) +
  geom_sf(aes(geometry = WKT, fill = as.factor(file)), alpha = .5) +
  geom_line(color = "black") +
  labs(x = "", y = "", color = "File", fill = "File", title = "Linkage for Positive Data Subset")


## Clusters restricted to "extreme" positive growth
filename_pos_ex <- paste0("./data/linked_data_pos_growth_ex.csv")
pos_ex_data_trans <- read_csv(filename_pos_ex)
pos_ex_data_trans <- pos_ex_data_trans |> mutate(WKT = st_as_sfc(WKT))

## Identify and plot a subset of the linkage
filter_index <- which(pos_ex_data_trans$x > 326800 & pos_ex_data_trans$x < 327000 & pos_ex_data_trans$y > 4311500 & pos_ex_data_trans$y < 4311700)                              
pos_ex_data_trans[filter_index,] |> ggplot(aes(x = x, y = y, group = id, color = as.factor(file))) +
  geom_sf(aes(geometry = WKT, fill = as.factor(file)), alpha = .5) +
  geom_line(color = "black") +
  labs(x = "", y = "", color = "File", fill = "File", title = "Linkage for Extreme Positive Data Subset")
