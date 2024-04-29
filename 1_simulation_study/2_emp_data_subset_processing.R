########################################################################################################
#### This script processes the raw 2015 data to identify the low, medium, and high density subsets. ####
########################################################################################################

## Load the necessary packages
library(tidyr) ## data manip
library(readr) ## load and save results
library(spatstat)
library(dplyr)


## Read in and filter the data
tree_data_2015 <- read_csv("./resources/empirical_data/UER_lidar_canopy_segmentation/crown_attributes_2015.csv")
td_2015_con <- tree_data_2015 %>%
  filter(LCmajority == 1)

## Create a point process object and obtain 100m^2 subets
pp15con <- ppp(td_2015_con$XTOP, td_2015_con$YTOP,c(326096, 328096), c(4309939, 4311939))
pp.15.grid <- quadrats(as.owin(pp15con), nx = 20, ny = 20)
pp.list.15 <- list()
for(i in 1:20){
  pp.list.i <- list()
  for(j in 1:20){
    pp.data <- tree_data_2015 %>%
      filter(LCmajority == 1 & YTOP > pp.15.grid$ygrid[j] & YTOP < pp.15.grid$ygrid[j+1] & XTOP > pp.15.grid$xgrid[i] & XTOP < pp.15.grid$xgrid[i+1])
    pp.list.i[[j]] <- ppp(pp.data$XTOP, pp.data$YTOP,
                          c(pp.15.grid$xgrid[i], pp.15.grid$xgrid[i+1]), c(pp.15.grid$ygrid[j], pp.15.grid$ygrid[j+1]),
                          marks = pp.data$CANVOL2015)
  }
  pp.list.15[[i]] <- pp.list.i
}
all.pat.15 <- unlist(pp.list.15, recursive = FALSE)


## Obtain and save the low density expanded area
low.win.exp <- c((all.pat.15[[96]]$window$xrange[1] - 10),
                 (all.pat.15[[96]]$window$xrange[2] + 10),
                 (all.pat.15[[96]]$window$yrange[1] - 10),
                 (all.pat.15[[96]]$window$yrange[2] + 10))
low.dens.exp.data <- tree_data_2015 %>%
  dplyr::filter(LCmajority == 1 & YTOP > low.win.exp[3] & YTOP < low.win.exp[4] & XTOP > low.win.exp[1] & XTOP < low.win.exp[2])
low.dens.exp.pp <- ppp(low.dens.exp.data$XTOP, low.dens.exp.data$YTOP,
                       c(low.win.exp[1], low.win.exp[2]), c(low.win.exp[3], low.win.exp[4]),
                       marks = low.dens.exp.data$CANVOL2015)
low.times <- list()
low.sorted <- rev(sort(low.dens.exp.pp$marks))
for(i in 1:(low.dens.exp.pp$n - 1)){
  low.times[[i]] <- low.sorted[1] - low.sorted[i+1]
}
low.times <- c(0, unlist(low.times))
low.dens <- data.frame(x = low.dens.exp.pp$x, y = low.dens.exp.pp$y, size = low.dens.exp.pp$marks)
low.dens <- low.dens[order(low.dens$size, decreasing = TRUE),]
low.dens <- data.frame(x = low.dens$x, y = low.dens$y, time = low.times, size = low.dens$size)
results_file <- paste0("./1_simulation_study/size_models/low_dens_data_exp.csv")
write_csv(low.dens, file = results_file)


## Obtain and save the medium density expanded area
med.win.exp <- c((all.pat.15[[393]]$window$xrange[1] - 20),
                 (all.pat.15[[393]]$window$xrange[2]),
                 (all.pat.15[[393]]$window$yrange[1] - 10),
                 (all.pat.15[[393]]$window$yrange[2] + 10))
med.dens.exp.data <- tree_data_2015 %>%
  dplyr::filter(LCmajority == 1 & YTOP > med.win.exp[3] & YTOP < med.win.exp[4] & XTOP > med.win.exp[1] & XTOP < med.win.exp[2])
med.dens.exp.pp <- ppp(med.dens.exp.data$XTOP, med.dens.exp.data$YTOP,
                       c(med.win.exp[1], med.win.exp[2]), c(med.win.exp[3], med.win.exp[4]),
                       marks = med.dens.exp.data$CANVOL2015)
plot(med.dens.exp.pp)
med.times <- list()
med.sorted <- rev(sort(med.dens.exp.pp$marks))
for(i in 1:(med.dens.exp.pp$n - 1)){
  med.times[[i]] <- med.sorted[1] - med.sorted[i+1]
}
med.times <- c(0, unlist(med.times))
med.dens <- data.frame(x = med.dens.exp.pp$x, y = med.dens.exp.pp$y, size = med.dens.exp.pp$marks)
med.dens <- med.dens[order(med.dens$size, decreasing = TRUE),]
med.dens <- data.frame(x = med.dens$x, y = med.dens$y, time = med.times, size = med.dens$size)
results_file <- paste0("./1_simulation_study/size_models/med_dens_data_exp.csv")
write_csv(med.dens, file = results_file)


## Obtain and save the high density expanded area
high.win.exp <- c((all.pat.15[[274]]$window$xrange[1] - 10),
                  (all.pat.15[[274]]$window$xrange[2] + 10),
                  (all.pat.15[[274]]$window$yrange[1] - 10),
                  (all.pat.15[[274]]$window$yrange[2] + 10))
high.dens.exp.data <- tree_data_2015 %>%
  filter(LCmajority == 1 & YTOP > high.win.exp[3] & YTOP < high.win.exp[4] & XTOP > high.win.exp[1] & XTOP < high.win.exp[2])
high.dens.exp.pp <- ppp(high.dens.exp.data$XTOP, high.dens.exp.data$YTOP,
                        c(high.win.exp[1], high.win.exp[2]), c(high.win.exp[3], high.win.exp[4]),
                        marks = high.dens.exp.data$CANVOL2015)
high.times <- list()
high.sorted <- rev(sort(high.dens.exp.pp$marks))
for(i in 1:(high.dens.exp.pp$n - 1)){
  high.times[[i]] <- high.sorted[1] - high.sorted[i+1]
}
high.times <- c(0, unlist(high.times))
high.dens <- data.frame(x = high.dens.exp.pp$x, y = high.dens.exp.pp$y, size = high.dens.exp.pp$marks)
high.dens <- high.dens[order(high.dens$size, decreasing = TRUE),]
high.dens <- data.frame(x = high.dens$x, y = high.dens$y, time = high.times, size = high.dens$size)
results_file <- paste0("./1_simulation_study/size_models/high_dens_data_exp.csv")
write_csv(high.dens, file = results_file)
