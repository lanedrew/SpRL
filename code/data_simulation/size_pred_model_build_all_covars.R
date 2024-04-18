library(tidyr) ## data manip
library(readr) ## load and save results
library(sf)
library(terra)
library(xgboost)
library(dials)
library(rsample)
library(recipes)
library(parsnip)
library(tune)
library(workflows)
library(yardstick)
library(MASS)
library(tmvtnorm) ## for generating tmvn random variables
library(mvtnorm) ## for generating mvn random variables
library(usedist) ## for working with distance objects
# speed up computation with parallel processing (optional)
library(doParallel)
all_cores <- parallel::detectCores(logical = FALSE)
registerDoParallel(cores = all_cores)


## Set seed for reproducibility ----
set.seed(90210)
data2015 <- as_tibble(read.csv("./data/UER_lidar_canopy_segmentation/crown_attributes_2015.csv"))

## Read in the raster data for the covariates of interest
southness.rast <- scale(rast('./data/Snodgrass_aspect_southness_1m.tif'))
wetness.rast <- scale(rast('./data/Snodgrass_wetness_index_1m.tif'))
spp.rast <- scale(rast('./data/Snodgrass_Snowpack_Persistence_DOY_2013_2019.tif'))
gdd.rast <- scale(rast('./data/Snodgrass_Degree_Days_2013_2019.tif'))

## Creates a pseudo-mortality variable
data2015$deltaPropCANVOL <- (data2015$CANVOL2019 - data2015$CANVOL2015) / (data2015$CANVOL2015)
data2015$mortCANVOL <- as.numeric(data2015$deltaPropCANVOL < -0.2)
data2015$estGROWTH <- (data2015$CANVOL2019 - data2015$CANVOL2015)/4


## High density model build
a_x <- 327096
a_y <- 4311239
b_x <- 327196
b_y <- 4311339

a_x_exp <- a_x - 15
a_y_exp <- a_y - 15
b_x_exp <- b_x + 15
b_y_exp <- b_y + 15

south.high <- crop(southness.rast, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp))
wetness.high <- crop(wetness.rast, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp))
gdd.high <- lapply(gdd.rast, crop, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp))
spp.high <- lapply(spp.rast, crop, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp))
raster_list <- c(list(south.high), list(wetness.high), gdd.high, spp.high)

high.dens.filter <- which(a_x_exp < data2015$XTOP & data2015$XTOP < b_x_exp & a_y_exp < data2015$YTOP & data2015$YTOP < b_y_exp & data2015$LCmajority == 1)
data2015_sample.high <- data2015[high.dens.filter,]


## Sample the covariate values for the latents and recruits
s.high <- cbind(data2015_sample.high$XTOP, data2015_sample.high$YTOP)
# covars.2015.high <- cbind(terra::extract(south.high, vect(s.high), method = "bilinear", ID = FALSE),
#                           terra::extract(slope.high, vect(s.high), method = "bilinear", ID = FALSE),
#                           terra::extract(wetness.high, vect(s.high), method = "bilinear", ID = FALSE),
#                           terra::extract(DEM.high, vect(s.high), method = "bilinear", ID = FALSE))
scaled.covars.2015.high <- update_covars_arma(s.high, raster_list)


#### Run the sampler and obtain the results ####
H <- data2015_sample.high$ZTOP
S.1 <- data2015_sample.high$CANVOL2015
X <- as.data.frame(scaled.covars.2015.high[,-1])
names(X) <- c("south", "wet", "gdd_2013", "gdd_2014", "gdd_2015", "gdd_2016", "gdd_2017", "gdd_2018", "gdd_2019",
              "spp_2013", "spp_2014", "spp_2015", "spp_2016", "spp_2017", "spp_2018", "spp_2019")
xy <- as.matrix(cbind(data2015_sample.high$XTOP, data2015_sample.high$YTOP))

## Obtain pseudo-spatial covariates
X$south.nbr <- rep(NA, nrow(X))
X$wet.nbr <- rep(NA, nrow(X))
X$gdd_2013.nbr <- rep(NA, nrow(X))
X$gdd_2014.nbr <- rep(NA, nrow(X))
X$gdd_2015.nbr <- rep(NA, nrow(X))
X$gdd_2016.nbr <- rep(NA, nrow(X))
X$gdd_2017.nbr <- rep(NA, nrow(X))
X$gdd_2018.nbr <- rep(NA, nrow(X))
X$gdd_2019.nbr <- rep(NA, nrow(X))
X$spp_2013.nbr <- rep(NA, nrow(X))
X$spp_2014.nbr <- rep(NA, nrow(X))
X$spp_2015.nbr <- rep(NA, nrow(X))
X$spp_2016.nbr <- rep(NA, nrow(X))
X$spp_2017.nbr <- rep(NA, nrow(X))
X$spp_2018.nbr <- rep(NA, nrow(X))
X$spp_2019.nbr <- rep(NA, nrow(X))
X$near.nbr.dist <- rep(NA, nrow(X))
X$near.nbr.num <- rep(NA, nrow(X))
X$avg.nbr.dist.15 <- rep(NA, nrow(X))
X$near.nbr.size <- rep(NA, nrow(X))
X$near.nbr.size.all <- rep(NA, nrow(X))
X$near.nbr.size.dist.ratio <- rep(NA, nrow(X))
X$x <- xy[,1]
X$y <- xy[,2]
X$age <- rep(NA, nrow(X))
colnames(xy) <- c("x", "y")
distance.matrix <- as.matrix(dist(xy, method = "euclidean"))

for(i in 1:nrow(X)){
  close.points.15 <- unique(which(distance.matrix[i,] < 15 & distance.matrix[i,] != 0))
  close.sizes.15 <- S.1[close.points.15]
  X$south.nbr[i] <- sum(X$south[close.points.15])
  X$wet.nbr[i] <- sum(X$wet[close.points.15])
  X$gdd_2013.nbr[i] <- sum(X$gdd_2013[close.points.15])
  X$gdd_2014.nbr[i] <- sum(X$gdd_2014[close.points.15])
  X$gdd_2015.nbr[i] <- sum(X$gdd_2015[close.points.15])
  X$gdd_2016.nbr[i] <- sum(X$gdd_2016[close.points.15])
  X$gdd_2017.nbr[i] <- sum(X$gdd_2017[close.points.15])
  X$gdd_2018.nbr[i] <- sum(X$gdd_2018[close.points.15])
  X$gdd_2019.nbr[i] <- sum(X$gdd_2019[close.points.15])
  X$spp_2013.nbr[i] <- sum(X$spp_2013[close.points.15])
  X$spp_2014.nbr[i] <- sum(X$spp_2014[close.points.15])
  X$spp_2015.nbr[i] <- sum(X$spp_2015[close.points.15])
  X$spp_2016.nbr[i] <- sum(X$spp_2016[close.points.15])
  X$spp_2017.nbr[i] <- sum(X$spp_2017[close.points.15])
  X$spp_2018.nbr[i] <- sum(X$spp_2018[close.points.15])
  X$spp_2019.nbr[i] <- sum(X$spp_2019[close.points.15])
  X$near.nbr.dist[i] <- min(distance.matrix[i,][-i])
  X$near.nbr.num[i] <- length(close.points.15)
  X$avg.nbr.dist.15[i] <- mean(distance.matrix[i,][close.points.15])
  if(length(close.points.15) == 0){
    X$avg.nbr.dist.15[i] <- min(distance.matrix[i,][-i])
  }
  X$near.nbr.size[i] <- S.1[unique(which(distance.matrix[i,] == X$near.nbr.dist[i]))]
  X$near.nbr.size.all[i] <- mean(close.sizes.15)
  if(length(close.points.15) == 0){
    X$near.nbr.size.all[i] <- S.1[unique(which(distance.matrix[i,] == X$near.nbr.dist[i]))]
  }
  X$near.nbr.size.dist.ratio[i] <- X$near.nbr.size[i]/X$near.nbr.dist[i]
  X$age[i] <- max(S.1) - S.1[i]
}

red.index <- which(a_x < X$x & b_x > X$x & a_y < X$y & b_y > X$y)
X.red <- X[red.index,]


## Fit the size model based on predicted heights
mod.data.size <- data.frame(size = S.1[red.index], X.red)

sprl_split <- rsample::initial_split(
  mod.data.size, 
  prop = 0.8, 
  strata = size
)

preprocessing_recipe <- 
  recipes::recipe(size ~ ., data = training(sprl_split)) %>%
  # convert categorical variables to factors
  recipes::step_string2factor(all_nominal()) %>%
  # combine low frequency factor levels
  recipes::step_other(all_nominal(), threshold = 0.01) %>%
  # remove no variance predictors which provide no predictive information 
  # recipes::step_nzv(all_nominal()) %>%
  prep()


sprl_cv_folds <- 
  recipes::bake(
    preprocessing_recipe, 
    new_data = training(sprl_split)
  ) %>%  
  rsample::vfold_cv(v = 5)

xgboost_model <- 
  parsnip::boost_tree(
    mode = "regression",
    trees = 1000,
    min_n = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune()
  ) %>%
  set_engine("xgboost", objective = "reg:squarederror")

xgboost_params <- 
  dials::parameters(
    min_n(),
    tree_depth(),
    learn_rate(),
    loss_reduction()
  )

xgboost_grid <- 
  dials::grid_max_entropy(
    xgboost_params, 
    size = 200
  )

xgboost_wf <- 
  workflows::workflow() %>%
  add_model(xgboost_model) %>% 
  add_formula(size ~ .)

xgboost_tuned <- tune::tune_grid(
  object = xgboost_wf,
  resamples = sprl_cv_folds,
  grid = xgboost_grid,
  metrics = yardstick::metric_set(rmse, rsq, mae),
  control = tune::control_grid(verbose = TRUE)
)

xgboost_best_params <- xgboost_tuned %>%
  tune::select_best("mae")

xgboost_model_final <- xgboost_model %>% 
  finalize_model(xgboost_best_params)


train_processed <- bake(preprocessing_recipe,  new_data = training(sprl_split))
train_prediction <- xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = size ~ ., 
    data    = train_processed
  ) %>%
  # predict the sale prices for the training data
  predict(new_data = train_processed) %>%
  bind_cols(training(sprl_split))

xgboost_score_train <- 
  train_prediction %>%
  yardstick::metrics(size, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
knitr::kable(xgboost_score_train)


test_processed  <- bake(preprocessing_recipe, new_data = testing(sprl_split))
test_prediction <- xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = size ~ ., 
    data    = train_processed
  ) %>%
  # use the training model fit to predict the test data
  predict(new_data = test_processed) %>%
  bind_cols(testing(sprl_split))

# measure the accuracy of our model using `yardstick`
xgboost_score <- 
  test_prediction %>%
  yardstick::metrics(size, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
knitr::kable(xgboost_score)

size.mod <- xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = size ~ ., 
    data    = mod.data.size
  )

all_prediction <- size.mod %>%
  # use the training model fit to predict the test data
  predict(new_data = mod.data.size) %>%
  bind_cols(mod.data.size)

# measure the accuracy of our model using `yardstick`
xgboost_score <- 
  all_prediction %>%
  yardstick::metrics(size, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
knitr::kable(xgboost_score)

save(size.mod, file = './code/data_simulation/size_models/high_dens_size_mod2.RData')





## Medium density model build

a_x <- 326996
a_y <- 4311239
b_x <- 327096
b_y <- 4311339

a_x_exp <- a_x - 15
a_y_exp <- a_y - 15
b_x_exp <- b_x + 15
b_y_exp <- b_y + 15

south.med <- crop(southness.rast, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp))
wetness.med <- crop(wetness.rast, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp))
gdd.med <- lapply(gdd.rast, crop, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp))
spp.med <- lapply(spp.rast, crop, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp))
raster_list <- c(list(south.med), list(wetness.med), gdd.med, spp.med)

med.dens.filter <- which(a_x_exp < data2015$XTOP & data2015$XTOP < b_x_exp & a_y_exp < data2015$YTOP & data2015$YTOP < b_y_exp & data2015$LCmajority == 1)
data2015_sample.med <- data2015[med.dens.filter,]

s.med <- cbind(data2015_sample.med$XTOP, data2015_sample.med$YTOP)
scaled.covars.2015.med <- update_covars_arma(s.med, raster_list)

H <- data2015_sample.med$ZTOP
S.1 <- data2015_sample.med$CANVOL2015
X <- as.data.frame(scaled.covars.2015.med[,-1])
names(X) <- c("south", "wet", "gdd_2013", "gdd_2014", "gdd_2015", "gdd_2016", "gdd_2017", "gdd_2018", "gdd_2019",
              "spp_2013", "spp_2014", "spp_2015", "spp_2016", "spp_2017", "spp_2018", "spp_2019")
xy <- as.matrix(cbind(data2015_sample.med$XTOP, data2015_sample.med$YTOP))

## Obtain pseudo-spatial covariates
X$south.nbr <- rep(NA, nrow(X))
X$wet.nbr <- rep(NA, nrow(X))
X$gdd_2013.nbr <- rep(NA, nrow(X))
X$gdd_2014.nbr <- rep(NA, nrow(X))
X$gdd_2015.nbr <- rep(NA, nrow(X))
X$gdd_2016.nbr <- rep(NA, nrow(X))
X$gdd_2017.nbr <- rep(NA, nrow(X))
X$gdd_2018.nbr <- rep(NA, nrow(X))
X$gdd_2019.nbr <- rep(NA, nrow(X))
X$spp_2013.nbr <- rep(NA, nrow(X))
X$spp_2014.nbr <- rep(NA, nrow(X))
X$spp_2015.nbr <- rep(NA, nrow(X))
X$spp_2016.nbr <- rep(NA, nrow(X))
X$spp_2017.nbr <- rep(NA, nrow(X))
X$spp_2018.nbr <- rep(NA, nrow(X))
X$spp_2019.nbr <- rep(NA, nrow(X))
X$near.nbr.dist <- rep(NA, nrow(X))
X$near.nbr.num <- rep(NA, nrow(X))
X$avg.nbr.dist.15 <- rep(NA, nrow(X))
X$near.nbr.size <- rep(NA, nrow(X))
X$near.nbr.size.all <- rep(NA, nrow(X))
X$near.nbr.size.dist.ratio <- rep(NA, nrow(X))
X$x <- xy[,1]
X$y <- xy[,2]
X$age <- rep(NA, nrow(X))
colnames(xy) <- c("x", "y")
distance.matrix <- as.matrix(dist(xy, method = "euclidean"))

for(i in 1:nrow(X)){
  close.points.15 <- unique(which(distance.matrix[i,] < 15 & distance.matrix[i,] != 0))
  close.sizes.15 <- S.1[close.points.15]
  X$south.nbr[i] <- sum(X$south[close.points.15])
  X$wet.nbr[i] <- sum(X$wet[close.points.15])
  X$gdd_2013.nbr[i] <- sum(X$gdd_2013[close.points.15])
  X$gdd_2014.nbr[i] <- sum(X$gdd_2014[close.points.15])
  X$gdd_2015.nbr[i] <- sum(X$gdd_2015[close.points.15])
  X$gdd_2016.nbr[i] <- sum(X$gdd_2016[close.points.15])
  X$gdd_2017.nbr[i] <- sum(X$gdd_2017[close.points.15])
  X$gdd_2018.nbr[i] <- sum(X$gdd_2018[close.points.15])
  X$gdd_2019.nbr[i] <- sum(X$gdd_2019[close.points.15])
  X$spp_2013.nbr[i] <- sum(X$spp_2013[close.points.15])
  X$spp_2014.nbr[i] <- sum(X$spp_2014[close.points.15])
  X$spp_2015.nbr[i] <- sum(X$spp_2015[close.points.15])
  X$spp_2016.nbr[i] <- sum(X$spp_2016[close.points.15])
  X$spp_2017.nbr[i] <- sum(X$spp_2017[close.points.15])
  X$spp_2018.nbr[i] <- sum(X$spp_2018[close.points.15])
  X$spp_2019.nbr[i] <- sum(X$spp_2019[close.points.15])
  X$near.nbr.dist[i] <- min(distance.matrix[i,][-i])
  X$near.nbr.num[i] <- length(close.points.15)
  X$avg.nbr.dist.15[i] <- mean(distance.matrix[i,][close.points.15])
  if(length(close.points.15) == 0){
    X$avg.nbr.dist.15[i] <- min(distance.matrix[i,][-i])
  }
  X$near.nbr.size[i] <- S.1[unique(which(distance.matrix[i,] == X$near.nbr.dist[i]))]
  X$near.nbr.size.all[i] <- mean(close.sizes.15)
  if(length(close.points.15) == 0){
    X$near.nbr.size.all[i] <- S.1[unique(which(distance.matrix[i,] == X$near.nbr.dist[i]))]
  }
  X$near.nbr.size.dist.ratio[i] <- X$near.nbr.size[i]/X$near.nbr.dist[i]
  X$age[i] <- max(S.1) - S.1[i]
}

red.index <- which(a_x < X$x & b_x > X$x & a_y < X$y & b_y > X$y)
X.red <- X[red.index,]


## Fit the size model based on predicted heights
mod.data.size <- data.frame(size = S.1[red.index], X.red)

sprl_split <- rsample::initial_split(
  mod.data.size, 
  prop = 0.8, 
  strata = size
)

preprocessing_recipe <- 
  recipes::recipe(size ~ ., data = training(sprl_split)) %>%
  # convert categorical variables to factors
  recipes::step_string2factor(all_nominal()) %>%
  # combine low frequency factor levels
  recipes::step_other(all_nominal(), threshold = 0.01) %>%
  # remove no variance predictors which provide no predictive information 
  # recipes::step_nzv(all_nominal()) %>%
  prep()


sprl_cv_folds <- 
  recipes::bake(
    preprocessing_recipe, 
    new_data = training(sprl_split)
  ) %>%  
  rsample::vfold_cv(v = 5)

xgboost_model <- 
  parsnip::boost_tree(
    mode = "regression",
    trees = 1000,
    min_n = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune()
  ) %>%
  set_engine("xgboost", objective = "reg:squarederror")

xgboost_params <- 
  dials::parameters(
    min_n(),
    tree_depth(),
    learn_rate(),
    loss_reduction()
  )

xgboost_grid <- 
  dials::grid_max_entropy(
    xgboost_params, 
    size = 200
  )

xgboost_wf <- 
  workflows::workflow() %>%
  add_model(xgboost_model) %>% 
  add_formula(size ~ .)

xgboost_tuned <- tune::tune_grid(
  object = xgboost_wf,
  resamples = sprl_cv_folds,
  grid = xgboost_grid,
  metrics = yardstick::metric_set(rmse, rsq, mae),
  control = tune::control_grid(verbose = TRUE)
)

xgboost_best_params <- xgboost_tuned %>%
  tune::select_best("mae")

xgboost_model_final <- xgboost_model %>% 
  finalize_model(xgboost_best_params)


train_processed <- bake(preprocessing_recipe,  new_data = training(sprl_split))
train_prediction <- xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = size ~ ., 
    data    = train_processed
  ) %>%
  # predict the sale prices for the training data
  predict(new_data = train_processed) %>%
  bind_cols(training(sprl_split))
xgboost_score_train <- 
  train_prediction %>%
  yardstick::metrics(size, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
knitr::kable(xgboost_score_train)


test_processed  <- bake(preprocessing_recipe, new_data = testing(sprl_split))
test_prediction <- xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = size ~ ., 
    data    = train_processed
  ) %>%
  # use the training model fit to predict the test data
  predict(new_data = test_processed) %>%
  bind_cols(testing(sprl_split))
# measure the accuracy of our model using `yardstick`
xgboost_score <- 
  test_prediction %>%
  yardstick::metrics(size, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
knitr::kable(xgboost_score)

size.mod <- xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = size ~ ., 
    data    = mod.data.size
  )

all_prediction <- size.mod %>%
  # use the training model fit to predict the test data
  predict(new_data = mod.data.size) %>%
  bind_cols(mod.data.size)
# measure the accuracy of our model using `yardstick`
xgboost_score <- 
  all_prediction %>%
  yardstick::metrics(size, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
knitr::kable(xgboost_score)

save(size.mod, file = './code/data_simulation/size_models/med_dens_size_mod2.RData')





## Low density model build

a_x <- 326496
a_y <- 4311439
b_x <- 326596
b_y <- 4311539

a_x_exp <- a_x - 15
a_y_exp <- a_y - 15
b_x_exp <- b_x + 15
b_y_exp <- b_y + 15

south.low <- crop(southness.rast, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp))
wetness.low <- crop(wetness.rast, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp))
gdd.low <- lapply(gdd.rast, crop, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp))
spp.low <- lapply(spp.rast, crop, ext(a_x_exp, b_x_exp, a_y_exp, b_y_exp))
raster_list <- c(list(south.low), list(wetness.low), gdd.low, spp.low)

low.dens.filter <- which(a_x_exp < data2015$XTOP & data2015$XTOP < b_x_exp & a_y_exp < data2015$YTOP & data2015$YTOP < b_y_exp & data2015$LCmajority == 1)
data2015_sample.low <- data2015[low.dens.filter,]

# Sample the covariate values for the latents and recruits
s.low <- cbind(data2015_sample.low$XTOP, data2015_sample.low$YTOP)
scaled.covars.2015.low <- update_covars_arma(s.low, raster_list)

H <- data2015_sample.low$ZTOP
S.1 <- data2015_sample.low$CANVOL2015
X <- as.data.frame(scaled.covars.2015.low[,-1])
names(X) <- c("south", "wet", "gdd_2013", "gdd_2014", "gdd_2015", "gdd_2016", "gdd_2017", "gdd_2018", "gdd_2019",
              "spp_2013", "spp_2014", "spp_2015", "spp_2016", "spp_2017", "spp_2018", "spp_2019")
xy <- as.matrix(cbind(data2015_sample.low$XTOP, data2015_sample.low$YTOP))

## Obtain pseudo-spatial covariates
X$south.nbr <- rep(NA, nrow(X))
X$wet.nbr <- rep(NA, nrow(X))
X$gdd_2013.nbr <- rep(NA, nrow(X))
X$gdd_2014.nbr <- rep(NA, nrow(X))
X$gdd_2015.nbr <- rep(NA, nrow(X))
X$gdd_2016.nbr <- rep(NA, nrow(X))
X$gdd_2017.nbr <- rep(NA, nrow(X))
X$gdd_2018.nbr <- rep(NA, nrow(X))
X$gdd_2019.nbr <- rep(NA, nrow(X))
X$spp_2013.nbr <- rep(NA, nrow(X))
X$spp_2014.nbr <- rep(NA, nrow(X))
X$spp_2015.nbr <- rep(NA, nrow(X))
X$spp_2016.nbr <- rep(NA, nrow(X))
X$spp_2017.nbr <- rep(NA, nrow(X))
X$spp_2018.nbr <- rep(NA, nrow(X))
X$spp_2019.nbr <- rep(NA, nrow(X))
X$near.nbr.dist <- rep(NA, nrow(X))
X$near.nbr.num <- rep(NA, nrow(X))
X$avg.nbr.dist.15 <- rep(NA, nrow(X))
X$near.nbr.size <- rep(NA, nrow(X))
X$near.nbr.size.all <- rep(NA, nrow(X))
X$near.nbr.size.dist.ratio <- rep(NA, nrow(X))
X$x <- xy[,1]
X$y <- xy[,2]
X$age <- rep(NA, nrow(X))
colnames(xy) <- c("x", "y")
distance.matrix <- as.matrix(dist(xy, method = "euclidean"))

for(i in 1:nrow(X)){
  close.points.15 <- unique(which(distance.matrix[i,] < 15 & distance.matrix[i,] != 0))
  close.sizes.15 <- S.1[close.points.15]
  X$south.nbr[i] <- sum(X$south[close.points.15])
  X$wet.nbr[i] <- sum(X$wet[close.points.15])
  X$gdd_2013.nbr[i] <- sum(X$gdd_2013[close.points.15])
  X$gdd_2014.nbr[i] <- sum(X$gdd_2014[close.points.15])
  X$gdd_2015.nbr[i] <- sum(X$gdd_2015[close.points.15])
  X$gdd_2016.nbr[i] <- sum(X$gdd_2016[close.points.15])
  X$gdd_2017.nbr[i] <- sum(X$gdd_2017[close.points.15])
  X$gdd_2018.nbr[i] <- sum(X$gdd_2018[close.points.15])
  X$gdd_2019.nbr[i] <- sum(X$gdd_2019[close.points.15])
  X$spp_2013.nbr[i] <- sum(X$spp_2013[close.points.15])
  X$spp_2014.nbr[i] <- sum(X$spp_2014[close.points.15])
  X$spp_2015.nbr[i] <- sum(X$spp_2015[close.points.15])
  X$spp_2016.nbr[i] <- sum(X$spp_2016[close.points.15])
  X$spp_2017.nbr[i] <- sum(X$spp_2017[close.points.15])
  X$spp_2018.nbr[i] <- sum(X$spp_2018[close.points.15])
  X$spp_2019.nbr[i] <- sum(X$spp_2019[close.points.15])
  X$near.nbr.dist[i] <- min(distance.matrix[i,][-i])
  X$near.nbr.num[i] <- length(close.points.15)
  X$avg.nbr.dist.15[i] <- mean(distance.matrix[i,][close.points.15])
  if(length(close.points.15) == 0){
    X$avg.nbr.dist.15[i] <- min(distance.matrix[i,][-i])
  }
  X$near.nbr.size[i] <- S.1[unique(which(distance.matrix[i,] == X$near.nbr.dist[i]))]
  X$near.nbr.size.all[i] <- mean(close.sizes.15)
  if(length(close.points.15) == 0){
    X$near.nbr.size.all[i] <- S.1[unique(which(distance.matrix[i,] == X$near.nbr.dist[i]))]
  }
  X$near.nbr.size.dist.ratio[i] <- X$near.nbr.size[i]/X$near.nbr.dist[i]
  X$age[i] <- max(S.1) - S.1[i]
}

red.index <- which(a_x < X$x & b_x > X$x & a_y < X$y & b_y > X$y)
X.red <- X[red.index,]


## Fit the size model based on predicted heights
mod.data.size <- data.frame(size = S.1[red.index], X.red)

sprl_split <- rsample::initial_split(
  mod.data.size, 
  prop = 0.8, 
  strata = size
)

preprocessing_recipe <- 
  recipes::recipe(size ~ ., data = training(sprl_split)) %>%
  # convert categorical variables to factors
  recipes::step_string2factor(all_nominal()) %>%
  # combine low frequency factor levels
  recipes::step_other(all_nominal(), threshold = 0.01) %>%
  # remove no variance predictors which provide no predictive information 
  # recipes::step_nzv(all_nominal()) %>%
  prep()


sprl_cv_folds <- 
  recipes::bake(
    preprocessing_recipe, 
    new_data = training(sprl_split)
  ) %>%  
  rsample::vfold_cv(v = 5)

xgboost_model <- 
  parsnip::boost_tree(
    mode = "regression",
    trees = 1000,
    min_n = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune()
  ) %>%
  set_engine("xgboost", objective = "reg:squarederror")

xgboost_params <- 
  dials::parameters(
    min_n(),
    tree_depth(),
    learn_rate(),
    loss_reduction()
  )

xgboost_grid <- 
  dials::grid_max_entropy(
    xgboost_params, 
    size = 200
  )

xgboost_wf <- 
  workflows::workflow() %>%
  add_model(xgboost_model) %>% 
  add_formula(size ~ .)

xgboost_tuned <- tune::tune_grid(
  object = xgboost_wf,
  resamples = sprl_cv_folds,
  grid = xgboost_grid,
  metrics = yardstick::metric_set(rmse, rsq, mae),
  control = tune::control_grid(verbose = TRUE)
)

xgboost_best_params <- xgboost_tuned %>%
  tune::select_best("mae")

xgboost_model_final <- xgboost_model %>% 
  finalize_model(xgboost_best_params)


train_processed <- bake(preprocessing_recipe,  new_data = training(sprl_split))
train_prediction <- xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = size ~ ., 
    data    = train_processed
  ) %>%
  # predict the sale prices for the training data
  predict(new_data = train_processed) %>%
  bind_cols(training(sprl_split))
xgboost_score_train <- 
  train_prediction %>%
  yardstick::metrics(size, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
knitr::kable(xgboost_score_train)


test_processed  <- bake(preprocessing_recipe, new_data = testing(sprl_split))
test_prediction <- xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = size ~ ., 
    data    = train_processed
  ) %>%
  # use the training model fit to predict the test data
  predict(new_data = test_processed) %>%
  bind_cols(testing(sprl_split))
# measure the accuracy of our model using `yardstick`
xgboost_score <- 
  test_prediction %>%
  yardstick::metrics(size, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
knitr::kable(xgboost_score)

size.mod <- xgboost_model_final %>%
  # fit the model on all the training data
  fit(
    formula = size ~ ., 
    data    = mod.data.size
  )

all_prediction <- size.mod %>%
  # use the training model fit to predict the test data
  predict(new_data = mod.data.size) %>%
  bind_cols(mod.data.size)
# measure the accuracy of our model using `yardstick`
xgboost_score <- 
  all_prediction %>%
  yardstick::metrics(size, .pred) %>%
  mutate(.estimate = format(round(.estimate, 2), big.mark = ","))
knitr::kable(xgboost_score)

save(size.mod, file = './code/data_simulation/size_models/low_dens_size_mod2.RData')


int.rad.mod <- lm(sqrt(data2015_sample.low$area/pi) ~ poly(CANVOL2015, 3), data = data2015_sample.low)
save(int.rad.mod, file = './code/data_simulation/size_models/low_dens_int_rad.RData')

int.rad.mod <- lm(sqrt(data2015_sample.med$area/pi) ~ poly(CANVOL2015, 3), data = data2015_sample.med)
save(int.rad.mod, file = './code/data_simulation/size_models/med_dens_int_rad.RData')

int.rad.mod <- lm(sqrt(data2015_sample.high$area/pi) ~ poly(CANVOL2015, 3), data = data2015_sample.high)
save(int.rad.mod, file = './code/data_simulation/size_models/high_dens_int_rad.RData')