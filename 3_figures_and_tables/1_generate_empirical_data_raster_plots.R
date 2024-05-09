########################################################
#### This script generates Figure 1 from the paper. ####
########################################################

## Load libraries
library(readr)
library(dplyr)
library(data.table)
library(ggplot2)
library(purrr)
library(sf)
library(terra)
library(patchwork)

## Set seed for reproducibility ----
set.seed(90210)

## Set a colorblind friendly palate to use for visualizations
cbbPalette <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Define the spatial domain
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

## Read in the raster data for the covariates of interest
southness.rast <- scale(rast('./resources/empirical_data/Snodgrass_aspect_southness_1m.tif'))
wetness.rast <- scale(rast('./resources/empirical_data/Snodgrass_wetness_index_1m.tif'))
GDD.rast <- scale(rast('./resources/empirical_data/Snodgrass_Degree_Days_2013_2019.tif'))
SPP.rast <- scale(rast('./resources/empirical_data/Snodgrass_Snowpack_Persistence_DOY_2013_2019.tif'))


## Collate and plot the full empirical dataset
all_file_data <- rbind(file2015, file2019) |> mutate(file = factor(file, levels = c("1", "2"), labels = c("2015", "2019")))
all_file_data |> ggplot(aes(x = x, y = y, color = file)) +
  geom_sf(aes(geometry = WKT, fill = file), alpha = .5) +
  scale_color_manual(values = cbbPalette[1:2]) + 
  scale_fill_manual(values = cbbPalette[1:2]) + 
  labs(x = "", y = "", color = "Scan", fill = "Scan", title = "Empirical Data Geometries") +
  theme_bw() +
  theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14), legend.position = "bottom") -> emp_data_plot


## Plot the subsetted empirical dataset
a_x <- 326996
a_y <- 4311239
b_x <- 327096
b_y <- 4311339

all_file_data |> filter(x > a_x & x < b_x & y > a_y & y < b_y) |> 
  ggplot(aes(x = x, y = y, color = file)) +
  geom_sf(aes(geometry = WKT, fill = file), alpha = .5) +
  scale_color_manual(values = cbbPalette[1:2]) + 
  scale_fill_manual(values = cbbPalette[1:2]) + 
  labs(x = "", y = "", color = "Scan", fill = "Scan", title = "Medium Density Inset") +
  theme_bw() +
  theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14), legend.position = "bottom") -> emp_data_inset_plot


## Extract the raster values and plot them
south_df <- as.data.frame(southness, xy = TRUE)
south_df |> ggplot(aes(x = x, y = y, fill = Snodgrass_aspect_southness_1m)) +
  geom_raster() +
  scale_fill_gradient(low = cbbPalette[3], high = cbbPalette[4]) +
  theme_bw() +
  labs(x = "", y = "", fill = "Folded Aspect", title = "Folded Aspect") +
  theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14), legend.position = "bottom") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = .5)) -> south_plot


TWI_df <- as.data.frame(wetness, xy = TRUE)
TWI_df |> ggplot(aes(x = x, y = y, fill = Snodgrass_wetness_index_1m)) +
  geom_raster() +
  scale_fill_gradient(low = cbbPalette[3], high = cbbPalette[4]) +
  theme_bw() +
  labs(x = "", y = "", fill = "TWI", title = "HAS Wetness Index") +
  theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14), legend.position = "bottom") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = .5))-> TWI_plot


SP_df <- as.data.frame(spp[[3]], xy = TRUE) |> rename(spp = "2015")
SP_df |> ggplot(aes(x = x, y = y, fill = spp)) +
  geom_raster() +
  scale_fill_gradient(low = cbbPalette[3], high = cbbPalette[4]) +
  theme_bw() +
  labs(x = "", y = "", fill = "Snowpack Persistence", title = "Snowpack Persistence") +
  theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14),
        legend.position = "bottom", legend.justification = "center") +
  guides(fill = guide_colorbar(title.position = "top")) -> SP_plot


GDD_df <- as.data.frame(gdd[[3]], xy = TRUE) |> rename(gdd = "2015")
GDD_df |> ggplot(aes(x = x, y = y, fill = gdd)) +
  geom_raster() +
  scale_fill_gradient(low = cbbPalette[3], high = cbbPalette[4]) +
  theme_bw() +
  labs(x = "", y = "", fill = "Growing Degree Days", title = "Growing Degree Days") +
  theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14),
        legend.position = "bottom", legend.justification = "center") +
  guides(fill = guide_colorbar(title.position = "top")) -> GDD_plot


## Specify the plot layout for patchwork and combine the plots
design <- "
  1122
  1122
  3456
"

emp_rast_plot <- emp_data_plot + emp_data_inset_plot + south_plot + TWI_plot + SP_plot + GDD_plot + plot_layout(design = design) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")")

## Save the generated plot
ggsave(filename = "1_empirical_data_raster_plot.png", plot = emp_rast_plot, path = "./3_figures_and_tables/",
       width = 30, height = 27, units = "cm", dpi = "retina")
