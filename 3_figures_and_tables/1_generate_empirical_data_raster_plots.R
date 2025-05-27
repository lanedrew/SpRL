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
  select(XTOP, YTOP, CANVOL2015, WKT) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2015) %>%
  mutate(file = 1)
file2019 <- read_csv("./resources/empirical_data/UER_lidar_canopy_segmentation/crown_attributes_2019.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(XTOP, YTOP, CANVOL2019, WKT) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2019) %>%
  mutate(file = 2)

## Read in the raster data for the covariates of interest
southness <- scale(rast('./resources/empirical_data/Snodgrass_aspect_southness_1m.tif'))
wetness <- scale(rast('./resources/empirical_data/Snodgrass_wetness_index_1m.tif'))
gdd <- scale(rast('./resources/empirical_data/Snodgrass_Degree_Days_2013_2019.tif'))
spp <- scale(rast('./resources/empirical_data/Snodgrass_Snowpack_Persistence_DOY_2013_2019.tif'))

## Collate and plot the full empirical dataset
all_file_data <- rbind(file2015, file2019) |> mutate(file = factor(file, levels = c("1", "2"), labels = c("2015", "2019"))) |> st_as_sf(coords = c("x", "y"), wkt = "WKT")

square_corners <- data.frame(x = c(326996, 327096, 327096, 326996, 326996),
                             y = c(4311239, 4311239, 4311339, 4311339, 4311239))
all_file_data |> ggplot(aes(x = x, y = y, color = file)) +
  geom_sf(aes(geometry = WKT, fill = file), alpha = .5) +
  scale_color_manual(values = cbbPalette[1:2]) + 
  scale_fill_manual(values = cbbPalette[1:2]) + 
  geom_polygon(data = square_corners, aes(x = x, y = y), fill = NA, color = "black", linewidth = .8) +
  labs(x = "", y = "", color = "Scan", fill = "Scan", title = "Empirical Data Geometries") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(family = "serif", size = 14), legend.position = "none") -> emp_data_plot

## Plot the subsetted empirical dataset
square_corners <- data.frame(x = c(326996 + 79, 327096 - 6, 327096 - 6, 326996 + 79, 326996 + 79),
                             y = c(4311239 + 84, 4311239 + 84, 4311339 - 1, 4311339 - 1, 4311239 + 84))

all_file_data |> filter(x > a_x & x < b_x & y > a_y & y < b_y) |> 
  ggplot(aes(x = x, y = y, color = file)) +
  geom_sf(aes(geometry = WKT, fill = file), alpha = .5) +
  geom_polygon(data = square_corners, aes(x = x, y = y), fill = NA, color = "black", linewidth = 1) +
  coord_sf(xlim = c(326996, 327096), ylim = c(4311239, 4311339)) +
  scale_color_manual(values = cbbPalette[1:2]) + 
  scale_fill_manual(values = cbbPalette[1:2]) + 
  labs(x = "", y = "", color = "Scan", fill = "Scan", title = "Medium Density Inset") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(family = "serif", size = 14),
        legend.position = "bottom") -> emp_data_inset_plot

## Plot the subsetted empirical dataset
square_corners <- data.frame(x = c(326996 + 79, 326996 + 87, 326996 + 87, 326996 + 79, 326996 + 79),
                             y = c(4311239 + 84, 4311239 + 84, 4311339 - 8, 4311339 - 8, 4311239 + 84))

all_file_data |> filter(x > a_x & x < b_x & y > a_y & y < b_y) |> 
  ggplot(aes(x = x, y = y, color = file)) +
  geom_sf(aes(geometry = WKT, fill = file), alpha = .5) +
  geom_polygon(data = square_corners, aes(x = x, y = y), fill = NA, color = "black", linewidth = 1.15) +
  coord_sf(xlim = c(326996 + 79, 327096 - 6), ylim = c(4311239 + 84, 4311339 - 1)) +
  scale_color_manual(values = cbbPalette[1:2]) + 
  scale_fill_manual(values = cbbPalette[1:2]) + 
  labs(x = "", y = "", color = "Scan", fill = "Scan", title = "Medium Density Inset Closeup") +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(family = "serif", size = 14),
        legend.position = "none") -> emp_data_closeup_plot

## Extract the raster values and plot them
south_df <- as.data.frame(southness, xy = TRUE)
south_df |> ggplot(aes(x = x, y = y, fill = Snodgrass_aspect_southness_1m)) +
  geom_raster() +
  scale_fill_gradient(low = cbbPalette[3], high = cbbPalette[4]) +
  theme_bw() +
  labs(x = "", y = "", fill = "Folded Aspect", title = "Folded Aspect") +
  theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14),
        legend.position = "bottom", legend.justification = "center",
        axis.ticks = element_blank()) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = .5)) -> south_plot


TWI_df <- as.data.frame(wetness, xy = TRUE)
TWI_df |> ggplot(aes(x = x, y = y, fill = Snodgrass_wetness_index_1m)) +
  geom_raster() +
  scale_fill_gradient(low = cbbPalette[3], high = cbbPalette[4]) +
  theme_bw() +
  labs(x = "", y = "", fill = "TWI", title = "HAS Wetness Index") +
  theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14),
        legend.position = "bottom", legend.justification = "center",
        axis.ticks = element_blank()) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = .5))-> TWI_plot


SP_df <- as.data.frame(spp[[3]], xy = TRUE) |> rename(spp = "2015")
SP_df |> ggplot(aes(x = x, y = y, fill = spp)) +
  geom_raster() +
  scale_fill_gradient(low = cbbPalette[3], high = cbbPalette[4]) +
  theme_bw() +
  labs(x = "", y = "", fill = "SP", title = "Snowpack Persistence") +
  theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14),
        legend.position = "bottom", legend.justification = "center",
        axis.ticks = element_blank()) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = .5)) -> SP_plot


GDD_df <- as.data.frame(gdd[[3]], xy = TRUE) |> rename(gdd = "2015")
GDD_df |> ggplot(aes(x = x, y = y, fill = gdd)) +
  geom_raster() +
  scale_fill_gradient(low = cbbPalette[3], high = cbbPalette[4]) +
  theme_bw() +
  labs(x = "", y = "", fill = "GDD", title = "Growing Degree Days") +
  theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14),
        legend.position = "bottom", legend.justification = "center",
        axis.ticks = element_blank()) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = .5)) -> GDD_plot

emp_rast_plot <- emp_data_plot + emp_data_inset_plot + south_plot + TWI_plot + SP_plot + GDD_plot + plot_layout(design = design) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")")

top_raw <- emp_data_plot + emp_data_inset_plot + emp_data_closeup_plot + plot_layout(ncol = 3)
bottom_raw <- south_plot + TWI_plot + SP_plot + GDD_plot + plot_layout(ncol = 4)
combined_plot <- top_raw / bottom_raw + plot_layout(heights = c(4, 2)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")")

## Save the generated plot
ggsave(filename = "1_empirical_data_raster_plot.png", plot = combined_plot, path = "./3_figures_and_tables/",
       width = 35, height = 30, units = "cm", dpi = "retina")

