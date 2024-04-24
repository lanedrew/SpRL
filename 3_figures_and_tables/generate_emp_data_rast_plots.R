library(readr) ## load and save results
library(dplyr)
library(data.table)
library(ggplot2)
library(purrr)
library(sf)
library(terra)
library(patchwork)

## Set seed for reproducibility ----
set.seed(90210)

cbbPalette <- c("#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Clusters restricted to "extreme" positive growth
filename_pos_ex <- paste0("./data/linked_data_pos_growth_ex.csv")
pos_ex_data_trans <- read_csv(filename_pos_ex)
pos_ex_data_trans <- pos_ex_data_trans |> mutate(WKT = st_as_sfc(WKT))

## Identify and plot a subset of the linkage
filter_index <- which(pos_ex_data_trans$x > 326600 & pos_ex_data_trans$x < 327000 & pos_ex_data_trans$y > 4311300 & pos_ex_data_trans$y < 4311700)                              
pos_ex_data_trans[filter_index,] |> ggplot(aes(x = x, y = y, group = id, color = as.factor(file))) +
  geom_sf(aes(geometry = WKT, fill = as.factor(file)), alpha = .5) +
  geom_line(color = "black") +
  labs(x = "", y = "", color = "File", fill = "File", title = "Linkage for Extreme Positive Data Subset")




## Specify the relevant indexes
a_x <- 326096
a_y <- 4309939
b_x <- 328096
b_y <- 4311939


southness.rast <- scale(rast('./data/Snodgrass_aspect_southness_1m.tif'))
wetness.rast <- scale(rast('./data/Snodgrass_wetness_index_1m.tif'))
GDD.rast <- scale(rast('./data/Snodgrass_Degree_Days_2013_2019.tif'))
SPP.rast <- scale(rast('./data/Snodgrass_Snowpack_Persistence_DOY_2013_2019.tif'))


## Crop the rasters to D* and discard the originals
southness <- crop(southness.rast, ext(a_x, b_x, a_y, b_y))
wetness <- crop(wetness.rast, ext(a_x, b_x, a_y, b_y))
spp <- crop(SPP.rast, ext(a_x, b_x, a_y, b_y))
gdd <- crop(GDD.rast, ext(a_x, b_x, a_y, b_y))
rm(southness.rast, wetness.rast, SPP.rast, GDD.rast)


## Read in and process the specified dataset
file2015 <- read_csv("./data/UER_lidar_canopy_segmentation/crown_attributes_2015.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(WKT, XTOP, YTOP, CANVOL2015) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2015) %>%
  mutate(file = 1, WKT = st_as_sfc(WKT))
file2019 <- read_csv("./data/UER_lidar_canopy_segmentation/crown_attributes_2019.csv", show_col_types = FALSE) %>%
  filter(LCmajority == 1) %>%
  select(WKT, XTOP, YTOP, CANVOL2019) %>%
  filter(XTOP > a_x & XTOP < b_x & YTOP > a_y & YTOP < b_y) %>%
  rename(x = XTOP, y = YTOP, size = CANVOL2019) %>%
  mutate(file = 2, WKT = st_as_sfc(WKT))


all_file_data <- rbind(file2015, file2019) |> mutate(file = factor(file, levels = c("1", "2"), labels = c("2015", "2019")))

all_file_data |> ggplot(aes(x = x, y = y, color = file)) +
  geom_sf(aes(geometry = WKT, fill = file), alpha = .5) +
  scale_color_manual(values = cbbPalette[1:2]) + 
  scale_fill_manual(values = cbbPalette[1:2]) + 
  labs(x = "", y = "", color = "Scan", fill = "Scan", title = "Empirical Data Geometries") +
  theme_bw() +
  theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14), legend.position = "bottom") -> emp_data_plot


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
  theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14), legend.position = "bottom", legend.justification = "center") +
  guides(fill = guide_colorbar(title.position = "top")) -> SP_plot

GDD_df <- as.data.frame(gdd[[3]], xy = TRUE) |> rename(gdd = "2015")

GDD_df |> ggplot(aes(x = x, y = y, fill = gdd)) +
  geom_raster() +
  scale_fill_gradient(low = cbbPalette[3], high = cbbPalette[4]) +
  theme_bw() +
  labs(x = "", y = "", fill = "Growing Degree Days", title = "Growing Degree Days") +
  theme(axis.text = element_blank(), text = element_text(family = "serif", size = 14), legend.position = "bottom", legend.justification = "center") +
  guides(fill = guide_colorbar(title.position = "top")) -> GDD_plot

# design <- "
#   1123
#   1145
# "
# 
# emp_rast_plot <- emp_data_plot + south_plot + TWI_plot + SP_plot + GDD_plot + plot_layout(design = design) + 
#   plot_annotation(tag_levels = 'a', tag_suffix = ")") 

# design <- "
#   11112222
#   11112222
#   11112222
#   11112222
#   33445566
#   33445566
#   33445566
# "
# 
# emp_rast_plot <- emp_data_plot + emp_data_inset_plot + south_plot + TWI_plot + SP_plot + GDD_plot + plot_layout(design = design) + 
#   plot_annotation(tag_levels = 'a', tag_suffix = ")") 

design <- "
  1122
  1122
  3456
"

emp_rast_plot <- emp_data_plot + emp_data_inset_plot + south_plot + TWI_plot + SP_plot + GDD_plot + plot_layout(design = design) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")")


# emp_rast_plot <- (emp_data_plot + emp_data_inset_plot + plot_layout(widths = 2, heights = 2))/(south_plot + TWI_plot + SP_plot + GDD_plot + plot_layout(nrow = 1, heights = 1, widths = 1)) + 
#   plot_annotation(tag_levels = 'a', tag_suffix = ")")


ggsave(filename = "emp_rast_plot.png", plot = emp_rast_plot, path = "./plots/F23/",
       width = 30, height = 27, units = "cm", dpi = "retina")

ggsave(filename = "emp_rast_plot.pdf", device = "pdf", plot = emp_rast_plot, path = "./plots/F23/",
       width = 30, height = 27, units = "cm", dpi = "retina", family = "serif")
