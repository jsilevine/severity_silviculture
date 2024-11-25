
library(speedglm)
library(data.table)
library(ggplot2)
library(terra)
library(raster)
library(sf)
library(tidyterra)
library(ggnewscale)
library(cowplot)
library(maps)
library(ggspatial)
library(ggnewscale)

source("code/utility/utility_functions.r")

template <- rast("data/templates/raster_template.tif")

combo_inv <- terra::mosaic(dixie_inv, north_complex_inv, sheep_inv, sugar_inv)
combo_sev <- terra::mosaic(north_complex_severity, dixie_severity, walker_severity, sheep_severity, sugar_severity)

perims_all <- st_read("data/perimeters/all_perimiters.shp")

dixie_perim <- readRDS("data/perimeters/dixie_perim.rds")
dixie_perim <- st_cast(dixie_perim, "POLYGON")
dixie_perim$area <- st_area(dixie_perim)
dixie_perim <- dixie_perim[as.numeric(dixie_perim$area) > 37756431,]

north_complex_perim <- readRDS("data/perimeters/northcomplex_perim.rds")
north_complex_perim <- st_cast(north_complex_perim, "POLYGON")
north_complex_perim$area <- st_area(north_complex_perim)
north_complex_perim <- north_complex_perim[as.numeric(north_complex_perim$area) > 127681018,]

sheep_perim <- readRDS("data/perimeters/sheep_perim.rds")
sheep_perim <- st_cast(sheep_perim, "POLYGON")
sheep_perim$area <- st_area(sheep_perim)
max(sheep_perim$area)
sheep_perim <- sheep_perim[as.numeric(sheep_perim$area) > 11940259,]

sugar_perim <- readRDS("data/perimeters/sugar_perim.rds")
sugar_perim <- st_cast(sugar_perim, "POLYGON")
sugar_perim$area <- st_area(sugar_perim)
max(sugar_perim$area)
sugar_perim <- sugar_perim[as.numeric(sugar_perim$area) > 42404068,]

walker_perim <- readRDS("data/perimeters/walker_perim.rds")
walker_perim <- st_cast(walker_perim, "POLYGON")
walker_perim$area <- st_area(walker_perim)
max(walker_perim$area)
walker_perim <- walker_perim[as.numeric(walker_perim$area) > 22028674,]

perims <- rbind(dixie_perim, north_complex_perim, sheep_perim, sugar_perim, walker_perim)
perims <- st_transform(perims, st_crs(3310))


data <- fread("data/complete_data.csv")

data <- data[,.(x,y,cbi,fire_name, slope,tpi,heat_load, em, cwd, datetime, datetime_formatted,
                prev_sev, hdw, avg_fuel_moisture, own_int, own_type, clust_30, mean_dens_30,
                mean_ht_30, mean_area_30, clust_180, mean_dens_180, mean_ht_180, mean_area_180)]
full_data_raster <- rast(data, type = "xyz", crs = "+init=epsg:3310")

plot(full_data_raster)

slope <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data, aes(x = x, y = y, fill = slope)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Slope (%)",
    palette = "viridi",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
slope

tpi <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data[,], aes(x = x, y = y, fill = tpi)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Topographic position index (m)",
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
tpi

heat_load <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data[,], aes(x = x, y = y, fill = heat_load)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Heat Load",
    palette = "viridi",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
heat_load

cwd <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data[,], aes(x = x, y = y, fill = cwd)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Climate Water Deficit",
    palette = "viridi",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
cwd

prev_sev <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data[,], aes(x = x, y = y, fill = prev_sev)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Previous Severity",
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
prev_sev


for (f in unique(data$fire_name)) {
  data[data$fire_name == f, "datetime_norm"] <- ((data[data$fire_name == f, "datetime"] - min(data[data$fire_name == f, "datetime"])) /
    (max(data[data$fire_name == f, "datetime"]) - min(data[data$fire_name == f, "datetime"]))) * 100
}

burn_time <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data, aes(x = x, y = y, fill = datetime_norm)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Time of burn (% of full fire length)",
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
burn_time

hs <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data, aes(x = x, y = y, fill = cbi)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "CBI",
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
hs

hdw <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data, aes(x = x, y = y, fill = hdw)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Hot-dry-windy index",
    palette = "viridi",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
hdw

fuel_moisture <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data, aes(x = x, y = y, fill = avg_fuel_moisture)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Fuel moisture (%)",
    palette = "viridi",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
fuel_moisture


dens_30 <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data[data$mean_dens_30 < 400,], aes(x = x, y = y, fill = mean_dens_30)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Mean Stem Density (ind. ha^-1)",
    palette = "viridi",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
dens_30

dens_180 <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data[,], aes(x = x, y = y, fill = mean_dens_180)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Mean Stem Density (ind. ha^-1)",
    palette = "viridi",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
dens_180

clust_30 <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data[data$clust_30 < 1.8 & data$clust_30 > 0.4,], aes(x = x, y = y, fill = clust_30)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Spatial Homogeneity",
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
clust_30

clust_180 <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data[,], aes(x = x, y = y, fill = clust_180)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Spatial Homogeneity",
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
clust_180



mean_ht_30 <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data[,], aes(x = x, y = y, fill = mean_ht_30)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Mean Stem Height (m)",
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
mean_ht_30

mean_ht_180 <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data[,], aes(x = x, y = y, fill = mean_ht_180)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Mean Stem Height (m)",
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 8,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
mean_ht_180



mean_area_30 <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data[,], aes(x = x, y = y, fill = mean_area_30)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Mean Gap Area (m^2)",
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 12,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent",
    trans = "log10"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
mean_area_30

mean_area_180 <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data[,], aes(x = x, y = y, fill = mean_area_180)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Mean Gap Area (m^2)",
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 12,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent",
    trans = "log10"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
mean_area_180


em <- ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data[,], aes(x = x, y = y, fill = em)) +
  scale_y_continuous(limits = c(NA, 260000)) +
  scale_fill_whitebox_c(
    name = "Ladder Fuels Index",
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 12,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.position = "bottom")
em



p1 <- plot_grid(cwd + theme(axis.title = element_blank()),
                slope + theme(axis.title = element_blank()),
                tpi + theme(axis.title = element_blank()),
                heat_load + theme(axis.title = element_blank()),
                align ="hv")
p1
ggsave("plots/figureS3a.png")

p1 <- plot_grid(cwd + theme(axis.title = element_blank(),
                            legend.position = "none"),
                slope + theme(axis.title = element_blank(),
                            legend.position = "none"),
                tpi + theme(axis.title = element_blank(),
                            legend.position = "none"),
                heat_load + theme(axis.title = element_blank(),
                            legend.position = "none"),
                align ="hv")
p1
ggsave("plots/figureS3_topography.png")


p2 <- plot_grid(hdw + theme(axis.title = element_blank()),
                fuel_moisture + theme(axis.title = element_blank()),
                prev_sev + theme(axis.title = element_blank()),
                burn_time + theme(axis.title = element_blank()),
                align ="hv")
p2


p2 <- plot_grid(hdw + theme(axis.title = element_blank(),
                            legend.position = "none"),
                fuel_moisture + theme(axis.title = element_blank(),
                            legend.position = "none"),
                prev_sev + theme(axis.title = element_blank(),
                            legend.position = "none"),
                burn_time + theme(axis.title = element_blank(),
                            legend.position = "none"),
                align ="hv")
p2

ggsave("plots/figureS3_weather.png")


p3 <- plot_grid(dens_30 + theme(axis.title = element_blank()),
                clust_30 + theme(axis.title = element_blank()),
                mean_ht_30 + theme(axis.title = element_blank()),
                mean_area_30 + theme(axis.title = element_blank()),
                align ="hv")
p3


p3 <- plot_grid(dens_30 + theme(axis.title = element_blank(),
                            legend.position = "none"),
                clust_30 + theme(axis.title = element_blank(),
                            legend.position = "none"),
                mean_ht_30 + theme(axis.title = element_blank(),
                            legend.position = "none"),
                mean_area_30 + theme(axis.title = element_blank(),
                            legend.position = "none"),
                align ="hv")
p3

ggsave("plots/figureS3_fs30.png")

p4 <- plot_grid(dens_180 + theme(axis.title = element_blank(),
                            legend.position = "none"),
                clust_180 + theme(axis.title = element_blank(),
                            legend.position = "none"),
                mean_ht_180 + theme(axis.title = element_blank(),
                            legend.position = "none"),
                mean_area_180 + theme(axis.title = element_blank(),
                            legend.position = "none"),
                align ="hv")
p4

ggsave("plots/figureS3_fs180.png")

p5 <- plot_grid(em + theme(axis.title = element_blank(),
                            legend.position = "none"),
                hs + theme(axis.title = element_blank(),
                            legend.position = "none"),
                align ="hv")
p5

mean_area_180

ggsave("plots/figureS3_misc.png")
