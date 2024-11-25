
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

data <- fread("data/complete_data.csv")
data <- data[complete.cases(data),]

template <- rast("data/templates/raster_template.tif")

## load severity data

dixie_s <- rast("data/cbi/DIXIE_2021_CBI_bc.tif")
dixie_inv <- snap_to_template(dixie_s, template, inverse = TRUE)
dixie_inv[!is.na(dixie_inv)] <- 1

north_complex_s <- rast("data/cbi/NORTH_COMPLEX_2020_CBI_bc.tif")
north_complex_inv <- snap_to_template(north_complex_s, template, inverse = TRUE)
north_complex_inv[!is.na(north_complex_inv)] <- 1

sheep_s <- rast("data/cbi/SHEEP_2020_CBI_bc.tif")
sheep_inv <- snap_to_template(sheep_s, template, inverse = TRUE)
sheep_inv[!is.na(sheep_inv)] <- 1

sugar_s <- rast("data/cbi/SUGAR_2021_CBI_bc.tif")
sugar_inv <- snap_to_template(sugar_s, template, inverse = TRUE)
sugar_inv[!is.na(sugar_inv)] <- 1

combo_inv <- terra::mosaic(dixie_inv, north_complex_inv, sheep_inv, sugar_inv)
combo_sev <- terra::mosaic(north_complex_severity, dixie_severity, walker_severity, sheep_severity, sugar_severity)

## load perimeter data

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

fig1a <- ggplot() +
  geom_spatraster(data = combo_inv, aes(fill = CBI_bc), alpha = 0.7) +
  scale_fill_gradient(na.value = "transparent", low = "gray", high = "gray") +
  new_scale_fill() +
  geom_spatraster(data = combo_sev, aes(fill = CBI_bc), maxcell = 5000000) +
  scale_fill_whitebox_c(
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 12,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  ) +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = NA, color = "black", linewidth = 0.4) +
  annotation_scale() +
  scale_y_continuous(limits = c(38.7, 40.4)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
fig1a

ggsave("plots/figure1a.png", fig1a)


ownership <- rast("data/ownership/ownership.tif")
ownership_df <- as.data.frame(ownership, xy = TRUE, na.rm = TRUE)
ownership_df$layer <- as.factor(ownership_df$layer)

colors <- met.brewer(name = "Egypt", n = 3)
colors
colors <- colors[c(1,3,4)]
colors

colors2 <- met.brewer(name = "Derain", n = 6)
colors2

ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = "black", linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = ownership_df, aes(x = x, y = y, fill = layer)) +
  scale_fill_manual(values = c(colors[3], colors2[1], colors[2])) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("plots/figure1e.png", bg = "transparent")

for (f in unique(data$fire_name)) {
  data[data$fire_name == f, "datetime_norm"] <- ((data[data$fire_name == f, "datetime"] - min(data[data$fire_name == f, "datetime"])) /
    (max(data[data$fire_name == f, "datetime"]) - min(data[data$fire_name == f, "datetime"]))) * 100
}


ggplot() +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "gray", color = NA, linewidth = 0.4, alpha = 0.5) +
  geom_raster(data = data, aes(x = x, y = y, fill = datetime_norm)) +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = NA, color = "black", linewidth = 0.4, alpha = 0.5) +
  scale_fill_whitebox_c(
    name = "Fire progression",
    palette = "muted",
    labels = scales::label_number(suffix = "%"),
    n.breaks = 12,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  )     +
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("plots/figure1f.png")

##---------------------------------------------------------------
## Lidar insets
##---------------------------------------------------------------

point_test <- st_read("data/lidar/points/TAO__-100000_178000__t4_0p75_lp3__HighPoints.shp")
point_test

box <- st_as_sfc(st_bbox(point_test))
box <- st_as_sf(box)
ext <- ext(box)

ext <- extent(c(-121.3, -121.1, 39.6, 39.7))

bb <- st_bbox(box)
bb <- st_transform(bb, st_crs(4326))
bb[1] <- -121.25
bb[2] <- 39.6
bb[3] <- -121.1
bb[4] <- 39.7
bb <- st_transform(bb, st_crs(3310))
ext <- ext(bb)

test <- crop(north_complex_severity, ext)
test[test$CBI_bc < 2.25] <- 0
test[test$CBI_bc >= 2.25] <- 1


ht <- rast(data[,.(x,y,mean_ht_30)], crs = "+init=epsg:3310")
ht <- crop(ht, ext(point_test))
ht_df <- as.data.frame(ht, xy = TRUE, na.rm = TRUE)

ownership <- rast("data/ownership/ownership.tif")
ownership_raster <- snap_to_template(ownership, template)
own_crop <- crop(ownership, ext)

own_crop_df <- as.data.frame(own_crop, xy = TRUE)
own_crop_df$own_class <- NA
own_crop_df[own_crop_df$layer == 3, "own_class"] <- "private"
own_crop_df[own_crop_df$layer == 2, "own_class"] <- "other"
own_crop_df[own_crop_df$layer == 1, "own_class"] <- "public"

own_crop_poly <- st_as_sf(as.polygons(own_crop))
own_crop_poly[own_crop_poly$layer == 1, "own_class"] <- "public"
own_crop_poly[own_crop_poly$layer == 2, "own_class"] <- "other"
own_crop_poly[own_crop_poly$layer == 3, "own_class"] <- "private"

t2 <- ggplot() +
  geom_sf(own_crop_poly, mapping = aes(geometry = geometry, fill = own_class)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  #scale_fill_manual(values = c(as.character(pal_unikn_light["pinky2"]), ## other
  #                             as.character(pal_unikn_light["seeblau3"]),                   ## private
  #                             as.character(pal_unikn_light["seegruen3"]))) +                ## public
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank())
t2

## st2 <- t2 +
##   new_scale_fill() +
##   scale_fill_manual(values = c("gray", "transparent")) +
##   geom_sf(data = test_poly, mapping = aes(geometry = geometry), fill = "gray", alpha = 0.4)
## t2

t2 <- t2 +
  geom_sf(data = box, fill = "transparent", color = "black", linewidth = 2)
t2


colors <- met.brewer(name = "Egypt", n = 6)
colors <- colors[3:4]

point_test <- point_test[order(point_test$height),]

t3 <- ggplot(point_test) +
  geom_sf(data = box, fill = "lightgray", alpha = 0.1) +
  geom_sf(data = st_crop(own_crop_poly, ext(point_test)), aes(fill = own_class), alpha = 0.3) +
  geom_sf(aes(color = height)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
   scale_color_whitebox_c(
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 12,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  )     +
  geom_sf(data = st_crop(own_crop_poly, ext(point_test)), fill = NA, color = "black",linewidth = 1) +
  scale_fill_manual(values = colors) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "white"))
t3

t4 <- ggplot(ht_df) +
  geom_sf(data = box, fill = "lightgray", alpha = 0.1) +
  geom_raster(aes(x = x, y = y, fill = mean_ht_30), alpha = 1.0) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_whitebox_c(
    palette = "muted",
    labels = scales::label_number(suffix = ""),
    n.breaks = 12,
    guide = guide_legend(reverse = TRUE),
    na.value = "transparent"
  )     +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "white"))
t4

fig1b <- plot_grid(t2, t3, align = "hv", rel_widths = c(0.5, 0.5))
fig1b


ggsave(filename = "plots/figure1b.svg", plot = fig1b)

fig1g <- plot_grid(t3, t4, align = "hv", rel_widths = c(0.5, 0.5))
fig1g

ggsave("plots/figure1g.svg", plot = fig1g)


##---------------------------------------------------------------
## California inset
##---------------------------------------------------------------

states <- map_data("state")
ca <- states[states$region == "california",]

fire_box <- st_as_sfc(st_bbox(point_test))
fire_box <- st_transform(box, crs = 4326)
fire_box <- st_as_sf(box)

perims <- st_make_valid(st_transform(perims, crs = 4326))


ggplot() +
  geom_polygon(data = ca, aes(x=long, y=lat, group = group), fill = "lightgray") +
  geom_sf(data = perims, mapping = aes(geometry = Shape), fill = "#3182bd", color = NA, linewidth = 0.4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("plots/figure1c.pdf")
