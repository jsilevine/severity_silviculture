##---------------------------------------------------------------
## ownership_cln.r
## BY: JACOB LEVINE - jacoblevine@princeton.edu
##---------------------------------------------------------------


##---------------------------------------------------------------
## 0. head
##---------------------------------------------------------------

library(sf)
library(raster)
library(terra)
library(fasterize)
library(ggplot2)
library(ggnewscale)
library(ggspatial)
library(mapdata)
library(ggmap)

source("code/utility/utility_functions.r")

##---------------------------------------------------------------
## 1. Body
##---------------------------------------------------------------

ownership_nonprivate <- st_read("data/ownership/ownership15_1.shp")
ownership_nonprivate <- st_transform(ownership_nonprivate, crs = st_crs(3310))
reclass <- data.frame(orig = unique(ownership_nonprivate$Own_Group),
                      mapping = c(1,1,1,1,1,1,1,1,1,1,2,1,1,1))

for (i in 1:nrow(ownership_nonprivate)) {
  ownership_nonprivate$own_int[i] <- reclass[reclass$orig == ownership_nonprivate$Own_Group[i], "mapping"]
}

ownership_private <- st_read("data/ownership//Forest_Industry_Owners_18_2.shp")
ownership_private <- st_transform(ownership_private, crs = st_crs(3310))
ownership_private$own_int <- 3

template <- rast("data/templates/isforest.tif")
ext <- terra::ext(template)

ownership_nonprivate_crop <- st_crop(ownership_nonprivate, st_bbox(ext))
ownership_private_crop <- st_crop(ownership_private, st_bbox(ext))

ownership_nonprivate_rast <- rast(fasterize(ownership_nonprivate_crop, raster(template), field = "own_int"))
ownership_private_rast <- rast(fasterize(ownership_private_crop, raster(template), field = "own_int"))
ownership <- merge(ownership_nonprivate_rast, ownership_private_rast)
plot(ownership)

ownership[is.na(ownership[])] <- 2
ownership <- snap_to_template(ownership, template)
plot(ownership)
writeRaster(ownership, "data/ownership/ownership.tif", overwrite = TRUE)

ownership_df <- as.data.frame(ownership, xy = TRUE, na.rm = TRUE)
colnames(ownership_df)[3] <- "own_int"
ownership_df[ownership_df$own_int == 1, "own_type"] <- "Federal"
ownership_df[ownership_df$own_int == 2, "own_type"] <- "Other"
ownership_df[ownership_df$own_int == 3, "own_type"] <- "Private Industrial"
write.csv(ownership_df, "data/ownership/ownership.csv")

## takes a long time to render
ggplot(ownership_df) +
  geom_tile(data = ownership_df, aes(x = x, y = y, fill = own_type)) +
  coord_sf(expand = TRUE) +
  scale_fill_manual(values = c("#bdbdbd", "#f0f0f0", "#636363")) +
  new_scale_fill() +
  geom_sf(data = perims, aes(color = FIRE_NA),
          fill = NA, linewidth = 0.65) +
  scale_color_brewer(type = "qual", palette = 3) +
  ##geom_sf_label(data = perims, aes(label = FIRE_NA), size = 3) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.title = element_blank(), legend.title = element_blank(),
        panel.grid = element_blank())


ggsave()
