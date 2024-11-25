
library(sf)
library(fasterize)
#library(raster)
library(ggplot2)
library(chron)
library(sp)
library(geosphere)
library(parallel)
library(terra)

source("code/utility/utility_functions.r")
template_latlon <- rast("data/templates/raster_template.tif")
template <- rast("data/templates/isforest.tif")

## load and merge dem files (from SRTM https://dwtkns.com/srtm30m/)
dem_files <- list.files("data/Topography/DEM")
for (d in dem_files) {
  if (exists("dem")) {
    dem <- merge(dem, rast(paste0("data/Topography/DEM/", d)))
  } else dem <- rast(paste0("data/Topography/DEM/", d))
}

#dem <- terra::project(dem, 'epsg:3310')
#plot(dem)

## load perimeter data for masking
perims <- st_read("data/perimeters/all_perimiters.shp")
perims <- st_transform(perims, st_crs(dem))
tomask <- st_buffer(perims, 1000)

## derive topographic vars
slope <- terrain(dem, v = "slope", unit = "degrees", neighbors = 8)
aspect <- terrain(dem, v = "aspect", unit = "degrees")
tpi <- focal(dem, w = 11, fun = function(x) x[72] - mean(x[-72]))


dem <- snap_to_template(dem, template_latlon)
slope <- snap_to_template(slope, template_latlon)
aspect <- snap_to_template(aspect, template_latlon)
tpi <- snap_to_template(tpi, template_latlon)

## convert to data frames
dem.df <- as.data.frame(dem, xy = TRUE, na.rm = TRUE)
slope.df <- as.data.frame(slope, xy = TRUE, na.rm = TRUE)
aspect.df <- as.data.frame(aspect, xy = TRUE, na.rm = TRUE)
tpi.df <- as.data.frame(tpi, xy = TRUE, na.rm = TRUE)

colnames(dem.df)[3] <- "elevation"
colnames(slope.df)[3] <- "slope"
colnames(aspect.df)[3] <- "aspect"
colnames(tpi.df)[3] <- "tpi"

## sort and merge data
dem.df <- dem.df[order(dem.df$x, dem.df$y),]
slope.df <- slope.df[order(slope.df$x, slope.df$y),]
aspect.df <- aspect.df[order(aspect.df$x, aspect.df$y),]
tpi.df <- tpi.df[order(tpi.df$x, tpi.df$y),]

topography <- cbind(dem.df, slope.df$slope, aspect.df$aspect, tpi.df$tpi)
colnames(topography)[4:6] <- c("slope", "aspect", "tpi")

## calculate heat load
topography$folded_aspect <- (180 - abs(topography$aspect - 180)) * pi / 180
topography$rad.lat <- topography$y * pi / 180
topography$rad.slope <- topography$slope * pi / 180
topography$heat_load <- -1.467 +
  (1.582 * cos(topography$rad.lat) * cos(topography$rad.slope)) +
  (-1.5 * cos(topography$folded_aspect) * sin(topography$rad.slope) * sin(topography$rad.lat)) +
  (-0.262 * sin(topography$rad.lat) * sin(topography$rad.slope)) +
  (0.607 * sin(topography$folded_aspect) * sin(topography$rad.slope))
topography$heat_load <- exp(topography$heat_load)

topography_raster <- rast(topography[,c(1:6,10)], type = "xyz", crs = "+init=epsg:4326")
plot(topography_raster)

topography_raster <- snap_to_template(topography_raster, template)
plot(topography_raster)

dem_snap <- topography_raster["elevation"]
slope_snap <- topography_raster["slope"]
aspect_snap <- topography_raster["aspect"]
tpi_snap <- topography_raster["tpi"]
heat_load_snap <- topography_raster["heat_load"]


terra::writeRaster(dem_snap, "data/Topography/derived/dem.tif", overwrite = TRUE)
terra::writeRaster(slope_snap, "data/Topography/derived/slope.tif", overwrite = TRUE)
terra::writeRaster(aspect_snap, "data/Topography/derived/aspect.tif", overwrite = TRUE)
terra::writeRaster(tpi_snap, "data/Topography/derived/tpi.tif", overwrite = TRUE)
terra::writeRaster(heat_load_snap, "data/Topography/derived/heat_load.tif", overwrite = TRUE)

saveRDS(dem_snap, "data/Topography/derived/dem.rds")
saveRDS(slope_snap, "data/Topography/derived/slope.rds")
saveRDS(aspect_snap, "data/Topography/derived/aspect.rds")
saveRDS(tpi_snap, "data/Topography/derived/tpi.rds")

topography_3310 <- as.data.frame(topography_raster, xy = TRUE, na.rm = TRUE)

write.csv(topography_3310, "data/Topography/derived/topography.csv")
