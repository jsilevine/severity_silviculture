library(rgdal)
library(sf)
library(fasterize)
library(raster)
library(ggplot2)
library(chron)
library(sp)
library(geosphere)
library(parallel)
library(terra)
setwd("../")

## load and merge dem files
dem_files <- list.files("data/Topography/DEM")
for (d in dem_files) {
  if (exists("dem")) {
    dem <- merge(dem, raster(paste0("data/Topography/DEM/", d)))
  } else dem <- raster(paste0("data/Topography/DEM/", d))
}

## load perimeter data for masking
perims <- st_read("data/perimeters/all_perimiters.shp")
perims <- st_transform(perims, st_crs(dem))
tomask <- st_buffer(perims, 1000)

## derive topographic vars
slope <- terrain(dem, v = "slope", unit = "degrees", neighbors = 8)
aspect <- terrain(dem, opt = "aspect", unit = "degrees")
tpi <- raster(focal(rast(dem), w = 11, fun = function(x) x[72] - mean(x[-72])))

## mask data and save
dem <- mask(dem, as_Spatial(tomask))
slope <- mask(slope, as_Spatial(tomask))
aspect <- mask(aspect, as_Spatial(tomask))
tpi <- mask(tpi, as_Spatial(tomask))

## resample to daily progression grid
togrid <- readRDS("data/fire_progression//dixie_daily_burned.rds")
togrid <- merge(togrid, readRDS("data/fire_progression//northcomplex_daily_burned.rds"))
togrid <- merge(togrid, readRDS("data/fire_progression//sheep_daily_burned.rds"))
togrid <- merge(togrid, readRDS("data/fire_progression//walker_daily_burned.rds"))
togrid <- merge(togrid, readRDS("data/fire_progression//sugar_daily_burned.rds"))

dem <- resample(dem, togrid)
slope <- resample(slope, togrid)
aspect <- resample(aspect, togrid)
tpi <- resample(tpi, togrid)

saveRDS(dem, "data/Topography/derived/dem.rds")
saveRDS(slope, "data/Topography/derived/slope.rds")
saveRDS(aspect, "data/Topography/derived/aspect.rds")
saveRDS(tpi, "data/Topography/derived/tpi.rds")

## convert to data frames
dem.df <- as.data.frame(dem, xy = TRUE, na.rm = TRUE)
slope.df <- as.data.frame(slope, xy = TRUE, na.rm = TRUE)
aspect.df <- as.data.frame(aspect, xy = TRUE, na.rm = TRUE)
tpi.df <- as.data.frame(tpi, xy = TRUE, na.rm = TRUE)

colnames(dem.df)[3] <- "elevation"
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
write.csv(topography, "data/Topography/derived/topography.csv")
