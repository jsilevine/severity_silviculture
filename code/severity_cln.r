
library(rgdal)
library(rvest)
library(sf)
library(rmapshaper)
library(fasterize)
library(raster)
library(ggplot2)
library(gganimate)
setwd("../")

##---------------------------------------------------------------
## Clean and save final fire perimeters from LANDFIRE/FRAP
##---------------------------------------------------------------
perimeters <- st_read("data/perimeters/FRAP_perims.gdb")
perimeters <- st_transform(perimeters, crs = st_crs(4326))
perimeters <- perimeters[!is.na(perimeters$FIRE_NAME),]
dixie_perim <- st_make_valid(perimeters[perimeters$FIRE_NAME == "DIXIE" & perimeters$YEAR_ == "2021",])
northcomplex_perim <- st_make_valid(perimeters[perimeters$FIRE_NAME == "NORTH COMPLEX",])
sheep_perim <- st_make_valid(perimeters[perimeters$FIRE_NAME == "SHEEP" & perimeters$YEAR_ == "2020" & perimeters$UNIT_ID == "PNF",])
walker_perim <- st_make_valid(perimeters[perimeters$FIRE_NAME == "WALKER" & perimeters$YEAR_ == "2019",])
sugar_perim <- st_make_valid(perimeters[perimeters$FIRE_NAME == "SUGAR" & perimeters$YEAR_ == "2021",])

saveRDS(dixie_perim, "data/perimeters/dixie_perim.rds")
saveRDS(northcomplex_perim, "data/perimeters/northcomplex_perim.rds")
saveRDS(sheep_perim, "data/perimeters/sheep_perim.rds")
saveRDS(walker_perim, "data/perimeters/walker_perim.rds")
saveRDS(sugar_perim, "data/perimeters/sugar_perim.rds")

## save all perimiters as single shapefile
full_perims <- rbind(dixie_perim, northcomplex_perim, sheep_perim, walker_perim, sugar_perim)
st_write(full_perims, "data/perimeters/all_perimiters.shp", append = FALSE)

ggplot(dixie_perim) +
  ylim(c(min(dixie_daily_severity$y), max(dixie_daily_severity$y))) +
  xlim(c(min(dixie_daily_severity$x), max(dixie_daily_severity$x))) +
  geom_tile(data = dixie_daily_severity, aes(x = x, y = y, fill = time)) +
  geom_sf(alpha = 0.1, color = "gray") +
  theme_bw()

ggplot(walker_perim) +
  geom_raster(data = as.data.frame(as(rdnbr, "SpatialPixelsDataFrame")), aes(x = x, y = y, fill = "walker_rdnbr")) +
  geom_sf(alpha = 0.1, color = "gray") +
  theme_bw()

## save perimeters in format usable by GEE script (Parks et al. 2019)
gee_perims <- full_perims[, c("YEAR_", "FIRE_NAME")]
colnames(gee_perims)  <- c("Fire_Year", "Fire_Name", "Shape")
gee_perims[, "Start_Day"] <- 152 ## see table 2 in Parks et al.
gee_perims[, "End_Day"] <- 305 ## need to use late date for sugar and north complex
for (i in 1:nrow(gee_perims)) {
  gee_perims$Fire_ID[i] <- paste0(gee_perims$Fire_Name[i], "_", gee_perims$Fire_Year[i])
}
st_write(gee_perims, "data/perimeters/gee_perimeters.shp", append = FALSE)

##---------------------------------------------------------------
## get severity data from standard perimeters, no daily data
##---------------------------------------------------------------

pull_full_severity <- function(fire_name) {
  rdnbr <- raster(paste0("data/rdnbr/", fire_name, "_rdnbr.tif"))
  perim_data <- readRDS(paste0("data/perimeters/", fire_name, "_perim.rds"))
  perim_data <- st_transform(perim_data, crs = crs(rdnbr))
  mask <- fasterize(perim_data, rdnbr)
  sev <- mask(rdnbr, mask)
  sev <- as.data.frame(as(sev, "SpatialPixelsDataFrame"))
  colnames(sev)[1] <- "rdnbr"
  return(sev)
}

dixie_severity <- pull_full_severity("dixie")
northcomplex_severity <- pull_full_severity("northcomplex")
sheep_severity <- pull_full_severity("sheep")
walker_severity <- pull_full_severity("walker")
sugar_severity <- pull_full_severity("sugar")

saveRDS(dixie_severity, "data/full_severity/dixie_severity.rds")
saveRDS(northcomplex_severity, "data/full_severity/northcomplex_severity.rds")
saveRDS(sheep_severity, "data/full_severity/sheep_severity.rds")
saveRDS(walker_severity, "data/full_severity/walker_severity.rds")
saveRDS(sugar_severity, "data/full_severity/sugar_severity.rds")
