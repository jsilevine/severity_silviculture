
library(rvest)
library(sf)
library(rmapshaper)
library(fasterize)
library(terra)
library(ggplot2)

source("code/utility/utility_functions.r")

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

## save perimeters in format usable by GEE script (Parks et al. 2019)
gee_perims <- full_perims[, c("YEAR_", "FIRE_NAME")]
colnames(gee_perims)  <- c("Fire_Year", "Fire_Name", "Shape")
gee_perims[, "Start_Day"] <- 152 ## see table 2 in Parks et al.
gee_perims[, "End_Day"] <- 305 ## need to use late date for sugar and north complex
for (i in 1:nrow(gee_perims)) {
  gee_perims$Fire_ID[i] <- paste0(gee_perims$Fire_Name[i], "_", gee_perims$Fire_Year[i])
}
st_write(gee_perims, "data/perimeters/gee_perimeters.shp", append = FALSE)
