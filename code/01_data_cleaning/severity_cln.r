
library(rgdal)
library(rvest)
library(sf)
library(rmapshaper)
library(fasterize)
library(raster)
library(ggplot2)
library(gganimate)
setwd("severity_and_silviculture/")

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


##---------------------------------------------------------------
## Severity
##---------------------------------------------------------------

template <- rast("data/templates/isforest.tif")

sevfiles <- list.files("data/cbi")

for (i in 1:length(sevfiles)) {

  fire_name <- strsplit(tolower(strsplit(sevfiles[i], "_CBI_")[[1]][1]), "_")[[1]][1]

  severity <- rast(paste0("data/cbi/", sevfiles[i]))
  severity <- snap_to_template(severity, template)
  severity_df <- as.data.frame(severity, xy = TRUE, na.rm = TRUE)
  writeRaster(severity, paste0("data/cbi/", fire_name, "_processed.tif"), overwrite = TRUE)
  if (i == 1) {
    severity_full_df <- severity_df
  } else {
    severity_full_df <- rbind(severity_full_df, severity_df)
  }
  write.csv(severity_df, paste0("data/cbi/", fire_name, "_processed.csv"))

}







##---------------------------------------------------------------
## OLD/SCRATCH
##---------------------------------------------------------------

##---------------------------------------------------------------
## get severity data from standard perimeters, no daily data
##---------------------------------------------------------------

## fire_name <- "dixie"
## pull_full_severity <- function(fire_name) {
##   perim_data <- readRDS(paste0("data/perimeters/", fire_name, "_perim.rds"))
##   cbi <- raster(paste0("data/cbi/ravg/ravg_", perim_data$YEAR_[1], "_cbi4.tif"))
##   perim_data <- st_transform(perim_data, crs = crs(rdnbr))
##   rdnbr <- crop(rdnbr, perim_data)
##   mask <- fasterize(perim_data, rdnbr)
##   sev <- mask(rdnbr, mask)
##   sev[sev > 4]  <- NA
##   return(sev)
## }


## dixie_severity <- pull_full_severity("dixie")
## northcomplex_severity <- pull_full_severity("northcomplex")
## sheep_severity <- pull_full_severity("sheep")
## walker_severity <- pull_full_severity("walker")
## sugar_severity <- pull_full_severity("sugar")

## saveRDS(dixie_severity, "data/cbi/ravg/dixie_cbi_ravg.rds")
## saveRDS(northcomplex_severity, "data/cbi/ravg/northcomplex_cbi_ravg.rds")
## saveRDS(sheep_severity, "data/cbi/ravg/sheep_cbi_ravg.rds")
## saveRDS(walker_severity, "data/cbi/ravg/walker_cbi_ravg.rds")
## saveRDS(sugar_severity, "data/cbi/ravg/sugar_cbi_ravg.rds")
