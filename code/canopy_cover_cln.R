
##---------------------------------------------------------------
## canopy_cover_cln.r
## SCRIPT TO CLEAN CANOPY COVER DATA FROM LIDAR
## BY: JACOB LEVINE - jacoblevine@princeton.edu
## 11/26/22
##---------------------------------------------------------------

##---------------------------------------------------------------
## 0.
##---------------------------------------------------------------
setwd("../")

library(rgdal)
library(sf)
library(raster)
library(rayshader)
library(ggplot2)
library(data.table)

##---------------------------------------------------------------
## 1. load canopy cover data
##---------------------------------------------------------------

cc_2.8 <- raster("data/canopy_cover/CanopyCover_2-8_PNF_1_BS30.tif")
cc_8.16 <- raster("data/canopy_cover/CanopyCover_8-16_PNF_1_BS30.tif")
cc_16.32 <- raster("data/canopy_cover/CanopyCover_16-32_PNF_1_BS30.tif")
cc_32.n <- raster("data/canopy_cover/CanopyCover_32-None_PNF_1_BS30.tif")
cc_full <- raster("data/canopy_cover/CanopyCover_2-None_PNF_1_BS30.tif")

##---------------------------------------------------------------
## 2. create raster template
##---------------------------------------------------------------

## resample to daily progression grid
togrid <- readRDS("data/fire_progression//dixie_daily_burned.rds")
togrid <- merge(togrid, readRDS("data/fire_progression//northcomplex_daily_burned.rds"))
togrid <- merge(togrid, readRDS("data/fire_progression//sheep_daily_burned.rds"))
togrid <- merge(togrid, readRDS("data/fire_progression//walker_daily_burned.rds"))
togrid <- merge(togrid, readRDS("data/fire_progression//sugar_daily_burned.rds"))

## snap to togrid, create template raster
cc_full <- projectRaster(cc_full, crs = crs(togrid))
cc_full <- resample(cc_full, togrid, "bilinear")
ext <- extent(togrid)
cc_full <- crop(cc_full, ext)
cc_full <- mask(cc_full, togrid)
togrid <- mask(togrid, cc_full)
togrid[!is.na(togrid[])] <- 1

## save template raster, for use in other cleaning scripts.
saveRDS(togrid, "data/templates/raster_template.rds")

cc_full.df <- as.data.frame(cc_full, xy = TRUE, na.rm = TRUE)
colnames(cc_full.df)[3] <- "cc"

saveRDS(cc_full, "data/canopy_cover/processed/cc_full.rds")
write.csv(cc_full.df, "data/canopy_cover/processed/cc_full.csv")

##---------------------------------------------------------------
## 3. use raster template to mask and resample cc data
##---------------------------------------------------------------

raster_template <- readRDS("data/templates/raster_template.rds")

## define function to snap to raster
snap_to_template <- function(raster, template) {

  if (as.character(crs(raster)) != as.character(crs(template))) {
    raster <- projectRaster(raster, crs = crs(template))
  }

  raster <- resample(raster, template, "bilinear")
  ext <- extent(template)
  raster <- crop(raster, ext)
  raster <- mask(raster, template)

  return(raster)

}

cc_2.8 <- snap_to_template(cc_2.8, raster_template)
cc_2.8.df <- as.data.frame(cc_2.8, xy = TRUE, na.rm = TRUE)
colnames(cc_2.8.df)[3] <- "cc"

saveRDS(cc_2.8, "data/canopy_cover/processed/cc_2-8.rds")
write.csv(cc_2.8.df, "data/canopy_cover/processed/cc_2-8.csv")

cc_8.16 <- snap_to_template(cc_8.16, raster_template)
cc_8.16.df <- as.data.frame(cc_8.16, xy = TRUE, na.rm = TRUE)
colnames(cc_8.16.df)[3] <- "cc"

saveRDS(cc_8.16, "data/canopy_cover/processed/cc_8-16.rds")
write.csv(cc_8.16.df, "data/canopy_cover/processed/cc_8-16.csv")

cc_16.32 <- snap_to_template(cc_16.32, raster_template)
cc_16.32.df <- as.data.frame(cc_16.32, xy = TRUE, na.rm = TRUE)
colnames(cc_16.32.df)[3] <- "cc"

saveRDS(cc_16.32, "data/canopy_cover/processed/cc_16-32.rds")
write.csv(cc_16.32.df, "data/canopy_cover/processed/cc_16-32.csv")

cc_32.n <- snap_to_template(cc_32.n, raster_template)
cc_32.n.df <- as.data.frame(cc_32.n, xy = TRUE, na.rm = TRUE)
colnames(cc_32.n.df)[3] <- "cc"

saveRDS(cc_32.n, "data/canopy_cover/processed/cc_32-n.rds")
write.csv(cc_32.n.df, "data/canopy_cover/processed/cc_32-n.csv")

##---------------------------------------------------------------
## 4. create complete dataset
##---------------------------------------------------------------

cc_2.8.df <- data.table(cc_2.8.df)
colnames(cc_2.8.df)[3] <- "cc2_8"
cc_8.16.df <- data.table(cc_8.16.df)
colnames(cc_8.16.df)[3] <- "cc8_16"
cc_16.32.df <- data.table(cc_16.32.df)
colnames(cc_16.32.df)[3] <- "cc16_32"
cc_32.n.df <- data.table(cc_32.n.df)
colnames(cc_32.n.df)[3] <- "cc32_n"



cc_complete <- cc_2.8.df[cc_8.16.df, on = .(x,y), nomatch = NULL]
cc_complete <- cc_complete[cc_16.32.df, on = .(x,y), nomatch = NULL]
cc_complete <- cc_complete[cc_32.n.df, on = .(x,y), nomatch = NULL]

fwrite(cc_complete, "data/canopy_cover/processed/cc_complete.csv")
