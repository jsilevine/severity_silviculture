
##---------------------------------------------------------------
## canopy_cover_cln.r
## SCRIPT TO CLEAN CANOPY COVER DATA FROM LIDAR
## BY: JACOB LEVINE - jacoblevine@princeton.edu
## 11/26/22
##---------------------------------------------------------------

## has to be run after delineate forest

##---------------------------------------------------------------
## 0.
##---------------------------------------------------------------

library(sf)
library(terra)
library(ggplot2)
library(data.table)

source("code/utility/utility_functions.r")

##---------------------------------------------------------------
## 1. load canopy cover data
##---------------------------------------------------------------

cc_2.8 <- rast("data/canopy_cover/CanopyCover_2-8_PNF_1_BS30.tif")
cc_8.16 <- rast("data/canopy_cover/CanopyCover_8-16_PNF_1_BS30.tif")
cc_16.32 <- rast("data/canopy_cover/CanopyCover_16-32_PNF_1_BS30.tif")
cc_32.n <- rast("data/canopy_cover/CanopyCover_32-None_PNF_1_BS30.tif")
cc_full <- rast("data/canopy_cover/CanopyCover_2-None_PNF_1_BS30.tif")

##---------------------------------------------------------------
## 2. create raster template
##---------------------------------------------------------------

## resample to template_grid
template <- rast("data/templates/isforest.tif")

cc_full <- snap_to_template(cc_full, template)
cc_full.df <- as.data.frame(cc_full, xy = TRUE, na.rm = TRUE)
colnames(cc_full.df)[3] <- "cc"

writeRaster(cc_full, "data/canopy_cover/processed/cc_full.tif", overwrite = TRUE)
write.csv(cc_full.df, "data/canopy_cover/processed/cc_full.csv")

cc_2.8 <- snap_to_template(cc_2.8, template)
cc_2.8.df <- as.data.frame(cc_2.8, xy = TRUE, na.rm = TRUE)
colnames(cc_2.8.df)[3] <- "cc"

writeRaster(cc_2.8, "data/canopy_cover/processed/cc_2-8.tif", overwrite = TRUE)
write.csv(cc_2.8.df, "data/canopy_cover/processed/cc_2-8.csv")

cc_8.16 <- snap_to_template(cc_8.16, template)
cc_8.16.df <- as.data.frame(cc_8.16, xy = TRUE, na.rm = TRUE)
colnames(cc_8.16.df)[3] <- "cc"

writeRaster(cc_8.16, "data/canopy_cover/processed/cc_8-16.tif", overwrite = TRUE)
write.csv(cc_8.16.df, "data/canopy_cover/processed/cc_8-16.csv")

cc_16.32 <- snap_to_template(cc_16.32, template)
cc_16.32.df <- as.data.frame(cc_16.32, xy = TRUE, na.rm = TRUE)
colnames(cc_16.32.df)[3] <- "cc"

writeRaster(cc_16.32, "data/canopy_cover/processed/cc_16-32.tif", overwrite = TRUE)
write.csv(cc_16.32.df, "data/canopy_cover/processed/cc_16-32.csv")

cc_32.n <- snap_to_template(cc_32.n, template)
cc_32.n.df <- as.data.frame(cc_32.n, xy = TRUE, na.rm = TRUE)
colnames(cc_32.n.df)[3] <- "cc"

writeRaster(cc_32.n, "data/canopy_cover/processed/cc_32-n.tif", overwrite = TRUE)
write.csv(cc_32.n.df, "data/canopy_cover/processed/cc_32-n.csv")

##---------------------------------------------------------------
## 4. create complete dataset
##---------------------------------------------------------------

names(cc_2.8) <- "cc_2.8"
names(cc_8.16) <- "cc_8.16"
names(cc_16.32) <- "cc_16.32"
names(cc_32.n) <- "cc_32.n"
names(cc_full) <- "cc_full"

cc_complete <- c(cc_2.8, cc_8.16, cc_16.32, cc_32.n, cc_full)

plot(cc_complete)

cc_complete_df <- as.data.frame(cc_complete, xy = TRUE, na.rm = TRUE)
write.csv(cc_complete_df, "data/canopy_cover/processed/cc_complete.csv", row.names = FALSE)

writeRaster(cc_complete, "data/canopy_cover/processed/cc_complete.tif")
