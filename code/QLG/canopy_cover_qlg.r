
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

## resample to daily progression grid
tp <- readRDS("data/templates/raster_template.rds")

qlg_data <- read.csv("data/qlg_lidar_data.csv", row.names = 1)
togrid <- rasterFromXYZ(read.csv("data/qlg_lidar_data.csv", row.names = 1)[,1:3], crs = CRS('+init=EPSG:3310'))
length(unique(qlg_data$x))

rl <- list()
for (i in 3:(ncol(qlg_data))) {
  rl[[i-2]] <- rasterFromXYZ(qlg_data[,c(1,2,i)], crs = CRS('+init=EPSG:3310'))
}
rs <- stack(rl[[1]], rl[[2]], rl[[3]], rl[[4]], rl[[5]], rl[[6]], rl[[7]])
rs <- projectRaster(rs, crs = crs(tp))
togrid <- projectRaster(rs, crs = crs(tp))

## snap to togrid, create template raster
cc_full <- projectRaster(cc_full, crs = crs(togrid))
cc_full <- resample(cc_full, togrid, "bilinear")
ext <- extent(togrid)
cc_full <- crop(cc_full, ext)
cc_full <- mask(cc_full, togrid)
togrid <- mask(togrid, cc_full)

saveRDS(togrid, "data/templates/qlg_template.rds")

##---------------------------------------------------------------
## 3. use raster template to mask and resample cc data
##---------------------------------------------------------------

raster_template <- readRDS("data/templates/qlg_template.rds")

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
full_data <- stack(raster_template, subset(cc_2.8,1))

cc_8.16 <- snap_to_template(cc_8.16, raster_template)
full_data <- stack(full_data, subset(cc_8.16,1))

cc_16.32 <- snap_to_template(cc_16.32, raster_template)
full_data <- stack(full_data, subset(cc_16.32,1))

cc_32.n <- snap_to_template(cc_32.n, raster_template)
full_data <- stack(full_data, subset(cc_32.n,1))

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

##---------------------------------------------------------------
## For QLG
##---------------------------------------------------------------

qlg_data <- as.data.frame(full_data, na.rm = TRUE, xy = TRUE)
qlg_data <- data.table(qlg_data)
colnames(qlg_data)[10:13] <- c("cc2_8", "cc8_16", "cc16_32", "cc32_n")

## calculate evenness as measure of ladder fuels.
calc_em <- function(i, data = qlg_data) {
  cc <- as.numeric(data[i, .(cc2_8, cc8_16, cc16_32, cc32_n)])
  cc <- cc[cc > 0.0]
  cc_p <- cc / sum(cc)
  iD <- sum(cc_p^2)
  D <- 1 / iD
  E <- D / length(cc)
  em <- E * mean(cc)
}

i <- 1:nrow(qlg_data)
qlg_data[,"em"] <- sapply(i, FUN = calc_em) ## takes a little while
qlg_data <- qlg_data[, c(1:9, 14)]
colnames(qlg_data)[10] <- "ladder_fuels"
fwrite(qlg_data, "data/qlg_data_full.csv")
