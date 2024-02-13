##---------------------------------------------------------------
## Script to generate daily fire progression maps from SUOMI VIIRS
##
## Author: Jacob Levine; jacoblevine@princeton.edu -- 310.754.6059
## This script is a draft, use at own risk (tests not implemented)
##---------------------------------------------------------------

## prerequisites
library(rgdal)
library(rvest)
library(sf)
library(rmapshaper)
library(fasterize)
library(raster)
library(ggplot2)
library(chron)
library(concaveman)
library(sp)
library(geosphere)
library(parallel)

source("function_def.r")

## TO USE THIS SCRIPT YOU WILL FIRST NEED TO FIRST REQUEST A MAPKEY FROM FIRMS
## TO DO SO GO TO BOTTOM OF FOLLOWING LINK
## https://firms.modaps.eosdis.nasa.gov/api/area/
mapkey <- "INSERT_YOUR_MAPKEY_HERE"

dixie_start <- as.Date("2021/07/13")
dixie_end <- as.Date("2021/10/24")
dixie_perim <- readRDS("dixie_perim.rds")
raster_template <- raster("data/rdnbr/dixie_rdnbr.tif")

viirs_data <- download_viirs_data(dixie_perim, dixie_start, dixie_end, mapkey = mapkey)

## this can take awhile depending on the size of the fire
## you might have to play around with the parameters in the model, particularly
## the cluster_dist paraterms and cluster_tol
dixie_daily_burned <- extract_daily_perims(start_date = dixie_start, ## must be in yyy/mm/dd format
                                           end_date = dixie_end,
                                           perim = dixie_perim, raster_template = raster_template,
                                           viirs = viirs_data, cluster_tol = 4,
                                           cluster_dist = 1125, large_cluster_dist = 20000)
plot(dixie_daily_burned)
