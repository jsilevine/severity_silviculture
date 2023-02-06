library(rgdal)
library(rvest)
library(sf)
library(rmapshaper)
library(fasterize)
library(raster)
library(ggplot2)
library(gganimate)
setwd("../")

ogrListLayers("data/FACTS/FACTS.gdb/")
silv <- st_read("data/FACTS/FACTS.gdb/", layer = "R5_Activity_SilvReforestation")
silv <- silv[silv$ADMIN_FOREST_NAME == "Plumas National Forest",]
silv <- silv[!is.na(silv$ADMIN_FOREST_NAME),]


silv <- st_read("data/FACTS/FACTS.gdb/", layer = "VegBurnSeverity_updated2021")


sev <- readOGR("data/FACTS/FACTS.gdb/", layer = "VegBurnSeverity_updated2021")
