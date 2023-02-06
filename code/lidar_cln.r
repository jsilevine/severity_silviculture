
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
library(FNN)

##---------------------------------------------------------------
## 1. point data
##---------------------------------------------------------------

point_test <- st_read("data/lidar/points/TAO__-4000_210000__t4_0p75_lp3__HighPoints.shp")
plot(point_test[,1])
head(point_test)

xy <- st_coordinates(point_test)
max(xy[,1])
min(xy[,1])
max(xy[,2])

ext <- extent(point_test)
area <- (ext[2] - ext[1]) * (
  ext[4] - ext[3])
re <- 0.5 * sqrt(area / nrow(point_test))

r0 <- mean(knn.dist(st_coordinates(point_test), k = 1))

r0 / re

plt <- xy[xy[,1] < min(xy[,1]) + 30,]
plot(plt[,1], plt[,2])
abline(max(xy[,2]) - 90, 0)

crs(point_test)

filelist <- list.files("data/lidar/points", ".shp")
filelist <- filelist[!grepl(".xml", filelist)]

poslist <- character()
gridpts <- data.frame(x = numeric(length(filelist)), y = numeric(length(filelist)))
for (i in 1:length(filelist)) {
  poslist[i] <- strsplit(filelist[i], "-")[[1]][2]
}
for (i in 1:length(filelist)) {
  poslist[i] <- strsplit(poslist[i], "__t")[[1]][1]
}
for (i in 1:length(filelist)) {
  coord <- strsplit(poslist[i], "_")[[1]]
  gridpts[i,1] <- as.numeric(coord[1])
  gridpts[i,2] <- as.numeric(coord[2])
}

gridpts <- gridpts[order(gridpts[,1], gridpts[,2]),]
gridpts


grid <- st_as_sf(gridpts, coords = c("x", "y"), crs = crs(point_test))
plot(grid)

template <- readRDS("data/templates/raster_template.rds")
template <- projectRaster(template, res = c(30,30), crs = crs(point_test))
targetcoords <- as.data.frame(template, xy = TRUE, na.rm = TRUE)
targetcoords <- as.matrix(targetcoords[,1:2])
crs(point_test)

crs(template)

## 1. load pane
create_filename <- function(x, y) {
 return(paste0("data/lidar/points/TAO__-", as.character(x), "_", as.character(y), "__t4_0p75_lp3__HighPoints.shp"))
}

as.character(y)
x = 4000
y = 10000

pane <- st_read(create_filename(gridpts[1,1], gridpts[1,2]))
ext <- extent(pane)

## 2. get all target coords within pane
tg <- targetcoords[targetcoords[,1] < ext[2] & targetcoords[,1] > ext[1] &
                   targetcoords[,2] < ext[4] & targetcoords[,2] > ext[3],]

## 3. also load adjacent panes




##---------------------------------------------------------------
## see if pasting together is feasible
##---------------------------------------------------------------



nrow(gridpts) / 50

subs <- split(unique(gridpts$x), ceiling(seq_along(unique(gridpts$x))/6))
subs[[1]]

gridpts_sub <- gridpts[gridpts[,1] %in% subs[[1]],]
panes <- st_read(create_filename(gridpts_sub[1,1], gridpts_sub[1,2]))
for (i in 2:nrow(gridpts_sub))) {
  newpane <- st_read(create_filename(gridpts_sub[i,1], gridpts_sub[i,2]))
  panes <- rbind(panes, newpane)
}


panes <- pane
for (i in 2:50) {
  newpane <- st_read(paste0("data/lidar/points/", filelist[i]))
  panes <- rbind(panes, newpane)
}

plot(panes[,1], pct = ".")

filelist
