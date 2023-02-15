
##---------------------------------------------------------------
## canopy_cover_cln.r
## SCRIPT TO CLEAN CANOPY COVER DATA FROM LIDAR
## BY: JACOB LEVINE - jacoblevine@princeton.edu
## 11/26/22
##---------------------------------------------------------------

##---------------------------------------------------------------
## 0. libraries
##---------------------------------------------------------------

library(rgdal)
library(sf)
library(raster)
library(ggplot2)
library(data.table)
library(FNN)
library(foreach)
library(doParallel)
library(doSNOW)
library(progress)

##---------------------------------------------------------------
## 1. point data
##---------------------------------------------------------------

point_test <- st_read("data/lidar/points/TAO__-42000_210000__t4_0p75_lp3__HighPoints.shp")

filelist <- list.files("data/lidar/points", ".shp")
filelist <- filelist[!grepl(".xml", filelist)]

poslist <- character()
gridpts <- data.frame(x = numeric(length(filelist)), y = numeric(length(filelist)))
for (i in 1:length(filelist)) {
  poslist[i] <- strsplit(filelist[i], "TAO__")[[1]][2]
}
for (i in 1:length(filelist)) {
  poslist[i] <- strsplit(poslist[i], "__t")[[1]][1]
}
for (i in 1:length(filelist)) {
  coord <- strsplit(poslist[i], "_")[[1]]
  gridpts[i,1] <- as.numeric(coord[1])
  gridpts[i,2] <- as.numeric(coord[2])
}

gridpts <- gridpts[order(gridpts$x, gridpts$y),]

grid <- st_as_sf(gridpts, coords = c("x", "y"), crs = crs(point_test))

template <- readRDS("data/templates/raster_template.rds")
template <- projectRaster(template, res = c(30,30), crs = crs(point_test))
targetcoords <- as.data.frame(template, xy = TRUE, na.rm = TRUE)
targetcoords <- as.matrix(targetcoords[,1:2])

create_filename <- function(x, y) {
 return(paste0("data/lidar/points/TAO__", format(x, scientific = FALSE), "_",
               format(y, scientific = FALSE), "__t4_0p75_lp3__HighPoints.shp"))
}

gridpts <- gridpts[gridpts$x > min(targetcoords[,1]) & gridpts$x < max(targetcoords[,1]) &
                   gridpts$y > min(targetcoords[,2]) & gridpts$y < max(targetcoords[,2]),]

subs <- split(unique(gridpts$x), ceiling(seq_along(unique(gridpts$x))/3))

gridpts_sub <- gridpts[gridpts[,1] %in% subs[[5]],]
panes <- st_read(create_filename(gridpts_sub[1,1], gridpts_sub[1,2]))
for (i in 2:nrow(gridpts_sub)) {
  newpane <- st_read(create_filename(gridpts_sub[i,1], gridpts_sub[i,2]))
  panes <- rbind(panes, newpane)
}

ext <- extent(panes)

## 2. get all target coords within pane
tg <- targetcoords[targetcoords[,1] < ext[2] & targetcoords[,1] > ext[1] &
                   targetcoords[,2] < ext[4] & targetcoords[,2] > ext[3],]


trees <- st_coordinates(panes)
t <- 1 ## for testing

cl <- makeSOCKcluster(4)
registerDoSNOW(cl)
iter <- nrow(tg)
pb <- progress_bar$new(format = "Calculating [:bar] :current/:total (:percent) | elapsed: :elapsed",
                       total = iter)
pb$tick(0)
progress <- function() pb$tick()
opts <- list(progress = progress)
tree_data <-
  foreach(t = 1:nrow(tg),
        .combine = rbind,
        .options.snow = opts) %dopar% {

            target <- tg[t,]
            target_group <- targetcoords[targetcoords[,1] >= target[1] - 30.0 & targetcoords[,1] <= target[1] + 30.0 &
                                         targetcoords[,2] >= target[2] - 30.0 & targetcoords[,2] <= target[2] + 30.0,]
            tree_clip <- trees[trees[,1] > target[1] - 45.0 & trees[,1] < target[1] + 45.0 &
                               trees[,2] > target[2] - 45.0 & trees[,2] < target[2] + 45.0, ]

            if (nrow(target_group) < 9 | nrow(tree_clip) == 0) {
              clust <- NA
              mean_dens <- NA
            }
            else {

              re <- 0.5 * sqrt(90*90 / nrow(tree_clip))
              r0 <- mean(knn.dist(tree_clip, k = 1))

              clust <- r0 / re
              mean_dens <- (nrow(tree_clip) / (90 * 90)) * 10000

            }

            data.frame(x = tg[t,1], y = tg[t,2],
                       clust = clust, mean_dens = mean_dens)
          }




test <- ext_gen(as.numeric(target))
x <- st_filter(panes, test)
nrow(x)
