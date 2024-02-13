
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

setwd("data/lidar_processing")

##---------------------------------------------------------------
## 1. point data
##---------------------------------------------------------------
point_test <- st_read("data/lidar/points/TAO__-42000_210000__t4_0p75_lp3__HighPoints.shp")

plot(point_test['height_ft'])
ggplot(point_test, aes(color = height_ft)) +
  geom_sf() +
  scale_color_gradient(low = "#e7e1ef", high = "#dd1c77") +
  theme_bw() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))
ggsave("../f1.pdf")

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
fullcoords <- as.data.frame(template, xy = TRUE, na.rm = FALSE)
fullcoords <- as.matrix(fullcoords[,1:2])

create_filename <- function(x, y) {
 return(paste0("lidar/points/TAO__", format(x, scientific = FALSE), "_",
               format(y, scientific = FALSE), "__t4_0p75_lp3__HighPoints.shp"))
}

gridpts <- gridpts[gridpts$x > min(targetcoords[,1]) & gridpts$x < max(targetcoords[,1]) &
                   gridpts$y > min(targetcoords[,2]) & gridpts$y < max(targetcoords[,2]),]

subs <- split(unique(gridpts$x), ceiling(seq_along(unique(gridpts$x))/2))

for (i in 2:length(subs)) {
  subs[[i]] <- c(subs[[i-1]][length(subs[[i-1]])], subs[[i]])
}

cl <- makeSOCKcluster(15)
registerDoSNOW(cl)

for (s in 3:length(subs)) {
  
  print(paste0("Starting grid subset ", s, " of ", length(subs)))
  
  gridpts_sub <- gridpts[gridpts[,1] %in% subs[[s]],]
  panes <- st_read(create_filename(gridpts_sub[1,1], gridpts_sub[1,2]))
  for (i in 2:nrow(gridpts_sub)) {
    newpane <- st_read(create_filename(gridpts_sub[i,1], gridpts_sub[i,2]))
    panes <- rbind(panes, newpane) 
  }

  ext <- extent(min(gridpts_sub[,1]), max(gridpts_sub[, 1])+2000, 
                min(gridpts_sub[,2]), max(gridpts_sub[, 2])+2000)

  if (s == 1) {
    ext_sub <- ext
  } else {
   ext_sub <- extent(unique(gridpts_sub[,1])[2]-90, max(gridpts_sub[, 1])+2000, 
                      min(gridpts_sub[,2]), max(gridpts_sub[, 2])+2000)
  }

  ## 2. get all target coords within pane
  fg <- fullcoords[fullcoords[,1] < ext[2] & fullcoords[,1] > ext[1] &
                     fullcoords[,2] < ext[4] & fullcoords[,2] > ext[3],]

  tg <- targetcoords[targetcoords[,1] < ext_sub[2] & targetcoords[,1] > ext_sub[1] &
                    targetcoords[,2] < ext_sub[4] & targetcoords[,2] > ext_sub[3],]

  trees <- st_coordinates(panes)

  iter <- nrow(tg)
  pb <- progress_bar$new(format = "Calculating [:bar] :current/:total (:percent) | elapsed: :elapsed",
                         total = iter)
  pb$tick(0)
  progress <- function() pb$tick()
  opts <- list(progress = progress)
  tree_data <-
   foreach(t = 1:nrow(tg),
        .combine = rbind,
        .packages = c("FNN"),
        .options.snow = opts) %dopar% {

            target <- tg[t,]
            target_group <- fg[fg[,1] >= target[1] - 30.0 & fg[,1] <= target[1] + 30.0 &
                               fg[,2] >= target[2] - 30.0 & fg[,2] <= target[2] + 30.0,]
            tree_clip <- trees[trees[,1] > target[1] - 45.0 & trees[,1] < target[1] + 45.0 &
                               trees[,2] > target[2] - 45.0 & trees[,2] < target[2] + 45.0, ]

            if (nrow(target_group) < 9 | nrow(matrix(tree_clip, ncol = 2)) <= 1) {
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

  write.csv(tree_data, paste0("processed/tree_data_", s, ".csv"))

}


clust_data <- data.frame(x = numeric(0), y = numeric(0), z = numeric(0))
for (i in 1:31) {

  nt <- read.csv(paste0("data/processed/tree_data_", i, ".csv"))
  clust_data <- rbind(clust_data, data.frame(x = nt$x, y = nt$y, z = nt$clust))

}

clust_raster <- rasterFromXYZ(clust_data)
plot(clust_raster)


dens_data <- data.frame(x = numeric(0), y = numeric(0), z = numeric(0))
for (i in 1:31) {

  nt <- read.csv(paste0("data/processed/tree_data_", i, ".csv"))
  dens_data <- rbind(dens_data, data.frame(x = nt$x, y = nt$y, z = nt$mean_dens))

}

dens_raster <- rasterFromXYZ(dens_data)
plot(dens_raster)





##---------------------------------------------------------------
## 2. basin data
##---------------------------------------------------------------

basin_test <- st_read("data/lidar/basins/TAO__-42000_210000__t4_0p75_lp3__Basin_Map.shp")

ggplot(basin_test, aes(fill = height_ft)) +
  geom_sf(color = NA) +
  scale_fill_gradient(low = "#fa9fb5", high = "#ae017e") +
  theme_bw() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))
ggsave("../f2.pdf")



bbox <- st_bbox(basin_test)
bbox[1] <- -42000
bbox[2] <- 210000
bbox[3] <- -40000
bbox[4] <- 212000
bbox <- st_as_sfc(bbox)

ggplot(basin_test) + 
  geom_sf(aes(fill = height)) + 
  geom_sf(data = bbox, alpha = 0.1)

dif_test <- st_difference(bbox, st_union(basin_test))
dif_test <- st_cast(dif_test, "POLYGON")
dif_test <- st_as_sf(dif_test)
dif_test$area <- st_area(dif_test)
plot(dif_test)

test <- st_intersection(basin_test, bbox)


filelist <- list.files("data/lidar/basins", ".shp")
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

grid <- st_as_sf(gridpts, coords = c("x", "y"), crs = crs(basin_test))

create_filename <- function(x, y) {
  return(paste0("data/lidar/basins/TAO__", format(x, scientific = FALSE), "_",
                format(y, scientific = FALSE), "__t4_0p75_lp3__Basin_Map.shp"))
}
create_filename2 <- function(x, y) {
  return(paste0("data/lidar/gaps/TAO__", format(x, scientific = FALSE), "_",
                format(y, scientific = FALSE), "__t4_0p75_lp3__GAP.shp"))
}

for (i in 1:nrow(gridpts)) {
  
  b <- st_read(create_filename(gridpts[i,"x"], gridpts[i,"y"]))
  bbox <- st_bbox(b)
  bbox[1] <- gridpts[i,"x"]
  bbox[2] <- gridpts[i,"y"]
  bbox[3] <- gridpts[i,"x"] + 2000
  bbox[4] <- gridpts[i,"y"] + 2000
  bbox <- st_as_sfc(bbox)
  
  dif <- st_difference(bbox, st_union(b))
  dif <- st_cast(dif, "POLYGON")
  dif <- st_as_sf(dif)
  st_write(dif, create_filename2(gridpts[i,"x"], gridpts[i,"y"]))
  
}

calc_frac <- function(p, a) {
  (2 * log(0.25 * p)) / log(a)
}

template <- readRDS("data/templates/raster_template.rds")
template <- projectRaster(template, res = c(30,30), crs = crs(basin_test))
targetcoords <- as.data.frame(template, xy = TRUE, na.rm = TRUE)
targetcoords <- as.matrix(targetcoords[,1:2])
fullcoords <- as.data.frame(template, xy = TRUE, na.rm = FALSE)
fullcoords <- as.matrix(fullcoords[,1:2])

gridpts <- gridpts[gridpts$x > min(targetcoords[,1]) & gridpts$x < max(targetcoords[,1]) &
                     gridpts$y > min(targetcoords[,2]) & gridpts$y < max(targetcoords[,2]),]

subs <- split(unique(gridpts$x), ceiling(seq_along(unique(gridpts$x))/2))

for (i in 2:length(subs)) {
  subs[[i]] <- c(subs[[i-1]][length(subs[[i-1]])], subs[[i]])
}

cl <- makeSOCKcluster(15)
registerDoSNOW(cl)

for (s in 1:length(subs)) {
  
  print(paste0("Starting grid subset ", s, " of ", length(subs)))
  
  gridpts_sub <- gridpts[gridpts[,1] %in% subs[[s]],]
  panes <- st_read(create_filename2(gridpts_sub[1,1], gridpts_sub[1,2]))
  for (i in 2:nrow(gridpts_sub)) {
    newpane <- st_read(create_filename2(gridpts_sub[i,1], gridpts_sub[i,2]))
    panes <- rbind(panes, newpane) 
  }

  ggplot(panes[1000:10000,]) +
    geom_sf(aes(fill = FID))

  
  ext <- extent(min(gridpts_sub[,1]), max(gridpts_sub[, 1])+2000, 
                min(gridpts_sub[,2]), max(gridpts_sub[, 2])+2000)
  
  if (s == 1) {
    ext_sub <- ext
  } else {
    ext_sub <- extent(unique(gridpts_sub[,1])[2]-90, max(gridpts_sub[, 1])+2000, 
                      min(gridpts_sub[,2]), max(gridpts_sub[, 2])+2000)
  }
  
  ## 2. get all target coords within pane
  fg <- fullcoords[fullcoords[,1] < ext[2] & fullcoords[,1] > ext[1] &
                     fullcoords[,2] < ext[4] & fullcoords[,2] > ext[3],]
  
  tg <- targetcoords[targetcoords[,1] < ext_sub[2] & targetcoords[,1] > ext_sub[1] &
                       targetcoords[,2] < ext_sub[4] & targetcoords[,2] > ext_sub[3],]
  
  trees <- st_coordinates(panes)
  
  iter <- nrow(tg)
  pb <- progress_bar$new(format = "Calculating [:bar] :current/:total (:percent) | elapsed: :elapsed",
                         total = iter)
  pb$tick(0)
  progress <- function() pb$tick()
  opts <- list(progress = progress)
  tree_data <-
    foreach(t = 1:nrow(tg),
            .combine = rbind,
            .packages = c("FNN",
                          "sf"),
            .options.snow = opts) %dopar% {
              
              target <- tg[t,]
              target_group <- fg[fg[,1] >= target[1] - 30.0 & fg[,1] <= target[1] + 30.0 &
                                   fg[,2] >= target[2] - 30.0 & fg[,2] <= target[2] + 30.0,]
              bbox <- st_bbox(b)
              bbox[1] <- target[1] - 30.0
              bbox[2] <- target[2] - 30.0
              bbox[3] <- target[1] + 30.0
              bbox[4] <- target[2] + 30.0
              bbox <- st_as_sfc(bbox)
              
              clip <- st_intersection(panes, bbox)
              clip$area <- st_area(clip)
              
              if (nrow(target_group) < 9 | nrow(clip) < 1) {
                
                mean_area <- NA
                median_area <- NA
                sd_area <- NA
                mean_frac <- NA
                
              } else {
                
                mean_area <- mean(clip$area)
                median_area <- median(clip$area)
                sd_area <- sd(clip$area)
                
                clip$perim <- st_length(st_cast(clip, "MULTILINESTRING"))
                mean_frac <- mean(mapply(calc_frac, clip$perim, clip$area))
                
              }
              
              data.frame(x = tg[t,1], y = tg[t,2],
                         mean_area = mean_area, 
                         median_area = median_area,
                         sd_area = sd_area,
                         mean_frac = mean_frac)
              
            }
  
  write.csv(tree_data, paste0("processed/gap_data_", s, ".csv"))
  
}


check_gap <- read.csv("data/processed/gap_data_4.csv")

ras <- rasterFromXYZ(data.frame(x = check_gap$x, y = check_gap$y, z = check_gap$median_area))
plot(ras)

head(check_gap)


mean_gap_area_data <- data.frame(x = numeric(0), y = numeric(0), z = numeric(0))
for (i in 1:21) {

  nt <- read.csv(paste0("data/processed/gap_data_", i, ".csv"))
  mean_gap_area_data <- rbind(mean_gap_area_data, data.frame(x = nt$x, y = nt$y, z = nt$mean_area))

}

mean_gap_area_raster <- rasterFromXYZ(mean_gap_area_data)
plot(mean_gap_area_raster)


dens_data <- data.frame(x = numeric(0), y = numeric(0), z = numeric(0))
for (i in 1:31) {

  nt <- read.csv(paste0("data/processed/tree_data_", i, ".csv"))
  dens_data <- rbind(dens_data, data.frame(x = nt$x, y = nt$y, z = nt$mean_dens))

}

dens_raster <- rasterFromXYZ(dens_data)
plot(dens_raster)



##---------------------------------------------------------------
## 2. basin data -- try new approach
##---------------------------------------------------------------

basin_test <- st_read("data/lidar/basins/TAO__-42000_210000__t4_0p75_lp3__Basin_Map.shp")
bbox <- st_bbox(basin_test)
bbox[1] <- -42000
bbox[2] <- 210000
bbox[3] <- -40000
bbox[4] <- 212000
bbox <- st_as_sfc(bbox)

ggplot(basin_test) +
  geom_sf(aes(fill = height)) +
  geom_sf(data = bbox, alpha = 0.1)

dif_test <- st_difference(bbox, st_union(basin_test))
dif_test <- st_cast(dif_test, "POLYGON")
dif_test <- st_as_sf(dif_test)
dif_test$area <- st_area(dif_test)
plot(dif_test)

test <- st_intersection(basin_test, bbox)


filelist <- list.files("data/lidar/basins", ".shp")
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

grid <- st_as_sf(gridpts, coords = c("x", "y"), crs = crs(basin_test))

create_filename <- function(x, y) {
  return(paste0("data/lidar/basins/TAO__", format(x, scientific = FALSE), "_",
                format(y, scientific = FALSE), "__t4_0p75_lp3__Basin_Map.shp"))
}
create_filename2 <- function(x, y) {
  return(paste0("data/lidar/gaps/TAO__", format(x, scientific = FALSE), "_",
                format(y, scientific = FALSE), "__t4_0p75_lp3__GAP.shp"))
}

for (i in 1:nrow(gridpts)) {



  i <- 40
  gridpts_sub <- gridpts[i,]
  gridpts_sub <- rbind(gridpts_sub, gridpts[gridpts[,"x"] == gridpts[i,"x"] & abs(gridpts[,"y"] - gridpts[i,"y"]) == 2000,])
  gridpts_sub <- rbind(gridpts_sub, gridpts[gridpts[,"x"] == gridpts[i,"x"] - 2000 & abs(gridpts[,"y"] - gridpts[i,"y"]) == 2000,])
  gridpts_sub <- rbind(gridpts_sub, gridpts[gridpts[,"x"] == gridpts[i,"x"] + 2000 & abs(gridpts[,"y"] - gridpts[i,"y"]) == 2000,])

  panes <- st_read(create_filename(gridpts_sub[1,1], gridpts_sub[1,2]))
  cpane <- panes
  for (i in 2:nrow(gridpts_sub)) {
    newpane <- st_read(create_filename(gridpts_sub[i,1], gridpts_sub[i,2]))
    panes <- rbind(panes, newpane)
  }

  bbox <- st_bbox(cpane)
  bbo <- bbox ## make a copy
  bbox[1] <- gridpts_sub[1,"x"] - 50
  bbox[2] <- gridpts_sub[1,"y"] - 50
  bbox[3] <- gridpts_sub[1,"x"] + 2050
  bbox[4] <- gridpts_sub[1,"y"] + 2050

  crop <- st_crop(panes, bbox)
  bbox <- st_as_sfc(st_bbox(crop))
  dif <- st_difference(bbox, st_union(crop))
  dif <- st_cast(dif, "POLYGON")
  dif <- st_as_sf(dif)

  bbox <- st_bbox(cpane)
  bbox[1] <- gridpts_sub[1,"x"]
  bbox[2] <- gridpts_sub[1,"y"]
  bbox[3] <- gridpts_sub[1,"x"] + 2000
  bbox[4] <- gridpts_sub[1,"y"] + 2000
  bbox <- st_as_sfc(bbox)
  dif <- st_crop(dif, bbox)

    ggplot(dif) +
    geom_sf(fill = "gray")

    ggplot(crop) +
    geom_sf(aes(fill = height))

  +
    geom_sf(data = bbox, fill = "gray", alpha = 0.2)


  bbox <- st_bbox(panes_test)
  bbo <- bbox ## make a copy
  bbox[1] <- gridpts_sub[1,"x"]
  bbox[2] <- gridpts_sub[1,"y"]
  bbox[3] <- gridpts_sub[1,"x"] + 1990
  bbox[4] <- gridpts_sub[1,"y"] + 1990

  for (i in 1:4) {
    if (abs(bbox[i] - bbo[i]) > 50) {
      if (i <= 2) {
        bbox[i] <- bbo[i] + 10
      } else {
        bbox[i] <- bbo[i] - 10
      }
    }
  }
  bbox <- st_as_sfc(bbox)

  dif <- st_difference(bbox, st_union(panes_test))
  dif <- st_cast(dif, "POLYGON")
  dif <- st_as_sf(dif)

  ggplot(dif) +
    geom_sf(fill = "gray")





  b <- st_read(create_filename(gridpts[i,"x"], gridpts[i,"y"]))
  bbox <- st_bbox(b)
  bbox[1] <- gridpts[i,"x"]
  bbox[2] <- gridpts[i,"y"]
  bbox[3] <- gridpts[i,"x"] + 2000
  bbox[4] <- gridpts[i,"y"] + 2000
  bbox <- st_as_sfc(bbox)

  dif <- st_difference(bbox, st_union(b))
  dif <- st_cast(dif, "POLYGON")
  dif <- st_as_sf(dif)
  st_write(dif, create_filename2(gridpts[i,"x"], gridpts[i,"y"]))

}

calc_frac <- function(p, a) {
  (2 * log(0.25 * p)) / log(a)
}

template <- readRDS("data/templates/raster_template.rds")
template <- projectRaster(template, res = c(30,30), crs = crs(basin_test))
targetcoords <- as.data.frame(template, xy = TRUE, na.rm = TRUE)
targetcoords <- as.matrix(targetcoords[,1:2])
fullcoords <- as.data.frame(template, xy = TRUE, na.rm = FALSE)
fullcoords <- as.matrix(fullcoords[,1:2])

gridpts <- gridpts[gridpts$x > min(targetcoords[,1]) & gridpts$x < max(targetcoords[,1]) &
                     gridpts$y > min(targetcoords[,2]) & gridpts$y < max(targetcoords[,2]),]

subs <- split(unique(gridpts$x), ceiling(seq_along(unique(gridpts$x))/2))

for (i in 2:length(subs)) {
  subs[[i]] <- c(subs[[i-1]][length(subs[[i-1]])], subs[[i]])
}

cl <- makeSOCKcluster(15)
registerDoSNOW(cl)

## function to create a bounding box of irregular shape
generate_bound <- function(grd, crs = 4326) {

  bound <- grd[grd$x == min(grd$x), ]
  for (i in 1:length(unique(grd$x))) {
    x <- unique(grd$x)[i]
    if (i > 1) {

      if (abs(max(grd[grd$x == x, "y"]) - max(grd[grd$x == unique(grd$x)[i-1], "y"])) >= 2000) {

        if (max(grd[grd$x == x, "y"]) - max(grd[grd$x == unique(grd$x)[i-1], "y"]) > 0) {

          bound <- rbind(bound, grd[grd$x == x & grd$y > max(grd[grd$x == unique(grd$x)[i-1], "y"]), ])
          bound <- rbind(bound, data.frame(x = unique(grd$x)[i],
                                           y = max(grd[grd$x == unique(grd$x)[i], "y"]) + 2000))

        } else {

          bound <- rbind(bound, grd[grd$x == unique(grd$x)[i-1] & grd$y > max(grd[grd$x == x, "y"]), ])
          bound <- rbind(bound, data.frame(x = unique(grd$x)[i],
                                         y = max(grd[grd$x == unique(grd$x)[i], "y"]) + 2000))

        }

      } else {
        bound <- rbind(bound, data.frame(x = unique(grd$x)[i],
                                         y = max(grd[grd$x == unique(grd$x)[i], "y"]) + 2000))
      }
    } else {
      bound <- rbind(bound, data.frame(x = unique(grd$x)[i],
                                       y = max(grd[grd$x == unique(grd$x)[i], "y"]) + 2000))
    }

  }
  rs <- grd[grd$x == max(grd$x),]
  rs$x <- rs$x + 2000
  rs <- rbind(rs, data.frame(x = rs$x[1], y = max(rs$y) + 2000))
  bound <- rbind(bound, rs[order(rs$y, decreasing = TRUE),])
  for (i in length(unique(grd$x)):1) {
    x <- unique(grd$x)[i]
    if (i < length(unique(grd$x))) {
      if (abs(min(grd[grd$x == x, "y"]) - min(grd[grd$x == unique(grd$x)[i+1], "y"])) >= 2000) {

        if (min(grd[grd$x == x, "y"]) - min(grd[grd$x == unique(grd$x)[i+1], "y"]) < 0) {

          bound <- rbind(bound, grd[grd$x == x & grd$y <= min(grd[grd$x == unique(grd$x)[i+1], "y"]), ])

        } else {

          bound <- rbind(bound, grd[grd$x == unique(grd$x)[i+1] & grd$y <= min(grd[grd$x == x, "y"]), ])
          bound <- rbind(bound, data.frame(x = unique(grd$x)[i],
                                         y = min(grd[grd$x == unique(grd$x)[i], "y"])))

        }

      } else {
        bound <- rbind(bound, data.frame(x = unique(grd$x)[i],
                                         y = min(grd[grd$x == unique(grd$x)[i], "y"])))
      }
    } else {
      bound <- rbind(bound, data.frame(x = unique(grd$x)[i],
                                       y = min(grd[grd$x == unique(grd$x)[i], "y"])))
    }

  }

  #bound <- st_as_sf(bound, coords = c("x", "y"), crs = 4326)
  #bound <- st_combine(bound)
  #bound <- st_cast(bound, "POLYGON")
  #st_crs(bound) <- crs
  return(bound)

}

plot(st_combine(st_as_sf(gridpts, coords = c("x", "y"))))

nper <- 1
for (s in 1:round((length(unique(gridpts$x))/nper))) {

  print(paste0("Starting grid subset ", s, " of ", length(subs)))

  x_subs <- unique(gridpts$x)[((s-1)*nper+1):(((s-1)*nper)+nper)]
  x_subs_plus <- c(x_subs[1]-2000, x_subs, x_subs[nper]+2000)

  gridpts_sub <- gridpts[gridpts[,1] %in% x_subs_plus,]

  bound <- generate_bound(gridpts_sub, st_crs(basin_test))

  plot(bound)
  plot(bound$x, bound$y)
  plot(gridpts_sub$x, gridpts_sub$y)

  mask <- st_bbox(bound)
  mask[1] <- x_subs[1] - 150
  mask[2] <- min(gridpts_sub[gridpts_sub$x == x_subs[1], "y"])
  mask[3] <- x_subs[nper] + 2150
  mask[4] <- max(gridpts_sub[gridpts_sub$x == x_subs[1], "y"] + 2000)
  bbox <- mask ## save copy
  mask <- st_as_sfc(mask)

  bound_crop <- st_crop(bound, mask)
  plot(bound_crop)

  sqrt(2000^2 + 2000^2)


  x <- as.matrix(dist(gridpts[1:10,], upper = TRUE, diag = TRUE))
  w <- as.numeric(which(x[1,] < 2829))
  gridpts_sub <- gridpts[w,]
  panes <- st_read(create_filename(gridpts_sub[1,1], gridpts_sub[1,2]))
  for (i in 2:nrow(gridpts_sub)) {
    newpane <- st_read(create_filename(gridpts_sub[i,1], gridpts_sub[i,2]))
    panes <- rbind(panes, newpane)
  }

  ggplot(panes) +
    geom_sf(aes(fill = height))

  ch <- st_convex_hull(st_combine(panes))
  plot(ch)


  plot(gridpts[w,"x"], gridpts[w, "y"])


  panes_crop <- st_crop(panes, bound_crop)

  panes_test <- panes_crop[1:75000,]

  panes <- st_read(create_filename(gridpts_sub[1,1], gridpts_sub[1,2]))
  for (i in 2:nrow(gridpts_sub)) {
    newpane <- st_read(create_filename(gridpts_sub[i,1], gridpts_sub[i,2]))
    panes <- rbind(panes, newpane)
  }

  ch <- st_convex_hull(st_combine(basin_test))
  ggplot(basin_test) +
    geom_sf(aes(fill = height)) +
    geom_sf(data = ch, fill = "gray", alpha = 0.2)

  dif <- st_difference(bound_crop, st_union(panes))
  dif <- st_cast(dif, "POLYGON")
  dif <- st_as_sf(dif)

  ggplot(dif) +
    geom_sf(fill = "gray")

  ext <- extent(min(gridpts_sub[,1]), max(gridpts_sub[, 1])+2000,
                min(gridpts_sub[,2]), max(gridpts_sub[, 2])+2000)

  if (s == 1) {
    ext_sub <- ext
  } else {
    ext_sub <- extent(unique(gridpts_sub[,1])[2]-90, max(gridpts_sub[, 1])+2000,
                      min(gridpts_sub[,2]), max(gridpts_sub[, 2])+2000)
  }

  ## 2. get all target coords within pane
  fg <- fullcoords[fullcoords[,1] < ext[2] & fullcoords[,1] > ext[1] &
                     fullcoords[,2] < ext[4] & fullcoords[,2] > ext[3],]

  tg <- targetcoords[targetcoords[,1] < ext_sub[2] & targetcoords[,1] > ext_sub[1] &
                       targetcoords[,2] < ext_sub[4] & targetcoords[,2] > ext_sub[3],]

  trees <- st_coordinates(panes)

  iter <- nrow(tg)
  pb <- progress_bar$new(format = "Calculating [:bar] :current/:total (:percent) | elapsed: :elapsed",
                         total = iter)
  pb$tick(0)
  progress <- function() pb$tick()
  opts <- list(progress = progress)
  tree_data <-
    foreach(t = 1:nrow(tg),
            .combine = rbind,
            .packages = c("FNN",
                          "sf"),
            .options.snow = opts) %dopar% {

              target <- tg[t,]
              target_group <- fg[fg[,1] >= target[1] - 30.0 & fg[,1] <= target[1] + 30.0 &
                                   fg[,2] >= target[2] - 30.0 & fg[,2] <= target[2] + 30.0,]
              bbox <- st_bbox(b)
              bbox[1] <- target[1] - 30.0
              bbox[2] <- target[2] - 30.0
              bbox[3] <- target[1] + 30.0
              bbox[4] <- target[2] + 30.0
              bbox <- st_as_sfc(bbox)

              clip <- st_intersection(panes, bbox)
              clip$area <- st_area(clip)

              if (nrow(target_group) < 9 | nrow(clip) < 1) {

                mean_area <- NA
                median_area <- NA
                sd_area <- NA
                mean_frac <- NA

              } else {

                mean_area <- mean(clip$area)
                median_area <- median(clip$area)
                sd_area <- sd(clip$area)

                clip$perim <- st_length(st_cast(clip, "MULTILINESTRING"))
                mean_frac <- mean(mapply(calc_frac, clip$perim, clip$area))

              }

              data.frame(x = tg[t,1], y = tg[t,2],
                         mean_area = mean_area,
                         median_area = median_area,
                         sd_area = sd_area,
                         mean_frac = mean_frac)

            }

  write.csv(tree_data, paste0("processed/gap_data_", s, ".csv"))

}


check_gap <- read.csv("data/processed/gap_data_4.csv")

ras <- rasterFromXYZ(data.frame(x = check_gap$x, y = check_gap$y, z = check_gap$median_area))
plot(ras)

head(check_gap)


mean_gap_area_data <- data.frame(x = numeric(0), y = numeric(0), z = numeric(0))
for (i in 1:21) {

  nt <- read.csv(paste0("data/processed/gap_data_", i, ".csv"))
  mean_gap_area_data <- rbind(mean_gap_area_data, data.frame(x = nt$x, y = nt$y, z = nt$mean_area))

}

mean_gap_area_raster <- rasterFromXYZ(mean_gap_area_data)
plot(mean_gap_area_raster)


dens_data <- data.frame(x = numeric(0), y = numeric(0), z = numeric(0))
for (i in 1:31) {

  nt <- read.csv(paste0("data/processed/tree_data_", i, ".csv"))
  dens_data <- rbind(dens_data, data.frame(x = nt$x, y = nt$y, z = nt$mean_dens))

}

dens_raster <- rasterFromXYZ(dens_data)
plot(dens_raster)
