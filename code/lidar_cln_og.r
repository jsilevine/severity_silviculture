
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

setwd("Documents/Jacob/lidar_processing")

##---------------------------------------------------------------
## 1. point data
##---------------------------------------------------------------

point_test <- st_read("data/lidar/points/TAO__-42000_210000__t4_0p75_lp3__HighPoints.shp")

filelist <- list.files("lidar/points", ".shp")
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

template <- readRDS("templates/raster_template.rds")
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

##---------------------------------------------------------------
## 2. basin data
##---------------------------------------------------------------

basin_test <- st_read("data/lidar/basins/TAO__-42000_210000__t4_0p75_lp3__Basin_Map.shp")
bbox <- st_bbox(basin_test)
bbox[1] <- -42000
bbox[2] <- 210000
bbox[3] <- -40000
bbox[4] <- 212000
bbox <- st_as_sfc(bbox)

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

## for (i in 1:nrow(gridpts)) {

##   b <- st_read(create_filename(gridpts[i,"x"], gridpts[i,"y"]))
##   bbox <- st_bbox(b)
##   bbox[1] <- gridpts[i,"x"]
##   bbox[2] <- gridpts[i,"y"]
##   bbox[3] <- gridpts[i,"x"] + 2000
##   bbox[4] <- gridpts[i,"y"] + 2000
##   bbox <- st_as_sfc(bbox)

##   dif <- st_difference(bbox, st_union(b))
##   dif <- st_cast(dif, "POLYGON")
##   dif <- st_as_sf(dif)
##   st_write(dif, create_filename2(gridpts[i,"x"], gridpts[i,"y"]))

## }

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

## subs <- split(unique(gridpts$x), ceiling(seq_along(unique(gridpts$x))/2))

## for (i in 2:length(subs)) {
##   subs[[i]] <- c(subs[[i-1]][length(subs[[i-1]])], subs[[i]])
## }

st_fast_intersection <- function(x,y,...){
#faster replacement for st_intersection(x, y,...)

  y_subset <-
    st_intersects(x, y) %>%
    unlist() %>%
    unique() %>%
    sort() %>%
    {y[.,]}

  st_intersection(x, y_subset,...)
}


cl <- makeSOCKcluster(8)
registerDoSNOW(cl)

sample_radii <- c(30, 90, 180)

for (r in sample_radii) {

  for (s in 580:582) {

    print(paste0("Starting grid subset ", s, " of ", length(subs)))

    x <- as.matrix(dist(gridpts, upper = TRUE, diag = TRUE))
    w <- as.numeric(which(x[s,] < 2829))
    gridpts_sub <- gridpts[w,]
    panes <- st_read(create_filename(gridpts_sub[1,1], gridpts_sub[1,2]))
    if (nrow(gridpts_sub) > 1) {
      for (i in 2:nrow(gridpts_sub)) {
        newpane <- st_read(create_filename(gridpts_sub[i,1], gridpts_sub[i,2]))
        panes <- rbind(panes, newpane)
      }
    }
    panes <- st_union(panes)

    ext <- extent(gridpts[s,1], gridpts[s,1]+2000,
                  gridpts[s,2], gridpts[s,2]+2000)

    if (any(gridpts_sub$x < gridpts[s,"x"])) {
      ext[1] <- ext[1] - (r*2)
    }
    if (any(gridpts_sub$x > gridpts[s,"x"])) {
      ext[2] <- ext[2] + (r*2)
    }
    if (any(gridpts_sub$y < gridpts[s,"y"])) {
      ext[3] <- ext[3] - (r*2)
    }
    if (any(gridpts_sub$y > gridpts[s,"y"])) {
      ext[4] <- ext[4] + (r*2)
    }

    ## 2. get all target coords within pane
    fg <- fullcoords[fullcoords[,1] < ext[2] & fullcoords[,1] > ext[1] &
                       fullcoords[,2] < ext[4] & fullcoords[,2] > ext[3],]

    tg <- targetcoords[targetcoords[,1] < ext[2] & targetcoords[,1] > ext[1] &
                         targetcoords[,2] < ext[4] & targetcoords[,2] > ext[3],]

    if (nrow(matrix(tg, ncol = 2)) > 0) {
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
                  target_group <- fg[fg[,1] >= target[1] - r & fg[,1] <= target[1] + r &
                                       fg[,2] >= target[2] - r & fg[,2] <= target[2] + r,]
                  bbox <- st_bbox(panes)
                  bbox[1] <- target[1] - (r + 15)
                  bbox[2] <- target[2] - (r + 15)
                  bbox[3] <- target[1] + (r + 15)
                  bbox[4] <- target[2] + (r + 15)
                  bbox <- st_as_sfc(bbox)

                  clip <- st_intersection(panes, bbox)
                  gaps <- st_difference(bbox, st_union(clip))
                  gaps <- st_cast(gaps, "POLYGON")
                  gaps <- st_as_sf(gaps)
                  gaps$area <- st_area(gaps)

                  if (nrow(target_group) < (r/30)^2 | nrow(gaps) < 1) {

                    mean_area <- NA
                    median_area <- NA
                    sd_area <- NA
                    mean_frac <- NA
                    percent_open <- NA

                  } else {

                    mean_area <- as.numeric(mean(gaps$area))
                    median_area <- as.numeric(median(gaps$area))
                    sd_area <- as.numeric(sd(gaps$area))
                    percent_open <- as.numeric(sum(gaps$area) / ((r*2 + 30)^2))

                    gaps$perim <- st_length(st_cast(gaps, "MULTILINESTRING"))
                    mean_frac <- mean(mapply(calc_frac, gaps$perim, gaps$area))

                  }

                  data.frame(x = tg[t,1], y = tg[t,2],
                             mean_area = mean_area,
                             median_area = median_area,
                             sd_area = sd_area,
                             percent_open = percent_open,
                             mean_frac = mean_frac)

                }

      write.csv(tree_data, paste0("data/processed/gap_data_", r, "m_", s, ".csv"))

    }

  }

}

df <- read.csv("data/processed/gap_data_30m_1.csv", row.names = 1)
for (i in c(2,3,15:18)) {
  df <- rbind(df, read.csv(paste0("data/processed/gap_data_30m_", i, ".csv"), row.names = 1))
}

plot(rasterFromXYZ(df[,1:3], res = c(30,30), crs = st_crs(panes)))
plot(rasterFromXYZ(read.csv("data/processed/gap_data_30m_500.csv", row.names = 1)[,1:3], res = c(30,30), crs = st_crs(panes)))

plot(rasterFromXYZ(tree_data[,1:3], res = c(30,30), crs = st_crs(panes)))


ggplot(panes) +
  geom_sf(aes(fill = height))

+
  geom_tile(data = tree_data, aes(x = x, y = y, fill = mean_area), alpha = 0.3)
