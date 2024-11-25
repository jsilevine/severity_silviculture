
##---------------------------------------------------------------
## canopy_cover_cln.r
## SCRIPT TO CLEAN CANOPY COVER DATA FROM LIDAR
## BY: JACOB LEVINE - jacoblevine@princeton.edu
## 11/26/22
##---------------------------------------------------------------

##---------------------------------------------------------------
## 0. libraries
##---------------------------------------------------------------

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

point_test <- st_read("data/lidar/points/TAO__-112000_176000__t4_0p75_lp3__HighPoints.shp")

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

#template <- readRDS("data/templates/raster_template.rds")
#template <- projectRaster(template, res = c(30,30), crs = crs(point_test))

template <- readRDS("data/templates/isforest_template.rds")

targetcoords <- as.data.frame(template, xy = TRUE, na.rm = TRUE)
targetcoords <- as.matrix(targetcoords[,1:2])
fullcoords <- as.data.frame(template, xy = TRUE, na.rm = FALSE)
fullcoords <- as.matrix(fullcoords[,1:2])

create_filename <- function(x, y) {
 return(paste0("data/lidar/points/TAO__", format(x, scientific = FALSE), "_",
               format(y, scientific = FALSE), "__t4_0p75_lp3__HighPoints.shp"))
}

gridpts <- gridpts[gridpts$x > min(targetcoords[,1]) & gridpts$x < max(targetcoords[,1]) &
                   gridpts$y > min(targetcoords[,2]) & gridpts$y < max(targetcoords[,2]),]

cl <- makeSOCKcluster(10)
registerDoSNOW(cl)


## need to rerun s = 393
sample_radii <- c(30, 180)

for (r in sample_radii) {

  for (s in 1:nrow(gridpts)) {

    print(paste0("Starting grid subset ", s, " of ", nrow(gridpts)))

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

    ext <- extent(gridpts[s,1], gridpts[s,1]+2000,
                  gridpts[s,2], gridpts[s,2]+2000)

    tg <- targetcoords[targetcoords[,1] < ext[2] & targetcoords[,1] > ext[1] &
                         targetcoords[,2] < ext[4] & targetcoords[,2] > ext[3],]

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

    trees <- st_coordinates(panes)

    if (nrow(matrix(tg, ncol = 2)) > 0) {
      if (nrow(matrix(tg, ncol = 2)) == 1) {
        iter <- 1
        tg <- matrix(tg, ncol = 2)
      } else {
        iter <- nrow(tg)
      }

      #pb <- progress_bar$new(format = "Calculating [:bar] :current/:total (:percent) | elapsed: :elapsed",
      #                       total = iter)
      #pb$tick(0)
      #progress <- function() pb$tick()
      #opts <- list(progress = progress)

      tree_data <-
        foreach(t = 1:iter,
              .combine = rbind,
              .packages = c("FNN")) %dopar% {

                target <- tg[t,]
                target_group <- fg[fg[,1] >= target[1] - r & fg[,1] <= target[1] + r &
                                   fg[,2] >= target[2] - r & fg[,2] <= target[2] + r,]
                tree_clip <- trees[trees[,1] > target[1] - (r+15) & trees[,1] < target[1] + (r+15) &
                                   trees[,2] > target[2] - (r+15) & trees[,2] < target[2] + (r+15), ]
                ht_clip <- panes[trees[,1] > target[1] - (r+15) & trees[,1] < target[1] + (r+15) &
                                 trees[,2] > target[2] - (r+15) & trees[,2] < target[2] + (r+15), ]

                tree_clip <- matrix(tree_clip, ncol = 2)
                if (nrow(target_group) < (r/30)^2 | nrow(tree_clip) < 1) {
                  clust <- NA
                  mean_dens <- NA
                  mean_ht <- NA
                  median_ht <- NA
                  max_ht <- NA
                } else {

                  if (nrow(tree_clip) == 1) {
                    clust <- NA
                  } else {
                    re <- 0.5 * sqrt((30+(2*r))*(30+(2*r)) / nrow(tree_clip))
                    r0 <- mean(knn.dist(tree_clip, k = 1))
                  }

                  clust <- r0 / re
                  mean_dens <- (nrow(tree_clip) / ((30+(2*r))*(30+(2*r)))) * 10000

                  mean_ht <- mean(ht_clip$height)
                  median_ht <- median(ht_clip$height)
                  max_ht <- max(ht_clip$height)

                  }

                  data.frame(x = tg[t,1], y = tg[t,2],
                             clust = clust, mean_dens = mean_dens,
                             mean_ht = mean_ht, median_ht = median_ht, max_ht = max_ht)
              }

      write.csv(tree_data, paste0("data/processed/tree_data_", r, "m_", s, ".csv"))

    }

  }

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


x <- as.matrix(dist(gridpts, upper = TRUE, diag = TRUE))
boundary <- list()
for (i in 1:nrow(x)) {
    ext <- extent(gridpts[i,1], gridpts[i,1]+2000,
                  gridpts[i,2], gridpts[i,2]+2000)
    nbox <- st_as_sf(as(ext, "SpatialPolygons"))
    st_crs(nbox) <- st_crs(3310)
    boundary[[i]] <- nbox
}
boundary <- st_union(dplyr::bind_rows(boundary))

plot(boundary)

bb <- st_bbox(boundary)
bb[1] <- bb[1] - 10000
bb[2] <- bb[2] - 10000
bb[3] <- bb[3] + 10000
bb[4] <- bb[4] + 10000
bb <- st_as_sfc(bb)

outline <- st_difference(bb, boundary)


create_filename <- function(x, y) {
  return(paste0("data/lidar/basins/TAO__", format(x, scientific = FALSE), "_",
                format(y, scientific = FALSE), "__t4_0p75_lp3__Basin_Map.shp"))
}
create_filename2 <- function(x, y) {
  return(paste0("data/lidar/gaps/TAO__", format(x, scientific = FALSE), "_",
                format(y, scientific = FALSE), "__t4_0p75_lp3__GAP.shp"))
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

cl <- makeSOCKcluster(10)
registerDoSNOW(cl)

sample_radii <- c(30)

for (r in sample_radii) {

  #for (s in 1:nrow(gridpts)) {
  for (s in 1:10) {

    print(paste0("Starting grid subset ", s, " of ", nrow(gridpts)))

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

    ext <- extent(gridpts[s,1], gridpts[s,1]+2000,
                  gridpts[s,2], gridpts[s,2]+2000)

    tg <- targetcoords[targetcoords[,1] < ext[2] & targetcoords[,1] > ext[1] &
                         targetcoords[,2] < ext[4] & targetcoords[,2] > ext[3],]

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

    fg <- fullcoords[fullcoords[,1] < ext[2] & fullcoords[,1] > ext[1] &
                        fullcoords[,2] < ext[4] & fullcoords[,2] > ext[3],]


    if (nrow(matrix(tg, ncol = 2)) == 1) {
      tgdf <- data.frame(x = tg["x"], y = tg["y"])
    } else {
      tgdf <- as.data.frame(tg)
    }
    cells <- st_as_sf(tgdf, coords = c("x", "y"), crs = st_crs(panes))

    d <- as.numeric(st_distance(cells, outline))
    tgdf <- tgdf[d > r+15,] ## remove all pixels too close to edge to properly calculate metrics


    if (nrow(tgdf) > 0) {
      panes <- st_union(panes)
      if (nrow(tgdf) == 1) {
        iter <- 1
        tg <- tgdf
      } else {
        tg <- tgdf
        iter <- nrow(tgdf)
      }
      #pb <- progress_bar$new(format = "Calculating [:bar] :current/:total (:percent) | elapsed: :elapsed",
      #                       total = iter)
      #pb$tick(0)
      #progress <- function() pb$tick()
      #opts <- list(progress = progress)
      tree_data <-
        foreach(t = 1:iter,
                .combine = rbind,
                .packages = c("FNN",
                              "sf")) %dopar% {

                  target <- tg[t,]
                  #target_group <- fg[fg[,1] >= target[1] - r & fg[,1] <= target[1] + r &
                  #                     fg[,2] >= target[2] - r & fg[,2] <= target[2] + r,]
                  bbox <- st_bbox(panes)
                  bbox[1] <- as.numeric(target[1]) - (r + 15)
                  bbox[2] <- as.numeric(target[2]) - (r + 15)
                  bbox[3] <- as.numeric(target[1]) + (r + 15)
                  bbox[4] <- as.numeric(target[2]) + (r + 15)
                  bbox <- st_as_sfc(bbox)

                  clip <- st_intersection(panes, bbox)
                  gaps <- st_difference(bbox, st_union(clip))
                  gaps <- st_cast(gaps, "POLYGON")
                  gaps <- st_as_sf(gaps)
                  gaps$area <- st_area(gaps)

                  if (nrow(gaps) < 1) {

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


print(head(data), digits = 15)

filelist <- list.files("data/processed", "tree_data_30m",)
data <- read.csv(paste0("data/processed/", filelist[1]), row.names = 1)
for (i in 2:length(filelist)) {
  data <- rbind(data, read.csv(paste0("data/processed/", filelist[i]), row.names = 1))
}

head(data)
plot(rasterFromXYZ(data[,c(1,2,3)], res = c(30, 30), crs = crs('+init=EPSG:3310')))
write.csv(data, "data/processed/tree_data_30m_full.csv")

filelist <- list.files("data/processed", "tree_data_180m",)
data <- read.csv(paste0("data/processed/", filelist[1]), row.names = 1)
for (i in 2:length(filelist)) {

  print(i)
  data <- rbind(data, read.csv(paste0("data/processed/", filelist[i]), row.names = 1 ))
}
plot(rasterFromXYZ(data[,c(1,2,5)], res = c(30, 30)))
write.csv(data, "data/processed/tree_data_180m_full.csv")


filelist <- list.files("data/processed", "gap_data_30m",)
data <- read.csv(paste0("data/processed/", filelist[1]), row.names = 1)
for (i in 2:10) {

  print(i)
  data <- rbind(data, read.csv(paste0("data/processed/", filelist[i]), row.names = 1 ))
}
head(data)
plot(rasterFromXYZ(data[,c(1,2,6)], res = c(30, 30)))
write.csv(data, "data/processed/gap_data_30m_full.csv")

filelist <- list.files("data/processed", "gap_data_180m",)
data <- read.csv(paste0("data/processed/", filelist[1]), row.names = 1)
for (i in 2:length(filelist)) {
  print(i)
  data <- rbind(data, read.csv(paste0("data/processed/", filelist[i]), row.names = 1 ))
}
head(data)
plot(rast(rasterFromXYZ(data[,c(1,2,6)], res = c(30, 30))))
write.csv(data, "data/processed/gap_data_180m_full.csv")

sum(!is.na(data$mean_area))
