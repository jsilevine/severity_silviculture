
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
library(terra)


##---------------------------------------------------------------
## 1. basin data
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

template <- readRDS("data/templates/raster_template.rds")
template <- projectRaster(template, res = c(30,30), crs = crs(basin_test))
targetcoords <- as.data.frame(template, xy = TRUE, na.rm = TRUE)
targetcoords <- as.matrix(targetcoords[,1:2])
fullcoords <- as.data.frame(template, xy = TRUE, na.rm = FALSE)
fullcoords <- as.matrix(fullcoords[,1:2])

bb_template <- st_bbox(basin_test)

gridpts <- gridpts[gridpts$x > min(targetcoords[,1]) & gridpts$x < max(targetcoords[,1]) &
                     gridpts$y > min(targetcoords[,2]) & gridpts$y < max(targetcoords[,2]),]

cl <- makeSOCKcluster(10)
registerDoSNOW(cl)

for (s in 1:nrow(gridpts)) {

  print(paste0("Starting grid subset ", s, " of ", nrow(gridpts)))

  print("Loading crown polygons")
  x <- as.matrix(dist(gridpts, upper = TRUE, diag = TRUE))
  w <- as.numeric(which(x[s,] < 2829))
  gridpts_sub <- gridpts[w,]

  panes <- invisible(st_read(create_filename(gridpts_sub[1,1], gridpts_sub[1,2])))
  if (nrow(gridpts_sub) > 1) {
    for (i in 2:nrow(gridpts_sub)) {
      newpane <- invisible(st_read(create_filename(gridpts_sub[i,1], gridpts_sub[i,2])))
      panes <- rbind(panes, newpane)
    }
  }

  ext <- extent(gridpts[s,1], gridpts[s,1] + 2000,
                gridpts[s,2], gridpts[s,2] + 2000)

  ext2 <- ext
  ext2[1] <- ext2[1] - 45
  ext2[2] <- ext2[2] + 45
  ext2[3] <- ext2[3] - 45
  ext2[4] <- ext2[4] + 45

  print("Cropping crown polygons")
  panes <- st_crop(panes, st_bbox(ext2))

  ## 2. get all target coords within pane
  tg <- targetcoords[targetcoords[,1] < ext[2] & targetcoords[,1] > ext[1] &
                       targetcoords[,2] < ext[4] & targetcoords[,2] > ext[3],]

  if (nrow(matrix(tg, ncol = 2)) == 1) {
    tgdf <- data.frame(x = tg["x"], y = tg["y"])
  } else {
    tgdf <- as.data.frame(tg)
  }
  cells <- st_as_sf(tgdf, coords = c("x", "y"), crs = st_crs(panes))

  ## ggplot() +
  ##   geom_sf(data = panes, aes(fill = height)) +
  ##   geom_sf(data = cells, color = "black") +
  ##   theme_bw()

  ## plot(panes)

  print("Calculating cover")
  if (nrow(matrix(tg, ncol = 2)) > 0) {
    panes <- st_union(panes)
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
    cover_data <-
      foreach(t = 1:iter,
              .combine = rbind,
              .packages = c("FNN",
                            "sf",
                            "raster")) %dopar% {

                cell_extent <- extent(tg[t,1]-15, tg[t,1]+15, tg[t,2]-15, tg[t,2]+15)

                cell <- st_as_sf(as(cell_extent, "SpatialPolygons"))
                st_crs(cell) <- st_crs(panes)
                check <- st_intersection(panes, cell)
                cover <- as.numeric(sum(st_area(check))) / 30^2

                data.frame(x = tg[t,1], y = tg[t,2],
                           cover = cover)

              }

    write.csv(cover_data, paste0("data/processed/cover_data_", s, ".csv"))

  }

}

filelist <- list.files("data/processed", "cover_data",)
data <- read.csv(paste0("data/processed/", filelist[1]), row.names = 1)
for (i in 2:length(filelist)) {
  data <- rbind(data, read.csv(paste0("data/processed/", filelist[i]), row.names = 1))
}
head(data)

plot(rasterFromXYZ(data[,c(1,2,3)], res = c(30, 30)))
plot(gridpts$x, gridpts$y, add = TRUE)
write.csv(data, "data/processed/cover_full.csv")

data <- read.csv("data/processed/cover_full.csv")

## create binary forest/not forest map
data$isforest <- 1
data[data$cover <= 0.1, "isforest"] <- 0

library(sp)

isforest_raster <- rasterFromXYZ(data[,c(1,2,4)], res = c(30,30), crs = crs(template))
saveRDS(isforest_raster, "data/templates/isforest.rds")
writeRaster(isforest_raster, "data/templates/isforest.tif", overwrite = TRUE)
isforest_df <- as.data.frame(isforest_raster, xy = TRUE, na.rm = TRUE)
write.csv(isforest_df, "data/templates/isforest.csv")

isforest_template <- isforest_raster
isforest_template$isforest[isforest_template$isforest == 0] <- NA
saveRDS(isforest_template, "data/templates/isforest_template.rds")
writeRaster(isforest_template, "data/templates/isforest_template.tif", overwrite = TRUE)

cover_raster <- rasterFromXYZ(data[,c(1,2,3)], res = c(30,30), crs = crs(template))
saveRDS(cover_raster, "data/processed/cover.rds")
writeRaster(cover_raster, "data/processed/cover.tif", overwrite = TRUE)
