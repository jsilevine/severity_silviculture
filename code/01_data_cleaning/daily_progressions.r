##---------------------------------------------------------------
## daily_progressions.r
##
##
##---------------------------------------------------------------

#library(rgdal)
library(rvest)
library(sf)
library(rmapshaper)
library(fasterize)
library(terra)
library(ggplot2)
library(chron)
library(concaveman)
library(sp)
library(geosphere)
library(parallel)
library(FNN)
library(doParallel)
library(doSNOW)
library(progress)

source("code/utility/utility_functions.r")

convert_to_time <- function(time) {
  out <- paste0(substring(time, 1,2), ":", substring(time, 3, 4), ":00")
  return(out)
}

convert_to_datetime <- function(date, time) {
  out <- chron(dates = as.character(date), times = time, format = c('y-m-d', 'h:m:s'))
  return(out)
}

small_clust <- function(clusters) {
  for (i in unique(clusters$cluster_id)) {
    if (nrow(clusters[clusters$cluster_id == i,]) < 4) {
      if (!exists("rem")) {
        rem <- i
      } else {
        rem <- append(rem, i)
      }
    }
  }
  if (!exists("rem")) rem <- 1e+50
  return(rem)
}

cluster <- function(data, h) {
  z <- st_coordinates(data)
  dist <- dist(z)
  cluster <- hclust(as.dist(dist), method = "complete")
  data$cluster_id <- cutree(cluster, h = h)
  data <- data[!(data$cluster_id %in% small_clust(data)),]
  return(data)
}

split_perims <- function(perims, groups = 7) {
  breaks <- ceiling(seq(0, nrow(perims), length.out = groups+1))
  split <- list()
  for (i in 1:groups) {
    split[[i]] <- perims[(breaks[i]+1):(breaks[i+1]),]
  }
  return(split)
}

extract_daily_perims <- function(start_date, end_date, perim,
                                 template, template_large, viirs,
                                 start_buffer = 0, end_buffer = 0,
                                 cores = 7, cluster_dist = 1125, large_cluster_dist = 8000,
                                 cluster_tol = 1, concavity = 3,
                                 is_dixie = FALSE, smooth_dist = 25) {

  print("Running VIIRS extraction")
  sf_use_s2(FALSE)
  viirs <- st_transform(viirs, crs = st_crs(perim))

  sub_viirs <- viirs[viirs$ACQ_DATE >= start_date-start_buffer &
                     viirs$ACQ_DATE <= end_date+end_buffer,]
  sub_viirs <- st_intersection(sub_viirs, perim)

  sub_viirs$ACQ_TIME <- sapply(sub_viirs$ACQ_TIME, convert_to_time)
  sub_viirs$datetime <- chron(mapply(convert_to_datetime, sub_viirs$ACQ_DATE, sub_viirs$ACQ_TIME))
  sub_viirs <- sub_viirs[order(sub_viirs$datetime),]

  print("Beginning perimeter delineation")

  ## estimate daily perims
  suppressMessages(for (t in unique(sub_viirs$datetime)) {
                     subdata <- sub_viirs[sub_viirs$datetime == t,]
                     if (nrow(subdata) < 4) next
                     subdata <- cluster(subdata, h = cluster_dist)
                     if (nrow(subdata) < 4) next
                     subdata$large_cluster <- cluster(subdata, h = large_cluster_dist)$cluster_id
                     if (length(unique(subdata$large_cluster)) < cluster_tol) {
                       for (c in unique(subdata$large_cluster)) {
                         hn <- concaveman(st_as_sf(st_union(subdata[subdata$large_cluster == c,])),
                                          concavity = concavity)
                         if (!exists("h")) {
                           h <- hn
                         } else h <- rbind(h, hn)
                       }
                       hull <- h
                       rm(h)
                     } else hull <- concaveman(st_as_sf(st_union(subdata)), concavity = concavity)
                     hull <- st_intersection(hull, perim)
                     hull <- st_as_sf(st_union(hull))
                     hull[,"datetime"] <- t
                     if (exists("hull_prev")) {
                       if(t - t_prev > 30) next
                       hull <- st_union(hull, hull_prev)
                       hull <- hull[, c("datetime")]
                       daily_perims <- rbind(daily_perims, hull)
                       hull_prev <- hull
                       t_prev <- t
                     } else {
                       hull <- hull[, c("datetime")]
                       daily_perims <- hull
                       hull_prev <- hull
                       t_prev <- t
                     }
                   })

  ## finish perimeter
  final_perim <- perim
  final_perim[, "datetime"] <- max(daily_perims$datetime)+1
  final_perim <- final_perim[, c("datetime", "Shape")]
  colnames(final_perim)[2] <- "polygons"
  st_geometry(final_perim) <- "polygons"
  colnames(daily_perims)[2] <- "polygons"
  st_geometry(daily_perims) <- "polygons"
  daily_perims <- rbind(daily_perims, final_perim)

  print("Finished perimeter delineation")
  print("Beginning perimeter differencing")
  ## get daily burned areas
  daily_burned <- daily_perims[1,]
  suppressMessages(for (i in 2:nrow(daily_perims)) {
                     n <- st_difference(daily_perims[i,], daily_perims[i-1,])
                     n <- n[, c("datetime")]
                     if(nrow(n) > 0) daily_burned[i,] <- n
                   })
  daily_burned <- daily_burned[!is.na(daily_burned$datetime),]
  print("Finished perimeter differencing")
  print("Rasterizing")
  perims_for_par <- split_perims(daily_burned, groups = cores)

  for (i in 1:length(perims_for_par)) {

    if (i == 1) {
      br <- terra::rasterize(st_collection_extract(perims_for_par[[i]], "POLYGON"),
                             template, field = "datetime")
    } else {
      br <- merge(br, terra::rasterize(st_collection_extract(perims_for_par[[i]], "POLYGON"),
                                       template, field = "datetime"))
    }
  }

  smallbr <- br

  smooth <- focal(smallbr, w = matrix(1,smooth_dist,smooth_dist), mean, na.rm = TRUE)

  smooth_round <- smooth
  values(smooth_round) <- round_to_nearest(values(smooth_round), 8)

  return(list(raw = br,
              smooth = smooth,
              smooth_round = smooth_round))

  print("Finished")
}

crs <- st_crs(3310)

viirs <- st_read("data/viirs/suomi_viirs/fire_archive_SV-C2_311481.shp")
viirs <- st_transform(viirs, crs)

raster_template <- rast(readRDS("data/templates/isforest.rds"))
raster_template_large <- aggregate(raster_template, fact = 12)

test <- resample(raster_template_large, raster_template, method = "near", na.rm = TRUE)


plot(raster_template)
plot(test)
plot(raster_template_large)

test <- snap_to_template(raster_template_large, raster_template, "near")
plot(test)

test.df <- as.data.frame(test, xy = TRUE, na.rm = TRUE)
raster_template.df <- as.data.frame(raster_template, xy = TRUE, na.rm = TRUE)

head(test.df)


nrow(raster_template.df)
nrow(test.df)

plot(snap_to_template(raster_template_large, raster_template, "near"))

t <- rast("data/fire_progression/smooth_round/dixie_progression.tif")


fires <- c("dixie", "northcomplex", "sheep", "sugar", "walker")

cl <- makeSOCKcluster(10)
registerDoSNOW(cl)

system.time({

  for (f in fires) {

    print(paste("Processing the", f, "fire"))

    perim <- readRDS(paste0("data/perimeters/", f, "_perim.rds"))
    perim <- st_transform(perim, crs)

    if (f == "dixie") {
      start <- as.Date("2021/07/13")
      end <- as.Date("2021/10/24")
    } else {
      start <- as.Date(perim$ALARM_DATE)
      end <- as.Date(perim$CONT_DATE)
    }

    progression <- extract_daily_perims(start_date = start, end_date = end,
                                        perim = perim, template = raster_template,
                                        template_large = raster_template_large,
                                        viirs = viirs, is_dixie = (f == "dixie"),
                                        cluster_tol = 4, cluster_dist = 1125,
                                        large_cluster_dist = 10000)

    for (i in 1:length(progression)) {

      saveRDS(progression[[i]], paste0("data/fire_progression/", names(progression)[i], "/", f, "_progression.rds"))
      writeRaster(progression[[i]], paste0("data/fire_progression/", names(progression)[i], "/", f, "_progression.tif"), overwrite = TRUE)
      df <- as.data.frame(progression[[i]], xy = TRUE, na.rm = TRUE)
      colnames(df)[3] <- "datetime"
      write.csv(df, paste0("data/fire_progression/", names(progression)[i], "/", f, "_progression.csv"))

    }

  }

})

perim <- readRDS("data/perimeters/walker_perim.rds")
perim <- st_transform(perim, 3310)
test <- rast("data/fire_progression/smooth_round/walker_progression.tif")
test <- crop(test, perim)
plot(test)
