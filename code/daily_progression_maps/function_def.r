##---------------------------------------------------------------
## Function definition for daily fire progression mapping
## Author: Jacob Levine; jacoblevine@princeton.edu
##---------------------------------------------------------------

download_viirs_data <- function(fire_perimeter, fire_start, fire_end,
                                mapkey) {

  bb <- st_bbox(fire_perimeter)

  dates <- seq(fire_start, fire_end, by = 1)
  date_splt <- split(dates, ceiling(seq_along(dates)/10))

  if (!dir.exists("viirs_data")) {
    dir.create("viirs_data")
  }

  for (i in 1:length(date_splt)) {

    dl_url <- paste0("https://firms.modaps.eosdis.nasa.gov/api/area/csv/",
                     mapkey,
                     "/VIIRS_SNPP_SP/",
                     bb[1], ",", bb[2], ",", bb[3], ",", bb[4],
                     "/",
                     length(date_splt[[i]]),
                     "/", date_splt[[i]][1])

    download.file(dl_url,
                  paste0("viirs_data/viirs_", i, ".csv"))

  }

  full_viirs <- read.csv("viirs_data/viirs_1.csv")

  for (i in 2:length(date_splt)) {

    full_viirs  <- rbind(full_viirs, read.csv(paste0("viirs_data/viirs_", i, ".csv")))

  }

  viirs_shp <- st_as_sf(full_viirs, coords = c(2,1), crs = st_crs(4326))

  write.csv(full_viirs, "viirs_data/full_viirs.csv")
  st_write(viirs_shp, "viirs_data/full_viirs.shp")

  return(viirs_shp)

}

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
  z <- as_Spatial(data)
  dist <- distm(z)
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

extract_daily_perims <- function(start_date, end_date, perim, viirs,
                                 raster_template,
                                 start_buffer = 0, end_buffer = 0,
                                 cores = 7, cluster_dist = 1125, large_cluster_dist = 8000,
                                 cluster_tol = 1, concavity = 3) {

  print("Running VIIRS extraction")
  sf_use_s2(FALSE)
  viirs <- st_transform(viirs, crs = st_crs(perim))

  sub_viirs <- viirs[viirs$acq_date >= start_date-start_buffer &
                     viirs$acq_date <= end_date+end_buffer,]
  sub_viirs <- st_intersection(sub_viirs, perim)

  sub_viirs$acq_time <- sapply(sub_viirs$acq_time, convert_to_time)
  sub_viirs$datetime <- chron(mapply(convert_to_datetime, sub_viirs$acq_date, sub_viirs$acq_time))
  sub_viirs <- sub_viirs[!is.na(sub_viirs$acq_date),]
  sub_viirs <- sub_viirs[order(sub_viirs$datetime),]

  print("Beginning perimeter delineation")


  ## estimate daily perims using convex hull algorithm.
  ## algorithm adapted from Briones-Herrera 2020, Remote Sensing
  ## https://doi.org/10.3390/rs12122061
  suppressMessages(for (t in unique(sub_viirs$datetime)) {
                   subdata <- sub_viirs[sub_viirs$datetime == t,]
                   subdata <- subdata[!is.na(subdata$acq_date),]
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
  ## assumes that the user supplied perimeter is the true final perimeter,
  ## assigns to last timestep in specified range +1 day
  final_perim <- perim
  final_perim[, "datetime"] <- max(daily_perims$datetime)+1
  final_perim <- final_perim[, c("datetime", "Shape")]
  colnames(final_perim)[2] <- "polygons"
  st_geometry(final_perim) <- "polygons"
  colnames(daily_perims)[2] <- "polygons"
  st_geometry(daily_perims) <- "polygons"
  daily_perims <- rbind(daily_perims, final_perim)
  print("Finished perimeter delineation")

  ## get daily burned areas by sequntially differencing daily perimeters
  print("Beginning perimeter differencing")
  daily_burned <- daily_perims[1,]
  suppressMessages(for (i in 2:nrow(daily_perims)) {
                     n <- st_difference(daily_perims[i,], daily_perims[i-1,])
                     n <- n[, c("datetime")]
                     if(nrow(n) > 0) daily_burned[i,] <- n
                   })
  daily_burned <- daily_burned[!is.na(daily_burned$datetime),]
  print("Finished perimeter differencing")

  ## if user supplies a template raster, rasterize the daily perimeters and return
  if (exists("raster_template")) {
    print("Rasterizing")
    perims_for_par <- split_perims(daily_burned, groups = cores)
    br_list <- mclapply(perims_for_par,
                        function(x) fasterize(st_collection_extract(x, "POLYGON"),
                                              raster_template, field = "datetime"),
                        mc.cores = cores)
    br <- merge(br_list[[1]], br_list[[2]],
                br_list[[3]], br_list[[4]],
                br_list[[5]], br_list[[6]],
                br_list[[7]])
    return(br)

  ## otherwise simply return the daily polygons
  } else {

    return(daily_burned)

  }

  print("Finished")

}
