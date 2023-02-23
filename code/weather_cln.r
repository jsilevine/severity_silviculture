##---------------------------------------------------------------
## weather_cln.r
## SCRIPT TO DOWNLOAD WEATHER DATA, CALCULATE WEIGHTED AVERAGES
## BY: JACOB LEVINE - jacoblevine@princeton.edu
## 11/26/22
##---------------------------------------------------------------

## NOTE: for script to work, need to sign up for API key and token from:
## https://developers.synopticdata.com/mesonet/
## api_token argument required for multiple functions, see below


##---------------------------------------------------------------
## 0. Load prereq libraries and set working directory
##---------------------------------------------------------------
library(rgdal)
library(sf)
library(fasterize)
library(raster)
library(ggplot2)
library(chron)
library(ncdf4)
library(RAWSmet)
library(rjson)
library(basemaps)
library(firebehavioR)
library(data.table)
library(foreach)
library(doParallel)
library(doSNOW)
library(progress)
library(ggnewscale)
setwd("../")

##---------------------------------------------------------------
## 1. Function definition
##---------------------------------------------------------------

distance_elev_weighted_avg <- function(vals, d, elevations, focal_elevation, max_elev) {
  el_dist <- abs(focal_elevation - elevations)
  el_dist <- el_dist / max_elev
  weights <- (1/(d*el_dist))
  avg <- sum(weights*vals, na.rm = TRUE) / sum(weights[!is.na(vals)])
  return(avg)
}

correct_missing <- function(x, expected_length) {
  for (i in 1:expected_length) {
    if(is.null(x[[i]])) x[[i]] <- NA
  }
  return(x)
}

jsontochron <- function(datetime) {
  datetime <- gsub("T", ",", datetime)
  datetime <- gsub("Z", "", datetime)
  date <- strsplit(datetime, ",")[[1]][1]
  time <- strsplit(datetime, ",")[[1]][2]
  out <- chron(dates = date, time = time, format = c(dates = "y-m-d", times = "h:m:s"))
  return(out)
}

chartochron <-function(datetime) {
  datetime <- gsub("\\(", "", datetime)
  datetime <- gsub("\\)", "", datetime)
  date <- strsplit(datetime, " ")[[1]][1]
  time <- strsplit(datetime, " ")[[1]][2]
  out <- chron(dates = date, time = time, format = c(dates = "m/d/y", times = "h:m:s"))
  return(out)
}

pull_station_list <- function(perims, search_distance = 120,
                              api_base_url = "https://api.synopticdata.com/v2/",
                              api_token = "8a99301c9d054a60b7f4711a7a6ac270") {

  ## get centroid for searching
  centroid <- st_centroid(st_union(st_make_valid(perims)))
  x <- st_coordinates(centroid)[,"X"]
  y <- st_coordinates(centroid)[,"Y"]

  ## convert to miles
  mi_rad <- search_distance*0.621371

  ## first get list of all stations within 120km of fire centroid
  full_url <- paste0(api_base_url,
                   "stations/metadata?token=", api_token,"&radius=",
                   as.character(y), ",",
                   as.character(x), ",",
                   as.character(mi_rad), "&limit=50&network=2&output=json")

  download.file(full_url, destfile = "data/weather/station_list.json")

  ## read json to dataframe
  station_list <- rjson::fromJSON(file = "data/weather/station_list.json")
  st_ids <- character()
  lat <- numeric()
  lon <- numeric()
  elevation <- numeric()
  for (i in 1:length(station_list$STATION)) {
    st_ids <- append(st_ids, station_list$STATION[[i]]$STID)
    lat <- append(lat, station_list$STATION[[i]]$LATITUDE)
    lon <- append(lon, station_list$STATION[[i]]$LONGITUDE)
    if (length(as.numeric(station_list$STATION[[i]]$ELEVATION)) == 0) {
      elevation <- append(elevation, NA)
    } else elevation <- append(elevation, as.numeric(station_list$STATION[[i]]$ELEVATION)*0.3048)
  }

  stations <- data.frame(st_id = st_ids, lat = lat, lon = lon, elevation = elevation)
  stations_sf <- st_as_sf(stations, coords = c("lon", "lat"), crs = 4326, agr = "constant")
  saveRDS(stations_sf, "data/weather/stations/stations_sf.rds")
  write.csv(stations, "data/weather/stations/stations.csv")
}

pull_fire_weather <- function(fire_name, perims, buffer_dist = 20000,
                              api_base_url = "https://api.synopticdata.com/v2/",
                              api_token = "8a99301c9d054a60b7f4711a7a6ac270") {

  ## get list of stations within buffer dist of fire perimeter
  fire_perim <- perims[perims$FIRE_NA == toupper(fire_name), ]
  fire_name <- gsub(" ", "", fire_name)
  stations <- read.csv("data/weather/stations/stations.csv")
  stations_sf <- readRDS("data/weather/stations/stations_sf.rds")

  b <- TRUE
  while (b) {
    perim_buffer <- st_transform(st_buffer(st_make_valid(fire_perim), buffer_dist), st_crs(stations_sf))
    stations_sub <- st_intersects(stations_sf, perim_buffer)
    if (nrow(stations[as.data.frame(stations_sub)[,1],]) > 2) {
      b <- FALSE
    } else buffer_dist <- buffer_dist + 10000
  }
  stations <- stations[as.data.frame(stations_sub)[,1],]
  station_list <- stations[1, "st_id"]
  for (i in 1:nrow(stations)) {
    station_list <- paste0(station_list, ",", stations[i,"st_id"])
  }

  ## get fire start and end times
  daily_burned <- as.data.frame(fread(paste0("data/fire_progression/", fire_name, "_daily_burned.csv")))
  ## deal with chron not being able to read itself when character
  daily_burned$datetime <- as.chron(daily_burned$datetime)
  endtime <- max(daily_burned$datetime, na.rm = TRUE)
  starttime <- min(daily_burned$datetime, na.rm = TRUE)

  ## if time period too long API wont pull data, halve if longer than 60 days
  if (endtime - starttime > 60) {
    ts <- c(starttime-0.25, chron(as.numeric(endtime) - ((endtime - starttime) / 2)), endtime)
  } else ts <- c(starttime-0.25, endtime)

  ## convert to units interpretable by API
  timefull <- character(length(ts))
  for (i in 1:length(ts)) {
    date <- gsub("-", "", as.character(as.Date(ts[i])))
    hour <- hours(ts[i])
    if (nchar(hour) < 2) {
      hour <- paste0("0", as.character(hour))
    } else hour <- as.character(hour)
    minute <- minutes(ts[i])
    if (nchar(minute) < 2) {
      minute <- paste0("0", as.character(minute))
    } else minute <- as.character(minute)
    timefull[i] <- paste0(date, hour, minute)
  }

  ## download data from synoptic mesonet (RAWS) API
  weather_dir <- paste0("data/weather/", fire_name, "_weather")
  dir.create(weather_dir, showWarnings = FALSE)
  for (i in 1:(length(timefull)-1)) {

    destfile <- paste0(weather_dir, "/", fire_name, "_weather_", i, ".json")
    full_url <- paste0(api_base_url,
                       "stations/timeseries?token=", api_token,
                       "&start=", timefull[i],
                       "&end=", timefull[i+1],
                       "&vars=air_temp,relative_humidity,wind_speed,fuel_temp,fuel_moisture",
                       "&stid=", station_list,
                       "&output=json")

    download.file(full_url, destfile = destfile)
  }

  ## read jsons and convert to dataframe
  files <- list.files(weather_dir)
  jsons <- files[grepl(".json", files)]
  air_temp <- numeric()
  relative_humidity <- numeric()
  wind_speed <- numeric()
  fuel_temp <- numeric()
  fuel_moisture <- numeric()
  datetime <- numeric()
  st_id <- numeric()
  for (i in 1:length(jsons)) {
    data <- rjson::fromJSON(file = paste0(weather_dir, "/", jsons[i]))
    for (s in 1:length(data$STATION)) {
      if (!is.null(data$STATION[[s]]$OBSERVATIONS$fuel_temp_set_1)) {
        datetime <- c(datetime, unlist(data$STATION[[s]]$OBSERVATIONS$date_time))
        air_temp <- c(air_temp, unlist(correct_missing(data$STATION[[s]]$OBSERVATIONS$air_temp_set_1, length(data$STATION[[s]]$OBSERVATIONS$date_time))))
        relative_humidity <- c(relative_humidity, unlist(correct_missing(data$STATION[[s]]$OBSERVATIONS$relative_humidity_set_1, length(data$STATION[[s]]$OBSERVATIONS$date_time))))
        wind_speed <- c(wind_speed, unlist(correct_missing(data$STATION[[s]]$OBSERVATIONS$wind_speed_set_1,length(data$STATION[[s]]$OBSERVATIONS$date_time))))
        fuel_moisture <- c(fuel_moisture, unlist(correct_missing(data$STATION[[s]]$OBSERVATIONS$fuel_moisture_set_1,length(data$STATION[[s]]$OBSERVATIONS$date_time))))
        fuel_temp <- c(fuel_temp, unlist(correct_missing(data$STATION[[s]]$OBSERVATIONS$fuel_temp_set_1, length(data$STATION[[s]]$OBSERVATIONS$date_time))))
        st_id <- c(st_id, rep(data$STATION[[s]]$STID, times = length(data$STATION[[s]]$OBSERVATIONS$date_time)))
      }
    }
  }

  weather <- data.frame(st_id = st_id, datetime = datetime,
                        air_temp = air_temp, relative_humidity = relative_humidity,
                        wind_speed = wind_speed, fuel_moisture = fuel_moisture, fuel_temp)
  weather$datetime <- sapply(weather$datetime, FUN = jsontochron)
  weather$datetime <- chron(weather$datetime)
  return(weather)
}

calc_time_avgd_weather <- function(fire_name, perims, return = FALSE, buffer_dist = 20000) {

  print("Pulling weather data from API")
  weather_data <- pull_fire_weather(fire_name, perims = perims,
                                    buffer_dist = buffer_dist)
  b <- TRUE
  while (b) {
    stations <- read.csv("data/weather/stations/stations.csv", row.names = 1)
    ## subset station data to only include stations with data pulled by API
    ## (other stations did not have data avail for requested time periods)
    stations <- stations[stations$st_id %in% unique(weather_data$st_id),]
    if (nrow(stations) > 2) {
      b <- FALSE
    } else {
      buffer_dist <- buffer_dist + 5000
      weather_data <- pull_fire_weather(fire, perims = perims, buffer_dist = buffer_dist)
    }
  }

  fire_name <- gsub(" ", "", fire_name)
  daily_burned <- as.data.frame(fread(paste0("data/fire_progression/", fire_name, "_daily_burned.csv")))
  daily_burned$datetime <- as.chron(daily_burned$datetime)
  step_times <- unique(daily_burned$datetime)[!is.na(unique(daily_burned$datetime))]
  step_times <- chron(step_times[order(step_times)])

  write.csv(stations, paste0("data/weather/", fire_name, "_weather/stations.csv"))

  avg_air_temp <- matrix(NA, nrow = length(step_times), ncol = nrow(stations))
  max_air_temp <- matrix(NA, nrow = length(step_times), ncol = nrow(stations))
  avg_wind_speed <- matrix(NA, nrow = length(step_times), ncol = nrow(stations))
  max_wind_speed <- matrix(NA, nrow = length(step_times), ncol = nrow(stations))
  avg_relative_humidity <- matrix(NA, nrow = length(step_times), ncol = nrow(stations))
  avg_fuel_moisture <- matrix(NA, nrow = length(step_times), ncol = nrow(stations))
  avg_fuel_temp <- matrix(NA, nrow = length(step_times), ncol = nrow(stations))

  ## loop over time intervals and average data for each station
  for (t in 1:length(step_times)) {
    if (t == 1) {
      interval <- c(step_times[t] - 0.25, step_times[t])
    } else interval <- c(step_times[t-1], step_times[t])
    dat <- data.frame(st_id = stations$st_id)

    ##avgd vals
    ag <- aggregate(. ~ st_id,
                    data = weather_data[weather_data$datetime > interval[1] & weather_data$datetime <= interval[2],],
                    FUN = mean)
    vals <- merge(dat, ag, by = "st_id", all.x = TRUE)
    avg_air_temp[t,] <- vals$air_temp
    avg_wind_speed[t,] <- vals$wind_speed
    avg_relative_humidity[t,] <- vals$relative_humidity
    avg_fuel_moisture[t,] <- vals$fuel_moisture
    avg_fuel_temp[t,] <- vals$fuel_temp

    ##maxd vals
    ag <- aggregate(. ~ st_id,
                    data = weather_data[weather_data$datetime > interval[1] & weather_data$datetime <= interval[2],],
                    FUN = max)
    vals <- merge(dat, ag, by = "st_id", all.x = TRUE)
    max_air_temp[t,] <- vals$air_temp
    max_wind_speed[t,] <- vals$wind_speed
  }

  ## collect into list
  weather_list <- list(avg_air_temp = avg_air_temp,
                       max_air_temp = max_air_temp,
                       avg_wind_speed = avg_wind_speed,
                       max_wind_speed = max_wind_speed,
                       avg_relative_humidity = avg_relative_humidity,
                       avg_fuel_moisture = avg_fuel_moisture,
                       avg_fuel_temp = avg_fuel_temp)
  step_times.df <- data.frame(datetime = step_times)

  for (i in 1:length(weather_list)) {
    ## rename columns
    colnames(weather_list[[i]]) <- stations$st_id
    ## add column for time data
    weather_list[[i]] <- cbind(step_times.df, weather_list[[i]])
    ## write to csv
    write.csv(weather_list[[i]], paste0("data/weather/", fire_name, "_weather/", names(weather_list)[i], ".csv"))
  }
  if (return) return(weather_list)
}


## function to calculate distance between each point in fire perimeter and each weather station
calc_station_distances <- function(fire_name, return = FALSE) {

  fire_name <- gsub(" ", "", fire_name)
  ## read station data
  stations <- read.csv(paste0("data/weather/", fire_name, "_weather/stations.csv"), row.names = 1)
  stations$lon <- as.numeric(stations$lon)
  stations$lat <- as.numeric(stations$lat)
  ## read daily burned data
  daily_burned <- as.data.frame(fread(paste0("data/fire_progression/", fire_name, "_daily_burned.csv")))

  distances <- pointDistance(as.matrix(daily_burned[,c("x", "y")]),
                             as.matrix(stations[, c("lon","lat")]), lonlat = TRUE)
  colnames(distances) <- stations$st_id
  saveRDS(distances, paste0("data/weather/", fire_name, "_weather/station_distances.rds"))
  if(return) return(distances)
}

## function to calculate weighted average weather data for each point in fire perimeter
calc_weighted_weather_avgs <- function(fire_name, ncores = 6, return = FALSE) {
  fire_name <- gsub(" ", "", fire_name)
  ## load weather data
  weather_dir <- paste0("data/weather/", fire_name, "_weather")
  avg_air_temp <- read.csv(paste0(weather_dir, "/avg_air_temp.csv"), row.names = 1)
  max_air_temp <- read.csv(paste0(weather_dir, "/max_air_temp.csv"), row.names = 1)
  avg_wind_speed <- read.csv(paste0(weather_dir, "/avg_wind_speed.csv"), row.names = 1)
  max_wind_speed <- read.csv(paste0(weather_dir, "/max_wind_speed.csv"), row.names = 1)
  avg_relative_humidity <- read.csv(paste0(weather_dir, "/avg_relative_humidity.csv"), row.names = 1)
  avg_fuel_moisture <- read.csv(paste0(weather_dir, "/avg_fuel_moisture.csv"), row.names = 1)
  avg_fuel_temp <-  read.csv(paste0(weather_dir, "/avg_fuel_temp.csv"), row.names = 1)

  ## load elevation data
  print("Reading topography data")
  topography <- fread("data/Topography/derived/topography.csv")
  topography <- topography[, c("x", "y", "elevation")]

  ## load daily burned data
  daily_burned <- as.data.frame(fread(paste0("data/fire_progression/", fire_name, "_daily_burned.csv")))
  daily_burned$datetime <- as.character(as.chron(daily_burned$datetime))
  dts <- unique(daily_burned$datetime)[order(unique(daily_burned$datetime))]

  ## merge topography and daily burned data
  ## need to round to account for small discrepancies in raster centerpoints
  ## (shouldn't affect anything but would love to fix)
  topography$y <- round(topography$y, 5)
  topography$x <- round(topography$x, 5)
  daily_burned <- as.data.table(daily_burned)
  daily_burned$y <- round(daily_burned$y, 5)
  daily_burned$x <- round(daily_burned$x, 5)
  daily_burned <- as.data.frame(merge(daily_burned, topography, by = c("x", "y"),
                                      all.x = TRUE, all.y = FALSE))
  max_elev <- max(topography$elevation)

  ## load distances
  distances <- readRDS(paste0(weather_dir, "/station_distances.rds"))
  ## load stations
  stations <- read.csv(paste0(weather_dir, "/stations.csv"))

  ## start parallel process
  print("Spinning up cluster")
  cl <- makeSOCKcluster(ncores)
  registerDoSNOW(cl)

  iter <- length(dts)
  pb <- progress_bar$new(format = "Calculating [:bar] :current/:total (:percent) | elapsed: :elapsed",
                         total = iter)
  pb$tick(0)
  progress <- function() pb$tick()
  opts <- list(progress=progress)

  ## perform calculations
  withweather <-
    foreach(i=1:length(dts),
            .combine = rbind,
            .options.snow = opts,
            .export = "distance_elev_weighted_avg") %dopar% {
              ## compile data
              indx <- avg_air_temp$datetime == dts[i]
              indy <- 2:(ncol(avg_air_temp))
              weather_list <- list(avg_air_temp[indx, indy],
                                   max_air_temp[indx, indy],
                                   avg_wind_speed[indx, indy],
                                   max_wind_speed[indx, indy],
                                   avg_relative_humidity[indx, indy],
                                   avg_fuel_moisture[indx, indy],
                                   avg_fuel_temp[indx, indy])
              data <- daily_burned[daily_burned$datetime == dts[i],]
              dst <- distances[daily_burned$datetime == dts[i],]
              dst <- dst / 15000
              if (!is.null(nrow(dst))) {
                dst_splt <- split(dst, rep(1:nrow(dst), each = ncol(dst)))
              } else dst_splt <- dst

              out_list <- list()
              ## perform calculation for each var
              for (j in 1:length(weather_list)) {
                out_list[[j]] <- mapply(dst_splt,
                                        data[,"elevation"],
                                        FUN = distance_elev_weighted_avg,
                                        MoreArgs = list(vals = weather_list[[j]],
                                                        elevations = stations$elevation,
                                                        max_elev = max_elev))
              }

              ## collate into data frame for collection by master thread
              data.frame(x = daily_burned[daily_burned$datetime == dts[i],"x"],
                         y = daily_burned[daily_burned$datetime == dts[i],"y"],
                         datetime = rep(dts[i], times = length(out_list[[1]])),
                         avg_air_temp = out_list[[1]],
                         max_air_temp = out_list[[2]],
                         avg_wind_speed = out_list[[3]],
                         max_wind_speed = out_list[[4]],
                         avg_relative_humidity = out_list[[5]],
                         avg_fuel_moisture = out_list[[6]],
                         avg_fuel_temp = out_list[[7]])

            }

  stopCluster(cl)
  print("Finished parallel processing")

  write.csv(withweather, paste0(weather_dir, "/complete_weather.csv"))
  if (return) return(withweather)
}

## write function for
plot_stations <- function(perims) {
  stations_sf <- readRDS("data/weather/stations/stations_sf.rds")
  stations_sf <- st_transform(stations_sf, st_crs(3857))
  perims <- st_transform(perims, st_crs(3857))
  perim_buffer <- st_transform(st_buffer(st_make_valid(perims), 20000), st_crs(stations_sf))
  stations_sub <- st_intersects(stations_sf, perim_buffer)
  stations_sf <- stations_sf[as.data.frame(stations_sub)[,1],]

  ggplot() +
    basemap_gglayer(ext = st_bbox(perim_buffer), map_type = "terrain_bg") +
    scale_fill_identity() +
    new_scale_fill() +
    geom_sf(data = perims, alpha = 1, aes(fill = as.factor(FIRE_NA)), color = NA) +
    geom_sf(data = stations_sf, color = "black", size = 3) +
    scale_fill_manual(values = c("#fb9a99","#33a02c", "#fdbf6f", "#1f78b4", "#e31a1c")) +
    theme_bw() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    ggtitle("RAWS station positions relative to fire perimeters") +
    theme(axis.title = element_blank(),
          legend.title = element_blank(),
          legend.position = "bottom")
}

##---------------------------------------------------------------
## 2. Body
##---------------------------------------------------------------

## processing
perims <- st_read("data/perimeters/all_perimiters.shp")

for (fire in fire_list) {
  calc_time_avgd_weather(fire, perims = perims)
  print("Calculating distance matrix")
  calc_station_distances(fire)
  calc_weighted_weather_avgs(fire, ncores = 4)
}

## plotting
station_locations <- plot_stations(perims = perims)
station_locations

test <- read.csv("data/weather/northcomplex_weather/complete_weather.csv", row.names = 1)
ggplot(test, aes(x = x, y = y, color = avg_air_temp)) +
  geom_tile() +
  theme_bw()

test2 <- read.csv("data/fire_progression//northcomplex_daily_burned.csv")
ggplot(test2, aes(x = x, y = y, color = datetime)) +
  geom_tile() +
  theme_bw()

