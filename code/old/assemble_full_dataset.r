## assemble full datasets::
library(raster)
library(data.table)
library(chron)
library(ggplot2)
library(ggnewscale)
library(cowplot)
library(parallel)

## Need to combine:
## 1. Topography
## 2. canopy_cover (ladder_fuels)
## 3. cbi (severity)
## 4. daily progression data
## 5. weather
## 6. forest structure (lidar)

## first read in template
template <- readRDS("data/templates/raster_template.rds")

## define function to snap to raster
snap_to_template <- function(raster, template) {

  if (as.character(crs(raster)) != as.character(crs(template))) {
    raster <- projectRaster(raster, crs = crs(template))
  }

  raster <- resample(raster, template, "bilinear")
  ext <- extent(template)
  raster <- crop(raster, ext)
  raster <- mask(raster, template)

  return(raster)

}

##---------------------------------------------------------------
## 1. topography
##---------------------------------------------------------------
topography <- fread("data/Topography/derived/topography.csv")
topography <- topography[,.(x, y, elevation, slope, tpi, heat_load)] ## select columns of interest
topography <- rasterFromXYZ(topography, crs = crs(template))
topography <- snap_to_template(topography, template) ## mask and resample
topography <- data.table(as.data.frame(topography, xy = TRUE, na.rm = TRUE)) ## convert back to data.table
topography$x_j <- round(topography$x, 5) ## round to avoid rounding errors in matching
topography$y_j <- round(topography$y, 5) ## round to avoid rounding errors in matching

##---------------------------------------------------------------
## 2. canopy cover
##---------------------------------------------------------------
canopy_cover <- fread("data/canopy_cover/processed/cc_complete.csv")
canopy_cover <- canopy_cover[, .(x,y,cc2_8, cc8_16, cc16_32, cc32_n)]
canopy_cover$x_j <- round(canopy_cover$x, 5)
canopy_cover$y_j <- round(canopy_cover$y, 5)

## join topography and canopy cover
full_data <- canopy_cover[topography, on = .(x_j,y_j), nomatch = NULL]
full_data <- full_data[,.(x,y,cc2_8,cc8_16,cc16_32,cc32_n,x_j,y_j,elevation,slope,tpi,heat_load)]

rm(topography, canopy_cover)

## calculate evenness as measure of ladder fuels.
calc_em <- function(i, data = full_data) {
  cc <- as.numeric(data[i, .(cc2_8, cc8_16, cc16_32, cc32_n)])
  cc <- cc[cc > 0.0]
  cc_p <- cc / sum(cc)
  iD <- sum(cc_p^2)
  D <- 1 / iD
  E <- D / length(cc)
  em <- E * mean(cc)
}

i <- 1:nrow(full_data)
em <- mclapply(i, FUN = calc_em, mc.cores = 8) ## takes a little while
full_data[,"em"] <- unlist(em)
rm(em)


##---------------------------------------------------------------
## 3. severity
##---------------------------------------------------------------
dixie_severity <- raster("data/cbi/DIXIE_2021_CBI_bc.tif")
walker_severity <- raster("data/cbi/WALKER_2019_CBI_bc.tif")
north_complex_severity <- raster("data/cbi/NORTH_COMPLEX_2020_CBI_bc.tif")
sheep_severity <- raster("data/cbi/SHEEP_2020_CBI_bc.tif")
sugar_severity <- raster("data/cbi/SUGAR_2021_CBI_bc.tif")

dixie_severity <- snap_to_template(dixie_severity, template)
walker_severity <- snap_to_template(walker_severity, template)
north_complex_severity <- snap_to_template(north_complex_severity, template)
sheep_severity <- snap_to_template(sheep_severity, template)
sugar_severity <- snap_to_template(sugar_severity, template)

dixie_severity <- data.table(as.data.frame(dixie_severity, xy = TRUE, na.rm = TRUE))
walker_severity <- data.table(as.data.frame(walker_severity, xy = TRUE, na.rm = TRUE))
north_complex_severity <- data.table(as.data.frame(north_complex_severity, xy = TRUE, na.rm = TRUE))
sheep_severity <- data.table(as.data.frame(sheep_severity, xy = TRUE, na.rm = TRUE))
sugar_severity <- data.table(as.data.frame(sugar_severity, xy = TRUE, na.rm = TRUE))

dixie_severity$fire_name <- "dixie"
north_complex_severity$fire_name <- "north_complex"
walker_severity$fire_name <- "walker"
sheep_severity$fire_name <- "sheep"
sugar_severity$fire_name <- "sugar"

full_severity <- rbind(dixie_severity, walker_severity, north_complex_severity, sheep_severity, sugar_severity)
full_severity$x_j <- round(full_severity$x, 5)
full_severity$y_j <- round(full_severity$y, 5)

full_data <- full_data[full_severity, on = .(x_j, y_j)]
full_data <- full_data[, .(x, y, CBI_bc, fire_name, cc2_8, cc8_16, cc16_32, cc32_n, em, elevation, slope, tpi, heat_load, x_j, y_j)]

## high severity cutoffs
full_data$hs <- 0
full_data[CBI_bc > 2.25, "hs"] <- 1
full_data$hs_class <- "not_high"
full_data[CBI_bc > 2.25, "hs_class"] <- "high"

rm(full_severity)

##---------------------------------------------------------------
## 4. Daily progressions
##---------------------------------------------------------------

dixie_daily_burned <- fread("data/fire_progression/dixie_daily_burned.csv")
dixie_daily_burned <- dixie_daily_burned[,.(x,y,datetime)]
dixie_daily_burned$datetime_formatted <- chron(dixie_daily_burned$datetime)

walker_daily_burned <- fread("data/fire_progression/walker_daily_burned.csv")
walker_daily_burned <- walker_daily_burned[,.(x,y,datetime)]
walker_daily_burned$datetime_formatted <- chron(walker_daily_burned$datetime)

north_complex_daily_burned <- fread("data/fire_progression/northcomplex_daily_burned.csv")
north_complex_daily_burned <- north_complex_daily_burned[,.(x,y,datetime)]
north_complex_daily_burned$datetime_formatted <- chron(north_complex_daily_burned$datetime)

sheep_daily_burned <- fread("data/fire_progression/sheep_daily_burned.csv")
sheep_daily_burned <- sheep_daily_burned[,.(x,y,datetime)]
sheep_daily_burned$datetime_formatted <- chron(sheep_daily_burned$datetime)

sugar_daily_burned <- fread("data/fire_progression/sugar_daily_burned.csv")
sugar_daily_burned <- sugar_daily_burned[,.(x,y,datetime)]
sugar_daily_burned$datetime_formatted <- chron(sugar_daily_burned$datetime)

full_daily_burned <- rbind(dixie_daily_burned, walker_daily_burned, north_complex_daily_burned, sheep_daily_burned, sugar_daily_burned)

full_daily_burned$x_j <- round(full_daily_burned$x, 5)
full_daily_burned$y_j <- round(full_daily_burned$y, 5)
full_data <- full_daily_burned[full_data, on = .(x_j, y_j)]
full_data <- full_data[, .(x,y,datetime_formatted,x_j,y_j,CBI_bc,fire_name,cc2_8,cc8_16,cc16_32,cc32_n,em,elevation,slope,tpi,heat_load,hs)]
rm(full_daily_burned)

dwm <- function(i, data1, data2) {

  d <- as.matrix(dist(as.matrix(rbind(data2[i,.(x,y)], data1[,.(x,y)])), method = "euclidean"))
  dwm <- sum((1/d[2:(nrow(data1)+1),1]) * data1[,.(CBI_bc)]) / (sum(1/d[2:(nrow(data1)+1),1]))
  return(dwm)

}


r <- 0.0002694946 * (1000 / 30) ## 1000m radius in decimal degrees
## Calculate weighted average severity on previous timestep
for (fire in unique(full_data$fire_name)) {
  dt <- unique(full_data[full_data$fire_name == fire, .(datetime_formatted)])
  dt  <- dt[order(dt)][[1]]
  for (i in 1:length(dt)) {
    if (i == 1) {
      full_data[full_data$fire_name == fire & full_data$datetime_formatted == dt[i], "prev_sev"] <- 0
    } else {

      data1 <- full_data[full_data$fire_name == fire & full_data$datetime_formatted == dt[i-1], ]
      data2 <- full_data[full_data$fire_name == fire & full_data$datetime_formatted == dt[i], ]

      dwm <- function(j) {
        x_c <- data2[j,x]
        y_c <- data2[j,y]
        data1_sub <-  data1[x > x_c - r & x < x_c + r &
                            y > y_c - r & y < y_c + r]
        if (nrow(data1_sub) > 0) {
          d <- as.matrix(dist(as.matrix(rbind(data2[j,.(x,y)], data1_sub[,.(x,y)])), method = "euclidean"))
          dwm <- sum((1/d[2:(nrow(data1_sub)+1),1]) * data1_sub[,.(CBI_bc)]) / (sum(1/d[2:(nrow(data1_sub)+1),1]))
          return(dwm)
        } else {
          return(0)
        }
      }

      full_data[full_data$fire_name == fire & full_data$datetime_formatted == dt[i], "prev_sev"] <-
        unlist(mclapply(1:nrow(full_data[full_data$fire_name == fire & full_data$datetime_formatted == dt[i], ]),
               dwm,
               mc.cores = 8))
    }
    print(paste0("completed: ", fire, ", iteration ", i))
  }
}

##---------------------------------------------------------------
## 5. weather
##---------------------------------------------------------------

dixie_weather <- fread("data/weather/dixie_weather/complete_weather.csv")
dixie_weather <- dixie_weather[,.(x,y,datetime,avg_air_temp,max_air_temp,avg_wind_speed,max_wind_speed,avg_relative_humidity,
                                  avg_fuel_moisture,avg_fuel_temp)]

walker_weather <- fread("data/weather/walker_weather/complete_weather.csv")
walker_weather <- walker_weather[,.(x,y,datetime,avg_air_temp,max_air_temp,avg_wind_speed,max_wind_speed,avg_relative_humidity,
                                  avg_fuel_moisture,avg_fuel_temp)]

north_complex_weather <- fread("data/weather/northcomplex_weather/complete_weather.csv")
north_complex_weather <- north_complex_weather[,.(x,y,datetime,avg_air_temp,max_air_temp,avg_wind_speed,max_wind_speed,avg_relative_humidity,
                                  avg_fuel_moisture,avg_fuel_temp)]

sheep_weather <- fread("data/weather/sheep_weather/complete_weather.csv")
sheep_weather <- sheep_weather[,.(x,y,datetime,avg_air_temp,max_air_temp,avg_wind_speed,max_wind_speed,avg_relative_humidity,
                                  avg_fuel_moisture,avg_fuel_temp)]

sugar_weather <- fread("data/weather/sugar_weather/complete_weather.csv")
sugar_weather <- sugar_weather[,.(x,y,datetime,avg_air_temp,max_air_temp,avg_wind_speed,max_wind_speed,avg_relative_humidity,
                                  avg_fuel_moisture,avg_fuel_temp)]

full_weather <- rbind(dixie_weather, walker_weather, north_complex_weather, sheep_weather, sugar_weather)


calc_es <- function(temp){
  es <- 6.11 * exp((2.5e6 / 461) * (1 / 273 - 1 / (273 + temp)))
  return(es)
}

calc_vpd <- function(rh, temp){
  es <- calc_es(temp)
  vpd <- ((100 - rh) / 100) * es
  return(vpd)
}

full_weather$vpd <- unlist(mapply(calc_vpd, full_weather$avg_fuel_moisture, full_weather$avg_air_temp))
full_weather$hdw <- full_weather$vpd * full_weather$avg_wind_speed

full_weather$x_j <- round(full_weather$x, 5)
full_weather$y_j <- round(full_weather$y, 5)
fwrite(full_weather, "data/weather/full_weather.csv")

full_data <- full_weather[full_data, on = .(x_j, y_j)]
full_data <- full_data[,.(x,y,hs,CBI_bc,fire_name,datetime_formatted,avg_air_temp,max_air_temp,avg_wind_speed,max_wind_speed,
                          avg_relative_humidity,vpd,hdw,
                          avg_fuel_moisture,avg_fuel_temp,cc2_8,cc8_16,cc16_32,cc32_n,em,prev_sev,
                          elevation,slope,tpi,heat_load,x_j,y_j)]

##---------------------------------------------------------------
## 6. forest structure
##---------------------------------------------------------------

## gap data 30m
gap_data_30 <- fread("data/processed/gap_data_30m_full.csv")
gap_data_30 <- gap_data_30[,.(x,y,mean_area,median_area,sd_area,percent_open,mean_frac)]

gap_data_raster_30 <- rasterFromXYZ(gap_data_30, crs = 3310)
gap_data_raster_30 <- snap_to_template(gap_data_raster_30, template)
gap_data_30 <- data.table(as.data.frame(gap_data_raster_30, xy = TRUE, na.rm = TRUE))
gap_data_30$x_j <- round(gap_data_30$x, 5)
gap_data_30$y_j <- round(gap_data_30$y, 5)
setnames(gap_data_30,
         c("mean_area","median_area","sd_area","percent_open","mean_frac"),
         c("mean_area_30","median_area_30","sd_area_30","percent_open_30","mean_frac_30"))

full_data <- gap_data_30[full_data, on = .(x_j, y_j)]
rm(gap_data_30)

## gap data 180m
gap_data_180 <- fread("data/processed/gap_data_180m_full.csv")
gap_data_180 <- gap_data_180[,.(x,y,mean_area,median_area,sd_area,percent_open,mean_frac)]

gap_data_raster_180 <- rasterFromXYZ(gap_data_180, crs = 3310)
gap_data_raster_180 <- snap_to_template(gap_data_raster_180, template)
gap_data_180 <- data.table(as.data.frame(gap_data_raster_180, xy = TRUE, na.rm = TRUE))
gap_data_180$x_j <- round(gap_data_180$x, 5)
gap_data_180$y_j <- round(gap_data_180$y, 5)
setnames(gap_data_180,
         c("mean_area","median_area","sd_area","percent_open","mean_frac"),
         c("mean_area_180","median_area_180","sd_area_180","percent_open_180","mean_frac_180"))

full_data <- gap_data_180[full_data, on = .(x_j, y_j)]
rm(gap_data_180)

## tree data 30m
tree_data_30 <- fread("data/processed/tree_data_30m_full.csv")

tree_data_30 <- tree_data_30[,.(x,y,clust,mean_dens,mean_ht)]
tree_data_raster_30 <- rasterFromXYZ(tree_data_30, crs = 3310)
tree_data_raster_30 <- snap_to_template(tree_data_raster_30, template)
tree_data_30 <- data.table(as.data.frame(tree_data_raster_30, xy = TRUE, na.rm = TRUE))
tree_data_30$x_j <- round(tree_data_30$x, 5)
tree_data_30$y_j <- round(tree_data_30$y, 5)
setnames(tree_data_30,
         c("clust","mean_dens", "mean_ht"),
         c("clust_30","mean_dens_30", "mean_ht_30"))

full_data <- tree_data_30[full_data, on = .(x_j, y_j)]
rm(tree_data_30)

## tree data 180m
tree_data_180 <- fread("data/processed/tree_data_180m_full.csv")

tree_data_180 <- tree_data_180[,.(x,y,clust,mean_dens,mean_ht)]
tree_data_raster_180 <- rasterFromXYZ(tree_data_180, crs = 3310)
tree_data_raster_180 <- snap_to_template(tree_data_raster_180, template)
tree_data_180 <- data.table(as.data.frame(tree_data_raster_180, xy = TRUE, na.rm = TRUE))
tree_data_180$x_j <- round(tree_data_180$x, 5)
tree_data_180$y_j <- round(tree_data_180$y, 5)
setnames(tree_data_180,
         c("clust","mean_dens", "mean_ht"),
         c("clust_180","mean_dens_180", "mean_ht_180"))

full_data <- tree_data_180[full_data, on = .(x_j, y_j)]
rm(tree_data_180)

##---------------------------------------------------------------
## 7. ownership
##---------------------------------------------------------------
ownership <- fread("data/ownership/ownership.csv")
ownership <- ownership[,.(x,y,own_int,own_type)]
ownership$x_j <- round(ownership$x, 5)
ownership$y_j <- round(ownership$y, 5)

full_data <- ownership[full_data, on = .(x_j, y_j)]
full_data <- full_data[, .(x,y,hs,CBI_bc,fire_name,datetime_formatted,own_int,own_type,clust_30,mean_dens_30,mean_ht_30,mean_area_30,median_area_30,
                           sd_area_30,percent_open_30,mean_frac_30,clust_180,mean_dens_180,mean_ht_180,mean_area_180,median_area_180,
                           sd_area_180,percent_open_180,mean_frac_180,avg_air_temp,max_air_temp,avg_wind_speed,max_wind_speed,avg_relative_humidity,vpd,hdw,
                           avg_fuel_moisture,avg_fuel_temp,elevation,slope,tpi,heat_load,cc2_8,cc8_16,cc16_32,cc32_n,em,prev_sev)]

full_data$hs_class <- "not_high"
full_data[CBI_bc > 2.25, "hs_class"] <- "high"


coords <- st_as_sf(full_data[,c(1,2,3)], coords = c("x", "y"), crs = st_crs(4326))
coords <- st_transform(coords, crs = 3310)
xy <- st_coordinates(coords)
full_data$x <- xy[,1]
full_data$y <- xy[,2]

fwrite(full_data, "data/complete_data.csv")

full_data <- fread("data/complete_data.csv")


##---------------------------------------------------------------
## 8. check out the data
##---------------------------------------------------------------
full_data <- fread("data/complete_data.csv")


summary(lm(CBI_bc ~ hdw * mean_dens_30, data = full_data))



full_data[is.na(percent_open_30) & !is.na(clust_30),"percent_open_30"] <- 1
plot(rasterFromXYZ(full_data[,.(x,y,percent_open_30)]))
full_data[is.na(percent_open_180) & !is.na(clust_180),"percent_open_180"] <- 1
plot(rasterFromXYZ(full_data[,.(x,y,percent_open_180)]))

plot(rasterFromXYZ(full_data[,.(x,y,own_int)]))


p1 <- ggplot() +
  geom_tile(full_data = full_data[full_data$clust < 2 & full_data$x < -121 & full_data$y < 39.8,],
            aes(x = x, y = y, fill = clust)) +
  scale_fill_distiller(type = "div", palette = "Spectral") +
  ggtitle("stem clustering") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()

p2 <- ggplot() +
  geom_tile(full_data = full_data[full_data$clust < 2 & full_data$x < -121 & full_data$y < 39.8,],
            aes(x = x, y = y, fill = mean_dens)) +
  scale_fill_distiller(type = "div", palette = "Spectral") +
  ggtitle("stem density") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()

p3 <- ggplot() +
  geom_tile(full_data = full_data[full_data$clust < 2 & full_data$x < -121 & full_data$y < 39.8,],
            aes(x = x, y = y, fill = CBI_bc)) +
  scale_fill_distiller(type = "div", palette = "Spectral") +
  ggtitle("mean gap size") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()
p3

p4 <- ggplot() +
  geom_tile(full_data = full_data[full_data$fire_name == "dixie",],
            aes(x = x, y = y, fill = CBI_bc)) +
  scale_fill_distiller(type = "div", palette = "Spectral") +
  ggtitle("fire severity") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw()
p4

plot_grid(p1, p2, p3, p4, nrow = 2)

ggplot(full_full_data) +
  geom_density(aes(x = CBI_bc, fill = own_type), alpha = 0.4) +
  theme_bw()

mean(as.numeric(full_full_data[full_full_data$own_type == "Private Industrial", ][["CBI_bc"]]), na.rm = TRUE)
mean(as.numeric(full_full_data[full_full_data$own_type == "Federal", ][["CBI_bc"]]), na.rm = TRUE)
