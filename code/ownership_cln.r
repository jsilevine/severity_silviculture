##---------------------------------------------------------------
## ownership_cln.r
## BY: JACOB LEVINE - jacoblevine@princeton.edu
##---------------------------------------------------------------


##---------------------------------------------------------------
## 0. head
##---------------------------------------------------------------

##source("~/Documents/Science/severity_and_silviculture/code/variable_def.r")
setwd("../")

library(sf)
library(raster)
library(fasterize)
library(ggplot2)
library(ggnewscale)
library(ggspatial)
library(mapdata)
library(ggmap)

##---------------------------------------------------------------
## 1. Body
##---------------------------------------------------------------

ownership_nonprivate <- st_read("data/ownership/ownership15_1.shp")
ownership_nonprivate <- st_transform(ownership_nonprivate, crs = st_crs(4326))
reclass <- data.frame(orig = unique(ownership_nonprivate$Own_Group),
                      mapping = c(1,1,1,1,1,1,1,1,1,1,2,1,1,1))

for (i in 1:nrow(ownership_nonprivate)) {
  ownership_nonprivate$own_int[i] <- reclass[reclass$orig == ownership_nonprivate$Own_Group[i], "mapping"]
}

ownership_private <- st_read("data/ownership//Forest_Industry_Owners_18_2.shp")
ownership_private <- st_transform(ownership_private, crs = st_crs(4326))
ownership_private$own_int <- 3

## resample to daily progression grid
togrid <- readRDS("data/fire_progression//dixie_daily_burned.rds")
togrid <- merge(togrid, readRDS("data/fire_progression//northcomplex_daily_burned.rds"))
togrid <- merge(togrid, readRDS("data/fire_progression//sheep_daily_burned.rds"))
togrid <- merge(togrid, readRDS("data/fire_progression//walker_daily_burned.rds"))
togrid <- merge(togrid, readRDS("data/fire_progression//sugar_daily_burned.rds"))

## convert to raster, merge,
ownership_nonprivate.rast <- fasterize(ownership_nonprivate, togrid, field = "own_int")
ownership_private.rast <- fasterize(ownership_private, togrid, field = "own_int")
ownership <- merge(ownership_nonprivate.rast, ownership_private.rast)
ownership[is.na(ownership[])] <- 2
perims <- st_read("data/perimeters/all_perimiters.shp")
tomask <- st_buffer(perims, 1000)
ownership <- mask(ownership, as_Spatial(perims))
saveRDS(ownership, "data/ownership/ownership.rds")

## convert to dataframe and write to csv
ownership.df <- as.data.frame(ownership, xy = TRUE, na.rm = TRUE)
colnames(ownership.df)[3] <- "own_int"
ownership.df[ownership.df$own_int == 1, "own_type"] <- "Federal"
ownership.df[ownership.df$own_int == 2, "own_type"] <- "Other"
ownership.df[ownership.df$own_int == 3, "own_type"] <- "Private Industrial"
write.csv(ownership.df, "data/ownership/ownership.csv")

## takes a long time to render

ggplot(ownership.df) +
  geom_tile(data = ownership.df, aes(x = x, y = y, fill = own_type)) +
  coord_sf(expand = TRUE) +
  scale_fill_manual(values = c("#bdbdbd", "#f0f0f0", "#636363")) +
  new_scale_fill() +
  geom_sf(data = perims, aes(color = FIRE_NA),
          fill = NA, linewidth = 0.65) +
  scale_color_brewer(type = "qual", palette = 3) +
  ##geom_sf_label(data = perims, aes(label = FIRE_NA), size = 3) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.title = element_blank(), legend.title = element_blank(),
        panel.grid = element_blank())


ggsave()
