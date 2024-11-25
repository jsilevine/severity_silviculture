library(terra)
library(fasterize)
library(sf)

sevfiles <- list.files("data/rdnbr/old")
severity <- rast(paste0("data/rdnbr/old/", sevfiles[1]))
for (i in 2:length(sevfiles)) {
  severity <- merge(severity, rast(paste0("data/rdnbr/old/", sevfiles[i])))
}

cc_full <- rast("data/canopy_cover/CanopyCover_2-None_PNF_1_BS30.tif")

perims <- st_read("data/perimeters/all_perimiters.shp")
perims_rast <- rast(fasterize(perims, raster(severity)))

template <- severity

## snap to template, create template raster
cc_full <- terra::project(cc_full, template)
cc_full <- resample(cc_full, template, "bilinear")
ext <- ext(template)
cc_full <- crop(cc_full, ext)
cc_full <- mask(cc_full, template)
template <- mask(template, cc_full)
template <- mask(template, perims_rast)
template[!is.na(template[])] <- 1

plot(template)

## save template raster, for use in other cleaning scripts.
saveRDS(raster(template), "data/templates/raster_template.rds")
writeRaster(template, "data/templates/raster_template.tif")
