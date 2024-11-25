## define function to snap to raster

round_to_nearest <- function(time, nearest = 8) {
  # Convert chron time to numeric hours
  hours <- as.numeric(time) * 24  # Convert to hours
  # Round to nearest 8
  rounded_hours <- round(hours / nearest) * nearest
  # Convert back to chron
  chron_time <- chron(time = rounded_hours / 24)  # Convert back to chron
  return(chron_time)
}

snap_to_template <- function(raster, template, method = "bilinear", inverse = FALSE) {

  if (as.character(crs(raster)) != as.character(crs(template))) {
    raster <- terra::project(raster, crs(template))
  }

  raster <- terra::resample(raster, template, method)
  ext <- terra::ext(template)
  raster <- terra::crop(raster, ext)
  raster <- terra::mask(raster, template, inverse = inverse)

  return(raster)

}
