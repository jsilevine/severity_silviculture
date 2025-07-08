
library(rvest)
library(sf)
library(rmapshaper)
library(fasterize)
library(terra)
library(ggplot2)

source("code/utility/utility_functions.r")

##---------------------------------------------------------------
## Severity
##---------------------------------------------------------------

template <- rast("data/templates/isforest.tif")

sevfiles <- list.files("data/cbi")

for (i in 1:length(sevfiles)) {

  fire_name <- strsplit(tolower(strsplit(sevfiles[i], "_CBI_")[[1]][1]), "_")[[1]][1]

  severity <- rast(paste0("data/cbi/", sevfiles[i]))
  severity <- snap_to_template(severity, template)
  severity_df <- as.data.frame(severity, xy = TRUE, na.rm = TRUE)
  writeRaster(severity, paste0("data/cbi/", fire_name, "_processed.tif"), overwrite = TRUE)
  if (i == 1) {
    severity_full_df <- severity_df
  } else {
    severity_full_df <- rbind(severity_full_df, severity_df)
  }
  write.csv(severity_df, paste0("data/cbi/", fire_name, "_processed.csv"))

}

dixie_cbi <- rast("data/cbi/SUGAR_2021_CBI_bc.tif")
dixie_cbi <- snap_to_template(dixie_cbi, template)
plot(dixie_cbi)

plot(dixie_cbi > 2.25)

dixie_dnbr <- rast("data/old/sugar_rdnbr.tif")
dixie_dnbr <- snap_to_template(dixie_dnbr, template)
plot(dixie_dnbr)
plot(dixie_dnbr > 641)


sevfiles <- list.files("data/dnbr")

for (i in 1:length(sevfiles)) {

  fire_name <- strsplit(tolower(strsplit(sevfiles[i], "_CBI_")[[1]][1]), "_")[[1]][1]

  severity <- rast(paste0("data/dnbr/", sevfiles[i]))
  severity <- snap_to_template(severity, template)
  severity_df <- as.data.frame(severity, xy = TRUE, na.rm = TRUE)
  colnames(severity_df)[3] <- "dnbr"
  writeRaster(severity, paste0("data/dnbr/", fire_name, "_processed.tif"), overwrite = TRUE)
  if (i == 1) {
    severity_full_df <- severity_df
  } else {
    severity_full_df <- rbind(severity_full_df, severity_df)
  }
  write.csv(severity_df, paste0("data/dnbr/", fire_name, "_processed.csv"))

}


sevfiles <- list.files("data/rdnbr")

for (i in 1:length(sevfiles)) {

  fire_name <- strsplit(tolower(strsplit(sevfiles[i], "_CBI_")[[1]][1]), "_")[[1]][1]

  severity <- rast(paste0("data/rdnbr/", sevfiles[i]))
  severity <- snap_to_template(severity, template)
  severity_df <- as.data.frame(severity, xy = TRUE, na.rm = TRUE)
  colnames(severity_df)[3] <- "rdnbr"
  writeRaster(severity, paste0("data/rdnbr/", fire_name, "_processed.tif"), overwrite = TRUE)
  if (i == 1) {
    severity_full_df <- severity_df
  } else {
    severity_full_df <- rbind(severity_full_df, severity_df)
  }
  write.csv(severity_df, paste0("data/rdnbr/", fire_name, "_processed.csv"))

}

sevfiles <- list.files("data/old")

for (i in 1:length(sevfiles)) {

  fire_name <- strsplit(tolower(strsplit(sevfiles[i], "_CBI_")[[1]][1]), "_")[[1]][1]

  severity <- rast(paste0("data/old/", sevfiles[i]))
  severity <- snap_to_template(severity, template)
  severity_df <- as.data.frame(severity, xy = TRUE, na.rm = TRUE)
  colnames(severity_df)[3] <- "rdnbr"
  writeRaster(severity, paste0("data/old/", fire_name, "_processed.tif"), overwrite = TRUE)
  if (i == 1) {
    severity_full_df <- severity_df
  } else {
    severity_full_df <- rbind(severity_full_df, severity_df)
  }
  write.csv(severity_df, paste0("data/old/", fire_name, "_processed.csv"))

}



