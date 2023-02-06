
get_ftp_filename <- function(url) {

  lnks <- read_html(url) |>
    html_nodes("a") |>
    html_text(trim = T)
  if (length(lnks) == 1) {
    ret <- lnks
  } else ret <- lnks[grepl("gdb.zip", lnks) | grepl("GDB.zip", lnks) | grepl("gbd.zip", lnks) | (grepl("vent", lnks) & grepl(".zip", lnks))] ## some typos, not a well maintained dataset
  if (length(ret > 1)) {
    ret <- ret[grepl(".gdb", ret)]
    if (length(ret > 1)) ret <- ret[1]
  }
  return(ret)
}

download_ftp_data <- function(base_url, start_date, end_date, destfolder, fire_name) {

  dates <- seq(as.Date(start_date), as.Date(end_date), by = 1)
  dates <- gsub("-", "", dates)

  for (date in dates) {
    folder <- paste0(base_url, date, "/")
    dfile <- paste0(destfolder, "/", fire_name, "_", date, ".gdb.zip")
    tryCatch({
      file <- get_ftp_filename(folder)
      download.file(paste0(folder, file), destfile = dfile)
      unzip(dfile, exdir = "tmp")
      file.rename(paste0("tmp/", list.files("tmp/")), gsub(".zip", "", dfile))
      if (length(list.files(gsub(".zip", "", dfile))) == 1) {
        file.rename(paste0(gsub(".zip", "", dfile), "/", list.files(gsub(".zip", "", dfile))),
                    paste0(gsub(".zip", "", dfile), "1"))
        unlink(gsub(".zip", "", dfile), recursive = TRUE)
        file.rename(paste0(gsub(".zip", "", dfile), "1"), gsub(".zip", "", dfile))
      }
      unlink("tmp", recursive = TRUE)
      file.remove(dfile)
      closeAllConnections()
    }, error = function(e){print(paste0("missing file for: ", date))})
  }

  ## remove empty folders:
  for (f in list.files(destfolder)) {
    if (length(list.files(paste0(destfolder, f))) == 0) {
      unlink(paste0(destfolder, f), recursive = TRUE)
    }
  }

}

Dixie_start <- "2021/07/17"
Dixie_end <- "2021/09/12"
Dixie_base_url <- "https://ftp.wildfire.gov/public/incident_specific_data/calif_n/!CALFIRE/2021_Incidents/CA-BTU-009205_Dixie/GIS/IncidentData/"

download_ftp_data(Dixie_base_url, Dixie_start, Dixie_end, "data/Dixie", "Dixie")

Walker_start <- "2019/09/06"
Walker_end <- "2019/10/03"
Walker_base_url <- "https://ftp.wildfire.gov/public/incident_specific_data/calif_n/2019_FEDERAL_Incidents/CA-PNF-001324_Walker/GIS/IncidentData/"

download_ftp_data(Walker_base_url, Walker_start, Walker_end, "data/Walker", "Walker")

NorthComplex_start <- "2020/08/19"
NorthComplex_end <- "2020/10/29"
NorthComplex_base_url <- "https://ftp.wildfire.gov/public/incident_specific_data/calif_n/2020_FEDERAL_Incidents/CA-PNF-001308_PNF_North_Complex/GIS/IncidentData/"

download_ftp_data(NorthComplex_base_url, NorthComplex_start, NorthComplex_end, "data/NorthComplex", "NorthComplex")

dir <- "data/NorthComplex/"
fire_name <- c("Sheep", "Claremont", "Bear", "Claremont-Bear")

process_ftp_data <- function(dir, fire_name) {

  files <- list.files(dir)
  files <- files[order(files)] ## ensure in correct order

   for (f in list.files(dir)) {
    if (length(list.files(paste0(dir, f))) == 0) {
      unlink(paste0(dir, f), recursive = TRUE)
    }
  }

  for (f in 2:length(files)) {
    print(f)
    curr <- readOGR(paste0(dir, files[f]), "Event_Polygon")
    curr <- curr[!is.na(curr$IncidentName),]
    curr <- st_make_valid(st_as_sf(
      curr[curr$IncidentName %in% fire_name &
           curr$FeatureCategory == "Wildfire Daily Fire Perimeter",
           c("DateCurrent")]))
    curr <- st_transform(curr, crs = 4326)
    prev <- readOGR(paste0(dir, files[f-1]), "Event_Polygon")
    if (fire_name == "Dixie" & f == 2) { ## some weirdness here
      prev <- st_transform(st_as_sf(prev[2, c("DateCurrent")]), crs = crs(curr))
    } else {
      prev <- prev[!is.na(prev$IncidentName),]
      prev <- st_make_valid(st_as_sf(
        prev[prev$IncidentName %in% fire_name &
             prev$FeatureCategory == "Wildfire Daily Fire Perimeter",
             c("DateCurrent")]))
    }
    prev <- st_transform(prev, crs = 4326)
    if (!(nrow(prev) == 0 | nrow(curr) == 0)) {
      dif <- ms_erase(curr, prev)[1,]
      if (!is.na(st_bbox(dif)[1])) {
        date <- strsplit(strsplit(files[f], "_")[[1]][2], ".g")[[1]][1]
        dif$file_date <- date
        if (!exists("outdata")) {
          prev$file_date <- date
          prev <- prev[, colnames(dif)]
          outdata <- rbind(prev, dif)
        } else {
          outdata <- rbind(outdata, dif)
          print(outdata)
        }
      rm(dif)
    }
    }
  }

  return(outdata)
}

dixie_data <- process_ftp_data("data/Dixie/", "Dixie")
saveRDS(dixie_data, "data/full_fires/dixie_data.rds")
dixie_data <- readRDS("data/full_fires/dixie_data.rds")

northcomplex_data <- process_ftp_data("data/NorthComplex/",  c("Sheep", "Claremont", "Bear", "Claremont-Bear"))
saveRDS(northcomplex_data, "data/full_fires/northcomplex_data.rds")
northcomplex_data <- readRDS("data/full_fires/northcomplex_data.rds")

## load severity data
sheep_rdnbr <- raster("data/rdnbr/sheep_rdnbr.tif")
northcomplex_rdnbr <- raster("data/rdnbr/northcomplex_rdnbr.tif")
northcomplex_rdnbr <- merge(northcomplex_rdnbr, sheep_rdnbr)
writeRaster(northcomplex_rdnbr, "data/rdnbr/northcomplex_rdnbr.tif", overwrite = TRUE)

pull_daily_severity <- function(fire_name, proj = 4326, perim_data) {
  rdnbr <- raster(paste0("data/rdnbr/", fire_name, "_rdnbr.tif"))
  perim_data[1,] <- st_transform(perim_data[1,], crs = crs(rdnbr))
  mask <- fasterize(perim_data[1,], rdnbr)
  sev <- mask(rdnbr, mask)
  sev <- as.data.frame(as(sev, "SpatialPixelsDataFrame"))
  sev$time <- 1
  for (i in 2:nrow(perim_data)) {
    print(paste0("finished ", i, "of ", nrow(perim_data)))
    perim_data[i,] <- st_transform(perim_data[i,], crs = crs(rdnbr))
    mask <- fasterize(perim_data[i,], rdnbr)
    if (all(is.na(values(mask)))) {}
    else {
    nsev <- mask(rdnbr, mask)
    nsev <- as.data.frame(as(nsev, "SpatialPixelsDataFrame"))
    nsev$time <- i
    sev <- rbind(sev, nsev)
    }
  }
  return(sev)
}

dixie_daily_severity <- pull_daily_severity("dixie", perim_data = dixie_data)
saveRDS(dixie_daily_severity, "data/daily_severity/dixie_daily_severity.RDS")
northcomplex_daily_severity <- pull_daily_severity("northcomplex", perim_data = northcomplex_data)
saveRDS(northcomplex_daily_severity, "data/daily_severity/northcomplex_daily_severity.RDS")

dixie_daily_severity <- readRDS("data/daily_severity/dixie_daily_severity.RDS")
