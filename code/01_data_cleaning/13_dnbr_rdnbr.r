library(data.table)
library(speedglm)
library(sf)
library(gstat)
library(parallel)

data <- fread("data/complete_data.csv")

dnbr_files <- list.files(path = "data/dnbr", pattern = "\\.csv$", full.names = TRUE)
data[, dnbr := NA_real_]

for (i in seq_along(dnbr_files)) {
  sev <- fread(dnbr_files[i])
  f <- strsplit(strsplit(tolower(strsplit(dnbr_files[i], "_processed")[[1]][1]), "_")[[1]][1], "/")[[1]][3]
  sdata <- data[data$fire_name == f,]
  if ("dnbr" %in% names(sdata)) sdata[, dnbr := NULL]
  sdata[, orig_row := .I]
  sdata <- merge(sdata, sev[,.(x,y,dnbr)], by = c("x", "y"), all.x = TRUE)
  setorder(sdata, orig_row)  # Restore original order
  sdata[, orig_row := NULL]  # Remove the temporary column
  data[fire_name == f, dnbr := sdata$dnbr]
}

data[, hs_dnbr := NA_integer_]
data[data$dnbr > 367, hs_dnbr := 1]
data[data$dnbr <= 367, hs_dnbr := 0]

rdnbr_files <- list.files(path = "data/rdnbr", pattern = "\\.csv$", full.names = TRUE)
data[, rdnbr := NA_real_]

for (i in seq_along(rdnbr_files)) {
  sev <- fread(rdnbr_files[i])
  f <- strsplit(strsplit(tolower(strsplit(rdnbr_files[i], "_processed")[[1]][1]), "_")[[1]][1], "/")[[1]][3]
  sdata <- data[data$fire_name == f,]
  if ("rdnbr" %in% names(sdata)) sdata[, rdnbr := NULL]
  sdata[, orig_row := .I]
  sdata <- merge(sdata, sev[,.(x,y,rdnbr)], by = c("x", "y"), all.x = TRUE)
  setorder(sdata, orig_row)  # Restore original order
  sdata[, orig_row := NULL]  # Remove the temporary column
  data[fire_name == f, rdnbr := sdata$rdnbr]
}

data[, hs_rdnbr := NA_integer_]
data[data$rdnbr > 641, hs_rdnbr := 1]
data[data$rdnbr <= 641, hs_rdnbr := 0]

fwrite(data, "data/complete_data.csv", row.names = FALSE)
