
library(terra)

source("code/utility/utility_functions.r")

cwd <- rast("data/cwd/cwd1981_2010_ave_HST_1678990308.tif")
template <- readRDS("data/templates/isforest_template.rds")

cwd <- snap_to_template(cwd, rast(template))

writeRaster(cwd, "data/cwd/cwd_raster.tif", overwrite = TRUE)
saveRDS(cwd, "data/cwd/cwd_raster.rds")

cwd_df <- as.data.frame(cwd, na.rm = TRUE, xy = TRUE)
colnames(cwd_df)[3] <- "cwd"

write.csv(cwd_df, "data/cwd/cwd.csv")
