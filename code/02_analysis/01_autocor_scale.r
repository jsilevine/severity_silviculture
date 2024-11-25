
library(data.table)
library(speedglm)
library(sf)
library(gstat)
library(parallel)

data <- fread("data/complete_data.csv")
data <- data[complete.cases(data),]

data$mean_dens_30_scaled <- scale(data$mean_dens_30)
data$clust_30_scaled <- scale(data$clust_30)
data$mean_ht_30_scaled <- scale(data$mean_ht_30)
data$mean_area_30_scaled <- scale(data$mean_area_30)
data$em_scaled <- scale(data$em)
data$hdw_scaled <- scale(data$hdw)
data$avg_fuel_moisture_scaled <- scale(data$avg_fuel_moisture)
data$cwd_scaled <- scale(data$cwd)
data$slope_scaled <- scale(data$slope)
data$tpi_scaled <- scale(data$tpi)
data$heat_load_scaled <- scale(data$heat_load)
data$prev_sev_scaled <- scale(data$prev_sev)

set.seed(1211)
s <- sample(1:nrow(data), round(0.3 * nrow(data)))

naive_model <- glm(hs ~ hdw_scaled + avg_fuel_moisture_scaled + prev_sev_scaled +
                     mean_dens_30_scaled + mean_dens_30_scaled:hdw_scaled +
                     clust_30_scaled +
                     mean_ht_30_scaled + mean_ht_30_scaled:hdw_scaled +
                     mean_area_30_scaled + mean_area_30_scaled:hdw_scaled +
                     em_scaled + em_scaled:hdw_scaled +
                     cwd_scaled + slope_scaled + tpi_scaled + heat_load_scaled +
                     fire_name,
                  data = data[s,], family = binomial(), model = FALSE, y = FALSE)
summary(naive_model)
saveRDS(naive_model, "data/model_objects/naive_model_structure.rds")


## create variogram to estimate maximum range of autocorrelation
resid.rm.glm <- residuals(naive_model)
df.resid <- data.frame(z = resid.rm.glm, x = data[s,"x"], y = data[s,"y"])
df.resid <- df.resid[complete.cases(df.resid),]
set.seed(3)
df.resid <- df.resid[sample(1:nrow(df.resid), 200000),]
v1 <- variogram(z~1, data = df.resid, locations = ~x+y, cutoff = 3000)
f1 <- fit.variogram(v1, vgm("Sph"))
max.dist <- f1$range[2]
bb <- ceiling(max.dist) ##round up to nearest whole number
saveRDS(bb, "data/autocor_scale_m.rds")
