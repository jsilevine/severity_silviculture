library(data.table)
library(ggplot2)
library(ggfortify)
library(biglm)
library(cowplot)
library(gstat)
library(sf)

data <- fread("data/complete_data.csv")

ggplot(data = data[data$fire_name == "north_complex",], aes(x = x, y = y, fill = CBI_bc)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_gradientn(
    colours = c("#7fbc41", "white", "#c51b7d"),
    values = c(0, 0.7, 1)) +
  theme(panel.background = element_rect(fill = "#292e3d"),
        plot.background = element_rect(fill = "#292e3d"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("../figtest.png")



data$fire_name <- as.factor(data$fire_name)
data$own_type <- as.factor(data$own_type)
data <- data[complete.cases(data),]

data$clust_s <- scale(data[,clust])
data$mean_dens_s <- scale(data[,mean_dens])
data$mean_area_s <- scale(data[,mean_area])
data$median_area_s <- scale(data[,median_area])
data$em_s <- scale(data[,em])
data$elevation_s <- scale(data[,elevation])
data$slope_s <- scale(data[,slope])
data$tpi_s <- scale(data[,tpi])
data$heat_load_s <- scale(data[,heat_load])
data$prev_sev_s <- scale(data[,prev_sev])
data$avg_fuel_moisture_s <- scale(data[,avg_fuel_moisture])

calc_es <- function(temp){
  es <- 6.11 * exp((2.5e6 / 461) * (1 / 273 - 1 / (273 + temp)))
  return(es)
}

calc_vpd <- function(rh, temp){
  es <- calc_es(temp)
  vpd <- ((100 - rh) / 100) * es
  return(vpd)
}

data$vpd <- unlist(mapply(calc_vpd, data$avg_fuel_moisture, data$avg_air_temp))
data$hdw <- data$vpd * data$avg_wind_speed
data$hdw_s <- scale(data[,hdw])

pc <- prcomp(data[,.(clust_s, mean_dens_s, median_area_s, em_s)], retx = TRUE)
summary(pc)
autoplot(pc, data = data[own_type == "Private Industrial"], color = 'own_type', loadings = TRUE, loadings.label = TRUE, alpha = 0.25)
autoplot(pc, data = data[own_type == "Federal"], color = 'own_type', loadings = TRUE, loadings.label = TRUE, alpha = 0.25)

pcx <- pc$x
pcx <- as.data.frame(pcx)
colnames(pcx)

data <- cbind(data, pcx)
data$PC1 <- scale(data[,PC1])
data$PC2 <- scale(data[,PC2])
data$PC3 <- scale(data[,PC3])
data$PC4 <- scale(data[,PC4])

# make density plots

pc1 <- seq(min(data$PC1), max(data$PC1), length = 100)
pc2 <- seq(min(data$PC2), max(data$PC2), length = 100)
pc3 <- seq(min(data$PC3), max(data$PC3), length = 100)

pc1.gs <- pc1[2] - pc1[1]
pc2.gs <- pc2[2] - pc2[1]
pc3.gs <- pc3[2] - pc3[1]

pc_grid <- expand.grid(pc1, pc3)
own <- "Private Industrial"
data_sub <- data[own_type == own]
for (i in 1:nrow(pc_grid)) {
  pc_grid[i,3] <- nrow(data_sub[PC1 > pc_grid[i,1] - pc1.gs & PC1 < pc_grid[i,1] + pc1.gs &
           PC3 > pc_grid[i,2] - pc3.gs & PC3 < pc_grid[i,2] + pc3.gs])
}



own <- "Federal"
data_sub <- data[own_type == own]
for (i in 1:nrow(pc_grid)) {
  pc_grid[i,4] <- nrow(data_sub[PC1 > pc_grid[i,1] - pc1.gs & PC1 < pc_grid[i,1] + pc1.gs &
           PC3 > pc_grid[i,2] - pc3.gs & PC3 < pc_grid[i,2] + pc3.gs])
}

pc_grid[,3] <- pc_grid[,3] / sum(pc_grid[,3])
pc_grid[,4] <- pc_grid[,4] / sum(pc_grid[,4])
pc_grid[,5] <- pc_grid[,3] - pc_grid[,4]

ggplot(pc_grid[abs(pc_grid[,5])>0,], aes(x = Var1, y = Var2, fill = V5)) +
  geom_tile() +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  ggtitle("Private")

pc_grid


p1 <- ggplot(data[own_type == "Private Industrial"], aes(x = PC1, y = PC2)) +
  geom_hex(bins = 70) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  ggtitle("Private")

p2 <- ggplot(data[own_type == "Federal"], aes(x = PC1, y = PC2)) +
  geom_hex(bins = 70) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  ggtitle("Federal")

p3 <- ggplot(data[own_type == "Other"], aes(x = PC1, y = PC2)) +
  geom_hex(bins = 70) +
  #stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  ggtitle("Other")

plot_grid(p1, p2, p3)


nrow(data[hs == 0]) / nrow(data)

##---------------------------------------------------------------
## Fit naive models (no S.E. corrections)
##---------------------------------------------------------------

full_model_gamma <- glm(CBI_bc ~ fire_name + PC1 + PC2 + PC3 + PC4 +
                       elevation_s + slope_s + tpi_s + heat_load_s + prev_sev_s +
                       hdw_s * avg_fuel_moisture_s + PC1:hdw_s + PC2:hdw_s +
                       PC3:hdw_s + PC4:hdw_s, data = data[data$CBI_bc > 0, ], family = Gamma(log),
                     maxit = 1000)
summary(full_model_gamma)

full_model_binomial <- bigglm(hs ~ fire_name + PC1 + PC2 + PC3 + PC4 +
                       elevation_s + slope_s + tpi_s + heat_load_s + prev_sev_s +
                       hdw_s + avg_fuel_moisture_s + PC1:hdw_s + PC2:hdw_s +
                       PC3:hdw_s + PC4:hdw_s, data = data, family = binomial(link = "logit"),
                     maxit = 1000)

summary(full_model_binomial)

summary(bigglm(CBI_bc ~ fire_name + own_type +
               elevation_s + slope_s + tpi_s + heat_load_s + prev_sev_s +
               hdw_s + avg_fuel_moisture_s + own_type:hdw_s, data = data, family = Gamma(log),
               maxit = 1000))

##---------------------------------------------------------------
## Determine scale of residual autocorrelation
##---------------------------------------------------------------


## calculate deviance residuals

## logodds <- predict(full_model_binomial, newdata = data)
## pr <- 1 / (1 + exp(-logodds))

## resid_naive <- numeric(nrow(data))
## resid_naive[data$hs == 1] <- sqrt(-2 * log(pr[data$hs == 1]))
## resid_naive[data$hs == 0] <- -1 * sqrt(-2 * log(1 - pr[data$hs == 0]))

resid_naive <- residuals(full_model_gamma)
df.resid <- data.frame(z = resid_naive, x = data[data$CBI_bc > 0, "x"], y = data[data$CBI_bc > 0, "y"])
df.sub <- df.resid[sample(1:nrow(df.resid), round(nrow(df.resid) / 10)),]
sp_df <- st_as_sf(df.sub, coords = c("x", "y"), crs = 4326)
sp_df <- st_transform(sp_df, crs = 32610)
df.sub <- data.frame(sp_df)
df.sub[, c("x", "y")] <- st_coordinates(sp_df)

v1 <- variogram(z~1, data = df.sub, locations = ~x+y, cutoff = 1500, width = 1500 / 30)
f1 <- fit.variogram(v1, vgm("Sph"))
plot(v1, model = f1)
max.dist <- f1$range[2]
round(max.dist, 3)
saveRDS(v1, "data/model_objects/variogram.rds")

write.csv(round((1000 / 30) * (data[3000, "x"] - data[2999, "x"]), 3), "data/autocor_scale.csv")


##---------------------------------------------------------------
## Implement block bootstrapping
##---------------------------------------------------------------

## first convert data to UTM zone 10N and save file
sp_data <- st_as_sf(data, coords = c("x", "y"), crs = 4326)
sp_data <- st_transform(sp_data, crs = 32610)
data_transform <- data.frame(sp_data)
data_transform <- data_transform[, 1:47]
data_transform[, c("x", "y")] <- st_coordinates(sp_data)
data_transform <- data.table(data_transform)
fwrite(data_transform, "data/complete_data_utm10.csv")

data_transform[4000, "x"] - data_transform[3999, "x"]

data$mean_area

data_for_QLG <- data[, c("x", "y", "fire_name", "clust", "mean_dens", "mean_area", "sd_area", "percent_open", "em")]
colnames(data_for_QLG) <- c("x", "y", "fire_name", "clustering_index", "mean_stem_density_per_ha", "mean_gap_area_m",
                            "sd_gap_area_m", "percent_open", "ladder_fuels_metric")
write.csv(data_for_QLG, "data/data_for_QLG.csv")
