setwd("~/Documents/Science/fire/severity_and_silviculture/")

library(data.table)
library(speedglm)
library(sf)
library(gstat)
library(parallel)
library(ggplot2)

##---------------------------------------------------------------
## Neighborhood Scale
##---------------------------------------------------------------

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

set.seed(1211) #bday
s <- sample(1:nrow(data), round(0.25 * nrow(data)))

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

sub_lev <- seq(0.02, 0.62, by = 0.04)

output <- matrix(0, length(naive_model$coefficients), length(sub_lev))
row.names(output) <- names(naive_model$coefficients)

for (i in 1:length(sub_lev)) {

  set.seed(1211) #bday
  s <- sample(1:nrow(data), round(sub_lev[i] * nrow(data)))

  nm <- glm(hs ~ hdw_scaled + avg_fuel_moisture_scaled + prev_sev_scaled +
              mean_dens_30_scaled + mean_dens_30_scaled:hdw_scaled +
              clust_30_scaled +
              mean_ht_30_scaled + mean_ht_30_scaled:hdw_scaled +
              mean_area_30_scaled + mean_area_30_scaled:hdw_scaled +
              em_scaled + em_scaled:hdw_scaled +
              cwd_scaled + slope_scaled + tpi_scaled + heat_load_scaled +
              fire_name,
            data = data[s,], family = binomial(), model = FALSE, y = FALSE)

  output[,i] <- nm$coefficients

}

df_long <- data.frame(
  coef = rep(rownames(output), each = ncol(output)),
  sl = rep(sub_lev, times = nrow(output)),
  value = as.vector(t(output))
)

df_long

ggplot(df_long, aes(x = sl, y = value, group = coef)) +
  geom_line() +
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "red") +
  facet_wrap(~ coef) +  # separate panels for each row
  theme_bw() +
  labs(x = "Subsample Level (prop.)", y = "Estimate", title = "Sensitivity to Subsample Level (Neighborhood)")

ggsave("plots/figureS7.pdf")



## sensitivity to CBI threshold

cbi_lev <- seq(2.0, 2.5, length.out = 20)

output <- matrix(0, length(naive_model$coefficients), length(cbi_lev))
row.names(output) <- names(naive_model$coefficients)

for (i in 1:length(cbi_lev)) {

  nd <- data
  nd$hs <- 0
  nd[nd$cbi > cbi_lev[i], "hs"] <- 1

  set.seed(1211) #bday
  s <- sample(1:nrow(nd), round(0.25 * nrow(nd)))

  nm <- glm(hs ~ hdw_scaled + avg_fuel_moisture_scaled + prev_sev_scaled +
              mean_dens_30_scaled + mean_dens_30_scaled:hdw_scaled +
              clust_30_scaled +
              mean_ht_30_scaled + mean_ht_30_scaled:hdw_scaled +
              mean_area_30_scaled + mean_area_30_scaled:hdw_scaled +
              em_scaled + em_scaled:hdw_scaled +
              cwd_scaled + slope_scaled + tpi_scaled + heat_load_scaled +
              fire_name,
            data = nd[s,], family = binomial(), model = FALSE, y = FALSE)

  output[,i] <- nm$coefficients

}

df_long <- data.frame(
  coef = rep(rownames(output), each = ncol(output)),
  cbi = rep(cbi_lev, times = nrow(output)),
  value = as.vector(t(output))
)

df_long

ggplot(df_long, aes(x = cbi, y = value, group = coef)) +
  geom_line() +
  geom_vline(xintercept = 2.25, linetype = "dashed", color = "red") +
  facet_wrap(~ coef) +  # separate panels for each row
  theme_bw() +
  labs(x = "CBI Threshold", y = "Estimate", title = "Sensitivity to CBI Threshold (Neighborhood)")

ggsave("plots/figureS8.pdf")

##---------------------------------------------------------------
## stand scale
##---------------------------------------------------------------

data <- fread("data/complete_data.csv")
data <- data[complete.cases(data),]

data$mean_dens_180_scaled <- scale(data$mean_dens_180)
data$clust_180_scaled <- scale(data$clust_180)
data$mean_ht_180_scaled <- scale(data$mean_ht_180)
data$mean_area_180_scaled <- scale(data$mean_area_180)
data$em_scaled <- scale(data$em)
data$hdw_scaled <- scale(data$hdw)
data$avg_fuel_moisture_scaled <- scale(data$avg_fuel_moisture)
data$cwd_scaled <- scale(data$cwd)
data$slope_scaled <- scale(data$slope)
data$tpi_scaled <- scale(data$tpi)
data$heat_load_scaled <- scale(data$heat_load)
data$prev_sev_scaled <- scale(data$prev_sev)

set.seed(1211) #bday
s <- sample(1:nrow(data), round(0.25 * nrow(data)))

naive_model <- glm(hs ~ hdw_scaled + avg_fuel_moisture_scaled + prev_sev_scaled +
                     mean_dens_180_scaled + mean_dens_180_scaled:hdw_scaled +
                     clust_180_scaled +
                     mean_ht_180_scaled + mean_ht_180_scaled:hdw_scaled +
                     mean_area_180_scaled + mean_area_180_scaled:hdw_scaled +
                     em_scaled + em_scaled:hdw_scaled +
                     cwd_scaled + slope_scaled + tpi_scaled + heat_load_scaled +
                     fire_name,
                     data = data[s,], family = binomial(), model = FALSE, y = FALSE)
summary(naive_model)

sub_lev <- seq(0.02, 0.62, by = 0.04)

output <- matrix(0, length(naive_model$coefficients), length(sub_lev))
row.names(output) <- names(naive_model$coefficients)

for (i in 1:length(sub_lev)) {

  set.seed(1211) #bday
  s <- sample(1:nrow(data), round(sub_lev[i] * nrow(data)))

  nm <- glm(hs ~ hdw_scaled + avg_fuel_moisture_scaled + prev_sev_scaled +
              mean_dens_180_scaled + mean_dens_180_scaled:hdw_scaled +
              clust_180_scaled +
              mean_ht_180_scaled + mean_ht_180_scaled:hdw_scaled +
              mean_area_180_scaled + mean_area_180_scaled:hdw_scaled +
              em_scaled + em_scaled:hdw_scaled +
              cwd_scaled + slope_scaled + tpi_scaled + heat_load_scaled +
              fire_name,
            data = data[s,], family = binomial(), model = FALSE, y = FALSE)

  output[,i] <- nm$coefficients

}

df_long <- data.frame(
  coef = rep(rownames(output), each = ncol(output)),
  sl = rep(sub_lev, times = nrow(output)),
  value = as.vector(t(output))
)

df_long

ggplot(df_long, aes(x = sl, y = value, group = coef)) +
  geom_line() +
  geom_vline(xintercept = 0.25, linetype = "dashed", color = "red") +
  facet_wrap(~ coef) +  # separate panels for each row
  theme_bw() +
  labs(x = "Subsample Level (prop.)", y = "Estimate", title = "Sensitivity to Subsample Level (Stand)")

ggsave("plots/figureS9.pdf")



## sensitivity to CBI threshold

cbi_lev <- seq(2.0, 2.5, length.out = 20)

output <- matrix(0, length(naive_model$coefficients), length(cbi_lev))
row.names(output) <- names(naive_model$coefficients)

for (i in 1:length(cbi_lev)) {

  nd <- data
  nd$hs <- 0
  nd[nd$cbi > cbi_lev[i], "hs"] <- 1

  set.seed(1211) #bday
  s <- sample(1:nrow(nd), round(0.25 * nrow(nd)))

  nm <- glm(hs ~ hdw_scaled + avg_fuel_moisture_scaled + prev_sev_scaled +
              mean_dens_180_scaled + mean_dens_180_scaled:hdw_scaled +
              clust_180_scaled +
              mean_ht_180_scaled + mean_ht_180_scaled:hdw_scaled +
              mean_area_180_scaled + mean_area_180_scaled:hdw_scaled +
              em_scaled + em_scaled:hdw_scaled +
              cwd_scaled + slope_scaled + tpi_scaled + heat_load_scaled +
              fire_name,
            data = nd[s,], family = binomial(), model = FALSE, y = FALSE)

  output[,i] <- nm$coefficients

}

df_long <- data.frame(
  coef = rep(rownames(output), each = ncol(output)),
  cbi = rep(cbi_lev, times = nrow(output)),
  value = as.vector(t(output))
)

df_long

ggplot(df_long, aes(x = cbi, y = value, group = coef)) +
  geom_line() +
  geom_vline(xintercept = 2.25, linetype = "dashed", color = "red") +
  facet_wrap(~ coef) +  # separate panels for each row
  theme_bw() +
  labs(x = "CBI Threshold", y = "Estimate", title = "Sensitivity to CBI Threshold (Stand)")

ggsave("plots/figureS10.pdf")



nrow(data[data$own_type == "Federal",]) / nrow(data)
