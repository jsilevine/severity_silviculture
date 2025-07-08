
# Load necessary libraries
library(ggplot2)
library(reshape2)
library(data.table)

##---------------------------------------------------------------
## Neighborhood scale
##---------------------------------------------------------------

# Use built-in dataset
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

data <- data[,c("hdw_scaled", "avg_fuel_moisture_scaled", "prev_sev_scaled",
                     "mean_dens_30_scaled",
                     "clust_30_scaled", "mean_ht_30_scaled",
                     "mean_area_30_scaled",
                     "em_scaled",
                     "cwd_scaled", "slope_scaled", "tpi_scaled", "heat_load_scaled")]

# Compute correlation matrix
cor_matrix <- round(cor(data), 2)

# Melt the correlation matrix into long format
cor_long <- reshape2::melt(cor_matrix)
for (i in 1:nrow(cor_long)) {
  if(cor_long[i,"Var1"] == cor_long[i, "Var2"]) {
    cor_long[i,"value"] <- NA
  }
}

write.csv(cor_long, "data/correlation_matrices/correlation_matrix_neighborhood.csv", row.names = FALSE)

# Plot using ggplot2
ggplot(cor_long, aes(x = Var1, y = Var2)) +
  geom_point(aes(size = abs(value), color = value)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Matrix", x = "", y = "", size = "Strength")
ggsave("plots/figureS11.pdf")

cor_long[abs(cor_long$value) > 0.5 & !is.na(cor_long$value), ]

##---------------------------------------------------------------
## Stand scale
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

data <- data[,c("hdw_scaled", "avg_fuel_moisture_scaled", "prev_sev_scaled",
                     "mean_dens_180_scaled",
                     "clust_180_scaled", "mean_ht_180_scaled",
                     "mean_area_180_scaled",
                     "em_scaled",
                     "cwd_scaled", "slope_scaled", "tpi_scaled", "heat_load_scaled")]

# Compute correlation matrix
cor_matrix <- round(cor(data), 2)

# Melt the correlation matrix into long format
cor_long <- reshape2::melt(cor_matrix)
for (i in 1:nrow(cor_long)) {
  if(cor_long[i,"Var1"] == cor_long[i, "Var2"]) {
    cor_long[i,"value"] <- NA
  }
}

write.csv(cor_long, "data/correlation_matrices/correlation_matrix_stand.csv", row.names = FALSE)

# Plot using ggplot2
ggplot(cor_long, aes(x = Var1, y = Var2)) +
  geom_point(aes(size = abs(value), color = value)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Matrix", x = "", y = "", size = "Strength")
ggsave("plots/figureS12.pdf")
