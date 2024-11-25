library(speedglm)
library(data.table)
library(ggplot2)
library(cowplot)
library(paletteer)
library(metBrewer)

data <- fread("data/complete_data.csv")
data <- data[complete.cases(data),]
bootstrap_results <- read.csv("data/bootstraps/bootstrap_coefs_30.csv")
naive_model <- readRDS("data/model_objects/naive_model_structure.rds")

## scale continuous predictors
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

head(bootstrap_results)
colnames(bootstrap_results)[1] <- "intercept"

results_summary <- data.frame(variable = colnames(bootstrap_results),
                              mean = numeric(ncol(bootstrap_results)),
                              ci_lower = numeric(ncol(bootstrap_results)),
                              ci_upper = numeric(ncol(bootstrap_results)))

for (i in 1:nrow(results_summary)) {

  results_summary[i,"mean"] <- mean(bootstrap_results[,results_summary[i,"variable"]])
  results_summary[i,"ci_lower"] <- quantile(bootstrap_results[,results_summary[i,"variable"]], 0.025)
  results_summary[i,"ci_upper"] <- quantile(bootstrap_results[,results_summary[i,"variable"]], 0.975)

}

results_summary

write.csv(results_summary, "data/model_objects/30m_results_summary.csv")

gen_pdata <- function(var1, var2 = NA, var2_levels = NA, npts = 100) {

  vars <- c("mean_dens_30_scaled", "hdw_scaled", "clust_30_scaled", "mean_ht_30_scaled",
            "mean_area_30_scaled", "em_scaled", "avg_fuel_moisture_scaled", "cwd_scaled", "tpi_scaled",
            "heat_load_scaled", "prev_sev_scaled", "slope_scaled")


  if (!is.na(var2)) {

    pdata <- data.frame(x = numeric(npts * length(var2_levels)))
    indx <- which(colnames(data) == var1)

    pdata[,var1] <- rep(seq(min(data[[indx]]), max(data[[indx]]), length.out = npts), times = length(var2_levels))
    pdata[,var2] <- rep(var2_levels, each = npts)

    ocols <- vars[!(vars %in% c(var1, var2))]
    ocols

    for (i in 1:length(ocols)) {
      pdata[,ocols[i]] <- 0
    }

  } else {

    pdata <- data.frame(x = numeric(npts))
    indx <- which(colnames(data) == var1)

    pdata[,var1] <- seq(min(data[[indx]]), max(data[[indx]]), length.out = npts)

    ocols <- vars[!(vars %in% c(var1))]
    ocols

    for (i in 1:length(ocols)) {
      pdata[,ocols[i]] <- 0
    }

  }

    pdata <- pdata[,2:ncol(pdata)]
    pdata$fire_name <- "dixie"

    pdata[,"p_mean"] <- predict(naive_model, newdata = pdata)

    head(pdata)

    for (i in 1:nrow(pdata)) {

      pvec <- numeric(nrow(bootstrap_results))
      for (j in 1:nrow(bootstrap_results)) {

        pvec[j] <-
          bootstrap_results[j,"intercept"] +
          bootstrap_results[j,"mean_dens_30_scaled"] * pdata[i, "mean_dens_30_scaled"] +
          bootstrap_results[j,"hdw_scaled"] * pdata[i, "hdw_scaled"] +
          bootstrap_results[j,"clust_30_scaled"] * pdata[i, "clust_30_scaled"] +
          bootstrap_results[j,"mean_ht_30_scaled"] * pdata[i, "mean_ht_30_scaled"] +
          bootstrap_results[j,"mean_area_30_scaled"] * pdata[i, "mean_area_30_scaled"] +
          bootstrap_results[j,"em_scaled"] * pdata[i, "em_scaled"] +
          bootstrap_results[j,"avg_fuel_moisture_scaled"] * pdata[i, "avg_fuel_moisture_scaled"] +
          bootstrap_results[j,"cwd_scaled"] * pdata[i, "cwd_scaled"] +
          bootstrap_results[j,"slope_scaled"] * pdata[i, "slope_scaled"] +
          bootstrap_results[j,"tpi_scaled"] * pdata[i, "tpi_scaled"] +
          bootstrap_results[j,"heat_load_scaled"] * pdata[i, "heat_load_scaled"] +
          bootstrap_results[j,"prev_sev_scaled"] * pdata[i, "prev_sev_scaled"] +
          bootstrap_results[j,"hdw_scaled.mean_dens_30_scaled"] * pdata[i, "mean_dens_30_scaled"] * pdata[i, "hdw_scaled"] +
          bootstrap_results[j,"hdw_scaled.mean_ht_30_scaled"] * pdata[i, "mean_ht_30_scaled"] * pdata[i, "hdw_scaled"] +
          bootstrap_results[j,"hdw_scaled.mean_area_30_scaled"] * pdata[i, "mean_area_30_scaled"] * pdata[i, "hdw_scaled"] +
          bootstrap_results[j,"hdw_scaled.em_scaled"] * pdata[i, "em_scaled"] * pdata[i, "hdw_scaled"]

      }

      pdata[i,"p_lower"] <- quantile(pvec, 0.025)
      pdata[i,"p_higher"] <- quantile(pvec, 0.975)
      #pdata[i, "p_mean"] <- mean(pvec)

    }


    pdata_descaled <- pdata
    for (i in 1:(ncol(pdata_descaled)-4)) {

      cname <- gsub(pattern = "_scaled", x = colnames(pdata_descaled)[i], replacement = "")
      colnames(pdata_descaled)[i] <- cname
      indx <- which(colnames(data) == cname)
      pdata_descaled[,i] <- pdata_descaled[,i] * sd(data[[indx]]) + mean(data[[indx]])

    }
    head(pdata_descaled)

    ## transform predictions from log odds to probability
    pdata_descaled$p_mean <- exp(pdata_descaled$p_mean) / (1 + exp(pdata_descaled$p_mean))
    pdata_descaled$p_lower <- exp(pdata_descaled$p_lower) / (1 + exp(pdata_descaled$p_lower))
    pdata_descaled$p_higher <- exp(pdata_descaled$p_higher) / (1 + exp(pdata_descaled$p_higher))

  return(pdata_descaled)

}

pd_clust <- gen_pdata("clust_30_scaled")
pd_clust$prev_sev <- as.factor(pd_clust$prev_sev)

clust_pred <- ggplot(data = pd_clust, aes(x = clust_30, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_higher), fill = "gray", alpha = 0.9) +
  geom_line(color = "black", linewidth = 2) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  xlab("Spatial homogeneity") +
  ylab("High severity probability") +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
clust_pred

## Hot dry windy
pd_hdw <- gen_pdata("hdw_scaled")

hdw_pred <- ggplot(data = pd_hdw, aes(x = hdw, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_higher), fill = "gray", alpha = 0.9) +
  geom_line(color = "black", size = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  xlab("Hot dry windy (HDW) index") +
  ylab("High severity probability") +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
hdw_pred

## mean gap area
pd_area <- gen_pdata("mean_area_30_scaled", npts = 1000)

area_pred <- ggplot(data = pd_area, aes(x = mean_area_30, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_higher), fill = "gray", alpha = 0.9) +
  geom_line(color = "black", linewidth = 2) +
  scale_x_continuous(expand = c(0,0), limits = c(NA, NA)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  xlab("Mean gap area") +
  ylab("High severity probability") +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
area_pred

## fuel moisture
pd_fuel <- gen_pdata("avg_fuel_moisture_scaled")

fuel_pred <- ggplot(data = pd_fuel, aes(x = avg_fuel_moisture, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_higher), fill = "gray", alpha = 0.9) +
  geom_line(color = "black", size = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  xlab("Mean fuel moisture") +
  ylab("High severity probability") +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
fuel_pred

## cwd
pd_cwd <- gen_pdata("cwd_scaled")

cwd_pred <- ggplot(data = pd_cwd, aes(x = cwd, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_higher), fill = "gray", alpha = 0.9) +
  geom_line(color = "black", size = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  guides(x.sec = "axis", y.sec = "axis") +
  xlab("cwd") +
  ylab("High severity probability") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
cwd_pred

## slope
pd_slope <- gen_pdata("slope_scaled")

slope_pred <- ggplot(data = pd_slope, aes(x = slope, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_higher), fill = "gray", alpha = 0.9) +
  geom_line(color = "black", size = 1, linetype = "dashed") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  xlab("Slope") +
  ylab("High severity probability") +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
slope_pred

## topographic position index
pd_tpi <- gen_pdata("tpi_scaled")

tpi_pred <- ggplot(data = pd_tpi, aes(x = tpi, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_higher), fill = "gray", alpha = 0.9) +
  geom_line(color = "black", size = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  xlab("Topographic position index") +
  ylab("High severity probability") +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
tpi_pred

## heat load
pd_heatload <- gen_pdata("heat_load_scaled")

heatload_pred <- ggplot(data = pd_heatload, aes(x = heat_load, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_higher), fill = "gray", alpha = 0.9) +
  geom_line(color = "black", size = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  xlab("Heat load") +
  ylab("High severity probability") +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
heatload_pred

## psev
pd_psev <- gen_pdata("prev_sev_scaled")

psev_pred <- ggplot(data = pd_psev, aes(x = prev_sev, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_higher), fill = "gray", alpha = 0.9) +
  geom_line(color = "black", size = 1) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +
  xlab("CBI (previous time interval)") +
  ylab("High severity probability") +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
psev_pred


##---------------------------------------------------------------
## Interactions with HDW
##---------------------------------------------------------------

x <- c(10, 75, 150, 250)
x <- (x - mean(data$hdw)) / sd(data$hdw)

pd_dens_hdw <- gen_pdata("mean_dens_30_scaled", "hdw_scaled", x, npts = 100)
pd_dens_hdw$hdw <- as.factor(pd_dens_hdw$hdw)

dens_hdw_pred <- ggplot(data = pd_dens_hdw, aes(x = mean_dens_30, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_higher, fill = hdw), alpha = 0.7) +
  geom_line(aes(color = hdw), size = 2) +
  scale_x_continuous(expand = c(0,0), limits = c(40, NA)) +
  scale_y_continuous(expand = c(0,0), limits = c(0.0, 1)) +
  scale_color_manual(values = c("#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d")) +
  scale_fill_manual(values = c("#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d")) +
  xlab("Mean stem density") +
  ylab("High severity probability") +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
dens_hdw_pred

pd_ht_hdw <- gen_pdata("mean_ht_30_scaled", "hdw_scaled", x, npts = 100)
pd_ht_hdw$hdw <- as.factor(pd_ht_hdw$hdw)

ht_hdw_pred <- ggplot(data = pd_ht_hdw, aes(x = mean_ht_30, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_higher, fill = hdw), alpha = 0.7) +
  geom_line(aes(color = hdw), size = 2) +
  scale_x_continuous(expand = c(0,0), limits = c(8.5, NA)) +
  scale_y_continuous(expand = c(0,0), limits = c(0.0, 1)) +
  scale_color_manual(values = c("#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d")) +
  scale_fill_manual(values = c("#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d")) +
  xlab("Mean stem height (m)") +
  ylab("High severity probability") +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
ht_hdw_pred

pd_area_hdw <- gen_pdata("mean_area_30_scaled", "hdw_scaled", x, npts = 100)
pd_area_hdw$hdw <- as.factor(pd_area_hdw$hdw)

area_hdw_pred <- ggplot(data = pd_area_hdw, aes(x = mean_area_30, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_higher, fill = hdw), alpha = 0.7) +
  geom_line(aes(color = hdw), size = 2) +
  scale_x_continuous(expand = c(0,0), limits = c(NA, NA)) +
  scale_y_continuous(expand = c(0,0), limits = c(0.0, 1)) +
  scale_color_manual(values = c("#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d")) +
  scale_fill_manual(values = c("#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d")) +
  xlab("Mean gap area (m^2)") +
  ylab("High severity probability") +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
area_hdw_pred

pd_em_hdw <- gen_pdata("em_scaled", "hdw_scaled", x, npts = 100)
pd_em_hdw$hdw <- as.factor(pd_em_hdw$hdw)

em_hdw_pred <- ggplot(data = pd_em_hdw, aes(x = em, y = p_mean)) +
  geom_ribbon(aes(ymin = p_lower, ymax = p_higher, fill = hdw), alpha = 0.7) +
  geom_line(aes(color = hdw), size = 2) +
  scale_x_continuous(expand = c(0,0), limits = c(NA, NA)) +
  scale_y_continuous(expand = c(0,0), limits = c(0.0, 1)) +
  scale_color_manual(values = c("#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d")) +
  scale_fill_manual(values = c("#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d")) +
  xlab("Ladder fuels index") +
  ylab("High severity probability") +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
em_hdw_pred

vars <- c("mean_dens_30_scaled", "hdw_scaled", "clust_30_scaled", "mean_ht_30_scaled",
          "mean_area_30_scaled", "em_scaled", "avg_fuel_moisture_scaled", "cwd_scaled", "tpi_scaled",
          "heat_load_scaled", "prev_sev_scaled", "slope_scaled")

ndata <- data.frame(dixie = numeric(1),
                    northcomplex = numeric(1),
                    sheep = numeric(1),
                    sugar = numeric(1),
                    walker = numeric(1))

for (f in unique(data$fire_name)) {
  nndata <- ndata
  nndata[,"fire_name"] <- f
  nndata[,f] <- 1
  if (f == unique(data$fire_name)[1]) {
    pdata <- nndata
  } else {
    pdata <- rbind(pdata, nndata)
  }

}

for (i in 1:nrow(pdata)) {

  pvec <- numeric(nrow(bootstrap_results))
  for (j in 1:nrow(bootstrap_results)) {

    pvec[j] <-
      bootstrap_results[j,"intercept"] +
      bootstrap_results[j,"fire_namenorthcomplex"] * pdata[i, "northcomplex"] +
      bootstrap_results[j,"fire_namesheep"] * pdata[i, "sheep"] +
      bootstrap_results[j,"fire_namesugar"] * pdata[i, "sugar"] +
      bootstrap_results[j,"fire_namewalker"] * pdata[i, "walker"]

  }
  pdata[i, "p_mean"] <- mean(pvec)
  pdata[i,"p_lower"] <- quantile(pvec, 0.025)
  pdata[i,"p_higher"] <- quantile(pvec, 0.975)

}

## transform predictions from log odds to probability
pdata$p_mean <- exp(pdata$p_mean) / (1 + exp(pdata$p_mean))
pdata$p_lower <- exp(pdata$p_lower) / (1 + exp(pdata$p_lower))
pdata$p_higher <- exp(pdata$p_higher) / (1 + exp(pdata$p_higher))

colors = met.brewer(name = "Archambault", n = 5)
met.brewer(name = "Egypt", n = 3)
ggsave("palette.pdf")

## fires
fireid_pred <- ggplot(data = pdata, aes(x = fire_name, y = p_mean)) +
  geom_errorbar(aes(ymin = p_lower, ymax = p_higher, color = fire_name), linewidth = 1, width = 0.5) +
  geom_point(aes(color = fire_name), size = 2) +
  scale_y_continuous(expand = c(0,0), limits = c(0.0, 1.0)) +
  scale_color_manual(name = "Fire Name", values = colors) +
  xlab("Fire Name") +
  ylab("High severity probability") +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.9, 0.1),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y.right = element_blank(),
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank())
fireid_pred




pg_hdw <- plot_grid(dens_hdw_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 12)),
                    ht_hdw_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 12)),
                    em_hdw_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 12)),
                    clust_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 12)),
                    area_hdw_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 12)),
                    ncol = 3)
pg_hdw

pg_wx_topo <- plot_grid(hdw_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                   fuel_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                   psev_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                   cwd_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                   tpi_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                   slope_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                   heatload_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                   nrow = 2)
pg_wx_topo

pg_topo <- plot_grid(fireid_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                     tpi_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                     slope_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                     heatload_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                     nrow = 1)
pg_topo

pg_wx <- plot_grid(hdw_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                   fuel_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                   psev_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                   cwd_pred + theme(legend.position = "none", axis.title.y = element_blank(), axis.text.y = element_blank(), axis.title = element_text(size = 10), axis.text = element_text(size = 8)),
                   nrow = 1)
pg_wx


pg_wx_topo <- plot_grid(pg_wx, pg_topo, ncol = 1)

plot_grid(pg_hdw, pg_wx_topo, ncol = 1, rel_heights = c(0.6, 0.4), align = "hv")

ggsave("plots/figureS2.pdf")
