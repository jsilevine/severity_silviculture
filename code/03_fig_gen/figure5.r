
##---------------------------------------------------------------
## Density plots comparing private and public land ownership
##---------------------------------------------------------------
library(data.table)
library(ggplot2)
library(cowplot)

data <- fread("data/complete_data.csv")
data <- data[complete.cases(data),]
set.seed(1)

dsub_private <- data[own_type == "Private Industrial"]
#dsub_private <- dsub_private[sample(1:nrow(dsub_private), 598409)]
dsub_public <- data[own_type == "Federal"]
#dsub_public <- dsub_public[sample(1:nrow(dsub_public), 598409)]

dprivate <- density(dsub_private$mean_dens_30, n = 10000, from = min(data$mean_dens_30), to = max(data$mean_dens_30))
dpublic <- density(dsub_public$mean_dens_30, n = 10000, from = min(data$mean_dens_30), to = max(data$mean_dens_30))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
dens_30 <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0 & pd$x > 100,], linewidth = 2, color = "#08519c") +
  geom_line(data = pd[pd$dif >= 0 & pd$x < 100,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 350)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Stem density (45m)") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
dens_30

dprivate <- density(dsub_private$clust_30, n = 10000, from = min(data$clust_30), to = max(data$clust_30))
dpublic <- density(dsub_public$clust_30, n = 10000, from = min(data$clust_30), to = max(data$clust_30))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
clust_30 <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0 & pd$x < 1,], linewidth = 2, color = "#08519c") +
  geom_line(data = pd[pd$dif >= 0 & pd$x > 1,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0), limits = c(0.5, 1.5)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Spatial regularity (45m)") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
clust_30


dprivate <- density(dsub_private$mean_ht_30, n = 10000, from = min(data$mean_ht_30), to = max(data$mean_ht_30))
dpublic <- density(dsub_public$mean_ht_30, n = 10000, from = min(data$mean_ht_30), to = max(data$mean_ht_30))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
ht_30 <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0 & pd$x < 20,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 20,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0), limits = c(5, 45)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Mean stem height (45m)") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
ht_30

dprivate <- density(dsub_private$mean_area_30, n = 40000, from = min(data$mean_area_30), to = max(data$mean_area_30))
dpublic <- density(dsub_public$mean_area_30, n = 40000, from = min(data$mean_area_30), to = max(data$mean_area_30))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
area_30 <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 2000)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Mean gap area (45m)") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
area_30

dprivate <- density(dsub_private$mean_dens_180, n = 10000, from = min(data$mean_dens_180), to = max(data$mean_dens_180))
dpublic <- density(dsub_public$mean_dens_180, n = 10000, from = min(data$mean_dens_180), to = max(data$mean_dens_180))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
dens_180 <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0 & pd$x < 300,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0 & pd$x < 100,], linewidth = 2, color = "#08519c") +
  geom_line(data = pd[pd$dif >= 0 & pd$x >= 100,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 350)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Stem density (390m)") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
dens_180

dprivate <- density(dsub_private$clust_180, n = 10000, from = min(data$clust_180), to = max(data$clust_180))
dpublic <- density(dsub_public$clust_180, n = 10000, from = min(data$clust_180), to = max(data$clust_180))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
clust_180 <- ggplot(data = pd, aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  #geom_line(data = pd, linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0,], linewidth = 2, color = "#08519c") +
  #geom_line(data = pd[pd$dif >= 0 & pd$x > 1,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0), limits = c(0.5, NA)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Spatial regularity (390m)") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
clust_180

dprivate <- density(dsub_private$mean_ht_180, n = 10000, from = min(data$mean_ht_180), to = max(data$mean_ht_180))
dpublic <- density(dsub_public$mean_ht_180, n = 10000, from = min(data$mean_ht_180), to = max(data$mean_ht_180))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
ht_180 <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0 & pd$x < 20,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 20,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0), limits = c(5, 45)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Mean stem height (390m)") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
ht_180

dprivate <- density(dsub_private$mean_area_30, n = 40000, from = min(data$mean_area_180), to = max(data$mean_area_180))
dpublic <- density(dsub_public$mean_area_30, n = 40000, from = min(data$mean_area_180), to = max(data$mean_area_180))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
area_180 <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 2000)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Mean gap area (390m)") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
area_180


dprivate <- density(dsub_private$em, n = 10000, from = min(data$em), to = max(data$em))
dpublic <- density(dsub_public$em, n = 10000, from = min(data$em), to = max(data$em))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
em <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0 & pd$x > 10,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Ladder fuels index") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
em


dprivate <- density(dsub_private$hdw, n = 10000, from = min(data$hdw), to = max(data$hdw))
dpublic <- density(dsub_public$hdw, n = 10000, from = min(data$hdw), to = max(data$hdw))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
hdw <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0 & pd$x < 100,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0 & pd$x > 10,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Hot-dry-windy index") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
hdw

dprivate <- density(dsub_private$avg_fuel_moisture, n = 30000, from = min(data$avg_fuel_moisture), to = max(data$avg_fuel_moisture))
dpublic <- density(dsub_public$avg_fuel_moisture, n = 30000, from = min(data$avg_fuel_moisture), to = max(data$avg_fuel_moisture))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
fuels <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0 & pd$x < 4,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 4 & pd$x < 5,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 5 & pd$x < 7,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 7,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0 & pd$x < 3,], linewidth = 2, color = "#08519c") +
  geom_line(data = pd[pd$dif >= 0 & pd$x > 3 & pd$x < 4.5,], linewidth = 2, color = "#08519c") +
  geom_line(data = pd[pd$dif >= 0 & pd$x > 4.5 & pd$x < 5.5,], linewidth = 2, color = "#08519c") +
  geom_line(data = pd[pd$dif >= 0 & pd$x > 5.5 & pd$x < 10,], linewidth = 2, color = "#08519c") +
  geom_line(data = pd[pd$dif >= 0 & pd$x > 10,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0), limits = c(0, 15)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Mean fuel moisture (%)") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
fuels


dprivate <- density(dsub_private$cwd, n = 10000, from = min(data$cwd), to = max(data$cwd))
dpublic <- density(dsub_public$cwd, n = 10000, from = min(data$cwd), to = max(data$cwd))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
cwd <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 300 & pd$x < 430,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 480 & pd$x < 520,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 520 & pd$x < 550,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 550 & pd$x < 600,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 600 & pd$x < 690,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 695,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif > 0 & pd$x < 400,], linewidth = 2, color = "#08519c") +
  geom_line(data = pd[pd$dif > 0 & pd$x > 400 & pd$x < 500,], linewidth = 2, color = "#08519c") +
  geom_line(data = pd[pd$dif > 0 & pd$x > 500 & pd$x < 530,], linewidth = 2, color = "#08519c") +
  geom_line(data = pd[pd$dif > 0 & pd$x > 530 & pd$x < 550,], linewidth = 2, color = "#08519c") +
  geom_line(data = pd[pd$dif > 0 & pd$x > 585 & pd$x < 650,], linewidth = 2, color = "#08519c") +
  geom_line(data = pd[pd$dif > 0 & pd$x > 650 & pd$x < 700,], linewidth = 2, color = "#08519c") +
  geom_line(data = pd[pd$dif > 0 & pd$x > 700,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Climate water deficit") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
cwd


dprivate <- density(dsub_private$tpi, n = 10000, from = min(data$tpi), to = max(data$tpi))
dpublic <- density(dsub_public$tpi, n = 10000, from = min(data$tpi), to = max(data$tpi))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
tpi <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0 & pd$x < 0,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 0,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0 & pd$x < 25,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Topographic position index") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
tpi


dprivate <- density(dsub_private$slope, n = 10000, from = min(data$slope), to = max(data$slope))
dpublic <- density(dsub_public$slope, n = 10000, from = min(data$slope), to = max(data$slope))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
slope <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0 & pd$x < 20,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 20,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Slope (%)") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
slope


dprivate <- density(dsub_private$heat_load, n = 10000, from = min(data$heat_load), to = max(data$heat_load))
dpublic <- density(dsub_public$heat_load, n = 10000, from = min(data$heat_load), to = max(data$heat_load))
pd <- data.frame(x = dprivate$x, dif = dprivate$y - dpublic$y)
heat_load <- ggplot(data = pd[pd$dif < 0,], aes(x = x, y = dif)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = pd[pd$dif < 0 & pd$x < 0.8,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif < 0 & pd$x > 0.8,], linewidth = 2, color = "#2ca25f") +
  geom_line(data = pd[pd$dif >= 0,], linewidth = 2, color = "#08519c") +
  scale_x_continuous(expand = c(0,0)) +
  #scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Heat Load") +
  ylab("Density (Private Industrial - Public)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank(),
    axis.title.y = element_blank())
heat_load

plot_grid(dens_30, clust_30, ht_30, area_30,
          dens_180, clust_180, ht_180, area_180,
          em, hdw, fuels,
          cwd, tpi, slope, heat_load, align = "hv")

ggsave("plots/figure5/figure5b.pdf")


##---------------------------------------------------------------
##
##---------------------------------------------------------------


dens30 <- ggplot(data = dsub_private, aes(x = mean_dens_30)) +
  geom_density(linewidth = 2, fill = "#6baed6", color = "#08519c", alpha =0.6) +
  geom_density(data = dsub_public, linewidth = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 350)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.011)) +
  xlab("Stem density (45m)") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())


clust30 <- ggplot(data = dsub_private, aes(x = clust_30)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha =0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0), limits = c(0.5, 1.5)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 6.2)) +
  xlab("Stem clustering (45m)") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())
clust30

ht30 <- ggplot(data = dsub_private, aes(x = mean_ht_30)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha =0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0), limits = c(5, 45)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.07)) +
  xlab("Mean tree height (45m)") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())

area30 <- ggplot(data = dsub_private, aes(x = mean_area_30)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha =0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 3500)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.005)) +
  xlab("Mean gap area (45m)") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())


dens180 <- ggplot(data = dsub_private, aes(x = mean_dens_180)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha =0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 350)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.014)) +
  xlab("Stem density (195m)") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())

clust180 <- ggplot(data = dsub_private, aes(x = clust_180)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha =0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0), limits = c(0.5, 1.5)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 7.5)) +
  xlab("Stem clustering (195m)") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())
clust180

ht180 <- ggplot(data = dsub_private, aes(x = mean_ht_180)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha =0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0), limits = c(5, 45)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.09)) +
  xlab("Mean tree height (195m)") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())

area180 <- ggplot(data = dsub_private, aes(x = mean_area_180)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha =0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 3500)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.0055)) +
  xlab("Mean gap area (195m)") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())


em <- ggplot(data = dsub_private, aes(x = em)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha =0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.05)) +
  xlab("Vertical fuel continuity") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())

hdw <- ggplot(data = dsub_private, aes(x = hdw)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha =0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.018)) +
  xlab("Hot dry windy index") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())

fuels <- ggplot(data = dsub_private, aes(x = avg_fuel_moisture)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha =0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 15)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 1.1)) +
  xlab("Mean fuel moisture") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())

elev <- ggplot(data = dsub_private, aes(x = elevation)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha = 0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.0018)) +
  xlab("Elevation") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())

tpi <- ggplot(data = dsub_private, aes(x = tpi)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha = 0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.045)) +
  xlab("Topographic position index") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())

slope <- ggplot(data = dsub_private, aes(x = slope)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha = 0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 0.05)) +
  xlab("Slope (%)") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())

heat_load <- ggplot(data = dsub_private, aes(x = heat_load)) +
  geom_density(size = 2, fill = "#6baed6", color = "#08519c", alpha = 0.6) +
  geom_density(data = dsub_public, size = 2, fill = "#66c2a4", color = "#2ca25f", alpha = 0.6) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0.0), limits = c(0, 3.5)) +
  xlab("Heat load") +
  ylab("") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 12),
    axis.line.x.top = element_line(size = 1),
    axis.line.x.bottom = element_line(size = 1),
    axis.line.y.left = element_line(size = 1),
    axis.line.y.right = element_line(size = 1),
    axis.text.y = element_blank())


plot_grid(dens30, clust30, ht30, area30,
          dens180, clust180, ht180, area180,
          em, hdw, fuels,
          elev, tpi, slope, heat_load, align = "hv")

ggsave("plots/figure5/figure5.pdf")
