##### QLG Project #####

### Packages
library(ggplot2)
library(raster)
library(sf)
library(chron)
library(ggnewscale)
library(cowplot)

sort(sapply(ls(),function(x){object.size(get(x))}))

baseline <- read.csv("data/qlg/baseline_treated.csv")

new_pred <- read.csv("data/qlg/results_prob.csv")
prediction <- read.csv("data/qlg/pred.csv")
prediction <- cbind(prediction, rowMeans(new_pred, na.rm = TRUE))
colnames(prediction)[9] <- "new_pred_mean"
prediction <- cbind(prediction, baseline[,"pr"])
colnames(prediction)[10] <- "baseline_pred"
prediction$sev_dif <- prediction$new_pred - prediction$baseline_pred

new_pred

p_dif <- numeric(ncol(new_pred)-1)
for (i in 2:ncol(new_pred)) {

  n_p <- mean(new_pred[prediction$qlg_treatment_type != "none", i])
  b_p <- mean(prediction[prediction$qlg_treatment_type != "none", "baseline_pred"])

  p_dif[i] <- n_p - b_p

}
p_dif <- p_dif[2:ncol(new_pred)]

p_dif

mean(p_dif)
quantile(p_dif, c(0.05, 0.95))

a_dif <- numeric(ncol(new_pred))
for (i in 1:ncol(new_pred)) {

  n_area <- length(new_pred[prediction$qlg_treatment_type != "none" & new_pred[,i] >= 0.5, i])
  b_area <- nrow(prediction[prediction$qlg_treatment_type != "none" & prediction$baseline_pred >= 0.5, ])

  a_dif[i] <- (b_area - n_area) * (30^2) / 10000

}
a_dif <- a_dif[2:ncol(new_pred)]

mean(a_dif)
quantile(a_dif, c(0.05, 0.95))

nrow(prediction[prediction$qlg_treatment_type != "none", ]) * (30^2) / 10000
nrow(prediction[prediction$qlg_treatment_type != "none" & prediction$sev_dif > 0, ]) * (30^2) / 10000
nrow(prediction[prediction$qlg_treatment_type != "none" & prediction$sev_dif < 0, ]) * (30^2) / 10000


mean(prediction[prediction$qlg_treatment_type != "none" & prediction$sev_dif > 0, "sev_dif"])
quantile(prediction[prediction$qlg_treatment_type != "none" & prediction$sev_dif > 0, "sev_dif"], c(0.05, 0.95))

mean(prediction[prediction$qlg_treatment_type != "none" & prediction$sev_dif < 0, "sev_dif"])
quantile(prediction[prediction$qlg_treatment_type != "none" & prediction$sev_dif < 0, "sev_dif"], c(0.05, 0.95))

mean(prediction[prediction$qlg_treatment_type != "none" & prediction$sev_dif > 0, "baseline_pred"])
quantile(prediction[prediction$qlg_treatment_type != "none" & prediction$sev_dif > 0, "baseline_pred"], c(0.05, 0.95))
mean(prediction[prediction$qlg_treatment_type != "none" & prediction$sev_dif < 0, "baseline_pred"])
quantile(prediction[prediction$qlg_treatment_type != "none" & prediction$sev_dif < 0, "baseline_pred"], c(0.05, 0.95))

nrow(prediction[prediction$qlg_treatment_type != "none" & prediction$new_pred >= 0.5 & prediction$baseline_pred <= 0.5, ]) * (30^2) / 10000





mean(prediction[prediction$qlg_treatment_type != "none", "new_pred"])
mean(prediction[prediction$qlg_treatment_type != "none", "baseline_pred"])

mean(prediction[prediction$qlg_treatment_type != "none", "sev_dif"], na.rm = TRUE)

ggplot() +
  geom_raster(data = prediction[prediction$qlg_treatment_type != "none" & prediction$sev_dif > 0,],
              aes(x = x, y = y, fill = sev_dif)) +
  #geom_sf(data = st_make_valid(existing_trt_all), color = "black", fill = NA, alpha = 1.0) +
  #geom_sf(data = st_make_valid(proposed_dfpz_all), color = "green", fill = NA, alpha = 1.0) +
  #geom_sf(data = st_make_valid(proposed_matrix_all), color = "pink", fill = NA, alpha = 1.0) +
  scale_fill_gradient2(low = "#fff5f0", mid = "#fb6a4a", high = "#a50f15") +
  new_scale_fill() +
  geom_raster(data = prediction[prediction$qlg_treatment_type != "none" & prediction$sev_dif < 0,],
              aes(x = x, y = y, fill = sev_dif)) +
  scale_fill_gradient2(low = "#08519c", mid = "#6baed6", high = "#f7fbff") +
  geom_raster(data = prediction[prediction$qlg_treatment_type == "none",],
              aes(x = x, y = y), fill = "lightgray", alpha = 0.4) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave("prob_plot.pdf")


p1 <- ggplot() +
  geom_raster(data = prediction[prediction$qlg_treatment_type != "none",],
              aes(x = x, y = y, fill = baseline_pred)) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0.5) +
  geom_raster(data = prediction[prediction$qlg_treatment_type == "none",],
              aes(x = x, y = y), fill = "lightgray", alpha = 0.4) +
  xlab("") +
  ylab("") +
  ggtitle("Baseline prediction") +
  labs(fill = "Probabaility of \n high severity fire") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
leg <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

p2 <- ggplot() +
  geom_raster(data = prediction[prediction$qlg_treatment_type != "none",],
              aes(x = x, y = y, fill = new_pred_mean)) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0.5) +
  geom_raster(data = prediction[prediction$qlg_treatment_type == "none",],
              aes(x = x, y = y), fill = "lightgray", alpha = 0.4) +
  ggtitle("Treatment prediction") +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

plot_grid(p1, p2, leg, nrow = 1, rel_widths = c(0.45, 0.45, 0.1))
ggsave("comparison_plot.pdf")



## create fake data boxes to show total area comparisons
n_decr <- nrow(prediction[prediction$qlg_treatment_type != "none" & prediction$new_pred <= 0.5 &
                          prediction$baseline_pred >= 0.5,])
n_incr <- nrow(prediction[prediction$qlg_treatment_type != "none" & prediction$new_pred >= 0.5 &
                          prediction$baseline_pred <= 0.5,])
n_unch <- nrow(prediction[prediction$qlg_treatment_type != "none" & prediction$new_pred >= 0.5 &
                          prediction$baseline_pred >= 0.5,]) +
  nrow(prediction[prediction$qlg_treatment_type != "none" & prediction$new_pred <= 0.5 &
                          prediction$baseline_pred <= 0.5,])

nrow <- 150
ncol_decr <- floor(n_decr / nrow)
ncol_incr <- floor(n_incr / nrow)
ncol_unch <- floor(n_unch / nrow)

extra_row_decr <- n_decr %% nrow
extra_row_incr <- n_incr %% nrow
extra_row_unch <- n_unch %% nrow

min(prediction$y)

row_y_incr <- unique(prediction$y)[order(unique(prediction$y))][1:nrow]
row_y_decr <- unique(prediction$y)[order(unique(prediction$y))][(nrow+51):((nrow*2)+50)]
row_y_unch <- unique(prediction$y)[order(unique(prediction$y))][((nrow*2)+101):((nrow*3)+100)]
col_options <- unique(prediction$x)[order(unique(prediction$x))][2598:(2598+ncol_unch)]

y_dif <- unique(prediction$y)[order(unique(prediction$y))][2] - unique(prediction$y)[order(unique(prediction$y))][1]

decr_data <- data.frame(x = numeric(0), y = numeric(0))
for (i in 1:ncol_decr) {
  ndata <- data.frame(x = rep(col_options[i], times = nrow), y = row_y_decr)
  decr_data <- rbind(decr_data, ndata)
}
ndata <- data.frame(x = rep(col_options[i], times = extra_row_decr), y = row_y_decr[1:extra_row_decr])
decr_data <- rbind(decr_data, ndata)


incr_data <- data.frame(x = numeric(0), y = numeric(0))
for (i in 1:ncol_incr) {
  ndata <- data.frame(x = rep(col_options[i], times = nrow), y = row_y_incr)
  incr_data <- rbind(incr_data, ndata)
}
ndata <- data.frame(x = rep(col_options[i], times = extra_row_incr), y = row_y_incr[1:extra_row_incr])
incr_data <- rbind(incr_data, ndata)


#unch_data <- data.frame(x = numeric(0), y = numeric(0))
#for (i in 1:ncol_unch) {
#  ndata <- data.frame(x = rep(col_options[i], times = nrow), y = row_y_unch)
#  unch_data <- rbind(unch_data, ndata)
#}
#ndata <- data.frame(x = rep(col_options[i], times = extra_row_unch), y = row_y_unch[1:extra_row_unch])
#unch_data <- rbind(unch_data, ndata)


ggplot() +
  geom_raster(data = prediction[prediction$qlg_treatment_type != "none" & prediction$new_pred >= 0.5 &
                                prediction$baseline_pred <= 0.5,],
              aes(x = x, y = y), fill = "#a50f15") +
  geom_raster(data = incr_data,
              aes(x = x, y = y), fill = "#a50f15") +
  geom_raster(data = prediction[prediction$qlg_treatment_type != "none" & prediction$new_pred <= 0.5 &
                                prediction$baseline_pred >= 0.5,],
              aes(x = x, y = y), fill = "#08519c") +
  geom_raster(data = decr_data,
              aes(x = x, y = y), fill = "#08519c") +
  geom_raster(data = prediction[prediction$qlg_treatment_type != "none" & prediction$new_pred >= 0.5 &
                                prediction$baseline_pred >= 0.5,],
              aes(x = x, y = y), fill = "black", alpha = 0.2) +
  geom_raster(data = prediction[prediction$qlg_treatment_type != "none" & prediction$new_pred <= 0.5 &
                                prediction$baseline_pred <= 0.5,],
              aes(x = x, y = y), fill = "black", alpha = 0.2) +
  geom_raster(data = prediction[prediction$qlg_treatment_type == "none",],
              aes(x = x, y = y), fill = "lightgray", alpha = 0.4) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave("area_plot.pdf")

prediction <- merge(prediction, baseline, by = c("x", "y"), all.x = TRUE)

prediction$sev_dif <- prediction$new_pred - prediction$pr

mean(prediction$sev_dif, na.rm = TRUE)

mean(prediction$pr, na.rm = TRUE)
mean(prediction[prediction$sev_dif != 0, "sev_dif"], na.rm = TRUE)

mean(baseline$pr)


mean(prediction[prediction$qlg_treatment_type != "none" & prediction$sev_dif != 0, "sev_dif"], na.rm = TRUE)









