
library(ggplot2)

shap_30_overall <- read.csv("data/shapley_total_30m.csv")
shap_30_private <- read.csv("data/shapley_private_30m.csv")
shap_30_public <- read.csv("data/shapley_public_30m.csv")
shap_30_balanced <- read.csv("data/shapley_balanced_30m.csv")

shap_180_overall <- read.csv("data/shapley_total_180m.csv")
shap_180_private <- read.csv("data/shapley_private_180m.csv")
shap_180_public <- read.csv("data/shapley_public_180m.csv")
shap_180_balanced <- read.csv("data/shapley_balanced_180m.csv")

plot_shapley <- function(data, title = "") {

  data[data$feature == "hdw_scaled", "feature"] <- "HDW index"
  data[data$feature == "mean_dens_30_scaled", "feature"] <- "Stem density"
  data[data$feature == "mean_dens_180_scaled", "feature"] <- "Stem density"
  data[data$feature == "cwd_scaled", "feature"] <- "cwd"
  data[data$feature == "em_scaled", "feature"] <- "Vertical fuel continuity"
  data[data$feature == "avg_fuel_moisture_scaled", "feature"] <- "Fuel moisture"
  data[data$feature == "prev_sev_scaled", "feature"] <- "CBI (previous timestep)"
  data[data$feature == "mean_area_30_scaled", "feature"] <- "Mean gap area"
  data[data$feature == "mean_area_180_scaled", "feature"] <- "Mean gap area"
  data[data$feature == "fire_name", "feature"] <- "Fire ID"
  data[data$feature == "slope_scaled", "feature"] <- "Slope"
  data[data$feature == "tpi_scaled", "feature"] <- "TPI"
  data[data$feature == "mean_ht_30_scaled", "feature"] <- "Mean stem height"
  data[data$feature == "mean_ht_180_scaled", "feature"] <- "Mean stem height"
  data[data$feature == "clust_30_scaled", "feature"] <- "Stem clustering"
  data[data$feature == "clust_180_scaled", "feature"] <- "Stem clustering"
  data[data$feature == "heat_load_scaled", "feature"] <- "Heat load"

  data$feature <- factor(data$feature, levels = data$feature[order(data$mean_high, decreasing = FALSE)])

  p <- ggplot(data) +
    geom_point(aes(x = mean_high, y = feature), size = 8) +
    theme_bw() +
    ggtitle(title) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(data$mean_high) + 0.05*max(data$mean_high)),
                       breaks = c(0, 0.05, 0.1)) +
    xlab("Shapley feature importance") +
    ylab("Variable") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          plot.title = element_blank(),
          axis.text = element_text(size = 20))
  for (i in 1:nrow(data)) {
    p <- p + geom_segment(x = 0,
                          xend = data[i,"mean_high"],
                          y = data[i,"feature"],
                          yend = data[i,"feature"],
                          linetype = "dashed", linewidth = 1)
  }
  return(p)

}

plot_shapley(shap_30_overall, "30m overall")
ggsave("plots/figure3/figure3a.pdf", width = 8, height = 5)
plot_shapley(shap_30_private, "30m private")
ggsave("plots/figure3/figure3b.pdf", width = 8, height = 5)
plot_shapley(shap_30_public, "30m public")
ggsave("plots/figure3/figure3c.pdf", width = 8, height = 5)
plot_shapley(shap_30_balanced, "30m balanced")
ggsave("plots/figure3/figure3d.pdf", width = 8, height = 5)

plot_shapley(shap_180_overall, "180m overall")
ggsave("plots/figure3/figure3e.pdf", width = 8, height = 5)
plot_shapley(shap_180_private, "180m private")
ggsave("plots/figure3/figure3f.pdf", width = 8, height = 5)
plot_shapley(shap_180_public, "180m public")
ggsave("plots/figure3/figure3g.pdf", width = 8, height = 5)
plot_shapley(shap_180_balanced, "180m balanced")
ggsave("plots/figure3/figure3h.pdf", width = 8, height = 5)
