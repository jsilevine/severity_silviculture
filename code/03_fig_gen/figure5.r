library(speedglm)
library(data.table)
library(ggplot2)
library(cowplot)
library(paletteer)
library(MetBrewer)

bootstrap_results <- read.csv("data/bootstraps/bootstrap_coefs_ownership.csv")
naive_model <- readRDS("data/model_objects/naive_model_ownership.rds")

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

write.csv(results_summary, "data/model_objects/ownership_results_summary.csv")


summary(naive_model)


pdata <- data.frame(ownership = c("public", "private_industrial", "other"),
                    mean = numeric(3), ci_lower = numeric(3), ci_upper = numeric(3))

pdata[pdata$ownership == "public", "mean"] <- exp(naive_model$coefficients[1]) / (1 + exp(naive_model$coefficients[1]))
pdata[pdata$ownership == "private_industrial", "mean"] <- exp(naive_model$coefficients[1] + naive_model$coefficients[3]) /
  (1 + exp(naive_model$coefficients[1] + naive_model$coefficients[3]))
pdata[pdata$ownership == "other", "mean"] <- exp(naive_model$coefficients[1] + naive_model$coefficients[2]) /
  (1 + exp(naive_model$coefficients[1] + naive_model$coefficients[2]))

pdata[pdata$ownership == "public", "ci_lower"] <- exp(results_summary[1,"ci_lower"]) / (1 + exp(results_summary[1,"ci_lower"]))
pdata[pdata$ownership == "private_industrial", "ci_lower"] <- exp(naive_model$coefficients[1] + results_summary[3,"ci_lower"]) /
  (1 + exp(naive_model$coefficients[1] + results_summary[3,"ci_lower"]))
pdata[pdata$ownership == "other", "ci_lower"] <- exp(naive_model$coefficients[1] + results_summary[2,"ci_lower"]) /
  (1 + exp(naive_model$coefficients[1] + results_summary[2,"ci_lower"]))

pdata[pdata$ownership == "public", "ci_upper"] <- exp(results_summary[1,"ci_upper"]) / (1 + exp(results_summary[1,"ci_upper"]))
pdata[pdata$ownership == "private_industrial", "ci_upper"] <- exp(naive_model$coefficients[1] + results_summary[3,"ci_upper"]) /
  (1 + exp(naive_model$coefficients[1] + results_summary[3,"ci_upper"]))
pdata[pdata$ownership == "other", "ci_upper"] <- exp(naive_model$coefficients[1] + results_summary[2,"ci_upper"]) /
  (1 + exp(naive_model$coefficients[1] + results_summary[2,"ci_upper"]))



pdata$ownership <- c("Public", "Private Industrial", "Other")

pdata$ownership <- factor(pdata$ownership, levels = c("Public", "Private Industrial", "Other"))

colors <- met.brewer(name = "Egypt", n = 6)
colors
colors <- colors[c(1,3,4)]
colors

colors2 <- met.brewer(name = "Derain", n = 6)
colors2

ggplot(pdata, aes(x = ownership)) +
  geom_point(aes(y = mean, color = ownership), size = 8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper, color = ownership), linewidth = 1.5, width = 0.3) +
  xlab("Ownership") +
  ylab("High-severity fire probability") +
  scale_color_manual(name = "Ownership", values = c(colors2[1], colors[2], colors[3])) +
  guides(x.sec = "axis", y.sec = "axis") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
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

ggsave("plots/figure5.png")
ggsave("plots/figure5.pdf")
