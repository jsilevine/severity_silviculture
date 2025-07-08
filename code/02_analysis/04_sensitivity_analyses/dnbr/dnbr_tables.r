library(speedglm)
library(data.table)
library(ggplot2)
library(cowplot)
library(paletteer)
library(MetBrewer)

data <- fread("data/complete_data.csv")
data <- data[complete.cases(data),]

br_own <- read.csv("data/bootstraps/bootstrap_coefs_ownership_dnbr.csv")

head(br_own)
colnames(br_own)[1] <- "intercept"

results_summary <- data.frame(variable = colnames(br_own),
                              mean = numeric(ncol(br_own)),
                              ci_lower = numeric(ncol(br_own)),
                              ci_upper = numeric(ncol(br_own)))

for (i in 1:nrow(results_summary)) {

  results_summary[i,"mean"] <- mean(br_own[,results_summary[i,"variable"]])
  results_summary[i,"ci_lower"] <- quantile(br_own[,results_summary[i,"variable"]], 0.025)
  results_summary[i,"ci_upper"] <- quantile(br_own[,results_summary[i,"variable"]], 0.975)

}

write.csv(results_summary, "data/model_objects/ownership_dnbr_summary.csv", row.names = FALSE)



br_30 <- read.csv("data/bootstraps/bootstrap_coefs_30_dnbr.csv")

colnames(br_30)[1] <- "intercept"

results_summary_30 <- data.frame(variable = colnames(br_30),
                                 mean = numeric(ncol(br_30)),
                                 ci_lower = numeric(ncol(br_30)),
                                 ci_upper = numeric(ncol(br_30)))

for (i in 1:nrow(results_summary_30)) {
    results_summary_30[i,"mean"] <- mean(br_30[,results_summary_30[i,"variable"]])
    results_summary_30[i,"ci_lower"] <- quantile(br_30[,results_summary_30[i,"variable"]], 0.025)
    results_summary_30[i,"ci_upper"] <- quantile(br_30[,results_summary_30[i,"variable"]], 0.975)
}

write.csv(results_summary_30, "data/model_objects/30_dnbr_summary.csv", row.names = FALSE)


br_180 <- read.csv("data/bootstraps/bootstrap_coefs_180_dnbr.csv")

colnames(br_180)[1] <- "intercept"

results_summary_180 <- data.frame(variable = colnames(br_180),
                                 mean = numeric(ncol(br_180)),
                                 ci_lower = numeric(ncol(br_180)),
                                 ci_upper = numeric(ncol(br_180)))

for (i in 1:nrow(results_summary_180)) {
    results_summary_180[i,"mean"] <- mean(br_180[,results_summary_180[i,"variable"]])
    results_summary_180[i,"ci_lower"] <- quantile(br_180[,results_summary_180[i,"variable"]], 0.025)
    results_summary_180[i,"ci_upper"] <- quantile(br_180[,results_summary_180[i,"variable"]], 0.975)
}

write.csv(results_summary_180, "data/model_objects/180_dnbr_summary.csv", row.names = FALSE)
