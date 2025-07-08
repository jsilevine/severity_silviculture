
library(rstudioapi)
library(parallel)
library(data.table)
library(speedglm)

print(Sys.time())

load("bsdownership_nocwd.Rdata")
bs.coefs <- read.csv("data/bootstraps/bootstrap_coefs_ownership_nocwd.csv")

print("starting bootstrap")
set.seed(NULL)
coef.list <- run_bootstrap(1)
print(coef.list)
colnames(bs.coefs) <- colnames(coef.list)
bs.coefs <- rbind(bs.coefs, coef.list)
write.csv(bs.coefs, file = "data/bootstraps/bootstrap_coefs_ownership_nocwd.csv", row.names = FALSE)

gc()

print(paste0("completed iteration: ", nrow(bs.coefs)))

print(Sys.time())

if (nrow(bs.coefs) < 100) {

  restartSession(command='source("code/02_analysis/04_sensitivity_analyses/nocwd/run_bootstrap_ownership_nocwd.R")')

}
