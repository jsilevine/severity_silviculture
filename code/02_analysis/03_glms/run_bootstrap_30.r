
library(rstudioapi)
library(parallel)
library(data.table)
library(speedglm)

print(Sys.time())

load("bsd.Rdata")
bs.coefs <- read.csv("data/bootstraps/bootstrap_coefs_30.csv")

print("starting bootstrap")
set.seed(NULL)
coef.list <- run_bootstrap(1)
print(coef.list)
colnames(bs.coefs) <- colnames(coef.list)
bs.coefs <- rbind(bs.coefs, coef.list)
write.csv(bs.coefs, file = "data/bootstraps/bootstrap_coefs_30.csv", row.names = FALSE)

gc()

print(paste0("completed iteration: ", nrow(bs.coefs)))

print(Sys.time())

if (nrow(bs.coefs) < 200) {

  restartSession(command='source("code/02_analysis/run_bootstrap_30.R")')

}
