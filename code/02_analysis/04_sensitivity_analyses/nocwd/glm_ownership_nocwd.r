library(data.table)
library(speedglm)
library(sf)
library(gstat)
library(parallel)

data <- fread("data/complete_data.csv")
data <- data[complete.cases(data),]

data$hdw_scaled <- scale(data$hdw)
data$avg_fuel_moisture_scaled <- scale(data$avg_fuel_moisture)
data$slope_scaled <- scale(data$slope)
data$tpi_scaled <- scale(data$tpi)
data$heat_load_scaled <- scale(data$heat_load)
data$prev_sev_scaled <- scale(data$prev_sev)

set.seed(1211)
s <- sample(1:nrow(data), round(0.25 * nrow(data)))

naive_model <- glm(hs ~ own_type + hdw_scaled + avg_fuel_moisture_scaled +
                          slope_scaled + tpi_scaled + heat_load_scaled + prev_sev_scaled + fire_name,
                  data = data[s,], family = binomial(), model = FALSE, y = FALSE)
summary(naive_model)
saveRDS(naive_model, "data/model_objects/naive_model_ownership_nocwd.rds")


## create variogram to estimate maximum range of autocorrelation
resid.rm.glm <- residuals(naive_model)
df.resid <- data.frame(z = resid.rm.glm, x = data[s,"x"], y = data[s,"y"])
df.resid <- df.resid[complete.cases(df.resid),]
set.seed(1211)
df.resid <- df.resid[sample(1:nrow(df.resid), 200000),]
v1 <- variogram(z~1, data = df.resid, locations = ~x+y, cutoff = 3000)
f1 <- fit.variogram(v1, vgm("Sph"))
max.dist <- f1$range[2]
bb <- ceiling(max.dist) ##round up to nearest whole number
saveRDS(bb, "data/autocor_scale/autocor_scale_m_ownership_nocwd.rds")

bb <- readRDS("data/autocor_scale/autocor_scale_m_ownership_nocwd.rds")

## get list of objectids to stratify over
objectids <- unique(data$fire_name)

## create empty dataframe for potential sample points
sample_points <- data.frame(objectid = objectids,
                            xmin = rep(0, times = length(objectids)),
                            xmax = rep(0, times = length(objectids)),
                            ymin = rep(0, times = length(objectids)),
                            ymax = rep(0, times = length(objectids)),
                            num_boxes = rep(0, times = length(objectids)))

## loop over individual fires:
for (i in 1:length(objectids)) {

  ## determine minimum rectangle
  fire_area <- list(xmin = min(data[fire_name == objectids[i], "x"]),
                    xmax = max(data[fire_name == objectids[i], "x"]),
                    ymin = min(data[fire_name == objectids[i], "y"]),
                    ymax = max(data[fire_name == objectids[i], "y"]))

  ## divide study area into boxes of size bb
  box_centers <- expand.grid(x = seq(fire_area$xmin, fire_area$xmax, by = bb) + bb/2,
                             y = seq(fire_area$ymin, fire_area$ymax, by = bb) + bb/2)


  num_boxes <- nrow(box_centers)

  sample_points[i, 2] <- min(box_centers$x)
  sample_points[i, 3] <- max(box_centers$x)
  sample_points[i, 4] <- min(box_centers$y)
  sample_points[i, 5] <- max(box_centers$y)
  sample_points[i, 6] <- num_boxes

}

## req fields for analysis:
fields <- c("x",
            "y",
            "hs",
            "own_type",
            "hdw_scaled",
            "avg_fuel_moisture_scaled",
            "slope_scaled",
            "tpi_scaled",
            "heat_load_scaled",
            "prev_sev_scaled",
            "fire_name")

## function to pull data given spatial extent, fire name, and fields
pull_data <- function(xcoord, ycoord, data, bb, objid) {

  p.data <- data[fire_name == objid]
  setkey(p.data, x)
  p.data <- p.data[x > (xcoord-(bb/2)) & x < (xcoord+(bb/2))]
  setkey(p.data, y)
  p.data <- p.data[y > (ycoord-(bb/2)) & y < (ycoord+(bb/2))]

  return(p.data)

}

## function to randomly choose boxes and pull data:
rpull <- function(nbox, objectid, sdata, bb, sample_points) {

  xcoords <- runif(nbox, sample_points[sample_points$objectid == objectid, "xmin"], sample_points[sample_points$objectid == objectid, "xmax"])
  ycoords <- runif(nbox, sample_points[sample_points$objectid == objectid, "ymin"], sample_points[sample_points$objectid == objectid, "ymax"])

  bs.data <- mcmapply(FUN = pull_data,
                      xcoords,
                      ycoords,
                      MoreArgs = list(data = sdata, bb = bb, objid = objectid),
                      SIMPLIFY = FALSE,
                      mc.cores = 7)
  out <- do.call(rbind, bs.data)

  return(out)

}

## set number of bootstraps
n <- 300

## create empty dataframe of fitted coefficients:
cols <- rownames(coef(summary(naive_model)))
bs.coefs <- data.frame(matrix(NA, nrow = 1, ncol = length(cols)))
colnames(bs.coefs) <- cols
bs.coefs <- bs.coefs[0,] 

## write function to create bootstrap dataset, run model, append estimates to dataframe
run_bootstrap <- function(iter) {

  system.time({
  outlist <- list()
  for (i in 1:nrow(sample_points)) {

    print(paste0("starting iteration: ", i))

    outlist[[i]] <- rpull(nbox = sample_points[i, "num_boxes"],
                          objectid = sample_points[i, "objectid"],
                          sdata = data[s,..fields], bb = bb, sample_points = sample_points)

  }
  bs.data <- do.call(rbind, outlist)
  })

  print("starting glm")

  bs.glm <- speedglm(hs ~ own_type + hdw_scaled + avg_fuel_moisture_scaled +
                          slope_scaled + tpi_scaled + heat_load_scaled + prev_sev_scaled + fire_name,
                      data = bs.data,
                      family = binomial(),
                      y = FALSE,
                      model = FALSE,
                      fitted = FALSE)

  coefs <- as.data.frame(matrix(coef(summary(bs.glm))[, "Estimate"],
                                nrow = 1))
  rm(bs.glm)
  gc()

  colnames(coefs) <- cols

  return(coefs)

}

## run bootstraps, can be done in parallel but my poor computer could not handle it
## a note - when doing this with speedglm I had memory clearing issues. gc() did not remove all the memory, and so in order to fully clear the memory between runs I had to restart R. I did so by looping on an alternative script which included calls to restart R after each iteration. (See below)
for(i in 1:1) {

  if (i == 1) print("starting bootstrap")

  coef.list <- mclapply(1, run_bootstrap, mc.cores = 1)
  n.coefs <- as.data.frame(do.call(rbind, coef.list))
  bs.coefs <- rbind(bs.coefs, n.coefs)

  write.csv(bs.coefs, file = "data/bootstraps/bootstrap_coefs_ownership_nocwd.csv")

  gc()

  print(paste0("completed iteration: ", i))

  print(Sys.time())
}


write.csv(bs.coefs, file = "data/bootstraps/bootstrap_coefs_ownership_nocwd.csv", row.names = FALSE)

## this is the janky way to force R to restart after each iteration
## first save workspace image
save.image(file = "bsdownership_nocwd.Rdata")
## then run external file (contents below)
source("code/02_analysis/04_sensitivity_analyses/nocwd/run_bootstrap_ownership_nocwd.R")
