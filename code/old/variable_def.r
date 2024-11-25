##---------------------------------------------------------------
## variable_def.r
## DEFINE GLOBAL VARIABLES AND PACKAGES HERE,
## TO BE USED IN ALL SCRIPTS
## BY: JACOB LEVINE - jacoblevine@princeton.edu
##---------------------------------------------------------------

##---------------------------------------------------------------
## 0. User supplied variables
##---------------------------------------------------------------

## i. working directory, change this to the location of repo on your machine
root_dir <- "~/Documents/Science/severity_and_silviculture/"

## ii. mesonet weather API token, download from ___
api_token <- "8a99301c9d054a60b7f4711a7a6ac270"

## iii. list of fires: include spaces, all lowercase
fire_list <- c("dixie",
               "north complex", ## this is just the bear fire
               "sheep",
               "walker",
               "sugar")

## iv. coordinate reference system for analyses (CRS): supply as EPSG - don't recommend changing this
study_epsg <- 4326 ## unprojected WGS84 in lat/lon

## v.

##---------------------------------------------------------------
## 1. Commonly used libraries
##---------------------------------------------------------------
library(rgdal)
library(sf)
library(fasterize)
library(raster)
library(ggplot2)
library(chron)
library(rjson)
library(basemaps)
library(data.table)
library(foreach)
library(doParallel)
library(doSNOW)
library(progress)

##---------------------------------------------------------------
## 2. Set root working directory
##---------------------------------------------------------------

setwd(root_dir)
