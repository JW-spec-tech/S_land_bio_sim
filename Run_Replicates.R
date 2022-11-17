#### Load packages ####
library(arrow)
library(NLMR)
library(sf)
library(raster)
library(fasterize)
library(sspm)
library(rgeos)
library(tidyr)
library(readr)
library(dplyr)
library(mgcv)
library(readr)
library(foreach)
library(doParallel)
library(parallelly)
library(ranger)
library(tidyverse)
library(kableExtra)




#### Load the Functions ####
source(file = "Simulation_Depth_Temp_Coordinates_Biomass_function_V2.R")
source(file = "function_patches.R")
source(file = "Make_local.domain_local.arena.R")
source(file = "Resampling_test.R")
source(file = "Run_GAM_predict_replicates.R")




print(paste("Start of Sim generation @",Sys.time()))
reps= as.numeric(Sys.getenv('REPS')) # Number of Replicates
sims= as.numeric(Sys.getenv('SIMS')) # Number of Sims
size= as.numeric(Sys.getenv('SIZE')) # Size of the Landscape
seed= as.numeric(Sys.getenv('SEED')) # Starting seed
var = as.numeric(Sys.getenv('VAR'))  # Variation in biomass field --> higher variation = increased biomass variation

percent = as.numeric(Sys.getenv('PERCENT')) # Sets sampling percentage of the sampling of the entire dataset
# List all sims

# getwd()  <-- need to be in Data folder

setwd("Data_2022-11-17 14:07:14/")

# Get files names
f_list <- list.files()

# Load sim data

for (i in f_list) {
  setwd(i)
  # #### 1. Make replicates of each rep (landscape) ####
  Resample_write_PBfall(200)
  setwd("../")
}
parallel::stopCluster(cl = main.cluster)
