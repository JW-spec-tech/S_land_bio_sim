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



#### Load the Functions ####
source(file = "Simulation_Depth_Temp_Coordinates_Biomass_function_V2.R")
source(file = "function_patches.R")
source(file = "Make_local.domain_local.arena.R")
source(file = "Make_PB_fall.dat.R")


memory.limit(40000)

print(paste("Start of Sim generation @",Sys.time()))
reps= as.numeric(Sys.getenv('REPS')) # Number of Replicates
sims= as.numeric(Sys.getenv('SIMS')) # Number of Sims
size= as.numeric(Sys.getenv('SIZE')) # Size of the Landscape
seed= as.numeric(Sys.getenv('SEED')) # Starting seed
var = as.numeric(Sys.getenv('VAR'))  # Variation in biomass field --> higher variation = increased biomass variation

percent = as.numeric(Sys.getenv('PERCENT')) # Sets sampling percentage of the sampling of the entire dataset

#### Loop to run replicates of simulations in individual folders ####
for (rep in 1:reps) {
  print(paste("Replicate #",rep))
  seeds = seed - 1 + rep
  set.seed(seeds)
  cwd <- getwd()          # CURRENT dir
  setwd("Data/")
  newdir <- paste("Run",rep,"Size",size,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
  dir.create(newdir)     # Create new directory
  setwd(newdir) 
  write.table(as.data.frame(newdir),"seed")
  #### 1. Run the sim ####
  results <- S_land_bio_sim(sims,size,variation = var) # higher variation = increased biomass variation
  
  #### 2. Write the ogmap files ####
  Make_patch_domain_arena_DAT(size,patches=results$patches_list$patches,the_stack=results$the_stack,percent=percent)
  
  #### 3. Slice data file ####
  Make_PB_fall.dat()
  setwd(cwd)
  gc()
}

print(paste("End of Sim generation @",Sys.time()))


# ################################################