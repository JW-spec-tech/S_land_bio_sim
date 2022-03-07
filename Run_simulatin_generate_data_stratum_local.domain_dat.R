#### Load packages ####
library(arrow)
library(NLMR)
library(sf)
library(dplyr)
library(raster)
library(fasterize)
library(sspm)
library(rgeos)
library(tidyr)
library(readr)


#### Load the Functions ####
source(file = "Simulation_Depth_Temp_Coordinates_Biomass_function.R")
source(file = "function_patches.R")
source(file = "Make_local.domain_local.arena.R")
source(file = "Make_PB_fall.dat.R")

#### 1. Run the sim ####
results <- S_land_bio_sim(5,100,100,50,5,0.2,0.5)

#### 2. Write the ogmap files ####
Make_patch_domain_arena_DAT(patches=results$patches_list$patches,the_stack=results$the_stack)

#### 3. Slice data file ####
Make_PB_fall.dat()



################################################
