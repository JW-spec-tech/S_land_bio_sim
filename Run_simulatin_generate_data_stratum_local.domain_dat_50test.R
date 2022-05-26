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

#### 1. Run the sim ####
memory.limit(40000)

Sys.time()
size=500
results <- S_land_bio_sim(500,size) # highr variation = increased biomass variation
Sys.time()
#### 2. Write the ogmap files ####
Make_patch_domain_arena_DAT(size,patches=results$patches_list$patches,the_stack=results$the_stack)

#### 3. Slice data file ####
Make_PB_fall.dat()



# ################################################
