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
library(doParallel)
library(erer)



#### Load the Functions ####
source(file = "Simulation_Depth_Temp_Coordinates_Biomass_function_V2.R")
source(file = "function_patches.R")
source(file = "Make_local.domain_local.arena.R")
source(file = "Resampling_test.R")

# memory.limit(40000)

reps=2
sims=10
size=500
seed=round(runif(1,0,1000000))
var=1.5
percent=2

print(paste("Start of Sim generation @",Sys.time()))
# reps= as.numeric(Sys.getenv('REPS')) # Number of Replicates
# sims= as.numeric(Sys.getenv('SIMS')) # Number of Sims
# size= as.numeric(Sys.getenv('SIZE')) # Size of the Landscape
# seed= as.numeric(Sys.getenv('SEED')) # Starting seed
# var = as.numeric(Sys.getenv('VAR'))  # Variation in biomass field --> higher variation = increased biomass variation

# percent = as.numeric(Sys.getenv('PERCENT')) # Sets sampling percentage of the sampling of the entire dataset

#### create the cluster ####

  n.cores <- 5
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  #check cluster definition (optional)
  print(my.cluster)
  
  
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  
  #how many workers are available? (optional)
  foreach::getDoParWorkers()

Sys.time()
#### Loop to run replicates of simulations in individual folders ####
  Sim_Loop <- foreach(
    rep = 1:reps,
    .packages = c('arrow','NLMR','sf','raster','fasterize','sspm','rgeos','tidyr','readr','dplyr','erer')
  ) %dopar% {
    #### 0. Set-up ####
  print(paste("Replicate #",rep,Sys.time()))
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
  
  # Save size of each strata
  patches=results$patches_list$patches
  patches=st_set_geometry(patches,NULL)
  write_parquet(patches,"patches")
  gc()
  
  #### 2. Write the ogmap files ####
  Make_patch_domain_arena_DAT(size,patches=results$patches_list$patches,the_stack=results$the_stack,percent=percent)
  gc()
  
  #### 3. Slice data file ####
  Make_PB_fall.dat()
  gc()
  
  #### 4. Run STRAP calculations ####
  STRAP()
  
  setwd(cwd)
  gc()
}

print(paste("End of Sim generation @",Sys.time()))

stopCluster(my.cluster)

# ################################################
