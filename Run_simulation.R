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
source(file = "Strap_calculations.R")
source(file = "Make_PB_fall.dat.R")
# source(file = "Resampling_test.R")
# source(file = "Run_GAM_predict_replicates.R")




print(paste("Start of Sim generation @",Sys.time(),"test"))
reps= as.numeric(Sys.getenv('REPS')) # Number of Replicates
sims= as.numeric(Sys.getenv('SIMS')) # Number of Sims
size= as.numeric(Sys.getenv('SIZE')) # Size of the Landscape
seed= as.numeric(Sys.getenv('SEED')) # Starting seed
var = as.numeric(Sys.getenv('VAR'))  # Variation in biomass field --> higher variation = increased biomass variation
percent = as.numeric(Sys.getenv('PERCENT')) # Sets sampling percentage of the sampling of the entire dataset

memory.limit(20000)

print(paste("Start of Sim generation @",Sys.time(),"test"))
reps= 5
sims= 20
size= 500
seed= sample(0:100000,1)
var = 1
percent = 2




#### 1. Create and Start Cluster ####
print(paste("Start of Sim generation @",Sys.time(),"test"))

#create the cluster
# n.cores <- parallelly::availableCores()/2   
# For windows
n.cores <- 8
main.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(main.cluster)


#register it to be used by %dopar%
doParallel::registerDoParallel(cl = main.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

dir_now <- paste0("~/Git projects/S_land_bio_sim/exp/Experiment_Variation_",var)

dir.create(dir_now)
setwd(paste0(dir_now,"/"))
# dir.create(dir_now)
#   setwd("Data_test/")
  print(dir_now)
print(paste("Start of Sim generation @",Sys.time(),"TEST"))

cwd <- getwd()
cwd

foreach(
  rep = 1:reps,
  .packages = c('mgcv','dplyr','purrr','NLMR','arrow','sspm','raster','foreach','doParallel','parallelly','readr','fasterize')
) %dopar% {
  setwd(paste0(dir_now,"/"))
  print(paste("Replicate #",rep))
  seeds = seed - 10 + rep
  set.seed(seeds)
  newdir <- paste("Run",rep,"Size",size,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
  dir.create(newdir)     # Create new directory
  setwd(paste0(newdir,"/"))
  cwd=getwd()
  cwd
  write.table(as.data.frame(newdir),"seed")
  #### 1. Run the sim ####
  results <- S_land_bio_sim(n=sims,size=size,variation = var) # higher variation = increased biomass variation
  
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
  
  gc()
  
}

#### Kill the cluster
parallel::stopCluster(cl = main.cluster)

setwd("~/Git projects/S_land_bio_sim")



print(paste("Start of Sim generation @",Sys.time(),"test"))
reps= 5
sims= 20
size= 500
seed= sample(0:100000,1)
var = 1
percent = 2




#### 1. Create and Start Cluster ####
print(paste("Start of Sim generation @",Sys.time(),"test"))

#create the cluster
# n.cores <- parallelly::availableCores()/2   
# For windows
n.cores <- 8
main.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(main.cluster)


#register it to be used by %dopar%
doParallel::registerDoParallel(cl = main.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

dir_now <- paste0("~/Git projects/S_land_bio_sim/exp/Experiment_Variation_",var)

dir.create(dir_now)
setwd(paste0(dir_now,"/"))
# dir.create(dir_now)
#   setwd("Data_test/")
print(dir_now)
print(paste("Start of Sim generation @",Sys.time(),"TEST"))

cwd <- getwd()
cwd

foreach(
  rep = 1:reps,
  .packages = c('mgcv','dplyr','purrr','NLMR','arrow','sspm','raster','foreach','doParallel','parallelly','readr','fasterize')
) %dopar% {
  setwd(paste0(dir_now,"/"))
  print(paste("Replicate #",rep))
  seeds = seed - 10 + rep
  set.seed(seeds)
  newdir <- paste("Run",rep,"Size",size,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
  dir.create(newdir)     # Create new directory
  setwd(paste0(newdir,"/"))
  cwd=getwd()
  cwd
  write.table(as.data.frame(newdir),"seed")
  #### 1. Run the sim ####
  results <- S_land_bio_sim(n=sims,size=size,variation = var) # higher variation = increased biomass variation
  
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
  
  gc()
  
}

#### Kill the cluster
parallel::stopCluster(cl = main.cluster)

setwd("~/Git projects/S_land_bio_sim")



print(paste("Start of Sim generation @",Sys.time(),"test"))
reps= 5
sims= 20
size= 500
seed= sample(0:100000,1)
var = 5
percent = 2




#### 1. Create and Start Cluster ####
print(paste("Start of Sim generation @",Sys.time(),"test"))

#create the cluster
# n.cores <- parallelly::availableCores()/2   
# For windows
n.cores <- 8
main.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(main.cluster)


#register it to be used by %dopar%
doParallel::registerDoParallel(cl = main.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

dir_now <- paste0("~/Git projects/S_land_bio_sim/exp/Experiment_Variation_",var)

dir.create(dir_now)
setwd(paste0(dir_now,"/"))
# dir.create(dir_now)
#   setwd("Data_test/")
print(dir_now)
print(paste("Start of Sim generation @",Sys.time(),"TEST"))

cwd <- getwd()
cwd

foreach(
  rep = 1:reps,
  .packages = c('mgcv','dplyr','purrr','NLMR','arrow','sspm','raster','foreach','doParallel','parallelly','readr','fasterize')
) %dopar% {
  setwd(paste0(dir_now,"/"))
  print(paste("Replicate #",rep))
  seeds = seed - 10 + rep
  set.seed(seeds)
  newdir <- paste("Run",rep,"Size",size,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
  dir.create(newdir)     # Create new directory
  setwd(paste0(newdir,"/"))
  cwd=getwd()
  cwd
  write.table(as.data.frame(newdir),"seed")
  #### 1. Run the sim ####
  results <- S_land_bio_sim(n=sims,size=size,variation = var) # higher variation = increased biomass variation
  
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
  
  gc()
  
}

#### Kill the cluster
parallel::stopCluster(cl = main.cluster)

setwd("~/Git projects/S_land_bio_sim")




print(paste("Start of Sim generation @",Sys.time(),"test"))
reps= 5
sims= 20
size= 500
seed= sample(0:100000,1)
var = 10
percent = 2




#### 1. Create and Start Cluster ####
print(paste("Start of Sim generation @",Sys.time(),"test"))

#create the cluster
# n.cores <- parallelly::availableCores()/2   
# For windows
n.cores <- 8
main.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(main.cluster)


#register it to be used by %dopar%
doParallel::registerDoParallel(cl = main.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

dir_now <- paste0("~/Git projects/S_land_bio_sim/exp/Experiment_Variation_",var)

dir.create(dir_now)
setwd(paste0(dir_now,"/"))
# dir.create(dir_now)
#   setwd("Data_test/")
print(dir_now)
print(paste("Start of Sim generation @",Sys.time(),"TEST"))

cwd <- getwd()
cwd

foreach(
  rep = 1:reps,
  .packages = c('mgcv','dplyr','purrr','NLMR','arrow','sspm','raster','foreach','doParallel','parallelly','readr','fasterize')
) %dopar% {
  setwd(paste0(dir_now,"/"))
  print(paste("Replicate #",rep))
  seeds = seed - 10 + rep
  set.seed(seeds)
  newdir <- paste("Run",rep,"Size",size,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
  dir.create(newdir)     # Create new directory
  setwd(paste0(newdir,"/"))
  cwd=getwd()
  cwd
  write.table(as.data.frame(newdir),"seed")
  #### 1. Run the sim ####
  results <- S_land_bio_sim(n=sims,size=size,variation = var) # higher variation = increased biomass variation
  
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
  
  gc()
  
}

#### Kill the cluster
parallel::stopCluster(cl = main.cluster)

setwd("~/Git projects/S_land_bio_sim")



print(paste("Start of Sim generation @",Sys.time(),"test"))
reps= 5
sims= 20
size= 500
seed= sample(0:100000,1)
var = 25
percent = 2




#### 1. Create and Start Cluster ####
print(paste("Start of Sim generation @",Sys.time(),"test"))

#create the cluster
# n.cores <- parallelly::availableCores()/2   
# For windows
n.cores <- 8
main.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(main.cluster)


#register it to be used by %dopar%
doParallel::registerDoParallel(cl = main.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

dir_now <- paste0("~/Git projects/S_land_bio_sim/exp/Experiment_Variation_",var)

dir.create(dir_now)
setwd(paste0(dir_now,"/"))
# dir.create(dir_now)
#   setwd("Data_test/")
print(dir_now)
print(paste("Start of Sim generation @",Sys.time(),"TEST"))

cwd <- getwd()
cwd

foreach(
  rep = 1:reps,
  .packages = c('mgcv','dplyr','purrr','NLMR','arrow','sspm','raster','foreach','doParallel','parallelly','readr','fasterize')
) %dopar% {
  setwd(paste0(dir_now,"/"))
  print(paste("Replicate #",rep))
  seeds = seed - 10 + rep
  set.seed(seeds)
  newdir <- paste("Run",rep,"Size",size,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
  dir.create(newdir)     # Create new directory
  setwd(paste0(newdir,"/"))
  cwd=getwd()
  cwd
  write.table(as.data.frame(newdir),"seed")
  #### 1. Run the sim ####
  results <- S_land_bio_sim(n=sims,size=size,variation = var) # higher variation = increased biomass variation
  
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
  
  gc()
  
}

#### Kill the cluster
parallel::stopCluster(cl = main.cluster)

setwd("~/Git projects/S_land_bio_sim")


print(paste("Start of Sim generation @",Sys.time(),"test"))
reps= 5
sims= 20
size= 500
seed= sample(0:100000,1)
var = 50
percent = 2




#### 1. Create and Start Cluster ####
print(paste("Start of Sim generation @",Sys.time(),"test"))

#create the cluster
# n.cores <- parallelly::availableCores()/2   
# For windows
n.cores <- 8
main.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(main.cluster)


#register it to be used by %dopar%
doParallel::registerDoParallel(cl = main.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

dir_now <- paste0("~/Git projects/S_land_bio_sim/exp/Experiment_Variation_",var)

dir.create(dir_now)
setwd(paste0(dir_now,"/"))
# dir.create(dir_now)
#   setwd("Data_test/")
print(dir_now)
print(paste("Start of Sim generation @",Sys.time(),"TEST"))

cwd <- getwd()
cwd

foreach(
  rep = 1:reps,
  .packages = c('mgcv','dplyr','purrr','NLMR','arrow','sspm','raster','foreach','doParallel','parallelly','readr','fasterize')
) %dopar% {
  setwd(paste0(dir_now,"/"))
  print(paste("Replicate #",rep))
  seeds = seed - 10 + rep
  set.seed(seeds)
  newdir <- paste("Run",rep,"Size",size,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
  dir.create(newdir)     # Create new directory
  setwd(paste0(newdir,"/"))
  cwd=getwd()
  cwd
  write.table(as.data.frame(newdir),"seed")
  #### 1. Run the sim ####
  results <- S_land_bio_sim(n=sims,size=size,variation = var) # higher variation = increased biomass variation
  
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
  
  gc()
  
}

#### Kill the cluster
parallel::stopCluster(cl = main.cluster)

setwd("~/Git projects/S_land_bio_sim")
#### Loop to run replicates of simulations in individual folders ####
# for (rep in 1:reps) {
#   print(paste("Replicate #",rep))
#   seeds = seed - 1 + rep
#   set.seed(seeds)
#   cwd <- getwd()          # CURRENT dir
#   setwd("Data/")
#   newdir <- paste("Run",rep,"Size",size,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
#   dir.create(newdir)     # Create new directory
#   setwd(newdir) 
#   write.table(as.data.frame(newdir),"seed")
#   #### 1. Run the sim ####
#   results <- S_land_bio_sim(sims,size,variation = var) # higher variation = increased biomass variation
#   
#   # Save size of each strata
#   patches=results$patches_list$patches
#   patches=st_set_geometry(patches,NULL)
#   write_parquet(patches,"patches")
#   
#   #### 2. Write the ogmap files ####
#   Make_patch_domain_arena_DAT(size,patches=results$patches_list$patches,the_stack=results$the_stack,percent=percent)
#   
#   # #### 3. Make replicates of each rep (landscape) ####
#   Resample_sims_random(200)
#   # 
#   # #### 4. Run the GAM ####
#   replicates_gam()
#   # 
#   # #### 5. Run Comparison + Print graphs ####
#   # Compare_Graph()
#   
#   #### 6. Return to Original WD ####
#   setwd(cwd)
#   gc()
# }


print(paste("Start of Sim generation @",Sys.time(),"test"))
reps= as.numeric(Sys.getenv('REPS')) # Number of Replicates
sims= as.numeric(Sys.getenv('SIMS')) # Number of Sims
size= as.numeric(Sys.getenv('SIZE')) # Size of the Landscape
seed= as.numeric(Sys.getenv('SEED')) # Starting seed
var = as.numeric(Sys.getenv('VAR'))  # Variation in biomass field --> higher variation = increased biomass variation

percent = as.numeric(Sys.getenv('PERCENT')) # Sets sampling percentage of the sampling of the entire dataset


#### 1. Create and Start Cluster ####
print(paste("Start of Sim generation @",Sys.time(),"test"))

#create the cluster
# n.cores <- parallelly::availableCores()/2   
# For windows
n.cores <- as.numeric(Sys.getenv('OMP_NUM_THREADS'))
main.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(main.cluster)


#register it to be used by %dopar%
doParallel::registerDoParallel(cl = main.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()


reps= 20
sims= 20
size= 500
seed= sample(1:1000000, 1, replace=TRUE)
var = 1

dir_now <- paste0("Data_test_",var)

dir.create(dir_now)
setwd(paste0(dir_now,"/"))
# dir.create(dir_now)
#   setwd("Data_test/")
  print(dir_now)
print(paste("Start of Sim generation @",Sys.time(),"TEST"))

cwd <- getwd()
cwd

foreach(
  rep = 1:reps,
  .packages = c('mgcv','dplyr','purrr','NLMR','arrow','sspm','raster','foreach','doParallel','parallelly','readr','fasterize')
) %dopar% {
  cwd <- getwd()
  print(paste("Replicate #",rep))
  seeds = seed - 10 + rep
  set.seed(seeds)
  newdir <- paste("Run",rep,"Size",size,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
  dir.create(newdir)     # Create new directory
  setwd(paste0(newdir,"/")) 
  write.table(as.data.frame(newdir),"seed")
  #### 1. Run the sim ####
  results <- S_land_bio_sim(n= sims,size = size,variation = var) # higher variation = increased biomass variation
  
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

setwd("~/Git projects/S_land_bio_sim")


reps= 20
sims= 20
size= 500
seed= sample(1:1000000, 1, replace=TRUE)
var = 5

dir_now <- paste0("Data_test_",var)

dir.create(dir_now)
setwd(paste0(dir_now,"/"))
# dir.create(dir_now)
#   setwd("Data_test/")
print(dir_now)
print(paste("Start of Sim generation @",Sys.time(),"TEST"))

cwd <- getwd()
cwd

foreach(
  rep = 1:reps,
  .packages = c('mgcv','dplyr','purrr','NLMR','arrow','sspm','raster','foreach','doParallel','parallelly','readr','fasterize')
) %dopar% {
  cwd <- getwd()
  print(paste("Replicate #",rep))
  seeds = seed - 10 + rep
  set.seed(seeds)
  newdir <- paste("Run",rep,"Size",size,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
  dir.create(newdir)     # Create new directory
  setwd(paste0(newdir,"/")) 
  write.table(as.data.frame(newdir),"seed")
  #### 1. Run the sim ####
  results <- S_land_bio_sim(n= sims,size = size,variation = var) # higher variation = increased biomass variation
  
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

setwd("~/Git projects/S_land_bio_sim")


reps= 20
sims= 20
size= 500
seed= sample(1:1000000, 1, replace=TRUE)
var = 10

dir_now <- paste0("Data_test_",var)

dir.create(dir_now)
setwd(paste0(dir_now,"/"))
# dir.create(dir_now)
#   setwd("Data_test/")
print(dir_now)
print(paste("Start of Sim generation @",Sys.time(),"TEST"))

cwd <- getwd()
cwd

foreach(
  rep = 1:reps,
  .packages = c('mgcv','dplyr','purrr','NLMR','arrow','sspm','raster','foreach','doParallel','parallelly','readr','fasterize')
) %dopar% {
  cwd <- getwd()
  print(paste("Replicate #",rep))
  seeds = seed - 10 + rep
  set.seed(seeds)
  newdir <- paste("Run",rep,"Size",size,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
  dir.create(newdir)     # Create new directory
  setwd(paste0(newdir,"/")) 
  write.table(as.data.frame(newdir),"seed")
  #### 1. Run the sim ####
  results <- S_land_bio_sim(n= sims,size = size,variation = var) # higher variation = increased biomass variation
  
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


setwd("~/Git projects/S_land_bio_sim")


reps= 20
sims= 20
size= 500
seed= sample(1:1000000, 1, replace=TRUE)
var = 25

dir_now <- paste0("Data_test_",var)

dir.create(dir_now)
setwd(paste0(dir_now,"/"))
# dir.create(dir_now)
#   setwd("Data_test/")
print(dir_now)
print(paste("Start of Sim generation @",Sys.time(),"TEST"))

cwd <- getwd()
cwd

foreach(
  rep = 1:reps,
  .packages = c('mgcv','dplyr','purrr','NLMR','arrow','sspm','raster','foreach','doParallel','parallelly','readr','fasterize')
) %dopar% {
  cwd <- getwd()
  print(paste("Replicate #",rep))
  seeds = seed - 10 + rep
  set.seed(seeds)
  newdir <- paste("Run",rep,"Size",size,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
  dir.create(newdir)     # Create new directory
  setwd(paste0(newdir,"/")) 
  write.table(as.data.frame(newdir),"seed")
  #### 1. Run the sim ####
  results <- S_land_bio_sim(n= sims,size = size,variation = var) # higher variation = increased biomass variation
  
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

setwd("~/Git projects/S_land_bio_sim")


reps= 20
sims= 20
size= 500
seed= sample(1:1000000, 1, replace=TRUE)
var = 50

dir_now <- paste0("Data_test_",var)

dir.create(dir_now)
setwd(paste0(dir_now,"/"))
# dir.create(dir_now)
#   setwd("Data_test/")
print(dir_now)
print(paste("Start of Sim generation @",Sys.time(),"TEST"))

cwd <- getwd()
cwd

foreach(
  rep = 1:reps,
  .packages = c('mgcv','dplyr','purrr','NLMR','arrow','sspm','raster','foreach','doParallel','parallelly','readr','fasterize')
) %dopar% {
  cwd <- getwd()
  print(paste("Replicate #",rep))
  seeds = seed - 10 + rep
  set.seed(seeds)
  newdir <- paste("Run",rep,"Size",size,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
  dir.create(newdir)     # Create new directory
  setwd(paste0(newdir,"/")) 
  write.table(as.data.frame(newdir),"seed")
  #### 1. Run the sim ####
  results <- S_land_bio_sim(n= sims,size = size,variation = var) # higher variation = increased biomass variation
  
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


#################################################