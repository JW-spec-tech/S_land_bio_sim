
library(parallel)
library(parallelly)
library(doParallel)
library(mgcv)
library(dplyr)
library(purrr)
library(raster)
library(geodiv)
library(arrow)


print(paste("Start of Surface analysis @",Sys.time()))
reps= as.numeric(Sys.getenv('REPS')) # Number of Replicates
sims= as.numeric(Sys.getenv('SIMS')) # Number of Sims
sizes= as.numeric(Sys.getenv('SIZE')) # Size of the Landscape
seed= as.numeric(Sys.getenv('SEED')) # Starting seed
var = as.numeric(Sys.getenv('VAR'))  # Variation in biomass field --> higher variation = increased biomass variation

percent = as.numeric(Sys.getenv('PERCENT')) # Sets sampling percentage of the sampling of the entire dataset

# #### create the cluster ####
# 
# n.cores <- as.numeric(Sys.getenv('OMP_NUM_THREADS'))
# my.cluster <- parallel::makeCluster(
#   n.cores, 
#   type = "PSOCK"
# )
# 
# #check cluster definition (optional)
# print(my.cluster)
# 
# 
# #register it to be used by %dopar%
# doParallel::registerDoParallel(cl = my.cluster)
# 
# #check if it is registered (optional)
# foreach::getDoParRegistered()
# 
# #how many workers are available? (optional)
# foreach::getDoParWorkers()


Sys.time()


Surface_analysis <- function() {
  
  
  for (cwd in list.dirs(full.names = T,recursive = F)) {
    
    print(cwd)
    setwd(cwd)
    cwd=getwd()
    setwd(cwd)
    
    #### 1. Find all sims and files ####
    dir.create(paste0(cwd,"/Result_Surface"), recursive = T)
    raster_sim <- arrow::read_parquet(paste0(cwd,"/sim/sim1"))
    raster_sim <- data.frame(raster_sim$coord.x,raster_sim$coord.y,raster_sim$biomass)
    raster_sim <- rasterFromXYZ(raster_sim)
    
    #### 2. Surface analysis ####
    # Roughness Calculation
    
    RC <- geodiv::sa(raster_sim)
    
    # Surface bearing index (peaks)
    
    SBI <- geodiv::sbi(raster_sim)
    
    #### 3. Write files ####
    surface_C <- data.frame(RC=RC,SBI=SBI)
    write.table(surface_C, paste0(cwd,"/Result_Surface/Surface.charaterization"),row.names=F)
    
    #### 4. GC ####
    gc()
  # 
  # #### 1. Find all sims and files ####
  # start_year =S_year
  # files <-  list.dirs(recursive = F, full.names = TRUE)
  # years = length(files)
  # size = 499
  # 
  # 
  # #### Loop to run Surface analysis on each sim individually ####
  # 
  #   Surface <- foreach(
  #     d = files,
  #     .packages = c('mgcv','dplyr','purrr','raster','geodiv','arrow')
  #   ) %dopar% {
  #     
  #     #### 1. Find all sims and files ####
  #     dir.create(paste0(d,"/Result_Surface"), recursive = T)
  #     raster_sim <- arrow::read_parquet(paste0(d,"/sim/sim1"))
  #     raster_sim <- data.frame(raster_sim$coord.x,raster_sim$coord.y,raster_sim$biomass)
  #     raster_sim <- rasterFromXYZ(raster_sim)
  #     
  #     #### 2. Surface analysis ####
  #       # Roughness Calculation
  #       
  #       RC <- geodiv::sa(raster_sim)
  #       
  #       # Surface bearing index (peaks)
  #       
  #       SBI <- geodiv::sbi(raster_sim)
  #       
  #     #### 3. Write files ####
  #       surface_C <- data.frame(RC=RC,SBI=SBI)
  #       write.table(surface_C, paste0(d,"/Result_Surface/Surface.charaterization"),row.names=F)
  #       
  #     #### 4. GC ####
  #       gc()
    
    setwd("../")
        
    }
  
  Surface_roughness_list <- list()
  counter = 1
  
  files <- list.dirs(full.names = T,recursive = F)
  
  for (file in files) {
    file <- substring(file, 2)
    file = paste0(getwd(),file)
    print(file)
    Surface_roughness_list[[counter]] <- readr::read_table(paste0(file,"/Result_Surface/Surface.charaterization"))
    counter = counter + 1
  }
  
  Surface_roughness_df <-  do.call(rbind.data.frame, Surface_roughness_list)
  
  Surface_roughness_df_mean <- colMeans(Surface_roughness_df)
  
  write.table(Surface_roughness_df_mean, "Surface_roughness_mean")
  
  setwd("../")

}

# Get files names
f_list <- paste0(getwd(),"/",list.dirs(path = "exp", full.names = TRUE, recursive = F))

for (i in f_list) {
  print(i)
  setwd(i)
  Surface_analysis()
  setwd("~/Git projects/S_land_bio_sim")
}

parallel::stopCluster(cl = my.cluster)

print(paste("End of Surface analysis @",Sys.time()))