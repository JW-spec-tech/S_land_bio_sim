

print(paste("Start of Surface analysis @",Sys.time()))
reps= as.numeric(Sys.getenv('REPS')) # Number of Replicates
sims= as.numeric(Sys.getenv('SIMS')) # Number of Sims
sizes= as.numeric(Sys.getenv('SIZE')) # Size of the Landscape
seed= as.numeric(Sys.getenv('SEED')) # Starting seed
var = as.numeric(Sys.getenv('VAR'))  # Variation in biomass field --> higher variation = increased biomass variation

percent = as.numeric(Sys.getenv('PERCENT')) # Sets sampling percentage of the sampling of the entire dataset

#### create the cluster ####

n.cores <- as.numeric(Sys.getenv('OMP_NUM_THREADS'))
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


Surface_analysis <- function(S_year =1991, c_wd=cwd) {
  
  
  #### 1. Find all sims and files ####
  start_year =S_year
  files <-  list.dirs("Data_test/", recursive = F, full.names = TRUE)
  years = length(files)
  size = 499
  
  
  #### Loop to run Surface analysis on each sim individually ####
  
    Surface <- foreach(
      d = files,
      .packages = c('mgcv','dplyr','purrr','raster','geodiv','arrow')
    ) %dopar% {
      
      #### 1. Find all sims and files ####
      dir.create(paste0(d,"/Result_Surface"), recursive = T)
      raster_sim <- arrow::read_parquet(paste0(d,"/sim/sim1"))
      raster_sim <- data.frame(raster_sim$coord.x,raster_sim$coord.y,raster_sim$biomass)
      raster_sim <- rasterFromXYZ(raster_sim)
      
      #### 2. Surface analysis ####
        # Roughness Calculation
        
        RC <- geodiv::sa(raster_sim)
        
        # Surface bearing index (peaks)
        
        SBI <- geodiv::sbi(raster_sim)
        
      #### 3. Write files ####
        surface_C <- data.frame(RC=RC,SBI=SBI)
        write.table(surface_C, paste0(d,"/Result_Surface/Surface.charaterization"),row.names=F)
        
      #### 4. GC ####
        gc()

    }

  }


parallel::stopCluster(cl = my.cluster)

print(paste("End of Surface analysis @",Sys.time()))