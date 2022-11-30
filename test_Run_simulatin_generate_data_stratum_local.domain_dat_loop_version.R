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
source(file = "Make_PB_fall.dat.R")
source(file = "Resampling_test.R")

# memory.limit(40000)

reps=20
sims=20
sizes=500
seed=round(runif(1,0,1000000))
var=1.5
percent=2

print(paste("Start of Sim generation @",Sys.time()))
# reps= as.numeric(Sys.getenv('REPS')) # Number of Replicates
# sims= as.numeric(Sys.getenv('SIMS')) # Number of Sims
# sizes= as.numeric(Sys.getenv('SIZE')) # Size of the Landscape
# seed= as.numeric(Sys.getenv('SEED')) # Starting seed
# var = as.numeric(Sys.getenv('VAR'))  # Variation in biomass field --> higher variation = increased biomass variation

# percent = as.numeric(Sys.getenv('PERCENT')) # Sets sampling percentage of the sampling of the entire dataset

#### create the cluster ####

  n.cores <- 10
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

  dir.create("Data_test")

  Sys.time()
#### Loop to run replicates of simulations in individual folders ####
  Sim_Loop <- foreach(
    rep = 1:reps,
    .packages = c('arrow','NLMR','sf','raster','fasterize','sspm','rgeos','tidyr','readr','erer',"mgcv", "plyr", "dplyr", "raster","landscapetools","devtools","openxlsx","sspm")
  ) %dopar% {
    #### 0. Set-up ####
  print(paste("Replicate #",rep,Sys.time()))
  seeds = seed + rep
  set.seed(seeds)
  cwd <- "Data_test"          # CURRENT dir
  # setwd("Data/")
  newdir <- paste("Run",rep,"Size",sizes,"seed",seeds,"nsim",sims,"Percent",percent,Sys.Date(),sep = "_")
  dir.create(paste0(cwd,"/",newdir))    # Create new directory
  cwd <- paste0(getwd(),"/",cwd,"/",newdir)
  # setwd(newdir) 
  write.table(as.data.frame(newdir),paste0(cwd,"/","seed"))
  
  #### 1. Run the sim ####
  
  results <- S_land_bio_sim(n=sims, size=sizes,variation = var) # higher variation = increased biomass variation
  
  # Save size of each strata
  patches=results$patches_list$patches
  patches=st_set_geometry(patches,NULL)
  write_parquet(patches,paste0(cwd,'/',"patches"))
  gc()
  # 
  #### 2. Write the ogmap files ####
  Make_patch_domain_arena_DAT(size=sizes,patches=results$patches_list$patches,the_stack=results$the_stack,percent=percent)
  gc()

  # #### 3. Slice data file ####
  Make_PB_fall.dat()
  gc()
  # 
  #### 4. Run STRAP Calculations ####
  STRAP()

  # setwd(cwd)
  gc()
  }
  
  replicates_gam <- function(S_year =1991) {
    
    start_year =S_year
    files <-  list.dirs("Data_test/", recursive = F, full.names = TRUE)
    years = length(files)
    size = 499
    
    for (d in files) {
      dir.create(paste0(d,"/Result"), recursive = T)
      trawl_data <- readr::read_table(paste0(d,"/PB_fall.dat"))
      print(d)
      trawl_data$year_f <- factor(trawl_data$year)
      #Filtering 
      # trawl_data <- trawl_data[trawl_data$year %in% 1991:2000,]
      
      dat_grid_x_y <- as.data.frame(expand_grid(long=seq(0.5,size,by=1),
                                                lat= seq(0.5,size,by=1)))
      
      dat_grid_year <- c(start_year:(start_year+years-1))
      
      #### testing ####
      # years_trawl <- unique(trawl_data$year)
      # 
      # trawl_data <- trawl_data[trawl_data$year %in% 1991:2000,]
      # 
      # diff <- trawl_data-trawl_data_1
      # sum(trawl_data$biomass)
      # sum(trawl_data_1$biomass)
      # nyears <- as.data.frame(table(trawl_data$year))
      # hist(trawl_data$year, breaks = 501)
      ####
      head(trawl_data)
      gc()
      trawl_data$year_f <- factor(trawl_data$year)
      # Sys.time()
      # simple_gam  <- bam((biomass/1000)~te(long, lat, year_f, bs= c("tp","re"), d = c(2,1)), family= "tw", data = trawl_data, method="REML")
      # Sys.time()
      #### Create models ####
      print("#### Create models ####")
      simple_gam <- list()
      
      #define number of data frames to split into
      split=5
      n_chunk <- length(dat_grid_year)/split
      #split data frame into groups per year
      split_data <- trawl_data %>% 
        group_split(year)
      
      # for (n_chunk in 1:n_chunk) {
      #   print(n_chunk)
      #   chunk <- split_data[(((n_chunk-1)*split)+1):(n_chunk*split)] %>% reduce(full_join)
      #   simple_gam[[n_chunk]]  <- bam((biomass/1000)~te(long, lat, year_f, bs= c("tp","re"), d = c(2,1)), family= "tw", data = chunk, method="REML")
      # }
      # 
      # Sys.time()
      
      # split_gams <- function(years_gam=trawl_data$year_f,dat_grid_year=dat_grid_year,split=10){
      
      #  #split data frame into n equal-sized data frames
      #   split_data <- trawl_data %>%
      #     group_split(year)
      #   for (n_chunk in 1:n_chunk) {
      #   print(n_chunk)
      #   chunk <- split_data[(((n_chunk-1)*split)+1):(n_chunk*split)] %>% reduce(full_join)
      #   simple_gam[[n_chunk]]  <- bam((biomass/1000)~te(long, lat, year_f, bs= c("tp","re"), d = c(2,1)), family= "tw", data = chunk, method="REML")
      # }
      # 
      
      
      # Sys.time()
      
      simple_gam <- foreach(
        n_chunk = 1:n_chunk,
        .packages = c('mgcv','dplyr','purrr')
      ) %dopar% {
        chunk <- split_data[(((n_chunk-1)*split)+1):(n_chunk*split)] %>% reduce(full_join)
        simple_gam[[n_chunk]]  <- bam((biomass/1000)~te(long, lat, year_f, bs= c("tp","re"), d = c(2,1)), family= "tw", data = chunk, method="REML")
        simple_gam[[n_chunk]]
      }
      # Sys.time()
      
      # save(simple_gam, file=paste0("resample_data/","rep",r,"/Result/gam_model.gam"))
      
      #### Get PRedicted biomass + CI function ####
      
      print("#### Get PRedicted biomass + CI function ####")
      Get_biomass_Ci_write <- function(fit,dat_per_year=dplyr::bind_cols(dat_grid_x_y,year_f=as.factor(year_f),year=year_f)){
        
        sims <- sspm:::produce_sims(fit, dat_per_year, 1000)
        sims <- exp(sims)
        
        sims_total <- apply(sims, MARGIN = 2, FUN = "sum")
        sims_point <- mean(sims_total)
        
        alpha = 0.05
        sims_CI <- quantile(sims_total, prob = c(alpha/2, 1-alpha/2))
        output <- data.frame(year=year_f,point_est = sims_point, lower = sims_CI[1], upper = sims_CI[2])
        dir.create(paste0(d,"/Result"))
        write_parquet(output,paste0(d,"/Result/model_", year_f))
        
      }
      
      #### 3. Run GAM Predictions ####
      print("#### 3. Run GAM Predictions ####")
      Sys.time()
      x <- foreach(
        year_f = dat_grid_year,
        .packages = c('mgcv','dplyr','sspm','arrow')
      ) %dopar% {
        # load(paste0("rep",r,"/Result/gam_model.gam"))
        for(gam in simple_gam){
          if(as.factor(year_f) %in% unique(gam$model$year_f)){
            Get_biomass_Ci_write(fit=gam)
            # return(list(year=year,years=unique(gam$model$year_f)))
          }
          gc()
        }
      }
      Sys.time()
      
      #### Stop Cluster ####
      
      
      
      
    }
    # parallel::stopCluster(cl = my.cluster)
  }
  
  replicates_gam()

print(paste("End of GAM predict @",Sys.time()))

stopCluster(my.cluster)

# END ################################################
