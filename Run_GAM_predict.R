#### Load packages ####
library(mgcv)
library(readr)
library(sspm)
library(tidyr)
library(arrow)
library(foreach)
library(doParallel)
library(parallelly)
library(ranger)
library(tidyverse)
library(kableExtra)

#### 1. Create and Start Cluster ####


#create the cluster
n.cores <- parallelly::availableCores() - 3
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
memory.limit(n.cores*4000)
Sys.time()

#### 2. Define grid size and number of years ####

start_year =1991
years = files <-  as.numeric(length(list.files('sim/', pattern = ".",all.files = FALSE, recursive = TRUE, full.names = TRUE)))
size = 500

trawl_data <- readr::read_table("PB_fall.dat")
head(trawl_data)

simple_gam <- bam((biomass/1000) ~
                    te(long, lat, year, bs= c("tp", "tp"), d = c(2,1)),
                  data = trawl_data, method="REML", family = "tw")
dat_grid_x_y <- as.data.frame(expand_grid(long=seq(0.5,size,by=1),
                                          lat= seq(0.5,size,by=1)))

dat_grid_year <- c(start_year:(start_year+years-1))


#### Get PRedicted biomass + CI function ####

Get_biomass_Ci_write <- function(fit=simple_gam,dat_per_year=dplyr::bind_cols(dat_grid_x_y,year=year)){
  year=year
  sims <- sspm:::produce_sims(fit, dat_per_year, 500)
  sims <- exp(sims)
  
  sims_total <- apply(sims, MARGIN = 2, FUN = "sum")
  sims_point <- mean(sims_total)
  
  alpha = 0.05
  sims_CI <- quantile(sims_total, prob = c(alpha/2, 1-alpha/2))
  output <- data.frame(year=year,point_est = sims_point, lower = sims_CI[1], upper = sims_CI[2])
  write_parquet(output,paste0('Result/model_', year))
}

#### 3. Run GAM Predictions ####
Sys.time()
x <- foreach(
  year = dat_grid_year,
  .packages = c('mgcv','dplyr','sspm','arrow')
) %dopar% {
  Get_biomass_Ci_write()
  gc()
}
Sys.time()

#### Stop Cluster ####

parallel::stopCluster(cl = my.cluster)