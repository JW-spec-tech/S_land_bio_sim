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


#### 3. Run GAM Predictions ####
Sys.time()
x <- foreach(
  year = dat_grid_year,
  .packages = c('mgcv','dplyr','sspm','arrow')
) %dopar% {
  dat_per_year <- dplyr::bind_cols(dat_grid_x_y,year=year) %>% 
    dplyr::mutate(fit_simple_gam = predict.bam(simple_gam,type = "response", newdata = .)) %>%
    dplyr::bind_cols(predict_intervals(object_fit = simple_gam, new_data = ., PI= F, n=500))
  # write.table(dat_per_year, file = paste0('Result/model_', year,'.csv'), sep = ",",row.names = F)
  write_parquet(dat_per_year,paste0('Result/model_', year))
  gc()
}
Sys.time()

#### Stop Cluster ####

parallel::stopCluster(cl = my.cluster)
