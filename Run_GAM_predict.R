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
library(dplyr)

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
# memory.limit(n.cores*4000)
memory.limit(300000)
Sys.time()

#### 2. Define grid size and number of years ####

start_year =1991
years = files <-  as.numeric(length(list.files('sim/', pattern = ".",all.files = FALSE, recursive = TRUE, full.names = TRUE)))
size = 500

trawl_data <- readr::read_table("PB_fall.dat")
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
simple_gam <<- list()

#define number of data frames to split into
split=10
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
  print(n_chunk)
  chunk <- split_data[(((n_chunk-1)*split)+1):(n_chunk*split)] %>% reduce(full_join)
  simple_gam[[n_chunk]]  <- bam((biomass/1000)~te(long, lat, year_f, bs= c("tp","re"), d = c(2,1)), family= "tw", data = chunk, method="REML")
  gc()
  return(simple_gam[[n_chunk]])
}
# Sys.time()


#### Get PRedicted biomass + CI function ####


Get_biomass_Ci_write <- function(fit,dat_per_year=dplyr::bind_cols(dat_grid_x_y,year_f=as.factor(year))){
  
  year=year
  
  sims <- sspm:::produce_sims(fit, dat_per_year, 1000)
  sims <- exp(sims)
  
  sims_total <- apply(sims, MARGIN = 2, FUN = "sum")
  sims_point <- mean(sims_total)
  
  alpha = 0.05
  sims_CI <- quantile(sims_total, prob = c(alpha/2, 1-alpha/2))
  output <- data.frame(year=year,point_est = sims_point, lower = sims_CI[1], upper = sims_CI[2])
  write_parquet(output,paste0('Result/model_', year))
  return(output)
}

#### 3. Run GAM Predictions ####
Sys.time()
x <- foreach(
  year = dat_grid_year,
  .packages = c('mgcv','dplyr','sspm','arrow')
) %dopar% {
    for(gam in simple_gam){
      if(as.factor(year) %in% unique(gam$model$year_f)){
        Get_biomass_Ci_write(fit=gam)
        # return(list(year=year,years=unique(gam$model$year_f)))
      }
      gc()
  }
}
Sys.time()

#### Stop Cluster ####


parallel::stopCluster(cl = my.cluster)
# 
# 
# dat_per_year <- dplyr::bind_cols(dat_grid_x_y,year=year) %>% 
#   dplyr::mutate(fit_simple_gam = predict.bam(simple_gam,type = "response", newdata = .)) 
# 
# sum(dat_per_year$fit_simple_gam)
# Predictions_summary[1]

read_parquet('Result/model')
