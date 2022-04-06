library(mgcv)
library(tidyr)
library(sspm)
# Gam_prediction <-function(path="PB.fall.dat",size=500)

# 1. Load the biomass sample data & generate new prediction grid
# x=size
# y=size
sample <- readr::read_table("PB_fall.dat")
# grid <- expand_grid(x=seq(0,1,length=500), y= seq(0,1,length=500))


# 2. Run the models
simple_gam <- bam(biomass ~ te(long, lat, year, bs= c("tp", "tp"), d = c(2,1)), data= sample, method="REML", family = "tw")

spatialonly_gam <- gam(biomass ~s(long, lat, bs="tp"), data= sample, method="REML", family = "tw")


# 3. Run the predictions
dat_grid <- as.data.frame(expand_grid(long=seq(0,1,length=50), lat= seq(0,1,length=50),year=c(1991:2000)))%>%
  dplyr::mutate(fit_simple_gam  = predict.gam(simple_gam,type = "response",newdata = .),
                fit_spatialonly_gam = predict.gam(spatialonly_gam,type = "response",newdata = .))




library(mgcv)
library(readr)
library(sspm)
library(tidyr)
library(arrow)
trawl_data <- readr::read_table("PB_fall.dat")
head(trawl_data)

simple_gam <- bam(biomass ~
                    te(long, lat, year, bs= c("tp", "tp"), d = c(2,1)),
                  data = trawl_data, method="REML", family = "tw")


# spatialonly_gam <- gam(biomass ~s(long, lat, bs="tp"),
#                        data = trawl_data, method="REML", family = "tw")

dat_grid <- as.data.frame(expand_grid(long=seq(0,5,by=0.5),
                                      lat= seq(0,5,by=0.5),
                                      year=c(1991:2040)))


# predict_intervals <- sspm:::predict_productivity_intervals
# produce_sims <- function (fit, new_data, n = 100) 
# {
#   checkmate::assert_class(fit, "gam")
#   coefs <- stats::coef(fit)
#   lp <- predict(fit, newdata = new_data, type = "lpmatrix")
#   vcv <- stats::vcov(fit)
#   coefs_sim <- t(rmvn(n = n, coefs, vcv))
#   sims <- lp %*% coefs_sim
#   return(sims)
# }
# 
# sspm:::pred
# find_quantiles <- sspm:::find_quantiles
# confidence_interval <- sspm:::confidence_interval
# prediction_interval <- sspm:::prediction_interval

# predict(simple_gam, newdata = dat_grid, type = "lpmatrix")
sspm::predict_intervals()
memory.limit(280000)
gc()
library(ff)
# require(parallel)  
# nc <- 2   ## cluster size, set for example portability
# if (detectCores()>1) { ## no point otherwise
#   cl <- makeCluster(nc) 
#   ## could also use makeForkCluster, but read warnings first!
# } else cl <- NULL
library(data.table)
gc()
Sys.time()

dat_grid_pred_1 <- dat_grid %>%
  dplyr::mutate(fit_simple_gam = predict.bam(simple_gam,type = "response",
                                             newdata = .)) %>%
  dplyr::bind_cols(predict_intervals(object_fit = simple_gam, new_data = ., PI= F, n=5))
Sys.time()

dat_grid_pred_1 <- dat_grid %>%
  dplyr::mutate(fit_simple_gam = predict.bam(simple_gam,type = "response",
                                             newdata = .)) %>%
  dplyr::bind_cols(predict_intervals(object_fit = simple_gam, new_data = ., PI= F, n=5)) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(fit_simple_gam = sum(fit_simple_gam),
                   CI_upper_gam = sum(CI_upper), CI_lower_gam = sum(CI_lower)) %>%
  dplyr::ungroup()

Sys.time()
library(mgcv)
library(readr)
library(sspm)
library(tidyr)
library(arrow)

trawl_data <- readr::read_table("PB_fall.dat")
head(trawl_data)

simple_gam <- bam(biomass ~
                    te(long, lat, year, bs= c("tp", "tp"), d = c(2,1)),
                  data = trawl_data, method="REML", family = "tw")
dat_grid_x_y <- as.data.frame(expand_grid(long=seq(0,500,by=0.5),
                                          lat= seq(0,500,by=0.5)))

dat_grid_year <- c(1991:2020)

file_list <- list.files("Result/", full.names=TRUE)
unlink(file_list)
# for (year in dat_grid_year) {
#   print(Sys.time())
#   dat_per_year <- dplyr::bind_cols(dat_grid_x_y,year=year) %>% 
#     dplyr::mutate(fit_simple_gam = predict.bam(simple_gam,type = "response",
#                                                newdata = .)) %>%
#     dplyr::bind_cols(predict_intervals(object_fit = simple_gam, new_data = ., PI= F, n=5))
#   col.names=ifelse(year==dat_grid_year[1],T,F)
#   write.table(dat_per_year, file = 'Result/Predictions.csv', sep = ",",
#               append = TRUE, quote = FALSE,
#               col.names = col.names, row.names = FALSE)
#   print(Sys.time())
# }
Sys.time()
for (year in dat_grid_year) {
  print(Sys.time())
  dat_per_year <- dplyr::bind_cols(dat_grid_x_y,year=year) %>% 
    dplyr::mutate(fit_simple_gam = predict.bam(simple_gam,type = "response", newdata = .)*1000) %>%
    dplyr::bind_cols(predict_intervals(object_fit = simple_gam, new_data = ., PI= F, n=5)*1000)
  # write_feather(dat_per_year,paste0('Result/model_', year))
  print(year)
}
Sys.time()
run_gam_predict <- function(dat_grid_x_y,year){
  print(Sys.time())
  dat_per_year <- dplyr::bind_cols(dat_grid_x_y,year=year) %>% 
  dplyr::mutate(fit_simple_gam = predict.bam(simple_gam,type = "response", newdata = .)*1000) %>%
  dplyr::bind_cols(predict_intervals(object_fit = simple_gam, new_data = ., PI= F, n=5)*1000)
  write_feather(dat_per_year,paste0('Result/model_', year))
  print(year)
}
data_prediction<- list()
  lapply(FUN=run_gam_predict,
         dat_grid_x_y=dat_grid_x_y,
         year=dat_grid_year
         )

file_list <- list.files("Result/", full.names=TRUE)
allData <- plyr::ldply(as.list(file_list), readr::read_csv)
Predictions <- data.table::as.data.table(allData)
write_parquet(Predictions,'Predictions')
Predictions_summary <- Predictions[, lapply(.SD,mean), by=.(year)]
write_parquet(Predictions_summary,'Predictions_summary')
Sys.time()

#automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "palmerpenguins",
  "tidyverse",
  "kableExtra"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}

parallel::detectCores()
n.cores <- parallel::detectCores() - 6
n.cores <- parallelly::availableCores() - 1


#create the cluster
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
memory.limit(200000)
Sys.time()
x <- foreach(
  year = dat_grid_year,
  .packages = c('mgcv','dplyr','sspm','arrow')
) %dopar% {
  dat_per_year <- dplyr::bind_cols(dat_grid_x_y,year=year) %>% 
    dplyr::mutate(fit_simple_gam = predict.bam(simple_gam,type = "response", newdata = .)*1000) %>%
    dplyr::bind_cols(predict_intervals(object_fit = simple_gam, new_data = ., PI= F, n=5)*1000)
  # write.table(dat_per_year, file = paste0('Result/model_', year,'.csv'), sep = ",",row.names = F)
  write_parquet(dat_per_year,paste0('Result/model_', year))
  gc()
}
Sys.time()
x
parallel::stopCluster(cl = my.cluster)
