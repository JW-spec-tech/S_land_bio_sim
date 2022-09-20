#  Goal of this is to resample landscapes and creat a distribution

library(dplyr)
library(arrow)


#  Read full simulation data
F_data <- arrow::read_parquet("PB_fall.dat.complete")

#  #function to create a single sample
Resample_sims <- function(n,percent_f=0.025,path="PB_fall.dat.complete"){
  

  
  propotion_strata <- F_data %>% 
    dplyr::group_by(stratum,year) %>% 
    dplyr::summarise(n=n()) %>% 
    dplyr::mutate(min=3)
  
  propotion_strata$prop = round(ifelse(propotion_strata$n*percent_f<=propotion_strata$min,propotion_strata$min,ifelse(propotion_strata$n*percent_f>=15,15,propotion_strata$n*percent_f)))
  
  propotion_strata$prop_random = propotion_strata$prop + sample(-2:2,length(propotion_strata$prop),replace = T)
  
  propotion_strata$prop_random= ifelse(propotion_strata$prop_random<propotion_strata$min,propotion_strata$min,propotion_strata$prop_random)
  
  # hist(propotion_strata$prop_random) # testing shows the histogram of # of trawls per strata
  
  F_data_prop <- dplyr::left_join(F_data,propotion_strata[,c("year","stratum","prop_random","n")])
  
  # Initialize list
  Resampling_list <- list()
  
  for (n in 1:n) {
    
    #  Sample the landscape
    single_sample <- F_data_prop %>% 
      group_by(stratum,year)  %>% 
      group_split() %>%
      purrr::map_dfr(~slice_sample(.x,n=.x$prop_random[1])) %>% 
      bind_rows() %>% 
      mutate(biomass=rTweedie((biomass*1000), p = 1.76, phi= 2))
    
    # store said sample
    Resampling_list[[n]] <- single_sample

  }
  return(Resampling_list)
}

Resample_sims_random <- function(n,percent_f=0.025,path="PB_fall.dat.complete"){
  
  # Initialize list
  Resampling_list <- list()
  
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
  
  n=2000
  
  Resampling_list <- foreach(
    n = 1:n,
    .packages = c('mgcv','dplyr','purrr')
  ) %dopar% {
    propotion_strata <- F_data %>% 
      dplyr::group_by(stratum,year) %>% 
      dplyr::summarise(n=n()) %>% 
      dplyr::mutate(min=3)
    
    propotion_strata$prop = round(ifelse(propotion_strata$n*percent_f<=propotion_strata$min,propotion_strata$min,ifelse(propotion_strata$n*percent_f>=15,15,propotion_strata$n*percent_f)))
    
    propotion_strata$prop_random = propotion_strata$prop + sample(-2:2,length(propotion_strata$prop),replace = T)
    
    propotion_strata$prop_random= ifelse(propotion_strata$prop_random<propotion_strata$min,propotion_strata$min,propotion_strata$prop_random)
    
    # hist(propotion_strata$prop_random) # testing shows the histogram of # of trawls per strata
    
    F_data_prop <- dplyr::left_join(F_data,propotion_strata[,c("year","stratum","prop_random","n")])
    
    
    
    #  Sample the landscape
    single_sample <- F_data_prop %>% 
      group_by(stratum,year)  %>% 
      group_split() %>%
      purrr::map_dfr(~slice_sample(.x,n=.x$prop_random[1])) %>% 
      bind_rows() %>% 
      mutate(biomass=rTweedie((biomass*1000), p = 1.76, phi= 2)) %>% 
      filter(year==1991)
    
    # store said sample
    Resampling_list[[n]] <- single_sample
  }

  return(Resampling_list)
  parallel::stopCluster(cl = my.cluster)
}

parallel::stopCluster(cl = my.cluster)


# Initialize list
Resampling_list <- list()

# Run resampling
Resample_list <- Resample_sims(10000)



# Random resampling
Random_Resample_list <- Resample_sims_random(10,percent_f=0.1)

par(mfrow=c(1,1))

hist(mean(Random_Resample_list[[1]][["biomass"]]))
hist(Random_Resample_list[[2]][["biomass"]])
hist(Random_Resample_list[[3]][["biomass"]])
hist(Random_Resample_list[[4]][["biomass"]])


mean_list <- list() 
for (n in 1:length(Random_Resample_list)) {

     mean <- mean(Random_Resample_list[[n]][["biomass"]])
     
      mean_list[n] <- mean
}

mean_list2 <- list() 
for (n in 1:length(Resampling_list_100)) {
  
  Biomass_mean <- Resampling_list_100[[n]] %>%
                  group_by(stratum,n) %>%
                  dplyr::summarize(mean_biomass= mean(biomass)*n) %>%
                  distinct(stratum, .keep_all = T) %>%
                  ungroup() %>% 
                  dplyr::summarize(mean_biomass= sum(mean_biomass)/(nrow(F_data_prop)/length(F_data_prop)))

  
  # Biomass_mean <- Resampling_list_100
  # 
  # Biomass_mean <- as.data.table(Biomass_mean)
  # 
  # Biomass_mean <- Biomass_mean[,sum(mean(biomass)*n)/sum(n), by = list(stratum,n)]
  # 
  # Biomass_mean <- Biomass_mean[,mean(V1)]

  mean_list2[n] <- Biomass_mean
}

dat <- mtcars[c(2:4,11)]



grp <- function(x) {
  group_by(Resampling_list_100,!!as.name(x)) %>%
    dplyr::summarize(mean_biomass= mean(biomass)*n) %>%
    distinct(stratum, .keep_all = T) %>%
    ungroup()
}

mean_list2 <- grp(stratum)


Resampling_list_100 %>% 
  map(~summarise(group_by(., stratum,n),
                 mean_biomass= mean(biomass)*n))

lapply(colnames(dat), grp)


df_list_2 <- lapply(Resampling_list_100, function(x) {
  aggregate(x$Frequency, by=list(stratum=x$Category), FUN=sum)
})


weighted_means <- lapply(Resampling_list_100, function(x) as.data.table(x)[,sum(mean(biomass)*n)/sum(n), by = list(stratum,n)])

weighted_means <- lapply(Resampling_list_100, function(x) as.data.table(x)[,mean(biomass)*n, by = list(stratum,n)])


Sum_weighted_means <- lapply(weighted_means, function(x) as.data.table(x)[,sum(V1)/(nrow(F_data_prop)/length(F_data_prop))])

sum(mean_biomass)/(nrow(F_data_prop)/length(F_data_prop))

mean_df <- t(as.data.frame(do.call(cbind, Sum_weighted_means)))

hist(mean_df, xlim = c(min(400),max(mean_df)))
abline(v=mean(bio_1991_F$biomass), col = "blue", lwd = 4, lty = 4)


Resample_hist <- function(data) {
  
  weighted_means <- lapply(data, function(x) as.data.table(x)[,mean(biomass)*n, by = list(stratum,n)])
  
  Sum_weighted_means <- lapply(weighted_means, function(x) as.data.table(x)[,sum(V1)/(nrow(F_data_prop)/length(F_data_prop))])
  
  mean_df <- t(as.data.frame(do.call(cbind, Sum_weighted_means)))
  
  hist(mean_df, xlim = c(min(400),max(mean_df)))
  abline(v=mean(bio_1991_F$biomass), col = "blue", lwd = 4, lty = 4)
  
  
}


  Resample_hist(Resampling_list)



for (n in 1:length(Resampling_list)) {
  
  mean <- mean(Resampling_list[[n]][["biomass"]])
  
  mean_list[n] <- mean
}


mean_strata_pre <- data.frame(biomass=Random_Resample_list[[1]][["biomass"]],stratum=Random_Resample_list[[1]][["stratum"]])

mean_test <- mean_strata_pre %>% 
  group_by(stratum) %>% 
  dplyr::summarize(mean_biomass= mean(biomass)) %>% 
  dplyr::ungroup()  %>% 
  dplyr::summarize(mean= mean(mean_biomass))
mean(mean_strata_pre$biomass)


for (n in 1:length(Random_Resample_list)) {
  
  
  for (year in unique(Random_Resample_list[[n]][["year"]])) {
    bio_year_s <- (Random_Resample_list %>% 
                     dplyr::filter(year==1991) %>% 
                     dplyr::select(biomass)) 
  }

  
}


mean_df <- t(as.data.frame(do.call(cbind, mean_list)))



hist(mean_df, xlim = c(min(mean_df),max(mean_df)))
abline(v=mean(bio_1991_F$biomass), col = "blue", lwd = 4, lty = 4)


mean(mean_df)  
bio_1991_F <- (F_data %>% 
       dplyr::filter(year==1991) %>% 
       dplyr::select(biomass)) %>% 
  mutate(biomass=rTweedie((biomass*1000), p = 1.76, phi= 2))

mean(bio_1991_F$biomass)











  