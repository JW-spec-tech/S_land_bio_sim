

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


Surface_analysis <- function(S_year =1991) {
  
  
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
    dir.create(paste0(d,"/Result_Graphs"), recursive = T)
    sim_files <- sort(list.files(paste0(d,"/sim")))
    
    # Load real biomass
    sim_biomass_predict_total <- list()
    
    # Load predicted biomass
    predict_files <- sort(list.files(paste0(d,"/Result")))
    
    sim_biomass_predict_total <- list()
    
    for (f in predict_files) {
      counter = as.numeric(substr(f, 7,nchar(f)))-1990
      sim_biomass_total <- read_parquet(paste0(d,"/sim/sim",counter)) %>% 
                 summarise(., real_sum=sum(biomass)) %>% 
                 mutate(year=(1990+counter))
      sim_biomass_predicted <- read_parquet(paste0(d,"/Result/",f))
      sim_biomass_predicted$year <- as.numeric(sim_biomass_predicted$year)
      
      sim_real_predict<- left_join(sim_biomass_predicted,sim_biomass_total, by='year')
      
      sim_biomass_predict_total[[f]] <- sim_real_predict
    }

    comparison <- list()
    for (n in 1:length(sim_biomass_predict_total)) {
      comparison[[n]] <- sim_biomass_predict_total[[n]] %>%
        dplyr::mutate(GAM_CI = real_sum >= lower & real_sum <= upper)
    }
    
    
    comparison <- do.call(rbind.data.frame, comparison)
    
    Gam_percent   <- sum(comparison$GAM_CI, na.rm = TRUE)/nrow(comparison)*100
    
    write.table(Gam_percent, paste0(d,"/GAM.percent"))
    
    # #### PErcentage of tim in the interval 
    # Ogmap_percent <- sum(comparison$Ogmap_CI, na.rm = TRUE)/nrow(comparison)*100
    # Gam_percent   <- sum(comparison$GAM_CI, na.rm = TRUE)/nrow(comparison)*100
    # Result_CI <- data.frame(
    #   model=c("Ogmap_percent","Gam_percent") ,  
    #   value=c(Ogmap_percent,Gam_percent)
    # )
    # # Barplot
    # plot <- ggplot(Result_CI, aes(x=model, y=value)) + 
    #   geom_bar(stat = "identity", fill='lightblue', color ='black')+
    #   geom_text(aes(label=paste0(value,'%')),  
    #             position = position_dodge(width = 1),
    #             vjust = 7)+
    #   ggtitle('Percentage of time the simulated biomass falls within the model CI
    #            20 simulations to run for parameters (predict_intervals)
    #            500 samples using 2% of the entire dataset')+
    #   theme(plot.title = element_text(hjust = 0.5))
    # plot
    #  
    # ggsave("plot_test.png",
    #        plot = plot,
    #        device = "png")
    
    #### 4. GC ####
    gc()
    
  }
  
}

Surface_analysis()


parallel::stopCluster(cl = my.cluster)



files <-  list.dirs("Data_test/", recursive = F, full.names = TRUE)

Gam_percent_results <- list()
for (f in files) {
  Gam_percent_results[[f]] <- (read.table(paste0(f,"/GAM.percent")))
}

Gam_percent_results <- do.call(rbind.data.frame, Gam_percent_results)

Surface_analysis_results <- list()
for (f in files) {
  Surface_analysis_results[[f]] <- (read.table(paste0(f,"/Result_Surface/Surface.charaterization")))
}

Surface_analysis_results <- do.call(rbind.data.frame, Surface_analysis_results)



colnames(Surface_analysis_results) <- c("RC","SBI")

graph1 <- data.frame(Surface_analysis_results$RC,Gam_percent_results$x)

graph1 <- graph1[order(graph1$Gam_percent_results.x),c(1,2)]

plot(graph1$Gam_percent_results.x,graph1$Surface_analysis_results.RC)


print(paste("End of Surface analysis @",Sys.time()))
