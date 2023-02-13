
library(arrow)
library(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
getwd()
source("Gam_analysis.R")


#  1. Compile Gam Data
Analysis_graph <- function(S_year =1991) {
  
  start_year =S_year
  years = files <-  20
  size = 499
  
  for (cwd in list.dirs(full.names = T,recursive = F)) {
    
    print(cwd)
    setwd(cwd)
    cwd=getwd()
    # Analyse_Gam()
    Graphing_1()
    setwd("../")
  }
}


# 2. Generate biomass stas
biomass_Stats <- function(S_year =1991) {
  
  start_year =S_year
  years = files <-  20
  size = 499
  
  for (cwd in list.dirs(full.names = T,recursive = F)) {
    
    print(cwd)
    setwd(cwd)
    cwd=getwd()
  
  #### 1. Load sim names
  file_list <- list.files("Results/", full.names=TRUE)
  
  #### 2. Merge sims
  biomass_data_GAM <- plyr::ldply(as.list(file_list), arrow::read_parquet)

# Read Ogmap data
# Define the path to the directory containing the files
dir_path <- cwd

# Get a list of all the files in the directory
files <- list.files(dir_path)

# Filter the list of files to only include those that start with "biomass"
biomass_files <- files[grep("^biomass", files)]

biomass_file <- biomass_files[1]
biomass_data_Ogmap <- readr::read_table(biomass_file, skip=2)

# Load strap Data

biomass_data_STRAP <- read_parquet("Strap_estimate")

# Actual stats
F_data <- arrow::read_parquet("PB_fall.dat.complete")


biomass_data_true <- F_data %>%
  dplyr::group_by(year) %>%
  dplyr::select(biomass,year) %>%
  dplyr::summarise(t_bio=sum(biomass)
  )

biomass_all <- data.frame(year=biomass_data_true$year,
                          true_bio=biomass_data_true$t_bio,
                          Ogmap_mean_bio=biomass_data_Ogmap$Estimate*1000,
                          Ogmap_Upper_bio=biomass_data_Ogmap$UpCIval*1000,
                          Ogmap_Lower_bio=biomass_data_Ogmap$LowCIval*1000,
                          Gam_mean_bio=biomass_data_GAM$point_est,
                          Gam_Upper_bio=biomass_data_GAM$upper,
                          Gam_Lower_bio=biomass_data_GAM$lower,
                          STRAP_mean_bio=biomass_data_STRAP$B_total,
                          STRAP_UPPER_bio=biomass_data_STRAP$upper,
                          STRAP_Lower_bio=biomass_data_STRAP$lower)

write.table(biomass_all,"Biomass_stats")

setwd("../")
  }

}

# 3. Generate Graph CI interval
Generate_CI_Graph <- function(S_year =1991) {
  
  start_year =S_year
  years = files <-  20
  size = 499
  
  for (cwd in list.dirs(full.names = T,recursive = F)) {
    
    print(cwd)
    setwd(cwd)
    cwd=getwd()
    
    # Read data 
    Data <- read.table("biomass_Stats")
    
    #graph Confidence Interval coverage per year
    # Create ggplot object
    data_Graph <- ggplot(Data, aes(year))
    
    # Create graph
    colors <- c("GAM" = "red", "Ogmap" = "lightblue", "Total Biomass" = "black")
    
    CI_plot_gam_ogmap <- data_Graph +
      geom_ribbon(aes(ymin = Gam_Lower_bio, ymax = Gam_Upper_bio), fill = "red", alpha=0.25)+
      geom_line(aes(y=Gam_mean_bio, color= "GAM"))+
      geom_ribbon(aes(ymin = Ogmap_Lower_bio, ymax = Ogmap_Upper_bio), fill = "lightblue", alpha=0.5)+
      geom_line(aes(y=Ogmap_mean_bio, color = "Ogmap"))+
      # geom_ribbon(aes(ymin = CI_STRAP_lower, ymax = CI_STRAP_upper), fill = "green", alpha=0.25)+
      # geom_line(aes(y=STRAP), color = "green")+
      geom_line(aes(y=true_bio, color = "Total Biomass"))+
      scale_color_manual(name= "Biomass", labels = c("GAM","OGmap", "Total Biomass"), values = colors)+
      labs(title = "Coverage of Confidence Intervals
GAM VS OGmap VS STRAP", x="Year (simulation #)", y="Biomass in kg")+
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_y_continuous(limits = c(80000, NA), expand = c(0,0))
    
    colors2 <- c("STRAP" = "orange", "Total Biomass" = "black")
    
    CI_plot_Strap <- data_Graph +
      # geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha=0.25)+
      # geom_line(aes(y=Gam_mean_bio), color= "red")+
      # geom_ribbon(aes(ymin = Ogmap_Lower_bio, ymax = Ogmap_Upper_bio), fill = "lightblue", alpha=0.5)+
      # geom_line(aes(y=ogmap), color = "lightblue")+
      geom_ribbon(aes(ymin = STRAP_Lower_bio, ymax = STRAP_UPPER_bio), fill = "orange", alpha=0.25)+
      geom_line(aes(y=STRAP_mean_bio, color = "STRAP"))+
      geom_line(aes(y=true_bio, color = "Total Biomass"))+
      scale_color_manual(name= "Biomass", labels = c("STRAP","Total biomass"),values = colors2)+
      labs(title = "Coverage of Confidence Intervals
       STRAP", x="Year (simulation #)", y="Biomass in kg")+
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_y_continuous(limits = c(80000, NA), expand = c(0,0))
    
    CI_plot <- grid.arrange(CI_plot_gam_ogmap,CI_plot_Strap)
    
    
    plot_name <- paste0(basename(getwd()),"_graph")
    
    ggsave(plot_name,
            plot = CI_plot,
            device = "png",
            width = 24000,
            height = 1835,
            units = "px",
            limitsize = F)
    
    setwd("../")
  }
    
    
  }


# 4. Compile biomass files
Compile_biomass <- function(S_year =1991) {
  
  start_year =S_year
  years = files <-  20
  size = 499
  
  
  DATA_list <- list()
  for (cwd in list.dirs(full.names = T,recursive = F)) {
    
    print(cwd)
    setwd(cwd)
    counter <- substr(cwd, 7, 7)
    cwd=getwd()
    # Read data 
    DATA_list[[counter]] <- read.table("biomass_Stats")
    
    setwd("../")
    
  }
  
  df_DATA_list <- as.data.frame(do.call(rbind, DATA_list))
  
  write_parquet(df_DATA_list,"Summary_data")
}

# 4. CI interval calculation
CI_calc_percent <- function(S_year =1991) {
  start_year =S_year
  years = files <-  20
  size = 499
  
  S_DATA <- read_parquet("Summary_data")
  
  head(S_DATA)

  S_DATA$result_Ogmap <- ifelse(S_DATA$true_bio >= S_DATA$Ogmap_Lower_bio & S_DATA$true_bio <= S_DATA$Ogmap_Upper_bio, TRUE, FALSE)
  head(S_DATA)
  S_DATA$result_GAM <- ifelse(S_DATA$true_bio >= S_DATA$Gam_Lower_bio & S_DATA$true_bio <= S_DATA$Gam_Upper_bio, TRUE, FALSE)
  head(S_DATA)
    S_DATA$result_STRAP <- ifelse(S_DATA$true_bio >= S_DATA$STRAP_Lower_bio & S_DATA$true_bio <= S_DATA$STRAP_UPPER_bio, TRUE, FALSE)
  head(S_DATA)
  
  
  percentage_of_true_Ogmap <- sum(S_DATA$result_Ogmap == TRUE) / nrow(S_DATA) * 100
  percentage_of_true_Ogmap
  
  percentage_of_true_GAM <- sum(S_DATA$result_GAM == TRUE) / nrow(S_DATA) * 100
  percentage_of_true_GAM
  
  percentage_of_true_STRAP <- sum(S_DATA$result_STRAP == TRUE) / nrow(S_DATA) * 100
  percentage_of_true_STRAP
  
  CI_results <- data.frame(percentage_of_true_Ogmap=percentage_of_true_Ogmap,percentage_of_true_GAM=percentage_of_true_GAM,percentage_of_true_STRAP=percentage_of_true_STRAP)
  write_parquet(CI_results,"CI_results")
  }



# Get files names
f_list <- paste0(getwd(),"/",list.dirs(path = "exp", full.names = TRUE, recursive = F))

# Load sim data

for (i in f_list) {
  print(i)
  setwd(i)
  Analysis_graph()
  setwd("~/Git projects/S_land_bio_sim")
}

for (i in f_list) {
  print(i)
  setwd(i)
  biomass_Stats()
  setwd("~/Git projects/S_land_bio_sim")
}

for (i in f_list) {
  print(i)
  setwd(i)
  Generate_CI_Graph()
  setwd("~/Git projects/S_land_bio_sim")
}


for (i in f_list) {
  print(i)
  setwd(i)
  Compile_biomass()
  setwd("~/Git projects/S_land_bio_sim")
}


for (i in f_list) {
  print(i)
  setwd(i)
  CI_calc_percent()
  setwd("~/Git projects/S_land_bio_sim")
}
