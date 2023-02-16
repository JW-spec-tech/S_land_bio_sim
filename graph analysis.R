
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
      scale_y_continuous(limits = c(95000, NA), expand = c(0,0))
    
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
      scale_y_continuous(limits = c(95000, NA), expand = c(0,0))
    
    CI_plot <- grid.arrange(CI_plot_gam_ogmap,CI_plot_Strap)
    
    
    plot_name <- paste0(basename(getwd()),"_graph.png")
    
    ggsave(plot_name,
            plot = CI_plot,
            device = "png",
            width = 3000,
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

# 3. Generate Graph CI interval vs variation
Generate_CI_Graph_vs_Var <- function(f_names=f_list) {
  
  CI_results_list <- list()
  counter=1
  
  for (folder in f_names) {
    CI_results_list[[counter]] <- read_parquet(paste0(folder,"/CI_results"))
    counter=counter+1
  }
  
  CI_results_df <-  do.call(rbind.data.frame, CI_results_list)
  
  f_names <- sub("^.*Experiment_", "", f_names)
  print(f_names)
  
  rownames(CI_results_df) <- f_names
  
  # Define desired order of row names
  desired_order <- c("Variation_1", "Variation_5", "Variation_10", "Variation_25")
  
  # Create factor with desired order
  row_order <- factor(rownames(CI_results_df), levels = desired_order)
  
  # Sort dataframe by factor
  CI_results_df <- CI_results_df[order(row_order),]
  
  # Extract numeric portion of row names
  variation <- as.numeric(gsub("Variation_", "", rownames(CI_results_df)))
  
  # Add variation column to dataframe
  CI_results_df$variation <- variation
  
  library(ggplot2)
  
  # Define plot data
  plot_data <- data.frame(
    variation = CI_results_df$variation,
    Ogmap = CI_results_df$percentage_of_true_Ogmap,
    GAM = CI_results_df$percentage_of_true_GAM,
    STRAP = CI_results_df$percentage_of_true_STRAP
  )
  
  # Create ggplot with percentage of true values as a function of variation, with lines between points and title
  ggplot(plot_data, aes(x = variation)) +
    geom_line(aes(y = Ogmap, color = "Ogmap")) +
    geom_point(aes(y = Ogmap, color = "Ogmap"), size = 3) +
    geom_line(aes(y = GAM, color = "GAM")) +
    geom_point(aes(y = GAM, color = "GAM"), size = 3) +
    geom_line(aes(y = STRAP, color = "STRAP")) +
    geom_point(aes(y = STRAP, color = "STRAP"), size = 3) +
    scale_color_manual(values = c("Ogmap" = "red", "GAM" = "green", "STRAP" = "blue")) +
    labs(x = "Variation", y = "Percentage of true values", 
         title = "Percentage of the time the Confidence interval of models 
         captures the true Biomass at different level of Landscape Variation")
  
}

# 5. Generate SBI and RC vs variation

Generate_CI_Graph_Surface_vs_CI <- function(f_names=f_list) {
  
  CI_results_list <- list()
  counter=1
  
  for (folder in f_names) {
    CI_results_list[[counter]] <- read_parquet(paste0(folder,"/CI_results"))
    counter=counter+1
  }
  
  CI_results_df <-  do.call(rbind.data.frame, CI_results_list)
  
  f_names_col <- sub("^.*Experiment_", "", f_names)
  print(f_names_col)
  
  
  Surface_results_list <- list()
  counter=1
  
  for (folder in f_names) {
    Surface_results_list[[counter]] <- read.table(paste0(folder,"/Surface_roughness_mean"))
    counter=counter+1
  }
  
  Surface_results_df <-  do.call(cbind.data.frame, Surface_results_list)
  
  f_names_col <- sub("^.*Experiment_", "", f_names)
  
  colnames(Surface_results_df) <- f_names_col
  
  
  # Extract the numbers from the column names
  col_nums <- as.numeric(gsub("Variation_", "", colnames(Surface_results_df)))
  
  # Order the column numbers and get the corresponding column names
  sorted_col_names <- colnames(Surface_results_df)[order(col_nums)]
  
  # Reorder the columns in the data frame using the sorted column names
  Surface_results_df <- Surface_results_df[, sorted_col_names]
  
  # Transpose
  Surface_results_df <- t(Surface_results_df)
  
  CI_results_df$RC <- Surface_results_df[,1]
  
  CI_results_df$SBI <- Surface_results_df[,2]
  
  CI_results_df$variation <- as.charCI_results_df$variation
  
  library(ggplot2)
  library(gridExtra)
  
  # add row names as a column in the data frame
  CI_results_df$name <- rownames(CI_results_df)
  
  # define the colors for each column
  colours <- c("percentage_of_true_Ogmap" = "blue", "percentage_of_true_GAM" = "red", "percentage_of_true_STRAP" = "green")
  

  
  # plot 1: percentage_of_true_Ogmap, percentage_of_true_GAM, and percentage_of_true_STRAP as a function of SBI
  p1 <- ggplot(CI_results_df, aes(x = log10(SBI), group = factor(variation))) +
    geom_point(aes(y = percentage_of_true_Ogmap, color = "Ogmap"), size = 3) +
    geom_point(aes(y = percentage_of_true_GAM, color = "GAM"), size = 3) +
    geom_point(aes(y = percentage_of_true_STRAP, color = "STRAP"), size = 3) +
    labs(x = "log10(SBI)", y = "Percentage of True Biomass Captured") +
    ggtitle("SBI vs. Percentage of True Biomass Captured") +
    scale_color_manual(values = c("blue", "red", "green")) +
    geom_text(aes(label = factor(variation), y = percentage_of_true_Ogmap),nudge_x = 0.006) +
    geom_text(aes(label = factor(variation), y = percentage_of_true_GAM),nudge_x = 0.006) +
    geom_text(aes(label = factor(variation), y = percentage_of_true_STRAP),nudge_x = 0.006)
  
  # plot 2: percentage_of_true_Ogmap, percentage_of_true_GAM, and percentage_of_true_STRAP as a function of log10(RC)
  p2 <- ggplot(CI_results_df, aes(x = log10(RC), group = factor(variation))) +
    geom_point(aes(y = percentage_of_true_Ogmap, color = "Ogmap"), size = 3) +
    geom_point(aes(y = percentage_of_true_GAM, color = "GAM"), size = 3) +
    geom_point(aes(y = percentage_of_true_STRAP, color = "STRAP"), size = 3) +
    labs(x = "log10(RC)", y = "Percentage of True Biomass Captured") +
    ggtitle("RC vs. Percentage of True Biomass Captured") +
    scale_color_manual(values = c("blue", "red", "green")) +
    geom_text(aes(label = factor(variation), y = percentage_of_true_Ogmap),nudge_x = 0.5) +
    geom_text(aes(label = factor(variation), y = percentage_of_true_GAM),nudge_x = 0.5) +
    geom_text(aes(label = factor(variation), y = percentage_of_true_STRAP),nudge_x = 0.5)
  
  # combine the two plots into one using grid.arrange
  grid.arrange(p1, p2, ncol = 2)
  
  
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
