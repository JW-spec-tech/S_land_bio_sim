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
library(mgcv)
library(readr)
library(foreach)
library(doParallel)
library(parallelly)
library(ranger)
library(tidyverse)
library(kableExtra)
library(ggplot2)
library(patchwork)
library(viridis)





#### Load the Functions ####
source(file = "function_patches.R")
source(file = "Make_local.domain_local.arena.R")
source(file = "Strap_calculations.R")
source(file = "Make_PB_fall.dat.R")


#### Load model - TODO ####
# have the depth model here

#### Set Simulation Parameters ####

# Landscape dimensions
x=y=size <- 500

# Number of years/simulations
n <- 1

# Vairation
v1 <- 1
v2 <- v1*2
v3 <- v1*3

#### 1. Generate base Landscape ####

message("1. Generate base Landscape")

Main_L <- nlm_mpd(ncol = size, nrow = size, roughness = roughness)*nlm_mpd(ncol = size, nrow = size, roughness = roughness)*nlm_mpd(ncol = size, nrow = size, roughness = roughness)

plot(Main_L)

#### 1.B Generate secondary landscapes which will be used to vary temperature.####

message("1.B Generate secondary landscapes which will be used to vary temperature.")

list_x <- rep(x-1, n)
list_y <- rep(y-1, n)

Sub_L <- mapply(FUN = nlm_gaussianfield,
                ncol = list_x, nrow = list_y, 
                resolution = 1,
                autocorr_range = 1,
                mag_var = 5,
                nug = 0.2,
                mean = 0.5)

plot(Sub_L[[1]])

#### 1.C Generate secondary depths which will be used to generate alternate temperatures for the different years ####
Sub_L_M <- list()
for(i in 1:n){
  Sub_L_M[[i]] <- (Main_L + (Sub_L[[i]]/10))/2
}
plot(Sub_L_M[[1]])

#### 1.D Vizualizing before and after Variation

df1 <- as.data.frame(Main_L, xy = TRUE)
df1$layer <- df1$layer/max(df1$layer)
df2 <- as.data.frame(Sub_L_M[[1]], xy = TRUE)
df2$layer <- df2$layer/max(df2$layer)

plot1 <- ggplot(df1, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  theme_minimal() +
  labs(title = "Before yearly variation") +
  viridis::scale_fill_viridis()

plot2 <- ggplot(df2, aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  theme_minimal() +
  labs(title = "After yearly variation") +
  viridis::scale_fill_viridis()

combined_plot <- plot1 + plot2
combined_plot

#### 2. Create Depth patches ####

message("2. Create Depth patches")

 # Add patches of depth variation

 depth_patch_variation <- nlm_randomcluster(ncol = size-1, nrow = size-1,

                                        p = 0.54,
                                        ai = c(0.7, 0.10, 0.10, 0.10))
 
  plot(depth_patch_variation)
  
 # Multiply the landscape depth at every location
  Sub_L_M_patch_v1 <- exp(depth_patch_variation*v1)*Sub_L_M[[1]]@data@values
  Sub_L_M_patch_v1@data@values <- Sub_L_M_patch_v1@data@values/max(Sub_L_M_patch_v1@data@values)
  Sub_L_M_patch_v2 <- exp(depth_patch_variation*v2)*Sub_L_M[[1]]@data@values
  Sub_L_M_patch_v2@data@values <- Sub_L_M_patch_v2@data@values/max(Sub_L_M_patch_v2@data@values)
  Sub_L_M_patch_v3 <- exp(depth_patch_variation*v3)*Sub_L_M[[1]]@data@values
  Sub_L_M_patch_v3@data@values <- Sub_L_M_patch_v3@data@values/max(Sub_L_M_patch_v3@data@values)
  
#### 3. Generate depth ####
  
  Sub_L_M_patch_v1@data@values <- (Sub_L_M_patch_v1@data@values*1126)+58
  Sub_L_M_patch_v2@data@values <- (Sub_L_M_patch_v2@data@values*1126)+58
  Sub_L_M_patch_v3@data@values <- (Sub_L_M_patch_v3@data@values*1126)+58

  Sub_L_M_patch_v1_df <- as.data.frame(Sub_L_M_patch_v1, xy = TRUE)
  Sub_L_M_patch_v2_df <- as.data.frame(Sub_L_M_patch_v2, xy = TRUE)
  Sub_L_M_patch_v3_df <- as.data.frame(Sub_L_M_patch_v3, xy = TRUE)
  
  # Rename the 3rd column
  names(Sub_L_M_patch_v1_df)[3] <- "depth"
  names(Sub_L_M_patch_v2_df)[3] <- "depth"
  names(Sub_L_M_patch_v3_df)[3] <- "depth"
  
#### 4. Generate Temperature ####
  
  message("4. Generate Temperature")
  
  real <- read.csv("/trawl_nl.csv")

  # Generate gam based on real depth and temp
  gam_depth_sim <- gam(temp_bottom ~ s(sqrt(depth), bs="ad"), data = real)

  
  # Predict the temperature from depth
  Sub_L_M_patch_v1_df$temperature <- predict(gam_depth_sim, newdata = Sub_L_M_patch_v1_df, se.fit = T)$fit
  Sub_L_M_patch_v2_df$temperature <- predict(gam_depth_sim, newdata = Sub_L_M_patch_v2_df, se.fit = T)$fit
  Sub_L_M_patch_v3_df$temperature <- predict(gam_depth_sim, newdata = Sub_L_M_patch_v3_df, se.fit = T)$fit

#### 5. Generate Biomass ####
  
  message("5. Generate Biomass")
  

  # Biomass parameters #
  depth_sd = 200
  temp_sd  = 2
  scale_depth <- dnorm(0,0,depth_sd)
  scale_temp  <- dnorm(0,0,2)
  
  Sub_L_M_patch_v1_df$biomass <- dnorm((Sub_L_M_patch_v1_df$depth - 312.5), 0, 100)/dnorm(0,0,100)*dnorm((Sub_L_M_patch_v1_df$temperature - 2.916), 0, 2)/dnorm(0,0,2)
  Sub_L_M_patch_v2_df$biomass <- dnorm((Sub_L_M_patch_v2_df$depth - 312.5), 0, 100)/dnorm(0,0,100)*dnorm((Sub_L_M_patch_v2_df$temperature - 2.916), 0, 2)/dnorm(0,0,2)
  Sub_L_M_patch_v3_df$biomass <- dnorm((Sub_L_M_patch_v3_df$depth - 312.5), 0, 100)/dnorm(0,0,100)*dnorm((Sub_L_M_patch_v3_df$temperature - 2.916), 0, 2)/dnorm(0,0,2)


#### 7. Plotting ####
  
  library(ggplot2)
  library(viridis)
  library(cowplot)
  
  # Function to create stacked plot for each data frame and variable
  create_stacked_plot <- function(df, variable, title) {
    plot <- ggplot(df, aes(x = x, y = y, fill = .data[[variable]])) +
      geom_tile() +
      scale_fill_viridis() +
      labs(title = title)
    
    return(plot)
  }
  
  # Create the stacked plots for each data frame and variable
  plot_v1_depth <- create_stacked_plot(Sub_L_M_patch_v1_df, "depth", "Patch v1 - Depth")
  plot_v1_temperature <- create_stacked_plot(Sub_L_M_patch_v1_df, "temperature", "Patch v1 - Temperature")
  plot_v1_biomass <- create_stacked_plot(Sub_L_M_patch_v1_df, "biomass", "Patch v1 - Biomass")
  
  plot_v2_depth <- create_stacked_plot(Sub_L_M_patch_v2_df, "depth", "Patch v2 - Depth")
  plot_v2_temperature <- create_stacked_plot(Sub_L_M_patch_v2_df, "temperature", "Patch v2 - Temperature")
  plot_v2_biomass <- create_stacked_plot(Sub_L_M_patch_v2_df, "biomass", "Patch v2 - Biomass")
  
  plot_v3_depth <- create_stacked_plot(Sub_L_M_patch_v3_df, "depth", "Patch v3 - Depth")
  plot_v3_temperature <- create_stacked_plot(Sub_L_M_patch_v3_df, "temperature", "Patch v3 - Temperature")
  plot_v3_biomass <- create_stacked_plot(Sub_L_M_patch_v3_df, "biomass", "Patch v3 - Biomass")
  
  # Arrange the plots in a 3x3 grid
  grid_plot <- plot_grid(
    plot_v1_depth, plot_v2_depth, plot_v3_depth,
    plot_v1_temperature, plot_v2_temperature, plot_v3_temperature,
    plot_v1_biomass, plot_v2_biomass, plot_v3_biomass,
    ncol = 3
  )
  
  # Display the grid plot
  print(grid_plot)
  
  
  
  
  

  