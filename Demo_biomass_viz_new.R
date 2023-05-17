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
x=y=size <- 100

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
 depth_patch_variation <- nlm_randomcluster(ncol = 500, nrow = 500,
                                        p = 0.40,
                                        ai = c(0.7, 0.10, 0.10, 0.10))
 
  plot(depth_patch_variation)
 # Multiply the landscape depth at every location
 for (i in 1:n) {
   Sub_L_M[[i]]@data@values <- exp(depth_patch_variation_2)*Sub_L_M[[i]]@data@values

 }
