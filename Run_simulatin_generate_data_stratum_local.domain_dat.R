#### Load packages ####
library(arrow)
library(NLMR)
library(sf)
library(dplyr)
library(raster)
library(fasterize)
library(sspm)
library(rgeos)
library(tidyr)
library(readr)


#### Load the Functions ####
source(file = "Simulation_Depth_Temp_Coordinates_Biomass_function.R")
source(file = "function_patches.R")
source(file = "Make_local.domain_local.arena.R")
source(file = "Make_PB_fall.dat.R")

#### 1. Run the sim ####
results <- S_land_bio_sim(50,500,500,50,5,0.2,0.5)

#### 2. Write the ogmap files ####
Make_patch_domain_arena_DAT(patches=results$patches_list$patches,the_stack=results$the_stack)

#### 3. Slice data file ####
Make_PB_fall.dat()



################################################
library(dplyr)
# Total biomass
F_data <- readr::read_table("PB_fall.dat.complete")
biomass_year <- F_data %>% 
  dplyr::group_by(year) %>% 
  dplyr::select(biomass,year) %>% 
  dplyr::summarise(tbio=sum(biomass))

View(biomass_year)

trawl_year <- data.frame(table(s_data$year)) 



hist_biomass <- s_data %>% 
           filter(year==1991) %>% 
            dplyr::select(biomass)
hist(hist_biomass$biomass)



gf_raster(lat~long, fill = ~log(biomass), data = F_data,
         show.legend = NA) %>% 
  gf_facet_wrap(~year)+
  scale_fill_viridis_c()



.D <- expand.grid(x = 0:5, y = 0:5)
D$z <- runif(nrow(D))
gf_tile(y ~ x, fill = ~z, data = D)
gf_tile(z ~ x + y, data = D)



setwd("C:/Users/jpw/OneDrive/Masters 2020-2022/Main Project/Ogmap2 - Demo")
ogmap_estimates <- readr::read_table("biomass__Ogmap2 - Demo.log")
names(ogmap_estimates)=ogmap_estimates[2,]
ogmap_estimates <- ogmap_estimates[-c(1:2),]

gam_biomass_year <- dat_grid %>% 
  dplyr::group_by(year) %>% 
  dplyr::select(fit_simple_gam,fit_spatialonly_gam,year) %>% 
  dplyr::summarise(tfit_simple_gam=sum(fit_simple_gam)/1000,tfit_spatialonly_gam=sum(fit_spatialonly_gam)/1000)



comparison <- data.frame(biomass_year, ogmap=(as.numeric(ogmap_estimates$Estimate)*1000),)

comparison <- comparison %>% 
              mutate(percent=tbio/ogmap)

summary(comparison$percent)


years <- c(1990:(1990+50))
year_commas <- paste(years, sep = "", collapse = ",")
write.table(year_commas,"years.txt")
