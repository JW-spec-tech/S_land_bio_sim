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


#### Functions needed ####

# Function: make_patches
# Purpose: Generate patches from a raster and convert them into polygons.
# Inputs: patch - Raster object representing the main landscape
# Outputs: List containing the rasterized patches and the polygon representation of patches

make_patches <- function(patch){ #,plot=F
  
  # patch <- raster("main_L.gri")
  
  Main_L_copy <- patch
  Main_L_copy_smt <- patch
  
  # Smooth the raster 
  #####Removed this step #######
  
  # Main_L_copy_smt <- raster::focal(Main_L_copy, w=matrix(1,9,9), pad = TRUE,
  # na.rm = TRUE)
  
  # Reclass the raster
  # ?reclassify
  brks <- 16
  rcl_matrix <- matrix(c(seq(0, 1200, length.out = brks),
                         c(seq(0, 1200, length.out = brks)[-1], 1400),
                         1:brks),
                       ncol=3, nrow=brks)
  Main_L_copy_rcl <- reclassify(Main_L_copy_smt, 
                                rcl_matrix)
  # plot(Main_L_copy_rcl)
  # Turn into polygons (sp object)
  Main_L_copy_rcl_clumped <- rasterToPolygons(Main_L_copy_rcl, dissolve = TRUE)
  # plot(Main_L_copy_rcl_clumped)
  
  # Turn the sp object into sf object for easier manipulation
  Main_L_copy_rcl_clumped_sf <- st_as_sf(Main_L_copy_rcl_clumped) |>
    st_cast("POLYGON")
  
  Main_L_copy_rcl_clumped_sf$area <- 0
  # Calculate the area amd cretae 
  Main_L_copy_rcl_clumped_sf$area <- st_area(Main_L_copy_rcl_clumped_sf)
  Main_L_copy_rcl_clumped_sf$bound_id <- 1:nrow(Main_L_copy_rcl_clumped_sf)
  # plot(Main_L_copy_rcl_clumped)
  
  # Identify small polygons and give them the id of their container polygon
  small_polygons <- Main_L_copy_rcl_clumped_sf$bound_id[which(Main_L_copy_rcl_clumped_sf$area < 2500)]
  Main_L_copy_rcl_clumped_sf_edges <- sf::st_intersects(Main_L_copy_rcl_clumped_sf)
  names(Main_L_copy_rcl_clumped_sf_edges) <- Main_L_copy_rcl_clumped_sf$bound_id
  
  # Loop to check each polygon for its size
  for (i in small_polygons) {
    current_polygons <- Main_L_copy_rcl_clumped_sf[Main_L_copy_rcl_clumped_sf_edges[[i]], ] %>%
      dplyr::filter(.data$area == max(.data$area))
    max_id <- current_polygons$bound_id
    Main_L_copy_rcl_clumped_sf$bound_id[Main_L_copy_rcl_clumped_sf$bound_id == i] <- max_id
  }
  
  # Combine small polygons
  Main_L_copy_rcl_clumped_sf <-
    Main_L_copy_rcl_clumped_sf %>%
    dplyr::group_by(.data$bound_id) %>%
    dplyr::summarize() %>%
    dplyr::ungroup()
  
  # plot(Main_L_copy_rcl_clumped_sf)
  
  # Make sure bound id is a factor (bug in sspm)
  Main_L_copy_rcl_clumped_sf <- Main_L_copy_rcl_clumped_sf %>% 
    mutate(bound_id = as.factor(bound_id))
  
  # plot(st_simplify(Main_L_copy_rcl_clumped_sf, dTolerance = 10)$geometry)
  
  
  # using sspm to tesselate
  sspm_boundary <- spm_as_boundary(boundaries = Main_L_copy_rcl_clumped_sf, 
                                   boundary = "bound_id")
  
  # set number of strata per area as a function of size
  nb_nodes_bound_id <- round(as.numeric(sspm_boundary@boundaries$area_bound_id/2500))
  
  # This makes sure we sample the surface of the polygons at random points 
  voronoi <- spm_discretize(sspm_boundary, method = "tesselate_voronoi", 
                            sample_surface = TRUE, nb_samples = nb_nodes_bound_id, min_size = 50)
  
  # Other bug in sspm: we make sure all poygons are of the type POLUYGON and NOT MULTIPOLYGON
  patches <- st_cast(st_make_valid(spm_patches(voronoi), "POLYGON"))
  sf::st_bbox(patches)
  
  sum(patches$patch_area)
  
  # Rasterize the sf object
  patches_raster <- fasterize(patches, raster = Main_L_copy, field = "patch_id")
  
  
  # # Stack the rasters and turn into df values (alignment  )
  # the_stack <- stack(Main_L_copy, patches_raster)
  # names(the_stack) <<- c("Main", "Patches")
  # test <- as.data.frame(the_stack)
  # plot(patches["patch_id"], main = "Strata Generation Simulation")
  # hist(patches$patch_area)
  return(list(patches_raster=patches_raster,patches=patches))
}


# Function: Make_patch_domain_arena_DAT
# Purpose: Generate local domain and arena files for the simulation.
# Inputs: size - Landscape dimensions
#         patches - Rasterized patches
#         the_stack - Stack of the main landscape and patches
#         percent - Percent sampling for local domain
#         c_wd - Current working directory
# Outputs: local.domain and local.arena files

Make_patch_domain_arena_DAT <- function(size=sizes,patches,the_stack,percent,c_wd=cwd){
  # Generating areas
  
  patches_area <- patches %>% 
    dplyr::select(-bound_id) %>% 
    dplyr::rename(bound_id=patch_id, bound_area=patch_area)
  
  # using sspm to tesselate
  sspm_boundary_areas <- spm_as_boundary(boundaries = patches_area, 
                                         boundary = "bound_id")
  # plot(sspm_boundary_areas)
  
  # This makes sure we sample the surface of the polygons at random points 
  # 10 points per polygon
  # It also does the same process than above with a min size poof 20k
  sample_number_per_area <- patches_area %>% 
    dplyr::select(bound_id,bound_area) %>% 
    mutate(samples=as.numeric(round(bound_area/(1/(percent/100))))) %>% # sets # of samples to % of original data
    mutate(samples = replace(samples, samples<3, 3)) # Make sure each strata has 3 or more samples.
  
  sample_vector <- sample_number_per_area$samples
  names(sample_vector)=sample_number_per_area$bound_id
  
  
  voronoi_areas <- spm_discretize(sspm_boundary_areas, method = "tesselate_voronoi", 
                                  sample_surface = TRUE, nb_samples = sample_vector, min_size = 1)
  
  # Other bug in sspm: we make sure all poygons are of the type POLUYGON and NOT MULTIPOLYGON
  patches_area <- st_cast(st_make_valid(spm_patches(voronoi_areas), "POLYGON"))
  
  # plot(patches_area$geometry)
  
  points_area <- spm_points(voronoi_areas)
  
  # plot(spm_points(voronoi_areas)$geometry)
  
  
  points_area_join <- st_join(points_area,patches_area)
  
  coord_area <- as.data.frame(st_coordinates(points_area_join))
  names(coord_area) <- c("x","y")
  
  
  points_area_coord <- points_area_join %>% 
    bind_cols(coord_area) %>%
    st_drop_geometry() %>% 
    dplyr::select(patch_id,bound_id.x,patch_area,x,y)
  
  points_area_coord$bound_id.x <- substring(points_area_coord$bound_id.x, 2)
  
  points_area_coord
  
  
  # Extracting depth for local.domain
  coord_areas <- data.frame(points_area_coord$x,points_area_coord$y)
  
  coord_areas_depth <- data.frame(coord_areas,raster::extract(the_stack,coord_areas))
  
  
  
  
  
  # Create local.domain file
  local.domain <- data.frame(lat =points_area_coord$x,
                             long = points_area_coord$y,
                             rootdepth = sqrt(coord_areas_depth$Main),
                             stratum = points_area_coord$bound_id.x,
                             depth = coord_areas_depth$Main,
                             NAFO = rep("3K", nrow(points_area_coord)),
                             SFA = rep(6,nrow(points_area_coord)),
                             area = points_area_coord$patch_area
  )
  
  
  round_df <- function(x, digits) {
    # round all numeric variables
    # x: data frame 
    # digits: number of digits to round
    numeric_columns <- sapply(x, mode) == 'numeric'
    x[numeric_columns] <-  round(x[numeric_columns], digits)
    x
  }
  
  local.domain <- round_df(local.domain,4)
  
  write.table(local.domain, file = paste0(c_wd,"/","local.domain"),sep = " ", quote = F, row.names = F )
  write.table(local.domain, file = paste0(c_wd,"/","local.arena"),sep = " ", quote = F, row.names = F )
  
  
  
  
  # Generate data file PB_fall.dat
  
  # 1. Get all simulation data
  
  # Get files names
  f_list <- list.files(paste0(c_wd,"/","sim"))
  
  # Load sim data
  output <- list()
  for (i in f_list) {
    
    output[[i]] <- read_parquet(file = paste0(c_wd,"/","sim/",i))
    read_parquet(paste0("sim/",i))
  }
  
  
  
  # Write sim data
  listOfDataFrames<- list()
  for (i in f_list) {
    
    listOfDataFrames[[i]] <- data.frame(
      year=   rep(1990+as.numeric(substring(i,4)),((size-1)^2)),
      lat =   output[[i]][["x"]],
      long =  output[[i]][["y"]],
      temp =  output[[i]][["temperature"]],
      depth = output[[i]][["depth"]],
      NAFO =  rep("3K", ((size-1)^2)),
      sfa =   rep(6,((size-1)^2)),
      stratum = output[[i]][["stratum"]],
      biomass=(output[[i]][["biomass"]])
    )
  }
  df <- do.call("rbind", listOfDataFrames)
  
  df <- round_df(df,4)
  
  #write.table(df, file = "PB_fall.dat.complete",sep = " ", quote = F, row.names = F )
  write_parquet(df, paste0(c_wd,"/","PB_fall.dat.complete"))
}

# Function: Make_PB_fall.dat
# Purpose: Generate the PB_fall.dat file for the simulation.
# Inputs: percent_f - Percent sampling for PB_fall.dat
#         path - Path to the complete PB_fall.dat file
#         fname - Name of the generated PB_fall.dat file
#         c_wd - Current working directory
# Outputs: PB_fall.dat file

Make_PB_fall.dat <- function(percent_f=0.025,path="PB_fall.dat.complete",fname="PB_fall.dat",c_wd=cwd.){
  
  F_data <- arrow::read_parquet(paste0(c_wd,"/","PB_fall.dat.complete"))
  
  propotion_strata <- df %>%
    dplyr::group_by(stratum, year) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::mutate(prop = round(ifelse(n * percent_f <= 3, 3, n * percent_f)))
  
  S_data <- df %>%
    dplyr::left_join(propotion_strata, by = c("stratum", "year")) %>%
    dplyr::group_by(stratum, year) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(~ dplyr::slice_sample(.x, n = min(.x$prop[1], nrow(.x)), replace = FALSE)) %>%
    dplyr::ungroup()
  
  # S_data$biomass <- abs(rTweedie((S_data$biomass), p = 1.76, phi = 2))
  S_data$biomass <- rTweedie((S_data$biomass), p = 1.76, phi = 2)
  
  
  # S_data <- return(S_data)
  
  write.table(S_data, file = paste0(c_wd,"/",fname),sep = " ", quote = F, row.names = F )
  gc()
  
}

# Function: s_land_bio_sim_V2
# Purpose: Simulate a landscape with depth, temperature, and biomass variations.
# Inputs: cwd - Current working directory
#         size - Landscape dimensions
#         n - Number of years/simulations
#         roughness - Main nlm_mpd roughness
#         V - Variation parameter
# Outputs: Simulations + List containing the stack of the main landscape and patches

s_land_bio_sim_V2 <- function(cwd=cwd.,size=size.,n=n.,roughness=roughness.,V=V.) {
  x=y=size
  # Vairation
  v1 <- V
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
  
  Sub_L_M <- mapply(FUN = nlm_gaussianfield,
                    ncol = list_x, nrow = list_y, 
                    resolution = 1,
                    autocorr_range = 1,
                    mag_var = 5,
                    nug = 0.2,
                    mean = 0.5)
  
  
  #### 1.C Generate secondary depths which will be used to generate alternate temperatures for the different years ####
  
  for(i in 1:n){
    Sub_L_M[[i]] <- (Main_L + (Sub_L_M[[i]]/10))/2
  }
  # plot(Sub_L_M[[1]])
  
  #### 2. Create Depth patches ####
  
  message("2. Create Depth patches")
  
  # Add patches of depth variation
  
  depth_patch_variation <- nlm_randomcluster(ncol = size-1, nrow = size-1,
                                             p = 0.57,
                                             ai = c(0.6, 0.13, 0.13, 0.13))
  
  # plot(depth_patch_variation)
  # show_landscape(depth_patch_variation)
  
  # Smooth depth
  depth_patch_variation <- (focal(depth_patch_variation, w=focalWeight(depth_patch_variation,15,type = "circle"), sum, pad=T, padValue=0))
  
  # show_landscape(depth_patch_variation)
  
  depth_patch_variation <- reclassify(depth_patch_variation, cbind(0, 0.3, 0))
  
  # show_landscape(depth_patch_variation)
  
  
  for(i in 1:n){
    Sub_L_M[[i]] <- (exp(depth_patch_variation*v1))*Sub_L_M[[i]]@data@values
  }
  
  
  # Adjusting the mean of "depth" in Sub_L_M
  desired_mean_depth <- 200
  Sub_L_M_df_list <- list()
  for(i in 1:n){
    print(paste0("Sub_L_M","_df_",i))
    Sub_L_M[[i]]@data@values <- (Sub_L_M[[i]]@data@values*1126)+58
    temp_df <- as.data.frame(Sub_L_M[[i]], xy = TRUE)
    
    # Rename the 3rd column
    names(temp_df)[3] <- "depth"
    original_mean_depth <- mean(temp_df$depth)
    
    temp_df$depth <- temp_df$depth + desired_mean_depth - original_mean_depth
    
    Sub_L_M_df_list[[i]] <- temp_df
  }
  
  
  #### 4. Generate Temperature ####
  
  message("4. Generate Temperature")
  
  # Predict the temperature from depth
  for(i in 1:n){
    
    Sub_L_M_df_list[[i]]$temperature <- predict(gam_depth_sim, newdata = Sub_L_M_df_list[[i]], se.fit = T)$fit
    
  }
  
  
  #### 5. Generate Biomass ####
  
  message("5. Generate Biomass")
  
  
  # Biomass parameters #
  depth_sd = 200
  temp_sd  = 2
  scale_depth <- dnorm(0,0,depth_sd)
  scale_temp  <- dnorm(0,0,2)
  
  for(i in 1:n){
    
    Sub_L_M_df_list[[i]]$biomass <- dnorm((Sub_L_M_df_list[[i]]$depth - 312.5), 0, 100)/dnorm(0,0,100)*dnorm((Sub_L_M_df_list[[i]]$temperature - 2.916), 0, 2)/dnorm(0,0,2)
    
  }
  
  
  #### 7. Generate stratums
  # Stack the rasters and turn into df values (alignment  )
  patches_list <- make_patches(patch=Main_L)
  
  
  # Stack the rasters and turn into df values (alignment  )
  the_stack <- stack(Main_L, patches_list$patches_raster)
  names(the_stack) <- c("Main", "Patches")
  
  
  stratum = values(patches_list$patches_raster)
  
  add_stratum_column <- function(df, stratum) {
    df$stratum <- stratum
    return(df)
  }
  
  for(i in 1:n){
    
    Sub_L_M_df_list[[i]] <- add_stratum_column(Sub_L_M_df_list[[i]], stratum)
    
  }
  
  #### 8. Assemble and store the data ####
  
  #Get all the necessary data into a single list
  message("8. Assemble and store the data")
  for (i in 1:n) {
    year = 1990+i
    Sub_L_M_df_list[[i]]$year <- year
  }
  
  #write the files individually
  dir.create(paste0(cwd,"/","sim"))
  for (i in 1:n) {
    write_parquet(Sub_L_M_df_list[[i]],paste0(cwd,"/","sim/sim",i))
  }
  
  stack_patch <- list(the_stack=the_stack,patches_list=patches_list)
  return(stack_patch)
}

#### Load model ####
gam_depth_sim <- readRDS("gam_depth_sim.rds")


#### Set Simulation Parameters ####

# Working Directory
cwd. = getwd()

# Landscape dimensions
size. <- 300

# Set Main nlm_mpd roughness
roughness.=0.6

# Number of years/simulations
n. <- 20

# Vairation
V. <- 1

# Number of Replicates
reps=5

# Percent Sampling
percent. = 0.025

# Seed
seed = sample((1:50000),1)

# #### 1. Run the sim ####
# results <- s_land_bio_sim_V2() 
# 
# # Save size of each strata
# patches=results$patches_list$patches
# patches=st_set_geometry(patches,NULL)
# write_parquet(patches,"patches")
# 
# #### 2. Write the ogmap files ####
# Make_patch_domain_arena_DAT(size=size.,patches=results$patches_list$patches,the_stack=results$the_stack,percent=percent.,c_wd = cwd.)
# 
# #### 3. Slice data file ####
# Make_PB_fall.dat()
# 
# #### 4. Return to Original WD ####
# setwd(cwd)
# gc()


#### 1. Create and Start Cluster ####


#create the cluster
# n.cores <- parallelly::availableCores()/2   
# For windows
n.cores <- as.numeric(Sys.getenv('OMP_NUM_THREADS'))
main.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

#check cluster definition (optional)
print(main.cluster)


#register it to be used by %dopar%
doParallel::registerDoParallel(cl = main.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

setwd(cwd.)
print(cwd.)

foreach(
  rep = 1:reps,
  .packages = c('mgcv','dplyr','purrr','NLMR','arrow','sspm','raster','foreach','doParallel','parallelly','readr','fasterize')
) %dopar% {
  print(paste("Replicate #",rep))
  seeds = seed - 1 + rep
  set.seed(seeds)
  cwd <- getwd()          # CURRENT dir
  newdir <- paste("Run",rep,"Size",size.,"seed",seeds,"nsim",n.,"Percent",percent.,Sys.Date(),sep = "_")
  dir.create(newdir)     # Create new directory
  setwd(newdir) 
  write.table(as.data.frame(newdir),"seed")
  #### 1. Run the sim ####
  results <- s_land_bio_sim_V2()
  
  # Save size of each strata
  patches=results$patches_list$patches
  patches=st_set_geometry(patches,NULL)
  write_parquet(patches,"patches")
  
  #### 2. Write the ogmap files ####
  Make_patch_domain_arena_DAT(size=size.,patches=results$patches_list$patches,the_stack=results$the_stack,percent=percent.,c_wd = newdir)
  
  #### 6. Return to Original WD ####
  setwd("..")
  gc()
  
  
  
}


#### Kill the cluster
parallel::stopCluster(cl = main.cluster)

