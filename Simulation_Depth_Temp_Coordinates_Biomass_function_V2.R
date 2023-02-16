S_land_bio_sim <- function(c_wd=cwd,n,size,roughness=0.6,variation=1.5){
  y=size
  x=size

  
  #### 1. Load Packages ####
  # Package names
  packages <- c("mgcv", "plyr", "dplyr", "raster","landscapetools","devtools","openxlsx","arrow","sspm")
  
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
  # # install.packages("devtools")                    ~If needed~
  # devtools::install_github("ropensci/NLMR")
  
  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
  
  
  #### 2. Generate Main Landscape using Gaussian Field ####
  
  Main_L <- nlm_mpd(ncol = size, nrow = size, roughness = roughness)*nlm_mpd(ncol = size, nrow = size, roughness = roughness)*nlm_mpd(ncol = size, nrow = size, roughness = roughness)
  
  message("Main_landscape_generation")
  #### 2.B Generate secondary landscapes which will be used to vary temperature.####
  # Sub_L <- list()
  # for(i in 1:n){
  #   
  #   
  #   Sub_L[[i]] <- nlm_gaussianfield(x,
  #                                   y,
                                    # resolution = 10,
                                    # autocorr_range = 1,
                                    # mag_var = 5,
                                    # nug = 0.2,
                                    # mean = 0.5)
  # }
  # 
  # x <- 1000
  # y <- 1000
  
  list_x <- rep(x-1, n)
  list_y <- rep(y-1, n)
  
  Sub_L <- mapply(FUN = nlm_gaussianfield,
                   ncol = list_x, nrow = list_y, 
                  resolution = 1,
                  autocorr_range = 1,
                  mag_var = 5,
                  nug = 0.2,
                  mean = 0.5)
  
  message("2.B Generate secondary landscapes which will be used to vary temperature.")
  #### 2.C Generate secondary depths which will be used to generate alternate temperatures for the different years ####
  Sub_L_M <- list()
  for(i in 1:n){
    Sub_L_M[[i]] <- (Main_L + (Sub_L[[i]]/10))/2
  }
  
  # #visualize
  # for(i in 1:20){
  #   show_landscape(Sub_L_M[[i]])
  # }
  
  message("2.C Generate secondary depths which will be used to generate alternate temperatures for the different years")
  #### 3. Generate Main Depth and temporary sub depths ####
  Main_L_copy <- Main_L #make a copy
  Main_L_copy@data@values <- (Main_L_copy@data@values*1126)+58 #make depth
  # landscapetools::show_landscape(Main_L_copy) # visualize main ladnscape
  value <- Main_L_copy@data@values # get main landscape depths
  
  writeRaster(Main_L_copy,paste0(c_wd,"/","main_L"),overwrite=TRUE)
  # Make sub landscapes depths
  for(i in 1:n){
    Sub_L_M[[i]]@data@values <- (Sub_L_M[[i]]@data@values*1126)+58
  }
  
  # #visualize
  # for(i in 1:n){
  # landscapetools::show_landscape(Sub_L_M[[i]])
  # }
  message("3. Generate Main Depth and temporary sub depths")
  
  #### 4. Generate temperature ####
  
  
  # read real data
  real <- read.csv("../../../trawl_nl.csv")
  
  # Generate gam based on real depth and temp
  gam_depth_sim <- gam(temp_bottom ~ s(sqrt(depth), bs="ad"), data = real)
  
  # extract main landscape depth
  depth <- Main_L_copy@data@values
  depth <- as.data.frame(depth)
  
  # predict main landscape temp from depth using gam model
  predict_temp_main <- predict(gam_depth_sim, newdata = depth, se.fit = T)
  
  
  # name list depth value to depth for the gam predict
  s_depth <- list()
  for(i in 1:n){
    s_depth[[i]] <- Sub_L_M[[i]]@data@values
  }
  
  #make individual dataframes for easier renaming
  for (i in 1:length(s_depth)) {
    assign(paste0("s_depth_DF",i), as.data.frame(s_depth[[i]]))
  }
  
  # Rename column names in each DF
  for (i in paste0("s_depth_DF",1:n)){
    d=get(i)
    colnames(d)[1]="depth"
    assign(i,d)
  }
  # predict sub landscape temp from depth using gam model
  temps_sub <- list()
  for (i in paste0("s_depth_DF",1:n)){
    d=get(i)
    temps_sub[[i]] <- predict(gam_depth_sim, newdata = d[1] , se.fit = T)
  }
  message("4. Generate temperature")
  
  #### 5. Generate Biomass ####
  
  # mean_biomass = dnorm(depth - depth_opt, 0, 100)/dnorm(0,0,100)*dnorm(temp - temp_opt, 0, 2)/dnorm(0,0,2)

  biomass_mean_sub <- data.frame(matrix(nrow=nrow(s_depth_DF1),ncol=n))
  
  # for (i in 1:nrow(s_depth_DF1)){
  #   for (j in 1:n){
  #     k=paste0("s_depth_DF",j)
  #     biomass_mean_sub[i,j] <- dnorm(((Sub_L_M[[j]]@data@values[i]) - 312.5), 0, 100)/dnorm(0,0,100)*dnorm(((temps_sub[[k]][["fit"]][[i]]) - 2.916), 0, 2)/dnorm(0,0,2)
  #   }}
  depth_sd = 200
  temp_sd  = 2
  scale_depth <-dnorm(0,0,depth_sd)
  scale_temp<- dnorm(0,0,2)
  for (j in 1:n){
    k=paste0("s_depth_DF",j)
    biomass_mean_sub[,j] <- (dnorm(((Sub_L_M[[j]]@data@values) - 312.5), 0, depth_sd)/scale_depth *dnorm(((temps_sub[[k]][["fit"]]) - 2.916), 0, temp_sd)/scale_temp) 
  }
  
  biomass_mean_sub
  
  message("5. Generate Biomass")
  
  #### 6. Add Variation in shrimp biomass
  
  # Generate a landscape for variation
  biomass_variation <- nlm_randomcluster(ncol = x-1, nrow = y-1,
                                      p = 0.57,
                                      ai = c(0.5, 0.25, 0.25))
  
  biomass_variation@data@values <- biomass_variation@data@values*variation
  # biomass_variation <- nlm_gaussianfield(x,
  #                                        y,
  #                                        resolution = 1,
  #                                        autocorr_range = 100*10^variation,
  #                                        mag_var = 50,
  #                                        nug = 0.2,
  #                                        mean = 0.5)
  # 
  
  # Multiply the landscape with mean biomass at every location
  for (i in 1:n) {
    biomass_mean_sub[,i] <- exp(biomass_variation@data@values)*biomass_mean_sub[,i]
    
  }
   
  message("6. Add Variation in shrimp biomass")
  #### 7. Generate stratums
  patches_list <- make_patches(patch=Main_L_copy)
  
  
  # Stack the rasters and turn into df values (alignment  )
  the_stack <- stack(Main_L_copy, patches_list$patches_raster)
  names(the_stack) <- c("Main", "Patches")
  
  message("6. Generate stratums")
  
  #### 8. Assemble and store the data ####
  
  #Get all the necessary data into a single list
  sim<-list()
  stratum = values(patches_list$patches_raster)
  for (i in 1:n) {
    year = 1990+i
    k=paste0("s_depth_DF",i)
    depth   = value  # get froim Main_L_COPY
    temp    = temps_sub[[k]][["fit"]]
    coord   = coordinates(Main_L)
    biomass = biomass_mean_sub[,i]
      #### NEEDS PATCH_RASTER CODE ABOVE ####
    name <- paste('item:',i,sep='')
    tmp <- list(depth=depth, temp=temp, coord=coord, biomass=biomass,stratum=stratum,year=year)
    sim[[name]] <- as.data.frame(tmp)
  }
  
  #write the files individually
  dir.create(paste0(c_wd,"/","sim"))
  for (i in 1:n) {
    k=paste0("item:",i)
    write_parquet(sim[[k]],paste0(c_wd,"/","sim/sim",i))
  }
  message("7. Assemble and store the data")
  return(list(the_stack=the_stack,patches_list=patches_list))
}




