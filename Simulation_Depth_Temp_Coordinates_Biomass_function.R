S_land_bio_sim <- function(n,x,y,res,auto,var,nug,mean){
  # x=100
  # y=100
  #### 1. Load Packages ####
  # Package names
  packages <- c("NLMR", "mgcv", "plyr", "dplyr", "raster","landscapetools","devtools","openxlsx","arrow")
  
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
  
  Main_L <- nlm_gaussianfield(x,
                              y,
                              resolution = res,
                              autocorr_range = auto,
                              mag_var = var,
                              nug = nug,
                              mean = mean)
  
  
  #### 2.B Generate secondary landscapes which will be used to vary temperature.####
  Sub_L <- list()
  for(i in 1:n){
    Sub_L[[i]] <- nlm_gaussianfield(x,
                                    y,
                                    resolution = 10,
                                    autocorr_range = 1,
                                    mag_var = 5,
                                    nug = 0.2,
                                    mean = 0.5)
  }
  
  #### 2.C Generate secondary depths which will be used to generate alternate temperatures for the different years ####
  Sub_L_M <- list()
  for(i in 1:n){
    Sub_L_M[[i]] <- (Main_L + Sub_L[[i]])/2
  }
  
  # #visualize
  # for(i in 1:20){
  #   show_landscape(Sub_L_M[[i]])
  # }
  #### 3. Generate Main Depth and temporary sub depths ####
  Main_L_copy <- Main_L #make a copy
  Main_L_copy@data@values <- (Main_L_copy@data@values*1126)+58 #make depth
  landscapetools::show_landscape(Main_L_copy) # visualize main ladnscape
  value <- Main_L_copy@data@values # get main landscape depths
  
  
  # Make sub landscapes depths
  for(i in 1:n){
    Sub_L_M[[i]]@data@values <- (Sub_L_M[[i]]@data@values*1126)+58
  }
  
  # #visualize
  # for(i in 1:n){
  # landscapetools::show_landscape(Sub_L_M[[i]])
  # }
  
  #### 4. Generate temperature ####
  
  
  # read real data
  real <- read.csv("trawl_nl.csv")
  
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
  
  #### 5. Generate Biomass ####
  
  # mean_biomass = dnorm(depth - depth_opt, 0, 100)/dnorm(0,0,100)*dnorm(temp - temp_opt, 0, 2)/dnorm(0,0,2)
  biomass_mean_sub <<- data.frame(matrix(nrow=nrow(s_depth_DF1),ncol=n))
  
  for (i in 1:nrow(s_depth_DF1)){
    for (j in 1:n){
      k=paste0("s_depth_DF",j)
      biomass_mean_sub[i,j] <- dnorm(((Sub_L_M[[j]]@data@values[i]) - 312.5), 0, 100)/dnorm(0,0,100)*dnorm(((temps_sub[[k]][["fit"]][[i]]) - 2.916), 0, 2)/dnorm(0,0,2)
    }}
  biomass_mean_sub
  #### 6. Assemble and store the data ####
  
  #Get all the necessary data into a single list
  sim<-list()
  for (i in 1:n) {
    k=paste0("s_depth_DF",i)
      depth   = Sub_L_M[[i]]@data@values
      temp    = temps_sub[[k]][["fit"]]
      coord   = coordinates(Main_L)
      biomass = biomass_mean_sub[,i]
      name <- paste('item:',i,sep='')
      tmp <- list(depth=depth, temp=temp, coord=coord, biomass=biomass)
      sim[[name]] <- as.data.frame(tmp)
  }
  
  #write the files individually
  for (i in 1:n) {
    k=paste0("item:",i)
    write_parquet(sim[[k]],paste0("sim",i))
  }

  
}


biomass <- S_land_bio_sim(25,1000,1000,1,50,5,0.2,0.5)
