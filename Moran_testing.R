
  
  #### Load simulations ####
  file_list_sim <- list.files("sim/", full.names=TRUE)
  
  #### Create raster from existing simulation data ####
  Moran_list <- list()
  counter=0
  for (fname in file_list_sim) {
    Global_Moran <- read_parquet(fname) %>% 
      dplyr::select(coord.x,coord.y,biomass) %>% 
      rasterFromXYZ()
    
      Queen=Moran(Global_Moran)
      
      Rook=Moran(Global_Moran,matrix(c(0,1,0,1,0,1,0,1,0),nrow = 3))
      
      
    counter = counter + 1
    
    #### Calculate Global Moran I ####
    Moran_list[[paste0("sim",counter)]] <- c(Queen=Queen,Rook=Rook)
  }

    
    
  
  
  ################################################################################
  # Testing stuff
  # library(spdep)
  # if (require(rgdal, quietly=TRUE)) {
  #   eire <- readOGR(system.file("shapes/eire.shp", package="spData")[1])
  # } else {
  #   require(maptools, quietly=TRUE)
  #   eire <- readShapeSpatial(system.file("shapes/eire.shp", package="spData")[1])
  # }
  # row.names(eire) <- as.character(eire$names)
  # proj4string(eire) <- CRS("+proj=utm +zone=30 +ellps=airy +units=km")
  # 
  # class(eire)
  # 
  # names(eire)
  # 
  # data("mafragh")
  # 
  # nb.patch <- spdep::poly2nb(results$patches_list$patches$geometry)
  # 
  # s.Spatial(as_Spatial(results$patches_list$patches$geometry), nb = nb.patch, plabel.cex = 0, pnb.edge.col = 'red')
  # 
  # one_year_sample_xy <- read_table("PB_fall.dat") %>% 
  #   dplyr::filter(year==1991) %>% 
  #   select(x=long,y=lat)
  # 
  # rownames(one_year_sample_xy) <- NULL
  # 
  # coords <- as.matrix(one_year_sample_xy)
  # 
  # 
  # nbtri <- tri2nb(coords)
  # 
  # s.label(one_year_sample_xy, nb = nbtri, pnb.edge.col = "red", main = "Delaunay")
  # 
  # 
  # # Defining spatial weight matrix
  # 
  # nbgab <- graph2nb(gabrielneigh(coords), sym = TRUE)
  # 
  # nb2listw(nbgab)
  # 
  # # More sophisticated forms of spatial weighting matrices can be defined.
  # distgab <- nbdists(nbgab, coords)
  # 
  # #  spatial weights are defined as a function of distance 
  # fdist <- lapply(distgab, function(x) 1 - x/max(dist(coords)))
  # 
  # # spatial weighting matrix is then created
  # listwgab <- nb2listw(nbgab, glist = fdist)
  # listwgab
  # 
  # mem.gab <- mem(listwgab)
  # mem.gab
  # 
  # 
  # barplot(attr(mem.gab, "values"), 
  #         main = "Eigenvalues of the spatial weighting matrix", cex.main = 0.7)
  # 
  # MC.env <- moran.randtest(one_year_sample$biomass, listwgab, nrepet = 999)
  # MC.env
  # 
  # r <- raster(nrows=10, ncols=10)
  # values(r) <- 1:ncell(r)
  # 
  # biomass_moran <- read_parquet("sim/sim1")
  # 
  # raster_biomass_moran <- Main_L_copy
  # 
  # values(raster_biomass_moran) <- biomass_moran$biomass
  # 
  # Moran(raster_biomass_moran)