
Make_patch_domain_arena_DAT <- function(patches,the_stack){
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
    mutate(samples=as.numeric(round(bound_area/60)))
  
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
  
  write.table(local.domain, file = "local.domain",sep = " ", quote = F, row.names = F )
  write.table(local.domain, file = "local.arena",sep = " ", quote = F, row.names = F )
  
  
  
  
  # Generate data file PB_fall.dat
  
  # 1. Get all simulation data
  
  # Get files names
  f_list <- list.files("sim/")
  
  # Load sim data
  output <- list()
  for (i in f_list) {
    
    output[[i]] <- read_parquet(file = paste0("sim/",i))
  }
  
  
  # Write sim data
  listOfDataFrames<- list()
  for (i in f_list) {
    
    listOfDataFrames[[i]] <- data.frame(
      year=   rep(1990+as.numeric(substring(i,4)),250000),
      lat =   output[[i]][["coord.x"]],
      long =  output[[i]][["coord.y"]],
      depth = output[[i]][["depth"]],
      NAFO =  rep("3K", 250000),
      sfa =   rep(6,250000),
      stratum = output[[i]][["stratum"]],
      biomass=(output[[i]][["biomass"]])
    )
  }
  df <- do.call("rbind", listOfDataFrames)
  
  df <- round_df(df,4)
  
  write.table(df, file = "PB_fall.dat.complete",sep = " ", quote = F, row.names = F )
}
