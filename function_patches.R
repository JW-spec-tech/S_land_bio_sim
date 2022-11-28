make_patches <- function(patch){ #,plot=F
   
  # patch <- raster("main_L.gri")
  
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


