library(arrow)
library(NLMR)
# data <- read_parquet("sim1")

# data$biomass = data$biomass*500
# data$trawl_sample <- mgcv::rTweedie(data$biomass,phi=2)

Main_L <- nlm_gaussianfield(300,
                            300,
                            resolution = 10,
                            autocorr_range = 50,
                            mag_var = 5,
                            nug = 0.2,
                            mean = 0.5)
Main_L_copy <- Main_L #make a copy
Main_L_copy@data@values <- (Main_L_copy@data@values*1126)+58

landscapetools::show_landscape(Main_L_copy)

# Smooth the raster

Main_L_copy_smt <- focal(Main_L_copy, w=matrix(1,9,9), fun=mean)

# Reclass thr raster
library(raster)
# ?reclassify
brks <- 16
rcl_matrix <- matrix(c(seq(0, 1200, length.out = brks),
                       c(seq(0, 1200, length.out = brks)[-1], 1400),
                       1:brks),
                       ncol=3, nrow=brks)
Main_L_copy_rcl <- reclassify(Main_L_copy_smt, 
                              rcl_matrix)
# plot(Main_L_copy_rcl)
Main_L_copy_rcl_clumped <- rasterToPolygons(Main_L_copy_rcl, dissolve = TRUE)
plot(Main_L_copy_rcl_clumped)











HL <- matrix(c(1,1,2,3,3,2,1, 1,2,4,2,1,0,1), ncol=2, byrow=FALSE)
require(rgdal)
require(sp)
L = Line(HL)
Ls = Lines(list(L), ID = "a")
SL = SpatialLines(list(Ls))
proj4string(SL) = CRS("+proj=utm +zone=46 +datum=WGS84 +units=m +no_defs")
plot(SL)
vp<-seq(0, 4, by= 1)
hp<-seq(0, 4, by= 1)
vpXmin<-cbind(rep(min(hp),length(vp)), vp)
vpXmax<-cbind(rep(max(hp),length(vp)), vp)
l <- lapply(1:nrow(vpXmin),function(i) rbind(vpXmax[i,],vpXmin[i,]))

