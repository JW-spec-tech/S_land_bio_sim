
# Gam_prediction <-function(path="PB.fall.dat",size=500)

# 1. Load the biomass sample data & generate new prediction grid
# x=size
# y=size
sample <- readr::read_table("PB_fall.dat")
# grid <- expand_grid(x=seq(0,1,length=500), y= seq(0,1,length=500))


# 2. Run the models
simple_gam <- gam(biomass ~ te(long, lat, year, bs= c("tp", "tp"), d = c(2,1)), data= sample, method="REML", family = "tw")

spatialonly_gam <- gam(biomass ~s(long, lat, bs="tp"), data= sample, method="REML", family = "tw")


# 3. Run the predictions
dat_grid <- as.data.frame(expand_grid(long=seq(0,1,length=500), lat= seq(0,1,length=500),year=c(1991:2000)))%>%
  dplyr::mutate(fit_simple_gam  = predict.gam(simple_gam,type = "response",newdata = .),
                fit_spatialonly_gam = predict.gam(spatialonly_gam,type = "response",newdata = .))
