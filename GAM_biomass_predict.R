library(mgcv)
library(tidyr)
library(sspm)
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
dat_grid <- as.data.frame(expand_grid(long=seq(0,1,length=50), lat= seq(0,1,length=50),year=c(1991:2000)))%>%
  dplyr::mutate(fit_simple_gam  = predict.gam(simple_gam,type = "response",newdata = .),
                fit_spatialonly_gam = predict.gam(spatialonly_gam,type = "response",newdata = .))




library(mgcv)
library(readr)
library(sspm)
library(tidyr)

trawl_data <- read_table("PB_fall.dat")
head(trawl_data)

simple_gam <- gam(biomass ~
                    te(long, lat, year, bs= c("tp", "tp"), d = c(2,1)),
                  data = trawl_data, method="REML", family = "tw")

predict(simple_gam)

# spatialonly_gam <- gam(biomass ~s(long, lat, bs="tp"),
#                        data = trawl_data, method="REML", family = "tw")

dat_grid <- as.data.frame(expand_grid(long=seq(0,1,length=10),
                                      lat= seq(0,1,length=10),
                                      year=c(1991:1996)))




# 
# predict_intervals <- sspm:::predict_productivity_intervals
# produce_sims <- function (fit, new_data, n = 100) 
# {
#   checkmate::assert_class(fit, "gam")
#   coefs <- stats::coef(fit)
#   lp <- predict(fit, newdata = new_data, type = "lpmatrix")
#   vcv <- stats::vcov(fit)
#   coefs_sim <- t(rmvn(n = n, coefs, vcv))
#   sims <- lp %*% coefs_sim
#   return(sims)
# }
# 
# sspm:::pred
# find_quantiles <- sspm:::find_quantiles
# confidence_interval <- sspm:::confidence_interval
# prediction_interval <- sspm:::prediction_interval

# predict(simple_gam, newdata = dat_grid, type = "lpmatrix")

dat_grid_pred <- dat_grid %>%
  dplyr::mutate(fit_simple_gam = predict.gam(simple_gam,type = "response",
                                             newdata = .)) %>%
  dplyr::bind_cols(predict_intervals(object_fit = simple_gam, new_data = ., PI= F, n=100)) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(fit_simple_gam = sum(fit_simple_gam),
                   CI_upper_gam = sum(CI_upper), CI_lower_gam = sum(CI_lower)) %>%
  dplyr::ungroup()

