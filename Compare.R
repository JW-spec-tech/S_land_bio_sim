library(dplyr)
library(mosaic)
# Simulation stats 

# Actual stats
F_data <- readr::read_table("PB_fall.dat.complete")

biomass_year <- F_data %>%
  dplyr::group_by(year) %>%
  dplyr::select(biomass,year) %>%
  dplyr::summarise(t_bio=sum(biomass),SD_bio=sd(biomass),CI_lower=quantile(biomass, prob = c(0.025)),CI_upper=quantile(biomass, prob = c(0.975)))
      # dont need CI here
#                       
# biomass_year <- F_data %>% 
#   dplyr::group_by(year) %>% 
#   dplyr::select(biomass,year) %>%
#   mean
# 
#                                                      
# quantile(F_data$biomass, prob = c(0.025, 0.975))
# 
# 
# 
# View(biomass_year)
# 
# trawl_year <- data.frame(table(s_data$year)) 
# 
# 
# 
# hist_biomass <- s_data %>% 
#   filter(year==1991) %>% 
#   dplyr::select(biomass)
# hist(hist_biomass$biomass)
# 
# 
# 
# gf_raster(lat~long, fill = ~log(biomass), data = F_data,
#           show.legend = NA) %>% 
#   gf_facet_wrap(~year)+
#   scale_fill_viridis_c()
# 
# 
# 
# D <- expand.grid(x = 0:5, y = 0:5)
# D$z <- runif(nrow(D))
# gf_tile(y ~ x, fill = ~z, data = D)
# gf_tile(z ~ x + y, data = D)


ogmap_estimates <- readr::read_table("biomass__Ogmap2 - Demo.log", skip = 2)

# gam_biomass_year <- dat_grid %>% 
#   dplyr::group_by(year) %>% 
#   dplyr::select(fit_simple_gam,fit_spatialonly_gam,year) %>% 
#   dplyr::summarise(tfit_simple_gam=sum(fit_simple_gam)/1000,tfit_spatialonly_gam=sum(fit_spatialonly_gam)/1000)



comparison <- data.frame(dat_grid_pred,
                         ogmap=(as.numeric(ogmap_estimates$Estimate[1:6])*1000), 
                         CI_upper_ogmap= ogmap_estimates$UpCIval[1:6]*1000, 
                         CI_lower_ogmap=ogmap_estimates$LowCIval[1:6]*1000)

comparison <- comparison %>% 
  mutate(Ogmap_CI = ogmap >= CI_lower_ogmap & value <= CI_upper_ogmap) %>% 
  mutate(GAM_CI = fit_simple_gam >= CI_lower_gam & value <= CI_upper_gam)

# comparison <- comparison %>% 
#   mutate(percent=tbio/ogmap)
# 
# summary(comparison$percent)
# 
# 
# years <- c(1990:(1990+50))
# year_commas <- paste(years, sep = "", collapse = ",")
# write.table(year_commas,"years.txt")





