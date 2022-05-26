library(dplyr)
library(mosaic)
library(arrow)
library(ggplot2)
# Simulation stats 

# Actual stats
F_data <- readr::read_table("PB_fall.dat.complete")


biomass_year <- F_data %>%
  dplyr::group_by(year) %>%
  dplyr::select(biomass,year) %>%
  dplyr::summarise(t_bio=sum(biomass)
                   )


#Load predictions
ogmap_estimates <- readr::read_table("biomass__Ogmap2 - Demo.log", skip = 2)
Predictions_summary <- read_parquet("Predictions_summary")

# gam_biomass_year <- dat_grid %>% 
#   dplyr::group_by(year) %>% 
#   dplyr::select(fit_simple_gam,fit_spatialonly_gam,year) %>% 
#   dplyr::summarise(tfit_simple_gam=sum(fit_simple_gam)/1000,tfit_spatialonly_gam=sum(fit_spatialonly_gam)/1000)


# in kg
comparison <- data.frame(biomass_year, Predictions_summary,
                         ogmap=(as.numeric(ogmap_estimates$Estimate)*1e3), 
                         CI_upper_ogmap= ogmap_estimates$UpCIval*1e3, 
                         CI_lower_ogmap=ogmap_estimates$LowCIval*1e3)
comparison <- comparison[ -c(3) ]
comparison <- comparison[c(1:6,8,7)]



comparison <- comparison %>% 
  dplyr::mutate(Ogmap_CI = t_bio >= CI_lower_ogmap & t_bio <= CI_upper_ogmap) %>% 
  dplyr::mutate(GAM_CI = t_bio >= lower & t_bio <= upper)

#### PErcentage of tim in the interval 
Ogmap_percent <- sum(comparison$Ogmap_CI, na.rm = TRUE)/nrow(comparison)*100
Gam_percent   <- sum(comparison$GAM_CI, na.rm = TRUE)/nrow(comparison)*100
Result_CI <- data.frame(
  model=c("Ogmap_percent","Gam_percent") ,  
  value=c(Ogmap_percent,Gam_percent)
)
# Barplot
plot <- ggplot(Result_CI, aes(x=model, y=value)) + 
  geom_bar(stat = "identity", fill='lightblue', color ='black')+
  geom_text(aes(label=paste0(value,'%')),  
            position = position_dodge(width = 1),
            vjust = 7)+
  ggtitle('Percentage of time the simulated biomass falls within the model CI
           20 simulations to run for parameters (predict_intervals)
           500 samples using 0.1% of the entire dataset')+
  theme(plot.title = element_text(hjust = 0.5))
plot
 
ggsave("plot_test.png",
       plot = plot,
       device = "png")


#graph Confidence Interval coverage per year
# Create ggplot object
data_Graph <- ggplot(comparison, aes(year))

# Create graph
CI_plot <- data_Graph +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha=0.25)+
  geom_line(aes(y=point_est, colour = "red"))+
  geom_ribbon(aes(ymin = CI_lower_ogmap, ymax = CI_upper_ogmap), fill = "lightblue", alpha=0.5)+
  geom_line(aes(y=ogmap, colour = "lightblue"))+
  geom_line(aes(y=t_bio, colour = "black"))+
  scale_color_manual(name= "Biomass", labels = c("Total biomass","GAM", "OGmap"),values = c("black", "red","blue"))+
  labs(title = "Coverage of Confidence Intervals
       GAM VS OGmap", x="Year (simulation #)", y="Biomass in kg")+
  theme(plot.title = element_text(hjust = 0.5))
  
ggsave("CI_plot.png",
       plot = CI_plot,
       device = "png")

# 
# dplyr::between()
# 
# between(1:12, 7, 9)

# comparison <- comparison %>% 
#   mutate(percent=tbio/ogmap)
# 
# summary(comparison$percent)
# 
# 
# years <- c(1990:(1990+50))
# year_commas <- paste(years, sep = "", collapse = ",")
# write.table(year_commas,"years.txt")


###### Graphing ######
# library(ggformula)
# library(ggpubr)
#   #### Histogram biomass ####
# 
# truth <- gf_histogram(~t_bio, data = biomass_year,
#              fill = "skyblue", 
#              color = "black"
# )
#   
# ogmap <-  gf_histogram(~Estimate, data = ogmap_estimates,
#                fill = "skyblue", 
#                color = "black"
#   )
# 
#   gf_histogram(~UpCIval, data = ogmap_estimates,
#                fill = "skyblue", 
#                color = "black"
#   )
#   
#   gf_histogram(~LowCIval, data = ogmap_estimates,
#                fill = "skyblue", 
#                color = "black"
#   )
# 
# 
# gam <-  gf_histogram(~fit_simple_gam, data = Predictions_summary,
#                fill = "skyblue", 
#                color = "black"
#   )
#   
#   gf_histogram(~CI_upper, data = Predictions_summary,
#                fill = "skyblue", 
#                color = "black"
#   )
#   
#   gf_histogram(~CI_lower, data = Predictions_summary,
#                fill = "skyblue", 
#                color = "black"
#   )
# 
# ggarrange(ogmap,truth,gam,ncol = 3)
#   
# #### Does the simulation fit the models interval? ####
#   Interval_gam <- list()
#   for (i in 1:nrow(Predictions_summary)) {
#     Interval_gam[[i]] <- dplyr::between(biomass_year$t_bio[i], Predictions_summary$CI_lower[i], Predictions_summary$CI_upper[i])
#   }
#   Interval_ogmap <- list()
#   for (i in 1:nrow(Predictions_summary)) {
#     Interval_ogmap[[i]] <- dplyr::between(biomass_year$t_bio[i]*1000, ogmap_estimates$LowCIval[i]*1e6, ogmap_estimates$UpCIval[i]*1e6)
#   }
# 
#   Compare_CI <- as.data.frame(Interval_gam)
# plot(biomass_year$t_bio[1:40],Predictions_summary$fit_simple_gam)
# plot(biomass_year$t_bio[1:40],ogmap_estimates$Estimate[1:40])
