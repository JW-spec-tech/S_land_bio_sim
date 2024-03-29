# Strap calculation to estimate biomass from survey data for comparison with ogmap/gam
STRAP <- function(fname="Strap_estimate",c_wd=cwd) {
  

#### Packages ####
library(dplyr)
library(readr)
library(ggplot2)
  library(arrow)

#### 1. Load data ####
  F_data <- read_parquet(paste0(c_wd,"/PB_fall.dat.complete"))
  # write.table(F_data, "test.table")
  # read_parquet("PB_fall.dat.complete")
  # readr::read_table("test.table")
  
  
  biomass_year <- F_data %>%
    dplyr::group_by(year) %>%
    dplyr::select(biomass,year) %>%
    dplyr::summarise(t_bio=sum(biomass)
    )
  
survey_raw_data <- readr::read_table(paste0(c_wd,"/","PB_fall.dat"))
patches <- read_parquet(paste0(c_wd,"/","patches"))


strata_area <- patches %>% 
  dplyr::mutate(stratum=as.numeric(patch_id)) %>% 
  units::drop_units() %>% 
  dplyr::select(stratum,patch_area) 

#### 2. Merge the area and survey data ####

Survey_W_Area <- dplyr::left_join(survey_raw_data,strata_area)

#### 2. Run strap calculation ####
Strap_estimate <- Survey_W_Area %>%
  dplyr::group_by(year, stratum) %>%
  dplyr::mutate(biomass = biomass / 1000) %>%
  dplyr::filter(!is.na(patch_area) & !is.na(biomass)) %>% # filter out missing values
  dplyr::summarize(Bj = patch_area * mean(biomass),
                   s2j = (patch_area^2) * var(biomass) / (n())) %>%
  distinct(year, stratum, .keep_all = TRUE) %>%
  dplyr::group_by(year) %>%
  dplyr::summarize(B_total = sum(Bj), B_se = sqrt(sum(s2j))) %>%
  dplyr::mutate(lower = B_total - 1.96*B_se,
                upper = B_total + 1.96*B_se)



#### 4. Save the data ####

# join the data 
Strap_estimate <- dplyr::left_join(Strap_estimate,biomass_year)
# write the data
write_parquet(Strap_estimate,paste0(c_wd,"/",fname))

}

# 
# data_Graph <- ggplot(Strap_estimate, aes(year))
# 
# CI_plot <- data_Graph +
#   geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha=0.25)+
#   geom_line(aes(y=B_total, colour = "blue"))+
#   geom_line(aes(y=t_bio, colour = "black"))+
#   scale_color_manual(name= "Biomass", labels = c("STRAP","Real"),values = c("blue","black"))+
#   labs(title = "Coverage of Confidence Intervals
#        GAM VS OGmap", x="Year (simulation #)", y="Biomass in kg")+
#   theme(plot.title = element_text(hjust = 0.5))
# 
# print(CI_plot)

