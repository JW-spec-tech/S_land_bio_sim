# Strap calculation to estimate biomass from survey data for comparison with ogmap/gam

#### Packages ####
library(dplyr)
library(readr)

#### 1. Load data ####
survey_raw_data <- readr::read_table("PB_fall.dat")

strata_area <- patches %>% 
  st_set_geometry(.,NULL) %>% 
  dplyr::mutate(stratum=as.numeric(patch_id)) %>% 
  units::drop_units() %>% 
  dplyr::select(stratum,patch_area) 

#### 2. Merge the area and survey data ####

Survey_W_Area <- dplyr::left_join(survey_raw_data,strata_area)

#### 2. Run strap calculation ####
Strap_estimate <- Survey_W_Area %>%
  
  dplyr::group_by(year,stratum)%>%
  
  dplyr::mutate(biomass=biomass/1000) %>% 

  dplyr::summarize(Bj = patch_area*mean(biomass),
            s2j = (patch_area^2)*var(biomass)/(n())) %>% 
  
  dplyr::group_by(year)%>%
  
  dplyr::summarize(B_total = sum(Bj),B_se= sqrt(sum(s2j))) %>%
  
  dplyr::mutate(lower = B_total - 1.96*B_se,
         upper = B_total + 1.96*B_se)
