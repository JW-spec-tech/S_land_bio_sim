# Strap calculation to estimate biomass from survey data for comparison with ogmap/gam

#### Packages ####
library(dplyr)
library(readr)

#### 1. Load data ####
survey_raw_data <- readr::read_table("Data/Run_1_Size_100_seed_1_n_sim_100_Percent_2_2022-06-16/PB_fall.dat")

strata_area <- patches_area %>% 
  select()

#### 2. Run strap calculation ####
Strap_estimate <- survey_raw_data %>%
  
  group_by(year,stratum)%>%

  summarize(Bj = area[1]*mean(biomass),
            s2j = area[1]^2*var(biomass)/n()) %>% 
  
  group_by(year)%>%
  
  summarize(B_total = sum(Bj),B_se= sqrt(sum(s2j))) %>%
  
  mutate(lower = B_total - 1.96*B_se,
         upper = B_total + 1.96*N_se)