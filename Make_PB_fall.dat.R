Make_PB_fall.dat <- function(percent=0.1,path="PB_fall.dat.complete",fname="PB_fall.dat"){
  
  F_data <- readr::read_table("PB_fall.dat.complete")
  
  # names(F_data) <- F_data[1,]
  # 
  # F_data = F_data[-1,]
  
  
  propotion_strata <- F_data %>% 
                      dplyr::group_by(stratum,year) %>% 
                      dplyr::summarise(n=n()) %>% 
                      dplyr::mutate(min=3)
  
  propotion_strata$prop = round(ifelse(propotion_strata$n*percent<=propotion_strata$min,propotion_strata$min,ifelse(propotion_strata$n*percent>=15,15,propotion_strata$n*percent)))

  propotion_strata$prop_random = propotion_strata$prop + sample(-2:2,length(propotion_strata$prop),replace = T)
    
  propotion_strata$prop_random= ifelse(propotion_strata$prop_random<propotion_strata$min,propotion_strata$min,propotion_strata$prop_random)
  
  F_data_prop <- dplyr::left_join(F_data,propotion_strata[,c("year","stratum","prop_random")])
  
  S_data <- F_data_prop %>% 
    group_by(stratum,year) %>% 
    dplyr::sample_n(size=prop_random) %>% 
   mutate(biomass=rTweedie((biomass*1000), p = 1.76, phi= 2))
    
 
# S_data <- return(S_data)
  
  write.table(S_data, file = fname,sep = " ", quote = F, row.names = F )
  
}
