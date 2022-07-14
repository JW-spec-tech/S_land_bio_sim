Make_PB_fall.dat <- function(percent_f=0.025,path="PB_fall.dat.complete",fname="PB_fall.dat"){
  
  F_data <- arrow::read_parquet("PB_fall.dat.complete")
  
  propotion_strata <- F_data %>% 
                      dplyr::group_by(stratum,year) %>% 
                      dplyr::summarise(n=n()) %>% 
                      dplyr::mutate(min=3)
  
  propotion_strata$prop = round(ifelse(propotion_strata$n*percent_f<=propotion_strata$min,propotion_strata$min,ifelse(propotion_strata$n*percent_f>=15,15,propotion_strata$n*percent_f)))

  propotion_strata$prop_random = propotion_strata$prop + sample(-2:2,length(propotion_strata$prop),replace = T)
    
  propotion_strata$prop_random= ifelse(propotion_strata$prop_random<propotion_strata$min,propotion_strata$min,propotion_strata$prop_random)
  
  hist(propotion_strata$prop_random)
  
  F_data_prop <- dplyr::left_join(F_data,propotion_strata[,c("year","stratum","prop_random")])
  
  S_data <- F_data_prop %>% 
    group_by(stratum,year)  %>% 
  group_split() %>%
    purrr::map_dfr(~slice_sample(.x,n=.x$prop_random[1])) %>% 
    bind_rows() %>% 
   mutate(biomass=rTweedie((biomass*1000), p = 1.76, phi= 2))
    
 
# S_data <- return(S_data)
  
  write.table(S_data, file = fname,sep = " ", quote = F, row.names = F )
  gc()
  
}
