Make_PB_fall.dat <- function(percent=0.1,path="PB_fall.dat.complete",fname="PB_fall.dat"){
  
  F_data <- readr::read_table("PB_fall.dat.complete")
  
  # names(F_data) <- F_data[1,]
  # 
  # F_data = F_data[-1,]

  S_data <- F_data %>%
         group_by(stratum,year) %>%
         slice_sample(prop = percent/100) %>% 
         mutate(biomass=rTweedie((biomass*1000), p = 1.76, phi= 2))
 
# S_data <- return(S_data)
  
  write.table(S_data, file = fname,sep = " ", quote = F, row.names = F )
  
}
