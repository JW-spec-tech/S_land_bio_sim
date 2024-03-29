Analyse_Gam <- function(file_name='Predictions_summary') {
  

library(arrow)
memory.limit(100000)
#### 1. Load sim names
file_list <- list.files("sim/", full.names=TRUE)

#### 2. Merge sims
allData <- plyr::ldply(as.list(file_list), arrow::read_parquet)
# Do a Dirty test
allData$year <- rep(1991:2010, each = 249001)

#### 3. Turn data into a data.table + Write a copy
# Predictions <- read_parquet("Predictions")
Predictions <- data.table::as.data.table(allData)
write_parquet(Predictions,'Predictions')




#### 4. Calculate mean on each column + Write a copy
Predictions_summary <- Predictions[, lapply(.SD,sum), by=.(year)]
write_parquet(Predictions_summary,file_name)
}



