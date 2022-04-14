library(arrow)
memory.limit(100000)
#### 1. Load sim names
file_list <- list.files("Result/", full.names=TRUE)

#### 2. Merge sims
allData <- plyr::ldply(as.list(file_list), arrow::read_parquet)

#### 3. Turn data into a data.table + Write a copy
# Predictions <- read_parquet("Predictions")
Predictions <- data.table::as.data.table(allData)
write_parquet(Predictions,'Predictions')

#### 4. Calculate mean on each column + Write a copy
Predictions_summary <- Predictions[, lapply(.SD,sum), by=.(year)]
write_parquet(Predictions_summary,'Predictions_summary')
