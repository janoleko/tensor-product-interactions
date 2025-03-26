
library(dplyr)
library(zoo)


## data preprocessing
data <- read.delim("./data/muskox_tracking_data_winter_20200417.txt", header = TRUE)
data$ID = data$burst_id
data <- data %>% select(
  ID, step, angle, x, y, datetime, tday, month, julian
)
data$datetime <- as.POSIXct(data$datetime, format = "%Y-%m-%d %H:%M:%S")

# linearly interpolate missing locations x,y
data[c("x", "y")] <- data.frame(lapply(data[c("x", "y")], na.approx, na.rm = FALSE))

## save data
saveRDS(data, "./data/muskox_winter.rds")

muskox = readRDS("./data/muskox_winter.rds")
