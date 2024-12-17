### Code that loads and reformats temperature data
##  Created by: Nicholas Bussberg
##  Date last updated: 2024-12-16
##  Filename: irfkrg7a-parse-data.R
##  Purpose: loads global temperature data and reformats it to be used in 
##  study

## Data downloaded from: https://www.nsstc.uah.edu/data/msu/t2lt/
## For Feb 1979 and Jan 1982 (not used in the present study/code), there were 
## several instances where there was no break between neighboring values. 
## A space should be manually added so that the following code would work. 
## README for data: https://www.nsstc.uah.edu/data/msu/docs/readme.msu 

### Setup ###
rm(list=ls())
options(width = 80)

filename <-"irfkrg7a-parse-data"

sink(paste(filename, ".log", sep = "")) # create log file

# R version
R.version.string


### 1 Load functions

source(file = "irf-fcn03-distance.R")
source(file = "irf-fcn09-parse-data.R")


### 2 Create matrix to store results

# convert lat/long degrees to radians
latitude <- (seq(-88.75,88.75,2.5))*(pi/180)
latitude <- latitude + pi/2 # switch range to 0 to pi

longitude <- (seq(-178.75,178.75,2.5))*(pi/180)

num.cols <- length(latitude)    # should be 72
num.rows <- length(longitude)   # should be 144


# create matrix of locations
Lmat <- matrix(c(kronecker(latitude,rep(1,num.rows)),
                    rep(longitude,num.cols)), num.cols * num.rows, 2) 
# In data, each lat value is repeated 144 times in order with corresponding long

# Store Lmat
saveRDS(Lmat, file = paste("data/",filename,"-Lmat.Rds", sep = ""))



## 3 Load data from directory, parse the data into new format, save as txt file

# create matrix to store temp anomalies 
# 1 col for each month; each row is a year
loc.data <- cbind(Lmat[,1],Lmat[,2], NA, NA, NA, NA, NA, NA,
                  NA, NA, NA, NA, NA, NA) 
colnames(loc.data) <- c("latitude", "longitude", "01", "02", "03", "04",
                        "05", "06", "07", "08", "09", "10", "11", "12")


# create vector of years of where there is data
# For this study, only 2015 was needed
# Note: datasets must start in January (month=01).
years <- seq(2015, 2015, 1)

for (i in 1:length(years)){
  data.name <- paste("data/tltmonamg.", years[i], "_5.6", sep = "")
  temp.data <- parse.data(data.name, loc.data)
  saveRDS(temp.data, 
              file = paste("data/",filename,"-",years[i], ".Rds", sep = ""))
}



### 4 Compute distance vector (Dvec)

P <- length(Lmat[,1])

Dvec <- Dfun(Lmat, P)

saveRDS(Dvec, file=paste("data/", filename, "-Dvec.Rds", sep = ""))



# close log
sink()
