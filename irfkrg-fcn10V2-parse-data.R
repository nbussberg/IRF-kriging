### Function to parse temperature anomaly data ###
##  Created by: Nicholas Bussberg
##  Date last updated: 2021-12-16
##  Filename: irf-fcn10V2-parse-data.R
##  Purpose: create a function that parses temperature anomaly data that are 
##      from NSSTC (https://www.nsstc.uah.edu/data/msu/t2lt/)


parse.data <- function(file, loc.data){
    # Parses temperature anomaly data (file) by month. 
    # Each month's data is put in a separate column with corresponding location
    #   data (loc.data).
    # The function takes a filename (file) in the form of a string and a matrix
    #   with location information and 12 empty columns starting at column 
    #   index 3 (loc.data)
    # The function output is an updated location matrix that includes the 
    #   temperature data. 
    
## 1 Import raw data
temp.data <- read.table(file, fill=T)


## 2 
for (i in 1:12){
    startvec <- (i-1)*649 + 2
    endvec   <- i*649
    temp.vec <- as.vector(t(temp.data[startvec:endvec,]))
    temp.vec <- as.numeric(temp.vec) / 100 # data are in 100ths of degrees
    
    loc.data[, i+2] <- temp.vec
}

return(loc.data)
}





