### Code that simulates IRF data for kriging and splits each
##  Created by: Nicholas Bussberg
##  Date last updated: 2022-12-13
##  Filename: irfkrg3V3-kriging-split.R
##  Purpose: simulates IRF2 and IRF3 data and splits each into training and 
    ## testing datasets to be used in kriging

### Setup ###
rm(list=ls())
options(width = 80)

filename <-"irfkrg3V3-kriging-split"

sink(paste(filename, ".log", sep = "")) # create log file

# R version
R.version.string


### 00 Load functions

source(file = "irfkrg-fcn01V3-sphere-harmonics.R")


###  01 Load data

Lmat <- readRDS("Data/irfkrg1V3-irf-data-Lmat.Rds")

ud.k2 <- readRDS("Data/irfkrg1V3-irf-data-ud-k2.Rds")
ud.k3 <- readRDS("Data/irfkrg1V3-irf-data-ud-k3.Rds")


P <- length(ud.k2[1,])



###  02 Simulate X.0 (data)

set.seed(2345)

X.0.k2 <- ud.k2 %*% matrix(rnorm(P,0,1),P,1)
X.0.k3 <- ud.k3 %*% matrix(rnorm(P,0,1),P,1)


# Combine X.0 and Lmat to form one data frame

X.Lmat.k2 <- data.frame(cbind(X.0.k2, Lmat))
colnames(X.Lmat.k2) <- c("X.0.k2", "Lmat1-lat", "Lmat2-lon")

X.Lmat.k3 <- data.frame(cbind(X.0.k3, Lmat))
colnames(X.Lmat.k3) <- c("X.0.k3", "Lmat1-lat", "Lmat2-lon")



###  03 Split data into train and test sets

set.seed(135)

train_seq.k2 <- sample(seq_len(nrow(X.Lmat.k2)), size=.9*length(X.Lmat.k2[,1]))
train_seq.k3 <- sample(seq_len(nrow(X.Lmat.k3)), size=.9*length(X.Lmat.k3[,1]))

data.train.k2 <- X.Lmat.k2[train_seq.k2, ]
data.train.k3 <- X.Lmat.k3[train_seq.k3, ]

data.test.k2  <- X.Lmat.k2[-train_seq.k2, ]
data.test.k3  <- X.Lmat.k3[-train_seq.k3, ]



###  04 Save datasets to files

saveRDS(X.Lmat.k2, file = paste("Data/", filename, "-k2-full.Rds", sep=""))
saveRDS(X.Lmat.k3, file = paste("Data/", filename, "-k3-full.Rds", sep=""))

saveRDS(data.train.k2, file = paste("Data/", filename, "-k2-train.Rds", sep=""))
saveRDS(data.train.k3, file = paste("Data/", filename, "-k3-train.Rds", sep=""))

saveRDS(data.test.k2, file = paste("Data/", filename, "-k2-test.Rds", sep=""))
saveRDS(data.test.k3, file = paste("Data/", filename, "-k3-test.Rds", sep=""))




# close log
sink()


