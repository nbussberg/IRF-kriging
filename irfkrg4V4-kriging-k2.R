### Code that tests IRF kriging on IRF2 data
##  Created by: Nicholas Bussberg
##  Date last updated: 2022-11-29
##  Filename: irfkrg4V4-kriging-k2.R
##  Purpose: test IRF kriging predictions and r estimates for IRF2 data

### Setup ###
rm(list=ls())
options(width = 80)

filename <-"irfkrg4V4-kriging-k2"

sink(paste(filename, ".log", sep = "")) # create log file

# R version
R.version.string


### 00 Load functions

source(file = "irfkrg-fcn01V3-sphere-harmonics.R")
source(file = "irfkrg-fcn02V3-legendre-polynomials.R")
source(file = "irfkrg-fcn03V3-distance.R")
source(file = "irfkrg-fcn04V3-icf.R")
source(file = "irfkrg-fcn05V3-MoM-estimator.R")
source(file = "irfkrg-fcn06V3-order-approx.R")
source(file = "irfkrg-fcn07V4-criteria.R")
source(file = "irfkrg-fcn08V3-r-estimation.R")
source(file = "irfkrg-fcn09V3-kriger.R")

library(MASS)


###  01 Load data

X.Lmat <- readRDS("Data/irfkrg3V3-kriging-split-k2-full.Rds")

train.data <- readRDS("Data/irfkrg3V3-kriging-split-k2-train.Rds")
test.data  <- readRDS("Data/irfkrg3V3-kriging-split-k2-test.Rds")

X0.train <- train.data$X.0.k2
Lmat.train <- train.data[, 2:3]
P.train <- length(X0.train)

# Load test data

X0.test <- test.data$X.0.k2
Lmat.test <- test.data[, 2:3]
P.test <- length(X0.test)


### 02 Set parameters

H <- 10         # Number of distances at which to estimate
Hvec <- as.vector(seq(0, 1, by = 1/H))[1:H]
eps <- 0.10     # use to bin distances in estimation

# Get distances
Dvec.train <- Dfun(Lmat.train, P.train)



### 03 Estimate kappa

GN.train <- order.approx(data = X0.train, Lmat.train, Dvec.train, 
                        Hvec, eps, P.train)

G.diff <- G.diff.fcn.noadjust(GN.train$G, H, Hvec)

crit.raw <- criterion.v01(G.diff)

# original units graph
pdf(file = paste("Graphs/", filename, "-crit-plot.pdf", sep = ""),
    width = 10, height = 7)

par(las=1, mar=c(5,6.2,2,1)+0.1)
plot(0:6, crit.raw, xlab = "k", ylab = "", cex = 1.2,
     cex.axis = 1.2, cex.lab = 1.2)
mtext("S(k)", 2, line=2, cex = 1.2, las=2, at=700)

dev.off()

# log units graph
pdf(file = paste("Graphs/", filename, "-crit-plot-log.pdf", sep = ""),
    width = 10, height = 7) 
par(las=1, mar=c(5,6.2,2,1)+0.1)
plot(0:6, log(crit.raw), xlab = "k", ylab = "", 
     cex = 1.2, cex.axis = 1.2, cex.lab = 1.2)
mtext("log(S(k))", 2, line=0.5, cex = 1.2, las=2, at=7)
dev.off()


# both version imply kappa-hat is 2



###  04 Estimate r using weighted least squares

# estimate r assuming k=1 (Ordinary kriging)
fit.k1 <- nls(GN.train$G[,2] ~ rfun.k(r, 1), data = GN.train, 
              weights = 1/GN.train$N[,2], start = list(r = .5)) 
b.icf.k1 <- coef(fit.k1)

print(paste("Estimate for r assuming k=1 is ", b.icf.k1, sep = ""))


# estimate r with kappa estimate (k=2)
fit.k2 <- nls(GN.train$G[,3] ~ rfun.k(r, 2), data = GN.train, 
              weights = 1/GN.train$N[,3], start = list(r = .5)) 
(b.icf.k2 <- coef(fit.k2))

print(paste("Estimate for r assuming k=2 is ", b.icf.k2, sep = ""))





###  05 Set parameters and estimates for kriging

r.true <- 0.75
r.est.k1 <- b.icf.k1
r.est.k2 <- b.icf.k2

sig <- 0

Dmat.train <- matrix(0,P.train,P.train)
Dmat.train[upper.tri(Dmat.train,diag = TRUE)] <- Dvec.train
Dmat.train <- Dmat.train + t(Dmat.train)



###  06 Set prelim matrices for IRFk method

PQinv.k1 <- PQinv.fcn(1, Dmat.train, r.est.k1, P.train, Lmat.train, sig)
PQinv.k2 <- PQinv.fcn(2, Dmat.train, r.est.k2, P.train, Lmat.train, sig)


###  07 Make predictions

# IRF kriging 

X.krig.k1 <- X.krig.k2 <- matrix(0, P.test, 1)

for (i in 1:P.test) {
    X.krig.k1[i] <- kriger.icfk(Lmat.test[i,], 1, Lmat.train, P.train, 
                                PQinv.k1, X0.train, sig, r.est.k1)
    X.krig.k2[i] <- kriger.icfk(Lmat.test[i,], 2, Lmat.train, P.train, 
                                PQinv.k2, X0.train, sig, r.est.k2)
}


# Calculate RMSE 

k1.rmse <- (sum((X0.test - X.krig.k1)^2) / P.test)^(1/2)

print(paste("RMSE for kriging using k=1: ", k1.rmse, sep = ""))


k2.rmse <- (sum((X0.test - X.krig.k2)^2) / P.test)^(1/2)

print(paste("RMSE for kriging using k=2: ", k2.rmse, sep = ""))





# close log
sink()


