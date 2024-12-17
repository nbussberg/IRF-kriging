### Code that tests IRF kriging on IRF3 data
##  Created by: Nicholas Bussberg
##  Date last updated: 2024-12-16
##  Filename: irfkrg5-kriging-k3.R
##  Purpose: test IRF kriging predictions and r estimates for IRF3 data

### Setup ###
rm(list=ls())
options(width = 80)

filename <-"irfkrg5-kriging-k3"

sink(paste(filename, ".log", sep = "")) # create log file

# R version
R.version.string


### 00 Load functions

source(file = "irf-fcn01-sphere-harmonics.R")
source(file = "irf-fcn02-legendre-polynomials.R")
source(file = "irf-fcn03-distance.R")
source(file = "irf-fcn04-covariances.R")
source(file = "irf-fcn05-MoM-estimator.R")
source(file = "irf-fcn06-kappa-approx.R")
source(file = "irf-fcn07-model-estimator.R")
source(file = "irf-fcn08-kriger.R")

library(MASS)


###  01 Load data

X.Lmat <- readRDS("data/irfkrg3-kriging-split-k3-sig1-full.Rds")

train.data <- readRDS("data/irfkrg3-kriging-split-k3-sig1-train.Rds")
test.data  <- readRDS("data/irfkrg3-kriging-split-k3-sig1-test.Rds")

X0.train <- train.data$X.0.k3
Lmat.train <- train.data[, 2:3]
P.train <- length(X0.train)

# Load test data

X0.test <- test.data$X.0.k3
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
pdf(file = paste("graphs/", filename, "-crit-plot.pdf", sep = ""),
    width = 10, height = 7)

par(las=1, mar=c(5,6.2,2,1)+0.1)
plot(0:6, crit.raw, xlab = "k", ylab = "", cex = 1.2,
     cex.axis = 1.2, cex.lab = 1.2)
mtext("S(k)", 2, line=2, cex = 1.2, las=2, at=700)

dev.off()

# log units graph
pdf(file = paste("graphs/", filename, "-crit-plot-log.pdf", sep = ""),
    width = 10, height = 7) 
par(las=1, mar=c(5,6.2,2,1)+0.1)
plot(0:6, log(crit.raw), xlab = "k", ylab = "", 
     cex = 1.2, cex.axis = 1.2, cex.lab = 1.2)
mtext("log(S(k))", 2, line=0.5, cex = 1.2, las=2, at=7)
dev.off()


# both version imply kappa-hat is 3



###  04 Estimate r and sig using iterative weighted least squares

# estimate r assuming k=1 (Ordinary kriging)
fit.k1 <- nls(GN.train$G[,2] ~ rfun.k(r, 1, sig = 1), data = GN.train, 
              weights = 1/GN.train$N[,2], start = list(r = .85)) 
b.icf.k1 <- coef(fit.k1)

print(paste("Estimate for r assuming k=1 is ", b.icf.k1, sep = ""))


# estimate r with kappa estimate (k=3)
fit.k3 <- nls(GN.train$G[,4] ~ rfun.k(r, 3, sig = 1), data = GN.train, 
              weights = 1/GN.train$N[,4], start = list(r = .5)) 
(b.icf.k3 <- coef(fit.k3))

print(paste("Estimate for r assuming k=3 is ", b.icf.k3, sep = ""))



###  05 Set parameters and estimates for kriging

r.true <- 0.75

Dmat.train <- matrix(0,P.train,P.train)
Dmat.train[upper.tri(Dmat.train,diag = TRUE)] <- Dvec.train
Dmat.train <- Dmat.train + t(Dmat.train)



###  06 Set prelim matrices for IRFk method

PQinv.k1 <- PQinv.fcn(1, Dmat.train, b.icf.k1, P.train, Lmat.train, 
                      1, 0)
PQinv.k3 <- PQinv.fcn(3, Dmat.train, b.icf.k3, P.train, Lmat.train, 
                      1, 0)


###  07 Make predictions

# IRF kriging 

X.krig.k1 <- X.krig.k3 <- matrix(0, P.test, 1)

for (i in 1:P.test) {
    X.krig.k1[i] <- kriger.icfk(Lmat.test[i,], 1, Lmat.train, P.train, 
                                PQinv.k1, X0.train, 0, 
                                b.icf.k1, 1)
    X.krig.k3[i] <- kriger.icfk(Lmat.test[i,], 3, Lmat.train, P.train, 
                                PQinv.k3, X0.train, 0, 
                                b.icf.k3, 1)
}


# Calculate RMSE 

k1.rmse <- (sum((X0.test - X.krig.k1)^2) / P.test)^(1/2)

print(paste("RMSE for kriging using k=1: ", k1.rmse, sep = ""))


k3.rmse <- (sum((X0.test - X.krig.k3)^2) / P.test)^(1/2)

print(paste("RMSE for kriging using k=3: ", k3.rmse, sep = ""))





# close log
sink()


