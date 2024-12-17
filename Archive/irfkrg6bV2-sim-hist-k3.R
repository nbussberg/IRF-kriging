### Code that simulates 1000 IRF3 datasets and plots kappa estimation
##  Created by: Nicholas Bussberg
##  Date last updated: 2022-12-14
##  Filename: irfkrg6bV2-sim-hist-k3.R
##  Purpose: simulates 1000 IRF3 datasets, estimates kappa for each, plot
## the kappa estimates on a histogram

### Setup ###
rm(list=ls())
options(width = 80)

filename <-"irfkrg6bV2-sim-hist-k3"

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


###  01 Load data

Lmat <- readRDS("Data/irfkrg1V3-irf-data-Lmat.Rds")

ud <- readRDS("Data/irfkrg1V3-irf-data-ud-k3.Rds")


# Parameters of ICF
H <- 10                   # Number Distances at Which to Estimate
eps <- .10                # Use to Bin Distances in Estimation
Hvec <- as.vector(seq(0,1,by = 1/H))[1:H]
P <- 1500


Dvec <- Dfun(Lmat, P)



### 02 Generate simulations, estimate kappa for each

nsim <- 1000

crit.1.nsim <- matrix(0,nsim,7)
crit.2.nsim <- matrix(0,nsim,7)

set.seed(384)

for(isim in 1:nsim){
    if(isim%%10 ==0) print(isim)
    
    data.X0 <- ud %*% matrix(rnorm(P,0,1),P,1)

    Gmat.GN <- order.approx(data.X0, Lmat, Dvec, Hvec, eps, P)
    
    G <- Gmat.GN$G
    
    G.diff.new <- G.diff.fcn.noadjust(G, H, Hvec)
    crit.1.nsim[isim,] <- criterion.v01(G.diff.new)
    
    G.diff.new2 <- G.diff.fcn(G, H, Hvec)
    crit.2.nsim[isim,] <- criterion.v01(G.diff.new2)
}


### 03 Plot kappa estimates on histogram

# Version 1: set alpha (penalty term) to 0.5
crit.3.nsim <- crit.2.nsim

crit.4.nsim <- rep(0,nsim)

for(isim in 1:nsim){
  value.max <- which.max(crit.2.nsim[isim,])
  crit.3.nsim[isim,] <- crit.2.nsim[isim,]+0.5*(0:6)
  crit.4.nsim[isim] <- which.min(crit.3.nsim[isim,value.max:7])+value.max-1	
}

table(crit.4.nsim)

pdf(file = paste("Graphs/", filename, "-alpha_half.pdf", sep = ""))
par(las=1, mar=c(5,4,0,0))
hist(crit.4.nsim, breaks = seq(0, 7),
     xlab='Order', main="",
     xaxt="n", ylim = c(0, 1000))
axis(side=1, c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5), 
     labels = c(0, 1, 2, 3, 4, 5, 6),
     line=0)
dev.off()


# Version 2: set alpha (penalty term) to 0.2
crit.5.nsim <- crit.2.nsim

crit.6.nsim <- rep(0,nsim)

for(isim in 1:nsim){
  value.max <- which.max(crit.2.nsim[isim,])
  crit.5.nsim[isim,] <- crit.2.nsim[isim,]+0.2*(0:6)
  crit.6.nsim[isim] <- which.min(crit.5.nsim[isim,value.max:7])+value.max-1	
}

table(crit.6.nsim)

pdf(file = paste("Graphs/", filename, "-alpha_fifth.pdf", sep = ""))
par(las=1, mar=c(5,4,0,0))
hist(crit.6.nsim, breaks = seq(0, 7),
     xlab='Order', main="",
     xaxt="n", ylim = c(0, 1000))
axis(side=1, c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5), 
     labels = c(0, 1, 2, 3, 4, 5, 6),
     line=0)
dev.off()



# Version 3: set alpha (penalty term) to 1.5
crit.7.nsim <- crit.2.nsim

crit.8.nsim <- rep(0,nsim)

for(isim in 1:nsim){
  value.max <- which.max(crit.2.nsim[isim,])
  crit.7.nsim[isim,] <- crit.2.nsim[isim,]+1.5*(0:6)
  crit.8.nsim[isim] <- which.min(crit.7.nsim[isim,value.max:7])+value.max-1	
}

table(crit.8.nsim)

pdf(file = paste("Graphs/", filename, "-alpha_one+half.pdf", sep = ""))
par(las=1, mar=c(5,4,0,0))
hist(crit.8.nsim, breaks = seq(0, 7),
     xlab='Order', main="",
     xaxt="n", ylim = c(0, 1000))
axis(side=1, c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5), 
     labels = c(0, 1, 2, 3, 4, 5, 6),
     line=0)
dev.off()


# close log
sink()




