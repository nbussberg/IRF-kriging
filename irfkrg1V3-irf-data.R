### Code that simulates IRFk data 
##  Created by: Nicholas Bussberg
##  Date last updated: 2022-12-13
##  Filename: irfkrg1V3-irf-data.R
##  Purpose: simulates IRF2 and IRF3 data

## WARNING: This code may take 3-4 hours to run ##

### Setup ###
rm(list=ls())
options(width = 80)

filename <-"irfkrg1V3-irf-data"

sink(paste(filename, ".log", sep = "")) # create log file to store code & output

# R version
R.version.string



### 00 Load functions

source(file = "irfkrg-fcn01V3-sphere-harmonics.R")
source(file = "irfkrg-fcn02V3-legendre-polynomials.R")
source(file = "irfkrg-fcn03V3-distance.R")
source(file = "irfkrg-fcn04V3-icf.R")
source(file = "irfkrg-fcn05V3-MoM-estimator.R")
source(file = "irfkrg-fcn06V3-order-approx.R")

library(fields)   # needed for image.plot()
library(maptools) # needed for wrld_simpl map - requires package sp
library(akima)    # for interp()



### 01 Initial setup

user.seed <- 1115    # random seed to investigate
r <- 0.75
P <- 1500

# set distance bins
H <- 10         # Number Distances at Which to Estimate
eps <- .10      # Use to Bin Distances in Estimation
Hvec <- as.vector(seq(0,1,by = 1/H))[1:H]


## Create tau matrices to test
# matrix values are lat, long

tau2 <- matrix(c(pi/9,  pi/3,
                 pi/3,  5*pi/6, 
                 2*pi/3, 6*pi/5,
                 8*pi/9, 5*pi/3),
               nrow = 4, ncol = 2, byrow = T)

tau3 <- matrix(c(pi/12, pi/6,
                 pi/9,  pi/3, # from k=2 tau1
                 pi/6,  2*pi/3,
                 pi/3,  5*pi/6, # from k=2 tau1
                 pi/2,  pi,
                 2*pi/3, 6*pi/5, # from k=2 tau1
                 5*pi/6, 3*pi/2,
                 8*pi/9, 5*pi/3, # from k=2 tau1
                 11*pi/12, 9*pi/5),
               nrow = 9, ncol = 2, byrow = T)



## Create corresponding A matrices (see Shields 2017)

A2 <- solve(matrix(c(apply(tau2,1,Y00),
                     apply(tau2,1,Y1n1),
                     apply(tau2,1,Y10),
                     apply(tau2,1,Y11)),
                   4, 4, byrow = T))

A3 <- solve(matrix(c(apply(tau3,1,Y00),
                     apply(tau3,1,Y1n1),
                     apply(tau3,1,Y10),
                     apply(tau3,1,Y11),
                     apply(tau3,1,Y2n2),
                     apply(tau3,1,Y2n1),
                     apply(tau3,1,Y20),
                     apply(tau3,1,Y21),
                     apply(tau3,1,Y22)),
                   9, 9, byrow = T))


## Create location matrix

set.seed(user.seed)

Lmat.1 <- matrix(runif(P*2,0,1),ncol=2)
Lmat.1[,1] <- Lmat.1[,1]*pi #lat
Lmat.1[,2] <- Lmat.1[,2]*2*pi #lon
Lmat <- Lmat.1

saveRDS(Lmat, paste("Data/", filename, "-Lmat.Rds", sep = ""))


## Create distance vector
Dvec <- Dfun(Lmat, P)



### 02 run loop to generate ud and Gmat

G.mat.k <- list()
G.index <- 1

kappa.vec <- c(2, 3)

for(k in 1:length(kappa.vec)){
    
    kappa <- kappa.vec[k]
    
    # Compute Cmat
    Cmat <- matrix(0, P, P)
    for (i in 1:P) {
        for (j in i:P) {
            if(kappa==2) Cmat[i,j] <- Cmat[j,i] <- CF2(Lmat[i,], Lmat[j,], 
                                                       A2, r, tau2)
            if(kappa==3) Cmat[i,j] <- Cmat[j,i] <- CF3(Lmat[i,], Lmat[j,], 
                                                       A3, r, tau3)
        }
    }
    
    # Compute svd and ud
    Cmat.svd <- svd(Cmat)
    ud <- Cmat.svd$u %*% diag(sqrt(Cmat.svd$d), nrow = P, ncol = P)
    
    saveRDS(ud, paste("Data/", filename, "-ud-k", kappa, ".Rds", sep = ""))
    
    # Compute X.0
    X0 <- ud %*% matrix(rnorm(P,0,1),P,1)
    
    # Compute G
    G <- order.approx(X0, Lmat, Dvec, Hvec, eps, P)
    
    G.mat.k[[G.index]] <- G$G
    
    G.index <- G.index + 1
    
    
    # Make plot
    loc.data <- cbind((Lmat[,2]-pi)*(180/pi),(Lmat[,1] - pi/2)*(180/pi), X0)
    
    val.int <- akima::interp(loc.data[,1], loc.data[,2], X0)
    
    X.lim <- max(max(X0), abs(min(X0)))
    
    pdf(file = paste("Graphs/", filename, kappa, ".pdf", sep = ""), 
        width = 10, height = 7)
    par(las=1)
    par(mar=c(2.5,3,1,1) + 0.1)
    
    fields::image.plot(val.int, cex.axis = 1.2, xlim = c(-180, 180), 
               ylim = c(-90, 90), zlim = c(-X.lim, X.lim))
    
    data(wrld_simpl)
    plot(wrld_simpl,border = "black",add = T)
    
    if(kappa==2){
        points((tau2[,2] - pi)*(180/pi),(tau2[,1] - pi/2)*(180/pi), 
               pch = 24, cex = 1.7, col = "black", bg = "white")
    }
    if(kappa==3){
        points((tau3[,2] - pi)*(180/pi),(tau3[,1] - pi/2)*(180/pi), 
               pch = 24, cex = 1.7, col = "black", bg = "white")
    }
    dev.off()
    
}



## 03 Save G matrices

saveRDS(G.mat.k, paste("Data/", filename, "-Gmat.Rds", sep = ""))



# close log
sink()


