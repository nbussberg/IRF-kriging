## Code that simulates IRFk data ##


# Script metadata --------------------------------------------------------------

# Created by: Nicholas Bussberg and Chunfeng Huang
# Date last updated: 2024-12-16
# Filename: irfkrg1-data-r-sig-model.R
# Purpose: simulates IRF2 and IRF3 data based on the model:
  # sig^2 * (1-r^2)/(4 pi) * (1 - 2r cos(h) +r^2)^(-3/2), 0<=r<1

## WARNING: This code may take 3-4 hours to run ##


# Prelim setup -----------------------------------------------------------------

rm(list=ls())
options(width = 80)

filename <-"irfkrg1-data-r-sig-model"

# create log file to store code & output
sink(paste(filename, ".log", sep = "")) 

# R version
R.version.string


# 00 Load functions and packages -----------------------------------------------

source(file = "irf-fcn01-sphere-harmonics.R")
source(file = "irf-fcn02-legendre-polynomials.R")
source(file = "irf-fcn03-distance.R")
source(file = "irf-fcn04-covariances.R")
source(file = "irf-fcn05-MoM-estimator.R")

library(fields)   # needed for image.plot()
library(maps)     # needed for world map
library(akima)    # for interp()


# 01 Initial setup -------------------------------------------------------------

user.seed <- 1115  # random seed
r <- 0.75   # parameter for model
sig <- 1.0  # parameter for model
P <- 1500   # number of points to simulate

# set bins for distances
H <- 10
Hvec <- c(0, 0.02, 0.04 ,0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18)
eps <- 0.001


## Create tau and A matrices ---------------------------------------------------
# See Shields 2017 for details on tau and A forms
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

# Create corresponding A matrices (see Shields 2017)

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


## Create location matrix and distance vector ----------------------------------

set.seed(user.seed)

# location matrix
Lmat.1 <- matrix(runif(P*2,0,1),ncol=2)
Lmat.1[,1] <- Lmat.1[,1]*pi #lat
Lmat.1[,2] <- Lmat.1[,2]*2*pi #lon
Lmat <- Lmat.1

saveRDS(Lmat, paste("data/", filename, "-Lmat.Rds", sep = ""))

# distance vector
Dvec <- Dfun(Lmat, P)



# 02 Generate ud and Gmat ------------------------------------------------------

G.mat.k <- list()
G.index <- 1

kappa.vec <- c(2, 3)

for(k in 1:length(kappa.vec)){
   
  kappa <- kappa.vec[k]
  
  # Compute Cmat
  Cmat <- matrix(0, P, P)
  for (i in 1:P) {
    for (j in i:P) {
      if(kappa==2) Cmat[i,j] <- Cmat[j,i] <- CF2.r(Lmat[i,], Lmat[j,], 
                                                 A2, r, sig, tau2)
      if(kappa==3) Cmat[i,j] <- Cmat[j,i] <- CF3.r(Lmat[i,], Lmat[j,], 
                                                 A3, r, sig, tau3)
    }
  }
  
  # Compute svd and ud
  Cmat.svd <- svd(Cmat)
  ud <- Cmat.svd$u %*% diag(sqrt(Cmat.svd$d), nrow = P, ncol = P)
  
  saveRDS(ud, paste("data/", filename, "-ud-k", kappa, ".Rds", sep = ""))
  
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
  
  pdf(file = paste("graphs/", filename, kappa, ".pdf", sep = ""), 
      width = 10, height = 7)
  par(las=1)
  par(mar=c(2.5,3,1,1) + 0.1)
  
  fields::image.plot(val.int, cex.axis = 1.2, xlim = c(-180, 180), 
             ylim = c(-90, 90), zlim = c(-X.lim, X.lim))
  
  maps::map("world", add = T)
  
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



# 03 Save G matrices -----------------------------------------------------------

saveRDS(G.mat.k, paste("data/", filename, "-Gmat.Rds", sep = ""))



# close log
sink()


