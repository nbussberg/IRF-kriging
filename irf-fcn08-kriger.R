## Function for kriging using IRFs ##

# Script metadata --------------------------------------------------------------

# Created by: Nicholas Bussberg, Chunfeng Huang
# Date last updated: 2024-08-15
# Filename: irf-fcn08-kriger.R
# Purpose: functions to implement IRF-based kriging


# Preliminary function for matrices needed for IRFk method ---------------------

PQinv.fcn <- function(k, Dmat, r, P, Lmat, sig, sigma){
  # must have spherical harmonic functions & ICF functions loaded  
  # k: order of IRF (must be greater than 0)
  # Dmat: matrix created with Dvec
  # r: spatial  parameter estimate
  # sig: spatial parameter estimate
  # P: number of locations
  # Lmat: location matrix
  # sigma: sigma value
  
  if(k==0) print("k must be greater than 0")
  
  Psi.irf <- ICFk.r(Dmat, r, sig, k) 
  Q.irf <- matrix(0, P, k^2) 
  for (i in 1:P) {
    Lmat.temp <- as.numeric(Lmat[i, ])
    if (k==1) Q.irf[i,] <- c(Y00(Lmat.temp))
    if (k==2){
        Q.irf[i,] <- c(Y00(Lmat.temp),Y1n1(Lmat.temp),
                       Y10(Lmat.temp),Y11(Lmat.temp))
    }
    if (k==3){
        Q.irf[i,] <- c(Y00(Lmat.temp),Y1n1(Lmat.temp),Y10(Lmat.temp),
                       Y11(Lmat.temp),Y2n2(Lmat.temp),Y2n1(Lmat.temp),
                       Y20(Lmat.temp),Y21(Lmat.temp),Y22(Lmat.temp))
    }
  }
  zero.irf <- matrix(0,k^2,k^2)
  PQ.irf <- rbind(cbind(Psi.irf + diag(sigma^2, P), Q.irf), 
                  cbind(t(Q.irf), zero.irf))
  PQinv.irf <- ginv(PQ.irf)
  
  return(PQinv.irf)
}


# Kriging function -------------------------------------------------------------

kriger.icfk <- function(x, k, Lmat, P, PQinv, X.0, sigma, r, sig){
  # must have spherical harmonic functions & ICF functions loaded
  # x: Lmat[1, ] 
  # k: order of IRF (must be greater than 0)
  # Lmat: location matrix
  # P: number of points
  # PQinv: prelim matrix from previous function
  # X.0: data
  # sigma: sigma value
  # r: spatial parameter estimate
  # sig: spatial parameter estimate
  
  xtemp <- as.numeric(x)
  
  psi <- matrix(0,P,1)
  for (i in 1:P) {
    psi[i,] <- ICFk.r(SphDist(Lmat[i,1], Lmat[i,2], xtemp[1], xtemp[2]), 
                      r, sig, k)
    }
  
  if (k==1) q <- matrix(c(Y00(xtemp)), k^2, 1)
  if (k==2){
    q <- matrix(c(Y00(xtemp), Y1n1(xtemp), Y10(xtemp), Y11(xtemp)),k^2,1)
  }
  if (k==3){
    q <- matrix(c(Y00(xtemp), Y1n1(xtemp), Y10(xtemp), Y11(xtemp), 
                  Y2n2(xtemp), Y2n1(xtemp), Y20(xtemp), Y21(xtemp),  Y22(xtemp)), 
                k^2, 1)
  }
  eta <- PQinv %*% rbind(psi, q)
  z <- t(eta[1:P, ]) %*% X.0
  
  return(c(z))
}


