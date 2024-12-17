## Functions for MoM estimators and MoM approximations ##

# Script metadata --------------------------------------------------------------

# Created by: Nicholas Bussberg, Chunfeng Huang, Jacob Shields
# Date last updated: 2024-08-15
# Filename: irf-fcn05-MoM-estimator.R
# Purpose: creates MoM estimator function for estimating the degree of 
  # non-homogeneity based on IRFk theory. 
  # Also creates function to approximate the MoM values for different degrees of 
  # non-homogeneity.


# Method of Moment Estimator ---------------------------------------------------

G_Hat <- function(Z, Dvec, Hvec, eps, P) {
  # Create Second Moment Matrix and Vector
  # Z: residuals obtained from regressing the full data process on the 
    # spherical harmonic functions
  # Dvec: vectors of distances between locations
  # Hvec: bins for distances
  # eps: used to bin distances
  # P: number of locations
  
  ZZmat <- matrix(0,P,P)
  ZZmat <- Z%*%t(Z)
  ZZvec <- ZZmat[upper.tri(ZZmat,diag = TRUE)]
  
  # Calculate Estimates using MOM Estimator
  G <- N <- numeric(H)
  s <- which(Dvec == 0)
  N[1] <- length(s)
  G[1] <- mean(ZZvec[s])
  for (i in 2:H) {
    s <- which(((Hvec[i] - eps) <  Dvec) & (Dvec <  (Hvec[i] + eps)))
    N[i] <- length(s)
    G[i] <- mean(ZZvec[s])
  }
  
  return(list(G = G, N = N))
  
}


# MoM Approximation Function ---------------------------------------------------

order.approx <- function(data, Lmat, Dvec, Hvec, eps, P){
  # data must be one time period
  # must have the Legendre polynomial functions loaded
  # Lmat: location matrix
  # Dvec: vectors of distances between points
  # Hvec: bins for distances
  # eps: used to bin distances
  # P: number of locations

  outcome <- "X.0"
  
  x1.vars <- "1"
  x2.vars <- paste(x1.vars,
                  "apply(Lmat,1,Y1n1) + apply(Lmat,1,Y10) + apply(Lmat,1,Y11)", 
                  sep = "+")
  x3.vars <- paste(x2.vars, 
                  "apply(Lmat,1,Y2n2) + apply(Lmat,1,Y2n1) + 
                  apply(Lmat,1,Y20) + apply(Lmat,1,Y21) + apply(Lmat,1,Y22)",
                  sep = "+")
  x4.vars <- paste(x3.vars, 
                  "apply(Lmat,1,Y3n3) + apply(Lmat,1,Y3n2) + apply(Lmat,1,Y3n1) + 
                  apply(Lmat,1,Y30) + apply(Lmat,1,Y31) + 
                  apply(Lmat,1,Y32) + apply(Lmat,1,Y33)", 
                  sep = "+")
  x5.vars <- paste(x4.vars, 
                  "apply(Lmat,1,Y4n4) + apply(Lmat,1,Y4n3) + 
                  apply(Lmat,1,Y4n2) + apply(Lmat,1,Y4n1) + 
                  apply(Lmat,1,Y40) + apply(Lmat,1,Y41) + apply(Lmat,1,Y42) + 
                  apply(Lmat,1,Y43) + apply(Lmat,1,Y44)", 
                  sep = "+")
  x6.vars <- paste(x5.vars, 
                  "apply(Lmat,1,Y5n5) + apply(Lmat,1,Y5n4) + 
                  apply(Lmat,1,Y5n3) + apply(Lmat,1,Y5n2) + apply(Lmat,1,Y5n1) + 
                  apply(Lmat,1,Y50) + apply(Lmat,1,Y51) + apply(Lmat,1,Y52) + 
                  apply(Lmat,1,Y53) + apply(Lmat,1,Y54) + apply(Lmat,1,Y55)", 
                  sep = "+")
  x7.vars <- paste(x6.vars, 
                  "apply(Lmat,1,Y6n6) + apply(Lmat,1,Y6n5) + apply(Lmat,1,Y6n4) + 
                  apply(Lmat,1,Y6n3) + apply(Lmat,1,Y6n2) + apply(Lmat,1,Y6n1) + 
                  apply(Lmat,1,Y60) + apply(Lmat,1,Y61) + apply(Lmat,1,Y62) + 
                  apply(Lmat,1,Y63) + apply(Lmat,1,Y64) + apply(Lmat,1,Y65) + 
                  apply(Lmat,1,Y66)", 
                  sep = "+")
  
  x1.form <- as.formula(paste(outcome, x1.vars, sep = "~"))
  x2.form <- as.formula(paste(outcome, x2.vars, sep = "~"))
  x3.form <- as.formula(paste(outcome, x3.vars, sep = "~"))
  x4.form <- as.formula(paste(outcome, x4.vars, sep = "~"))
  x5.form <- as.formula(paste(outcome, x5.vars, sep = "~"))
  x6.form <- as.formula(paste(outcome, x6.vars, sep = "~"))
  x7.form <- as.formula(paste(outcome, x7.vars, sep = "~"))
  
  X.0 <- data
  X.1 <- lm(x1.form)$res # order 1
  X.2 <- lm(x2.form)$res # order 2
  X.3 <- lm(x3.form)$res # order 3
  X.4 <- lm(x4.form)$res # order 4
  X.5 <- lm(x5.form)$res # order 5
  X.6 <- lm(x6.form)$res # order 6
  X.7 <- lm(x7.form)$res # order 7
  
  G.0 <- G_Hat(X.0, Dvec, Hvec, eps, P)
  G.1 <- G_Hat(X.1, Dvec, Hvec, eps, P)
  G.2 <- G_Hat(X.2, Dvec, Hvec, eps, P)
  G.3 <- G_Hat(X.3, Dvec, Hvec, eps, P)
  G.4 <- G_Hat(X.4, Dvec, Hvec, eps, P)
  G.5 <- G_Hat(X.5, Dvec, Hvec, eps, P)
  G.6 <- G_Hat(X.6, Dvec, Hvec, eps, P)
  G.7 <- G_Hat(X.7, Dvec, Hvec, eps, P)
  
  G <- t(rbind(G.0$G,G.1$G,G.2$G,G.3$G,G.4$G,G.5$G,G.6$G,G.7$G))
  N <- t(rbind(G.0$N,G.1$N,G.2$N,G.3$N,G.4$N,G.5$N,G.6$N,G.7$N))
  
  return(list(G = G, N = N))
}



