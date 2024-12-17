### Function to create MoM estimator
##  Created by: Nicholas Bussberg and Jacob Shields
##  Date last updated: 2022-12-12
##  Filename: irf-fcn05V3-MoM-estimator.R
##  Purpose: create the MoM estimator function that will be used in estimating
    ## the degree of non-homogeneity, kappa, based on IRFk theory


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
