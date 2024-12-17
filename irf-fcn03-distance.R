## Functions to calculate spherical distances ##

# Script metadata --------------------------------------------------------------

# Created by: Nicholas Bussberg, Jacob Shields, Chunfeng Huang
# Date last updated: 2024-08-15
# Filename: irf-fcn03-distance.R
# Purpose: functions to calculate Great Circle distance and distances 
  # between points in a location matrix


# Great Circle Distance between two locations ----------------------------------

SphDist <- function(lat1, long1, lat2, long2) {
  # Finds Great Circle distance between two locations

  d <- sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(abs(long2 - long1))
  
  # Vincenty Formula
  n <- sqrt((cos(lat2)*sin(abs(long2 - long1)))^2 + 
            (cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(abs(long2 - long1)))^2)
  
  atan2(n,d)
}


# Distance function to determine distances between all pairs of points ---------

Dfun <- function(Lmat, P) {
  # returns vector of all distance pairs
  # Lmat: location matrix
  # P: number of locations
  
  Dmat <- matrix(0,P,P)
  for (i in 1:P) {
    for (j in i:P) {
      Dmat[i,j] <- SphDist(Lmat[i,1],Lmat[i,2],Lmat[j,1],Lmat[j,2])
    }
  } 

  Dvec <- Dmat[upper.tri(Dmat,diag = TRUE)] 

  return(Dvec)
}

