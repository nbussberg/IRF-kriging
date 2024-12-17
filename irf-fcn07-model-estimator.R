## Functions to estimate spatial models ##

# Script metadata --------------------------------------------------------------

# Created by: Nicholas Bussberg, Chunfeng Huang, Jacob Shields
# Date last updated: 2024-08-15
# Filename: irf-fcn07-model-estimator.R
# Purpose: functions to estimate spatial models


# Spatial model with r and sigma -----------------------------------------------

rfun.k <- function(r, k, sig){
  # r: spatial model parameter
  # k: order of IRF
  # sig: sigma model parameter 
    
  est <- (((1 - 2*r*cos(Hvec) + r^2)^(-3/2)) * (1 - r^2) / (4*pi))*sig^2
  if(k>=1) est <- est - sig^2*1/(4*pi) * P.0(Hvec)
  if(k>=2) est <- est - sig^2*r*3/(4*pi) * P.1(cos(Hvec))
  if(k>=3) est <- est - sig^2*5*(r^2) / (4*pi) * P.2(cos(Hvec))
    
  return(est)
}



