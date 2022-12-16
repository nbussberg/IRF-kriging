### Functions to estimate the spatial parameter r ###
##  Created by: Nicholas Bussberg and Jacob Shields
##  Date last updated: 2022-12-12
##  Filename: irf-fcn08V3-r-estimation.R
##  Purpose: create the function to estimate the spatial parameter r


rfun.k <- function(r, k){
    ## Estimates the spatial parameter r (estimate is r.est)
    # r: spatial parameter
    # k: order of IRF
    
    r.est <- (((1 - 2*r*cos(Hvec) + r^2)^(-3/2)) * 
                ((2*r*cos(Hvec) - 2*r^2) + (1 - 2*r*cos(Hvec) + r^2)) / 
                (4*pi))
    if(k>=1) r.est <- r.est - 1/(4*pi) * P.0(Hvec)
    if(k>=2) r.est <- r.est - r*3/(4*pi) * P.1(cos(Hvec))
    if(k>=3) r.est <- r.est - 5*(r^2) / (4*pi) * P.2(cos(Hvec))
    
    return(r.est)
}


    


