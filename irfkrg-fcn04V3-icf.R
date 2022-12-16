### Function to create the intrinsic covariance functions ###
##  Created by: Nicholas Bussberg and Jacob Shields
##  Date last updated: 2022-12-12
##  Filename: irf-fcn04V3-icf.R
##  Purpose: create the intrinsic covariance functions (ICFs) for 
    ## simulating the covariance of IRF data


ICFk <- function(x, r, k){
    # x: value (spherical distance)
    # r: spatial parameter
    # k: order of IRF
    
    icf <- ((1 - 2*r*cos(x) + r*r)^(-3/2))*(1-r^2)/(4*pi)
    if(k>=1) icf <- icf - 1/(4*pi)
    if(k>=2) icf <- icf - 3*r*cos(x)/(4*pi)
    if(k>=3) icf <- icf - 5*(r^2)*P.2(cos(x))/(4*pi)
    if(k>=4) icf <- icf - 7*(r^3)*P.3(cos(x))/(4*pi)
    if(k>=5) icf <- icf - 9*(r^4)*P.4(cos(x))/(4*pi)
    if(k>=6) icf <- icf - 11*(r^5)*P.5(cos(x))/(4*pi)
    if(k>=7) icf <- icf - 13*(r^6)*P.6(cos(x))/(4*pi)
    return(icf)
}



# CF0 - only dependent on distance b/w points (homogeneous)
# k=0
CF0 <- function(x, y, r) {
    ICFk(SphDist(x[1],x[2],y[1],y[2]), r, 0)
}



# CF1 - intrinsically homogeneous, k=1
p1 <- function(x, A) {
    A[1,1] * Y00(x)
}

CF1 <- function(x, y, A, r, tau) {
    rk1 <- rk2 <- rk3 <- 0
    v <- u <- 1
    
    rk1 <- rk1 + 
        ICFk(SphDist(x[1],x[2],tau[v,1],tau[v,2]), r, 1) * p1(y, A) + 
        ICFk(SphDist(y[1],y[2],tau[v,1],tau[v,2]), r, 1) * p1(x, A)
    rk3 <- rk3 + p1(x, A) * p1(y, A)
    rk2 <- rk2 + ICFk(SphDist(tau[v,1],tau[v,2],tau[u,1],tau[u,2]), r, 1) * 
        p1(x, A) * p1(y, A)
    ICFk(SphDist(x[1],x[2],y[1],y[2]), r, 1) - rk1 + rk2 + rk3
}



# CF2, k=2

p2 <- function(x, i, A){
    Y <- c(Y00(x), Y1n1(x), Y10(x), Y11(x))
    as.numeric(A[i,1:4] %*% Y)
}

CF2 <- function(x, y, A, r, tau) {
    rk1 <- rk2 <- rk3 <- 0
    for (v in 1:4) {
        rk1 <-  rk1 + 
            ICFk(SphDist(x[1],x[2],tau[v,1],tau[v,2]), r, 2) * p2(y, v, A) + 
            ICFk(SphDist(y[1],y[2],tau[v,1],tau[v,2]), r, 2) * p2(x, v, A)
        rk3 <-  rk3 + p2(x, v, A) * p2(y, v, A)
        for (u in 1:4) {
            rk2 <- rk2 + ICFk(SphDist(tau[v,1],tau[v,2],tau[u,1],tau[u,2]), r, 2) * 
                    p2(x, v, A) * p2(y, u, A)
        }
    }
    ICFk(SphDist(x[1],x[2],y[1],y[2]), r, 2) - rk1 + rk2 + rk3
}




# CF3, k=3

p3 <- function(x, i, A){
    Y <- c(Y00(x), Y1n1(x), Y10(x), Y11(x), Y2n2(x), Y2n1(x), 
           Y20(x), Y21(x), Y22(x))
    as.numeric(A[i,1:9] %*% Y)
}

CF3 <- function(x, y, A, r, tau) {
    rk1 <- rk2 <- rk3 <- 0
    for (v in 1:9) {
        rk1 <- rk1 + 
            ICFk(SphDist(x[1],x[2],tau[v,1],tau[v,2]), r, 3) * p3(y, v, A) + 
            ICFk(SphDist(y[1],y[2],tau[v,1],tau[v,2]), r, 3) * p3(x, v, A) 
        rk3 <- rk3 + p3(x, v, A) * p3(y, v, A)
        for (u in 1:9) {
            rk2 <- rk2 + ICFk(SphDist(tau[v,1],tau[v,2],tau[u,1],tau[u,2]), r, 3) * 
                p3(x, v, A) * p3(y, u, A)
        }
    }
    ICFk(SphDist(x[1],x[2],y[1],y[2]), r, 3) - rk1 + rk2 + rk3
}










