## Functions to create intrinsic covariance functions ##

# Script metadata --------------------------------------------------------------

# Created by: Nicholas Bussberg, Jacob Shields, Chunfeng Huang
# Date last updated: 2024-08-15
# Filename: irf-fcn04-covariances.R
# Purpose: create the intrinsic covariance functions (ICFs) for 
  # simulating the covariance of IRF data given parametric models


# Model 1 ----------------------------------------------------------------------
## Model: sig^2 * (1-r^2)/(4 pi) * (1 - 2r cos(h) +r^2)^(-3/2), 0<=r<1

ICFk.r <- function(x, r, sig, k){
  # x: spherical distance
  # r: spatial parameter
  # sig: scaling parameter
  # k: order of IRF
  
  icf <- ((1 - 2*r*cos(x) + r^2)^(-3/2))*(1-r^2)/(4*pi)*sig^2
  if(k>=1) icf <- icf - sig^2*1/(4*pi)
  if(k>=2) icf <- icf - sig^2*3*r*cos(x)/(4*pi)
  if(k>=3) icf <- icf - sig^2*5*(r^2)*P.2(cos(x))/(4*pi)
  if(k>=4) icf <- icf - sig^2*7*(r^3)*P.3(cos(x))/(4*pi)
  if(k>=5) icf <- icf - sig^2*9*(r^4)*P.4(cos(x))/(4*pi)
  if(k>=6) icf <- icf - sig^2*11*(r^5)*P.5(cos(x))/(4*pi)
  if(k>=7) icf <- icf - sig^2*13*(r^6)*P.6(cos(x))/(4*pi)
  return(icf)
}


# CF0.r - only dependent on distance b/w points (homogeneous), k=0
CF0.r <- function(x, y, r, sig) {
  ICFk.r(SphDist(x[1], x[2], y[1], y[2]), r, sig, 0)
}


# CF1.r - intrinsically homogeneous, k=1
p1 <- function(x, A) {
  A[1,1] * Y00(x)
}

CF1.r <- function(x, y, A, r, sig, tau) {
  rk1 <- rk2 <- rk3 <- 0
  v <- u <- 1
  
  rk1 <- ICFk.r(SphDist(x[1],x[2],tau[v,1],tau[v,2]), r, sig, 1) * p1(y, A) + 
         ICFk.r(SphDist(y[1],y[2],tau[v,1],tau[v,2]), r, sig, 1) * p1(x, A)
  rk3 <- sig^2 * p1(x, A) * p1(y, A)
  rk2 <- ICFk.r(SphDist(tau[v,1],tau[v,2],tau[u,1],tau[u,2]), r, sig, 1) * 
         p1(x, A) * p1(y, A)
  
  ICFk.r(SphDist(x[1],x[2],y[1],y[2]), r, sig, 1) - rk1 + rk2 + rk3
}


# CF2.r, k=2

p2 <- function(x, i, A){
  Y <- c(Y00(x), Y1n1(x), Y10(x), Y11(x))
  as.numeric(A[i,1:4] %*% Y)
}

CF2.r <- function(x, y, A, r, sig, tau) {
  rk1 <- rk2 <- rk3 <- 0
  
  for(v in 1:4){
    rk1 <- rk1 + 
        ICFk.r(SphDist(x[1],x[2],tau[v,1],tau[v,2]), r, sig, 2) * p2(y, v, A) + 
        ICFk.r(SphDist(y[1],y[2],tau[v,1],tau[v,2]), r, sig, 2) * p2(x, v, A)
    rk3 <- rk3 + sig^2 * p2(x, v, A) * p2(y, v, A)
    for(u in 1:4){
      rk2 <- rk2 + 
             ICFk.r(SphDist(tau[v,1],tau[v,2],tau[u,1],tau[u,2]), r, sig, 2) * 
             p2(x, v, A) * p2(y, u, A)
    }
  }
  
  ICFk.r(SphDist(x[1],x[2],y[1],y[2]), r, sig, 2) - rk1 + rk2 + rk3
}


# CF3.r, k=3

p3 <- function(x, i, A){
  Y <- c(Y00(x), Y1n1(x), Y10(x), Y11(x), Y2n2(x), Y2n1(x), 
         Y20(x), Y21(x), Y22(x))
  as.numeric(A[i,1:9] %*% Y)
}

CF3.r <- function(x, y, A, r, sig, tau) {
  rk1 <- rk2 <- rk3 <- 0
  
  for(v in 1:9){
    rk1 <- rk1 + 
        ICFk.r(SphDist(x[1],x[2],tau[v,1],tau[v,2]), r, sig, 3) * p3(y, v, A) + 
        ICFk.r(SphDist(y[1],y[2],tau[v,1],tau[v,2]), r, sig, 3) * p3(x, v, A) 
    rk3 <- rk3 + sig^2 * p3(x, v, A) * p3(y, v, A)
    for (u in 1:9) {
      rk2 <- rk2 + 
             ICFk.r(SphDist(tau[v,1],tau[v,2],tau[u,1],tau[u,2]), r, sig, 3) * 
             p3(x, v, A) * p3(y, u, A)
    }
  }
  
  ICFk.r(SphDist(x[1],x[2],y[1],y[2]), r, sig, 3) - rk1 + rk2 + rk3
}



# Model 2 ----------------------------------------------------------------------
## Model: e^(x / a)

ICFk.exp <- function(x, a, k){
  # x: spherical distance
  # a: scaling parameter
  # k: order of IRF
  
  icf <- exp(x/a)
  if(k>=1) icf <- icf - 2*pi*a*(exp(-pi/a) + 1) / (a^2 + 1)
  if(k>=2) icf <- icf - 2*pi*a^2*(1 - exp(-pi/a)) / (4*a^2 + 1)
  if(k>=3) icf <- icf - 2*pi*a^2*(exp(-pi/a) + 1) / (9*a^4 + 10*a^2 + 1)
  if(k>=4) icf <- icf - 2*pi*a^2*(a^2 + 1)*(1 - exp(-pi/a)) / (64*a^4 + 20*a^2 + 1)
  if(k>=5) icf <- icf - 2*pi*a^2*(exp(-pi/a) + 1)*(4*a^2 + 1) / 
      (225*a^6 + 259*a^4 + 35*a^2 + 1)
  if(k>=6) icf <- icf - 2*pi*a^2*(9*a^4 + 10*a^2 + 1)*(1 - exp(-pi/a)) / 
      (2304*a^6 + 784*a^4 + 56*a^2 + 1)
  if(k>=7) icf <- icf - 2*pi*a^2*(exp(-pi/a) + 1)*(64*a^4 + 20*a^2 + 1) / 
      (11025*a^8 + 12916*a^6 + 1974*a^4 + 84*a^2 + 1)
  return(icf)
}


# CF0.exp - only dependent on distance b/w points (homogeneous), k=0
CF0.exp <- function(x, y, a) {
  ICFk.exp(SphDist(x[1],x[2],y[1],y[2]), a, 0)
}


# CF1.exp - intrinsically homogeneous, k=1
p1 <- function(x, A) {
  A[1,1] * Y00(x)
}

CF1.exp <- function(x, y, A, a, tau) {
  rk1 <- rk2 <- rk3 <- 0
  v <- u <- 1
  
  rk1 <- ICFk.exp(SphDist(x[1],x[2],tau[v,1],tau[v,2]), a, 1) * p1(y, A) + 
    ICFk.exp(SphDist(y[1],y[2],tau[v,1],tau[v,2]), a, 1) * p1(x, A)
  rk3 <- sig^2 * p1(x, A) * p1(y, A)
  rk2 <- ICFk.exp(SphDist(tau[v,1],tau[v,2],tau[u,1],tau[u,2]), a, 1) * 
    p1(x, A) * p1(y, A)
  
  ICFk.r(SphDist(x[1],x[2],y[1],y[2]), r, sig, 1) - rk1 + rk2 + rk3
}


# CF2.exp, k=2

p2 <- function(x, i, A){
  Y <- c(Y00(x), Y1n1(x), Y10(x), Y11(x))
  as.numeric(A[i,1:4] %*% Y)
}

CF2.exp <- function(x, y, A, a, tau) {
  rk1 <- rk2 <- rk3 <- 0
  
  for(v in 1:4){
    rk1 <- rk1 + 
      ICFk.exp(SphDist(x[1],x[2],tau[v,1],tau[v,2]), a, 2) * p2(y, v, A) + 
      ICFk.exp(SphDist(y[1],y[2],tau[v,1],tau[v,2]), a, 2) * p2(x, v, A)
    rk3 <- rk3 + sig^2 * p2(x, v, A) * p2(y, v, A)
    for(u in 1:4){
      rk2 <- rk2 + 
        ICFk.exp(SphDist(tau[v,1],tau[v,2],tau[u,1],tau[u,2]), a, 2) * 
        p2(x, v, A) * p2(y, u, A)
    }
  }
  
  ICFk.exp(SphDist(x[1],x[2],y[1],y[2]), a, 2) - rk1 + rk2 + rk3
}


# CF3.exp, k=3

p3 <- function(x, i, A){
  Y <- c(Y00(x), Y1n1(x), Y10(x), Y11(x), Y2n2(x), Y2n1(x), 
         Y20(x), Y21(x), Y22(x))
  as.numeric(A[i,1:9] %*% Y)
}

CF3.exp <- function(x, y, A, a, tau) {
  rk1 <- rk2 <- rk3 <- 0
  
  for(v in 1:9){
    rk1 <- rk1 + 
      ICFk.exp(SphDist(x[1],x[2],tau[v,1],tau[v,2]), a, 3) * p3(y, v, A) + 
      ICFk.exp(SphDist(y[1],y[2],tau[v,1],tau[v,2]), a, 3) * p3(x, v, A) 
    rk3 <- rk3 + sig^2 * p3(x, v, A) * p3(y, v, A)
    for (u in 1:9) {
      rk2 <- rk2 + 
        ICFk.exp(SphDist(tau[v,1],tau[v,2],tau[u,1],tau[u,2]), a, 3) * 
        p3(x, v, A) * p3(y, u, A)
    }
  }
  
  ICFk.exp(SphDist(x[1],x[2],y[1],y[2]), a, 3) - rk1 + rk2 + rk3
}


