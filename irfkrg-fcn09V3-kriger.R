### Functions to implement kriging using IRFs ###
##  Created by: Nicholas Bussberg and Chunfeng Huang
##  Date last updated: 2022-12-12
##  Filename: irf-fcn09V3-kriger.R
##  Purpose: create the functions needed to implement IRF-based kriging


PQinv.fcn <- function(k, Dmat, r, P, Lmat, sig){
    ##  Sets prelim matrices for IRFk method
    # k: order of IRF (must be greater than 0)
    # Dmat: matrix created with Dvec
    # r: spatial correlation parameter estimate
    # P: number of locations
    # Lmat: location matrix
    # sig: sigma value
    # must have spherical harmonic functions & ICF functions loaded
    
    if(k==0) print("k must be greater than 0")
    
    Psi.irf <- ICFk(Dmat, r, k) 
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
    PQ.irf <- rbind(cbind(Psi.irf + diag(sig^2, P), Q.irf), 
                    cbind(t(Q.irf), zero.irf))
    PQinv.irf <- ginv(PQ.irf)
    
    return(PQinv.irf)
}



kriger.icfk <- function(x, k, Lmat, P, PQinv, X.0, sig, r){
    ## Kriging function
    # x: Lmat[1, ] 
    # k: order of IRF (must be greater than 0)
    # Lmat: location matrix
    # P: number of points
    # PQinv: prelim matrix from previous function
    # X.0: data
    # sig: sigma value
    # r: spatial correlation parameter estimate
    # must have spherical harmonic functions & ICF functions loaded
    
    xtemp <- as.numeric(x)
    
    psi <- matrix(0,P,1)
    for (i in 1:P) {
        psi[i,] <- ICFk(SphDist(Lmat[i,1],Lmat[i,2],xtemp[1],xtemp[2]), r, k)
        }
    
    if (k==1) q <- matrix(c(Y00(xtemp)),k^2,1)
    if (k==2){
        q <- matrix(c(Y00(xtemp), Y1n1(xtemp), Y10(xtemp), Y11(xtemp)),k^2,1)
    }
    if (k==3){
        q <- matrix(c(Y00(xtemp), Y1n1(xtemp), Y10(xtemp), 
                      Y11(xtemp), Y2n2(xtemp), Y2n1(xtemp),
                      Y20(xtemp), Y21(xtemp),  Y22(xtemp)),k^2,1)
    }
    eta <- PQinv %*% rbind(psi,q)
    z <- t(eta[1:P,]) %*% X.0
    
    return(c(z))
}


