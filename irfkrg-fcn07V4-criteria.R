### Functions to estimate the degree of non-homogeneity, kappa ###
##  Created by: Nicholas Bussberg and Jacob Shields
##  Date last updated: 2022-12-12
##  Filename: irf-fcn07V4-criteria.R
##  Purpose: create the functions to estimate the degree of non-homogeneity,
##      kappa, based on IRFk theory



# Functions to calculate criterion for selecting kappa

G.diff.fcn <- function(G.mat, H, Hvec){
    ## MoM estimation that standardizes results by dividing by smallest value
    # G.mat: G matrix from each MoM estimation
    # H: for creating dimensions of G.diff.mat
    # Hvec: for evaluation in G.diff.mat
    
    G.diff.mat <- matrix(0, 7, H)
    for (i in 1:7) { # kappa
        for (j in 2:H) { # Hvec
            cos.Hvec.j <- cos(Hvec[j])
            G.diff.mat[i, j-1] <- (G.mat[j,i] - G.mat[j, i+1]) - 
                                        eval(parse(text = paste("P.",i - 1, "(", 
                                                                cos.Hvec.j, ")", 
                                                                sep = ""))) * 
                                  (G.mat[1,i] - G.mat[1,i + 1]) 
            
            G.diff.mat[i, j-1] <- G.diff.mat[i, j-1] / G.mat[1,8] 
                # divide by smallest value
        }
    }
    return(G.diff.mat)
}



G.diff.fcn.noadjust <- function(G.mat, H, Hvec){
    ## MoM estimation that does not standardize results (divide by smallest value)
    # G.mat: G matrix from each MoM estimation
    # H: for creating dimensions of G.diff.mat
    # Hvec: for evaluation in G.diff.mat
    
    G.diff.mat <- matrix(0, 7, H)
    for (i in 1:7) { # kappa
        for (j in 2:H) { # Hvec
            cos.Hvec.j <- cos(Hvec[j])
            G.diff.mat[i, j-1] <- (G.mat[j,i] - G.mat[j, i+1]) - 
                eval(parse(text = paste("P.",i - 1, "(", 
                                        cos.Hvec.j, ")", 
                                        sep = ""))) * 
                (G.mat[1,i] - G.mat[1,i + 1]) #)*(1/Hvec[j])
            
            G.diff.mat[i, j-1] <- G.diff.mat[i, j-1]
        }
    }
    return(G.diff.mat)
}



criterion.v01 <- function(G.diff.mat){
    ## Returns S(k) values that are used in estimating kappa
    # G.diff.mat: G difference matrix
    
    criterion <- numeric(7)
    for (i in 1:7) {
        criterion[i] <- t(G.diff.mat[i,]) %*% (G.diff.mat[i,])
    }
    
    return(criterion)
}  









