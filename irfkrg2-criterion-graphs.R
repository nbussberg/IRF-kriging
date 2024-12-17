### Code that creates criterion plots
##  Created by: Nicholas Bussberg
##  Date last updated: 2024-12-16
##  Filename: irfkrg2-criterion-graphs.R
##  Purpose: create criterion plots to estimate kappa for IRFk processes

### Setup ###
rm(list=ls())
options(width = 80)

filename <-"irfkrg2-criterion-graphs"

sink(paste(filename, ".log", sep = "")) # create log file

# R version
R.version.string



### 00 Load functions

source(file = "irf-fcn01-sphere-harmonics.R")
source(file = "irf-fcn02-legendre-polynomials.R")
source(file = "irf-fcn03-distance.R")
source(file = "irf-fcn04-covariances.R")
source(file = "irf-fcn05-MoM-estimator.R")
source(file = "irf-fcn06-kappa-approx.R")


### 01 set inputs and distance bins

H <- 10             # Number Distances at Which to Estimate
eps <- .10          # Use to Bin Distances in Estimation
Hvec <- as.vector(seq(0,1,by = 1/H))[1:H]


###  02 Run loop generating criterion

Lmat <- readRDS("data/irfkrg1-data-r-sig-model-Lmat.Rds")
ud.k2   <- readRDS("data/irfkrg1-data-r-sig-model-ud-k2.Rds")
ud.k3   <- readRDS("data/irfkrg1-data-r-sig-model-ud-k3.Rds")

P <- length(Lmat[,1])
Dvec <- Dfun(Lmat, P)

set.seed(2474)

X0.k2 <- ud.k2 %*% matrix(rnorm(P,0,1),P,1)
X0.k3 <- ud.k3 %*% matrix(rnorm(P,0,1),P,1)

G.k2 <- order.approx(X0.k2, Lmat, Dvec, Hvec, eps, P)$G
G.k3 <- order.approx(X0.k3, Lmat, Dvec, Hvec, eps, P)$G

G.diff.k2 <- G.diff.fcn.noadjust(G.k2, H, Hvec)
G.diff.k3 <- G.diff.fcn.noadjust(G.k3, H, Hvec)
            
(crit.k2 <- criterion.v01(G.diff.k2))
(crit.k3 <- criterion.v01(G.diff.k3))            

# Graphs with original units
pdf(file = paste("graphs/", filename, "-k2.pdf", sep = ""),
    width = 10, height = 7)

par(las=1, mar=c(5,5,2,1)+0.1)
plot(0:6, crit.k2, xlab = "k", ylab = "", cex = 1.5,
     cex.axis = 1.5, cex.lab = 1.5, pch = 16)
mtext("S(k)", 2, line=2.2, cex = 1.5, las=2, at=1070)

dev.off()


pdf(file = paste("graphs/", filename, "-k3.pdf", sep = ""),
    width = 10, height = 7)

par(las=1, mar=c(5,6,2,1)+0.1)
plot(0:6, crit.k3 / 1000000, xlab = "k", ylab = "", cex = 1.5,
     cex.axis = 1.5, cex.lab = 1.5, pch = 16)
mtext("S(k)", 2, line=2.2, cex = 1.5, las=2, at=4)
mtext("(x 10^6)", 2, line=1, cex = 1.5, las=2, at=3.5)

dev.off()
    
        
# Graphs with log units
pdf(file = paste("graphs/", filename, "-k2-log.pdf", sep = ""),
    width = 10, height = 7)

par(las=1, mar=c(5,6,2,1)+0.1)
plot(0:6, log(crit.k2), xlab = "k", ylab = "", cex = 1.5,
     cex.axis = 1.5, cex.lab = 1.5, pch = 16)
mtext("log(S(k))", 2, line=0.5, cex = 1.5, las=2, at=7)

dev.off()


pdf(file = paste("graphs/", filename, "-k3-log.pdf", sep = ""),
    width = 10, height = 7)

par(las=1, mar=c(5,6,2,1)+0.1)
plot(0:6, log(crit.k3), xlab = "k", ylab = "", cex = 1.5,
     cex.axis = 1.5, cex.lab = 1.5, pch = 16)
mtext("log(S(k))", 2, line=0.5, cex = 1.5, las=2, at=16.2)

dev.off()



# close log
sink()


